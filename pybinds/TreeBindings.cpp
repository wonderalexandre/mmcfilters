#include "ModuleBindings.hpp"

#include "PybindConversions.hpp"

#include "../mmcfilters/attributes/AttributeComputation.hpp"
#include "../mmcfilters/attributes/AttributeNames.hpp"
#include "../mmcfilters/trees/saliency/HierarchySaliencyMapValidation.hpp"
#include "../mmcfilters/trees/saliency/HierarchySaliencyMapProjection.hpp"
#include "../mmcfilters/trees/saliency/HierarchySaliencyMap.hpp"
#include "../mmcfilters/trees/saliency/ShapeSpaceSaliency.hpp"
#include "../mmcfilters/trees/MorphologicalTreeFactory.hpp"
#include "../mmcfilters/trees/TreeAltitudeAlgorithms.hpp"
#include "../mmcfilters/trees/WeightedMorphologicalTree.hpp"
#include "../mmcfilters/utils/RegularGridAdjacency2D.hpp"
#include "../mmcfilters/utils/Common.hpp"
#include "../mmcfilters/utils/Image.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <concepts>
#include <cstdint>
#include <memory>
#include <optional>
#include <set>
#include <span>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace mmcfilters::pybindings {

namespace py = pybind11;
using namespace pybind11::literals;

using UInt8InputArray = py::array;

namespace {

template <class Range> std::vector<NodeId> collectNodeIds(const Range& range) {
    std::vector<NodeId> ids;
    for (NodeId id : range) {
        ids.push_back(id);
    }
    return ids;
}

py::array_t<std::uint8_t> reconstructNodeMask(const MorphologicalTree& tree, NodeId nodeId) {
    if (!tree.isNode(nodeId) || !tree.isAlive(nodeId)) {
        throw std::invalid_argument("invalid NodeId for reconstruction");
    }

    ImageUInt8Ptr output = ImageUInt8::create(tree.getNumRowsOfGridDomain2D(), tree.getNumColsOfGridDomain2D());
    output->fill(0);
    for (NodeId properPart : tree.getConnectedComponent(nodeId)) {
        (*output)[properPart] = 255;
    }
    return pybind_utils::toNumpy(output);
}

ImageUInt8Ptr imageFromArray(const UInt8InputArray& input) {
    if (!input.dtype().is(py::dtype::of<uint8_t>())) {
        throw std::invalid_argument("input must be a 2D uint8 array");
    }
    auto buf = input.request();
    if (buf.ndim != 2) {
        throw std::invalid_argument("input must be a 2D uint8 array");
    }
    const int rows = static_cast<int>(buf.shape[0]);
    const int cols = static_cast<int>(buf.shape[1]);
    if (buf.strides[1] != static_cast<py::ssize_t>(sizeof(uint8_t)) || buf.strides[0] != static_cast<py::ssize_t>(cols * sizeof(uint8_t))) {
        throw std::invalid_argument("input must be a C-contiguous 2D uint8 array");
    }
    return ImageUInt8::fromExternal(static_cast<uint8_t*>(buf.ptr), rows, cols);
}

std::uint8_t pythonUInt8AltitudeValue(py::handle value, NodeId nodeId, const char* context) {
    if (py::isinstance<py::bool_>(value) || !PyIndex_Check(value.ptr())) {
        throw std::invalid_argument(std::string(context) + " must be an integer altitude in [0, 255].");
    }

    py::object indexValue = py::reinterpret_steal<py::object>(PyNumber_Index(value.ptr()));
    const long long rawAltitude = py::cast<long long>(indexValue);
    if (rawAltitude < 0 || rawAltitude > 255) {
        std::ostringstream oss;
        oss << context << " requires node altitudes in the uint8 domain [0, 255]; node " << nodeId << " has altitude " << rawAltitude << ".";
        throw std::invalid_argument(oss.str());
    }
    const auto altitude = static_cast<std::uint8_t>(rawAltitude);
    (void)TreeAltitudeAlgorithms::requireUInt8AltitudeValue(altitude, nodeId, context);
    return altitude;
}

std::vector<std::uint8_t> pythonUInt8AltitudeVector(py::handle values, const char* context) {
    if (py::isinstance<py::array>(values)) {
        py::array array = py::reinterpret_borrow<py::array>(values);
        if (!array.dtype().is(py::dtype::of<uint8_t>())) {
            throw std::invalid_argument(std::string(context) + " numpy arrays must have dtype uint8.");
        }
        auto buffer = array.request();
        if (buffer.ndim != 1) {
            throw std::invalid_argument(std::string(context) + " must be a 1D uint8 array.");
        }
        if (buffer.strides[0] != static_cast<py::ssize_t>(sizeof(uint8_t))) {
            throw std::invalid_argument(std::string(context) + " must be a C-contiguous 1D uint8 array.");
        }
        const auto* data = static_cast<const uint8_t*>(buffer.ptr);
        std::vector<std::uint8_t> altitude(static_cast<std::size_t>(buffer.shape[0]));
        for (std::size_t i = 0; i < altitude.size(); ++i) {
            altitude[i] = static_cast<std::uint8_t>(data[i]);
        }
        return altitude;
    }

    if (!py::isinstance<py::sequence>(values)) {
        throw std::invalid_argument(std::string(context) + " must be a sequence of integer altitudes in [0, 255].");
    }

    py::sequence sequence = py::reinterpret_borrow<py::sequence>(values);
    std::vector<std::uint8_t> altitude;
    altitude.reserve(static_cast<std::size_t>(py::len(sequence)));
    NodeId nodeId = 0;
    for (py::handle item : sequence) {
        altitude.push_back(pythonUInt8AltitudeValue(item, nodeId, context));
        ++nodeId;
    }
    return altitude;
}

template <AltitudeValue T> const MorphologicalTree& topology(const WeightedMorphologicalTree<T>& weighted) { return weighted.topology(); }

template <AltitudeValue T> T getAltitudeValue(const WeightedMorphologicalTree<T>& weighted, NodeId nodeId) { return weighted.getAltitude(nodeId); }

template <AltitudeValue T> AltitudeDiff<T> getResidueValue(const WeightedMorphologicalTree<T>& weighted, NodeId nodeId) {
    return weighted.getNodeResidue(nodeId);
}

template <AltitudeValue T> py::array_t<T> reconstructionImageOf(const WeightedMorphologicalTree<T>& weighted) {
    return pybind_utils::toNumpy(weighted.reconstructionImage());
}

template <AltitudeValue T> std::pair<std::vector<NodeId>, std::vector<T>> exportHigraHierarchyOf(const WeightedMorphologicalTree<T>& weighted) {
    return weighted.exportHigraHierarchy();
}

std::vector<Attribute> parseProjectionAttributes(py::handle attributes, int valuesPerNode) {
    std::vector<Attribute> parsed;

    try {
        parsed.push_back(py::cast<Attribute>(attributes));
    } catch (const py::cast_error&) {
        if (!py::isinstance<py::sequence>(attributes)) {
            throw std::invalid_argument("attributes must be an Attribute or a sequence of Attribute values");
        }
        py::sequence seq = py::reinterpret_borrow<py::sequence>(attributes);
        parsed.reserve(static_cast<size_t>(py::len(seq)));
        for (py::handle item : seq) {
            parsed.push_back(py::cast<Attribute>(item));
        }
    }

    if (static_cast<int>(parsed.size()) != valuesPerNode) {
        std::ostringstream oss;
        oss << "attributes must contain " << valuesPerNode << " item(s), got " << parsed.size();
        throw std::invalid_argument(oss.str());
    }
    if (std::set<Attribute>(parsed.begin(), parsed.end()).size() != parsed.size()) {
        throw std::invalid_argument("attributes must not contain duplicates");
    }
    return parsed;
}

AttributeNames makeProjectionAttributeNames(const std::vector<Attribute>& attributes) {
    std::unordered_map<Attribute, int> offsets;
    for (int i = 0; i < static_cast<int>(attributes.size()); ++i) {
        offsets[attributes[static_cast<size_t>(i)]] = i;
    }
    return AttributeNames(std::move(offsets));
}

template <std::floating_point Real>
py::array projectNodeValuesToExportedHigraTyped(const WeightedMorphologicalTree<std::uint8_t>& weighted, const py::array& nodeValues, py::object attributes) {
    const auto buffer = nodeValues.request();
    const int numNodeSlots = weighted.topology().getNumInternalNodeSlots();
    const int numHigraVertices = weighted.topology().getNumTotalProperParts() + weighted.topology().getNumNodes();

    if (buffer.ndim == 1) {
        pybind_utils::require1DArray(buffer, numNodeSlots, "nodeValues");
        const auto parsedAttributes = parseProjectionAttributes(attributes, 1);
        const AttributeNames attrNames = makeProjectionAttributeNames(parsedAttributes);
        auto projected = AttributeComputation::projectNodeValuesToExportedHigra<Real>(
            weighted, attrNames, std::span<const Real>(static_cast<const Real*>(buffer.ptr), static_cast<size_t>(numNodeSlots)));
        return pybind_utils::toNumpyOwned(std::move(projected), numHigraVertices);
    }

    if (buffer.ndim == 2) {
        if (buffer.shape[0] != numNodeSlots) {
            std::ostringstream oss;
            oss << "nodeValues must have " << numNodeSlots << " rows, got " << buffer.shape[0];
            throw std::invalid_argument(oss.str());
        }
        const int valuesPerNode = static_cast<int>(buffer.shape[1]);
        const auto parsedAttributes = parseProjectionAttributes(attributes, valuesPerNode);
        const AttributeNames attrNames = makeProjectionAttributeNames(parsedAttributes);
        auto projected = AttributeComputation::projectNodeValuesToExportedHigra<Real>(
            weighted, attrNames,
            std::span<const Real>(static_cast<const Real*>(buffer.ptr), static_cast<size_t>(numNodeSlots) * static_cast<size_t>(valuesPerNode)));
        return pybind_utils::toNumpyOwned2D(std::move(projected), numHigraVertices, valuesPerNode);
    }

    throw std::invalid_argument("nodeValues must be a 1D or 2D float32 or float64 array");
}

py::array projectNodeValuesToExportedHigraOf(const WeightedMorphologicalTree<std::uint8_t>& weighted, const py::array& nodeValues, py::object attributes) {
    py::array contiguous = py::array::ensure(nodeValues, py::array::c_style);
    if (!contiguous) {
        throw std::invalid_argument("nodeValues must be a 1D or 2D C-contiguous float32 or float64 array");
    }

    if (contiguous.dtype().is(py::dtype::of<float>())) {
        return projectNodeValuesToExportedHigraTyped<float>(weighted, contiguous, std::move(attributes));
    }
    if (contiguous.dtype().is(py::dtype::of<double>())) {
        return projectNodeValuesToExportedHigraTyped<double>(weighted, contiguous, std::move(attributes));
    }

    throw std::invalid_argument("nodeValues must be a 1D or 2D C-contiguous float32 or float64 array");
}

template <class Value> py::dict edgeSaliencyMapToDict(EdgeSaliencyMap<Value>&& edgeMap) {
    const int numEdges = static_cast<int>(edgeMap.size());
    py::dict out;
    out["numRows"] = edgeMap.numRows;
    out["numCols"] = edgeMap.numCols;
    out["adjacencyRadius"] = edgeMap.adjacencyRadius;
    out["sources"] = pybind_utils::toNumpyOwned(std::move(edgeMap.sources), numEdges);
    out["targets"] = pybind_utils::toNumpyOwned(std::move(edgeMap.targets), numEdges);
    out["values"] = pybind_utils::toNumpyOwned(std::move(edgeMap.values), numEdges);
    return out;
}

template <std::floating_point Real> py::dict shapeSpaceExtinctionResultToDict(ShapeSpaceExtinctionResult<Real>&& result) {
    py::list extrema;
    for (const ShapeSpaceExtremum<Real>& extremum : result.extrema) {
        py::dict item;
        item["representative"] = extremum.representative;
        item["birthLevel"] = extremum.birthLevel;
        item["deathLevel"] = extremum.deathLevel;
        item["extinction"] = extremum.extinction;
        extrema.append(std::move(item));
    }

    const int numNodeScores = static_cast<int>(result.nodeScores.size());
    py::dict out;
    out["extrema"] = std::move(extrema);
    out["nodeScores"] = pybind_utils::toNumpyOwned(std::move(result.nodeScores), numNodeScores);
    return out;
}

template <std::floating_point Real> py::dict shapeSpaceSaliencyResultToDict(ShapeSpaceSaliencyResult<Real>&& result) {
    py::list extrema;
    for (const ShapeSpaceExtremum<Real>& extremum : result.extrema) {
        py::dict item;
        item["representative"] = extremum.representative;
        item["birthLevel"] = extremum.birthLevel;
        item["deathLevel"] = extremum.deathLevel;
        item["extinction"] = extremum.extinction;
        extrema.append(std::move(item));
    }

    const int numNodeScores = static_cast<int>(result.nodeScores.size());
    py::dict out;
    out["extrema"] = std::move(extrema);
    out["nodeScores"] = pybind_utils::toNumpyOwned(std::move(result.nodeScores), numNodeScores);
    out["edgeMap"] = edgeSaliencyMapToDict(std::move(result.edgeMap));
    return out;
}

py::dict edgeContourMapToDict(EdgeContourMap&& contourMap) {
    const int numEdges = static_cast<int>(contourMap.size());
    py::dict out;
    out["numRows"] = contourMap.numRows;
    out["numCols"] = contourMap.numCols;
    out["adjacencyRadius"] = contourMap.adjacencyRadius;
    out["sources"] = pybind_utils::toNumpyOwned(std::move(contourMap.sources), numEdges);
    out["targets"] = pybind_utils::toNumpyOwned(std::move(contourMap.targets), numEdges);
    return out;
}

py::dict nodeContourEdgeMapToDict(NodeContourEdgeMap&& contourMap) {
    const int numEdges = static_cast<int>(contourMap.size());
    py::dict out;
    out["numRows"] = contourMap.numRows;
    out["numCols"] = contourMap.numCols;
    out["adjacencyRadius"] = contourMap.adjacencyRadius;
    out["sources"] = pybind_utils::toNumpyOwned(std::move(contourMap.sources), numEdges);
    out["targets"] = pybind_utils::toNumpyOwned(std::move(contourMap.targets), numEdges);
    out["nodes"] = pybind_utils::toNumpyOwned(std::move(contourMap.nodes), numEdges);
    return out;
}

py::dict incrementalNodeContourMapToDict(IncrementalNodeContourMap&& contourMap) {
    const int numEdges = static_cast<int>(contourMap.size());
    const int numOffsets = static_cast<int>(contourMap.offsets.size());
    py::dict out;
    out["numRows"] = contourMap.numRows;
    out["numCols"] = contourMap.numCols;
    out["numNodeSlots"] = contourMap.numNodeSlots;
    out["adjacencyRadius"] = contourMap.adjacencyRadius;
    out["offsets"] = pybind_utils::toNumpyOwned(std::move(contourMap.offsets), numOffsets);
    out["sources"] = pybind_utils::toNumpyOwned(std::move(contourMap.sources), numEdges);
    out["targets"] = pybind_utils::toNumpyOwned(std::move(contourMap.targets), numEdges);
    return out;
}

EdgeSaliencyMap<double> edgeSaliencyMapFromDictPy(py::dict edgeMap, const char* context) {
    for (const char* key : {"numRows", "numCols", "adjacencyRadius", "sources", "targets", "values"}) {
        if (!edgeMap.contains(key)) {
            throw std::invalid_argument(std::string(context) + " edgeMap is missing key '" + key + "'.");
        }
    }

    auto sources = py::array_t<NodeId, py::array::c_style | py::array::forcecast>::ensure(edgeMap["sources"]);
    auto targets = py::array_t<NodeId, py::array::c_style | py::array::forcecast>::ensure(edgeMap["targets"]);
    auto values = py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(edgeMap["values"]);
    if (!sources || !targets || !values) {
        throw std::invalid_argument(std::string(context) + " requires 1D numeric sources, targets, and values arrays.");
    }

    const py::buffer_info sourceInfo = sources.request();
    const py::buffer_info targetInfo = targets.request();
    const py::buffer_info valueInfo = values.request();
    if (sourceInfo.ndim != 1 || targetInfo.ndim != 1 || valueInfo.ndim != 1) {
        throw std::invalid_argument(std::string(context) + " requires 1D sources, targets, and values arrays.");
    }
    if (sourceInfo.shape[0] != targetInfo.shape[0] || sourceInfo.shape[0] != valueInfo.shape[0]) {
        throw std::invalid_argument(std::string(context) + " requires sources, targets, and values arrays with the same length.");
    }

    const auto numEdges = static_cast<std::size_t>(sourceInfo.shape[0]);
    const auto* sourceData = static_cast<const NodeId*>(sourceInfo.ptr);
    const auto* targetData = static_cast<const NodeId*>(targetInfo.ptr);
    const auto* valueData = static_cast<const double*>(valueInfo.ptr);

    EdgeSaliencyMap<double> parsed;
    parsed.numRows = py::cast<int>(edgeMap["numRows"]);
    parsed.numCols = py::cast<int>(edgeMap["numCols"]);
    parsed.adjacencyRadius = py::cast<double>(edgeMap["adjacencyRadius"]);
    parsed.sources.assign(sourceData, sourceData + numEdges);
    parsed.targets.assign(targetData, targetData + numEdges);
    parsed.values.assign(valueData, valueData + numEdges);
    return parsed;
}

HierarchyValuationPolicy saliencyPolicyFromStrict(bool strict) {
    return strict ? HierarchyValuationPolicy::RequireStrictHierarchy : HierarchyValuationPolicy::AllowLevelCollapse;
}

HierarchyValuationRangePolicy saliencyRangePolicyFromNonnegative(bool nonnegative) {
    return nonnegative ? HierarchyValuationRangePolicy::RequireNonNegative : HierarchyValuationRangePolicy::AllowAnyFinite;
}

template <class Real>
void validateHierarchyValuationTypedPy(const WeightedMorphologicalTree<std::uint8_t>& weighted, py::array valuation, bool strict, bool nonnegative) {
    auto scores = pybind_utils::requireNodeAttributeArray<Real>(std::move(valuation), weighted.topology(), "valuation");
    const py::buffer_info buffer = scores.request();
    HierarchySaliencyMapValidation::validateHierarchyValuation(
        weighted.topology(), std::span<const Real>(static_cast<const Real*>(buffer.ptr), static_cast<std::size_t>(buffer.shape[0])),
        saliencyPolicyFromStrict(strict), saliencyRangePolicyFromNonnegative(nonnegative), "HierarchySaliencyMapValidation.validateHierarchyValuation");
}

void validateHierarchyValuationPy(const WeightedMorphologicalTree<std::uint8_t>& weighted, py::array valuation, bool strict, bool nonnegative) {
    py::array contiguous = py::array::ensure(valuation, py::array::c_style);
    if (!contiguous) {
        throw std::invalid_argument("valuation must be a 1D C-contiguous int32, float32, or float64 array");
    }

    if (contiguous.dtype().is(py::dtype::of<int>())) {
        validateHierarchyValuationTypedPy<int>(weighted, contiguous, strict, nonnegative);
        return;
    }
    if (contiguous.dtype().is(py::dtype::of<float>())) {
        validateHierarchyValuationTypedPy<float>(weighted, contiguous, strict, nonnegative);
        return;
    }
    if (contiguous.dtype().is(py::dtype::of<double>())) {
        validateHierarchyValuationTypedPy<double>(weighted, contiguous, strict, nonnegative);
        return;
    }

    throw std::invalid_argument("valuation must be a 1D C-contiguous int32, float32, or float64 array");
}

template <class Real> py::array rankHierarchyValuationTypedPy(const WeightedMorphologicalTree<std::uint8_t>& weighted, py::array valuation, bool strict) {
    auto scores = pybind_utils::requireNodeAttributeArray<Real>(std::move(valuation), weighted.topology(), "valuation");
    const py::buffer_info buffer = scores.request();
    std::vector<int> ranks = HierarchySaliencyMapValidation::rankHierarchyValuation(
        weighted.topology(), std::span<const Real>(static_cast<const Real*>(buffer.ptr), static_cast<std::size_t>(buffer.shape[0])),
        saliencyPolicyFromStrict(strict));
    const int numRanks = static_cast<int>(ranks.size());
    return pybind_utils::toNumpyOwned(std::move(ranks), numRanks);
}

py::array rankHierarchyValuationPy(const WeightedMorphologicalTree<std::uint8_t>& weighted, py::array valuation, bool strict) {
    py::array contiguous = py::array::ensure(valuation, py::array::c_style);
    if (!contiguous) {
        throw std::invalid_argument("valuation must be a 1D C-contiguous int32, float32, or float64 array");
    }

    if (contiguous.dtype().is(py::dtype::of<int>())) {
        return rankHierarchyValuationTypedPy<int>(weighted, contiguous, strict);
    }
    if (contiguous.dtype().is(py::dtype::of<float>())) {
        return rankHierarchyValuationTypedPy<float>(weighted, contiguous, strict);
    }
    if (contiguous.dtype().is(py::dtype::of<double>())) {
        return rankHierarchyValuationTypedPy<double>(weighted, contiguous, strict);
    }

    throw std::invalid_argument("valuation must be a 1D C-contiguous int32, float32, or float64 array");
}

template <class Real>
py::array computeNormalizedScoresTypedPy(const WeightedMorphologicalTree<std::uint8_t>& weighted, py::array valuation, bool strict, bool nonnegative) {
    auto scores = pybind_utils::requireNodeAttributeArray<Real>(std::move(valuation), weighted.topology(), "valuation");
    const py::buffer_info buffer = scores.request();
    std::vector<double> normalized = HierarchySaliencyMapValidation::computeNormalizedScores(
        weighted.topology(), std::span<const Real>(static_cast<const Real*>(buffer.ptr), static_cast<std::size_t>(buffer.shape[0])),
        saliencyPolicyFromStrict(strict), saliencyRangePolicyFromNonnegative(nonnegative));
    const int numScores = static_cast<int>(normalized.size());
    return pybind_utils::toNumpyOwned(std::move(normalized), numScores);
}

py::array computeNormalizedScoresPy(const WeightedMorphologicalTree<std::uint8_t>& weighted, py::array valuation, bool strict, bool nonnegative) {
    py::array contiguous = py::array::ensure(valuation, py::array::c_style);
    if (!contiguous) {
        throw std::invalid_argument("valuation must be a 1D C-contiguous int32, float32, or float64 array");
    }

    if (contiguous.dtype().is(py::dtype::of<int>())) {
        return computeNormalizedScoresTypedPy<int>(weighted, contiguous, strict, nonnegative);
    }
    if (contiguous.dtype().is(py::dtype::of<float>())) {
        return computeNormalizedScoresTypedPy<float>(weighted, contiguous, strict, nonnegative);
    }
    if (contiguous.dtype().is(py::dtype::of<double>())) {
        return computeNormalizedScoresTypedPy<double>(weighted, contiguous, strict, nonnegative);
    }

    throw std::invalid_argument("valuation must be a 1D C-contiguous int32, float32, or float64 array");
}

IncrementalNodeContourMap incrementalNodeContourMapFromDictPy(py::dict contours) {
    constexpr const char* context = "HierarchySaliencyMapProjection incremental contour operation";
    for (const char* key : {"numRows", "numCols", "numNodeSlots", "adjacencyRadius", "offsets", "sources", "targets"}) {
        if (!contours.contains(key)) {
            throw std::invalid_argument(std::string(context) + " contours is missing key '" + key + "'.");
        }
    }

    auto offsets = py::array_t<std::size_t, py::array::c_style | py::array::forcecast>::ensure(contours["offsets"]);
    auto sources = py::array_t<NodeId, py::array::c_style | py::array::forcecast>::ensure(contours["sources"]);
    auto targets = py::array_t<NodeId, py::array::c_style | py::array::forcecast>::ensure(contours["targets"]);
    if (!offsets || !sources || !targets) {
        throw std::invalid_argument(std::string(context) + " requires 1D numeric offsets, sources, and targets arrays.");
    }

    const py::buffer_info offsetInfo = offsets.request();
    const py::buffer_info sourceInfo = sources.request();
    const py::buffer_info targetInfo = targets.request();
    if (offsetInfo.ndim != 1 || sourceInfo.ndim != 1 || targetInfo.ndim != 1) {
        throw std::invalid_argument(std::string(context) + " requires 1D offsets, sources, and targets arrays.");
    }
    if (sourceInfo.shape[0] != targetInfo.shape[0]) {
        throw std::invalid_argument(std::string(context) + " requires sources and targets arrays with the same length.");
    }

    const auto numOffsets = static_cast<std::size_t>(offsetInfo.shape[0]);
    const auto numEdges = static_cast<std::size_t>(sourceInfo.shape[0]);
    const auto* offsetData = static_cast<const std::size_t*>(offsetInfo.ptr);
    const auto* sourceData = static_cast<const NodeId*>(sourceInfo.ptr);
    const auto* targetData = static_cast<const NodeId*>(targetInfo.ptr);

    IncrementalNodeContourMap parsed;
    parsed.numRows = py::cast<int>(contours["numRows"]);
    parsed.numCols = py::cast<int>(contours["numCols"]);
    parsed.numNodeSlots = py::cast<int>(contours["numNodeSlots"]);
    parsed.adjacencyRadius = py::cast<double>(contours["adjacencyRadius"]);
    parsed.offsets.assign(offsetData, offsetData + numOffsets);
    parsed.sources.assign(sourceData, sourceData + numEdges);
    parsed.targets.assign(targetData, targetData + numEdges);
    return parsed;
}

py::dict computeTopologicalLevelEdgeMapPy(const WeightedMorphologicalTree<std::uint8_t>& weighted, std::optional<double> radius) {
    if (radius.has_value()) {
        RegularGridAdjacency2D adjacency =
            pybind_utils::makeRegularGridAdjacency2D(weighted.topology().getNumRowsOfGridDomain2D(), weighted.topology().getNumColsOfGridDomain2D(), *radius,
                                                     "HierarchySaliencyMap.computeTopologicalLevelEdgeMap");
        return edgeSaliencyMapToDict(HierarchySaliencyMap::computeTopologicalLevelEdgeMap(weighted.topology(), adjacency));
    }
    return edgeSaliencyMapToDict(HierarchySaliencyMap::computeTopologicalLevelEdgeMap(weighted.topology()));
}

py::dict computeNormalizedAltitudeEdgeMapPy(const WeightedMorphologicalTree<std::uint8_t>& weighted, std::optional<double> radius) {
    if (radius.has_value()) {
        RegularGridAdjacency2D adjacency =
            pybind_utils::makeRegularGridAdjacency2D(weighted.topology().getNumRowsOfGridDomain2D(), weighted.topology().getNumColsOfGridDomain2D(), *radius,
                                                     "HierarchySaliencyMap.computeNormalizedAltitudeEdgeMap");
        return edgeSaliencyMapToDict(HierarchySaliencyMap::computeNormalizedAltitudeEdgeMap(weighted, adjacency));
    }
    return edgeSaliencyMapToDict(HierarchySaliencyMap::computeNormalizedAltitudeEdgeMap(weighted));
}

RegularGridAdjacency2D saliencyAdjacencyPy(const WeightedMorphologicalTree<std::uint8_t>& weighted, std::optional<double> radius, const char* context) {
    if (radius.has_value()) {
        return pybind_utils::makeRegularGridAdjacency2D(weighted.topology().getNumRowsOfGridDomain2D(), weighted.topology().getNumColsOfGridDomain2D(), *radius,
                                                        context);
    }
    return HierarchySaliencyMap::requireProjectionAdjacency(weighted.topology(), context);
}

template <class Real>
py::dict computeSaliencyEdgeMapTypedPy(const WeightedMorphologicalTree<std::uint8_t>& weighted, py::array valuation, std::optional<double> radius, bool strict,
                                       HierarchyLevelConvention levelConvention, bool validateConnectivity, const char* context) {
    auto scores = pybind_utils::requireNodeAttributeArray<Real>(std::move(valuation), weighted.topology(), "valuation");
    const py::buffer_info buffer = scores.request();
    const auto* scoreData = static_cast<const Real*>(buffer.ptr);
    const auto valuationSpan = std::span<const Real>(scoreData, static_cast<std::size_t>(buffer.shape[0]));

    if (radius.has_value()) {
        RegularGridAdjacency2D adjacency = pybind_utils::makeRegularGridAdjacency2D(weighted.topology().getNumRowsOfGridDomain2D(),
                                                                                    weighted.topology().getNumColsOfGridDomain2D(), *radius, context);
        return edgeSaliencyMapToDict(HierarchySaliencyMap::computeSaliencyEdgeMap(
            weighted.topology(), adjacency, valuationSpan, saliencyPolicyFromStrict(strict), levelConvention,
            validateConnectivity ? HierarchyConnectivityPolicy::ValidateConnected : HierarchyConnectivityPolicy::AssumeConnected));
    }
    return edgeSaliencyMapToDict(HierarchySaliencyMap::computeSaliencyEdgeMap(
        weighted.topology(), valuationSpan, saliencyPolicyFromStrict(strict), levelConvention,
        validateConnectivity ? HierarchyConnectivityPolicy::ValidateConnected : HierarchyConnectivityPolicy::AssumeConnected));
}

py::dict computeSaliencyEdgeMapPy(const WeightedMorphologicalTree<std::uint8_t>& weighted, py::array valuation, std::optional<double> radius, bool strict,
                                  HierarchyLevelConvention levelConvention, bool validateConnectivity, const char* context) {
    py::array contiguous = py::array::ensure(valuation, py::array::c_style);
    if (!contiguous) {
        throw std::invalid_argument("valuation must be a 1D C-contiguous int32, float32, or float64 array");
    }

    if (contiguous.dtype().is(py::dtype::of<int>())) {
        return computeSaliencyEdgeMapTypedPy<int>(weighted, contiguous, radius, strict, levelConvention, validateConnectivity, context);
    }
    if (contiguous.dtype().is(py::dtype::of<float>())) {
        return computeSaliencyEdgeMapTypedPy<float>(weighted, contiguous, radius, strict, levelConvention, validateConnectivity, context);
    }
    if (contiguous.dtype().is(py::dtype::of<double>())) {
        return computeSaliencyEdgeMapTypedPy<double>(weighted, contiguous, radius, strict, levelConvention, validateConnectivity, context);
    }

    throw std::invalid_argument("valuation must be a 1D C-contiguous int32, float32, or float64 array");
}

template <class Real>
py::dict computeCanonicalRankedSaliencyEdgeMapTypedPy(const WeightedMorphologicalTree<std::uint8_t>& weighted, py::array valuation,
                                                      std::optional<double> radius, bool strict, bool validateConnectivity) {
    auto scores = pybind_utils::requireNodeAttributeArray<Real>(std::move(valuation), weighted.topology(), "valuation");
    const py::buffer_info buffer = scores.request();
    const auto valuationSpan = std::span<const Real>(static_cast<const Real*>(buffer.ptr), static_cast<std::size_t>(buffer.shape[0]));
    const RegularGridAdjacency2D adjacency = saliencyAdjacencyPy(weighted, radius, "HierarchySaliencyMap.computeCanonicalRankedSaliencyEdgeMap");
    return edgeSaliencyMapToDict(HierarchySaliencyMap::computeCanonicalRankedSaliencyEdgeMap(
        weighted.topology(), adjacency, valuationSpan, saliencyPolicyFromStrict(strict),
        validateConnectivity ? HierarchyConnectivityPolicy::ValidateConnected : HierarchyConnectivityPolicy::AssumeConnected));
}

py::dict computeCanonicalRankedSaliencyEdgeMapPy(const WeightedMorphologicalTree<std::uint8_t>& weighted, py::array valuation, std::optional<double> radius,
                                                 bool strict, bool validateConnectivity) {
    py::array contiguous = py::array::ensure(valuation, py::array::c_style);
    if (!contiguous) {
        throw std::invalid_argument("valuation must be a 1D C-contiguous int32, float32, or float64 array");
    }
    if (contiguous.dtype().is(py::dtype::of<int>())) {
        return computeCanonicalRankedSaliencyEdgeMapTypedPy<int>(weighted, contiguous, radius, strict, validateConnectivity);
    }
    if (contiguous.dtype().is(py::dtype::of<float>())) {
        return computeCanonicalRankedSaliencyEdgeMapTypedPy<float>(weighted, contiguous, radius, strict, validateConnectivity);
    }
    if (contiguous.dtype().is(py::dtype::of<double>())) {
        return computeCanonicalRankedSaliencyEdgeMapTypedPy<double>(weighted, contiguous, radius, strict, validateConnectivity);
    }
    throw std::invalid_argument("valuation must be a 1D C-contiguous int32, float32, or float64 array");
}

template <std::floating_point Real>
std::span<const Real> shapeSpaceDenseNodeSpan(py::array values, const MorphologicalTree& tree, py::array_t<Real, py::array::c_style>& owner,
                                              const char* argumentName) {
    owner = pybind_utils::requireNodeAttributeArray<Real>(std::move(values), tree, argumentName);
    const py::buffer_info buffer = owner.request();
    return std::span<const Real>(static_cast<const Real*>(buffer.ptr), static_cast<std::size_t>(buffer.shape[0]));
}

template <std::floating_point Real>
py::dict computeShapeSpaceExtinctionValuesTypedPy(const WeightedMorphologicalTree<std::uint8_t>& weighted, py::array attribute,
                                                  ShapeSpaceExtremaPolarity polarity) {
    py::array_t<Real, py::array::c_style> attributeOwner;
    const std::span<const Real> attributeSpan = shapeSpaceDenseNodeSpan<Real>(std::move(attribute), weighted.topology(), attributeOwner, "attribute");
    return shapeSpaceExtinctionResultToDict(ShapeSpaceSaliency::computeExtinctionValues(weighted.topology(), attributeSpan, polarity));
}

py::dict computeShapeSpaceExtinctionValuesPy(const WeightedMorphologicalTree<std::uint8_t>& weighted, py::array attribute, ShapeSpaceExtremaPolarity polarity) {
    py::array contiguous = py::array::ensure(attribute, py::array::c_style);
    if (!contiguous) {
        throw std::invalid_argument("attribute must be a 1D C-contiguous float32 or float64 array");
    }
    if (contiguous.dtype().is(py::dtype::of<float>())) {
        return computeShapeSpaceExtinctionValuesTypedPy<float>(weighted, contiguous, polarity);
    }
    if (contiguous.dtype().is(py::dtype::of<double>())) {
        return computeShapeSpaceExtinctionValuesTypedPy<double>(weighted, contiguous, polarity);
    }
    throw std::invalid_argument("attribute must be a 1D C-contiguous float32 or float64 array");
}

template <std::floating_point Real>
py::dict projectShapeSpaceContourScoresTypedPy(const WeightedMorphologicalTree<std::uint8_t>& weighted, py::array scores, std::optional<double> radius) {
    py::array_t<Real, py::array::c_style> scoreOwner;
    const std::span<const Real> scoreSpan = shapeSpaceDenseNodeSpan<Real>(std::move(scores), weighted.topology(), scoreOwner, "nodeScores");

    if (radius.has_value()) {
        RegularGridAdjacency2D adjacency = pybind_utils::makeRegularGridAdjacency2D(
            weighted.topology().getNumRowsOfGridDomain2D(), weighted.topology().getNumColsOfGridDomain2D(), *radius, "ShapeSpaceSaliency.projectContourScores");
        return edgeSaliencyMapToDict(ShapeSpaceSaliency::projectContourScores(weighted.topology(), scoreSpan, adjacency));
    }
    return edgeSaliencyMapToDict(ShapeSpaceSaliency::projectContourScores(weighted.topology(), scoreSpan));
}

py::dict projectShapeSpaceContourScoresPy(const WeightedMorphologicalTree<std::uint8_t>& weighted, py::array scores, std::optional<double> radius) {
    py::array contiguous = py::array::ensure(scores, py::array::c_style);
    if (!contiguous) {
        throw std::invalid_argument("nodeScores must be a 1D C-contiguous float32 or float64 array");
    }
    if (contiguous.dtype().is(py::dtype::of<float>())) {
        return projectShapeSpaceContourScoresTypedPy<float>(weighted, contiguous, radius);
    }
    if (contiguous.dtype().is(py::dtype::of<double>())) {
        return projectShapeSpaceContourScoresTypedPy<double>(weighted, contiguous, radius);
    }
    throw std::invalid_argument("nodeScores must be a 1D C-contiguous float32 or float64 array");
}

template <std::floating_point Real>
py::dict computeShapeSpaceSaliencyTypedPy(const WeightedMorphologicalTree<std::uint8_t>& weighted, py::array attribute, ShapeSpaceExtremaPolarity polarity,
                                          std::optional<double> radius) {
    py::array_t<Real, py::array::c_style> attributeOwner;
    const std::span<const Real> attributeSpan = shapeSpaceDenseNodeSpan<Real>(std::move(attribute), weighted.topology(), attributeOwner, "attribute");

    if (radius.has_value()) {
        RegularGridAdjacency2D adjacency = pybind_utils::makeRegularGridAdjacency2D(
            weighted.topology().getNumRowsOfGridDomain2D(), weighted.topology().getNumColsOfGridDomain2D(), *radius, "ShapeSpaceSaliency.compute");
        return shapeSpaceSaliencyResultToDict(ShapeSpaceSaliency::compute(weighted.topology(), attributeSpan, polarity, adjacency));
    }
    return shapeSpaceSaliencyResultToDict(ShapeSpaceSaliency::compute(weighted.topology(), attributeSpan, polarity));
}

py::dict computeShapeSpaceSaliencyPy(const WeightedMorphologicalTree<std::uint8_t>& weighted, py::array attribute, ShapeSpaceExtremaPolarity polarity,
                                     std::optional<double> radius) {
    py::array contiguous = py::array::ensure(attribute, py::array::c_style);
    if (!contiguous) {
        throw std::invalid_argument("attribute must be a 1D C-contiguous float32 or float64 array");
    }
    if (contiguous.dtype().is(py::dtype::of<float>())) {
        return computeShapeSpaceSaliencyTypedPy<float>(weighted, contiguous, polarity, radius);
    }
    if (contiguous.dtype().is(py::dtype::of<double>())) {
        return computeShapeSpaceSaliencyTypedPy<double>(weighted, contiguous, polarity, radius);
    }
    throw std::invalid_argument("attribute must be a 1D C-contiguous float32 or float64 array");
}

py::dict thresholdCutPy(py::dict edgeMap, double threshold) {
    EdgeSaliencyMap<double> parsed = edgeSaliencyMapFromDictPy(std::move(edgeMap), "HierarchySaliencyMapProjection.thresholdCut");
    return edgeContourMapToDict(HierarchySaliencyMapProjection::thresholdCut(parsed, threshold));
}

py::array edgeMapToPixelImagePy(py::dict edgeMap, EdgeToPixelReducer reducer) {
    EdgeSaliencyMap<double> parsed = edgeSaliencyMapFromDictPy(std::move(edgeMap), "HierarchySaliencyMapProjection.edgeMapToPixelImage");
    return pybind_utils::toNumpy(HierarchySaliencyMapProjection::edgeMapToPixelImage(parsed, reducer));
}

py::dict nodeContourEdgesPy(const WeightedMorphologicalTree<std::uint8_t>& weighted, std::optional<double> radius) {
    if (radius.has_value()) {
        RegularGridAdjacency2D adjacency =
            pybind_utils::makeRegularGridAdjacency2D(weighted.topology().getNumRowsOfGridDomain2D(), weighted.topology().getNumColsOfGridDomain2D(), *radius,
                                                     "HierarchySaliencyMapProjection.nodeContourEdges");
        return nodeContourEdgeMapToDict(HierarchySaliencyMapProjection::nodeContourEdges(weighted, adjacency));
    }
    return nodeContourEdgeMapToDict(HierarchySaliencyMapProjection::nodeContourEdges(weighted));
}

py::dict computeIncrementalNodeContoursPy(const WeightedMorphologicalTree<std::uint8_t>& weighted, std::optional<double> radius) {
    if (radius.has_value()) {
        RegularGridAdjacency2D adjacency =
            pybind_utils::makeRegularGridAdjacency2D(weighted.topology().getNumRowsOfGridDomain2D(), weighted.topology().getNumColsOfGridDomain2D(), *radius,
                                                     "HierarchySaliencyMapProjection.computeIncrementalNodeContours");
        return incrementalNodeContourMapToDict(HierarchySaliencyMapProjection::computeIncrementalNodeContours(weighted, adjacency));
    }
    return incrementalNodeContourMapToDict(HierarchySaliencyMapProjection::computeIncrementalNodeContours(weighted));
}

template <std::floating_point Real> py::dict projectNodeValuationOnIncrementalContoursTypedPy(py::dict contours, py::array nodeValuation) {
    IncrementalNodeContourMap parsed = incrementalNodeContourMapFromDictPy(std::move(contours));
    auto valuation = py::array_t<Real, py::array::c_style>::ensure(nodeValuation);
    if (!valuation) {
        throw std::invalid_argument("nodeValuation must be a 1D C-contiguous float32 or float64 array.");
    }
    const py::buffer_info info = valuation.request();
    pybind_utils::require1DArray(info, parsed.numNodeSlots, "nodeValuation");
    if (info.strides[0] != static_cast<py::ssize_t>(sizeof(Real))) {
        throw std::invalid_argument("nodeValuation must be C-contiguous.");
    }
    return edgeSaliencyMapToDict(HierarchySaliencyMapProjection::projectNodeValuation(
        parsed, std::span<const Real>(static_cast<const Real*>(info.ptr), static_cast<std::size_t>(info.shape[0]))));
}

py::dict projectNodeValuationOnIncrementalContoursPy(py::dict contours, py::array nodeValuation) {
    py::array contiguous = py::array::ensure(nodeValuation, py::array::c_style);
    if (!contiguous) {
        throw std::invalid_argument("nodeValuation must be a 1D C-contiguous float32 or float64 array");
    }
    if (contiguous.dtype().is(py::dtype::of<float>())) {
        return projectNodeValuationOnIncrementalContoursTypedPy<float>(std::move(contours), contiguous);
    }
    if (contiguous.dtype().is(py::dtype::of<double>())) {
        return projectNodeValuationOnIncrementalContoursTypedPy<double>(std::move(contours), contiguous);
    }
    throw std::invalid_argument("nodeValuation must be a 1D C-contiguous float32 or float64 array");
}

template <std::floating_point Real> py::dict thresholdIncrementalContoursByNodeValuationTypedPy(py::dict contours, py::array nodeValuation, double threshold) {
    IncrementalNodeContourMap parsed = incrementalNodeContourMapFromDictPy(std::move(contours));
    auto valuation = py::array_t<Real, py::array::c_style>::ensure(nodeValuation);
    if (!valuation) {
        throw std::invalid_argument("nodeValuation must be a 1D C-contiguous float32 or float64 array.");
    }
    const py::buffer_info info = valuation.request();
    pybind_utils::require1DArray(info, parsed.numNodeSlots, "nodeValuation");
    if (info.strides[0] != static_cast<py::ssize_t>(sizeof(Real))) {
        throw std::invalid_argument("nodeValuation must be C-contiguous.");
    }
    return edgeContourMapToDict(HierarchySaliencyMapProjection::thresholdByNodeValuation(
        parsed, std::span<const Real>(static_cast<const Real*>(info.ptr), static_cast<std::size_t>(info.shape[0])), threshold));
}

py::dict thresholdIncrementalContoursByNodeValuationPy(py::dict contours, py::array nodeValuation, double threshold) {
    py::array contiguous = py::array::ensure(nodeValuation, py::array::c_style);
    if (!contiguous) {
        throw std::invalid_argument("nodeValuation must be a 1D C-contiguous float32 or float64 array");
    }
    if (contiguous.dtype().is(py::dtype::of<float>())) {
        return thresholdIncrementalContoursByNodeValuationTypedPy<float>(std::move(contours), contiguous, threshold);
    }
    if (contiguous.dtype().is(py::dtype::of<double>())) {
        return thresholdIncrementalContoursByNodeValuationTypedPy<double>(std::move(contours), contiguous, threshold);
    }
    throw std::invalid_argument("nodeValuation must be a 1D C-contiguous float32 or float64 array");
}

template <class TreeLike, class PyClass> void bindTreeQueryApi(PyClass& cls) {
    cls.def_property_readonly(
           "numInternalNodeSlots", [](TreeLike& self) { return topology(self).getNumInternalNodeSlots(); }, "Size of the dense internal NodeId slot domain.")
        .def_property_readonly(
            "numTotalProperParts", [](TreeLike& self) { return topology(self).getNumTotalProperParts(); },
            "Number of proper parts in the finite support domain.")
        .def_property_readonly(
            "numHigraNodes", [](TreeLike& self) { return topology(self).getNumHigraNodes(); }, "Size of the preserved imported Higra node-id domain.")
        .def(
            "getRoot", [](TreeLike& self) { return topology(self).getRoot(); }, "Return the current root node id.")
        .def_property_readonly(
            "root", [](TreeLike& self) { return topology(self).getRoot(); }, "Current root node id.")
        .def_property_readonly(
            "numFreeNodeSlots", [](TreeLike& self) { return topology(self).getNumFreeNodeSlots(); }, "Number of currently free internal-node slots.")
        .def_property_readonly(
            "numLeafNodes", [](TreeLike& self) { return topology(self).getNumLeafNodes(); }, "Number of alive leaf nodes.")
        .def(
            "getAliveNodeIds", [](TreeLike& self) { return collectNodeIds(topology(self).getAliveNodeIds()); },
            "Return all alive internal-node ids in the dense node-id domain.")
        .def_property_readonly(
            "aliveNodeIds", [](TreeLike& self) { return collectNodeIds(topology(self).getAliveNodeIds()); },
            "Alive internal-node ids in the dense node-id domain.")
        .def_property_readonly(
            "alive_node_ids", [](TreeLike& self) { return collectNodeIds(topology(self).getAliveNodeIds()); },
            "Alive internal-node ids in the dense node-id domain.")
        .def(
            "getLeafNodeIds", [](TreeLike& self) { return topology(self).getLeaves(); }, "Return alive leaf node ids in the dense node-id domain.")
        .def_property_readonly(
            "leafNodeIds", [](TreeLike& self) { return topology(self).getLeaves(); }, "Alive leaf node ids in the dense node-id domain.")
        .def_property_readonly(
            "leaf_node_ids", [](TreeLike& self) { return topology(self).getLeaves(); }, "Alive leaf node ids in the dense node-id domain.")
        .def(
            "getChildren", [](TreeLike& self, NodeId nodeId) { return collectNodeIds(topology(self).getChildren(nodeId)); }, "nodeId"_a,
            "Return the direct children of a node in the dense node-id domain.")
        .def(
            "getNodeNumDescendants", [](TreeLike& self, NodeId nodeId) { return topology(self).getNodeNumDescendants(nodeId); }, "nodeId"_a,
            "Return the number of descendants of nodeId.")
        .def(
            "getNodeNumSiblings", [](TreeLike& self, NodeId nodeId) { return topology(self).getNodeNumSiblings(nodeId); }, "nodeId"_a,
            "Return the number of siblings of nodeId.")
        .def(
            "getNumProperParts", [](TreeLike& self, NodeId nodeId) { return topology(self).getNumProperParts(nodeId); }, "nodeId"_a,
            "Return the number of direct proper parts owned by nodeId.")
        .def(
            "getNodeTimePreOrder", [](TreeLike& self, NodeId nodeId) { return topology(self).getNodeTimePreOrder(nodeId); }, "nodeId"_a,
            "Return the preorder timestamp of nodeId.")
        .def(
            "getNodeTimePostOrder", [](TreeLike& self, NodeId nodeId) { return topology(self).getNodeTimePostOrder(nodeId); }, "nodeId"_a,
            "Return the postorder timestamp of nodeId.")
        .def(
            "getProperParts", [](TreeLike& self, NodeId nodeId) { return collectNodeIds(topology(self).getProperParts(nodeId)); }, "nodeId"_a,
            "Return the proper parts owned directly by a node.")
        .def(
            "getConnectedComponent",
            [](TreeLike& self, NodeId nodeId) {
                auto range = topology(self).getConnectedComponent(nodeId);
                return py::make_iterator(range.begin(), range.end());
            },
            py::keep_alive<0, 1>(), "nodeId"_a, "Iterate over all proper parts in the connected component represented by nodeId.")
        .def(
            "reconstructNode", [](TreeLike& self, NodeId nodeId) { return reconstructNodeMask(topology(self), nodeId); }, "nodeId"_a,
            "Reconstruct a binary mask for the connected component represented by nodeId.")
        .def(
            "getPostOrderNodes",
            [](TreeLike& self, std::optional<NodeId> rootNodeId) {
                return rootNodeId.has_value() ? collectNodeIds(topology(self).getPostOrderNodes(*rootNodeId))
                                              : collectNodeIds(topology(self).getPostOrderNodes());
            },
            "rootNodeId"_a = std::nullopt, "Return post-order traversal node ids under `rootNodeId`, or under the tree root.")
        .def(
            "getIteratorBreadthFirstTraversal",
            [](TreeLike& self, std::optional<NodeId> rootNodeId) {
                return rootNodeId.has_value() ? collectNodeIds(topology(self).getIteratorBreadthFirstTraversal(*rootNodeId))
                                              : collectNodeIds(topology(self).getIteratorBreadthFirstTraversal());
            },
            "rootNodeId"_a = std::nullopt, "Return breadth-first traversal node ids under `rootNodeId`, or under the tree root.")
        .def(
            "getPathToRootNodes", [](TreeLike& self, NodeId nodeId) { return collectNodeIds(topology(self).getPathToRootNodes(nodeId)); }, "nodeId"_a,
            "Return the path from `nodeId` to the current root.")
        .def(
            "getPathBetweenNodes",
            [](TreeLike& self, NodeId sourceNodeId, NodeId targetNodeId) {
                return collectNodeIds(topology(self).getPathBetweenNodes(sourceNodeId, targetNodeId));
            },
            "sourceNodeId"_a, "targetNodeId"_a, "Return the upward path from `sourceNodeId` toward `targetNodeId`.")
        .def(
            "getNodeSubtree", [](TreeLike& self, NodeId nodeId) { return collectNodeIds(topology(self).getNodeSubtree(nodeId)); }, "nodeId"_a,
            "Return all alive node ids in the subtree rooted at `nodeId`.")
        .def(
            "getDescendants", [](TreeLike& self, NodeId nodeId) { return collectNodeIds(topology(self).getDescendants(nodeId)); }, "nodeId"_a,
            "Return all strict descendant node ids of `nodeId`.")
        .def(
            "getNodeParent", [](TreeLike& self, NodeId nodeId) { return topology(self).getNodeParent(nodeId); }, "nodeId"_a,
            "Return the parent node id of `nodeId`, or `InvalidNode` when absent.")
        .def(
            "getProperPartOwner", [](TreeLike& self, NodeId properPartId) { return topology(self).getProperPartOwner(properPartId); }, "properPartId"_a,
            "Return the internal node that directly owns `properPartId`.")
        .def(
            "getHigraNodeId", [](TreeLike& self, NodeId nodeId) { return topology(self).getHigraNodeId(nodeId); }, "nodeId"_a,
            "Return the preserved imported Higra node id for a live internal NodeId, or InvalidNode.")
        .def(
            "getNumChildren", [](TreeLike& self, NodeId nodeId) { return topology(self).getNumChildren(nodeId); }, "nodeId"_a,
            "Return the number of direct children of `nodeId`.")
        .def(
            "getFirstChild", [](TreeLike& self, NodeId nodeId) { return topology(self).getFirstChild(nodeId); }, "nodeId"_a,
            "Return the first child of `nodeId`, or `InvalidNode`.")
        .def(
            "getNextSibling", [](TreeLike& self, NodeId nodeId) { return topology(self).getNextSibling(nodeId); }, "nodeId"_a,
            "Return the next sibling of `nodeId`, or `InvalidNode`.")
        .def(
            "isNode", [](TreeLike& self, NodeId nodeId) { return topology(self).isNode(nodeId); }, "nodeId"_a,
            "Return true when `nodeId` is in the internal-node slot domain.")
        .def(
            "isProperPart", [](TreeLike& self, NodeId nodeId) { return topology(self).isProperPart(nodeId); }, "nodeId"_a,
            "Return true when the id is in the proper-part domain.")
        .def(
            "isAlive", [](TreeLike& self, NodeId nodeId) { return topology(self).isAlive(nodeId); }, "nodeId"_a,
            "Return true when `nodeId` is an alive internal node.")
        .def(
            "isRoot", [](TreeLike& self, NodeId nodeId) { return topology(self).isRoot(nodeId); }, "nodeId"_a, "Return true when `nodeId` is the current root.")
        .def(
            "isLeaf", [](TreeLike& self, NodeId nodeId) { return topology(self).isLeaf(nodeId); }, "nodeId"_a,
            "Return true when `nodeId` is an alive leaf node.")
        .def(
            "isStructuralNode", [](TreeLike& self, NodeId nodeId) { return topology(self).isStructuralNode(nodeId); }, "nodeId"_a,
            "Return true when `nodeId` owns no direct proper parts.")
        .def(
            "hasChild", [](TreeLike& self, NodeId parentId, NodeId childId) { return topology(self).hasChild(parentId, childId); }, "parentId"_a, "childId"_a,
            "Return true when `childId` is a direct child of `parentId`.")
        .def(
            "pruneNode", [](TreeLike& self, NodeId nodeId) { self.pruneNode(nodeId); }, "nodeId"_a,
            "Prune one node from the topology, preserving a valid rooted tree.")
        .def(
            "mergeNodeIntoParent", [](TreeLike& self, NodeId nodeId) { self.mergeNodeIntoParent(nodeId); }, "nodeId"_a,
            "Merge `nodeId` into its parent and reattach descendants/proper parts.")
        .def_property_readonly(
            "descriptiveKind", [](TreeLike& self) { return topology(self).getDescriptiveKind(); }, "Optional descriptive hierarchy-family label.")
        .def_property_readonly(
            "altitudeOrder", [](TreeLike& self) { return topology(self).getAltitudeOrder(); }, "Global parent-to-child altitude ordering capability.")
        .def_property_readonly(
            "adjacencyMode", [](TreeLike& self) { return topology(self).getAdjacencyMode(); }, "Shape of the hierarchy's adjacency context.")
        .def_property_readonly(
            "hasUniformGridAdjacency2D", [](TreeLike& self) { return topology(self).hasUniformGridAdjacency2D(); },
            "Whether one immutable regular-grid 2D adjacency is available.")
        .def_property_readonly(
            "hasDirectionalGridAdjacency2D", [](TreeLike& self) { return topology(self).hasDirectionalGridAdjacency2D(); },
            "Whether decreasing/increasing regular-grid 2D relations are available.")
        .def_property_readonly(
            "hasGridDomain2D", [](TreeLike& self) { return topology(self).hasGridDomain2D(); }, "Whether proper-part ids have a regular row/column layout.")
        .def_property_readonly(
            "gridDomain2D", [](TreeLike& self) { return topology(self).getGridDomain2D(); }, "Optional regular 2D proper-part domain.")
        .def(
            "getUniformGridAdjacency2D",
            [](TreeLike& self) -> const RegularGridAdjacency2D& {
                const RegularGridAdjacency2D* adjacency = topology(self).getUniformGridAdjacency2D();
                if (!adjacency) {
                    throw std::runtime_error("Uniform regular-grid 2D adjacency context is not available.");
                }
                return *adjacency;
            },
            py::return_value_policy::reference_internal, "Return the borrowed immutable regular-grid 2D adjacency.")
        .def(
            "getDirectionalGridAdjacency2D",
            [](TreeLike& self) -> const DirectionalGridAdjacency2D& {
                const DirectionalGridAdjacency2D* adjacency = topology(self).getDirectionalGridAdjacency2D();
                if (!adjacency) {
                    throw std::runtime_error("Directional regular-grid 2D adjacency context is not available.");
                }
                return *adjacency;
            },
            py::return_value_policy::reference_internal, "Return the borrowed directional regular-grid 2D adjacency context.")
        .def(
            "getDecreasingGridAdjacency2D",
            [](TreeLike& self) -> const RegularGridAdjacency2D& {
                const RegularGridAdjacency2D* adjacency = topology(self).getDecreasingGridAdjacency2D();
                if (!adjacency) {
                    throw std::runtime_error("Decreasing regular-grid 2D adjacency context is not available.");
                }
                return *adjacency;
            },
            py::return_value_policy::reference_internal, "Return the borrowed adjacency used by decreasing branches.")
        .def(
            "getIncreasingGridAdjacency2D",
            [](TreeLike& self) -> const RegularGridAdjacency2D& {
                const RegularGridAdjacency2D* adjacency = topology(self).getIncreasingGridAdjacency2D();
                if (!adjacency) {
                    throw std::runtime_error("Increasing regular-grid 2D adjacency context is not available.");
                }
                return *adjacency;
            },
            py::return_value_policy::reference_internal, "Return the borrowed adjacency used by increasing branches.")
        .def_property_readonly(
            "numRows", [](TreeLike& self) { return topology(self).getNumRowsOfGridDomain2D(); }, "Number of rows in the active proper-part domain.")
        .def_property_readonly(
            "numCols", [](TreeLike& self) { return topology(self).getNumColsOfGridDomain2D(); }, "Number of columns in the active proper-part domain.")
        .def_property_readonly("numNodes", [](TreeLike& self) { return topology(self).getNumNodes(); }, "Number of currently alive internal nodes.");

    if constexpr (requires(TreeLike& self, NodeId nodeId) {
                      typename TreeLike::altitude_type;
                      self.getAltitude(nodeId);
                      self.getNodeResidue(nodeId);
                      self.reconstructionImage();
                      self.exportHigraHierarchy();
                  }) {
        cls.def(
               "getAltitude", [](TreeLike& self, NodeId nodeId) { return getAltitudeValue(self, nodeId); }, "nodeId"_a,
               "Return the altitude associated with nodeId.")
            .def(
                "getNodeResidue", [](TreeLike& self, NodeId nodeId) { return getResidueValue(self, nodeId); }, "nodeId"_a,
                "Return the residue between nodeId and its parent.")
            .def(
                "reconstructionImage", [](TreeLike& self) { return reconstructionImageOf(self); },
                "Reconstruct the current tree into a 2D image using the attached GridDomain2D.")
            .def(
                "exportHigraHierarchy", [](TreeLike& self) { return exportHigraHierarchyOf(self); },
                "Export the current rooted tree to a new compact Higra (parent, altitude) representation.");
    }
}

} // namespace

/**
 * @brief Registers morphological tree bindings in the Python module.
 *
 * @param m Python module receiving the bindings.
 */
void initMorphologicalTree(py::module_& m) {
    py::enum_<ToSInterpolation>(m, "ToSInterpolation", py::module_local(false))
        .value("SelfDual", ToSInterpolation::SelfDual)
        .value("Min4cMax8c", ToSInterpolation::Min4cMax8c)
        .value("Min8cMax4c", ToSInterpolation::Min8cMax4c)
        .export_values();

    py::enum_<ToSPaddingPolicy>(m, "ToSPaddingPolicy", py::module_local(false))
        .value("Exterior", ToSPaddingPolicy::Exterior)
        .value("NoPadding", ToSPaddingPolicy::NoPadding);

    py::class_<TreeOfShapesProducerOptions>(m, "TreeOfShapesProducerOptions", py::module_local(false),
                                            "Tree-of-Shapes producer options. These settings affect construction "
                                            "only; the resulting tree remains on the original image domain.")
        .def(py::init([](ToSInterpolation interpolation, ToSPaddingPolicy padding, int infinitySeedRow, int infinitySeedCol) {
                 return TreeOfShapesProducerOptions{interpolation, padding, infinitySeedRow, infinitySeedCol};
             }),
             "interpolation"_a = ToSInterpolation::SelfDual, "padding"_a = ToSPaddingPolicy::Exterior, "infinitySeedRow"_a = ToSDefaultInfinityRow,
             "infinitySeedCol"_a = ToSDefaultInfinityCol)
        .def_readwrite("interpolation", &TreeOfShapesProducerOptions::interpolation)
        .def_readwrite("padding", &TreeOfShapesProducerOptions::padding)
        .def_readwrite("infinitySeedRow", &TreeOfShapesProducerOptions::infinitySeedRow)
        .def_readwrite("infinitySeedCol", &TreeOfShapesProducerOptions::infinitySeedCol);

    py::enum_<NodeIdSpace>(m, "NodeIdSpace", py::module_local(false))
        .value("MORPHOLOGICAL_TREE", NodeIdSpace::MORPHOLOGICAL_TREE)
        .value("HIGRA", NodeIdSpace::HIGRA)
        .export_values();

    py::enum_<MorphologicalTreeKind>(m, "MorphologicalTreeKind", py::module_local(false))
        .value("GENERIC", MorphologicalTreeKind::GENERIC)
        .value("MAX_TREE", MorphologicalTreeKind::MAX_TREE)
        .value("MIN_TREE", MorphologicalTreeKind::MIN_TREE)
        .value("TREE_OF_SHAPES", MorphologicalTreeKind::TREE_OF_SHAPES)
        .value("SELF_DUAL_RESIDUAL_TREE", MorphologicalTreeKind::SELF_DUAL_RESIDUAL_TREE)
        .export_values();

    py::enum_<sdrt::SdrtTiePolicy>(m, "SdrtTiePolicy", py::module_local(false))
        .value("MAX_BEFORE_MIN_THEN_SPATIAL", sdrt::SdrtTiePolicy::MaxBeforeMinThenSpatial)
        .value("CONTRAST_INVARIANT_SPATIAL", sdrt::SdrtTiePolicy::ContrastInvariantSpatial)
        .export_values();

    py::enum_<AltitudeOrder>(m, "AltitudeOrder", py::module_local(false))
        .value("INCREASING_FROM_ROOT", AltitudeOrder::INCREASING_FROM_ROOT, "Every parent altitude is strictly smaller than each child altitude.")
        .value("DECREASING_FROM_ROOT", AltitudeOrder::DECREASING_FROM_ROOT, "Every parent altitude is strictly greater than each child altitude.")
        .value("UNCONSTRAINED", AltitudeOrder::UNCONSTRAINED, "No global parent-child altitude order is declared.")
        .export_values();

    py::enum_<AdjacencyMode>(m, "AdjacencyMode", py::module_local(false))
        .value("NONE", AdjacencyMode::NONE)
        .value("UNIFORM", AdjacencyMode::UNIFORM)
        .value("DIRECTIONAL", AdjacencyMode::DIRECTIONAL)
        .export_values();

    py::class_<GridDomain2D>(m, "GridDomain2D", py::module_local(false), "Optional row/column layout attached to a proper-part domain.")
        .def(py::init<int, int>(), "rows"_a, "cols"_a)
        .def_readonly("rows", &GridDomain2D::rows)
        .def_readonly("cols", &GridDomain2D::cols);

    py::class_<DirectionalGridAdjacency2D>(m, "DirectionalGridAdjacency2D", py::module_local(false),
                                           "Generic adjacency context for decreasing and increasing hierarchy branches.")
        .def(py::init([](RegularGridAdjacency2D decreasing, RegularGridAdjacency2D increasing) {
                 return DirectionalGridAdjacency2D{std::move(decreasing), std::move(increasing)};
             }),
             "decreasing"_a, "increasing"_a)
        .def_readonly("decreasing", &DirectionalGridAdjacency2D::decreasing, "Adjacency used by decreasing branches.")
        .def_readonly("increasing", &DirectionalGridAdjacency2D::increasing, "Adjacency used by increasing branches.");
    py::class_<HierarchySemantics>(m, "HierarchySemantics", py::module_local(false),
                                   "Generic interpretation capabilities attached to a morphological hierarchy.")
        .def(py::init([](AltitudeOrder altitudeOrder, MorphologicalTreeKind descriptiveKind, std::optional<RegularGridAdjacency2D> uniformAdjacency,
                         std::optional<DirectionalGridAdjacency2D> directionalAdjacency) {
                 if (uniformAdjacency && directionalAdjacency) {
                     throw std::invalid_argument("HierarchySemantics accepts either uniform or directional adjacency, not both.");
                 }
                 AdjacencyContext adjacency = NoAdjacency{};
                 if (directionalAdjacency) {
                     adjacency = std::move(*directionalAdjacency);
                 } else if (uniformAdjacency) {
                     adjacency = UniformGridAdjacency2D{std::move(*uniformAdjacency)};
                 }
                 return HierarchySemantics{descriptiveKind, altitudeOrder, std::move(adjacency)};
             }),
             "altitudeOrder"_a = AltitudeOrder::UNCONSTRAINED, "descriptiveKind"_a = MorphologicalTreeKind::GENERIC, "uniformAdjacency"_a = py::none(),
             "directionalAdjacency"_a = py::none())
        .def_readonly("descriptiveKind", &HierarchySemantics::descriptiveKind)
        .def_readonly("altitudeOrder", &HierarchySemantics::altitudeOrder)
        .def_property_readonly("adjacencyMode", &HierarchySemantics::adjacencyMode);

    py::enum_<EdgeToPixelReducer>(m, "EdgeToPixelReducer", py::module_local(false))
        .value("Max", EdgeToPixelReducer::Max)
        .value("Mean", EdgeToPixelReducer::Mean)
        .export_values();

    py::enum_<HierarchyLevelConvention>(m, "HierarchyLevelConvention", py::module_local(false))
        .value("EdgeSaliencyValue", HierarchyLevelConvention::EdgeSaliencyValue)
        .value("PartitionAppearanceLevel", HierarchyLevelConvention::PartitionAppearanceLevel)
        .export_values();

    py::enum_<ShapeSpaceExtremaPolarity>(m, "ShapeSpaceExtremaPolarity", py::module_local(false))
        .value("Minima", ShapeSpaceExtremaPolarity::Minima)
        .value("Maxima", ShapeSpaceExtremaPolarity::Maxima)
        .export_values();

    auto weightedCls = py::class_<WeightedMorphologicalTree<std::uint8_t>, std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>>>(
        m, "WeightedMorphologicalTree", py::module_local(false),
        "Wrapper pairing MorphologicalTree topology with an external dense altitude buffer. "
        "Imports can preserve an original Higra node-id domain; exports always create a new compact Higra domain.");

    py::class_<MorphologicalTreeFactory>(m, "MorphologicalTreeFactory", py::module_local(false))
        .def_static(
            "createMaxTree",
            [](UInt8InputArray input, double radius) {
                const double validatedRadius = pybind_utils::requireAdjacencyRadius(radius, "MorphologicalTreeFactory.createMaxTree");
                return std::make_shared<WeightedMorphologicalTree<std::uint8_t>>(
                    MorphologicalTreeFactory::createMaxTree(imageFromArray(input), validatedRadius));
            },
            "input"_a, "radius"_a = 1.5, "Create a max-tree from a 2D C-contiguous `np.uint8` image.")
        .def_static(
            "createMaxTree",
            [](UInt8InputArray input, RegularGridAdjacency2D adjacency) {
                return std::make_shared<WeightedMorphologicalTree<std::uint8_t>>(
                    MorphologicalTreeFactory::createMaxTree(imageFromArray(input), std::move(adjacency)));
            },
            "input"_a, "adjacency"_a, "Create a max-tree with an explicit regular-grid 2D adjacency.")
        .def_static(
            "createMinTree",
            [](UInt8InputArray input, double radius) {
                const double validatedRadius = pybind_utils::requireAdjacencyRadius(radius, "MorphologicalTreeFactory.createMinTree");
                return std::make_shared<WeightedMorphologicalTree<std::uint8_t>>(
                    MorphologicalTreeFactory::createMinTree(imageFromArray(input), validatedRadius));
            },
            "input"_a, "radius"_a = 1.5, "Create a min-tree from a 2D C-contiguous `np.uint8` image.")
        .def_static(
            "createMinTree",
            [](UInt8InputArray input, RegularGridAdjacency2D adjacency) {
                return std::make_shared<WeightedMorphologicalTree<std::uint8_t>>(
                    MorphologicalTreeFactory::createMinTree(imageFromArray(input), std::move(adjacency)));
            },
            "input"_a, "adjacency"_a, "Create a min-tree with an explicit regular-grid 2D adjacency.")
        .def_static(
            "createSelfDualResidualTree",
            [](UInt8InputArray input, double radius, sdrt::SdrtTiePolicy tiePolicy) {
                const double validatedRadius = pybind_utils::requireAdjacencyRadius(radius, "MorphologicalTreeFactory.createSelfDualResidualTree");
                return std::make_shared<WeightedMorphologicalTree<std::uint8_t>>(
                    MorphologicalTreeFactory::createSelfDualResidualTree(
                        imageFromArray(input), validatedRadius, sdrt::UnrestrictedResidualTreeOptions{tiePolicy}));
            },
            "input"_a, "radius"_a = 1.5, "tiePolicy"_a = sdrt::SdrtTiePolicy::ContrastInvariantSpatial,
            "Create the unrestricted residual tree with synchronized max-tree and min-tree states.")
        .def_static(
            "createSelfDualResidualTree",
            [](UInt8InputArray input, RegularGridAdjacency2D adjacency, sdrt::SdrtTiePolicy tiePolicy) {
                return std::make_shared<WeightedMorphologicalTree<std::uint8_t>>(
                    MorphologicalTreeFactory::createSelfDualResidualTree(
                        imageFromArray(input), std::move(adjacency), sdrt::UnrestrictedResidualTreeOptions{tiePolicy}));
            },
            "input"_a, "adjacency"_a, "tiePolicy"_a = sdrt::SdrtTiePolicy::ContrastInvariantSpatial,
            "Create the unrestricted residual tree with an explicit shared symmetric adjacency.")
        .def_static(
            "createSaturatedSelfDualResidualTree",
            [](UInt8InputArray input, NodeId infinityPixel, double radius, sdrt::SdrtTiePolicy tiePolicy) {
                const double validatedRadius = pybind_utils::requireAdjacencyRadius(radius, "MorphologicalTreeFactory.createSaturatedSelfDualResidualTree");
                return std::make_shared<WeightedMorphologicalTree<std::uint8_t>>(
                    MorphologicalTreeFactory::createSaturatedSelfDualResidualTree(
                        imageFromArray(input), infinityPixel, validatedRadius, sdrt::SaturatedResidualTreeOptions{tiePolicy}));
            },
            "input"_a, "infinityPixel"_a = NodeId{0}, "radius"_a = 1.5, "tiePolicy"_a = sdrt::SdrtTiePolicy::ContrastInvariantSpatial,
            "Create the saturated residual tree with synchronized max-tree and min-tree states.")
        .def_static(
            "createSaturatedSelfDualResidualTree",
            [](UInt8InputArray input, RegularGridAdjacency2D adjacency, NodeId infinityPixel, sdrt::SdrtTiePolicy tiePolicy) {
                return std::make_shared<WeightedMorphologicalTree<std::uint8_t>>(
                    MorphologicalTreeFactory::createSaturatedSelfDualResidualTree(
                        imageFromArray(input), std::move(adjacency), infinityPixel, sdrt::SaturatedResidualTreeOptions{tiePolicy}));
            },
            "input"_a, "adjacency"_a, "infinityPixel"_a = NodeId{0}, "tiePolicy"_a = sdrt::SdrtTiePolicy::ContrastInvariantSpatial,
            "Create the saturated residual tree with an explicit shared symmetric adjacency.")
        .def_static(
            "createTreeOfShapes",
            [](UInt8InputArray input, ToSInterpolation interpolation, int infinitySeedRow, int infinitySeedCol) {
                return std::make_shared<WeightedMorphologicalTree<std::uint8_t>>(
                    MorphologicalTreeFactory::createTreeOfShapes(imageFromArray(input), interpolation, infinitySeedRow, infinitySeedCol));
            },
            "input"_a, "interpolation"_a = ToSInterpolation::SelfDual, "infinitySeedRow"_a = ToSDefaultInfinityRow, "infinitySeedCol"_a = ToSDefaultInfinityCol,
            "Create a tree of shapes on the original image domain from a 2D C-contiguous `np.uint8` image.")
        .def_static(
            "createTreeOfShapes",
            [](UInt8InputArray input, TreeOfShapesProducerOptions options) {
                return std::make_shared<WeightedMorphologicalTree<std::uint8_t>>(MorphologicalTreeFactory::createTreeOfShapes(imageFromArray(input), options));
            },
            "input"_a, "options"_a,
            "Create a tree of shapes with producer-local interpolation, "
            "padding, and infinity-seed options. The returned tree uses the "
            "original image domain.")
        .def_static(
            "createFromNativeTopology",
            [](const std::vector<NodeId>& nodeParent, const std::vector<NodeId>& properPartOwner, py::object altitudeInput, NodeId root,
               HierarchySemantics semantics) {
                const std::vector<std::uint8_t> altitude =
                    pythonUInt8AltitudeVector(altitudeInput, "MorphologicalTreeFactory.createFromNativeTopology altitude");
                return std::make_shared<WeightedMorphologicalTree<std::uint8_t>>(
                    MorphologicalTreeFactory::createFromNativeTopology(std::span<const NodeId>(nodeParent), std::span<const NodeId>(properPartOwner),
                                                                       std::span<const std::uint8_t>(altitude), root, std::move(semantics)));
            },
            "nodeParent"_a, "properPartOwner"_a, "altitude"_a, "root"_a, "semantics"_a,
            R"doc(Create a morphological tree over an abstract finite proper-part domain.

No row/column interpretation is attached. Topological and support-based
algorithms remain available; image reconstruction and geometric attributes
require a separate regular 2D domain.)doc")
        .def_static(
            "createFromNativeTopology",
            [](const std::vector<NodeId>& nodeParent, const std::vector<NodeId>& properPartOwner, py::object altitudeInput, NodeId root, int rows, int cols,
               HierarchySemantics semantics) {
                const std::vector<std::uint8_t> altitude =
                    pythonUInt8AltitudeVector(altitudeInput, "MorphologicalTreeFactory.createFromNativeTopology altitude");
                return std::make_shared<WeightedMorphologicalTree<std::uint8_t>>(
                    MorphologicalTreeFactory::createFromNativeTopology(std::span<const NodeId>(nodeParent), std::span<const NodeId>(properPartOwner),
                                                                       std::span<const std::uint8_t>(altitude), root, rows, cols, std::move(semantics)));
            },
            "nodeParent"_a, "properPartOwner"_a, "altitude"_a, "root"_a, "rows"_a, "cols"_a, "semantics"_a,
            R"doc(Create a morphological tree from native partial-partition buffers and explicit generic capabilities.

`nodeParent` and `altitude` use the dense internal-node domain.
`properPartOwner` uses the row-major proper-part domain and must contain
`rows * cols` entries. Every committed node must have non-empty subtree
support, although structural nodes may own no direct proper parts.)doc")
        .def_static(
            "createFromHigraParent",
            [](const std::vector<NodeId>& parent, py::object altitudeInput, int rows, int cols, MorphologicalTreeKind kind, std::optional<double> radius) {
                const std::vector<std::uint8_t> altitude = pythonUInt8AltitudeVector(altitudeInput, "MorphologicalTreeFactory.createFromHigraParent altitude");
                if (parent.size() != altitude.size()) {
                    throw std::invalid_argument("parent and altitude must have the same size");
                }

                std::optional<RegularGridAdjacency2D> adjacency;
                if (radius.has_value()) {
                    adjacency.emplace(pybind_utils::makeRegularGridAdjacency2D(rows, cols, *radius, "MorphologicalTreeFactory.createFromHigraParent"));
                }

                return std::make_shared<WeightedMorphologicalTree<std::uint8_t>>(MorphologicalTreeFactory::createFromHigraParent(
                    std::span<const NodeId>(parent), std::span<const std::uint8_t>(altitude), rows, cols, kind, std::move(adjacency)));
            },
            "parent"_a, "altitude"_a, "rows"_a, "cols"_a, "kind"_a, "radius"_a = py::none(),
            "Create a weighted tree from an imported static Higra parent/altitude representation [leaves | internal nodes]. "
            "The imported Higra node-id domain is preserved until the tree is edited.");

    py::class_<HierarchySaliencyMapValidation>(m, "HierarchySaliencyMapValidation", py::module_local(false),
                                               "Utilities for validating and transforming hierarchy valuations used by saliency maps.")
        .def_static(
            "validateHierarchyConnectivity",
            [](std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weighted, std::optional<double> radius) {
                if (!weighted) {
                    throw std::invalid_argument("weighted tree must not be null");
                }
                const RegularGridAdjacency2D adjacency = saliencyAdjacencyPy(*weighted, radius, "HierarchySaliencyMapValidation.validateHierarchyConnectivity");
                HierarchySaliencyMapValidation::validateHierarchyConnectivity(weighted->topology(), adjacency,
                                                                              "HierarchySaliencyMapValidation.validateHierarchyConnectivity");
            },
            "tree"_a, "radius"_a = py::none(),
            R"doc(Validate that every hierarchy support is connected in the projection graph.

This is the graph-connectivity hypothesis required by the Cousty correspondence;
it is distinct from checking that the parent array forms one rooted tree.)doc")
        .def_static(
            "validateHierarchyValuation",
            [](std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weighted, py::array valuation, bool strict, bool nonnegative) {
                if (!weighted) {
                    throw std::invalid_argument("weighted tree must not be null");
                }
                validateHierarchyValuationPy(*weighted, std::move(valuation), strict, nonnegative);
            },
            "tree"_a, "valuation"_a, "strict"_a = false, "nonnegative"_a = false,
            R"doc(Validate that a dense node-indexed valuation is compatible with the hierarchy.

The valuation must be a 1D int32, float32, or float64 array with one value per
dense internal NodeId slot. With `strict=False`, parent values must be greater
than or equal to child values, allowing equal-valued tree levels to collapse in
the induced QFZ hierarchy. With `strict=True`, every live parent-child relation
must be strictly increasing toward the root. Set `nonnegative=True` to enforce
the paper's literal non-negative edge-weight domain.)doc")
        .def_static(
            "rankHierarchyValuation",
            [](std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weighted, py::array valuation, bool strict) {
                if (!weighted) {
                    throw std::invalid_argument("weighted tree must not be null");
                }
                return rankHierarchyValuationPy(*weighted, std::move(valuation), strict);
            },
            "tree"_a, "valuation"_a, "strict"_a = false,
            R"doc(Convert a compatible hierarchy valuation to dense non-negative integer levels.

The input valuation must be a 1D int32, float32, or float64 array. Distinct
live-node valuation values are ranked as `0..k-1` while preserving order. Equal
values receive the same rank, so level collapse is preserved. Use `strict=True`
to reject equal parent-child valuation levels before ranking. This ranks all
live nodes; use `HierarchySaliencyMap.computeCanonicalRankedSaliencyEdgeMap`
when only levels that occur on graph edges should define the dense scale.)doc")
        .def_static(
            "computeNormalizedScores",
            [](std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weighted, py::array valuation, bool strict, bool nonnegative) {
                if (!weighted) {
                    throw std::invalid_argument("weighted tree must not be null");
                }
                return computeNormalizedScoresPy(*weighted, std::move(valuation), strict, nonnegative);
            },
            "tree"_a, "valuation"_a, "strict"_a = false, "nonnegative"_a = false,
            R"doc(Normalize a compatible hierarchy valuation to double values in [0, 1].

The input valuation must be a 1D int32, float32, or float64 array. The valuation
is validated before normalization. The output preserves the valuation order with
an increasing affine transform over live-node values; equal levels remain equal.
Set `strict=True` to reject equal parent-child levels before normalization and
`nonnegative=True` to require the input valuation itself to be non-negative.)doc");

    py::class_<HierarchySaliencyMap>(m, "HierarchySaliencyMap", py::module_local(false),
                                     R"doc(Utilities for projecting a morphological hierarchy onto an image adjacency graph as an edge saliency map.

Primary reference: Jean Cousty, Laurent Najman, Yukiko Kenmochi, and Silvio
Guimarães, "Hierarchical segmentations with graphs: quasi-flat zones, minimum
spanning trees, and saliency maps," Journal of Mathematical Imaging and Vision
60(4):479-502, 2018, https://doi.org/10.1007/s10851-017-0768-7. The direct
projection corresponds to Section 4, Equations (5)-(6), and Section 7,
Algorithm 1.)doc")
        .def_static(
            "computeSaliencyEdgeMap",
            [](std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weighted, py::array valuation, std::optional<double> radius, bool strict,
               HierarchyLevelConvention levelConvention, bool validateConnectivity) {
                if (!weighted) {
                    throw std::invalid_argument("weighted tree must not be null");
                }
                return computeSaliencyEdgeMapPy(*weighted, std::move(valuation), radius, strict, levelConvention, validateConnectivity,
                                                "HierarchySaliencyMap.computeSaliencyEdgeMap");
            },
            "tree"_a, "valuation"_a, "radius"_a = py::none(), "strict"_a = false, "levelConvention"_a = HierarchyLevelConvention::EdgeSaliencyValue,
            "validateConnectivity"_a = true,
            R"doc(Compute the formal edge-indexed saliency map induced by a hierarchy valuation.

This method implements `Phi(H)` for a connected hierarchy under the
quasi-flat-zone saliency convention: each adjacency edge receives
`valuation[LCA(owner(source), owner(target))]` when the endpoint owners differ,
and value 0 when both endpoints already belong to the same finest represented
region under `EdgeSaliencyValue`. With `PartitionAppearanceLevel`, transition
edges receive `valuation[LCA] - 1`, exactly as in Algorithm 1 of Cousty et al.
It does not build the full `Psi(w) = Phi(QFZ(G, w))` pipeline from an
arbitrary edge-weighted graph. The valuation is validated before projection. Use
`strict=True` when the map must preserve every explicit parent-child level of
the tree; leave it false when equal-valued adjacent levels may be collapsed.
Formal saliency valuations must be non-negative. Graph connectivity is validated
by default; disable it only for a trusted producer-internal path.)doc")
        .def_static(
            "computeCanonicalRankedSaliencyEdgeMap",
            [](std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weighted, py::array valuation, std::optional<double> radius, bool strict,
               bool validateConnectivity) {
                if (!weighted) {
                    throw std::invalid_argument("weighted tree must not be null");
                }
                return computeCanonicalRankedSaliencyEdgeMapPy(*weighted, std::move(valuation), radius, strict, validateConnectivity);
            },
            "tree"_a, "valuation"_a, "radius"_a = py::none(), "strict"_a = false, "validateConnectivity"_a = true,
            R"doc(Project and rank only hierarchy levels that actually occur on graph edges.

Same-owner edges form the base rank when present. Unlike node-wise ranking,
unused leaf or incomparable-node values cannot introduce gaps in the canonical
edge scale. The input may contain negative values because only its order is
retained.)doc")
        .def_static(
            "computeTopologicalLevelEdgeMap",
            [](std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weighted, std::optional<double> radius) {
                if (!weighted) {
                    throw std::invalid_argument("weighted tree must not be null");
                }
                return computeTopologicalLevelEdgeMapPy(*weighted, radius);
            },
            "tree"_a, "radius"_a = py::none(),
            R"doc(Compute a topological-level hierarchy saliency map on adjacency edges.

Leaf internal nodes receive level 0, and each ancestor receives one plus the
maximum level of its children.)doc")
        .def_static(
            "computeNormalizedAltitudeEdgeMap",
            [](std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weighted, std::optional<double> radius) {
                if (!weighted) {
                    throw std::invalid_argument("weighted tree must not be null");
                }
                return computeNormalizedAltitudeEdgeMapPy(*weighted, radius);
            },
            "tree"_a, "radius"_a = py::none(),
            R"doc(Compute a normalized-altitude hierarchy saliency map on adjacency edges.

Values are doubles in [0, 1], oriented so that ancestors are greater than or
equal to descendants for max-trees and min-trees.)doc");

    py::class_<ComponentTreePartitionHierarchyAdapter>(m, "ComponentTreePartitionHierarchyAdapter", py::module_local(false),
                                                       "Explicit proper-part completion of a component tree into connected graph partitions.")
        .def_static(
            "validate",
            [](std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weighted, std::optional<double> radius) {
                if (!weighted) {
                    throw std::invalid_argument("weighted tree must not be null");
                }
                const RegularGridAdjacency2D adjacency = saliencyAdjacencyPy(*weighted, radius, "ComponentTreePartitionHierarchyAdapter.validate");
                ComponentTreePartitionHierarchyAdapter::validate(weighted->topology(), adjacency);
            },
            "tree"_a, "radius"_a = py::none(), "Validate the connected proper-part completion in the selected graph.")
        .def_static(
            "computePartitionAppearanceLevels",
            [](std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weighted) {
                if (!weighted) {
                    throw std::invalid_argument("weighted tree must not be null");
                }
                std::vector<int> levels = ComponentTreePartitionHierarchyAdapter::computePartitionAppearanceLevels(weighted->topology());
                const int size = static_cast<int>(levels.size());
                return pybind_utils::toNumpyOwned(std::move(levels), size);
            },
            "tree"_a,
            R"doc(Return positive partition-appearance levels for the component-tree completion.

Use these levels with `HierarchyLevelConvention.PartitionAppearanceLevel` to
apply the exact `level(LCA)-1` convention.)doc")
        .def_static(
            "computeSaliencyEdgeMap",
            [](std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weighted, py::array partitionAppearanceLevels, std::optional<double> radius,
               bool strict) {
                if (!weighted) {
                    throw std::invalid_argument("weighted tree must not be null");
                }
                return computeSaliencyEdgeMapPy(*weighted, std::move(partitionAppearanceLevels), radius, strict,
                                                HierarchyLevelConvention::PartitionAppearanceLevel, true,
                                                "ComponentTreePartitionHierarchyAdapter.computeSaliencyEdgeMap");
            },
            "tree"_a, "partitionAppearanceLevels"_a, "radius"_a = py::none(), "strict"_a = false,
            R"doc(Project positive partition-appearance levels using `level(LCA)-1`.

The component-tree support completion is validated for connectivity in the
selected graph before projection.)doc");

    py::class_<ShapeSpaceSaliency>(m, "ShapeSpaceSaliency", py::module_local(false),
                                   R"doc(Xu shaping extinction values on the tree-node graph and their maximum-on-contours image projection.

Primary reference: Yongchao Xu, Edwin Carlinet, Thierry Géraud, and Laurent
Najman, "Hierarchical Segmentation Using Tree-Based Shape Spaces," IEEE
Transactions on Pattern Analysis and Machine Intelligence 39(3):457-469, 2017,
https://doi.org/10.1109/TPAMI.2016.2554550, Section 4.3. The minima construction
is generalized here to minima or maxima and regular-grid graph edges.)doc")
        .def_static(
            "computeExtinctionValues",
            [](std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weighted, py::array attribute, ShapeSpaceExtremaPolarity polarity) {
                if (!weighted) {
                    throw std::invalid_argument("weighted tree must not be null");
                }
                return computeShapeSpaceExtinctionValuesPy(*weighted, std::move(attribute), polarity);
            },
            "tree"_a, "attribute"_a, "polarity"_a,
            R"doc(Compute shaping extinction values on the original tree-node graph.

`attribute` must be a dense float32 or float64 array with one value per
internal NodeId slot. The values need not be monotone on the original tree.
Parent-child relations form the graph for a second component-tree computation;
`polarity` selects its local minima or maxima. The returned dictionary contains
an `extrema` list and the dense extinction array `nodeScores`, preserving the
input dtype.)doc")
        .def_static(
            "projectContourScores",
            [](std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weighted, py::array nodeScores, std::optional<double> radius) {
                if (!weighted) {
                    throw std::invalid_argument("weighted tree must not be null");
                }
                return projectShapeSpaceContourScoresPy(*weighted, std::move(nodeScores), radius);
            },
            "tree"_a, "nodeScores"_a, "radius"_a = py::none(),
            R"doc(Project dense original-node scores onto their full region contours.

Every image adjacency edge receives the maximum score among original regions
whose contour contains that edge. The return value is the standard edge-map
dictionary. If `radius` is omitted, the tree's stored construction adjacency is
used. This operation is distinct from projecting an LCA valuation.)doc")
        .def_static(
            "compute",
            [](std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weighted, py::array attribute, ShapeSpaceExtremaPolarity polarity,
               std::optional<double> radius) {
                if (!weighted) {
                    throw std::invalid_argument("weighted tree must not be null");
                }
                return computeShapeSpaceSaliencyPy(*weighted, std::move(attribute), polarity, radius);
            },
            "tree"_a, "attribute"_a, "polarity"_a, "radius"_a = py::none(),
            R"doc(Compute Xu shaping extinctions and their maximum-on-contours map.

The returned dictionary contains `extrema`, dense extinction `nodeScores`, and
the standard edge-map dictionary `edgeMap`. If `radius` is omitted, the tree's
stored construction adjacency is used.)doc");

    py::class_<HierarchySaliencyMapProjection>(m, "HierarchySaliencyMapProjection", py::module_local(false),
                                               "Derived projections and contour materializations of hierarchy edge saliency maps.")
        .def_static(
            "edgeMapToPixelImage", [](py::dict edgeMap, EdgeToPixelReducer reducer) { return edgeMapToPixelImagePy(std::move(edgeMap), reducer); }, "edgeMap"_a,
            "reducer"_a = EdgeToPixelReducer::Max,
            R"doc(Rasterize an edge-indexed map into a pixel image for display.

This helper is a visualization projection only. Each edge value contributes to
both endpoint pixels. `EdgeToPixelReducer.Max` writes the maximum incident edge
value per pixel; `EdgeToPixelReducer.Mean` writes the arithmetic mean of
incident edge values. The formal saliency representation remains the edge map
dictionary.)doc")
        .def_static(
            "thresholdCut", [](py::dict edgeMap, double threshold) { return thresholdCutPy(std::move(edgeMap), threshold); }, "edgeMap"_a, "threshold"_a,
            R"doc(Threshold an edge saliency map into a contour edge set.

The input is the dictionary returned by `HierarchySaliencyMap`. The returned
dictionary contains `numRows`, `numCols`, `adjacencyRadius`, `sources`, and
`targets`, selecting edges whose saliency is greater than or equal to
`threshold`.)doc")
        .def_static(
            "nodeContourEdges",
            [](std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weighted, std::optional<double> radius) {
                if (!weighted) {
                    throw std::invalid_argument("weighted tree must not be null");
                }
                return nodeContourEdgesPy(*weighted, radius);
            },
            "tree"_a, "radius"_a = py::none(),
            R"doc(Project transition adjacency edges onto their hierarchy owner node.

Each returned edge carries the node id `LCA(owner(source), owner(target))` in the
`nodes` array. Edges whose endpoints have the same owner are omitted and have
implicit value 0 in the full formal saliency map. If `radius` is not provided,
the tree's stored construction adjacency is used.)doc")
        .def_static(
            "computeIncrementalNodeContours",
            [](std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weighted, std::optional<double> radius) {
                if (!weighted) {
                    throw std::invalid_argument("weighted tree must not be null");
                }
                return computeIncrementalNodeContoursPy(*weighted, radius);
            },
            "tree"_a, "radius"_a = py::none(),
            R"doc(Compute per-node incremental contour edges.

The returned dictionary contains `offsets`, `sources`, and `targets`.
For node `u`, the slice `offsets[u]:offsets[u + 1]` stores all transition
adjacency edges whose endpoint owners differ and whose endpoint-owner LCA is
`u`.)doc")
        .def_static(
            "projectNodeValuation",
            [](py::dict contours, py::array nodeValuation) {
                return projectNodeValuationOnIncrementalContoursPy(std::move(contours), std::move(nodeValuation));
            },
            "contours"_a, "nodeValuation"_a,
            R"doc(Project a dense node valuation onto per-node transition contour edges.

`contours` must be the dictionary returned by `computeIncrementalNodeContours`.
`nodeValuation` must be a 1D float32 or float64 array with one value per dense
internal node slot. The result is sparse over transition edges; same-owner graph
edges are omitted and have implicit value 0 in the full formal saliency map.)doc")
        .def_static(
            "thresholdByNodeValuation",
            [](py::dict contours, py::array nodeValuation, double threshold) {
                return thresholdIncrementalContoursByNodeValuationPy(std::move(contours), std::move(nodeValuation), threshold);
            },
            "contours"_a, "nodeValuation"_a, "threshold"_a,
            R"doc(Threshold per-node transition contours by dense node valuation.

All edges in a node slice are selected when `nodeValuation[node] >= threshold`.)doc");

    weightedCls
        .def(
            "setAltitude",
            [](WeightedMorphologicalTree<std::uint8_t>& tree, NodeId nodeId, py::object altitudeInput) {
                if (!tree.topology().isNode(nodeId)) {
                    throw std::invalid_argument("invalid NodeId for altitude update");
                }
                const std::uint8_t altitude = pythonUInt8AltitudeValue(altitudeInput, nodeId, "WeightedMorphologicalTree.setAltitude");
                tree.setAltitude(nodeId, altitude);
            },
            "nodeId"_a, "altitude"_a, "Set one node altitude while preserving the declared strict parent-child altitude order.")
        .def(
            "setAltitudeBuffer",
            [](WeightedMorphologicalTree<std::uint8_t>& tree, py::object altitudeInput) {
                const std::vector<std::uint8_t> altitude = pythonUInt8AltitudeVector(altitudeInput, "WeightedMorphologicalTree.setAltitudeBuffer");
                pybind_utils::requireVectorSize(altitude, static_cast<std::size_t>(tree.topology().getNumInternalNodeSlots()), "altitude");
                tree.setAltitudeBuffer(altitude);
            },
            "altitude"_a, "Replace the dense altitude buffer indexed by internal NodeId after validating the declared strict altitude order.")
        .def_property(
            "altitude", [](const WeightedMorphologicalTree<std::uint8_t>& tree) { return AltitudeBuffer<std::uint8_t>(tree.getAltitudeBuffer()); },
            [](WeightedMorphologicalTree<std::uint8_t>& tree, py::object altitudeInput) {
                const std::vector<std::uint8_t> altitude = pythonUInt8AltitudeVector(altitudeInput, "WeightedMorphologicalTree.altitude");
                pybind_utils::requireVectorSize(altitude, static_cast<std::size_t>(tree.topology().getNumInternalNodeSlots()), "altitude");
                tree.setAltitudeBuffer(altitude);
            },
            "Dense altitude buffer indexed by internal NodeId. Assignments validate the declared strict altitude order.")
        .def("validateAltitudeBufferShape",
             static_cast<void (WeightedMorphologicalTree<std::uint8_t>::*)() const>(&WeightedMorphologicalTree<std::uint8_t>::validateAltitudeBufferShape),
             "Validate that the dense altitude buffer covers every internal node slot.")
        .def("validateMonotoneAltitude",
             static_cast<void (WeightedMorphologicalTree<std::uint8_t>::*)() const>(&WeightedMorphologicalTree<std::uint8_t>::validateMonotoneAltitude),
             "Validate the declared strict parent-child altitude order where applicable.")
        .def("projectNodeValuesToExportedHigra", &projectNodeValuesToExportedHigraOf, "nodeValues"_a, "attributes"_a,
             "Project a node-indexed scalar or 2D attribute buffer to the compact Higra layout produced by exportHigraHierarchy().")
        .def("project_node_values_to_exported_higra", &projectNodeValuesToExportedHigraOf, "nodeValues"_a, "attributes"_a,
             "Project a node-indexed scalar or 2D attribute buffer to the compact Higra layout produced by exportHigraHierarchy().");

    bindTreeQueryApi<WeightedMorphologicalTree<std::uint8_t>>(weightedCls);
}

} // namespace mmcfilters::pybindings
