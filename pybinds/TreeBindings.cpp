#include "ModuleBindings.hpp"
#include "PythonValuedMorphologicalTree.hpp"

#include "PybindConversions.hpp"

#include "../mmcfilters/attributes/AttributeComputation.hpp"
#include "../mmcfilters/attributes/AttributeNames.hpp"
#include "../mmcfilters/trees/saliency/HierarchySaliencyMapValidation.hpp"
#include "../mmcfilters/trees/saliency/HierarchySaliencyMapProjection.hpp"
#include "../mmcfilters/trees/saliency/HierarchySaliencyMap.hpp"
#include "../mmcfilters/trees/saliency/ShapeSpaceSaliency.hpp"
#include "../mmcfilters/trees/MorphologicalTreeFactory.hpp"
#include "../mmcfilters/trees/TreeAltitudeAlgorithms.hpp"
#include "../mmcfilters/trees/ValuedMorphologicalTree.hpp"
#include "../mmcfilters/utils/RegularGridAdjacency2D.hpp"
#include "../mmcfilters/utils/Common.hpp"
#include "../mmcfilters/utils/Image.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <concepts>
#include <cstdint>
#include <memory>
#include <limits>
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

template <class Range> std::vector<PixelId> collectPixelIds(const Range& range) {
    std::vector<PixelId> ids;
    for (PixelId id : range) {
        ids.push_back(id);
    }
    return ids;
}

py::array_t<std::uint8_t> reconstructNodeMask(const MorphologicalTree& tree, NodeId nodeId) {
    if (!tree.isNode(nodeId) || !tree.isAlive(nodeId)) {
        throw std::invalid_argument("invalid NodeId for reconstruction");
    }

    ImageUInt8Ptr output = ImageUInt8::create(tree.numRows(), tree.numColumns());
    output->fill(0);
    for (PixelId pixel : tree.nodeSupport(nodeId)) {
        (*output)[pixel] = 255;
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
    const int columns = static_cast<int>(buf.shape[1]);
    if (buf.strides[1] != static_cast<py::ssize_t>(sizeof(uint8_t)) || buf.strides[0] != static_cast<py::ssize_t>(columns * sizeof(uint8_t))) {
        throw std::invalid_argument("input must be a C-contiguous 2D uint8 array");
    }
    return ImageUInt8::fromExternal(static_cast<uint8_t*>(buf.ptr), rows, columns);
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

template <std::unsigned_integral T> T pythonUnsignedAltitudeValue(py::handle value, NodeId nodeId, const char* context) {
    if (py::isinstance<py::bool_>(value) || !PyIndex_Check(value.ptr())) {
        throw std::invalid_argument(std::string(context) + " must be a non-negative integer altitude.");
    }
    py::object indexValue = py::reinterpret_steal<py::object>(PyNumber_Index(value.ptr()));
    const long long rawAltitude = py::cast<long long>(indexValue);
    if (rawAltitude < 0 || static_cast<unsigned long long>(rawAltitude) > static_cast<unsigned long long>(std::numeric_limits<T>::max())) {
        std::ostringstream oss;
        oss << context << " requires node altitudes in [0, " << std::numeric_limits<T>::max() << "]; node " << nodeId << " has altitude "
            << rawAltitude << ".";
        throw std::invalid_argument(oss.str());
    }
    return static_cast<T>(rawAltitude);
}

template <std::unsigned_integral T> std::vector<T> pythonUnsignedAltitudeVector(py::handle values, const char* context) {
    if (py::isinstance<py::array>(values)) {
        py::array array = py::reinterpret_borrow<py::array>(values);
        if (!array.dtype().is(py::dtype::of<T>())) {
            throw std::invalid_argument(std::string(context) + " numpy array dtype must match the tree altitude dtype.");
        }
        const auto buffer = array.request();
        if (buffer.ndim != 1 || buffer.strides[0] != static_cast<py::ssize_t>(sizeof(T))) {
            throw std::invalid_argument(std::string(context) + " must be a C-contiguous 1D altitude array.");
        }
        const auto* data = static_cast<const T*>(buffer.ptr);
        return std::vector<T>(data, data + buffer.shape[0]);
    }
    if (!py::isinstance<py::sequence>(values)) {
        throw std::invalid_argument(std::string(context) + " must be a sequence of non-negative integer altitudes.");
    }
    py::sequence sequence = py::reinterpret_borrow<py::sequence>(values);
    std::vector<T> altitude;
    altitude.reserve(static_cast<std::size_t>(py::len(sequence)));
    NodeId nodeId = 0;
    for (py::handle item : sequence) {
        altitude.push_back(pythonUnsignedAltitudeValue<T>(item, nodeId++, context));
    }
    return altitude;
}

template <AltitudeValue T> const MorphologicalTree& topology(const ValuedMorphologicalTree<T>& valuedTree) { return valuedTree.topology(); }

const MorphologicalTree& topology(const PythonValuedMorphologicalTree& valuedTree) { return valuedTree.topology(); }

template <AltitudeValue T>
std::shared_ptr<PythonValuedMorphologicalTree> wrapPythonValuedTree(ValuedMorphologicalTree<T>&& tree) {
    return std::make_shared<PythonValuedMorphologicalTree>(std::make_shared<ValuedMorphologicalTree<T>>(std::move(tree)));
}

template <AltitudeValue T> T nodeAltitudeValue(const ValuedMorphologicalTree<T>& valuedTree, NodeId nodeId) { return valuedTree.nodeAltitude(nodeId); }

template <AltitudeValue T> AltitudeDifference<T> nodeResidueValue(const ValuedMorphologicalTree<T>& valuedTree, NodeId nodeId) {
    return valuedTree.nodeResidue(nodeId);
}

template <AltitudeValue T> py::array_t<T> reconstructFromNodeAltitudesOf(const ValuedMorphologicalTree<T>& valuedTree) {
    return pybind_utils::toNumpy(valuedTree.reconstructFromNodeAltitudes());
}

template <std::floating_point Real>
py::array reconstructFromNodeContributionsTyped(const PythonValuedMorphologicalTree& valuedTree, py::array nodeContributions) {
    auto typed = pybind_utils::requireNodeAttributeArray<Real>(std::move(nodeContributions), valuedTree.topology(), "node_contributions");
    const auto* values = static_cast<const Real*>(typed.request().ptr);
    const std::span<const Real> contributions(values, static_cast<std::size_t>(valuedTree.topology().numInternalNodeSlots()));
    return pybind_utils::toNumpy(TreeAltitudeAlgorithms::reconstructFromNodeContributions(
        valuedTree.topology(), contributions, "ValuedMorphologicalTree.reconstruct_from_node_contributions"));
}

py::array reconstructFromNodeContributions(const PythonValuedMorphologicalTree& valuedTree, py::array nodeContributions) {
    if (pybind_utils::parseFloatingArrayDType(nodeContributions, "node_contributions") == pybind_utils::FloatingDType::Float64) {
        return reconstructFromNodeContributionsTyped<double>(valuedTree, std::move(nodeContributions));
    }
    return reconstructFromNodeContributionsTyped<float>(valuedTree, std::move(nodeContributions));
}

template <AltitudeValue T> std::pair<std::vector<NodeId>, std::vector<T>> exportHigraHierarchyOf(const ValuedMorphologicalTree<T>& valuedTree) {
    return valuedTree.exportHigraHierarchy();
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

template <AltitudeValue Altitude, std::floating_point Real>
py::array projectNodeValuesToExportedHigraTyped(const ValuedMorphologicalTree<Altitude>& valuedTree, const py::array& nodeValues, py::object attributes) {
    const auto buffer = nodeValues.request();
    const int numNodeSlots = valuedTree.topology().numInternalNodeSlots();
    const int numHigraVertices = valuedTree.topology().numPixels() + valuedTree.topology().numNodes();

    if (buffer.ndim == 1) {
        pybind_utils::require1DArray(buffer, numNodeSlots, "node_values");
        const auto parsedAttributes = parseProjectionAttributes(attributes, 1);
        const AttributeNames attrNames = makeProjectionAttributeNames(parsedAttributes);
        auto projected = AttributeComputation::projectNodeValuesToExportedHigra<Real>(
            valuedTree, attrNames, std::span<const Real>(static_cast<const Real*>(buffer.ptr), static_cast<size_t>(numNodeSlots)));
        return pybind_utils::toNumpyOwned(std::move(projected), numHigraVertices);
    }

    if (buffer.ndim == 2) {
        if (buffer.shape[0] != numNodeSlots) {
            std::ostringstream oss;
            oss << "node_values must have " << numNodeSlots << " rows, got " << buffer.shape[0];
            throw std::invalid_argument(oss.str());
        }
        const int valuesPerNode = static_cast<int>(buffer.shape[1]);
        const auto parsedAttributes = parseProjectionAttributes(attributes, valuesPerNode);
        const AttributeNames attrNames = makeProjectionAttributeNames(parsedAttributes);
        auto projected = AttributeComputation::projectNodeValuesToExportedHigra<Real>(
            valuedTree, attrNames,
            std::span<const Real>(static_cast<const Real*>(buffer.ptr), static_cast<size_t>(numNodeSlots) * static_cast<size_t>(valuesPerNode)));
        return pybind_utils::toNumpyOwned2D(std::move(projected), numHigraVertices, valuesPerNode);
    }

    throw std::invalid_argument("node_values must be a 1D or 2D float32 or float64 array");
}

template <AltitudeValue Altitude>
py::array projectNodeValuesToExportedHigraOf(const ValuedMorphologicalTree<Altitude>& valuedTree, const py::array& nodeValues, py::object attributes) {
    py::array contiguous = py::array::ensure(nodeValues, py::array::c_style);
    if (!contiguous) {
        throw std::invalid_argument("node_values must be a 1D or 2D C-contiguous float32 or float64 array");
    }

    if (contiguous.dtype().is(py::dtype::of<float>())) {
        return projectNodeValuesToExportedHigraTyped<Altitude, float>(valuedTree, contiguous, std::move(attributes));
    }
    if (contiguous.dtype().is(py::dtype::of<double>())) {
        return projectNodeValuesToExportedHigraTyped<Altitude, double>(valuedTree, contiguous, std::move(attributes));
    }

    throw std::invalid_argument("node_values must be a 1D or 2D C-contiguous float32 or float64 array");
}

py::array projectNodeValuesToExportedHigraOf(const PythonValuedMorphologicalTree& valuedTree, const py::array& nodeValues, py::object attributes) {
    return valuedTree.visit([&](const auto& native) { return projectNodeValuesToExportedHigraOf(*native, nodeValues, std::move(attributes)); });
}

template <class Value> py::dict edgeSaliencyMapToDict(EdgeSaliencyMap<Value>&& edgeMap) {
    const int numEdges = static_cast<int>(edgeMap.size());
    py::dict out;
    out["num_rows"] = edgeMap.numRows;
    out["num_columns"] = edgeMap.numColumns;
    out["adjacency_radius"] = edgeMap.adjacencyRadius;
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
        item["birth_level"] = extremum.birthLevel;
        item["death_level"] = extremum.deathLevel;
        item["extinction"] = extremum.extinction;
        extrema.append(std::move(item));
    }

    const int numNodeScores = static_cast<int>(result.nodeScores.size());
    py::dict out;
    out["extrema"] = std::move(extrema);
    out["node_scores"] = pybind_utils::toNumpyOwned(std::move(result.nodeScores), numNodeScores);
    return out;
}

template <std::floating_point Real> py::dict shapeSpaceSaliencyResultToDict(ShapeSpaceSaliencyResult<Real>&& result) {
    py::list extrema;
    for (const ShapeSpaceExtremum<Real>& extremum : result.extrema) {
        py::dict item;
        item["representative"] = extremum.representative;
        item["birth_level"] = extremum.birthLevel;
        item["death_level"] = extremum.deathLevel;
        item["extinction"] = extremum.extinction;
        extrema.append(std::move(item));
    }

    const int numNodeScores = static_cast<int>(result.nodeScores.size());
    py::dict out;
    out["extrema"] = std::move(extrema);
    out["node_scores"] = pybind_utils::toNumpyOwned(std::move(result.nodeScores), numNodeScores);
    out["edge_map"] = edgeSaliencyMapToDict(std::move(result.edgeMap));
    return out;
}

py::dict edgeContourMapToDict(EdgeContourMap&& contourMap) {
    const int numEdges = static_cast<int>(contourMap.size());
    py::dict out;
    out["num_rows"] = contourMap.numRows;
    out["num_columns"] = contourMap.numColumns;
    out["adjacency_radius"] = contourMap.adjacencyRadius;
    out["sources"] = pybind_utils::toNumpyOwned(std::move(contourMap.sources), numEdges);
    out["targets"] = pybind_utils::toNumpyOwned(std::move(contourMap.targets), numEdges);
    return out;
}

py::dict nodeContourEdgeMapToDict(NodeContourEdgeMap&& contourMap) {
    const int numEdges = static_cast<int>(contourMap.size());
    py::dict out;
    out["num_rows"] = contourMap.numRows;
    out["num_columns"] = contourMap.numColumns;
    out["adjacency_radius"] = contourMap.adjacencyRadius;
    out["sources"] = pybind_utils::toNumpyOwned(std::move(contourMap.sources), numEdges);
    out["targets"] = pybind_utils::toNumpyOwned(std::move(contourMap.targets), numEdges);
    out["nodes"] = pybind_utils::toNumpyOwned(std::move(contourMap.nodes), numEdges);
    return out;
}

py::dict incrementalNodeContourMapToDict(IncrementalNodeContourMap&& contourMap) {
    const int numEdges = static_cast<int>(contourMap.size());
    const int numOffsets = static_cast<int>(contourMap.offsets.size());
    py::dict out;
    out["num_rows"] = contourMap.numRows;
    out["num_columns"] = contourMap.numColumns;
    out["num_node_slots"] = contourMap.numNodeSlots;
    out["adjacency_radius"] = contourMap.adjacencyRadius;
    out["offsets"] = pybind_utils::toNumpyOwned(std::move(contourMap.offsets), numOffsets);
    out["sources"] = pybind_utils::toNumpyOwned(std::move(contourMap.sources), numEdges);
    out["targets"] = pybind_utils::toNumpyOwned(std::move(contourMap.targets), numEdges);
    return out;
}

EdgeSaliencyMap<double> edgeSaliencyMapFromDictPy(py::dict edgeMap, const char* context) {
    for (const char* key : {"num_rows", "num_columns", "adjacency_radius", "sources", "targets", "values"}) {
        if (!edgeMap.contains(key)) {
            throw std::invalid_argument(std::string(context) + " edge_map is missing key '" + key + "'.");
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
    parsed.numRows = py::cast<int>(edgeMap["num_rows"]);
    parsed.numColumns = py::cast<int>(edgeMap["num_columns"]);
    parsed.adjacencyRadius = py::cast<double>(edgeMap["adjacency_radius"]);
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
void validateHierarchyValuationTypedPy(const PythonValuedMorphologicalTree& valuedTree, py::array valuation, bool strict, bool nonnegative) {
    auto scores = pybind_utils::requireNodeAttributeArray<Real>(std::move(valuation), valuedTree.topology(), "valuation");
    const py::buffer_info buffer = scores.request();
    HierarchySaliencyMapValidation::validateHierarchyValuation(
        valuedTree.topology(), std::span<const Real>(static_cast<const Real*>(buffer.ptr), static_cast<std::size_t>(buffer.shape[0])),
        saliencyPolicyFromStrict(strict), saliencyRangePolicyFromNonnegative(nonnegative), "HierarchySaliencyMapValidation.validate_hierarchy_valuation");
}

void validateHierarchyValuationPy(const PythonValuedMorphologicalTree& valuedTree, py::array valuation, bool strict, bool nonnegative) {
    py::array contiguous = py::array::ensure(valuation, py::array::c_style);
    if (!contiguous) {
        throw std::invalid_argument("valuation must be a 1D C-contiguous int32, float32, or float64 array");
    }

    if (contiguous.dtype().is(py::dtype::of<int>())) {
        validateHierarchyValuationTypedPy<int>(valuedTree, contiguous, strict, nonnegative);
        return;
    }
    if (contiguous.dtype().is(py::dtype::of<float>())) {
        validateHierarchyValuationTypedPy<float>(valuedTree, contiguous, strict, nonnegative);
        return;
    }
    if (contiguous.dtype().is(py::dtype::of<double>())) {
        validateHierarchyValuationTypedPy<double>(valuedTree, contiguous, strict, nonnegative);
        return;
    }

    throw std::invalid_argument("valuation must be a 1D C-contiguous int32, float32, or float64 array");
}

template <class Real> py::array rankHierarchyValuationTypedPy(const PythonValuedMorphologicalTree& valuedTree, py::array valuation, bool strict) {
    auto scores = pybind_utils::requireNodeAttributeArray<Real>(std::move(valuation), valuedTree.topology(), "valuation");
    const py::buffer_info buffer = scores.request();
    std::vector<int> ranks = HierarchySaliencyMapValidation::rankHierarchyValuation(
        valuedTree.topology(), std::span<const Real>(static_cast<const Real*>(buffer.ptr), static_cast<std::size_t>(buffer.shape[0])),
        saliencyPolicyFromStrict(strict));
    const int numRanks = static_cast<int>(ranks.size());
    return pybind_utils::toNumpyOwned(std::move(ranks), numRanks);
}

py::array rankHierarchyValuationPy(const PythonValuedMorphologicalTree& valuedTree, py::array valuation, bool strict) {
    py::array contiguous = py::array::ensure(valuation, py::array::c_style);
    if (!contiguous) {
        throw std::invalid_argument("valuation must be a 1D C-contiguous int32, float32, or float64 array");
    }

    if (contiguous.dtype().is(py::dtype::of<int>())) {
        return rankHierarchyValuationTypedPy<int>(valuedTree, contiguous, strict);
    }
    if (contiguous.dtype().is(py::dtype::of<float>())) {
        return rankHierarchyValuationTypedPy<float>(valuedTree, contiguous, strict);
    }
    if (contiguous.dtype().is(py::dtype::of<double>())) {
        return rankHierarchyValuationTypedPy<double>(valuedTree, contiguous, strict);
    }

    throw std::invalid_argument("valuation must be a 1D C-contiguous int32, float32, or float64 array");
}

template <class Real>
py::array computeNormalizedScoresTypedPy(const PythonValuedMorphologicalTree& valuedTree, py::array valuation, bool strict, bool nonnegative) {
    auto scores = pybind_utils::requireNodeAttributeArray<Real>(std::move(valuation), valuedTree.topology(), "valuation");
    const py::buffer_info buffer = scores.request();
    std::vector<double> normalized = HierarchySaliencyMapValidation::computeNormalizedScores(
        valuedTree.topology(), std::span<const Real>(static_cast<const Real*>(buffer.ptr), static_cast<std::size_t>(buffer.shape[0])),
        saliencyPolicyFromStrict(strict), saliencyRangePolicyFromNonnegative(nonnegative));
    const int numScores = static_cast<int>(normalized.size());
    return pybind_utils::toNumpyOwned(std::move(normalized), numScores);
}

py::array computeNormalizedScoresPy(const PythonValuedMorphologicalTree& valuedTree, py::array valuation, bool strict, bool nonnegative) {
    py::array contiguous = py::array::ensure(valuation, py::array::c_style);
    if (!contiguous) {
        throw std::invalid_argument("valuation must be a 1D C-contiguous int32, float32, or float64 array");
    }

    if (contiguous.dtype().is(py::dtype::of<int>())) {
        return computeNormalizedScoresTypedPy<int>(valuedTree, contiguous, strict, nonnegative);
    }
    if (contiguous.dtype().is(py::dtype::of<float>())) {
        return computeNormalizedScoresTypedPy<float>(valuedTree, contiguous, strict, nonnegative);
    }
    if (contiguous.dtype().is(py::dtype::of<double>())) {
        return computeNormalizedScoresTypedPy<double>(valuedTree, contiguous, strict, nonnegative);
    }

    throw std::invalid_argument("valuation must be a 1D C-contiguous int32, float32, or float64 array");
}

IncrementalNodeContourMap incrementalNodeContourMapFromDictPy(py::dict contours) {
    constexpr const char* context = "HierarchySaliencyMapProjection incremental contour operation";
    for (const char* key : {"num_rows", "num_columns", "num_node_slots", "adjacency_radius", "offsets", "sources", "targets"}) {
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
    parsed.numRows = py::cast<int>(contours["num_rows"]);
    parsed.numColumns = py::cast<int>(contours["num_columns"]);
    parsed.numNodeSlots = py::cast<int>(contours["num_node_slots"]);
    parsed.adjacencyRadius = py::cast<double>(contours["adjacency_radius"]);
    parsed.offsets.assign(offsetData, offsetData + numOffsets);
    parsed.sources.assign(sourceData, sourceData + numEdges);
    parsed.targets.assign(targetData, targetData + numEdges);
    return parsed;
}

py::dict computeTopologicalLevelEdgeMapPy(const PythonValuedMorphologicalTree& valuedTree, std::optional<double> radius) {
    if (radius.has_value()) {
        RegularGridAdjacency2D adjacency =
            pybind_utils::makeRegularGridAdjacency2D(valuedTree.topology().numRows(), valuedTree.topology().numColumns(), *radius,
                                                     "HierarchySaliencyMap.compute_topological_level_edge_map");
        return edgeSaliencyMapToDict(HierarchySaliencyMap::computeTopologicalLevelEdgeMap(valuedTree.topology(), adjacency));
    }
    return edgeSaliencyMapToDict(HierarchySaliencyMap::computeTopologicalLevelEdgeMap(valuedTree.topology()));
}

py::dict computeNormalizedAltitudeEdgeMapPy(const PythonValuedMorphologicalTree& valuedTree, std::optional<double> radius) {
    if (radius.has_value()) {
        RegularGridAdjacency2D adjacency =
            pybind_utils::makeRegularGridAdjacency2D(valuedTree.topology().numRows(), valuedTree.topology().numColumns(), *radius,
                                                     "HierarchySaliencyMap.compute_normalized_altitude_edge_map");
        return valuedTree.visit([&](const auto& concreteTree) {
            return edgeSaliencyMapToDict(HierarchySaliencyMap::computeNormalizedAltitudeEdgeMap(*concreteTree, adjacency));
        });
    }
    return valuedTree.visit([](const auto& concreteTree) {
        return edgeSaliencyMapToDict(HierarchySaliencyMap::computeNormalizedAltitudeEdgeMap(*concreteTree));
    });
}

RegularGridAdjacency2D saliencyAdjacencyPy(const PythonValuedMorphologicalTree& valuedTree, std::optional<double> radius, const char* context) {
    if (radius.has_value()) {
        return pybind_utils::makeRegularGridAdjacency2D(valuedTree.topology().numRows(), valuedTree.topology().numColumns(), *radius,
                                                        context);
    }
    return HierarchySaliencyMap::requireProjectionAdjacency(valuedTree.topology(), context);
}

template <class Real>
py::dict computeSaliencyEdgeMapTypedPy(const PythonValuedMorphologicalTree& valuedTree, py::array valuation, std::optional<double> radius, bool strict,
                                       HierarchyLevelConvention levelConvention, bool validateConnectivity, const char* context) {
    auto scores = pybind_utils::requireNodeAttributeArray<Real>(std::move(valuation), valuedTree.topology(), "valuation");
    const py::buffer_info buffer = scores.request();
    const auto* scoreData = static_cast<const Real*>(buffer.ptr);
    const auto valuationSpan = std::span<const Real>(scoreData, static_cast<std::size_t>(buffer.shape[0]));

    if (radius.has_value()) {
        RegularGridAdjacency2D adjacency = pybind_utils::makeRegularGridAdjacency2D(valuedTree.topology().numRows(),
                                                                                    valuedTree.topology().numColumns(), *radius, context);
        return edgeSaliencyMapToDict(HierarchySaliencyMap::computeSaliencyEdgeMap(
            valuedTree.topology(), adjacency, valuationSpan, saliencyPolicyFromStrict(strict), levelConvention,
            validateConnectivity ? HierarchyConnectivityPolicy::ValidateConnected : HierarchyConnectivityPolicy::AssumeConnected));
    }
    return edgeSaliencyMapToDict(HierarchySaliencyMap::computeSaliencyEdgeMap(
        valuedTree.topology(), valuationSpan, saliencyPolicyFromStrict(strict), levelConvention,
        validateConnectivity ? HierarchyConnectivityPolicy::ValidateConnected : HierarchyConnectivityPolicy::AssumeConnected));
}

py::dict computeSaliencyEdgeMapPy(const PythonValuedMorphologicalTree& valuedTree, py::array valuation, std::optional<double> radius, bool strict,
                                  HierarchyLevelConvention levelConvention, bool validateConnectivity, const char* context) {
    py::array contiguous = py::array::ensure(valuation, py::array::c_style);
    if (!contiguous) {
        throw std::invalid_argument("valuation must be a 1D C-contiguous int32, float32, or float64 array");
    }

    if (contiguous.dtype().is(py::dtype::of<int>())) {
        return computeSaliencyEdgeMapTypedPy<int>(valuedTree, contiguous, radius, strict, levelConvention, validateConnectivity, context);
    }
    if (contiguous.dtype().is(py::dtype::of<float>())) {
        return computeSaliencyEdgeMapTypedPy<float>(valuedTree, contiguous, radius, strict, levelConvention, validateConnectivity, context);
    }
    if (contiguous.dtype().is(py::dtype::of<double>())) {
        return computeSaliencyEdgeMapTypedPy<double>(valuedTree, contiguous, radius, strict, levelConvention, validateConnectivity, context);
    }

    throw std::invalid_argument("valuation must be a 1D C-contiguous int32, float32, or float64 array");
}

template <class Real>
py::dict computeCanonicalRankedSaliencyEdgeMapTypedPy(const PythonValuedMorphologicalTree& valuedTree, py::array valuation,
                                                      std::optional<double> radius, bool strict, bool validateConnectivity) {
    auto scores = pybind_utils::requireNodeAttributeArray<Real>(std::move(valuation), valuedTree.topology(), "valuation");
    const py::buffer_info buffer = scores.request();
    const auto valuationSpan = std::span<const Real>(static_cast<const Real*>(buffer.ptr), static_cast<std::size_t>(buffer.shape[0]));
    const RegularGridAdjacency2D adjacency = saliencyAdjacencyPy(valuedTree, radius, "HierarchySaliencyMap.compute_canonical_ranked_saliency_edge_map");
    return edgeSaliencyMapToDict(HierarchySaliencyMap::computeCanonicalRankedSaliencyEdgeMap(
        valuedTree.topology(), adjacency, valuationSpan, saliencyPolicyFromStrict(strict),
        validateConnectivity ? HierarchyConnectivityPolicy::ValidateConnected : HierarchyConnectivityPolicy::AssumeConnected));
}

py::dict computeCanonicalRankedSaliencyEdgeMapPy(const PythonValuedMorphologicalTree& valuedTree, py::array valuation, std::optional<double> radius,
                                                 bool strict, bool validateConnectivity) {
    py::array contiguous = py::array::ensure(valuation, py::array::c_style);
    if (!contiguous) {
        throw std::invalid_argument("valuation must be a 1D C-contiguous int32, float32, or float64 array");
    }
    if (contiguous.dtype().is(py::dtype::of<int>())) {
        return computeCanonicalRankedSaliencyEdgeMapTypedPy<int>(valuedTree, contiguous, radius, strict, validateConnectivity);
    }
    if (contiguous.dtype().is(py::dtype::of<float>())) {
        return computeCanonicalRankedSaliencyEdgeMapTypedPy<float>(valuedTree, contiguous, radius, strict, validateConnectivity);
    }
    if (contiguous.dtype().is(py::dtype::of<double>())) {
        return computeCanonicalRankedSaliencyEdgeMapTypedPy<double>(valuedTree, contiguous, radius, strict, validateConnectivity);
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
py::dict computeShapeSpaceExtinctionValuesTypedPy(const PythonValuedMorphologicalTree& valuedTree, py::array attribute,
                                                  ShapeSpaceExtremaPolarity polarity) {
    py::array_t<Real, py::array::c_style> attributeOwner;
    const std::span<const Real> attributeSpan = shapeSpaceDenseNodeSpan<Real>(std::move(attribute), valuedTree.topology(), attributeOwner, "attribute");
    return shapeSpaceExtinctionResultToDict(ShapeSpaceSaliency::computeExtinctionValues(valuedTree.topology(), attributeSpan, polarity));
}

py::dict computeShapeSpaceExtinctionValuesPy(const PythonValuedMorphologicalTree& valuedTree, py::array attribute, ShapeSpaceExtremaPolarity polarity) {
    py::array contiguous = py::array::ensure(attribute, py::array::c_style);
    if (!contiguous) {
        throw std::invalid_argument("attribute must be a 1D C-contiguous float32 or float64 array");
    }
    if (contiguous.dtype().is(py::dtype::of<float>())) {
        return computeShapeSpaceExtinctionValuesTypedPy<float>(valuedTree, contiguous, polarity);
    }
    if (contiguous.dtype().is(py::dtype::of<double>())) {
        return computeShapeSpaceExtinctionValuesTypedPy<double>(valuedTree, contiguous, polarity);
    }
    throw std::invalid_argument("attribute must be a 1D C-contiguous float32 or float64 array");
}

template <std::floating_point Real>
py::dict projectShapeSpaceContourScoresTypedPy(const PythonValuedMorphologicalTree& valuedTree, py::array scores, std::optional<double> radius) {
    py::array_t<Real, py::array::c_style> scoreOwner;
    const std::span<const Real> scoreSpan = shapeSpaceDenseNodeSpan<Real>(std::move(scores), valuedTree.topology(), scoreOwner, "node_scores");

    if (radius.has_value()) {
        RegularGridAdjacency2D adjacency = pybind_utils::makeRegularGridAdjacency2D(
            valuedTree.topology().numRows(), valuedTree.topology().numColumns(), *radius, "ShapeSpaceSaliency.project_contour_scores");
        return edgeSaliencyMapToDict(ShapeSpaceSaliency::projectContourScores(valuedTree.topology(), scoreSpan, adjacency));
    }
    return edgeSaliencyMapToDict(ShapeSpaceSaliency::projectContourScores(valuedTree.topology(), scoreSpan));
}

py::dict projectShapeSpaceContourScoresPy(const PythonValuedMorphologicalTree& valuedTree, py::array scores, std::optional<double> radius) {
    py::array contiguous = py::array::ensure(scores, py::array::c_style);
    if (!contiguous) {
        throw std::invalid_argument("node_scores must be a 1D C-contiguous float32 or float64 array");
    }
    if (contiguous.dtype().is(py::dtype::of<float>())) {
        return projectShapeSpaceContourScoresTypedPy<float>(valuedTree, contiguous, radius);
    }
    if (contiguous.dtype().is(py::dtype::of<double>())) {
        return projectShapeSpaceContourScoresTypedPy<double>(valuedTree, contiguous, radius);
    }
    throw std::invalid_argument("node_scores must be a 1D C-contiguous float32 or float64 array");
}

template <std::floating_point Real>
py::dict computeShapeSpaceSaliencyTypedPy(const PythonValuedMorphologicalTree& valuedTree, py::array attribute, ShapeSpaceExtremaPolarity polarity,
                                          std::optional<double> radius) {
    py::array_t<Real, py::array::c_style> attributeOwner;
    const std::span<const Real> attributeSpan = shapeSpaceDenseNodeSpan<Real>(std::move(attribute), valuedTree.topology(), attributeOwner, "attribute");

    if (radius.has_value()) {
        RegularGridAdjacency2D adjacency = pybind_utils::makeRegularGridAdjacency2D(
            valuedTree.topology().numRows(), valuedTree.topology().numColumns(), *radius, "ShapeSpaceSaliency.compute");
        return shapeSpaceSaliencyResultToDict(ShapeSpaceSaliency::compute(valuedTree.topology(), attributeSpan, polarity, adjacency));
    }
    return shapeSpaceSaliencyResultToDict(ShapeSpaceSaliency::compute(valuedTree.topology(), attributeSpan, polarity));
}

py::dict computeShapeSpaceSaliencyPy(const PythonValuedMorphologicalTree& valuedTree, py::array attribute, ShapeSpaceExtremaPolarity polarity,
                                     std::optional<double> radius) {
    py::array contiguous = py::array::ensure(attribute, py::array::c_style);
    if (!contiguous) {
        throw std::invalid_argument("attribute must be a 1D C-contiguous float32 or float64 array");
    }
    if (contiguous.dtype().is(py::dtype::of<float>())) {
        return computeShapeSpaceSaliencyTypedPy<float>(valuedTree, contiguous, polarity, radius);
    }
    if (contiguous.dtype().is(py::dtype::of<double>())) {
        return computeShapeSpaceSaliencyTypedPy<double>(valuedTree, contiguous, polarity, radius);
    }
    throw std::invalid_argument("attribute must be a 1D C-contiguous float32 or float64 array");
}

py::dict thresholdCutPy(py::dict edgeMap, double threshold) {
    EdgeSaliencyMap<double> parsed = edgeSaliencyMapFromDictPy(std::move(edgeMap), "HierarchySaliencyMapProjection.threshold_cut");
    return edgeContourMapToDict(HierarchySaliencyMapProjection::thresholdCut(parsed, threshold));
}

py::array edgeMapToPixelImagePy(py::dict edgeMap, EdgeToPixelReducer reducer) {
    EdgeSaliencyMap<double> parsed = edgeSaliencyMapFromDictPy(std::move(edgeMap), "HierarchySaliencyMapProjection.edge_map_to_pixel_image");
    return pybind_utils::toNumpy(HierarchySaliencyMapProjection::edgeMapToPixelImage(parsed, reducer));
}

py::dict nodeContourEdgesPy(const PythonValuedMorphologicalTree& valuedTree, std::optional<double> radius) {
    if (radius.has_value()) {
        RegularGridAdjacency2D adjacency =
            pybind_utils::makeRegularGridAdjacency2D(valuedTree.topology().numRows(), valuedTree.topology().numColumns(), *radius,
                                                     "HierarchySaliencyMapProjection.node_contour_edges");
        return valuedTree.visit([&](const auto& concreteTree) {
            return nodeContourEdgeMapToDict(HierarchySaliencyMapProjection::nodeContourEdges(*concreteTree, adjacency));
        });
    }
    return valuedTree.visit([](const auto& concreteTree) {
        return nodeContourEdgeMapToDict(HierarchySaliencyMapProjection::nodeContourEdges(*concreteTree));
    });
}

py::dict computeIncrementalNodeContoursPy(const PythonValuedMorphologicalTree& valuedTree, std::optional<double> radius) {
    if (radius.has_value()) {
        RegularGridAdjacency2D adjacency =
            pybind_utils::makeRegularGridAdjacency2D(valuedTree.topology().numRows(), valuedTree.topology().numColumns(), *radius,
                                                     "HierarchySaliencyMapProjection.compute_incremental_node_contours");
        return valuedTree.visit([&](const auto& concreteTree) {
            return incrementalNodeContourMapToDict(HierarchySaliencyMapProjection::computeIncrementalNodeContours(*concreteTree, adjacency));
        });
    }
    return valuedTree.visit([](const auto& concreteTree) {
        return incrementalNodeContourMapToDict(HierarchySaliencyMapProjection::computeIncrementalNodeContours(*concreteTree));
    });
}

template <std::floating_point Real> py::dict projectNodeValuationOnIncrementalContoursTypedPy(py::dict contours, py::array nodeValuation) {
    IncrementalNodeContourMap parsed = incrementalNodeContourMapFromDictPy(std::move(contours));
    auto valuation = py::array_t<Real, py::array::c_style>::ensure(nodeValuation);
    if (!valuation) {
        throw std::invalid_argument("node_valuation must be a 1D C-contiguous float32 or float64 array.");
    }
    const py::buffer_info info = valuation.request();
    pybind_utils::require1DArray(info, parsed.numNodeSlots, "node_valuation");
    if (info.strides[0] != static_cast<py::ssize_t>(sizeof(Real))) {
        throw std::invalid_argument("node_valuation must be C-contiguous.");
    }
    return edgeSaliencyMapToDict(HierarchySaliencyMapProjection::projectNodeValuation(
        parsed, std::span<const Real>(static_cast<const Real*>(info.ptr), static_cast<std::size_t>(info.shape[0]))));
}

py::dict projectNodeValuationOnIncrementalContoursPy(py::dict contours, py::array nodeValuation) {
    py::array contiguous = py::array::ensure(nodeValuation, py::array::c_style);
    if (!contiguous) {
        throw std::invalid_argument("node_valuation must be a 1D C-contiguous float32 or float64 array");
    }
    if (contiguous.dtype().is(py::dtype::of<float>())) {
        return projectNodeValuationOnIncrementalContoursTypedPy<float>(std::move(contours), contiguous);
    }
    if (contiguous.dtype().is(py::dtype::of<double>())) {
        return projectNodeValuationOnIncrementalContoursTypedPy<double>(std::move(contours), contiguous);
    }
    throw std::invalid_argument("node_valuation must be a 1D C-contiguous float32 or float64 array");
}

template <std::floating_point Real> py::dict thresholdIncrementalContoursByNodeValuationTypedPy(py::dict contours, py::array nodeValuation, double threshold) {
    IncrementalNodeContourMap parsed = incrementalNodeContourMapFromDictPy(std::move(contours));
    auto valuation = py::array_t<Real, py::array::c_style>::ensure(nodeValuation);
    if (!valuation) {
        throw std::invalid_argument("node_valuation must be a 1D C-contiguous float32 or float64 array.");
    }
    const py::buffer_info info = valuation.request();
    pybind_utils::require1DArray(info, parsed.numNodeSlots, "node_valuation");
    if (info.strides[0] != static_cast<py::ssize_t>(sizeof(Real))) {
        throw std::invalid_argument("node_valuation must be C-contiguous.");
    }
    return edgeContourMapToDict(HierarchySaliencyMapProjection::thresholdByNodeValuation(
        parsed, std::span<const Real>(static_cast<const Real*>(info.ptr), static_cast<std::size_t>(info.shape[0])), threshold));
}

py::dict thresholdIncrementalContoursByNodeValuationPy(py::dict contours, py::array nodeValuation, double threshold) {
    py::array contiguous = py::array::ensure(nodeValuation, py::array::c_style);
    if (!contiguous) {
        throw std::invalid_argument("node_valuation must be a 1D C-contiguous float32 or float64 array");
    }
    if (contiguous.dtype().is(py::dtype::of<float>())) {
        return thresholdIncrementalContoursByNodeValuationTypedPy<float>(std::move(contours), contiguous, threshold);
    }
    if (contiguous.dtype().is(py::dtype::of<double>())) {
        return thresholdIncrementalContoursByNodeValuationTypedPy<double>(std::move(contours), contiguous, threshold);
    }
    throw std::invalid_argument("node_valuation must be a 1D C-contiguous float32 or float64 array");
}

template <class TreeLike, class PyClass> void bindTreeQueryApi(PyClass& cls) {
    cls.def_property_readonly(
           "num_internal_node_slots", [](TreeLike& self) { return topology(self).numInternalNodeSlots(); }, "Size of the dense internal NodeId slot domain.")
        .def_property_readonly(
            "num_pixels", [](TreeLike& self) { return topology(self).numPixels(); }, "Number of pixels in the finite tree domain.")
        .def_property_readonly(
            "num_higra_nodes", [](TreeLike& self) { return topology(self).getNumHigraNodes(); }, "Size of the preserved imported Higra node-id domain.")
        .def_property_readonly(
            "root", [](TreeLike& self) { return topology(self).root(); }, "Current root node id.")
        .def_property_readonly(
            "num_free_node_slots", [](TreeLike& self) { return topology(self).getNumFreeNodeSlots(); }, "Number of currently free internal-node slots.")
        .def_property_readonly(
            "num_leaf_nodes", [](TreeLike& self) { return topology(self).numLeafNodes(); }, "Number of alive leaf nodes.")
        .def_property_readonly(
            "alive_node_ids", [](TreeLike& self) { return collectNodeIds(topology(self).aliveNodeIds()); },
            "Alive internal-node ids in the dense node-id domain.")
        .def_property_readonly(
            "leaves", [](TreeLike& self) { return topology(self).leaves(); }, "Alive leaf node ids in the dense node-id domain.")
        .def(
            "children", [](TreeLike& self, NodeId nodeId) { return collectNodeIds(topology(self).children(nodeId)); }, "node_id"_a,
            "Return the direct children of a node in the dense node-id domain.")
        .def(
            "num_descendants", [](TreeLike& self, NodeId nodeId) { return topology(self).numDescendants(nodeId); }, "node_id"_a,
            "Return the number of strict descendants of node_id.")
        .def(
            "num_siblings", [](TreeLike& self, NodeId nodeId) { return topology(self).numSiblings(nodeId); }, "node_id"_a,
            "Return the number of siblings of node_id.")
        .def(
            "proper_part_cardinality", [](TreeLike& self, NodeId nodeId) { return topology(self).properPartCardinality(nodeId); }, "node_id"_a,
            "Return the number of pixels in the proper part of node_id.")
        .def(
            "dfs_entry_index", [](TreeLike& self, NodeId nodeId) { return topology(self).dfsEntryIndex(nodeId); }, "node_id"_a,
            "Return the zero-based entry-event index of node_id in the interleaved DFS event sequence.")
        .def(
            "dfs_exit_index", [](TreeLike& self, NodeId nodeId) { return topology(self).dfsExitIndex(nodeId); }, "node_id"_a,
            "Return the zero-based exit-event index of node_id in the interleaved DFS event sequence.")
        .def(
            "proper_part", [](TreeLike& self, NodeId nodeId) { return collectPixelIds(topology(self).properPart(nodeId)); }, "node_id"_a,
            "Return the pixels in the proper part of a node.")
        .def(
            "node_support",
            [](TreeLike& self, NodeId nodeId) {
                auto range = topology(self).nodeSupport(nodeId);
                return py::make_iterator(range.begin(), range.end());
            },
            py::keep_alive<0, 1>(), "node_id"_a, "Iterate over all pixels in the support represented by node_id.")
        .def(
            "reconstruct_node", [](TreeLike& self, NodeId nodeId) { return reconstructNodeMask(topology(self), nodeId); }, "node_id"_a,
            "Reconstruct a binary mask for the support represented by node_id.")
        .def(
            "post_order",
            [](TreeLike& self, std::optional<NodeId> rootNodeId) {
                return rootNodeId.has_value() ? collectNodeIds(topology(self).postOrder(*rootNodeId))
                                              : collectNodeIds(topology(self).postOrder());
            },
            "root_node_id"_a = std::nullopt, "Return a post-order traversal schedule under root_node_id, or under the tree root.")
        .def(
            "breadth_first_traversal",
            [](TreeLike& self, std::optional<NodeId> rootNodeId) {
                return rootNodeId.has_value() ? collectNodeIds(topology(self).breadthFirstTraversal(*rootNodeId))
                                              : collectNodeIds(topology(self).breadthFirstTraversal());
            },
            "root_node_id"_a = std::nullopt, "Return breadth-first traversal node ids under root_node_id, or under the tree root.")
        .def(
            "ancestors", [](TreeLike& self, NodeId nodeId) { return collectNodeIds(topology(self).ancestors(nodeId)); }, "node_id"_a,
            "Return node_id followed by its proper ancestors up to and including the root.")
        .def(
            "lowest_common_ancestor", [](TreeLike& self, NodeId u, NodeId v) { return topology(self).lowestCommonAncestor(u, v); }, "u"_a, "v"_a,
            "Return the inclusion-lowest common ancestor of two live nodes.")
        .def(
            "path_between_nodes",
            [](TreeLike& self, NodeId sourceNodeId, NodeId targetNodeId) {
                return collectNodeIds(topology(self).getPathBetweenNodes(sourceNodeId, targetNodeId));
            },
            "source_node_id"_a, "target_node_id"_a, "Return the upward path from source_node_id toward target_node_id.")
        .def(
            "subtree_nodes", [](TreeLike& self, NodeId nodeId) { return collectNodeIds(topology(self).subtreeNodes(nodeId)); }, "node_id"_a,
            "Return node_id and all of its descendants in pre-order.")
        .def(
            "descendants", [](TreeLike& self, NodeId nodeId) { return collectNodeIds(topology(self).descendants(nodeId)); }, "node_id"_a,
            "Return all strict descendant node ids of node_id in pre-order.")
        .def(
            "parent", [](TreeLike& self, NodeId nodeId) { return topology(self).parent(nodeId); }, "node_id"_a,
            "Return the parent of node_id; the root is its own parent.")
        .def(
            "smallest_node", [](TreeLike& self, PixelId pixel) { return topology(self).smallestNode(pixel); }, "pixel"_a,
            "Return the inclusion-smallest node containing pixel.")
        .def_property_readonly(
            "smallest_node_map",
            [](TreeLike& self) {
                const auto map = topology(self).smallestNodeMap();
                return std::vector<NodeId>(map.begin(), map.end());
            },
            "Pixel-indexed smallest-node map.")
        .def(
            "higra_node_id", [](TreeLike& self, NodeId nodeId) { return topology(self).getHigraNodeId(nodeId); }, "node_id"_a,
            "Return the preserved imported Higra node id for a live internal NodeId, or InvalidNode.")
        .def(
            "num_children", [](TreeLike& self, NodeId nodeId) { return topology(self).numChildren(nodeId); }, "node_id"_a,
            "Return the number of direct children of node_id.")
        .def(
            "first_child", [](TreeLike& self, NodeId nodeId) { return topology(self).getFirstChild(nodeId); }, "node_id"_a,
            "Return the first child of `node_id`, or `InvalidNode`.")
        .def(
            "next_sibling", [](TreeLike& self, NodeId nodeId) { return topology(self).getNextSibling(nodeId); }, "node_id"_a,
            "Return the next sibling of `node_id`, or `InvalidNode`.")
        .def(
            "is_node", [](TreeLike& self, NodeId nodeId) { return topology(self).isNode(nodeId); }, "node_id"_a,
            "Return true when `node_id` is in the internal-node slot domain.")
        .def(
            "is_pixel", [](TreeLike& self, PixelId pixel) { return topology(self).isPixel(pixel); }, "pixel"_a,
            "Return true when pixel is in the finite tree domain.")
        .def(
            "is_alive", [](TreeLike& self, NodeId nodeId) { return topology(self).isAlive(nodeId); }, "node_id"_a,
            "Return true when `node_id` is an alive internal node.")
        .def(
            "is_root", [](TreeLike& self, NodeId nodeId) { return topology(self).isRoot(nodeId); }, "node_id"_a, "Return true when node_id is the current root.")
        .def(
            "is_leaf", [](TreeLike& self, NodeId nodeId) { return topology(self).isLeaf(nodeId); }, "node_id"_a,
            "Return true when `node_id` is an alive leaf node.")
        .def(
            "has_empty_proper_part", [](TreeLike& self, NodeId nodeId) { return topology(self).hasEmptyProperPart(nodeId); }, "node_id"_a,
            "Return true when node_id has an empty proper part.")
        .def(
            "is_tree_of_partial_partitions", [](TreeLike& self) { return topology(self).isTreeOfPartialPartitions(); },
            "Whether every live node has a non-empty proper part.")
        .def(
            "validate_tree_of_partial_partitions", [](TreeLike& self) { topology(self).validateTreeOfPartialPartitions(); },
            "Validate the connected-subset tree and require every live node to have a non-empty proper part.")
        .def(
            "has_child", [](TreeLike& self, NodeId parentId, NodeId childId) { return topology(self).hasChild(parentId, childId); }, "parent_id"_a, "child_id"_a,
            "Return true when child_id is a direct child of parent_id.")
        .def(
            "prune_node", [](TreeLike& self, NodeId nodeId) { self.pruneNode(nodeId); }, "node_id"_a,
            "Prune one node from the topology, preserving a valid rooted tree.")
        .def(
            "merge_node_into_parent", [](TreeLike& self, NodeId nodeId) { self.mergeNodeIntoParent(nodeId); }, "node_id"_a,
            "Merge `nodeId` into its parent and reattach descendants/proper parts.")
        .def_property_readonly(
            "kind", [](TreeLike& self) { return topology(self).kind(); }, "Optional descriptive hierarchy-family label.")
        .def_property_readonly(
            "node_altitude_order", [](TreeLike& self) { return topology(self).nodeAltitudeOrder(); }, "Global parent-to-child altitude ordering capability.")
        .def_property_readonly(
            "semantics", [](TreeLike& self) -> const MorphologicalTreeSemantics& { return topology(self).semantics(); },
            py::return_value_policy::reference_internal, "Immutable scientific semantics retained by the tree.")
        .def_property_readonly(
            "construction_context", [](TreeLike& self) -> const MorphologicalTreeConstructionContext& { return topology(self).constructionContext(); },
            py::return_value_policy::reference_internal, "Typed construction context retained by the tree.")
        .def_property_readonly(
            "shared_adjacency_context", [](TreeLike& self) { return topology(self).sharedAdjacencyContext(); },
            py::return_value_policy::reference_internal, "Shared-adjacency context, or None.")
        .def_property_readonly(
            "saturated_residual_context", [](TreeLike& self) { return topology(self).saturatedResidualContext(); },
            py::return_value_policy::reference_internal, "Saturated-residual context, or None.")
        .def_property_readonly(
            "topographic_convention", [](TreeLike& self) { return topology(self).topographicConvention(); },
            py::return_value_policy::reference_internal, "Topographic convention, or None.")
        .def_property_readonly(
            "has_grid_domain_2d", [](TreeLike& self) { return topology(self).hasGridDomain2D(); }, "Whether pixel ids have a regular row/column layout.")
        .def_property_readonly(
            "grid_domain_2d", [](TreeLike& self) { return topology(self).gridDomain2D(); }, "Optional regular 2D pixel domain.")
        .def_property_readonly(
            "num_rows", [](TreeLike& self) { return topology(self).numRows(); }, "Number of rows in the active pixel domain.")
        .def_property_readonly(
            "num_columns", [](TreeLike& self) { return topology(self).numColumns(); }, "Number of columns in the active pixel domain.")
        .def_property_readonly("num_nodes", [](TreeLike& self) { return topology(self).numNodes(); }, "Number of currently alive internal nodes.");

    if constexpr (requires(TreeLike& self, NodeId nodeId) {
                      typename TreeLike::AltitudeType;
                      self.nodeAltitude(nodeId);
                      self.nodeResidue(nodeId);
                      self.reconstructFromNodeAltitudes();
                      self.exportHigraHierarchy();
                  }) {
        cls.def(
               "node_altitude", [](TreeLike& self, NodeId nodeId) { return nodeAltitudeValue(self, nodeId); }, "node_id"_a,
               "Return the altitude associated with node_id.")
            .def(
                "node_residue", [](TreeLike& self, NodeId nodeId) { return nodeResidueValue(self, nodeId); }, "node_id"_a,
                "Return the residue between node_id and its parent.")
            .def(
                "reconstruct_from_node_altitudes", [](TreeLike& self) { return reconstructFromNodeAltitudesOf(self); },
                "Reconstruct the current tree from node altitudes into a 2D image using the attached GridDomain2D.")
            .def(
                "export_higra_hierarchy", [](TreeLike& self) { return exportHigraHierarchyOf(self); },
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
    py::enum_<NodeIdSpace>(m, "NodeIdSpace", py::module_local(false))
        .value("MORPHOLOGICAL_TREE", NodeIdSpace::MorphologicalTree)
        .value("HIGRA", NodeIdSpace::Higra)
        .export_values();

    py::enum_<MorphologicalTreeKind>(m, "MorphologicalTreeKind", py::module_local(false))
        .value("GENERIC", MorphologicalTreeKind::Generic)
        .value("MAX_TREE", MorphologicalTreeKind::MaxTree)
        .value("MIN_TREE", MorphologicalTreeKind::MinTree)
        .value("TREE_OF_SHAPES", MorphologicalTreeKind::TreeOfShapes)
        .value("UNRESTRICTED_RESIDUAL_TREE", MorphologicalTreeKind::UnrestrictedResidualTree)
        .value("SATURATED_RESIDUAL_TREE", MorphologicalTreeKind::SaturatedResidualTree)
        .export_values();

    py::enum_<sdrt::Polarity>(m, "Polarity", py::module_local(false))
        .value("MAXIMUM", sdrt::Polarity::Maximum)
        .value("MINIMUM", sdrt::Polarity::Minimum)
        .export_values();

    py::class_<sdrt::SpatialOrder>(m, "SpatialOrder", py::module_local(false),
                                   "A total order over one dense pixel domain.")
        .def(py::init<std::vector<PixelId>>(), "pixels_in_order"_a)
        .def_property_readonly("is_row_major", &sdrt::SpatialOrder::isRowMajor)
        .def("precedes", &sdrt::SpatialOrder::precedes, "lhs"_a, "rhs"_a)
        .def("spatial_minimum", [](const sdrt::SpatialOrder& self, const std::vector<PixelId>& support) {
            return self.spatialMinimum(support);
        }, "support"_a);

    py::class_<sdrt::RowMajorSpatialOrder, sdrt::SpatialOrder>(m, "RowMajorSpatialOrder", py::module_local(false),
                                                                "The default row-major pixel order.")
        .def(py::init<>());

    py::class_<sdrt::SelfDualResidualKey>(m, "SelfDualResidualKey", py::module_local(false))
        .def(py::init<std::size_t, PixelId>(), "support_cardinality"_a, "spatial_minimum"_a)
        .def_readwrite("support_cardinality", &sdrt::SelfDualResidualKey::supportCardinality)
        .def_readwrite("spatial_minimum", &sdrt::SelfDualResidualKey::spatialMinimum);

    py::class_<sdrt::SelfDualResidualOrder>(m, "SelfDualResidualOrder", py::module_local(false))
        .def(py::init<sdrt::SpatialOrder>(), "spatial_order"_a = sdrt::RowMajorSpatialOrder{})
        .def("compare_residual_candidates", &sdrt::SelfDualResidualOrder::compareResidualCandidates, "lhs"_a, "rhs"_a);

    py::class_<sdrt::SelfDualResidualSchedule>(m, "SelfDualResidualSchedule", py::module_local(false))
        .def(py::init<sdrt::SpatialOrder>(), "spatial_order"_a = sdrt::RowMajorSpatialOrder{})
        .def("select_residual_candidate", [](const sdrt::SelfDualResidualSchedule& self,
                                              const std::vector<sdrt::SelfDualResidualKey>& keys) {
            return self.selectResidualCandidate(keys);
        }, "residual_keys"_a);

    py::enum_<NodeAltitudeOrder>(m, "NodeAltitudeOrder", py::module_local(false))
        .value("INCREASING", NodeAltitudeOrder::Increasing, "Every child altitude is strictly greater than its parent altitude.")
        .value("DECREASING", NodeAltitudeOrder::Decreasing, "Every child altitude is strictly smaller than its parent altitude.")
        .value("UNCONSTRAINED", NodeAltitudeOrder::Unconstrained, "No global parent-child altitude order is declared.")
        .export_values();

    py::class_<GridDomain2D>(m, "GridDomain2D", py::module_local(false), "Optional row/column layout attached to a pixel domain.")
        .def(py::init<int, int>(), "rows"_a, "columns"_a)
        .def_readonly("rows", &GridDomain2D::rows)
        .def_readonly("columns", &GridDomain2D::columns);

    py::class_<NoConstructionContext>(m, "NoConstructionContext", py::module_local(false), "Explicitly unavailable construction provenance.")
        .def(py::init<>());

    py::class_<SharedAdjacencyContext>(m, "SharedAdjacencyContext", py::module_local(false),
                                       "One adjacency shared by both construction polarities.")
        .def(py::init<RegularGridAdjacency2D>(), "adjacency"_a)
        .def_readonly("adjacency", &SharedAdjacencyContext::adjacency);

    py::class_<SaturatedResidualContext>(m, "SaturatedResidualContext", py::module_local(false),
                                         "Adjacency and infinity pixel of a saturated residual construction.")
        .def(py::init<RegularGridAdjacency2D, PixelId>(), "adjacency"_a, "infinity_pixel"_a = PixelId{0})
        .def_readonly("adjacency", &SaturatedResidualContext::adjacency)
        .def_readonly("infinity_pixel", &SaturatedResidualContext::infinityPixel);

    py::class_<ComplementaryAdjacencies>(m, "ComplementaryAdjacencies", py::module_local(false),
                                         "Minimum and maximum adjacencies of a complementary-grid immersion.")
        .def(py::init<RegularGridAdjacency2D, RegularGridAdjacency2D>(), "min_adjacency"_a, "max_adjacency"_a)
        .def_readonly("min_adjacency", &ComplementaryAdjacencies::minAdjacency)
        .def_readonly("max_adjacency", &ComplementaryAdjacencies::maxAdjacency);

    py::class_<SelfDualSpanImmersion>(m, "SelfDualSpanImmersion", py::module_local(false), "Self-dual span-valued immersion.")
        .def(py::init<>());

    py::class_<ComplementaryGridImmersion>(m, "ComplementaryGridImmersion", py::module_local(false), "Complementary-grid immersion.")
        .def(py::init<ComplementaryAdjacencies>(), "complementary_adjacencies"_a)
        .def_readonly("complementary_adjacencies", &ComplementaryGridImmersion::complementaryAdjacencies);

    py::enum_<TopographicDomainExtension>(m, "TopographicDomainExtension", py::module_local(false))
        .value("EXTERIOR_RING", TopographicDomainExtension::ExteriorRing)
        .value("NONE", TopographicDomainExtension::None)
        .export_values();

    py::class_<TopographicConvention>(m, "TopographicConvention", py::module_local(false),
                                      "Complete discrete convention retained by a tree of shapes.")
        .def(py::init<TreeOfShapesImmersion, TopographicDomainExtension, PixelId>(), "immersion"_a = SelfDualSpanImmersion{},
             "domain_extension"_a = TopographicDomainExtension::ExteriorRing, "infinity_pixel"_a = PixelId{0})
        .def_readonly("immersion", &TopographicConvention::immersion)
        .def_readonly("domain_extension", &TopographicConvention::domainExtension)
        .def_readonly("infinity_pixel", &TopographicConvention::infinityPixel);

    py::class_<MorphologicalTreeSemantics>(m, "MorphologicalTreeSemantics", py::module_local(false),
                                   "Immutable scientific semantics attached to a morphological tree.")
        .def(py::init([](MorphologicalTreeKind kind, NodeAltitudeOrder nodeAltitudeOrder,
                         MorphologicalTreeConstructionContext constructionContext) {
                 MorphologicalTreeSemantics semantics{kind, nodeAltitudeOrder, std::move(constructionContext)};
                 validateMorphologicalTreeSemantics(semantics);
                 return semantics;
             }),
             "kind"_a = MorphologicalTreeKind::Generic, "node_altitude_order"_a = NodeAltitudeOrder::Unconstrained,
             "construction_context"_a = NoConstructionContext{})
        .def_readonly("kind", &MorphologicalTreeSemantics::kind)
        .def_readonly("node_altitude_order", &MorphologicalTreeSemantics::nodeAltitudeOrder)
        .def_readonly("construction_context", &MorphologicalTreeSemantics::constructionContext);

    py::enum_<EdgeToPixelReducer>(m, "EdgeToPixelReducer", py::module_local(false))
        .value("MAX", EdgeToPixelReducer::Max)
        .value("MEAN", EdgeToPixelReducer::Mean)
        .export_values();

    py::enum_<HierarchyLevelConvention>(m, "HierarchyLevelConvention", py::module_local(false))
        .value("EDGE_SALIENCY_VALUE", HierarchyLevelConvention::EdgeSaliencyValue)
        .value("PARTITION_APPEARANCE_LEVEL", HierarchyLevelConvention::PartitionAppearanceLevel)
        .export_values();

    py::enum_<ShapeSpaceExtremaPolarity>(m, "ShapeSpaceExtremaPolarity", py::module_local(false))
        .value("MINIMA", ShapeSpaceExtremaPolarity::Minima)
        .value("MAXIMA", ShapeSpaceExtremaPolarity::Maxima)
        .export_values();

    auto valuedTreeCls = py::class_<PythonValuedMorphologicalTree, std::shared_ptr<PythonValuedMorphologicalTree>>(
        m, "ValuedMorphologicalTree", py::module_local(false),
        "Wrapper pairing MorphologicalTree topology with an external dense altitude buffer. "
        "Imports can preserve an original Higra node-id domain; exports always create a new compact Higra domain.");

    py::class_<MorphologicalTreeFactory>(m, "MorphologicalTreeFactory", py::module_local(false))
        .def_static(
            "create_max_tree",
            [](UInt8InputArray input, double radius) {
                const double validatedRadius = pybind_utils::requireAdjacencyRadius(radius, "MorphologicalTreeFactory.create_max_tree");
                return wrapPythonValuedTree(MorphologicalTreeFactory::createMaxTree(imageFromArray(input), validatedRadius));
            },
            "input"_a, "radius"_a = 1.5, "Create a max-tree from a 2D C-contiguous `np.uint8` image.")
        .def_static(
            "create_max_tree",
            [](UInt8InputArray input, RegularGridAdjacency2D adjacency) {
                return wrapPythonValuedTree(MorphologicalTreeFactory::createMaxTree(imageFromArray(input), std::move(adjacency)));
            },
            "input"_a, "adjacency"_a, "Create a max-tree with an explicit regular-grid 2D adjacency.")
        .def_static(
            "create_min_tree",
            [](UInt8InputArray input, double radius) {
                const double validatedRadius = pybind_utils::requireAdjacencyRadius(radius, "MorphologicalTreeFactory.create_min_tree");
                return wrapPythonValuedTree(MorphologicalTreeFactory::createMinTree(imageFromArray(input), validatedRadius));
            },
            "input"_a, "radius"_a = 1.5, "Create a min-tree from a 2D C-contiguous `np.uint8` image.")
        .def_static(
            "create_min_tree",
            [](UInt8InputArray input, RegularGridAdjacency2D adjacency) {
                return wrapPythonValuedTree(MorphologicalTreeFactory::createMinTree(imageFromArray(input), std::move(adjacency)));
            },
            "input"_a, "adjacency"_a, "Create a min-tree with an explicit regular-grid 2D adjacency.")
        .def_static(
            "create_unrestricted_residual_tree",
            [](UInt8InputArray input, double radius, sdrt::SpatialOrder spatialOrder) {
                const double validatedRadius =
                    pybind_utils::requireAdjacencyRadius(radius, "MorphologicalTreeFactory.create_unrestricted_residual_tree");
                return wrapPythonValuedTree(MorphologicalTreeFactory::createUnrestrictedResidualTree(
                    imageFromArray(input), validatedRadius, sdrt::UnrestrictedResidualTreeOptions{std::move(spatialOrder)}));
            },
            "input"_a, "radius"_a = 1.5, "spatial_order"_a = sdrt::RowMajorSpatialOrder{},
            "Create the unrestricted residual tree with synchronized max-tree and min-tree states.")
        .def_static(
            "create_unrestricted_residual_tree",
            [](UInt8InputArray input, RegularGridAdjacency2D adjacency, sdrt::SpatialOrder spatialOrder) {
                return wrapPythonValuedTree(MorphologicalTreeFactory::createUnrestrictedResidualTree(
                    imageFromArray(input), std::move(adjacency), sdrt::UnrestrictedResidualTreeOptions{std::move(spatialOrder)}));
            },
            "input"_a, "adjacency"_a, "spatial_order"_a = sdrt::RowMajorSpatialOrder{},
            "Create the unrestricted residual tree with an explicit shared symmetric adjacency.")
        .def_static(
            "create_saturated_residual_tree",
            [](UInt8InputArray input, PixelId infinityPixel, double radius, sdrt::SpatialOrder spatialOrder) {
                const double validatedRadius =
                    pybind_utils::requireAdjacencyRadius(radius, "MorphologicalTreeFactory.create_saturated_residual_tree");
                sdrt::SaturatedResidualTreeOptions options;
                options.spatialOrder = std::move(spatialOrder);
                return wrapPythonValuedTree(MorphologicalTreeFactory::createSaturatedResidualTree(
                    imageFromArray(input), infinityPixel, validatedRadius, std::move(options)));
            },
            "input"_a, "infinity_pixel"_a = PixelId{0}, "radius"_a = 1.5, "spatial_order"_a = sdrt::RowMajorSpatialOrder{},
            "Create the saturated residual tree with synchronized max-tree and min-tree states.")
        .def_static(
            "create_saturated_residual_tree",
            [](UInt8InputArray input, RegularGridAdjacency2D adjacency, PixelId infinityPixel, sdrt::SpatialOrder spatialOrder) {
                sdrt::SaturatedResidualTreeOptions options;
                options.spatialOrder = std::move(spatialOrder);
                return wrapPythonValuedTree(MorphologicalTreeFactory::createSaturatedResidualTree(
                    imageFromArray(input), std::move(adjacency), infinityPixel, std::move(options)));
            },
            "input"_a, "adjacency"_a, "infinity_pixel"_a = PixelId{0}, "spatial_order"_a = sdrt::RowMajorSpatialOrder{},
            "Create the saturated residual tree with an explicit shared symmetric adjacency.")
        .def_static(
            "create_tree_of_shapes",
            [](UInt8InputArray input, TopographicConvention convention) {
                return wrapPythonValuedTree(MorphologicalTreeFactory::createTreeOfShapes(imageFromArray(input), std::move(convention)));
            },
            "input"_a, "convention"_a = TopographicConvention{},
            "Create a tree of shapes using a complete topographic convention retained by the result.")
        .def_static(
            "create_from_native_topology",
            [](const std::vector<NodeId>& parent, const std::vector<NodeId>& smallestNodeMap, py::object nodeAltitudesInput, NodeId root,
               MorphologicalTreeSemantics semantics) {
                const std::vector<std::uint8_t> nodeAltitudes =
                    pythonUInt8AltitudeVector(nodeAltitudesInput, "MorphologicalTreeFactory.create_from_native_topology node_altitudes");
                return wrapPythonValuedTree(MorphologicalTreeFactory::createFromNativeTopology(
                    std::span<const NodeId>(parent), std::span<const NodeId>(smallestNodeMap), std::span<const std::uint8_t>(nodeAltitudes), root,
                    std::move(semantics)));
            },
            "parent"_a, "smallest_node_map"_a, "node_altitudes"_a, "root"_a, "semantics"_a,
            R"doc(Create a morphological tree over an abstract finite pixel domain.

No row/column interpretation is attached. Topological and support-based
algorithms remain available; image reconstruction and geometric attributes
require a separate regular 2D domain.)doc")
        .def_static(
            "create_from_native_topology",
            [](const std::vector<NodeId>& parent, const std::vector<NodeId>& smallestNodeMap, py::object nodeAltitudesInput, NodeId root, int rows, int columns,
               MorphologicalTreeSemantics semantics) {
                const std::vector<std::uint8_t> nodeAltitudes =
                    pythonUInt8AltitudeVector(nodeAltitudesInput, "MorphologicalTreeFactory.create_from_native_topology node_altitudes");
                return wrapPythonValuedTree(MorphologicalTreeFactory::createFromNativeTopology(
                    std::span<const NodeId>(parent), std::span<const NodeId>(smallestNodeMap), std::span<const std::uint8_t>(nodeAltitudes), root, rows,
                    columns, std::move(semantics)));
            },
            "parent"_a, "smallest_node_map"_a, "node_altitudes"_a, "root"_a, "rows"_a, "columns"_a, "semantics"_a,
            R"doc(Create a morphological tree from native connected-subset buffers and explicit generic capabilities.

`parent` and `node_altitudes` use the dense internal-node domain.
`smallest_node_map` uses the row-major pixel domain and must contain
`rows * columns` entries. Every committed node must have non-empty node
support, although a node may have an empty proper part.)doc")
        .def_static(
            "create_from_higra_parent",
            [](const std::vector<NodeId>& parent, py::object nodeAltitudesInput, int rows, int columns, MorphologicalTreeKind kind,
               std::optional<double> radius) {
                const std::vector<std::uint8_t> nodeAltitudes =
                    pythonUInt8AltitudeVector(nodeAltitudesInput, "MorphologicalTreeFactory.create_from_higra_parent node_altitudes");
                if (parent.size() != nodeAltitudes.size()) {
                    throw std::invalid_argument("parent and node_altitudes must have the same size");
                }

                std::optional<RegularGridAdjacency2D> adjacency;
                if (radius.has_value()) {
                    adjacency.emplace(
                        pybind_utils::makeRegularGridAdjacency2D(rows, columns, *radius, "MorphologicalTreeFactory.create_from_higra_parent"));
                }

                return wrapPythonValuedTree(MorphologicalTreeFactory::createFromHigraParent(
                    std::span<const NodeId>(parent), std::span<const std::uint8_t>(nodeAltitudes), rows, columns, kind, std::move(adjacency)));
            },
            "parent"_a, "node_altitudes"_a, "rows"_a, "columns"_a, "kind"_a, "radius"_a = py::none(),
            "Create a valued tree from an imported static Higra parent/node-altitude representation [leaves | internal nodes]. "
            "The imported Higra node-id domain is preserved until the tree is edited.");

    py::class_<HierarchySaliencyMapValidation>(m, "HierarchySaliencyMapValidation", py::module_local(false),
                                               "Utilities for validating and transforming hierarchy valuations used by saliency maps.")
        .def_static(
            "validate_hierarchy_connectivity",
            [](std::shared_ptr<PythonValuedMorphologicalTree> valuedTree, std::optional<double> radius) {
                if (!valuedTree) {
                    throw std::invalid_argument("valued tree must not be null");
                }
                const RegularGridAdjacency2D adjacency = saliencyAdjacencyPy(*valuedTree, radius, "HierarchySaliencyMapValidation.validate_hierarchy_connectivity");
                HierarchySaliencyMapValidation::validateHierarchyConnectivity(valuedTree->topology(), adjacency,
                                                                              "HierarchySaliencyMapValidation.validate_hierarchy_connectivity");
            },
            "tree"_a, "radius"_a = py::none(),
            R"doc(Validate that every hierarchy support is connected in the projection graph.

This is the graph-connectivity hypothesis required by the Cousty correspondence;
it is distinct from checking that the parent array forms one rooted tree.)doc")
        .def_static(
            "validate_hierarchy_valuation",
            [](std::shared_ptr<PythonValuedMorphologicalTree> valuedTree, py::array valuation, bool strict, bool nonnegative) {
                if (!valuedTree) {
                    throw std::invalid_argument("valued tree must not be null");
                }
                validateHierarchyValuationPy(*valuedTree, std::move(valuation), strict, nonnegative);
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
            "rank_hierarchy_valuation",
            [](std::shared_ptr<PythonValuedMorphologicalTree> valuedTree, py::array valuation, bool strict) {
                if (!valuedTree) {
                    throw std::invalid_argument("valued tree must not be null");
                }
                return rankHierarchyValuationPy(*valuedTree, std::move(valuation), strict);
            },
            "tree"_a, "valuation"_a, "strict"_a = false,
            R"doc(Convert a compatible hierarchy valuation to dense non-negative integer levels.

The input valuation must be a 1D int32, float32, or float64 array. Distinct
live-node valuation values are ranked as `0..k-1` while preserving order. Equal
values receive the same rank, so level collapse is preserved. Use `strict=True`
to reject equal parent-child valuation levels before ranking. This ranks all
live nodes; use `HierarchySaliencyMap.compute_canonical_ranked_saliency_edge_map`
when only levels that occur on graph edges should define the dense scale.)doc")
        .def_static(
            "compute_normalized_scores",
            [](std::shared_ptr<PythonValuedMorphologicalTree> valuedTree, py::array valuation, bool strict, bool nonnegative) {
                if (!valuedTree) {
                    throw std::invalid_argument("valued tree must not be null");
                }
                return computeNormalizedScoresPy(*valuedTree, std::move(valuation), strict, nonnegative);
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
            "compute_saliency_edge_map",
            [](std::shared_ptr<PythonValuedMorphologicalTree> valuedTree, py::array valuation, std::optional<double> radius, bool strict,
               HierarchyLevelConvention levelConvention, bool validateConnectivity) {
                if (!valuedTree) {
                    throw std::invalid_argument("valued tree must not be null");
                }
                return computeSaliencyEdgeMapPy(*valuedTree, std::move(valuation), radius, strict, levelConvention, validateConnectivity,
                                                "HierarchySaliencyMap.compute_saliency_edge_map");
            },
            "tree"_a, "valuation"_a, "radius"_a = py::none(), "strict"_a = false, "level_convention"_a = HierarchyLevelConvention::EdgeSaliencyValue,
            "validate_connectivity"_a = true,
            R"doc(Compute the formal edge-indexed saliency map induced by a hierarchy valuation.

This method implements `Phi(H)` for a connected hierarchy under the
quasi-flat-zone saliency convention: each adjacency edge receives
`valuation[LCA(smallest_node(source), smallest_node(target))]` when the endpoint smallest nodes differ,
and value 0 when both endpoints already belong to the same finest represented
region under `HierarchyLevelConvention.EDGE_SALIENCY_VALUE`. With
`HierarchyLevelConvention.PARTITION_APPEARANCE_LEVEL`, transition
edges receive `valuation[LCA] - 1`, exactly as in Algorithm 1 of Cousty et al.
It does not build the full `Psi(w) = Phi(QFZ(G, w))` pipeline from an
arbitrary edge-valuedTree graph. The valuation is validated before projection. Use
`strict=True` when the map must preserve every explicit parent-child level of
the tree; leave it false when equal-valued adjacent levels may be collapsed.
Formal saliency valuations must be non-negative. Graph connectivity is validated
by default; disable it only for a trusted producer-internal path.)doc")
        .def_static(
            "compute_canonical_ranked_saliency_edge_map",
            [](std::shared_ptr<PythonValuedMorphologicalTree> valuedTree, py::array valuation, std::optional<double> radius, bool strict,
               bool validateConnectivity) {
                if (!valuedTree) {
                    throw std::invalid_argument("valued tree must not be null");
                }
                return computeCanonicalRankedSaliencyEdgeMapPy(*valuedTree, std::move(valuation), radius, strict, validateConnectivity);
            },
            "tree"_a, "valuation"_a, "radius"_a = py::none(), "strict"_a = false, "validate_connectivity"_a = true,
            R"doc(Project and rank only hierarchy levels that actually occur on graph edges.

Same-smallest-node edges form the base rank when present. Unlike node-wise ranking,
unused leaf or incomparable-node values cannot introduce gaps in the canonical
edge scale. The input may contain negative values because only its order is
retained.)doc")
        .def_static(
            "compute_topological_level_edge_map",
            [](std::shared_ptr<PythonValuedMorphologicalTree> valuedTree, std::optional<double> radius) {
                if (!valuedTree) {
                    throw std::invalid_argument("valued tree must not be null");
                }
                return computeTopologicalLevelEdgeMapPy(*valuedTree, radius);
            },
            "tree"_a, "radius"_a = py::none(),
            R"doc(Compute a topological-level hierarchy saliency map on adjacency edges.

Leaf internal nodes receive level 0, and each ancestor receives one plus the
maximum level of its children.)doc")
        .def_static(
            "compute_normalized_altitude_edge_map",
            [](std::shared_ptr<PythonValuedMorphologicalTree> valuedTree, std::optional<double> radius) {
                if (!valuedTree) {
                    throw std::invalid_argument("valued tree must not be null");
                }
                return computeNormalizedAltitudeEdgeMapPy(*valuedTree, radius);
            },
            "tree"_a, "radius"_a = py::none(),
            R"doc(Compute a normalized-altitude hierarchy saliency map on adjacency edges.

Values are doubles in [0, 1], oriented so that ancestors are greater than or
equal to descendants for max-trees and min-trees.)doc");

    py::class_<ComponentTreePartitionHierarchyAdapter>(m, "ComponentTreePartitionHierarchyAdapter", py::module_local(false),
                                                       "Explicit proper-part completion of a component tree into connected graph partitions.")
        .def_static(
            "validate",
            [](std::shared_ptr<PythonValuedMorphologicalTree> valuedTree, std::optional<double> radius) {
                if (!valuedTree) {
                    throw std::invalid_argument("valued tree must not be null");
                }
                const RegularGridAdjacency2D adjacency = saliencyAdjacencyPy(*valuedTree, radius, "ComponentTreePartitionHierarchyAdapter.validate");
                ComponentTreePartitionHierarchyAdapter::validate(valuedTree->topology(), adjacency);
            },
            "tree"_a, "radius"_a = py::none(), "Validate the connected proper-part completion in the selected graph.")
        .def_static(
            "compute_partition_appearance_levels",
            [](std::shared_ptr<PythonValuedMorphologicalTree> valuedTree) {
                if (!valuedTree) {
                    throw std::invalid_argument("valued tree must not be null");
                }
                std::vector<int> levels = ComponentTreePartitionHierarchyAdapter::computePartitionAppearanceLevels(valuedTree->topology());
                const int size = static_cast<int>(levels.size());
                return pybind_utils::toNumpyOwned(std::move(levels), size);
            },
            "tree"_a,
            R"doc(Return positive partition-appearance levels for the component-tree completion.

Use these levels with `HierarchyLevelConvention.PARTITION_APPEARANCE_LEVEL` to
apply the exact `level(LCA)-1` convention.)doc")
        .def_static(
            "compute_saliency_edge_map",
            [](std::shared_ptr<PythonValuedMorphologicalTree> valuedTree, py::array partitionAppearanceLevels, std::optional<double> radius,
               bool strict) {
                if (!valuedTree) {
                    throw std::invalid_argument("valued tree must not be null");
                }
                return computeSaliencyEdgeMapPy(*valuedTree, std::move(partitionAppearanceLevels), radius, strict,
                                                HierarchyLevelConvention::PartitionAppearanceLevel, true,
                                                "ComponentTreePartitionHierarchyAdapter.compute_saliency_edge_map");
            },
            "tree"_a, "partition_appearance_levels"_a, "radius"_a = py::none(), "strict"_a = false,
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
            "compute_extinction_values",
            [](std::shared_ptr<PythonValuedMorphologicalTree> valuedTree, py::array attribute, ShapeSpaceExtremaPolarity polarity) {
                if (!valuedTree) {
                    throw std::invalid_argument("valued tree must not be null");
                }
                return computeShapeSpaceExtinctionValuesPy(*valuedTree, std::move(attribute), polarity);
            },
            "tree"_a, "attribute"_a, "polarity"_a,
            R"doc(Compute shaping extinction values on the original tree-node graph.

`attribute` must be a dense float32 or float64 array with one value per
internal NodeId slot. The values need not be monotone on the original tree.
Parent-child relations form the graph for a second component-tree computation;
`polarity` selects its local minima or maxima. The returned dictionary contains
an `extrema` list and the dense extinction array `node_scores`, preserving the
input dtype.)doc")
        .def_static(
            "project_contour_scores",
            [](std::shared_ptr<PythonValuedMorphologicalTree> valuedTree, py::array nodeScores, std::optional<double> radius) {
                if (!valuedTree) {
                    throw std::invalid_argument("valued tree must not be null");
                }
                return projectShapeSpaceContourScoresPy(*valuedTree, std::move(nodeScores), radius);
            },
            "tree"_a, "node_scores"_a, "radius"_a = py::none(),
            R"doc(Project dense original-node scores onto their full region contours.

Every image adjacency edge receives the maximum score among original regions
whose contour contains that edge. The return value is the standard edge-map
dictionary. If `radius` is omitted, the tree's stored construction adjacency is
used. This operation is distinct from projecting an LCA valuation.)doc")
        .def_static(
            "compute",
            [](std::shared_ptr<PythonValuedMorphologicalTree> valuedTree, py::array attribute, ShapeSpaceExtremaPolarity polarity,
               std::optional<double> radius) {
                if (!valuedTree) {
                    throw std::invalid_argument("valued tree must not be null");
                }
                return computeShapeSpaceSaliencyPy(*valuedTree, std::move(attribute), polarity, radius);
            },
            "tree"_a, "attribute"_a, "polarity"_a, "radius"_a = py::none(),
            R"doc(Compute Xu shaping extinctions and their maximum-on-contours map.

The returned dictionary contains `extrema`, dense extinction `node_scores`, and
the standard edge-map dictionary `edge_map`. If `radius` is omitted, the tree's
stored construction adjacency is used.)doc");

    py::class_<HierarchySaliencyMapProjection>(m, "HierarchySaliencyMapProjection", py::module_local(false),
                                               "Derived projections and contour materializations of hierarchy edge saliency maps.")
        .def_static(
            "edge_map_to_pixel_image", [](py::dict edgeMap, EdgeToPixelReducer reducer) { return edgeMapToPixelImagePy(std::move(edgeMap), reducer); }, "edge_map"_a,
            "reducer"_a = EdgeToPixelReducer::Max,
            R"doc(Rasterize an edge-indexed map into a pixel image for display.

This helper is a visualization projection only. Each edge value contributes to
both endpoint pixels. `EdgeToPixelReducer.MAX` writes the maximum incident edge
value per pixel; `EdgeToPixelReducer.MEAN` writes the arithmetic mean of
incident edge values. The formal saliency representation remains the edge map
dictionary.)doc")
        .def_static(
            "threshold_cut", [](py::dict edgeMap, double threshold) { return thresholdCutPy(std::move(edgeMap), threshold); }, "edge_map"_a, "threshold"_a,
            R"doc(Threshold an edge saliency map into a contour edge set.

The input is the dictionary returned by `HierarchySaliencyMap`. The returned
dictionary contains `num_rows`, `num_columns`, `adjacency_radius`, `sources`, and
`targets`, selecting edges whose saliency is greater than or equal to
`threshold`.)doc")
        .def_static(
            "node_contour_edges",
            [](std::shared_ptr<PythonValuedMorphologicalTree> valuedTree, std::optional<double> radius) {
                if (!valuedTree) {
                    throw std::invalid_argument("valued tree must not be null");
                }
                return nodeContourEdgesPy(*valuedTree, radius);
            },
            "tree"_a, "radius"_a = py::none(),
            R"doc(Project transition adjacency edges onto their hierarchy node.

Each returned edge carries the node id `LCA(smallest_node(source), smallest_node(target))` in the
`nodes` array. Edges whose endpoints have the same smallest node are omitted and have
implicit value 0 in the full formal saliency map. If `radius` is not provided,
the tree's stored construction adjacency is used.)doc")
        .def_static(
            "compute_incremental_node_contours",
            [](std::shared_ptr<PythonValuedMorphologicalTree> valuedTree, std::optional<double> radius) {
                if (!valuedTree) {
                    throw std::invalid_argument("valued tree must not be null");
                }
                return computeIncrementalNodeContoursPy(*valuedTree, radius);
            },
            "tree"_a, "radius"_a = py::none(),
            R"doc(Compute per-node incremental contour edges.

The returned dictionary contains `offsets`, `sources`, and `targets`.
For node `u`, the slice `offsets[u]:offsets[u + 1]` stores all transition
adjacency edges whose endpoint smallest nodes differ and whose endpoint-smallest-node LCA is
`u`.)doc")
        .def_static(
            "project_node_valuation",
            [](py::dict contours, py::array nodeValuation) {
                return projectNodeValuationOnIncrementalContoursPy(std::move(contours), std::move(nodeValuation));
            },
            "contours"_a, "node_valuation"_a,
            R"doc(Project a dense node valuation onto per-node transition contour edges.

`contours` must be the dictionary returned by `compute_incremental_node_contours`.
`node_valuation` must be a 1D float32 or float64 array with one value per dense
internal node slot. The result is sparse over transition edges; same-smallest-node graph
edges are omitted and have implicit value 0 in the full formal saliency map.)doc")
        .def_static(
            "threshold_by_node_valuation",
            [](py::dict contours, py::array nodeValuation, double threshold) {
                return thresholdIncrementalContoursByNodeValuationPy(std::move(contours), std::move(nodeValuation), threshold);
            },
            "contours"_a, "node_valuation"_a, "threshold"_a,
            R"doc(Threshold per-node transition contours by dense node valuation.

All edges in a node slice are selected when `node_valuation[node] >= threshold`.)doc");

    valuedTreeCls
        .def(
            "set_node_altitude",
            [](PythonValuedMorphologicalTree& tree, NodeId nodeId, py::object altitudeInput) {
                if (!tree.topology().isNode(nodeId)) {
                    throw std::invalid_argument("invalid NodeId for altitude update");
                }
                tree.visit([&](auto& native) {
                    using Altitude = typename std::remove_reference_t<decltype(*native)>::AltitudeType;
                    native->setNodeAltitude(nodeId, pythonUnsignedAltitudeValue<Altitude>(altitudeInput, nodeId, "ValuedMorphologicalTree.set_node_altitude"));
                });
            },
            "node_id"_a, "altitude"_a, "Set one node altitude while preserving the declared strict parent-child altitude order.")
        .def_property(
            "node_altitudes",
            [](const PythonValuedMorphologicalTree& tree) -> py::array {
                return tree.visit([](const auto& native) -> py::array {
                    using Altitude = typename std::remove_reference_t<decltype(*native)>::AltitudeType;
                    const auto altitude = native->nodeAltitudes();
                    return pybind_utils::toNumpyOwned(std::vector<Altitude>(altitude.begin(), altitude.end()), static_cast<int>(altitude.size()));
                });
            },
            [](PythonValuedMorphologicalTree& tree, py::object altitudeInput) {
                tree.visit([&](auto& native) {
                    using Altitude = typename std::remove_reference_t<decltype(*native)>::AltitudeType;
                    std::vector<Altitude> altitude = pythonUnsignedAltitudeVector<Altitude>(altitudeInput, "ValuedMorphologicalTree.node_altitudes");
                    pybind_utils::requireVectorSize(altitude, static_cast<std::size_t>(native->topology().numInternalNodeSlots()), "altitude");
                    native->setNodeAltitudes(std::move(altitude));
                });
            },
            "Dense altitude buffer indexed by internal NodeId. Assignments validate the declared strict altitude order.")
        .def(
            "node_altitude",
            [](const PythonValuedMorphologicalTree& tree, NodeId nodeId) {
                return tree.visit([&](const auto& native) -> std::uint64_t { return static_cast<std::uint64_t>(native->nodeAltitude(nodeId)); });
            },
            "node_id"_a, "Return the altitude associated with node_id.")
        .def(
            "node_residue",
            [](const PythonValuedMorphologicalTree& tree, NodeId nodeId) {
                return tree.visit([&](const auto& native) -> std::int64_t { return static_cast<std::int64_t>(native->nodeResidue(nodeId)); });
            },
            "node_id"_a, "Return the residue between node_id and its parent.")
        .def(
            "reconstruct_from_node_altitudes",
            [](const PythonValuedMorphologicalTree& tree) -> py::array {
                return tree.visit([](const auto& native) -> py::array { return reconstructFromNodeAltitudesOf(*native); });
            },
            "Reconstruct the current tree from node altitudes using the concrete altitude dtype.")
        .def("reconstruct_from_node_contributions", &reconstructFromNodeContributions, "node_contributions"_a,
             "Reconstruct floating-point node contributions by zero-baseline ancestor summation.")
        .def(
            "export_higra_hierarchy",
            [](const PythonValuedMorphologicalTree& tree) -> py::object {
                return tree.visit([](const auto& native) -> py::object { return py::cast(exportHigraHierarchyOf(*native)); });
            },
            "Export the current rooted tree to a compact Higra parent/altitude representation.")
        .def("validate_node_altitude_buffer_shape",
             [](const PythonValuedMorphologicalTree& tree) { tree.visit([](const auto& native) { native->validateNodeAltitudeBufferShape(); }); },
             "Validate that the dense altitude buffer covers every internal node slot.")
        .def("validate_monotone_node_altitudes",
             [](const PythonValuedMorphologicalTree& tree) { tree.visit([](const auto& native) { native->validateMonotoneNodeAltitudes(); }); },
             "Validate the declared strict parent-child altitude order where applicable.")
        .def("project_node_values_to_exported_higra",
             [](const PythonValuedMorphologicalTree& tree, const py::array& nodeValues, py::object attributes) {
                 return projectNodeValuesToExportedHigraOf(tree, nodeValues, std::move(attributes));
             },
             "node_values"_a, "attributes"_a,
             "Project a node-indexed scalar or 2D attribute buffer to the compact Higra layout produced by export_higra_hierarchy().");

    bindTreeQueryApi<PythonValuedMorphologicalTree>(valuedTreeCls);
}

} // namespace mmcfilters::pybindings
