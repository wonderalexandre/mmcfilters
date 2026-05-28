#pragma once

#include "../AttributeNames.hpp"
#include "../AttributeResultTypes.hpp"
#include "../computers/AreaComputer.hpp"
#include "../computers/BitquadAttributeComputer.hpp"
#include "../computers/BoundingBoxComputer.hpp"
#include "../computers/ContourSideAttributeComputer.hpp"
#include "../computers/GrayLevelStatsComputer.hpp"
#include "../computers/MaxDistComputer.hpp"
#include "../computers/MomentBasedAttributeComputer.hpp"
#include "../computers/TreeTopologyComputer.hpp"
#include "../computers/VolumeComputer.hpp"
#include "../../trees/MorphologicalTree.hpp"
#include "../../trees/detail/HigraExportLayoutDetail.hpp"
#include "../../utils/Altitude.hpp"
#include "../../utils/Common.hpp"
#include "../../utils/Image.hpp"

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcfilters::detail {

inline bool isGenericAltitudeUnitAttribute(Attribute attribute) noexcept {
    switch (attribute) {
        case VOLUME:
        case RELATIVE_VOLUME:
        case LEVEL:
        case MEAN_LEVEL:
        case VARIANCE_LEVEL:
        case GRAY_HEIGHT:
        case MAX_DIST:
            return true;
        default:
            return false;
    }
}

inline std::vector<NodeId> makeRowMajorProperPartOrder(const MorphologicalTree& tree) {
    std::vector<NodeId> properParts;
    properParts.reserve(static_cast<size_t>(tree.getNumTotalProperParts()));
    for (NodeId properPart = 0; properPart < tree.getNumTotalProperParts(); ++properPart) {
        properParts.push_back(properPart);
    }
    return properParts;
}

template<class Computer, std::floating_point Real>
inline void computeUnitAttributeWithComputer(const MorphologicalTree& tree, std::span<const NodeId> unitProperParts, std::span<Real> unitValues, const AttributeNames& attrNames, Attribute attribute) {
    const std::array<Attribute, 1> requested{attribute};
    Computer::computeUnitRows(UnitAttributeComputeContext<Real>{
        tree,
        unitProperParts,
        unitValues,
        attrNames,
        std::span<const Attribute>(requested)});
}

template <std::floating_point Real>
inline void computeTopologyUnitAttribute(const MorphologicalTree& tree, std::span<const NodeId> unitProperParts, std::span<Real> unitValues, const AttributeNames& attrNames, Attribute attribute) {
    using namespace ::mmcfilters::attributes::computers;

    switch (attribute) {
        case AREA:
            computeUnitAttributeWithComputer<AreaComputer>(tree, unitProperParts, unitValues, attrNames, attribute);
            return;

        case BOX_WIDTH:
        case BOX_HEIGHT:
        case DIAGONAL_LENGTH:
        case RECTANGULARITY:
        case RATIO_WH:
        case BOX_COL_MIN:
        case BOX_COL_MAX:
        case BOX_ROW_MIN:
        case BOX_ROW_MAX:
            computeUnitAttributeWithComputer<BoundingBoxComputer>(tree, unitProperParts, unitValues, attrNames, attribute);
            return;

        case CENTRAL_MOMENT_20:
        case CENTRAL_MOMENT_02:
        case CENTRAL_MOMENT_11:
        case CENTRAL_MOMENT_30:
        case CENTRAL_MOMENT_03:
        case CENTRAL_MOMENT_21:
        case CENTRAL_MOMENT_12:
            computeUnitAttributeWithComputer<CentralMomentsComputer>(tree, unitProperParts, unitValues, attrNames, attribute);
            return;

        case HU_MOMENT_1:
        case HU_MOMENT_2:
        case HU_MOMENT_3:
        case HU_MOMENT_4:
        case HU_MOMENT_5:
        case HU_MOMENT_6:
        case HU_MOMENT_7:
            computeUnitAttributeWithComputer<HuMomentsComputer>(tree, unitProperParts, unitValues, attrNames, attribute);
            return;

        case INERTIA:
        case COMPACTNESS:
        case ECCENTRICITY:
        case LENGTH_MAJOR_AXIS:
        case LENGTH_MINOR_AXIS:
        case AXIS_ORIENTATION:
        case CIRCULARITY:
            computeUnitAttributeWithComputer<MomentBasedAttributeComputer>(tree, unitProperParts, unitValues, attrNames, attribute);
            return;

        case BITQUADS_AREA:
        case BITQUADS_NUMBER_EULER:
        case BITQUADS_NUMBER_HOLES:
        case BITQUADS_PERIMETER:
        case BITQUADS_PERIMETER_CONTINUOUS:
        case BITQUADS_CIRCULARITY:
        case BITQUADS_PERIMETER_AVERAGE:
        case BITQUADS_LENGTH_AVERAGE:
        case BITQUADS_WIDTH_AVERAGE:
            computeUnitAttributeWithComputer<BitquadAttributeComputer>(tree, unitProperParts, unitValues, attrNames, attribute);
            return;

        case HEIGHT_NODE:
        case DEPTH_NODE:
        case IS_LEAF_NODE:
        case IS_ROOT_NODE:
        case NUM_CHILDREN_NODE:
        case NUM_SIBLINGS_NODE:
        case NUM_DESCENDANTS_NODE:
        case NUM_LEAF_DESCENDANTS_NODE:
        case LEAF_RATIO_NODE:
        case BALANCE_NODE:
        case AVG_CHILD_HEIGHT_NODE:
            computeUnitAttributeWithComputer<TreeTopologyComputer>(tree, unitProperParts, unitValues, attrNames, attribute);
            return;

        case CONTOUR_PIXELS:
        case CONTOUR_PERIMETER:
        case CONTOUR_SIDE_NORTH:
        case CONTOUR_SIDE_WEST:
        case CONTOUR_SIDE_EAST:
        case CONTOUR_SIDE_SOUTH:
            computeUnitAttributeWithComputer<ContourSideAttributeComputer>(tree, unitProperParts, unitValues, attrNames, attribute);
            return;

        default:
            throw std::runtime_error("No exported-Higra unit projection is registered for the requested attribute.");
    }
}

template <std::floating_point Real = float>
inline std::vector<Real> computeTopologyUnitAttributeRows(const MorphologicalTree& tree, std::span<const NodeId> unitProperParts, const AttributeNames& attrNames) {
    const size_t numColumns = static_cast<size_t>(attrNames.NUM_ATTRIBUTES);
    std::vector<Real> unitValues(unitProperParts.size() * numColumns, std::numeric_limits<Real>::quiet_NaN());

    for (const auto& [attribute, _] : attrNames.indexMap) {
        computeTopologyUnitAttribute<Real>(tree, unitProperParts, unitValues, attrNames, attribute);
    }

    return unitValues;
}

template<std::floating_point Real, AltitudeValue T>
inline void computeGenericAltitudeUnitAttribute(const MorphologicalTree& tree, std::span<const T> altitude, std::span<const NodeId> unitProperParts, std::span<Real> unitValues, const AttributeNames& attrNames, Attribute attribute) {
    if (!isGenericAltitudeUnitAttribute(attribute)) {
        return;
    }

    const std::array<Attribute, 1> requested{attribute};
    const AltitudeUnitAttributeComputeContext<Real, T> context{
        tree,
        altitude,
        unitProperParts,
        unitValues,
        attrNames,
        std::span<const Attribute>(requested)};

    using namespace ::mmcfilters::attributes::computers;
    switch (attribute) {
        case VOLUME:
        case RELATIVE_VOLUME:
            VolumeComputer::computeUnitRows(context);
            return;

        case LEVEL:
        case MEAN_LEVEL:
        case VARIANCE_LEVEL:
        case GRAY_HEIGHT:
            GrayLevelStatsComputer::computeUnitRows(context);
            return;

        case MAX_DIST:
            MaxDistComputer::computeUnitRows(context);
            return;

        default:
            return;
    }
}

template<std::floating_point Real = float, AltitudeValue T>
inline std::vector<Real> computeUnitAttributeRows(const MorphologicalTree& tree, std::span<const T> altitude, std::span<const NodeId> unitProperParts, const AttributeNames& attrNames) {
    const size_t numColumns = static_cast<size_t>(attrNames.NUM_ATTRIBUTES);
    std::vector<Real> unitValues(unitProperParts.size() * numColumns, std::numeric_limits<Real>::quiet_NaN());

    for (const auto& [attribute, _] : attrNames.indexMap) {
        if (isGenericAltitudeUnitAttribute(attribute)) {
            computeGenericAltitudeUnitAttribute<Real>(tree, altitude, unitProperParts, unitValues, attrNames, attribute);
            continue;
        }

        computeTopologyUnitAttribute<Real>(tree, unitProperParts, unitValues, attrNames, attribute);
    }

    return unitValues;
}

template <std::floating_point Real>
inline void copyScalarUnitRowsToPreservedHigra(
    const MorphologicalTree& tree,
    const AttributeNames& attrNames,
    std::span<const Real> unitValues,
    std::span<Real> projected)
{
    const int numProperParts = tree.getNumTotalProperParts();
    const int numColumns = attrNames.NUM_ATTRIBUTES;
    const size_t expectedUnitSize = static_cast<size_t>(numProperParts) * static_cast<size_t>(numColumns);
    if (unitValues.size() != expectedUnitSize || projected.size() < expectedUnitSize) {
        throw std::logic_error("Preserved Higra unit-row projection received incompatible buffer shapes.");
    }

    for (NodeId leafIndex = 0; leafIndex < numProperParts; ++leafIndex) {
        for (int column = 0; column < numColumns; ++column) {
            projected[static_cast<size_t>(leafIndex) * static_cast<size_t>(numColumns) + static_cast<size_t>(column)] =
                unitValues[static_cast<size_t>(leafIndex) * static_cast<size_t>(numColumns) + static_cast<size_t>(column)];
        }
    }
}

template <std::floating_point Real>
inline void copyScalarInternalRowsToPreservedHigra(
    const MorphologicalTree& tree,
    const AttributeNames& attrNames,
    std::span<const Real> nodeValues,
    std::span<Real> projected,
    int outputSize)
{
    const int numColumns = attrNames.NUM_ATTRIBUTES;
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        const NodeId outputNodeId = tree.getHigraNodeId(nodeId);
        if (outputNodeId == InvalidNode || outputNodeId < 0 || outputNodeId >= outputSize) {
            throw std::runtime_error("Cannot project attributes to Higra node-id space: a live node has no valid Higra id.");
        }
        for (int column = 0; column < numColumns; ++column) {
            projected[static_cast<size_t>(outputNodeId) * static_cast<size_t>(numColumns) + static_cast<size_t>(column)] =
                nodeValues[static_cast<size_t>(nodeId) * static_cast<size_t>(numColumns) + static_cast<size_t>(column)];
        }
    }
}

inline AttributeNames deltaUnitAttributeNames(const AttributeNamesWithDelta& deltaNames) {
    std::vector<Attribute> attributes;
    attributes.reserve(deltaNames.indexMap.size());
    for (const auto& [attributeKey, _] : deltaNames.indexMap) {
        if (std::find(attributes.begin(), attributes.end(), attributeKey.attr) == attributes.end()) {
            attributes.push_back(attributeKey.attr);
        }
    }
    return AttributeNames::fromList(attributes);
}

template <std::floating_point Real>
inline void copyDeltaUnitRowsToPreservedHigra(
    const MorphologicalTree& tree,
    const AttributeNamesWithDelta& deltaNames,
    const AttributeNames& unitNames,
    std::span<const Real> unitValues,
    std::span<Real> projected)
{
    const int numProperParts = tree.getNumTotalProperParts();
    const int numColumns = deltaNames.NUM_ATTRIBUTES;
    const size_t expectedUnitSize = static_cast<size_t>(numProperParts) * static_cast<size_t>(unitNames.NUM_ATTRIBUTES);
    const size_t expectedProjectedUnitSize = static_cast<size_t>(numProperParts) * static_cast<size_t>(numColumns);
    if (unitValues.size() != expectedUnitSize || projected.size() < expectedProjectedUnitSize) {
        throw std::logic_error("Preserved Higra delta unit-row projection received incompatible buffer shapes.");
    }

    // Proper parts have no ancestor/descendant neighbourhood in this domain, so
    // every delta column receives the unit value for its scalar attribute.
    for (NodeId leafIndex = 0; leafIndex < numProperParts; ++leafIndex) {
        for (const auto& [attributeKey, _] : deltaNames.indexMap) {
            projected[static_cast<size_t>(deltaNames.linearIndex(leafIndex, attributeKey))] =
                unitValues[static_cast<size_t>(unitNames.linearIndex(leafIndex, attributeKey.attr))];
        }
    }
}

template <std::floating_point Real>
inline void copyDeltaInternalRowsToPreservedHigra(
    const MorphologicalTree& tree,
    const AttributeNamesWithDelta& deltaNames,
    std::span<const Real> nodeValues,
    std::span<Real> projected,
    int outputSize)
{
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        const NodeId outputNodeId = tree.getHigraNodeId(nodeId);
        if (outputNodeId == InvalidNode || outputNodeId < 0 || outputNodeId >= outputSize) {
            throw std::runtime_error("Cannot project delta attributes to Higra node-id space: a live node has no valid Higra id.");
        }
        for (const auto& [attributeKey, _] : deltaNames.indexMap) {
            projected[static_cast<size_t>(deltaNames.linearIndex(outputNodeId, attributeKey.attr, attributeKey.delta))] =
                nodeValues[static_cast<size_t>(deltaNames.linearIndex(nodeId, attributeKey.attr, attributeKey.delta))];
        }
    }
}

/**
 * @brief Projects a scalar result from the internal node-id space to another
 * public node-id space.
 *
 * @details
 * The internal pipeline always computes in `MORPHOLOGICAL_TREE` space. This
 * projection step is only applied at the API boundary, for example to expose
 * results in the preserved imported Higra convention. Proper-part/leaf ids in
 * that target space receive the same unit-component values used by compact
 * Higra export; every live internal node must have a valid target id.
 */
template <std::floating_point Real>
inline ComputedAttributeData<Real> projectComputedDataToNodeIdSpace(const MorphologicalTree& tree, ComputedAttributeData<Real> computed, NodeIdSpace outputSpace){
    if (outputSpace == NodeIdSpace::MORPHOLOGICAL_TREE) {
        computed.nodeIdSpace = NodeIdSpace::MORPHOLOGICAL_TREE;
        return computed;
    }

    const int numAttributes = computed.first.NUM_ATTRIBUTES;
    const int outputSize = tree.getNodeIdSpaceSize(outputSpace);
    std::vector<Real> projected(static_cast<size_t>(outputSize) * static_cast<size_t>(numAttributes), std::numeric_limits<Real>::quiet_NaN());

    if (outputSpace != NodeIdSpace::HIGRA) {
        throw std::runtime_error("Unsupported node-id projection target.");
    }

    const std::vector<NodeId> properParts = makeRowMajorProperPartOrder(tree);
    const std::vector<Real> unitValues = computeTopologyUnitAttributeRows<Real>(tree, properParts, computed.first);
    copyScalarUnitRowsToPreservedHigra<Real>(tree, computed.first, unitValues, projected);
    copyScalarInternalRowsToPreservedHigra<Real>(tree, computed.first, computed.second, projected, outputSize);

    return {std::move(computed.first), std::move(projected), outputSpace};
}

template<std::floating_point Real, AltitudeValue T>
inline ComputedAttributeData<Real> projectComputedDataToNodeIdSpace(const MorphologicalTree& tree, std::span<const T> altitude, ComputedAttributeData<Real> computed, NodeIdSpace outputSpace){
    if (outputSpace == NodeIdSpace::MORPHOLOGICAL_TREE) {
        computed.nodeIdSpace = NodeIdSpace::MORPHOLOGICAL_TREE;
        return computed;
    }

    const int numAttributes = computed.first.NUM_ATTRIBUTES;
    const int outputSize = tree.getNodeIdSpaceSize(outputSpace);
    std::vector<Real> projected(static_cast<size_t>(outputSize) * static_cast<size_t>(numAttributes), std::numeric_limits<Real>::quiet_NaN());

    if (outputSpace != NodeIdSpace::HIGRA) {
        throw std::runtime_error("Unsupported node-id projection target.");
    }

    const std::vector<NodeId> properParts = makeRowMajorProperPartOrder(tree);
    const std::vector<Real> unitValues = computeUnitAttributeRows<Real>(tree, altitude, properParts, computed.first);
    copyScalarUnitRowsToPreservedHigra<Real>(tree, computed.first, unitValues, projected);
    copyScalarInternalRowsToPreservedHigra<Real>(tree, computed.first, computed.second, projected, outputSize);

    return {std::move(computed.first), std::move(projected), outputSpace};
}

/**
 * @brief Delta-aware counterpart of the scalar node-id-space projection.
 */
template <std::floating_point Real>
inline ComputedAttributeDataWithDelta<Real> projectComputedDataToNodeIdSpace(const MorphologicalTree& tree, ComputedAttributeDataWithDelta<Real> computed, NodeIdSpace outputSpace) {
    if (outputSpace == NodeIdSpace::MORPHOLOGICAL_TREE) {
        computed.nodeIdSpace = NodeIdSpace::MORPHOLOGICAL_TREE;
        return computed;
    }

    const int numAttributes = computed.first.NUM_ATTRIBUTES;
    const int outputSize = tree.getNodeIdSpaceSize(outputSpace);
    std::vector<Real> projected(
        static_cast<size_t>(outputSize) * static_cast<size_t>(numAttributes),
        std::numeric_limits<Real>::quiet_NaN());

    if (outputSpace != NodeIdSpace::HIGRA) {
        throw std::runtime_error("Unsupported node-id projection target.");
    }

    const std::vector<NodeId> properParts = makeRowMajorProperPartOrder(tree);
    const AttributeNames unitNames = deltaUnitAttributeNames(computed.first);
    const std::vector<Real> unitValues = computeTopologyUnitAttributeRows<Real>(tree, properParts, unitNames);
    copyDeltaUnitRowsToPreservedHigra<Real>(tree, computed.first, unitNames, unitValues, projected);
    copyDeltaInternalRowsToPreservedHigra<Real>(tree, computed.first, computed.second, projected, outputSize);

    return {std::move(computed.first), std::move(projected), outputSpace};
}

template<std::floating_point Real, AltitudeValue T>
inline ComputedAttributeDataWithDelta<Real> projectComputedDataToNodeIdSpace(const MorphologicalTree& tree, std::span<const T> altitude, ComputedAttributeDataWithDelta<Real> computed, NodeIdSpace outputSpace) {
    if (outputSpace == NodeIdSpace::MORPHOLOGICAL_TREE) {
        computed.nodeIdSpace = NodeIdSpace::MORPHOLOGICAL_TREE;
        return computed;
    }

    const int numAttributes = computed.first.NUM_ATTRIBUTES;
    const int outputSize = tree.getNodeIdSpaceSize(outputSpace);
    std::vector<Real> projected(
        static_cast<size_t>(outputSize) * static_cast<size_t>(numAttributes),
        std::numeric_limits<Real>::quiet_NaN());

    if (outputSpace != NodeIdSpace::HIGRA) {
        throw std::runtime_error("Unsupported node-id projection target.");
    }

    const std::vector<NodeId> properParts = makeRowMajorProperPartOrder(tree);
    const AttributeNames unitNames = deltaUnitAttributeNames(computed.first);
    const std::vector<Real> unitValues = computeUnitAttributeRows<Real>(tree, altitude, properParts, unitNames);
    copyDeltaUnitRowsToPreservedHigra<Real>(tree, computed.first, unitNames, unitValues, projected);
    copyDeltaInternalRowsToPreservedHigra<Real>(tree, computed.first, computed.second, projected, outputSize);

    return {std::move(computed.first), std::move(projected), outputSpace};
}

template<std::floating_point Real, AltitudeValue T>
inline std::vector<Real> projectNodeValuesToExportedHigraTyped(const MorphologicalTree& topology, std::span<const T> altitude, const AttributeNames& attrNames, std::span<const Real> nodeValues) {
    const int numColumns = attrNames.NUM_ATTRIBUTES;
    const size_t expectedSize = static_cast<size_t>(topology.getNumInternalNodeSlots()) * static_cast<size_t>(numColumns);
    if (nodeValues.size() != expectedSize) {
        throw std::invalid_argument("Node-value buffer size must match the dense internal-node domain and requested attributes.");
    }

    const auto layout = computeExportedHigraLayout(topology, altitude);
    const size_t numColumnValues = static_cast<size_t>(numColumns);
    const std::vector<Real> unitValues = computeUnitAttributeRows<Real>(topology, altitude, layout.properParts, attrNames);

    // Compact Higra buffers place unit proper-part rows before internal-node rows.
    std::vector<Real> projected(static_cast<size_t>(layout.numVertices) * numColumnValues, Real{0});

    for (NodeId leafIndex = 0; leafIndex < layout.numLeaves; ++leafIndex) {
        for (int column = 0; column < numColumns; ++column) {
            const size_t columnIndex = static_cast<size_t>(column);
            projected[static_cast<size_t>(leafIndex) * numColumnValues + columnIndex] = unitValues[static_cast<size_t>(leafIndex) * numColumnValues + columnIndex];
        }
    }

    for (NodeId nodeId : layout.sortedNodes) {
        const NodeId higraNodeId = layout.nodeToHigra[static_cast<size_t>(nodeId)];
        for (int column = 0; column < numColumns; ++column) {
            const size_t columnIndex = static_cast<size_t>(column);
            projected[static_cast<size_t>(higraNodeId) * numColumnValues + columnIndex] = nodeValues[static_cast<size_t>(nodeId) * numColumnValues + columnIndex];
        }
    }

    return projected;
}

/**
 * @brief Projects one scalar node attribute back to the original image domain.
 *
 * Each proper part receives the value stored at its owner component. The input
 * buffer must use the canonical dense internal-node layout described by
 * `attrNames`.
 */
template <std::floating_point Real>
inline ImagePtr<Real> mapNodeAttributeToImage(const MorphologicalTree& tree, const AttributeNames& attrNames, std::span<const Real> nodeValues, Attribute attribute) {
    ImagePtr<Real> imgPtr = Image<Real>::create(tree.getNumRowsOfImage(), tree.getNumColsOfImage());
    Real* img = imgPtr->rawData();
    for (int p = 0; p < imgPtr->getSize(); ++p) {
        const NodeId nodeId = tree.getProperPartOwner(p);
        img[p] = nodeValues[attrNames.linearIndex(nodeId, attribute)];
    }
    return imgPtr;
}

} // namespace mmcfilters::detail
