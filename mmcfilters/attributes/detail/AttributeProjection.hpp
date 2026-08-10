#pragma once

#include "AttributeFamilyScheduler.hpp"
#include "../AttributeNames.hpp"
#include "../AttributeResultTypes.hpp"
#include "../../trees/MorphologicalTree.hpp"
#include "../../trees/detail/CommittedTreeAccess.hpp"
#include "../../trees/detail/HigraExportLayoutDetail.hpp"
#include "../../utils/Altitude.hpp"
#include "../../utils/Common.hpp"
#include "../../utils/Contract.hpp"
#include "../../utils/Image.hpp"

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <span>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace mmcfilters::detail {

/**
 * @brief Tests whether generic altitude unit attribute holds.
 *
 * @param attribute Attribute requested by the operation.
 * @return True when generic altitude unit attribute; otherwise false.
 */
inline bool isGenericAltitudeUnitAttribute(Attribute attribute) noexcept {
    return attributeHasComputerDomain<attributes::computers::AttributeComputerDomain::Altitude>(attribute);
}

/**
 * @brief Creates the row-major ordering of direct proper-part pixels.
 *
 * @param tree Tree topology used by the operation.
 * @return Direct proper-part pixel identifiers in row-major order.
 */
inline std::vector<NodeId> makeRowMajorProperPartOrder(const MorphologicalTree& tree) {
    std::vector<NodeId> properParts;
    properParts.reserve(static_cast<size_t>(tree.getNumTotalProperParts()));
    for (NodeId properPart = 0; properPart < tree.getNumTotalProperParts(); ++properPart) {
        properParts.push_back(properPart);
    }
    return properParts;
}

/**
 * @brief Computes unit attribute with computer.
 *
 * @param tree Tree topology used by the operation.
 * @param unitProperParts Proper parts treated as unit supports.
 * @param unitValues Unit-support values read or written by the operation.
 * @param attrNames Layout mapping attributes to buffer columns.
 * @param attribute Attribute requested by the operation.
 */
template <class Computer, std::floating_point Real>
inline void computeUnitAttributeWithComputer(const MorphologicalTree& tree, std::span<const NodeId> unitProperParts, std::span<Real> unitValues,
                                             const AttributeNames& attrNames, Attribute attribute) {
    const std::array<Attribute, 1> requested{attribute};
    Computer::computeUnitRows(UnitAttributeComputeContext<Real>{tree, unitProperParts, unitValues, attrNames, std::span<const Attribute>(requested)});
}

/**
 * @brief Attempts to compute topology unit attribute with computer.
 *
 * @param tree Tree topology used by the operation.
 * @param unitProperParts Proper parts treated as unit supports.
 * @param unitValues Unit-support values read or written by the operation.
 * @param attrNames Layout mapping attributes to buffer columns.
 * @param attribute Attribute requested by the operation.
 * @return True when the documented condition holds; otherwise false.
 */
template <class Computer, std::floating_point Real>
inline bool tryComputeTopologyUnitAttributeWithComputer(const MorphologicalTree& tree, std::span<const NodeId> unitProperParts, std::span<Real> unitValues,
                                                        const AttributeNames& attrNames, Attribute attribute) {
    if constexpr (Computer::domain == attributes::computers::AttributeComputerDomain::Topology) {
        if (attributes::computers::producesAttribute<Computer>(attribute)) {
            computeUnitAttributeWithComputer<Computer>(tree, unitProperParts, unitValues, attrNames, attribute);
            return true;
        }
    }
    return false;
}

/**
 * @brief Computes topology unit attribute with computers.
 *
 * @param tree Tree topology used by the operation.
 * @param unitProperParts Proper parts treated as unit supports.
 * @param unitValues Unit-support values read or written by the operation.
 * @param attrNames Layout mapping attributes to buffer columns.
 * @param attribute Attribute requested by the operation.
 * @return Computed topology unit attribute with computers.
 */
template <std::floating_point Real, class... Computers>
inline bool computeTopologyUnitAttributeWithComputers(std::tuple<Computers...>, const MorphologicalTree& tree, std::span<const NodeId> unitProperParts,
                                                      std::span<Real> unitValues, const AttributeNames& attrNames, Attribute attribute) {
    bool computed = false;
    (void)std::initializer_list<int>{
        ((!computed && tryComputeTopologyUnitAttributeWithComputer<Computers>(tree, unitProperParts, unitValues, attrNames, attribute)) ? (computed = true, 0)
                                                                                                                                        : 0)...};
    return computed;
}

/**
 * @brief Computes topology unit attribute.
 *
 * @param tree Tree topology used by the operation.
 * @param unitProperParts Proper parts treated as unit supports.
 * @param unitValues Unit-support values read or written by the operation.
 * @param attrNames Layout mapping attributes to buffer columns.
 * @param attribute Attribute requested by the operation.
 */
template <std::floating_point Real>
inline void computeTopologyUnitAttribute(const MorphologicalTree& tree, std::span<const NodeId> unitProperParts, std::span<Real> unitValues,
                                         const AttributeNames& attrNames, Attribute attribute) {
    if (!computeTopologyUnitAttributeWithComputers(attributes::computers::RegisteredAttributeComputers{}, tree, unitProperParts, unitValues, attrNames,
                                                   attribute)) {
        throw std::runtime_error("No exported-Higra unit projection is registered for the requested attribute.");
    }
}

/**
 * @brief Computes topology unit attribute rows.
 *
 * @param tree Tree topology used by the operation.
 * @param unitProperParts Proper parts treated as unit supports.
 * @param attrNames Layout mapping attributes to buffer columns.
 * @return Computed topology unit attribute rows.
 */
template <std::floating_point Real = float>
inline std::vector<Real> computeTopologyUnitAttributeRows(const MorphologicalTree& tree, std::span<const NodeId> unitProperParts,
                                                          const AttributeNames& attrNames) {
    const size_t numColumns = static_cast<size_t>(attrNames.NUM_ATTRIBUTES);
    std::vector<Real> unitValues(unitProperParts.size() * numColumns, std::numeric_limits<Real>::quiet_NaN());

    for (const auto& [attribute, _] : attrNames.indexMap) {
        computeTopologyUnitAttribute<Real>(tree, unitProperParts, unitValues, attrNames, attribute);
    }

    return unitValues;
}

/**
 * @brief Attempts to compute altitude unit attribute with computer.
 *
 * @param context Operation name used in diagnostics.
 * @param attribute Attribute requested by the operation.
 * @return True when the documented condition holds; otherwise false.
 */
template <class Computer, std::floating_point Real, AltitudeValue T>
inline bool tryComputeAltitudeUnitAttributeWithComputer(const AltitudeUnitAttributeComputeContext<Real, T>& context, Attribute attribute) {
    if constexpr (Computer::domain == attributes::computers::AttributeComputerDomain::Altitude) {
        if (attributes::computers::producesAttribute<Computer>(attribute)) {
            Computer::computeUnitRows(context);
            return true;
        }
    }
    return false;
}

/**
 * @brief Computes altitude unit attribute with computers.
 *
 * @param context Operation name used in diagnostics.
 * @param attribute Attribute requested by the operation.
 * @return Computed altitude unit attribute with computers.
 */
template <std::floating_point Real, AltitudeValue T, class... Computers>
inline bool computeAltitudeUnitAttributeWithComputers(std::tuple<Computers...>, const AltitudeUnitAttributeComputeContext<Real, T>& context,
                                                      Attribute attribute) {
    bool computed = false;
    (void)std::initializer_list<int>{((!computed && tryComputeAltitudeUnitAttributeWithComputer<Computers>(context, attribute)) ? (computed = true, 0) : 0)...};
    return computed;
}

/**
 * @brief Computes generic altitude unit attribute.
 *
 * @param tree Tree topology used by the operation.
 * @param altitude Altitude values indexed by internal node identifier.
 * @param unitProperParts Proper parts treated as unit supports.
 * @param unitValues Unit-support values read or written by the operation.
 * @param attrNames Layout mapping attributes to buffer columns.
 * @param attribute Attribute requested by the operation.
 */
template <std::floating_point Real, AltitudeValue T>
inline void computeGenericAltitudeUnitAttribute(const MorphologicalTree& tree, std::span<const T> altitude, std::span<const NodeId> unitProperParts,
                                                std::span<Real> unitValues, const AttributeNames& attrNames, Attribute attribute) {
    if (!isGenericAltitudeUnitAttribute(attribute)) {
        return;
    }

    const std::array<Attribute, 1> requested{attribute};
    const AltitudeUnitAttributeComputeContext<Real, T> context{tree, altitude, unitProperParts, unitValues, attrNames, std::span<const Attribute>(requested)};

    if (!computeAltitudeUnitAttributeWithComputers(attributes::computers::RegisteredAttributeComputers{}, context, attribute)) {
        throw std::runtime_error("No altitude unit projection is registered for the requested attribute.");
    }
}

/**
 * @brief Computes unit attribute rows.
 *
 * @param tree Tree topology used by the operation.
 * @param altitude Altitude values indexed by internal node identifier.
 * @param unitProperParts Proper parts treated as unit supports.
 * @param attrNames Layout mapping attributes to buffer columns.
 * @return Computed unit attribute rows.
 */
template <std::floating_point Real = float, AltitudeValue T>
inline std::vector<Real> computeUnitAttributeRows(const MorphologicalTree& tree, std::span<const T> altitude, std::span<const NodeId> unitProperParts,
                                                  const AttributeNames& attrNames) {
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

/**
 * @brief Copies scalar unit rows to preserved higra.
 *
 * @param tree Tree topology used by the operation.
 * @param attrNames Layout mapping attributes to buffer columns.
 * @param unitValues Unit-support values read or written by the operation.
 * @param projected Destination buffer that receives the projected values.
 */
template <std::floating_point Real>
inline void copyScalarUnitRowsToPreservedHigra(const MorphologicalTree& tree, const AttributeNames& attrNames, std::span<const Real> unitValues,
                                               std::span<Real> projected) {
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

/**
 * @brief Copies scalar internal rows to preserved higra.
 *
 * @param tree Tree topology used by the operation.
 * @param attrNames Layout mapping attributes to buffer columns.
 * @param nodeValues Values indexed by internal node identifier.
 * @param projected Destination buffer that receives the projected values.
 * @param outputSize Count represented by `outputSize`.
 */
template <std::floating_point Real>
inline void copyScalarInternalRowsToPreservedHigra(const MorphologicalTree& tree, const AttributeNames& attrNames, std::span<const Real> nodeValues,
                                                   std::span<Real> projected, int outputSize) {
    const int numColumns = attrNames.NUM_ATTRIBUTES;
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        const NodeId outputNodeId = CommittedTreeAccess::higraNodeId(tree, nodeId);
        if (outputNodeId == InvalidNode || outputNodeId < 0 || outputNodeId >= outputSize) {
            throw std::runtime_error("Cannot project attributes to Higra node-id space: a live node has no valid Higra id.");
        }
        for (int column = 0; column < numColumns; ++column) {
            projected[static_cast<size_t>(outputNodeId) * static_cast<size_t>(numColumns) + static_cast<size_t>(column)] =
                nodeValues[static_cast<size_t>(nodeId) * static_cast<size_t>(numColumns) + static_cast<size_t>(column)];
        }
    }
}

/**
 * @brief Builds unit attribute names.
 *
 * @param deltaNames Names assigned to delta-dependent output fields.
 * @return Resulting unit attribute names.
 */
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

/**
 * @brief Copies delta unit rows to preserved higra.
 *
 * @param tree Tree topology used by the operation.
 * @param deltaNames Names assigned to delta-dependent output fields.
 * @param unitNames Unit labels assigned to projected output fields.
 * @param unitValues Unit-support values read or written by the operation.
 * @param projected Destination buffer that receives the projected values.
 */
template <std::floating_point Real>
inline void copyDeltaUnitRowsToPreservedHigra(const MorphologicalTree& tree, const AttributeNamesWithDelta& deltaNames, const AttributeNames& unitNames,
                                              std::span<const Real> unitValues, std::span<Real> projected) {
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

/**
 * @brief Copies delta internal rows to preserved higra.
 *
 * @param tree Tree topology used by the operation.
 * @param deltaNames Names assigned to delta-dependent output fields.
 * @param nodeValues Values indexed by internal node identifier.
 * @param projected Destination buffer that receives the projected values.
 * @param outputSize Count represented by `outputSize`.
 */
template <std::floating_point Real>
inline void copyDeltaInternalRowsToPreservedHigra(const MorphologicalTree& tree, const AttributeNamesWithDelta& deltaNames, std::span<const Real> nodeValues,
                                                  std::span<Real> projected, int outputSize) {
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        const NodeId outputNodeId = CommittedTreeAccess::higraNodeId(tree, nodeId);
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
 *
 * @param tree Tree topology used by the operation.
 * @param computed Flag controlling computed.
 * @param outputSpace Node-id domain used to index the output.
 * @return The projected scalar result from the internal node-id space to another public node-id space.
 */
template <std::floating_point Real>
inline ComputedAttributeData<Real> projectComputedDataToNodeIdSpace(const MorphologicalTree& tree, ComputedAttributeData<Real> computed,
                                                                    NodeIdSpace outputSpace) {
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

/**
 * @brief Projects computed data to node id space.
 *
 * @param tree Tree topology used by the operation.
 * @param altitude Altitude values indexed by internal node identifier.
 * @param computed Computed attribute result to transform.
 * @param outputSpace Node-id domain used to index the output.
 * @return Projected computed data to node id space.
 */
template <std::floating_point Real, AltitudeValue T>
inline ComputedAttributeData<Real> projectComputedDataToNodeIdSpace(const MorphologicalTree& tree, std::span<const T> altitude,
                                                                    ComputedAttributeData<Real> computed, NodeIdSpace outputSpace) {
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
 *
 * @param tree Tree topology used by the operation.
 * @param computed Flag controlling computed.
 * @param outputSpace Node-id domain used to index the output.
 * @return Projected delta-aware attribute data.
 */
template <std::floating_point Real>
inline ComputedAttributeDataWithDelta<Real> projectComputedDataToNodeIdSpace(const MorphologicalTree& tree, ComputedAttributeDataWithDelta<Real> computed,
                                                                             NodeIdSpace outputSpace) {
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
    const AttributeNames unitNames = deltaUnitAttributeNames(computed.first);
    const std::vector<Real> unitValues = computeTopologyUnitAttributeRows<Real>(tree, properParts, unitNames);
    copyDeltaUnitRowsToPreservedHigra<Real>(tree, computed.first, unitNames, unitValues, projected);
    copyDeltaInternalRowsToPreservedHigra<Real>(tree, computed.first, computed.second, projected, outputSize);

    return {std::move(computed.first), std::move(projected), outputSpace};
}

/**
 * @brief Projects computed data to node id space.
 *
 * @param tree Tree topology used by the operation.
 * @param altitude Altitude values indexed by internal node identifier.
 * @param computed Computed attribute result to transform.
 * @param outputSpace Node-id domain used to index the output.
 * @return Projected computed data to node id space.
 */
template <std::floating_point Real, AltitudeValue T>
inline ComputedAttributeDataWithDelta<Real> projectComputedDataToNodeIdSpace(const MorphologicalTree& tree, std::span<const T> altitude,
                                                                             ComputedAttributeDataWithDelta<Real> computed, NodeIdSpace outputSpace) {
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
    const AttributeNames unitNames = deltaUnitAttributeNames(computed.first);
    const std::vector<Real> unitValues = computeUnitAttributeRows<Real>(tree, altitude, properParts, unitNames);
    copyDeltaUnitRowsToPreservedHigra<Real>(tree, computed.first, unitNames, unitValues, projected);
    copyDeltaInternalRowsToPreservedHigra<Real>(tree, computed.first, computed.second, projected, outputSize);

    return {std::move(computed.first), std::move(projected), outputSpace};
}

/**
 * @brief Projects node values to exported higra typed.
 *
 * @param topology Tree topology used to propagate projected values.
 * @param altitude Altitude values indexed by internal node identifier.
 * @param attrNames Layout mapping attributes to buffer columns.
 * @param nodeValues Values indexed by internal node identifier.
 * @return Projected node values to exported higra typed.
 */
template <std::floating_point Real, AltitudeValue T>
inline std::vector<Real> projectNodeValuesToExportedHigraTyped(const MorphologicalTree& topology, std::span<const T> altitude, const AttributeNames& attrNames,
                                                               std::span<const Real> nodeValues) {
    const int numColumns = attrNames.NUM_ATTRIBUTES;
    const size_t expectedSize = static_cast<size_t>(topology.getNumInternalNodeSlots()) * static_cast<size_t>(numColumns);
    MMCFILTERS_CONTRACT_REQUIRE(
        nodeValues.size() == expectedSize,
        throw std::invalid_argument("Node-value buffer size must match the dense internal-node domain and requested attributes."));

    const auto layout = computeExportedHigraLayout(topology, altitude);
    const size_t numColumnValues = static_cast<size_t>(numColumns);
    const std::vector<Real> unitValues = computeUnitAttributeRows<Real>(topology, altitude, layout.properParts, attrNames);

    // Compact Higra buffers place unit proper-part rows before internal-node rows.
    std::vector<Real> projected(static_cast<size_t>(layout.numVertices) * numColumnValues, Real{0});

    for (NodeId leafIndex = 0; leafIndex < layout.numLeaves; ++leafIndex) {
        for (int column = 0; column < numColumns; ++column) {
            const size_t columnIndex = static_cast<size_t>(column);
            projected[static_cast<size_t>(leafIndex) * numColumnValues + columnIndex] =
                unitValues[static_cast<size_t>(leafIndex) * numColumnValues + columnIndex];
        }
    }

    for (NodeId nodeId : layout.sortedNodes) {
        const NodeId higraNodeId = layout.nodeToHigra[static_cast<size_t>(nodeId)];
        for (int column = 0; column < numColumns; ++column) {
            const size_t columnIndex = static_cast<size_t>(column);
            projected[static_cast<size_t>(higraNodeId) * numColumnValues + columnIndex] =
                nodeValues[static_cast<size_t>(nodeId) * numColumnValues + columnIndex];
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
 *
 * @param tree Tree topology used by the operation.
 * @param attrNames Layout that maps attributes to buffer columns.
 * @param nodeValues Dense values indexed by internal node identifier.
 * @param attribute Attribute requested by the operation.
 * @return The projected one scalar node attribute back to the original image domain.
 */
template <std::floating_point Real>
inline ImagePtr<Real> mapNodeAttributeToImage(const MorphologicalTree& tree, const AttributeNames& attrNames, std::span<const Real> nodeValues,
                                              Attribute attribute) {
    ImagePtr<Real> imgPtr = Image<Real>::create(tree.getNumRowsOfGridDomain2D(), tree.getNumColsOfGridDomain2D());
    Real* img = imgPtr->rawData();
    for (int p = 0; p < imgPtr->getSize(); ++p) {
        const NodeId nodeId = CommittedTreeAccess::properPartOwner(tree, p);
        img[p] = nodeValues[attrNames.linearIndex(nodeId, attribute)];
    }
    return imgPtr;
}

} // namespace mmcfilters::detail
