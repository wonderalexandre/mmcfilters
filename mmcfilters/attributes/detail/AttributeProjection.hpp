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
 * @param tree Tree topology.
 * @return Direct proper-part pixel identifiers in row-major order.
 */
inline std::vector<PixelId> makeRowMajorProperPartOrder(const MorphologicalTree& tree) {
    std::vector<PixelId> properParts;
    properParts.reserve(static_cast<size_t>(tree.numPixels()));
    for (PixelId pixel = 0; pixel < tree.numPixels(); ++pixel) {
        properParts.push_back(pixel);
    }
    return properParts;
}

/**
 * @brief Computes unit attribute with computer.
 *
 * @param tree Tree topology.
 * @param unitPixels Proper parts treated as unit supports.
 * @param unitValues Unit-support values read or written by the operation.
 * @param attrNames Layout mapping attributes to buffer columns.
 * @param attribute Attribute requested by the operation.
 */
template <class Computer, std::floating_point Real>
inline void computeUnitAttributeWithComputer(const MorphologicalTree& tree, std::span<const PixelId> unitPixels, std::span<Real> unitValues,
                                             const AttributeNames& attrNames, Attribute attribute) {
    const std::array<Attribute, 1> requested{attribute};
    Computer::computeUnitRows(UnitAttributeComputeContext<Real>{tree, unitPixels, unitValues, attrNames, std::span<const Attribute>(requested)});
}

/**
 * @brief Attempts to compute topology unit attribute with computer.
 *
 * @param tree Tree topology.
 * @param unitPixels Proper parts treated as unit supports.
 * @param unitValues Unit-support values read or written by the operation.
 * @param attrNames Layout mapping attributes to buffer columns.
 * @param attribute Attribute requested by the operation.
 * @return True when the documented condition holds; otherwise false.
 */
template <class Computer, std::floating_point Real>
inline bool tryComputeTopologyUnitAttributeWithComputer(const MorphologicalTree& tree, std::span<const PixelId> unitPixels, std::span<Real> unitValues,
                                                        const AttributeNames& attrNames, Attribute attribute) {
    if constexpr (Computer::domain == attributes::computers::AttributeComputerDomain::Topology) {
        if (attributes::computers::producesAttribute<Computer>(attribute)) {
            computeUnitAttributeWithComputer<Computer>(tree, unitPixels, unitValues, attrNames, attribute);
            return true;
        }
    }
    return false;
}

/**
 * @brief Computes topology unit attribute with computers.
 *
 * @param tree Tree topology.
 * @param unitPixels Proper parts treated as unit supports.
 * @param unitValues Unit-support values read or written by the operation.
 * @param attrNames Layout mapping attributes to buffer columns.
 * @param attribute Attribute requested by the operation.
 * @return Computed topology unit attribute with computers.
 */
template <std::floating_point Real, class... Computers>
inline bool computeTopologyUnitAttributeWithComputers(std::tuple<Computers...>, const MorphologicalTree& tree, std::span<const PixelId> unitPixels,
                                                      std::span<Real> unitValues, const AttributeNames& attrNames, Attribute attribute) {
    bool computed = false;
    (void)std::initializer_list<int>{
        ((!computed && tryComputeTopologyUnitAttributeWithComputer<Computers>(tree, unitPixels, unitValues, attrNames, attribute)) ? (computed = true, 0)
                                                                                                                                        : 0)...};
    return computed;
}

/**
 * @brief Computes topology unit attribute.
 *
 * @param tree Tree topology.
 * @param unitPixels Proper parts treated as unit supports.
 * @param unitValues Unit-support values read or written by the operation.
 * @param attrNames Layout mapping attributes to buffer columns.
 * @param attribute Attribute requested by the operation.
 */
template <std::floating_point Real>
inline void computeTopologyUnitAttribute(const MorphologicalTree& tree, std::span<const PixelId> unitPixels, std::span<Real> unitValues,
                                         const AttributeNames& attrNames, Attribute attribute) {
    if (!computeTopologyUnitAttributeWithComputers(attributes::computers::RegisteredAttributeComputers{}, tree, unitPixels, unitValues, attrNames,
                                                   attribute)) {
        throw std::runtime_error("No exported-Higra unit projection is registered for the requested attribute.");
    }
}

/**
 * @brief Computes topology unit attribute rows.
 *
 * @param tree Tree topology.
 * @param unitPixels Proper parts treated as unit supports.
 * @param attrNames Layout mapping attributes to buffer columns.
 * @return Computed topology unit attribute rows.
 */
template <std::floating_point Real = float>
inline std::vector<Real> computeTopologyUnitAttributeRows(const MorphologicalTree& tree, std::span<const PixelId> unitPixels,
                                                          const AttributeNames& attrNames) {
    const size_t numColumns = static_cast<size_t>(attrNames.NUM_ATTRIBUTES);
    std::vector<Real> unitValues(unitPixels.size() * numColumns, std::numeric_limits<Real>::quiet_NaN());

    for (const auto& [attribute, _] : attrNames.indexMap) {
        computeTopologyUnitAttribute<Real>(tree, unitPixels, unitValues, attrNames, attribute);
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
 * @param tree Tree topology.
 * @param altitude Altitude values indexed by internal node identifier.
 * @param unitPixels Proper parts treated as unit supports.
 * @param unitValues Unit-support values read or written by the operation.
 * @param attrNames Layout mapping attributes to buffer columns.
 * @param attribute Attribute requested by the operation.
 */
template <std::floating_point Real, AltitudeValue T>
inline void computeGenericAltitudeUnitAttribute(const MorphologicalTree& tree, std::span<const T> altitude, std::span<const PixelId> unitPixels,
                                                std::span<Real> unitValues, const AttributeNames& attrNames, Attribute attribute) {
    if (!isGenericAltitudeUnitAttribute(attribute)) {
        return;
    }

    const std::array<Attribute, 1> requested{attribute};
    const AltitudeUnitAttributeComputeContext<Real, T> context{tree, altitude, unitPixels, unitValues, attrNames, std::span<const Attribute>(requested)};

    if (!computeAltitudeUnitAttributeWithComputers(attributes::computers::RegisteredAttributeComputers{}, context, attribute)) {
        throw std::runtime_error("No altitude unit projection is registered for the requested attribute.");
    }
}

/**
 * @brief Computes unit attribute rows.
 *
 * @param tree Tree topology.
 * @param altitude Altitude values indexed by internal node identifier.
 * @param unitPixels Proper parts treated as unit supports.
 * @param attrNames Layout mapping attributes to buffer columns.
 * @return Computed unit attribute rows.
 */
template <std::floating_point Real = float, AltitudeValue T>
inline std::vector<Real> computeUnitAttributeRows(const MorphologicalTree& tree, std::span<const T> altitude, std::span<const PixelId> unitPixels,
                                                  const AttributeNames& attrNames) {
    const size_t numColumns = static_cast<size_t>(attrNames.NUM_ATTRIBUTES);
    std::vector<Real> unitValues(unitPixels.size() * numColumns, std::numeric_limits<Real>::quiet_NaN());

    for (const auto& [attribute, _] : attrNames.indexMap) {
        if (isGenericAltitudeUnitAttribute(attribute)) {
            computeGenericAltitudeUnitAttribute<Real>(tree, altitude, unitPixels, unitValues, attrNames, attribute);
            continue;
        }

        computeTopologyUnitAttribute<Real>(tree, unitPixels, unitValues, attrNames, attribute);
    }

    return unitValues;
}

/**
 * @brief Copies scalar unit rows to preserved higra.
 *
 * @param tree Tree topology.
 * @param attrNames Layout mapping attributes to buffer columns.
 * @param unitValues Unit-support values read or written by the operation.
 * @param projected Destination buffer that receives the projected values.
 */
template <std::floating_point Real>
inline void copyScalarUnitRowsToPreservedHigra(const MorphologicalTree& tree, const AttributeNames& attrNames, std::span<const Real> unitValues,
                                               std::span<Real> projected) {
    const int numPixels = tree.numPixels();
    const int numColumns = attrNames.NUM_ATTRIBUTES;
    const size_t expectedUnitSize = static_cast<size_t>(numPixels) * static_cast<size_t>(numColumns);
    if (unitValues.size() != expectedUnitSize || projected.size() < expectedUnitSize) {
        throw std::logic_error("Preserved Higra unit-row projection received incompatible buffer shapes.");
    }

    for (NodeId leafIndex = 0; leafIndex < numPixels; ++leafIndex) {
        for (int column = 0; column < numColumns; ++column) {
            projected[static_cast<size_t>(leafIndex) * static_cast<size_t>(numColumns) + static_cast<size_t>(column)] =
                unitValues[static_cast<size_t>(leafIndex) * static_cast<size_t>(numColumns) + static_cast<size_t>(column)];
        }
    }
}

/**
 * @brief Copies scalar internal rows to preserved higra.
 *
 * @param tree Tree topology.
 * @param attrNames Layout mapping attributes to buffer columns.
 * @param nodeValues Values indexed by internal node identifier.
 * @param projected Destination buffer that receives the projected values.
 * @param outputSize Count.
 */
template <std::floating_point Real>
inline void copyScalarInternalRowsToPreservedHigra(const MorphologicalTree& tree, const AttributeNames& attrNames, std::span<const Real> nodeValues,
                                                   std::span<Real> projected, int outputSize) {
    const int numColumns = attrNames.NUM_ATTRIBUTES;
    for (NodeId nodeId : tree.aliveNodeIds()) {
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
 * @param sampleLayout Layout assigned to sampled output fields.
 * @return Resulting unit attribute names.
 */
inline AttributeNames sampledUnitAttributeNames(const NodeAttributeSampleLayout& sampleLayout) {
    std::vector<Attribute> attributes;
    attributes.reserve(sampleLayout.indexMap.size());
    for (const auto& [sampleKey, _] : sampleLayout.indexMap) {
        if (std::find(attributes.begin(), attributes.end(), sampleKey.attribute) == attributes.end()) {
            attributes.push_back(sampleKey.attribute);
        }
    }
    return AttributeNames::fromList(attributes);
}

/**
 * @brief Copies sampled unit rows to preserved Higra space.
 *
 * @param tree Tree topology.
 * @param sampleLayout Layout assigned to sampled output fields.
 * @param unitNames Unit labels assigned to projected output fields.
 * @param unitValues Unit-support values read or written by the operation.
 * @param projected Destination buffer that receives the projected values.
 */
template <std::floating_point Real>
inline void copySampledUnitRowsToPreservedHigra(const MorphologicalTree& tree, const NodeAttributeSampleLayout& sampleLayout,
                                                const AttributeNames& unitNames, std::span<const Real> unitValues, std::span<Real> projected) {
    const int numPixels = tree.numPixels();
    const int numColumns = sampleLayout.NUM_ATTRIBUTES;
    const size_t expectedUnitSize = static_cast<size_t>(numPixels) * static_cast<size_t>(unitNames.NUM_ATTRIBUTES);
    const size_t expectedProjectedUnitSize = static_cast<size_t>(numPixels) * static_cast<size_t>(numColumns);
    if (unitValues.size() != expectedUnitSize || projected.size() < expectedProjectedUnitSize) {
        throw std::logic_error("Preserved Higra sampled unit-row projection received incompatible buffer shapes.");
    }

    // Proper parts have no ancestor/descendant neighbourhood in this domain, so
    // every sampled column receives the unit value for its scalar attribute.
    for (NodeId leafIndex = 0; leafIndex < numPixels; ++leafIndex) {
        for (const auto& [sampleKey, _] : sampleLayout.indexMap) {
            projected[static_cast<size_t>(sampleLayout.linearIndex(leafIndex, sampleKey))] =
                unitValues[static_cast<size_t>(unitNames.linearIndex(leafIndex, sampleKey.attribute))];
        }
    }
}

/**
 * @brief Copies sampled internal rows to preserved Higra space.
 *
 * @param tree Tree topology.
 * @param sampleLayout Layout assigned to sampled output fields.
 * @param nodeValues Values indexed by internal node identifier.
 * @param projected Destination buffer that receives the projected values.
 * @param outputSize Count.
 */
template <std::floating_point Real>
inline void copySampledInternalRowsToPreservedHigra(const MorphologicalTree& tree, const NodeAttributeSampleLayout& sampleLayout,
                                                    std::span<const Real> nodeValues, std::span<Real> projected, int outputSize) {
    for (NodeId nodeId : tree.aliveNodeIds()) {
        const NodeId outputNodeId = CommittedTreeAccess::higraNodeId(tree, nodeId);
        if (outputNodeId == InvalidNode || outputNodeId < 0 || outputNodeId >= outputSize) {
            throw std::runtime_error("Cannot project sampled attributes to Higra node-id space: a live node has no valid Higra id.");
        }
        for (const auto& [sampleKey, _] : sampleLayout.indexMap) {
            projected[static_cast<size_t>(sampleLayout.linearIndex(outputNodeId, sampleKey.attribute, sampleKey.sampleOffset))] =
                nodeValues[static_cast<size_t>(sampleLayout.linearIndex(nodeId, sampleKey.attribute, sampleKey.sampleOffset))];
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
 * @param tree Tree topology.
 * @param computed Flag controlling computed.
 * @param outputSpace Node-id domain used to index the output.
 * @return The projected scalar result from the internal node-id space to another public node-id space.
 */
template <std::floating_point Real>
inline ComputedAttributeData<Real> projectComputedDataToNodeIdSpace(const MorphologicalTree& tree, ComputedAttributeData<Real> computed,
                                                                    NodeIdSpace outputSpace) {
    if (outputSpace == NodeIdSpace::MorphologicalTree) {
        computed.nodeIdSpace = NodeIdSpace::MorphologicalTree;
        return computed;
    }

    const int numAttributes = computed.first.NUM_ATTRIBUTES;
    const int outputSize = tree.getNodeIdSpaceSize(outputSpace);
    std::vector<Real> projected(static_cast<size_t>(outputSize) * static_cast<size_t>(numAttributes), std::numeric_limits<Real>::quiet_NaN());

    if (outputSpace != NodeIdSpace::Higra) {
        throw std::runtime_error("Unsupported node-id projection target.");
    }

    const std::vector<PixelId> properParts = makeRowMajorProperPartOrder(tree);
    const std::vector<Real> unitValues = computeTopologyUnitAttributeRows<Real>(tree, properParts, computed.first);
    copyScalarUnitRowsToPreservedHigra<Real>(tree, computed.first, unitValues, projected);
    copyScalarInternalRowsToPreservedHigra<Real>(tree, computed.first, computed.second, projected, outputSize);

    return {std::move(computed.first), std::move(projected), outputSpace};
}

/**
 * @brief Projects computed data to node id space.
 *
 * @param tree Tree topology.
 * @param altitude Altitude values indexed by internal node identifier.
 * @param computed Computed attribute result to transform.
 * @param outputSpace Node-id domain used to index the output.
 * @return Projected computed data to node id space.
 */
template <std::floating_point Real, AltitudeValue T>
inline ComputedAttributeData<Real> projectComputedDataToNodeIdSpace(const MorphologicalTree& tree, std::span<const T> altitude,
                                                                    ComputedAttributeData<Real> computed, NodeIdSpace outputSpace) {
    if (outputSpace == NodeIdSpace::MorphologicalTree) {
        computed.nodeIdSpace = NodeIdSpace::MorphologicalTree;
        return computed;
    }

    const int numAttributes = computed.first.NUM_ATTRIBUTES;
    const int outputSize = tree.getNodeIdSpaceSize(outputSpace);
    std::vector<Real> projected(static_cast<size_t>(outputSize) * static_cast<size_t>(numAttributes), std::numeric_limits<Real>::quiet_NaN());

    if (outputSpace != NodeIdSpace::Higra) {
        throw std::runtime_error("Unsupported node-id projection target.");
    }

    const std::vector<PixelId> properParts = makeRowMajorProperPartOrder(tree);
    const std::vector<Real> unitValues = computeUnitAttributeRows<Real>(tree, altitude, properParts, computed.first);
    copyScalarUnitRowsToPreservedHigra<Real>(tree, computed.first, unitValues, projected);
    copyScalarInternalRowsToPreservedHigra<Real>(tree, computed.first, computed.second, projected, outputSize);

    return {std::move(computed.first), std::move(projected), outputSpace};
}

/**
 * @brief Sampled-attribute counterpart of the scalar node-id-space projection.
 *
 * @param tree Tree topology.
 * @param computed Flag controlling computed.
 * @param outputSpace Node-id domain used to index the output.
 * @return Projected sampled node-attribute data.
 */
template <std::floating_point Real>
inline SampledNodeAttributeData<Real> projectComputedDataToNodeIdSpace(const MorphologicalTree& tree, SampledNodeAttributeData<Real> computed,
                                                                             NodeIdSpace outputSpace) {
    if (outputSpace == NodeIdSpace::MorphologicalTree) {
        computed.nodeIdSpace = NodeIdSpace::MorphologicalTree;
        return computed;
    }

    const int numAttributes = computed.first.NUM_ATTRIBUTES;
    const int outputSize = tree.getNodeIdSpaceSize(outputSpace);
    std::vector<Real> projected(static_cast<size_t>(outputSize) * static_cast<size_t>(numAttributes), std::numeric_limits<Real>::quiet_NaN());

    if (outputSpace != NodeIdSpace::Higra) {
        throw std::runtime_error("Unsupported node-id projection target.");
    }

    const std::vector<PixelId> properParts = makeRowMajorProperPartOrder(tree);
    const AttributeNames unitNames = sampledUnitAttributeNames(computed.first);
    const std::vector<Real> unitValues = computeTopologyUnitAttributeRows<Real>(tree, properParts, unitNames);
    copySampledUnitRowsToPreservedHigra<Real>(tree, computed.first, unitNames, unitValues, projected);
    copySampledInternalRowsToPreservedHigra<Real>(tree, computed.first, computed.second, projected, outputSize);

    return {std::move(computed.first), std::move(projected), outputSpace};
}

/**
 * @brief Projects computed data to node id space.
 *
 * @param tree Tree topology.
 * @param altitude Altitude values indexed by internal node identifier.
 * @param computed Computed attribute result to transform.
 * @param outputSpace Node-id domain used to index the output.
 * @return Projected computed data to node id space.
 */
template <std::floating_point Real, AltitudeValue T>
inline SampledNodeAttributeData<Real> projectComputedDataToNodeIdSpace(const MorphologicalTree& tree, std::span<const T> altitude,
                                                                             SampledNodeAttributeData<Real> computed, NodeIdSpace outputSpace) {
    if (outputSpace == NodeIdSpace::MorphologicalTree) {
        computed.nodeIdSpace = NodeIdSpace::MorphologicalTree;
        return computed;
    }

    const int numAttributes = computed.first.NUM_ATTRIBUTES;
    const int outputSize = tree.getNodeIdSpaceSize(outputSpace);
    std::vector<Real> projected(static_cast<size_t>(outputSize) * static_cast<size_t>(numAttributes), std::numeric_limits<Real>::quiet_NaN());

    if (outputSpace != NodeIdSpace::Higra) {
        throw std::runtime_error("Unsupported node-id projection target.");
    }

    const std::vector<PixelId> properParts = makeRowMajorProperPartOrder(tree);
    const AttributeNames unitNames = sampledUnitAttributeNames(computed.first);
    const std::vector<Real> unitValues = computeUnitAttributeRows<Real>(tree, altitude, properParts, unitNames);
    copySampledUnitRowsToPreservedHigra<Real>(tree, computed.first, unitNames, unitValues, projected);
    copySampledInternalRowsToPreservedHigra<Real>(tree, computed.first, computed.second, projected, outputSize);

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
    const size_t expectedSize = static_cast<size_t>(topology.numInternalNodeSlots()) * static_cast<size_t>(numColumns);
    MMCFILTERS_CONTRACT_REQUIRE(
        nodeValues.size() == expectedSize,
        throw std::invalid_argument("Node-value buffer size must match the dense internal-node domain and requested attributes."));

    const auto layout = computeExportedHigraLayout(topology, altitude);
    const size_t numColumnValues = static_cast<size_t>(numColumns);
    const std::vector<Real> unitValues = computeUnitAttributeRows<Real>(topology, altitude, layout.properParts, attrNames);

    // Compact Higra buffers place rows for pixel leaves before internal-node rows.
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
 * Each proper part receives the value stored at its smallest node. The input
 * buffer must use the canonical dense internal-node layout described by
 * `attrNames`.
 *
 * @param tree Tree topology.
 * @param attrNames Layout that maps attributes to buffer columns.
 * @param nodeValues Dense values indexed by internal node identifier.
 * @param attribute Attribute requested by the operation.
 * @return The projected one scalar node attribute back to the original image domain.
 */
template <std::floating_point Real>
inline ImagePtr<Real> mapNodeAttributeToImage(const MorphologicalTree& tree, const AttributeNames& attrNames, std::span<const Real> nodeValues,
                                              Attribute attribute) {
    ImagePtr<Real> imgPtr = Image<Real>::create(tree.numRows(), tree.numColumns());
    Real* img = imgPtr->rawData();
    for (int p = 0; p < imgPtr->getSize(); ++p) {
        const NodeId nodeId = CommittedTreeAccess::smallestNodeMap(tree, p);
        img[p] = nodeValues[attrNames.linearIndex(nodeId, attribute)];
    }
    return imgPtr;
}

} // namespace mmcfilters::detail
