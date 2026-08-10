#pragma once

/**
 * @file AttributeDeltaMaterialization.hpp
 * @internal
 * @brief Internal helpers for delta-augmented attribute materialization.
 *
 * This header is installed only because the C++ core is header-only. It is not
 * part of the stable public API. Public callers should use
 * `AttributeComputation::computeSingleAttributeWithDelta(...)`.
 */

#include "../../trees/detail/TreeAltitudeDeltaNeighborhood.hpp"
#include "AttributeProjection.hpp"
#include "../AttributeResultTypes.hpp"
#include "../../trees/MorphologicalTree.hpp"
#include "../../trees/TreeAltitudeAlgorithms.hpp"
#include "../../utils/Altitude.hpp"
#include "../../utils/Common.hpp"
#include "../../utils/Contract.hpp"

#include <cmath>
#include <concepts>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mmcfilters::detail {

/**
 * @internal
 * @brief Validates the public padding strategy used by delta materialization.
 *
 * Supported strategies:
 * - `"last-padding"` repeats the nearest available center/previous value;
 * - `"nan-padding"` leaves missing samples as `NaN`;
 * - `"zero-padding"` fills missing samples with `0`;
 * - `"null-padding"` leaves the initially materialized `NaN` values unchanged.
 *
 * @param padding Padding strategy used outside the available hierarchy.
 *
 * @throws std::invalid_argument If `padding` is not one of the supported names.
 *
 */
inline void validateDeltaPaddingStrategy(const std::string& padding) {
    MMCFILTERS_CONTRACT_REQUIRE(padding == "last-padding" || padding == "nan-padding" || padding == "null-padding" || padding == "zero-padding",
                                throw std::invalid_argument("Unknown padding strategy."));
}

/**
 * @internal
 * @brief Rejects negative ancestor/descendant materialization radii.
 *
 * @param radius Maximum signed offset around the center attribute.
 * @param context Public API name used in diagnostics.
 *
 * @throws std::invalid_argument If `radius` is negative.
 */
inline void validateDeltaRadius(int radius, const char* context) {
    MMCFILTERS_CONTRACT_REQUIRE(radius >= 0,
                                throw std::invalid_argument(std::string(context) + " requires a non-negative radius."));
}

/**
 * @internal
 * @brief Materializes a delta layout with only the center value.
 *
 * This helper keeps the output shape compatible with
 * `ComputedAttributeDataWithDelta<float>` when no ancestor/descendant offsets are
 * requested. `base` must be expressed in the dense internal `MorphologicalTree`
 * node-id space because projection is applied only after the delta layout has
 * been assembled.
 *
 * @param tree Tree that defines the live dense internal node-id domain.
 * @param base Owning base attribute result in `NodeIdSpace::MORPHOLOGICAL_TREE`.
 * @param attribute Scalar attribute being materialized.
 * @param padding Padding strategy name validated for consistency with the typed
 * delta overload.
 * @param outputSpace Node-id space requested by the public caller.
 * @return Owning delta result projected to `outputSpace`.
 *
 * @throws std::invalid_argument If `padding` is unknown.
 * @throws std::logic_error If `base` is not in internal node-id space.
 */
template <std::floating_point Real>
inline ComputedAttributeDataWithDelta<Real> materializeSingleAttributeCenterOnly(const MorphologicalTree& tree, ComputedAttributeData<Real> base,
                                                                                 Attribute attribute, const std::string& padding, NodeIdSpace outputSpace) {
    validateDeltaPaddingStrategy(padding);
    if (base.nodeIdSpace != NodeIdSpace::MORPHOLOGICAL_TREE) {
        throw std::logic_error("Delta attribute materialization requires base attributes in internal node-id space.");
    }

    const AttributeNames& attributeNamesBase = base.first;
    const std::vector<Real>& attrsBase = base.second;
    const int n = tree.getNumInternalNodeSlots();
    const std::vector<Attribute> attrVec = {attribute};
    AttributeNamesWithDelta attributeNamesDelta(AttributeNamesWithDelta::create(0, attrVec));
    std::vector<Real> attrsDelta(static_cast<size_t>(n) * static_cast<size_t>(attributeNamesDelta.NUM_ATTRIBUTES), std::numeric_limits<Real>::quiet_NaN());

    for (NodeId nodeId : tree.getAliveNodeIds()) {
        const int outIdx = attributeNamesDelta.linearIndex(nodeId, attribute, 0);
        const int baseIdx = attributeNamesBase.linearIndex(nodeId, attribute);
        attrsDelta[static_cast<size_t>(outIdx)] = attrsBase[static_cast<size_t>(baseIdx)];
    }

    return projectComputedDataToNodeIdSpace(
        tree, ComputedAttributeDataWithDelta<Real>{std::move(attributeNamesDelta), std::move(attrsDelta), NodeIdSpace::MORPHOLOGICAL_TREE}, outputSpace);
}

/**
 * @internal
 * @brief Materializes ancestor/descendant samples for one scalar attribute.
 *
 * The output layout contains the center value at delta `0`, ancestor samples at
 * negative deltas, and descendant samples at positive deltas. For each distance
 * `d`, the altitude neighbourhood is selected with
 * `computeAscendantsAndDescendantsByAltitude(tree, altitude, deltaStep * d)`.
 * Missing samples are finalized according to `padding`.
 *
 * `base` must already contain `attribute` in internal node-id space. Dead node
 * slots keep their initialized `NaN` values; live-node iteration controls all
 * materialized writes.
 *
 * @tparam T Altitude type satisfying `AltitudeValue`.
 * @param tree Tree whose dense internal node-id domain indexes `base`.
 * @param altitude Dense altitude span indexed by internal `NodeId`.
 * @param base Owning base attribute result in `NodeIdSpace::MORPHOLOGICAL_TREE`.
 * @param attribute Scalar attribute being sampled.
 * @param deltaStep Altitude distance between consecutive sampled offsets.
 * @param radius Maximum positive and negative offset to materialize.
 * @param padding Strategy for missing ancestor/descendant samples.
 * @param outputSpace Node-id space requested by the public caller.
 * @return Owning delta result projected to `outputSpace`.
 *
 * @throws std::invalid_argument If radius, padding, altitude shape, or delta step
 * are invalid.
 * @throws std::logic_error If `base` is not in internal node-id space.
 *
 * Complexity: O(radius * (tree traversal + live nodes)) plus output projection.
 */
template <std::floating_point Real, AltitudeValue T>
inline ComputedAttributeDataWithDelta<Real>
materializeSingleAttributeWithTypedDelta(const MorphologicalTree& tree, std::span<const T> altitude, ComputedAttributeData<Real> base, Attribute attribute,
                                         AltitudeDiff<T> deltaStep, int radius, const std::string& padding, NodeIdSpace outputSpace) {
    validateDeltaPaddingStrategy(padding);
    validateDeltaRadius(radius, "computeSingleAttributeWithDelta");
    validateAltitudeDelta<T>(deltaStep, "computeSingleAttributeWithDelta");
    TreeAltitudeAlgorithms::validateAltitudeBufferShape(tree, altitude);
    if (base.nodeIdSpace != NodeIdSpace::MORPHOLOGICAL_TREE) {
        throw std::logic_error("Delta attribute materialization requires base attributes in internal node-id space.");
    }

    const AttributeNames& attributeNamesBase = base.first;
    const std::vector<Real>& attrsBase = base.second;
    const int n = tree.getNumInternalNodeSlots();
    const std::vector<Attribute> attrVec = {attribute};
    AttributeNamesWithDelta attributeNamesDelta(AttributeNamesWithDelta::create(radius, attrVec));

    std::vector<Real> attrsDelta(static_cast<size_t>(n) * static_cast<size_t>(attributeNamesDelta.NUM_ATTRIBUTES), std::numeric_limits<Real>::quiet_NaN());

    for (NodeId nodeId : tree.getAliveNodeIds()) {
        const NodeId nodeIndex = nodeId;
        const int outIdx = attributeNamesDelta.linearIndex(nodeIndex, attribute, 0);
        const int baseIdx = attributeNamesBase.linearIndex(nodeIndex, attribute);
        attrsDelta[static_cast<size_t>(outIdx)] = attrsBase[static_cast<size_t>(baseIdx)];
    }

    for (int d = 1; d <= radius; ++d) {
        const AltitudeDiff<T> distance = deltaStep * static_cast<AltitudeDiff<T>>(d);
        auto [ascendants, descendants] = computeAscendantsAndDescendantsByAltitude(tree, altitude, distance);

        for (NodeId nodeId : tree.getAliveNodeIds()) {
            const NodeId nodeIndex = nodeId;

            const NodeId ascIndex = (ascendants[static_cast<size_t>(nodeIndex)] != InvalidNode ? ascendants[static_cast<size_t>(nodeIndex)] : nodeIndex);
            if (ascIndex != nodeIndex) {
                const int outIdxAsc = attributeNamesDelta.linearIndex(nodeIndex, attribute, -d);
                const int baseIdxAsc = attributeNamesBase.linearIndex(ascIndex, attribute);
                attrsDelta[static_cast<size_t>(outIdxAsc)] = attrsBase[static_cast<size_t>(baseIdxAsc)];
            }

            const NodeId descIndex = (descendants[static_cast<size_t>(nodeIndex)] != InvalidNode ? descendants[static_cast<size_t>(nodeIndex)] : nodeIndex);
            if (descIndex != nodeIndex) {
                const int outIdxDesc = attributeNamesDelta.linearIndex(nodeIndex, attribute, +d);
                const int baseIdxDesc = attributeNamesBase.linearIndex(descIndex, attribute);
                attrsDelta[static_cast<size_t>(outIdxDesc)] = attrsBase[static_cast<size_t>(baseIdxDesc)];
            }
        }
    }

    if (padding == "last-padding" || padding == "nan-padding") {
        for (NodeId nodeId : tree.getAliveNodeIds()) {
            const NodeId nodeIndex = nodeId;

            for (int d = 1; d <= radius; ++d) {
                const int outIdx = attributeNamesDelta.linearIndex(nodeIndex, attribute, -d);
                const int refIdx = attributeNamesDelta.linearIndex(nodeIndex, attribute, -(d - 1));

                if (std::isnan(attrsDelta[static_cast<size_t>(outIdx)])) {
                    if (padding == "last-padding") {
                        attrsDelta[static_cast<size_t>(outIdx)] = attrsDelta[static_cast<size_t>(refIdx)];
                    } else {
                        attrsDelta[static_cast<size_t>(outIdx)] = std::numeric_limits<Real>::quiet_NaN();
                    }
                }
            }

            for (int d = 1; d <= radius; ++d) {
                const int outIdx = attributeNamesDelta.linearIndex(nodeIndex, attribute, +d);
                const int refIdx = attributeNamesDelta.linearIndex(nodeIndex, attribute, +(d - 1));

                if (CommittedTreeAccess::numChildren(tree, nodeIndex) == 0 || std::isnan(attrsDelta[static_cast<size_t>(outIdx)])) {
                    if (padding == "last-padding") {
                        attrsDelta[static_cast<size_t>(outIdx)] = attrsDelta[static_cast<size_t>(refIdx)];
                    } else {
                        attrsDelta[static_cast<size_t>(outIdx)] = std::numeric_limits<Real>::quiet_NaN();
                    }
                }
            }
        }
    } else if (padding == "zero-padding") {
        for (NodeId nodeId : tree.getAliveNodeIds()) {
            const NodeId nodeIndex = nodeId;
            for (int d = 1; d <= radius; ++d) {
                const int ascIdx = attributeNamesDelta.linearIndex(nodeIndex, attribute, -d);
                const int descIdx = attributeNamesDelta.linearIndex(nodeIndex, attribute, +d);
                if (std::isnan(attrsDelta[static_cast<size_t>(ascIdx)])) {
                    attrsDelta[static_cast<size_t>(ascIdx)] = Real{0};
                }
                if (std::isnan(attrsDelta[static_cast<size_t>(descIdx)])) {
                    attrsDelta[static_cast<size_t>(descIdx)] = Real{0};
                }
            }
        }
    }

    return projectComputedDataToNodeIdSpace(
        tree, altitude, ComputedAttributeDataWithDelta<Real>{std::move(attributeNamesDelta), std::move(attrsDelta), NodeIdSpace::MORPHOLOGICAL_TREE},
        outputSpace);
}

} // namespace mmcfilters::detail
