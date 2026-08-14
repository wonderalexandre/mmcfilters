#pragma once

#include "DepthStableRegionComputer.hpp"
#include "MSERComputer.hpp"
#include "NodeDecisionMasks.hpp"
#include "../trees/MorphologicalTree.hpp"
#include "../trees/ValuedMorphologicalTree.hpp"
#include "../utils/Altitude.hpp"
#include "../utils/Contract.hpp"

#include <concepts>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace mmcfilters {

/** @brief Policy used when a stability window cannot be evaluated. */
enum class IncompleteStabilityWindowPolicy {
    /** @brief Retain the input preservation decision at the rejected node. */
    PreserveInputDecision
};

/**
 * @brief Thresholds scalar node attributes into preservation decisions.
 *
 * A node is preserved exactly when `nodeAttributes[node] >= threshold`.
 * Thresholding is deliberately separate from stability adjustment.
 */
template <class AttributeValue>
    requires(std::is_arithmetic_v<AttributeValue> && !std::is_same_v<std::remove_cv_t<AttributeValue>, bool>)
[[nodiscard]] NodePreservationMask computeNodePreservationMask(std::span<const AttributeValue> nodeAttributes, AttributeValue threshold) {
    std::vector<bool> nodePreservationDecisions(nodeAttributes.size(), false);
    for (std::size_t index = 0; index < nodeAttributes.size(); ++index) {
        nodePreservationDecisions[index] = nodeAttributes[index] >= threshold;
    }
    return NodePreservationMask(std::move(nodePreservationDecisions));
}

namespace detail::node_preservation_stability {

inline void requireNodePreservationMaskShape(const MorphologicalTree& tree, const NodePreservationMask& nodePreservationMask,
                                             const char* context) {
    MMCFILTERS_CONTRACT_REQUIRE(
        nodePreservationMask.size() == static_cast<std::size_t>(tree.numInternalNodeSlots()),
        throw std::invalid_argument(std::string(context) + " nodePreservationMask size must match the internal node slot count."));
}

template <std::floating_point Real, class StabilityComputer>
[[nodiscard]] NodePreservationMask adjustNodePreservationMask(const MorphologicalTree& tree,
                                                              const NodePreservationMask& nodePreservationMask,
                                                              StabilityComputer& stabilityComputer,
                                                              IncompleteStabilityWindowPolicy incompleteWindowPolicy,
                                                              const char* context) {
    requireNodePreservationMaskShape(tree, nodePreservationMask, context);
    const std::vector<Real>& variations = stabilityComputer.getVariations();
    std::vector<bool> adjustedDecisions(static_cast<std::size_t>(tree.numInternalNodeSlots()), true);

    for (NodeId nodeId : tree.aliveNodeIds()) {
        if (nodePreservationMask[static_cast<std::size_t>(nodeId)]) {
            continue;
        }

        if (detail::isFiniteVariation(variations[static_cast<std::size_t>(nodeId)])) {
            const NodeId adjustedRejection = stabilityComputer.nodeWithMinimumVariationInWindow(nodeId);
            adjustedDecisions[static_cast<std::size_t>(adjustedRejection)] = false;
            continue;
        }

        switch (incompleteWindowPolicy) {
        case IncompleteStabilityWindowPolicy::PreserveInputDecision:
            adjustedDecisions[static_cast<std::size_t>(nodeId)] = false;
            break;
        }
    }

    return NodePreservationMask(std::move(adjustedDecisions));
}

} // namespace detail::node_preservation_stability

/**
 * @brief Relocates input rejections using an altitude-distance stability window.
 *
 * The valued tree must declare a globally monotone altitude order.
 */
template <AltitudeValue T>
[[nodiscard]] NodePreservationMask adjustNodePreservationMaskByAltitudeStability(
    const ValuedMorphologicalTree<T>& valuedTree, const NodePreservationMask& nodePreservationMask,
    AltitudeDifference<T> altitudeWindowRadius,
    IncompleteStabilityWindowPolicy incompleteWindowPolicy = IncompleteStabilityWindowPolicy::PreserveInputDecision) {
    MSERComputer<T> stabilityComputer(valuedTree);
    (void)stabilityComputer.computeMSER(altitudeWindowRadius);
    return detail::node_preservation_stability::adjustNodePreservationMask<float>(
        valuedTree.topology(), nodePreservationMask, stabilityComputer, incompleteWindowPolicy,
        "adjustNodePreservationMaskByAltitudeStability");
}

/** @brief Relocates input rejections using an edge-count stability window. */
[[nodiscard]] inline NodePreservationMask adjustNodePreservationMaskByDepthStability(
    const MorphologicalTree& tree, const NodePreservationMask& nodePreservationMask, int depthWindowRadius,
    IncompleteStabilityWindowPolicy incompleteWindowPolicy = IncompleteStabilityWindowPolicy::PreserveInputDecision) {
    DepthStableRegionComputer<float> stabilityComputer(tree);
    (void)stabilityComputer.computeByDepth(depthWindowRadius);
    return detail::node_preservation_stability::adjustNodePreservationMask<float>(
        tree, nodePreservationMask, stabilityComputer, incompleteWindowPolicy,
        "adjustNodePreservationMaskByDepthStability");
}

/** @brief Valued-tree convenience overload for topology-only depth stability. */
template <AltitudeValue T>
[[nodiscard]] inline NodePreservationMask adjustNodePreservationMaskByDepthStability(
    const ValuedMorphologicalTree<T>& valuedTree, const NodePreservationMask& nodePreservationMask, int depthWindowRadius,
    IncompleteStabilityWindowPolicy incompleteWindowPolicy = IncompleteStabilityWindowPolicy::PreserveInputDecision) {
    return adjustNodePreservationMaskByDepthStability(valuedTree.topology(), nodePreservationMask, depthWindowRadius,
                                                      incompleteWindowPolicy);
}

} // namespace mmcfilters
