#pragma once

#include "../../trees/MorphologicalTree.hpp"
#include "../../utils/Common.hpp"

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace mmcfilters::detail {

/**
 * @brief Tie-breaking policy used when preserving and removing a node have the same Viterbi cost.
 *
 * The default used by `AttributeFilters::filteringByViterbiRule` is
 * `PreferRemove`. With the threshold cost model in this file, equal costs happen
 * exactly at `attribute[node] == threshold`; preferring remove keeps the public
 * behavior aligned with the strict `attribute > threshold` convention already
 * used by the pruning attribute overloads.
 */
enum class ViterbiTieBreak {
    PreferRemove,
    PreferPreserve,
};

/**
 * @brief Options for the connected Viterbi decision on a morphological tree.
 */
struct ViterbiDecisionOptions {
    /** @brief Stores the tie break. */
    ViterbiTieBreak tieBreak = ViterbiTieBreak::PreferRemove;
};

/**
 * @brief Per-node unary costs for the two Viterbi states.
 *
 * Costs are indexed by the dense internal `NodeId` slot domain. Only alive
 * nodes reachable from the current root participate in the dynamic program,
 * but the vectors must still match `tree.getNumInternalNodeSlots()` so the
 * caller cannot accidentally mix Higra ids, compacted alive ids, and internal
 * node slots.
 */
template <std::floating_point Real> struct ViterbiNodeCosts {
    /** @brief Stores the preserve. */
    std::vector<Real> preserve;
    /** @brief Stores the remove. */
    std::vector<Real> remove;
};

namespace viterbi_decision_detail {

/**
 * @brief Validates finite non negative cost.
 *
 * @param value Value used by the operation.
 * @param nodeId Identifier of the node used by the operation.
 * @param costName Name of the cost term included in diagnostics.
 */
template <std::floating_point Real> void requireFiniteNonNegativeCost(Real value, NodeId nodeId, std::string_view costName) {
    if (!std::isfinite(value) || value < Real{}) {
        std::ostringstream oss;
        oss << "Viterbi " << costName << " cost must be finite and non-negative at node " << nodeId << ".";
        throw std::invalid_argument(oss.str());
    }
}

/**
 * @brief Validates finite attribute.
 *
 * @param value Value used by the operation.
 * @param nodeId Identifier of the node used by the operation.
 */
template <std::floating_point Real> void requireFiniteAttribute(Real value, NodeId nodeId) {
    if (!std::isfinite(value)) {
        std::ostringstream oss;
        oss << "Viterbi attribute value must be finite at node " << nodeId << ".";
        throw std::invalid_argument(oss.str());
    }
}

/**
 * @brief Adds saturating.
 *
 * @param lhs Left-hand operand.
 * @param rhs Right-hand operand.
 * @return Finite sum, or positive infinity when the addition saturates.
 */
template <std::floating_point Real> Real addSaturating(Real lhs, Real rhs) {
    if (lhs == std::numeric_limits<Real>::infinity() || rhs == std::numeric_limits<Real>::infinity()) {
        return std::numeric_limits<Real>::infinity();
    }
    if (lhs > std::numeric_limits<Real>::max() - rhs) {
        return std::numeric_limits<Real>::infinity();
    }
    return lhs + rhs;
}

/**
 * @brief Chooses preserve.
 *
 * @param preserveCost Cost of preserving the current node.
 * @param removeCost Cost of removing the current node.
 * @param options Options controlling the operation.
 * @return Selected preserve.
 */
template <std::floating_point Real> bool choosePreserve(Real preserveCost, Real removeCost, const ViterbiDecisionOptions& options) {
    if (preserveCost < removeCost) {
        return true;
    }
    if (removeCost < preserveCost) {
        return false;
    }
    return options.tieBreak == ViterbiTieBreak::PreferPreserve;
}

/**
 * @brief Validates cost shape and values.
 *
 * @param tree Tree topology used by the operation.
 * @param costs Per-node preservation and removal costs.
 */
template <std::floating_point Real> void requireCostShapeAndValues(const MorphologicalTree& tree, const ViterbiNodeCosts<Real>& costs) {
    const auto expectedSize = static_cast<std::size_t>(tree.getNumInternalNodeSlots());
    if (costs.preserve.size() != expectedSize || costs.remove.size() != expectedSize) {
        throw std::invalid_argument("Viterbi cost buffers must match the internal node slot count.");
    }

    for (NodeId nodeId : tree.getAliveNodeIds()) {
        requireFiniteNonNegativeCost(costs.preserve[static_cast<std::size_t>(nodeId)], nodeId, "preserve");
        requireFiniteNonNegativeCost(costs.remove[static_cast<std::size_t>(nodeId)], nodeId, "remove");
    }
}

/**
 * @brief Collects top down order.
 *
 * @param tree Tree topology used by the operation.
 * @return Values produced by the operation.
 */
inline std::vector<NodeId> collectTopDownOrder(const MorphologicalTree& tree) {
    std::vector<NodeId> order;
    order.reserve(static_cast<std::size_t>(tree.getNumNodes()));

    std::vector<NodeId> stack;
    stack.push_back(tree.getRoot());
    while (!stack.empty()) {
        const NodeId nodeId = stack.back();
        stack.pop_back();
        order.push_back(nodeId);

        for (NodeId childNodeId : tree.getChildren(nodeId)) {
            stack.push_back(childNodeId);
        }
    }

    return order;
}

} // namespace viterbi_decision_detail

/**
 * @brief Builds threshold-based Viterbi costs from an increasing node attribute.
 *
 * The cost model is intentionally simple and symmetric around `threshold`:
 *
 * - preserving a node costs `max(0, threshold - attribute[node])`;
 * - removing a node costs `max(0, attribute[node] - threshold)`.
 *
 * Large attribute values therefore favor preservation, small values favor
 * removal, and exact threshold ties are left to `ViterbiDecisionOptions`.
 *
 * @param tree Tree topology used by the operation.
 * @param attribute Attribute requested by the operation.
 * @param threshold Threshold applied by the operation.
 * @return The resulting threshold-based Viterbi costs from an increasing node attribute.
 */
template <std::floating_point Real>
[[nodiscard]] ViterbiNodeCosts<Real> makeThresholdViterbiCosts(const MorphologicalTree& tree, const Real* attribute, Real threshold) {
    if (attribute == nullptr) {
        throw std::invalid_argument("Viterbi threshold costs require a non-null attribute buffer.");
    }
    if (!std::isfinite(threshold)) {
        throw std::invalid_argument("Viterbi threshold must be finite.");
    }

    const auto numNodeSlots = static_cast<std::size_t>(tree.getNumInternalNodeSlots());
    ViterbiNodeCosts<Real> costs{
        std::vector<Real>(numNodeSlots, Real{}),
        std::vector<Real>(numNodeSlots, Real{}),
    };

    for (NodeId nodeId : tree.getAliveNodeIds()) {
        const auto index = static_cast<std::size_t>(nodeId);
        const Real value = attribute[index];
        viterbi_decision_detail::requireFiniteAttribute(value, nodeId);
        costs.preserve[index] = std::max(Real{}, threshold - value);
        costs.remove[index] = std::max(Real{}, value - threshold);
    }

    return costs;
}

/**
 * @brief Computes the optimal connected keep mask by dynamic programming on the tree.
 *
 * The two states are `preserve` and `remove`. The root is forced to
 * `preserve`, and connectivity is enforced by the transition constraint: once a
 * node is removed, all descendants are also removed. If a parent is preserved,
 * each child independently chooses the cheaper of its preserve/remove
 * subproblems. The returned vector is a dense internal-node keep criterion that
 * can be consumed directly by the direct reconstruction rule.
 *
 * @param tree Tree topology used by the operation.
 * @param costs Preserve/remove costs indexed by node.
 * @param options Policy options controlling the operation.
 * @return The computed optimal connected keep mask by dynamic programming on the tree.
 */
template <std::floating_point Real>
[[nodiscard]] std::vector<bool> computeViterbiKeepCriterion(const MorphologicalTree& tree, const ViterbiNodeCosts<Real>& costs,
                                                            const ViterbiDecisionOptions& options = {}) {
    viterbi_decision_detail::requireCostShapeAndValues(tree, costs);

    const auto numNodeSlots = static_cast<std::size_t>(tree.getNumInternalNodeSlots());
    std::vector<Real> preserveCost(numNodeSlots, Real{});
    std::vector<Real> removeCost(numNodeSlots, Real{});
    std::vector<std::uint8_t> childPreservedWhenParentPreserved(numNodeSlots, false);

    const std::vector<NodeId> order = viterbi_decision_detail::collectTopDownOrder(tree);
    for (auto it = order.rbegin(); it != order.rend(); ++it) {
        const NodeId nodeId = *it;
        const auto index = static_cast<std::size_t>(nodeId);
        Real preserve = costs.preserve[index];
        Real remove = costs.remove[index];

        for (NodeId childNodeId : tree.getChildren(nodeId)) {
            const auto childIndex = static_cast<std::size_t>(childNodeId);
            const bool preserveChild = viterbi_decision_detail::choosePreserve(preserveCost[childIndex], removeCost[childIndex], options);
            childPreservedWhenParentPreserved[childIndex] = static_cast<std::uint8_t>(preserveChild);
            preserve = viterbi_decision_detail::addSaturating(preserve, preserveChild ? preserveCost[childIndex] : removeCost[childIndex]);
            remove = viterbi_decision_detail::addSaturating(remove, removeCost[childIndex]);
        }

        preserveCost[index] = preserve;
        removeCost[index] = remove;
    }

    std::vector<bool> keep(numNodeSlots, false);
    std::vector<std::pair<NodeId, bool>> stack;
    stack.emplace_back(tree.getRoot(), true);
    while (!stack.empty()) {
        const auto [nodeId, preserve] = stack.back();
        stack.pop_back();
        keep[static_cast<std::size_t>(nodeId)] = preserve;

        for (NodeId childNodeId : tree.getChildren(nodeId)) {
            const auto childIndex = static_cast<std::size_t>(childNodeId);
            const bool preserveChild = preserve && (childPreservedWhenParentPreserved[childIndex] != 0);
            stack.emplace_back(childNodeId, preserveChild);
        }
    }

    return keep;
}

} // namespace mmcfilters::detail
