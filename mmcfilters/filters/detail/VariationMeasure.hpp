#pragma once

#include "../../trees/MorphologicalTree.hpp"
#include "../../utils/Common.hpp"
#include "../../utils/Contract.hpp"

#include <cmath>
#include <concepts>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace mmcfilters::detail {

/**
 * @brief Checks that a stability neighbourhood covers the dense node slot domain.
 *
 * @param tree Tree topology.
 * @param ancestors Ancestor-node data.
 * @param descendants Descendant-node data.
 * @param context Operation context or diagnostic label.
 */
inline void validateStabilityNeighborhoodShape(const MorphologicalTree& tree, const std::vector<NodeId>& ancestors, const std::vector<NodeId>& descendants,
                                               const char* context) {
    const auto expected = static_cast<std::size_t>(tree.numInternalNodeSlots());
    MMCFILTERS_CONTRACT_REQUIRE(ancestors.size() == expected && descendants.size() == expected,
                                throw std::invalid_argument(std::string(context) + " neighbourhood size must match the internal node slot count."));
}

/**
 * @brief Tests whether finite variation holds.
 *
 * @param value Value.
 * @return True when finite variation; otherwise false.
 */
template <std::floating_point Real> [[nodiscard]] inline bool isFiniteVariation(Real value) { return std::isfinite(static_cast<long double>(value)); }

/**
 * @brief Computes the area variation for one node.
 *
 * Missing neighbours, non-positive center attributes, and non-finite arithmetic
 * produce `NaN`, which downstream selection treats as an invalid variation.
 *
 * @param node Node identifier.
 * @param ancestors Ancestor-node data.
 * @param descendants Descendant-node data.
 * @param attrAt Attribute information.
 * @return The computed area variation for one node.
 */
template <std::floating_point Real, class AttrGetter>
[[nodiscard]] inline Real computeVariationValue(NodeId node, const std::vector<NodeId>& ancestors, const std::vector<NodeId>& descendants,
                                                AttrGetter& attrAt) {
    const NodeId ancestor = ancestors[static_cast<std::size_t>(node)];
    const NodeId descendant = descendants[static_cast<std::size_t>(node)];
    if (ancestor == InvalidNode || descendant == InvalidNode) {
        return std::numeric_limits<Real>::quiet_NaN();
    }

    const Real center = attrAt(node);
    if (!(center > Real{0}) || !isFiniteVariation(center)) {
        return std::numeric_limits<Real>::quiet_NaN();
    }

    const Real value = (attrAt(ancestor) - attrAt(descendant)) / center;
    return isFiniteVariation(value) ? value : std::numeric_limits<Real>::quiet_NaN();
}

/**
 * @brief Computes node-wise variation from a prepared ancestor/descendant window.
 *
 * @param tree Tree topology.
 * @param ancestors Ancestor-node data.
 * @param descendants Descendant-node data.
 * @param attrAt Attribute information.
 * @return The computed node-wise variation from a prepared ancestor/descendant window.
 */
template <std::floating_point Real, class AttrGetter>
[[nodiscard]] inline std::vector<Real> computeVariationsFromNeighborhood(const MorphologicalTree& tree, const std::vector<NodeId>& ancestors,
                                                                         const std::vector<NodeId>& descendants, AttrGetter& attrAt) {
    validateStabilityNeighborhoodShape(tree, ancestors, descendants, "computeVariationsFromNeighborhood");
    std::vector<Real> variation(static_cast<std::size_t>(tree.numInternalNodeSlots()), std::numeric_limits<Real>::quiet_NaN());

    for (NodeId nodeId : tree.aliveNodeIds()) {
        variation[static_cast<std::size_t>(nodeId)] = computeVariationValue<Real>(nodeId, ancestors, descendants, attrAt);
    }
    return variation;
}

/**
 * @brief Returns true when a node is a strict local minimum of variation.
 *
 * @param node Node identifier.
 * @param variation Per-node variation values.
 * @param ancestors Ancestor-node data.
 * @param descendants Descendant-node data.
 * @return True when a node is a strict local minimum of variation.
 */
template <std::floating_point Real>
[[nodiscard]] inline bool isStrictVariationMinimum(NodeId node, const std::vector<Real>& variation, const std::vector<NodeId>& ancestors,
                                                   const std::vector<NodeId>& descendants) {
    const NodeId ancestor = ancestors[static_cast<std::size_t>(node)];
    const NodeId descendant = descendants[static_cast<std::size_t>(node)];
    if (ancestor == InvalidNode || descendant == InvalidNode) {
        return false;
    }

    const Real center = variation[static_cast<std::size_t>(node)];
    const Real ancestorValue = variation[static_cast<std::size_t>(ancestor)];
    const Real descendantValue = variation[static_cast<std::size_t>(descendant)];
    return isFiniteVariation(center) && isFiniteVariation(ancestorValue) && isFiniteVariation(descendantValue) && center < ancestorValue &&
           center < descendantValue;
}

/**
 * @brief Selects nodes that are strict local minima and pass variation/attribute bounds.
 *
 * @param tree Tree topology.
 * @param variation Per-node variation values.
 * @param ancestors Ancestor-node data.
 * @param descendants Descendant-node data.
 * @param attrAt Attribute information.
 * @param maxVariation Maximum accepted stability variation.
 * @param minAttr Minimum accepted attribute value.
 * @param maxAttr Maximum accepted attribute value.
 * @param count Number.
 * @return The selected nodes that are strict local minima and pass variation/attribute bounds.
 */
template <std::floating_point Real, class AttrGetter>
[[nodiscard]] inline std::vector<uint8_t> selectStrictVariationMinima(const MorphologicalTree& tree, const std::vector<Real>& variation,
                                                                      const std::vector<NodeId>& ancestors, const std::vector<NodeId>& descendants,
                                                                      AttrGetter& attrAt, Real maxVariation, Real minAttr, Real maxAttr, int& count) {
    validateStabilityNeighborhoodShape(tree, ancestors, descendants, "selectStrictVariationMinima");
    MMCFILTERS_CONTRACT_REQUIRE(variation.size() == static_cast<std::size_t>(tree.numInternalNodeSlots()),
                                throw std::invalid_argument("selectStrictVariationMinima variation size must match the internal node slot count."));

    count = 0;
    std::vector<uint8_t> selected(static_cast<std::size_t>(tree.numInternalNodeSlots()), false);
    for (NodeId nodeId : tree.aliveNodeIds()) {
        if (!isStrictVariationMinimum(nodeId, variation, ancestors, descendants)) {
            continue;
        }
        const Real attribute = attrAt(nodeId);
        if (variation[static_cast<std::size_t>(nodeId)] < maxVariation && attribute >= minAttr && attribute <= maxAttr) {
            selected[static_cast<std::size_t>(nodeId)] = true;
            ++count;
        }
    }
    return selected;
}

/**
 * @brief Chooses the finite node with minimum variation in a node/descendant/ancestor window.
 *
 * Equal finite variations preserve the semantic candidate order: current node,
 * representative descendant, then ancestor. This keeps the decision independent
 * of the internal node-identifier assignment and avoids relocating a rejection
 * when the current node is equally stable.
 *
 * @param node Node identifier.
 * @param variation Per-node variation values.
 * @param ancestors Ancestor-node data.
 * @param descendants Descendant-node data.
 * @return The selected finite node with minimum variation in a node/descendant/ancestor window.
 */
template <std::floating_point Real>
[[nodiscard]] inline NodeId nodeWithMinimumVariationInWindow(NodeId node, const std::vector<Real>& variation, const std::vector<NodeId>& ancestors,
                                                             const std::vector<NodeId>& descendants) {
    NodeId bestNode = InvalidNode;
    Real bestValue = std::numeric_limits<Real>::infinity();
    const NodeId candidates[3] = {node, descendants[static_cast<std::size_t>(node)], ancestors[static_cast<std::size_t>(node)]};

    for (const NodeId candidate : candidates) {
        if (candidate == InvalidNode) {
            continue;
        }
        const Real value = variation[static_cast<std::size_t>(candidate)];
        if (!isFiniteVariation(value)) {
            continue;
        }
        if (bestNode == InvalidNode || value < bestValue) {
            bestNode = candidate;
            bestValue = value;
        }
    }

    return bestNode == InvalidNode ? node : bestNode;
}

} // namespace mmcfilters::detail
