#pragma once

#include "../../trees/MorphologicalTree.hpp"
#include "../../utils/Common.hpp"

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
 * @param tree Tree topology used by the operation.
 * @param ascendants Ancestor-node data represented by `ascendants`.
 * @param descendants Descendant-node data represented by `descendants`.
 * @param context Operation context or diagnostic label.
 */
inline void validateStabilityNeighborhoodShape(const MorphologicalTree& tree, const std::vector<NodeId>& ascendants, const std::vector<NodeId>& descendants,
                                               const char* context) {
    const auto expected = static_cast<std::size_t>(tree.getNumInternalNodeSlots());
    if (ascendants.size() != expected || descendants.size() != expected) {
        throw std::invalid_argument(std::string(context) + " neighbourhood size must match the internal node slot count.");
    }
}

/**
 * @brief Tests whether finite variation holds.
 *
 * @param value Value used by the operation.
 * @return True when finite variation; otherwise false.
 */
template <std::floating_point Real> [[nodiscard]] inline bool isFiniteVariation(Real value) { return std::isfinite(static_cast<long double>(value)); }

/**
 * @brief Computes the area variation for one node.
 *
 * Missing neighbours, non-positive center attributes, and non-finite arithmetic
 * produce `NaN`, which downstream selection treats as an invalid variation.
 *
 * @param node Node identifier used by the operation.
 * @param ascendants Ancestor-node data represented by `ascendants`.
 * @param descendants Descendant-node data represented by `descendants`.
 * @param attrAt Attribute information represented by `attrAt`.
 * @return The computed area variation for one node.
 */
template <std::floating_point Real, class AttrGetter>
[[nodiscard]] inline Real computeVariationValue(NodeId node, const std::vector<NodeId>& ascendants, const std::vector<NodeId>& descendants,
                                                AttrGetter& attrAt) {
    const NodeId ascendant = ascendants[static_cast<std::size_t>(node)];
    const NodeId descendant = descendants[static_cast<std::size_t>(node)];
    if (ascendant == InvalidNode || descendant == InvalidNode) {
        return std::numeric_limits<Real>::quiet_NaN();
    }

    const Real center = attrAt(node);
    if (!(center > Real{0}) || !isFiniteVariation(center)) {
        return std::numeric_limits<Real>::quiet_NaN();
    }

    const Real value = (attrAt(ascendant) - attrAt(descendant)) / center;
    return isFiniteVariation(value) ? value : std::numeric_limits<Real>::quiet_NaN();
}

/**
 * @brief Computes node-wise variation from a prepared ancestor/descendant window.
 *
 * @param tree Tree topology used by the operation.
 * @param ascendants Ancestor-node data represented by `ascendants`.
 * @param descendants Descendant-node data represented by `descendants`.
 * @param attrAt Attribute information represented by `attrAt`.
 * @return The computed node-wise variation from a prepared ancestor/descendant window.
 */
template <std::floating_point Real, class AttrGetter>
[[nodiscard]] inline std::vector<Real> computeVariationsFromNeighborhood(const MorphologicalTree& tree, const std::vector<NodeId>& ascendants,
                                                                         const std::vector<NodeId>& descendants, AttrGetter& attrAt) {
    validateStabilityNeighborhoodShape(tree, ascendants, descendants, "computeVariationsFromNeighborhood");
    std::vector<Real> variation(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), std::numeric_limits<Real>::quiet_NaN());

    for (NodeId nodeId : tree.getAliveNodeIds()) {
        variation[static_cast<std::size_t>(nodeId)] = computeVariationValue<Real>(nodeId, ascendants, descendants, attrAt);
    }
    return variation;
}

/**
 * @brief Returns true when a node is a strict local minimum of variation.
 *
 * @param node Node identifier used by the operation.
 * @param variation Per-node variation values.
 * @param ascendants Ancestor-node data represented by `ascendants`.
 * @param descendants Descendant-node data represented by `descendants`.
 * @return True when a node is a strict local minimum of variation.
 */
template <std::floating_point Real>
[[nodiscard]] inline bool isStrictVariationMinimum(NodeId node, const std::vector<Real>& variation, const std::vector<NodeId>& ascendants,
                                                   const std::vector<NodeId>& descendants) {
    const NodeId ascendant = ascendants[static_cast<std::size_t>(node)];
    const NodeId descendant = descendants[static_cast<std::size_t>(node)];
    if (ascendant == InvalidNode || descendant == InvalidNode) {
        return false;
    }

    const Real center = variation[static_cast<std::size_t>(node)];
    const Real ascendantValue = variation[static_cast<std::size_t>(ascendant)];
    const Real descendantValue = variation[static_cast<std::size_t>(descendant)];
    return isFiniteVariation(center) && isFiniteVariation(ascendantValue) && isFiniteVariation(descendantValue) && center < ascendantValue &&
           center < descendantValue;
}

/**
 * @brief Selects nodes that are strict local minima and pass variation/attribute bounds.
 *
 * @param tree Tree topology used by the operation.
 * @param variation Per-node variation values.
 * @param ascendants Ancestor-node data represented by `ascendants`.
 * @param descendants Descendant-node data represented by `descendants`.
 * @param attrAt Attribute information represented by `attrAt`.
 * @param maxVariation Maximum accepted stability variation.
 * @param minAttr Minimum accepted attribute value.
 * @param maxAttr Maximum accepted attribute value.
 * @param count Number represented by `count`.
 * @return The selected nodes that are strict local minima and pass variation/attribute bounds.
 */
template <std::floating_point Real, class AttrGetter>
[[nodiscard]] inline std::vector<uint8_t> selectStrictVariationMinima(const MorphologicalTree& tree, const std::vector<Real>& variation,
                                                                      const std::vector<NodeId>& ascendants, const std::vector<NodeId>& descendants,
                                                                      AttrGetter& attrAt, Real maxVariation, Real minAttr, Real maxAttr, int& count) {
    validateStabilityNeighborhoodShape(tree, ascendants, descendants, "selectStrictVariationMinima");
    if (variation.size() != static_cast<std::size_t>(tree.getNumInternalNodeSlots())) {
        throw std::invalid_argument("selectStrictVariationMinima variation size must match the internal node slot count.");
    }

    count = 0;
    std::vector<uint8_t> selected(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), false);
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        if (!isStrictVariationMinimum(nodeId, variation, ascendants, descendants)) {
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
 * @brief Chooses the finite node with minimum variation in a node/descendant/ascendant window.
 *
 * @param node Node identifier used by the operation.
 * @param variation Per-node variation values.
 * @param ascendants Ancestor-node data represented by `ascendants`.
 * @param descendants Descendant-node data represented by `descendants`.
 * @return The selected finite node with minimum variation in a node/descendant/ascendant window.
 */
template <std::floating_point Real>
[[nodiscard]] inline NodeId nodeWithMinimumVariationInWindow(NodeId node, const std::vector<Real>& variation, const std::vector<NodeId>& ascendants,
                                                             const std::vector<NodeId>& descendants) {
    NodeId bestNode = InvalidNode;
    Real bestValue = std::numeric_limits<Real>::infinity();
    const NodeId candidates[3] = {node, descendants[static_cast<std::size_t>(node)], ascendants[static_cast<std::size_t>(node)]};

    for (const NodeId candidate : candidates) {
        if (candidate == InvalidNode) {
            continue;
        }
        const Real value = variation[static_cast<std::size_t>(candidate)];
        if (!isFiniteVariation(value)) {
            continue;
        }
        if (bestNode == InvalidNode || value < bestValue || (value == bestValue && candidate < bestNode)) {
            bestNode = candidate;
            bestValue = value;
        }
    }

    return bestNode == InvalidNode ? node : bestNode;
}

} // namespace mmcfilters::detail
