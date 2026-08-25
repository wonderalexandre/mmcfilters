#pragma once

#include "ConnectedSubsetTreeLocalizer.hpp"
#include "FiniteWindowLocalEventCompiler.hpp"
#include "LocalEventModel.hpp"
#include "../trees/MorphologicalTree.hpp"
#include "../trees/detail/TreeTraversalDetail.hpp"
#include "../utils/Common.hpp"
#include "../utils/Contract.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mmcfilters::local_attributes {

namespace detail::kernel {

template <class Value, class Algebra>
    requires EventAlgebra<Algebra, Value>
inline void aggregateLocalAttributeIncrementValues(const MorphologicalTree& tree, std::span<Value> increments, const Algebra& algebra) {
    ::mmcfilters::detail::kernel::traversePostOrder(
        tree, tree.root(), [](NodeId) {},
        [&](NodeId parent, NodeId child) { algebra.addAssign(increments[static_cast<std::size_t>(parent)], increments[static_cast<std::size_t>(child)]); },
        [](NodeId) {});
}

template <class Value, class Algebra>
    requires EventAlgebra<Algebra, Value>
inline std::vector<NodeAttribute<Value>> aggregateLocalAttributeIncrements(const MorphologicalTree& tree,
                                                                          std::span<const LocalAttributeIncrement<Value>> increments,
                                                                          const Algebra& algebra) {
    std::vector<NodeAttribute<Value>> attributes;
    attributes.reserve(increments.size());
    for (const LocalAttributeIncrement<Value>& increment : increments) {
        attributes.push_back({increment.node, increment.value});
    }

    ::mmcfilters::detail::kernel::traversePostOrder(
        tree, tree.root(), [](NodeId) {},
        [&](NodeId parent, NodeId child) {
            algebra.addAssign(attributes[static_cast<std::size_t>(parent)].value, attributes[static_cast<std::size_t>(child)].value);
        },
        [](NodeId) {});
    return attributes;
}

template <LocalDecision Decision, class Algebra>
    requires EventAlgebra<Algebra, typename Decision::Value>
inline std::vector<NodeAttribute<typename Decision::Value>> computeFiniteWindowLocalAttribute(const MorphologicalTree& tree,
                                                                                              const ObservationWindow& window,
                                                                                              const Decision& decision, const Algebra& algebra) {
    using Value = typename Decision::Value;
    std::vector<Value> values = computeLocalAttributeIncrementValues(tree, window, decision, algebra);
    aggregateLocalAttributeIncrementValues(tree, std::span<Value>(values), algebra);

    std::vector<NodeAttribute<Value>> attributes;
    attributes.reserve(values.size());
    for (NodeId node = 0; node < tree.numInternalNodeSlots(); ++node) {
        attributes.push_back({node, std::move(values[static_cast<std::size_t>(node)])});
    }
    return attributes;
}

} // namespace detail::kernel

/** @brief Generic finite-window local-attribute computation over a connected-subset tree model. */
class FiniteWindowLocalAttributeComputer {
  public:
    /** @brief Locates a valid absolute sample on the anchor branch. @param tree Tree topology. @param anchorPixel Anchor pixel. @param samplePixel Sample
     * pixel. @return Anchored entry or no value for an invalid/unmapped pixel. */
    [[nodiscard]] static std::optional<NodeId> anchoredEntry(const MorphologicalTree& tree, PixelId anchorPixel, PixelId samplePixel) {
        if (anchorPixel < 0 || anchorPixel >= tree.numPixels() || samplePixel < 0 || samplePixel >= tree.numPixels()) {
            return std::nullopt;
        }
        const NodeId anchorSmallestNode = tree.smallestNode(anchorPixel);
        const NodeId sampleSmallestNode = tree.smallestNode(samplePixel);
        if (!tree.isAlive(anchorSmallestNode) || !tree.isAlive(sampleSmallestNode)) {
            return std::nullopt;
        }
        const NodeId entry = detail::kernel::anchoredEntry(tree, anchorPixel, samplePixel);
        return entry == InvalidNode ? std::nullopt : std::optional<NodeId>{entry};
    }

    /** @brief Locates one translated sample on the anchor branch. @param tree Tree with a 2D domain. @param anchorPixel Anchor pixel. @param windowOffset
     * Translated-sample offset. @return Anchored entry or no value outside the domain. */
    [[nodiscard]] static std::optional<NodeId> anchoredEntry(const MorphologicalTree& tree, PixelId anchorPixel, WindowOffset windowOffset) {
        const GridDomain2D& domain = tree.requireGridDomain2D("FiniteWindowLocalAttributeComputer::anchoredEntry");
        if (anchorPixel < 0 || anchorPixel >= domain.rows * domain.columns) {
            return std::nullopt;
        }
        const NodeId entry = detail::kernel::anchoredEntry(tree, anchorPixel, windowOffset);
        return entry == InvalidNode ? std::nullopt : std::optional<NodeId>{entry};
    }

    /** @brief Materializes the branch from the smallest node to the root. @param tree Tree topology. @param anchorPixel Anchor pixel. @return
     * Increasing-inclusion branch, or an empty vector for an invalid anchor. */
    [[nodiscard]] static AnchorBranch anchorBranch(const MorphologicalTree& tree, PixelId anchorPixel) {
        if (anchorPixel < 0 || anchorPixel >= tree.numPixels()) {
            return {};
        }
        NodeId node = tree.smallestNode(anchorPixel);
        if (!tree.isAlive(node)) {
            return {};
        }

        AnchorBranch branch;
        while (true) {
            branch.push_back(node);
            if (node == tree.root()) {
                break;
            }
            node = tree.parent(node);
        }
        return branch;
    }

    /** @brief Materializes anchored entries in window-coordinate order. @param tree Tree with a 2D domain. @param anchorPixel Anchor pixel. @param
     * observationWindow Valid observation window. @return Coordinate-preserving entry map. */
    [[nodiscard]] static AnchoredEntryMap anchoredEntryMap(const MorphologicalTree& tree, PixelId anchorPixel, const ObservationWindow& observationWindow) {
        detail::validateFiniteWindowLocalAttributeInput(tree);
        if (anchorPixel < 0 || anchorPixel >= tree.numPixels()) {
            return AnchoredEntryMap(anchorPixel, std::vector<std::optional<NodeId>>(observationWindow.size()));
        }
        return detail::kernel::anchoredEntryMap(tree, anchorPixel, observationWindow);
    }

    /** @brief Materializes the distinct anchored entries. @param tree Tree with a 2D domain. @param anchorPixel Anchor pixel. @param observationWindow Valid
     * observation window. @return Mathematical anchored-entry set. */
    [[nodiscard]] static AnchoredEntrySet anchoredEntrySet(const MorphologicalTree& tree, PixelId anchorPixel, const ObservationWindow& observationWindow) {
        AnchoredEntrySet entries;
        const AnchoredEntryMap entryMap = anchoredEntryMap(tree, anchorPixel, observationWindow);
        for (const std::optional<NodeId>& entry : entryMap.entries()) {
            if (entry.has_value()) {
                entries.insert(*entry);
            }
        }
        return entries;
    }

    /** @brief Groups and orders entries from the smallest node toward the root. @param tree Tree with a 2D domain. @param anchorPixel Anchor pixel. @param
     * observationWindow Valid observation window. @return Inclusion-ordered entry masks. */
    [[nodiscard]] static OrderedAnchoredEntries orderedAnchoredEntries(const MorphologicalTree& tree, PixelId anchorPixel,
                                                                       const ObservationWindow& observationWindow) {
        const AnchoredEntryMap entryMap = anchoredEntryMap(tree, anchorPixel, observationWindow);
        return detail::kernel::orderedAnchoredEntries(tree, entryMap);
    }

    /** @brief Groups and orders an existing entry map by increasing inclusion. @param tree Tree topology. @param entryMap Valid map on one anchor branch.
     * @return Inclusion-ordered entry masks. */
    [[nodiscard]] static OrderedAnchoredEntries orderAnchoredEntriesByInclusion(const MorphologicalTree& tree, const AnchoredEntryMap& entryMap) {
        detail::validateAnchoredEntryMap(tree, entryMap);
        return detail::kernel::orderedAnchoredEntries(tree, entryMap);
    }

    /** @brief Evaluates visibility at one node of the anchor branch. @param tree Tree topology. @param entryMap Valid anchored-entry map. @param node Node on
     * the map's anchor branch. @return Binary visibility state in window-coordinate order. */
    [[nodiscard]] static BinaryVisibilityState binaryVisibilityState(const MorphologicalTree& tree, const AnchoredEntryMap& entryMap, NodeId node) {
        std::uint32_t stateBits = 0;
        if (entryMap.size() == 0) {
            throw std::invalid_argument("Binary visibility state requires a non-empty anchored-entry map.");
        }
        detail::validateAnchoredEntryMap(tree, entryMap);
        if (!tree.isAlive(node)) {
            throw std::invalid_argument("Binary visibility state requires a live node.");
        }
        if (entryMap.anchorPixel() < 0 || entryMap.anchorPixel() >= tree.numPixels()) {
            throw std::invalid_argument("Binary visibility state requires a valid anchor pixel.");
        }
        const NodeId anchorSmallestNode = tree.smallestNode(entryMap.anchorPixel());
        if (!tree.isAlive(anchorSmallestNode) || !tree.isAncestor(node, anchorSmallestNode)) {
            throw std::invalid_argument("Binary visibility state is defined only on the anchor branch.");
        }
        for (std::size_t coordinate = 0; coordinate < entryMap.size(); ++coordinate) {
            if (entryMap[coordinate].has_value() && tree.isAncestor(node, *entryMap[coordinate])) {
                stateBits |= std::uint32_t{1} << coordinate;
            }
        }
        return BinaryVisibilityState(stateBits, entryMap.size());
    }

    /**
     * @brief Computes event deltas from a decision and a separately supplied additive algebra.
     * @tparam Decision Model of `LocalDecision`.
     * @tparam Algebra Model of `EventAlgebra` for the decision value.
     * @param tree Stable tree with a finite 2D domain.
     * @param anchorPixel Valid anchor pixel.
     * @param observationWindow Observation window.
     * @param decision Pure state-to-value decision.
     * @param algebra Additive difference algebra.
     * @return Inclusion-ordered event deltas for the anchor.
     */
    template <LocalDecision Decision, class Algebra>
        requires EventAlgebra<Algebra, typename Decision::Value>
    [[nodiscard]] static std::vector<EventDelta<typename Decision::Value>> computeEventDeltas(const MorphologicalTree& tree, PixelId anchorPixel,
                                                                                              const ObservationWindow& observationWindow,
                                                                                              const Decision& decision, const Algebra& algebra) {
        return FiniteWindowLocalEventCompiler::computeEventDeltas(tree, anchorPixel, observationWindow, decision, algebra);
    }

    /**
     * @brief Computes dense increments from a decision and an independent additive algebra.
     * @tparam Decision Model of `LocalDecision`.
     * @tparam Algebra Model of `EventAlgebra` for the decision value.
     * @param tree Stable tree with a finite 2D domain.
     * @param observationWindow Observation window.
     * @param decision Pure state-to-value decision.
     * @param algebra Additive difference algebra.
     * @return One non-aggregated increment per dense node slot.
     */
    template <LocalDecision Decision, class Algebra>
        requires EventAlgebra<Algebra, typename Decision::Value>
    [[nodiscard]] static std::vector<LocalAttributeIncrement<typename Decision::Value>>
    computeLocalAttributeIncrements(const MorphologicalTree& tree, const ObservationWindow& observationWindow, const Decision& decision,
                                    const Algebra& algebra) {
        return FiniteWindowLocalEventCompiler::computeLocalAttributeIncrements(tree, observationWindow, decision, algebra);
    }

    /**
     * @brief Aggregates previously compiled increments with an independent event algebra.
     * @tparam Value Increment value type.
     * @tparam Algebra Model of `EventAlgebra` for `Value`.
     * @param tree Stable tree topology.
     * @param increments Dense node-ordered increments.
     * @param algebra Additive aggregation algebra.
     * @return One aggregated attribute per dense node slot.
     */
    template <class Value, class Algebra>
        requires EventAlgebra<Algebra, Value>
    [[nodiscard]] static std::vector<NodeAttribute<Value>> aggregateEventIncrements(const MorphologicalTree& tree,
                                                                                   std::span<const LocalAttributeIncrement<Value>> increments,
                                                                                   const Algebra& algebra) {
        detail::validateLocalAttributeIncrements(tree, increments);
        return detail::kernel::aggregateLocalAttributeIncrements(tree, increments, algebra);
    }

    /**
     * @brief Computes a finite-window attribute from independent decision and algebra policies.
     * @tparam Decision Model of `LocalDecision`.
     * @tparam Algebra Model of `EventAlgebra` for the decision value.
     * @param tree Stable tree with a finite 2D domain.
     * @param observationWindow Observation window.
     * @param decision Pure state-to-value decision.
     * @param algebra Additive event and aggregation algebra.
     * @return One aggregated attribute per dense node slot.
     */
    template <LocalDecision Decision, class Algebra>
        requires EventAlgebra<Algebra, typename Decision::Value>
    [[nodiscard]] static std::vector<NodeAttribute<typename Decision::Value>> compute(const MorphologicalTree& tree,
                                                                                      const ObservationWindow& observationWindow,
                                                                                      const Decision& decision, const Algebra& algebra) {
        MMCFILTERS_CONTRACT_CHECKED_ONLY(detail::validateFiniteWindowLocalAttributeInput(tree));
        return detail::kernel::computeFiniteWindowLocalAttribute(tree, observationWindow, decision, algebra);
    }
};

} // namespace mmcfilters::local_attributes
