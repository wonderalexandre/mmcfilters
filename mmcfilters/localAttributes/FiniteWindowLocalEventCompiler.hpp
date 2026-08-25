#pragma once

#include "ConnectedSubsetTreeLocalizer.hpp"
#include "LocalEventModel.hpp"
#include "../trees/MorphologicalTree.hpp"
#include "../trees/detail/CommittedTreeAccess.hpp"
#include "../utils/Contract.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcfilters::local_attributes {

namespace detail {

/** @brief Trusted state operations used only after observation-window validation. */
struct BinaryVisibilityStateAccess {
    /** @brief Builds the trusted all-hidden state. @param coordinateCount Number of coordinates. @return Zero-bit visibility state. */
    [[nodiscard]] static BinaryVisibilityState zero(std::size_t coordinateCount) noexcept {
        return BinaryVisibilityState(BinaryVisibilityState::UncheckedConstructionTag{}, 0, coordinateCount);
    }

    /** @brief Marks coordinates visible without revalidation. @param state State to mutate. @param enteringOffsetMask Entering coordinate bits. */
    static void addEnteringOffsets(BinaryVisibilityState& state, std::uint32_t enteringOffsetMask) noexcept {
        state.bits_ |= enteringOffsetMask;
    }
};

inline void validateAnchoredEntryMap(const MorphologicalTree& tree, const AnchoredEntryMap& entryMap) {
    if (entryMap.size() == 0 || entryMap.size() > ObservationWindow::maxNumOffsets) {
        throw std::invalid_argument("AnchoredEntryMap size must match a valid observation-window state domain.");
    }
    if (entryMap.anchorPixel() < 0 || entryMap.anchorPixel() >= tree.numPixels()) {
        throw std::invalid_argument("AnchoredEntryMap requires a valid anchor pixel.");
    }
    const NodeId anchorSmallestNode = tree.smallestNode(entryMap.anchorPixel());
    if (!tree.isAlive(anchorSmallestNode)) {
        throw std::invalid_argument("AnchoredEntryMap requires an anchor pixel with a live smallest node.");
    }
    for (const std::optional<NodeId>& entry : entryMap.entries()) {
        if (entry.has_value() && (!tree.isAlive(*entry) || !tree.isAncestor(*entry, anchorSmallestNode))) {
            throw std::invalid_argument("Every anchored entry must be a live node on the anchor branch.");
        }
    }
}

inline void validateFiniteWindowLocalAttributeInput(const MorphologicalTree& tree) {
    const GridDomain2D& domain = tree.requireGridDomain2D("FiniteWindowLocalAttributeComputer");
    if (domain.rows <= 0 || domain.columns <= 0) {
        throw std::invalid_argument("FiniteWindowLocalAttributeComputer requires a non-empty 2D image domain.");
    }
    if (!tree.isAlive(tree.root())) {
        throw std::invalid_argument("FiniteWindowLocalAttributeComputer requires a live tree root.");
    }
}

template <class Value>
inline void validateLocalAttributeIncrements(const MorphologicalTree& tree, std::span<const LocalAttributeIncrement<Value>> increments) {
    if (increments.size() != static_cast<std::size_t>(tree.numInternalNodeSlots())) {
        throw std::invalid_argument("Local-attribute increments must cover every dense internal node slot.");
    }
    for (std::size_t slot = 0; slot < increments.size(); ++slot) {
        if (increments[slot].node != static_cast<NodeId>(slot)) {
            throw std::invalid_argument("Local-attribute increments must be ordered by dense node slot.");
        }
    }
}

namespace kernel {

inline AnchoredEntryMap anchoredEntryMap(const MorphologicalTree& tree, PixelId anchorPixel, const ObservationWindow& window) {
    std::vector<std::optional<NodeId>> entries;
    entries.reserve(window.size());
    for (WindowOffset offset : window) {
        const NodeId entry = anchoredEntry(tree, anchorPixel, offset);
        entries.push_back(entry == InvalidNode ? std::nullopt : std::optional<NodeId>{entry});
    }
    return AnchoredEntryMap(anchorPixel, std::move(entries));
}

inline OrderedAnchoredEntries orderedAnchoredEntries(const MorphologicalTree& tree, const AnchoredEntryMap& entryMap) {
    OrderedAnchoredEntries entries;
    entries.reserve(entryMap.size());
    for (std::size_t coordinate = 0; coordinate < entryMap.size(); ++coordinate) {
        if (entryMap[coordinate].has_value()) {
            entries.push_back({*entryMap[coordinate], std::uint32_t{1} << coordinate});
        }
    }

    std::sort(entries.begin(), entries.end(), [&](const AnchoredEntryMask& lhs, const AnchoredEntryMask& rhs) {
        if (lhs.node == rhs.node) {
            return false;
        }
        return ::mmcfilters::detail::CommittedTreeAccess::isAncestor(tree, rhs.node, lhs.node);
    });

    OrderedAnchoredEntries grouped;
    grouped.reserve(entries.size());
    for (const AnchoredEntryMask& entry : entries) {
        if (!grouped.empty() && grouped.back().node == entry.node) {
            grouped.back().enteringOffsetMask |= entry.enteringOffsetMask;
        } else {
            grouped.push_back(entry);
        }
    }
    return grouped;
}

using AnchoredEntryScratch = std::array<AnchoredEntryMask, ObservationWindow::maxNumOffsets>;

inline std::size_t fillOrderedAnchoredEntries(const MorphologicalTree& tree, const GridDomain2D& domain, PixelId anchorPixel,
                                              const ObservationWindow& window, AnchoredEntryScratch& scratch) {
    const int anchorRow = anchorPixel / domain.columns;
    const int anchorColumn = anchorPixel % domain.columns;
    const NodeId anchorSmallestNode = ::mmcfilters::detail::CommittedTreeAccess::smallestNodeMap(tree, anchorPixel);

    std::size_t numEntries = 0;
    for (std::size_t coordinate = 0; coordinate < window.size(); ++coordinate) {
        const WindowOffset offset = window[coordinate];
        const int sampleRow = anchorRow + offset.rowOffset;
        const int sampleColumn = anchorColumn + offset.columnOffset;
        if (sampleRow < 0 || sampleRow >= domain.rows || sampleColumn < 0 || sampleColumn >= domain.columns) {
            continue;
        }
        const PixelId samplePixel = static_cast<PixelId>(sampleRow * domain.columns + sampleColumn);
        const NodeId sampleSmallestNode = ::mmcfilters::detail::CommittedTreeAccess::smallestNodeMap(tree, samplePixel);
        scratch[numEntries++] = {connectedSubsetJoin(tree, anchorSmallestNode, sampleSmallestNode), std::uint32_t{1} << coordinate};
    }

    std::sort(scratch.begin(), scratch.begin() + static_cast<std::ptrdiff_t>(numEntries), [&](const AnchoredEntryMask& lhs, const AnchoredEntryMask& rhs) {
        if (lhs.node == rhs.node) {
            return false;
        }
        return ::mmcfilters::detail::CommittedTreeAccess::isAncestor(tree, rhs.node, lhs.node);
    });

    std::size_t numGroupedEntries = 0;
    for (std::size_t index = 0; index < numEntries; ++index) {
        const AnchoredEntryMask entry = scratch[index];
        if (numGroupedEntries > 0 && scratch[numGroupedEntries - 1].node == entry.node) {
            scratch[numGroupedEntries - 1].enteringOffsetMask |= entry.enteringOffsetMask;
        } else {
            scratch[numGroupedEntries++] = entry;
        }
    }
    return numGroupedEntries;
}

template <LocalDecision Decision, class Algebra, class Consumer>
    requires EventAlgebra<Algebra, typename Decision::Value>
inline void visitEventDeltas(std::span<const AnchoredEntryMask> entries, std::size_t coordinateCount, const Decision& decision, const Algebra& algebra,
                             Consumer&& consumer) {
    using Value = typename Decision::Value;
    BinaryVisibilityState state = BinaryVisibilityStateAccess::zero(coordinateCount);
    Value previousValue = algebra.additiveIdentity();
    bool hasPreviousValue = false;
    for (const AnchoredEntryMask& entry : entries) {
        BinaryVisibilityStateAccess::addEnteringOffsets(state, entry.enteringOffsetMask);
        Value currentValue = decision.evaluateLocalDecision(state);
        Value delta = currentValue;
        if (hasPreviousValue) {
            algebra.subtractAssign(delta, previousValue);
        }
        consumer(entry.node, std::move(delta));
        previousValue = std::move(currentValue);
        hasPreviousValue = true;
    }
}

template <LocalDecision Decision, class Algebra>
    requires EventAlgebra<Algebra, typename Decision::Value>
inline std::vector<EventDelta<typename Decision::Value>> computeEventDeltas(const MorphologicalTree& tree, PixelId anchorPixel,
                                                                            const ObservationWindow& window, const Decision& decision,
                                                                            const Algebra& algebra) {
    using Value = typename Decision::Value;
    const GridDomain2D& domain = ::mmcfilters::detail::CommittedTreeAccess::gridDomain2D(tree);
    AnchoredEntryScratch scratch;
    const std::size_t numEntries = fillOrderedAnchoredEntries(tree, domain, anchorPixel, window, scratch);
    const std::span<const AnchoredEntryMask> entries(scratch.data(), numEntries);

    std::vector<EventDelta<Value>> deltas;
    deltas.reserve(entries.size());
    visitEventDeltas(entries, window.size(), decision, algebra,
                     [&](NodeId entry, Value&& delta) { deltas.push_back({anchorPixel, entry, std::move(delta)}); });
    return deltas;
}

template <LocalDecision Decision, class Algebra>
    requires EventAlgebra<Algebra, typename Decision::Value>
inline void accumulateLocalAttributeIncrementValues(const MorphologicalTree& tree, const ObservationWindow& window, const Decision& decision,
                                                    const Algebra& algebra, std::span<typename Decision::Value> increments) {
    using Value = typename Decision::Value;
    const GridDomain2D& domain = ::mmcfilters::detail::CommittedTreeAccess::gridDomain2D(tree);
    const int totalPixels = domain.rows * domain.columns;
    AnchoredEntryScratch scratch;
    for (PixelId anchorPixel = 0; anchorPixel < totalPixels; ++anchorPixel) {
        const std::size_t numEntries = fillOrderedAnchoredEntries(tree, domain, anchorPixel, window, scratch);
        const std::span<const AnchoredEntryMask> entries(scratch.data(), numEntries);
        visitEventDeltas(entries, window.size(), decision, algebra,
                         [&](NodeId entry, Value&& delta) { algebra.addAssign(increments[static_cast<std::size_t>(entry)], delta); });
    }
}

template <LocalDecision Decision, class Algebra>
    requires EventAlgebra<Algebra, typename Decision::Value>
inline std::vector<typename Decision::Value> computeLocalAttributeIncrementValues(const MorphologicalTree& tree, const ObservationWindow& window,
                                                                                  const Decision& decision, const Algebra& algebra) {
    using Value = typename Decision::Value;
    std::vector<Value> increments;
    increments.reserve(static_cast<std::size_t>(tree.numInternalNodeSlots()));
    for (NodeId node = 0; node < tree.numInternalNodeSlots(); ++node) {
        increments.push_back(algebra.additiveIdentity());
    }
    accumulateLocalAttributeIncrementValues(tree, window, decision, algebra, increments);
    return increments;
}

template <LocalDecision Decision, class Algebra>
    requires EventAlgebra<Algebra, typename Decision::Value>
inline std::vector<LocalAttributeIncrement<typename Decision::Value>> computeLocalAttributeIncrements(const MorphologicalTree& tree,
                                                                                                      const ObservationWindow& window,
                                                                                                      const Decision& decision,
                                                                                                      const Algebra& algebra) {
    using Value = typename Decision::Value;
    std::vector<Value> values = computeLocalAttributeIncrementValues(tree, window, decision, algebra);
    std::vector<LocalAttributeIncrement<Value>> increments;
    increments.reserve(values.size());
    for (NodeId node = 0; node < tree.numInternalNodeSlots(); ++node) {
        increments.push_back({node, std::move(values[static_cast<std::size_t>(node)])});
    }
    return increments;
}

} // namespace kernel
} // namespace detail

/** @brief Compiles pure finite-window decisions into hierarchy-attached event deltas. */
class FiniteWindowLocalEventCompiler {
  public:
    /**
     * @brief Compiles the state changes of one anchor into sparse node events.
     * @tparam Decision Model of `LocalDecision`.
     * @tparam Algebra Model of `EventAlgebra` for the decision value.
     * @param tree Stable tree with a finite 2D domain.
     * @param anchorPixel Valid anchor pixel.
     * @param window Observation window.
     * @param decision Pure state-to-value decision.
     * @param algebra Additive difference algebra.
     * @return Inclusion-ordered event deltas for the anchor.
     */
    template <LocalDecision Decision, class Algebra>
        requires EventAlgebra<Algebra, typename Decision::Value>
    [[nodiscard]] static std::vector<EventDelta<typename Decision::Value>> computeEventDeltas(const MorphologicalTree& tree, PixelId anchorPixel,
                                                                                              const ObservationWindow& window,
                                                                                              const Decision& decision, const Algebra& algebra) {
        detail::validateFiniteWindowLocalAttributeInput(tree);
        if (anchorPixel < 0 || anchorPixel >= tree.numPixels()) {
            throw std::invalid_argument("Event-delta computation requires a valid anchor pixel.");
        }
        return detail::kernel::computeEventDeltas(tree, anchorPixel, window, decision, algebra);
    }

    /**
     * @brief Compiles every anchor into dense node-local increments.
     * @tparam Decision Model of `LocalDecision`.
     * @tparam Algebra Model of `EventAlgebra` for the decision value.
     * @param tree Stable tree with a finite 2D domain.
     * @param window Observation window.
     * @param decision Pure state-to-value decision.
     * @param algebra Additive difference algebra.
     * @return One non-aggregated increment per dense node slot.
     */
    template <LocalDecision Decision, class Algebra>
        requires EventAlgebra<Algebra, typename Decision::Value>
    [[nodiscard]] static std::vector<LocalAttributeIncrement<typename Decision::Value>>
    computeLocalAttributeIncrements(const MorphologicalTree& tree, const ObservationWindow& window, const Decision& decision, const Algebra& algebra) {
        MMCFILTERS_CONTRACT_CHECKED_ONLY(detail::validateFiniteWindowLocalAttributeInput(tree));
        return detail::kernel::computeLocalAttributeIncrements(tree, window, decision, algebra);
    }
};

} // namespace mmcfilters::local_attributes
