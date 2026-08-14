#pragma once

#include "../trees/MorphologicalTree.hpp"
#include "../trees/detail/CommittedTreeAccess.hpp"
#include "../trees/detail/TreeTraversalDetail.hpp"
#include "../utils/Common.hpp"
#include "../utils/Contract.hpp"

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <optional>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mmcfilters::local_attributes {

namespace detail {
struct BinaryVisibilityStateAccess;
}

/** @brief Relative row-column offset of one observation-window sample. */
struct WindowOffset {
    int rowOffset = 0;    ///< Positive values move down.
    int columnOffset = 0; ///< Positive values move right.

    /** @brief Compares both coordinate displacements. @param lhs Left offset. @param rhs Right offset. @return True when both fields match. */
    friend bool operator==(const WindowOffset& lhs, const WindowOffset& rhs) = default;
};

/**
 * @brief Indexed finite set of translated-sample offsets.
 *
 * @details
 * The offset order is part of the contract: offset `j` permanently identifies
 * coordinate `j` of every `BinaryVisibilityState`. A valid window contains
 * `(0, 0)` exactly once, contains no duplicate offsets, and has at most 32
 * coordinates because visibility is stored in a `std::uint32_t` mask.
 */
class ObservationWindow {
  public:
    /// Maximum coordinate count supported by the visibility mask.
    static constexpr std::size_t maxNumOffsets = 32;

    /** @brief Creates and validates an ordered window. @param offsets Ordered distinct offsets containing zero exactly once. */
    explicit ObservationWindow(std::vector<WindowOffset> offsets) : offsets_(std::move(offsets)) { validate(); }

    /** @brief Creates and validates an ordered window. @param offsets Ordered distinct offsets containing zero exactly once. */
    ObservationWindow(std::initializer_list<WindowOffset> offsets) : ObservationWindow(std::vector<WindowOffset>(offsets)) {}

    /** @brief Creates and validates an ordered fixed-size window. @tparam N Number of offsets. @param offsets Ordered offset array. */
    template <std::size_t N>
    explicit ObservationWindow(const std::array<WindowOffset, N>& offsets) : ObservationWindow(std::vector<WindowOffset>(offsets.begin(), offsets.end())) {}

    /** @brief Returns the number of state coordinates. @return Number of offsets. */
    [[nodiscard]] std::size_t size() const noexcept { return offsets_.size(); }
    /** @brief Returns one offset in semantic coordinate order. @param index Coordinate index. @return Offset at `index`. */
    [[nodiscard]] const WindowOffset& operator[](std::size_t index) const noexcept { return offsets_[index]; }
    /** @brief Returns all offsets in semantic coordinate order. @return Read-only offset span. */
    [[nodiscard]] std::span<const WindowOffset> offsets() const noexcept { return offsets_; }
    /** @brief Starts ordered offset iteration. @return Constant iterator to the first offset. */
    [[nodiscard]] auto begin() const noexcept { return offsets_.begin(); }
    /** @brief Ends ordered offset iteration. @return Constant past-the-end iterator. */
    [[nodiscard]] auto end() const noexcept { return offsets_.end(); }

  private:
    std::vector<WindowOffset> offsets_; ///< Offsets in permanent visibility-coordinate order.

    /** @brief Enforces the nonempty, unique, zero-containing, and width contracts. */
    void validate() const {
        if (offsets_.empty()) {
            throw std::invalid_argument("ObservationWindow requires at least one offset.");
        }
        if (offsets_.size() > maxNumOffsets) {
            throw std::invalid_argument("ObservationWindow supports at most 32 offsets.");
        }

        std::size_t numZeroOffsets = 0;
        for (std::size_t i = 0; i < offsets_.size(); ++i) {
            if (offsets_[i] == WindowOffset{}) {
                ++numZeroOffsets;
            }
            for (std::size_t j = 0; j < i; ++j) {
                if (offsets_[i] == offsets_[j]) {
                    throw std::invalid_argument("ObservationWindow does not permit duplicate offsets.");
                }
            }
        }
        if (numZeroOffsets != 1) {
            throw std::invalid_argument("ObservationWindow requires the zero offset exactly once.");
        }
    }
};

/** @brief Binary sample-visibility vector encoded in observation-window order. */
class BinaryVisibilityState {
  public:
    /** @brief Creates a state over a fixed coordinate domain. @param bits Visible coordinates. @param coordinateCount Domain size in `[1, 32]`. */
    BinaryVisibilityState(std::uint32_t bits, std::size_t coordinateCount) : bits_(bits), coordinateCount_(coordinateCount) {
        if (coordinateCount_ == 0) {
            throw std::invalid_argument("BinaryVisibilityState requires at least one coordinate.");
        }
        if (coordinateCount_ > ObservationWindow::maxNumOffsets) {
            throw std::invalid_argument("BinaryVisibilityState supports at most 32 coordinates.");
        }
        const std::uint32_t validMask = coordinateCount_ == 32 ? std::numeric_limits<std::uint32_t>::max()
                                                               : (coordinateCount_ == 0 ? 0 : (std::uint32_t{1} << coordinateCount_) - std::uint32_t{1});
        if ((bits_ & ~validMask) != 0) {
            throw std::invalid_argument("BinaryVisibilityState contains bits outside its coordinate domain.");
        }
    }

    /** @brief Returns the packed visibility bits. @return Visibility mask. */
    [[nodiscard]] std::uint32_t bits() const noexcept { return bits_; }
    /** @brief Returns the state-coordinate count. @return Number of valid bits. */
    [[nodiscard]] std::size_t coordinateCount() const noexcept { return coordinateCount_; }

    /** @brief Tests one visibility coordinate. @param coordinate Coordinate index. @return True when the coordinate is visible. */
    [[nodiscard]] bool isVisible(std::size_t coordinate) const {
        if (coordinate >= coordinateCount_) {
            throw std::out_of_range("BinaryVisibilityState coordinate is outside the observation window.");
        }
        return (bits_ & (std::uint32_t{1} << coordinate)) != 0;
    }

    /** @brief Adds all coordinates entering at one node. @param enteringOffsetMask Coordinates that become visible. @return Updated state. */
    [[nodiscard]] BinaryVisibilityState withEnteringOffsets(std::uint32_t enteringOffsetMask) const {
        return BinaryVisibilityState(bits_ | enteringOffsetMask, coordinateCount_);
    }

    /** @brief Compares state bits and coordinate domains. @param lhs Left state. @param rhs Right state. @return True when both states match. */
    friend bool operator==(const BinaryVisibilityState& lhs, const BinaryVisibilityState& rhs) = default;

  private:
    struct UncheckedConstructionTag {};

    /**
     * @brief Creates a state whose domain and bits were already validated by the finite-window kernel.
     * @param bits Packed visibility values.
     * @param coordinateCount Number of meaningful low bits.
     */
    BinaryVisibilityState(UncheckedConstructionTag, std::uint32_t bits, std::size_t coordinateCount) noexcept
        : bits_(bits), coordinateCount_(coordinateCount) {}

    std::uint32_t bits_ = 0;          ///< Packed visibility values.
    std::size_t coordinateCount_ = 0; ///< Number of meaningful low bits.

    friend struct detail::BinaryVisibilityStateAccess;
};

/** @brief One inclusion node and the observation coordinates that enter there. */
struct AnchoredEntryMask {
    NodeId node = InvalidNode;            ///< Anchored entry shared by the grouped coordinates.
    std::uint32_t enteringOffsetMask = 0; ///< Coordinates that first become visible at `node`.

    /** @brief Compares entry node and coordinate mask. @param lhs Left record. @param rhs Right record. @return True when both fields match. */
    friend bool operator==(const AnchoredEntryMask& lhs, const AnchoredEntryMask& rhs) = default;
};

/** @brief Window-coordinate map whose missing entries represent out-of-domain samples. */
class AnchoredEntryMap {
  public:
    /** @brief Creates a coordinate-preserving map. @param anchorPixel Anchor pixel. @param entries Optional entry for each window coordinate. */
    AnchoredEntryMap(PixelId anchorPixel, std::vector<std::optional<NodeId>> entries) : anchorPixel_(anchorPixel), entries_(std::move(entries)) {}

    /** @brief Returns the anchor pixel. @return Row-major anchor pixel. */
    [[nodiscard]] PixelId anchorPixel() const noexcept { return anchorPixel_; }
    /** @brief Returns the coordinate count. @return Number of mapped window coordinates. */
    [[nodiscard]] std::size_t size() const noexcept { return entries_.size(); }
    /** @brief Returns one optional entry. @param coordinate Window coordinate. @return Entry or no value for an out-of-domain sample. */
    [[nodiscard]] const std::optional<NodeId>& operator[](std::size_t coordinate) const noexcept { return entries_[coordinate]; }
    /** @brief Returns all optional entries in coordinate order. @return Read-only entry span. */
    [[nodiscard]] std::span<const std::optional<NodeId>> entries() const noexcept { return entries_; }

  private:
    PixelId anchorPixel_ = InvalidPixel;          ///< Pixel that fixes the anchor branch.
    std::vector<std::optional<NodeId>> entries_; ///< Optional entries in observation-window order.
};

using AnchorBranch = std::vector<NodeId>;
using AnchoredEntrySet = std::set<NodeId>;
using OrderedAnchoredEntries = std::vector<AnchoredEntryMask>;

/**
 * @brief One anchor-specific local-rule change attached to an anchored entry.
 * @tparam Value Additive-group value type of the local rule.
 */
template <class Value> struct EventDelta {
    PixelId anchorPixel = InvalidPixel;  ///< Anchor whose visibility transition produced the difference.
    NodeId anchoredEntry = InvalidNode; ///< Entry node to which the difference is attached.
    Value value{};                      ///< Signed local-rule difference at the entry.
};

/**
 * @brief Sum of all finite-window event deltas attached to one node.
 * @tparam Value Additive-group value type of the local rule.
 */
template <class Value> struct LocalAttributeIncrement {
    NodeId node = InvalidNode; ///< Dense node slot represented by this record.
    Value value{};             ///< Pre-aggregation contribution introduced at the node.
};

/**
 * @brief Final value of one finite-window node attribute after bottom-up aggregation.
 * @tparam Value Additive-group value type of the local rule.
 */
template <class Value> struct NodeAttribute {
    NodeId node = InvalidNode; ///< Dense node slot represented by this record.
    Value value{};             ///< Fully aggregated node-attribute value.
};

/**
 * @brief Compile-time contract for a rule valued in an additive Abelian group.
 *
 * @details
 * C++ can verify the required operations but not their algebraic laws. A model
 * of this concept must therefore ensure that `additiveIdentity`, `addAssign`,
 * and `subtractAssign` implement the identity, commutative addition, and
 * additive inverse of one Abelian group.
 */
template <class Rule>
concept LocalRule = requires(const Rule& rule, BinaryVisibilityState state, typename Rule::Value& target, const typename Rule::Value& source) {
    typename Rule::Value;
    { rule.additiveIdentity() } -> std::same_as<typename Rule::Value>;
    { rule.evaluateLocalRule(state) } -> std::same_as<typename Rule::Value>;
    { rule.addAssign(target, source) } -> std::same_as<void>;
    { rule.subtractAssign(target, source) } -> std::same_as<void>;
};

namespace detail {

/** @brief Trusted state operations used only after observation-window validation. */
struct BinaryVisibilityStateAccess {
    /**
     * @brief Creates the zero state for an already validated coordinate domain.
     * @param coordinateCount Number of observation-window coordinates.
     * @return Visibility state with every coordinate cleared.
     */
    [[nodiscard]] static BinaryVisibilityState zero(std::size_t coordinateCount) noexcept {
        return BinaryVisibilityState(BinaryVisibilityState::UncheckedConstructionTag{}, 0, coordinateCount);
    }

    /**
     * @brief Activates a mask produced from coordinates of the validated observation window.
     * @param state Visibility state updated in place.
     * @param enteringOffsetMask Coordinates that become visible at the current anchored entry.
     */
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
inline void validateLocalAttributeIncrements(const MorphologicalTree& tree, std::span<const LocalAttributeIncrement<Value>> localAttributeIncrements) {
    if (localAttributeIncrements.size() != static_cast<std::size_t>(tree.numInternalNodeSlots())) {
        throw std::invalid_argument("Local-attribute increments must cover every dense internal node slot.");
    }
    for (std::size_t slot = 0; slot < localAttributeIncrements.size(); ++slot) {
        if (localAttributeIncrements[slot].node != static_cast<NodeId>(slot)) {
            throw std::invalid_argument("Local-attribute increments must be ordered by dense node slot.");
        }
    }
}

namespace kernel {

inline NodeId anchoredEntryFromSmallestNodes(const MorphologicalTree& tree, NodeId anchorSmallestNode, NodeId sampleSmallestNode) {
    if (::mmcfilters::detail::CommittedTreeAccess::isAncestor(tree, anchorSmallestNode, sampleSmallestNode)) {
        return anchorSmallestNode;
    }
    if (::mmcfilters::detail::CommittedTreeAccess::isAncestor(tree, sampleSmallestNode, anchorSmallestNode)) {
        return sampleSmallestNode;
    }
    return ::mmcfilters::detail::CommittedTreeAccess::lowestCommonAncestor(tree, anchorSmallestNode, sampleSmallestNode);
}

inline NodeId anchoredEntry(const MorphologicalTree& tree, PixelId anchorPixel, PixelId samplePixel) {
    const NodeId anchorSmallestNode = ::mmcfilters::detail::CommittedTreeAccess::smallestNodeMap(tree, anchorPixel);
    const NodeId sampleSmallestNode = ::mmcfilters::detail::CommittedTreeAccess::smallestNodeMap(tree, samplePixel);
    return anchoredEntryFromSmallestNodes(tree, anchorSmallestNode, sampleSmallestNode);
}

inline NodeId anchoredEntry(const MorphologicalTree& tree, PixelId anchorPixel, WindowOffset offset) {
    const GridDomain2D& domain = ::mmcfilters::detail::CommittedTreeAccess::gridDomain2D(tree);
    const int anchorRow = anchorPixel / domain.columns;
    const int anchorColumn = anchorPixel % domain.columns;
    const int sampleRow = anchorRow + offset.rowOffset;
    const int sampleColumn = anchorColumn + offset.columnOffset;
    if (sampleRow < 0 || sampleRow >= domain.rows || sampleColumn < 0 || sampleColumn >= domain.columns) {
        return InvalidNode;
    }
    return anchoredEntry(tree, anchorPixel, static_cast<PixelId>(sampleRow * domain.columns + sampleColumn));
}

inline AnchoredEntryMap anchoredEntryMap(const MorphologicalTree& tree, PixelId anchorPixel, const ObservationWindow& observationWindow) {
    std::vector<std::optional<NodeId>> entries;
    entries.reserve(observationWindow.size());
    for (WindowOffset windowOffset : observationWindow) {
        const NodeId entry = anchoredEntry(tree, anchorPixel, windowOffset);
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

/**
 * @brief Builds grouped inclusion-ordered entries in reusable fixed-capacity storage.
 * @param tree Established connected-subset tree.
 * @param domain Established finite grid domain.
 * @param anchorPixel Valid row-major anchor pixel.
 * @param observationWindow Valid finite observation window.
 * @param scratch Caller-owned storage reused across anchors.
 * @return Number of grouped entries written at the beginning of `scratch`.
 */
inline std::size_t fillOrderedAnchoredEntries(const MorphologicalTree& tree, const GridDomain2D& domain, PixelId anchorPixel,
                                              const ObservationWindow& observationWindow, AnchoredEntryScratch& scratch) {
    const int anchorRow = anchorPixel / domain.columns;
    const int anchorColumn = anchorPixel % domain.columns;
    const NodeId anchorSmallestNode = ::mmcfilters::detail::CommittedTreeAccess::smallestNodeMap(tree, anchorPixel);

    std::size_t numEntries = 0;
    for (std::size_t coordinate = 0; coordinate < observationWindow.size(); ++coordinate) {
        const WindowOffset offset = observationWindow[coordinate];
        const int sampleRow = anchorRow + offset.rowOffset;
        const int sampleColumn = anchorColumn + offset.columnOffset;
        if (sampleRow < 0 || sampleRow >= domain.rows || sampleColumn < 0 || sampleColumn >= domain.columns) {
            continue;
        }
        const PixelId samplePixel = static_cast<PixelId>(sampleRow * domain.columns + sampleColumn);
        const NodeId sampleSmallestNode = ::mmcfilters::detail::CommittedTreeAccess::smallestNodeMap(tree, samplePixel);
        const NodeId entry = anchoredEntryFromSmallestNodes(tree, anchorSmallestNode, sampleSmallestNode);
        scratch[numEntries++] = {entry, std::uint32_t{1} << coordinate};
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

template <LocalRule Rule, class Consumer>
inline void visitEventDeltas(std::span<const AnchoredEntryMask> entries, std::size_t coordinateCount, const Rule& localRule, Consumer&& consumer) {
    using Value = typename Rule::Value;
    BinaryVisibilityState visibilityState = BinaryVisibilityStateAccess::zero(coordinateCount);
    Value previousRuleValue = localRule.additiveIdentity();
    bool hasPreviousRuleValue = false;
    for (const AnchoredEntryMask& entry : entries) {
        BinaryVisibilityStateAccess::addEnteringOffsets(visibilityState, entry.enteringOffsetMask);
        Value currentRuleValue = localRule.evaluateLocalRule(visibilityState);
        Value eventDelta = currentRuleValue;
        if (hasPreviousRuleValue) {
            localRule.subtractAssign(eventDelta, previousRuleValue);
        }
        consumer(entry.node, std::move(eventDelta));
        previousRuleValue = std::move(currentRuleValue);
        hasPreviousRuleValue = true;
    }
}

template <LocalRule Rule>
inline std::vector<EventDelta<typename Rule::Value>> computeEventDeltas(const MorphologicalTree& tree, PixelId anchorPixel,
                                                                        const ObservationWindow& observationWindow, const Rule& localRule) {
    using Value = typename Rule::Value;
    const GridDomain2D& domain = ::mmcfilters::detail::CommittedTreeAccess::gridDomain2D(tree);
    AnchoredEntryScratch entryScratch;
    const std::size_t numEntries = fillOrderedAnchoredEntries(tree, domain, anchorPixel, observationWindow, entryScratch);
    const std::span<const AnchoredEntryMask> entries(entryScratch.data(), numEntries);

    std::vector<EventDelta<Value>> eventDeltas;
    eventDeltas.reserve(entries.size());
    visitEventDeltas(entries, observationWindow.size(), localRule,
                     [&](NodeId entry, Value&& eventDelta) { eventDeltas.push_back({anchorPixel, entry, std::move(eventDelta)}); });
    return eventDeltas;
}

template <LocalRule Rule>
inline void accumulateLocalAttributeIncrementValues(const MorphologicalTree& tree, const ObservationWindow& observationWindow, const Rule& localRule,
                                                    std::span<typename Rule::Value> localAttributeIncrementValues) {
    using Value = typename Rule::Value;
    const GridDomain2D& domain = ::mmcfilters::detail::CommittedTreeAccess::gridDomain2D(tree);
    const int totalPixels = domain.rows * domain.columns;
    AnchoredEntryScratch entryScratch;
    for (PixelId anchorPixel = 0; anchorPixel < totalPixels; ++anchorPixel) {
        const std::size_t numEntries = fillOrderedAnchoredEntries(tree, domain, anchorPixel, observationWindow, entryScratch);
        const std::span<const AnchoredEntryMask> entries(entryScratch.data(), numEntries);
        visitEventDeltas(entries, observationWindow.size(), localRule, [&](NodeId entry, Value&& eventDelta) {
            localRule.addAssign(localAttributeIncrementValues[static_cast<std::size_t>(entry)], eventDelta);
        });
    }
}

template <LocalRule Rule>
inline std::vector<typename Rule::Value> computeLocalAttributeIncrementValues(const MorphologicalTree& tree, const ObservationWindow& observationWindow,
                                                                              const Rule& localRule) {
    using Value = typename Rule::Value;
    std::vector<Value> localAttributeIncrementValues;
    localAttributeIncrementValues.reserve(static_cast<std::size_t>(tree.numInternalNodeSlots()));
    for (NodeId node = 0; node < tree.numInternalNodeSlots(); ++node) {
        localAttributeIncrementValues.push_back(localRule.additiveIdentity());
    }
    accumulateLocalAttributeIncrementValues(tree, observationWindow, localRule, localAttributeIncrementValues);
    return localAttributeIncrementValues;
}

template <LocalRule Rule>
inline std::vector<LocalAttributeIncrement<typename Rule::Value>>
computeLocalAttributeIncrements(const MorphologicalTree& tree, const ObservationWindow& observationWindow, const Rule& localRule) {
    using Value = typename Rule::Value;
    std::vector<Value> localAttributeIncrementValues = computeLocalAttributeIncrementValues(tree, observationWindow, localRule);
    std::vector<LocalAttributeIncrement<Value>> localAttributeIncrements;
    localAttributeIncrements.reserve(localAttributeIncrementValues.size());
    for (NodeId node = 0; node < tree.numInternalNodeSlots(); ++node) {
        localAttributeIncrements.push_back({node, std::move(localAttributeIncrementValues[static_cast<std::size_t>(node)])});
    }
    return localAttributeIncrements;
}

template <LocalRule Rule>
inline void aggregateLocalAttributeIncrementValues(const MorphologicalTree& tree, std::span<typename Rule::Value> localAttributeIncrementValues,
                                                   const Rule& localRule) {
    ::mmcfilters::detail::kernel::traversePostOrder(
        tree, tree.root(), [](NodeId) {},
        [&](NodeId parent, NodeId child) {
            localRule.addAssign(localAttributeIncrementValues[static_cast<std::size_t>(parent)],
                                localAttributeIncrementValues[static_cast<std::size_t>(child)]);
        },
        [](NodeId) {});
}

template <LocalRule Rule>
inline std::vector<NodeAttribute<typename Rule::Value>>
aggregateLocalAttributeIncrements(const MorphologicalTree& tree, std::span<const LocalAttributeIncrement<typename Rule::Value>> localAttributeIncrements,
                                  const Rule& localRule) {
    using Value = typename Rule::Value;
    std::vector<NodeAttribute<Value>> nodeAttributes;
    nodeAttributes.reserve(localAttributeIncrements.size());
    for (const LocalAttributeIncrement<Value>& localAttributeIncrement : localAttributeIncrements) {
        nodeAttributes.push_back({localAttributeIncrement.node, localAttributeIncrement.value});
    }

    ::mmcfilters::detail::kernel::traversePostOrder(
        tree, tree.root(), [](NodeId) {},
        [&](NodeId parent, NodeId child) {
            localRule.addAssign(nodeAttributes[static_cast<std::size_t>(parent)].value, nodeAttributes[static_cast<std::size_t>(child)].value);
        },
        [](NodeId) {});
    return nodeAttributes;
}

template <LocalRule Rule>
inline std::vector<NodeAttribute<typename Rule::Value>> computeFiniteWindowLocalAttribute(const MorphologicalTree& tree,
                                                                                          const ObservationWindow& observationWindow, const Rule& localRule) {
    using Value = typename Rule::Value;
    std::vector<Value> localAttributeIncrementValues = computeLocalAttributeIncrementValues(tree, observationWindow, localRule);
    aggregateLocalAttributeIncrementValues(tree, localAttributeIncrementValues, localRule);

    std::vector<NodeAttribute<Value>> nodeAttributes;
    nodeAttributes.reserve(localAttributeIncrementValues.size());
    for (NodeId node = 0; node < tree.numInternalNodeSlots(); ++node) {
        nodeAttributes.push_back({node, std::move(localAttributeIncrementValues[static_cast<std::size_t>(node)])});
    }
    return nodeAttributes;
}

} // namespace kernel
} // namespace detail

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

    /** @brief Computes the rule changes for one anchor. @tparam Rule Model of `LocalRule`. @param tree Tree with a finite 2D domain. @param anchorPixel Anchor
     * pixel. @param observationWindow Valid observation window. @param localRule Rule valued in an additive Abelian group. @return Inclusion-ordered,
     * anchor-specific event deltas. */
    template <LocalRule Rule>
    [[nodiscard]] static std::vector<EventDelta<typename Rule::Value>> computeEventDeltas(const MorphologicalTree& tree, PixelId anchorPixel,
                                                                                          const ObservationWindow& observationWindow, const Rule& localRule) {
        detail::validateFiniteWindowLocalAttributeInput(tree);
        if (anchorPixel < 0 || anchorPixel >= tree.numPixels()) {
            throw std::invalid_argument("Event-delta computation requires a valid anchor pixel.");
        }
        return detail::kernel::computeEventDeltas(tree, anchorPixel, observationWindow, localRule);
    }

    /** @brief Sums all anchor-specific event deltas into dense node increments. @tparam Rule Model of `LocalRule`. @param tree Tree with a finite 2D domain.
     * @param observationWindow Valid observation window. @param localRule Rule valued in an additive Abelian group. @return One local-attribute increment for
     * every dense node slot. */
    template <LocalRule Rule>
    [[nodiscard]] static std::vector<LocalAttributeIncrement<typename Rule::Value>>
    computeLocalAttributeIncrements(const MorphologicalTree& tree, const ObservationWindow& observationWindow, const Rule& localRule) {
        MMCFILTERS_CONTRACT_CHECKED_ONLY(detail::validateFiniteWindowLocalAttributeInput(tree));
        return detail::kernel::computeLocalAttributeIncrements(tree, observationWindow, localRule);
    }

    /** @brief Aggregates node increments from children into parents. @tparam Rule Model of `LocalRule`. @param tree Tree topology. @param
     * localAttributeIncrements Dense node-ordered increments. @param localRule Rule supplying additive accumulation. @return One final node attribute for every
     * dense node slot. */
    template <LocalRule Rule>
    [[nodiscard]] static std::vector<NodeAttribute<typename Rule::Value>>
    aggregateLocalAttributeIncrements(const MorphologicalTree& tree, std::span<const LocalAttributeIncrement<typename Rule::Value>> localAttributeIncrements,
                                      const Rule& localRule) {
        detail::validateLocalAttributeIncrements(tree, localAttributeIncrements);
        return detail::kernel::aggregateLocalAttributeIncrements(tree, localAttributeIncrements, localRule);
    }

    /** @brief Computes the final node attribute induced by a finite window and additive rule. @tparam Rule Model of `LocalRule`. @param tree Tree with a finite
     * 2D domain. @param observationWindow Valid observation window. @param localRule Rule valued in an additive Abelian group. @return Dense node-indexed
     * attribute values after bottom-up aggregation. */
    template <LocalRule Rule>
    [[nodiscard]] static std::vector<NodeAttribute<typename Rule::Value>> compute(const MorphologicalTree& tree, const ObservationWindow& observationWindow,
                                                                                  const Rule& localRule) {
        MMCFILTERS_CONTRACT_CHECKED_ONLY(detail::validateFiniteWindowLocalAttributeInput(tree));
        return detail::kernel::computeFiniteWindowLocalAttribute(tree, observationWindow, localRule);
    }
};

} // namespace mmcfilters::local_attributes
