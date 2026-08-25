#pragma once

#include "ConnectedSubsetTreeLocalizer.hpp"
#include "../utils/Common.hpp"

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
#include <utility>
#include <vector>

namespace mmcfilters::local_attributes {

namespace detail {
struct BinaryVisibilityStateAccess;
}

/**
 * @brief Indexed finite set of translated-sample offsets.
 *
 * The offset order permanently identifies visibility-state coordinates.  A
 * valid window contains `(0, 0)` exactly once, has no duplicates, and contains
 * at most 32 offsets because the production state is a `std::uint32_t` mask.
 */
class ObservationWindow {
  public:
    /// Maximum number of offsets representable by one visibility bit mask.
    static constexpr std::size_t maxNumOffsets = 32;

    /** @brief Builds and validates a window. @param offsets Ordered translated-sample offsets. */
    explicit ObservationWindow(std::vector<WindowOffset> offsets) : offsets_(std::move(offsets)) { validate(); }
    /** @brief Builds and validates a window. @param offsets Ordered translated-sample offsets. */
    ObservationWindow(std::initializer_list<WindowOffset> offsets) : ObservationWindow(std::vector<WindowOffset>(offsets)) {}

    /** @brief Builds and validates a fixed-size window. @tparam N Number of offsets. @param offsets Ordered translated-sample offsets. */
    template <std::size_t N>
    explicit ObservationWindow(const std::array<WindowOffset, N>& offsets) : ObservationWindow(std::vector<WindowOffset>(offsets.begin(), offsets.end())) {}

    /** @brief Returns the number of offsets. @return Number of window coordinates. */
    [[nodiscard]] std::size_t size() const noexcept { return offsets_.size(); }
    /** @brief Returns one offset without bounds checking. @param index Window coordinate. @return Offset at `index`. */
    [[nodiscard]] const WindowOffset& operator[](std::size_t index) const noexcept { return offsets_[index]; }
    /** @brief Returns all ordered offsets. @return Borrowed offset span. */
    [[nodiscard]] std::span<const WindowOffset> offsets() const noexcept { return offsets_; }
    /** @brief Returns the first offset iterator. @return Beginning iterator. */
    [[nodiscard]] auto begin() const noexcept { return offsets_.begin(); }
    /** @brief Returns the past-the-end offset iterator. @return Ending iterator. */
    [[nodiscard]] auto end() const noexcept { return offsets_.end(); }

  private:
    /** @brief Verifies the finite-window size, uniqueness, and anchor invariants. */
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

    std::vector<WindowOffset> offsets_; ///< Ordered coordinates of the observation window.
};

/** @brief Binary sample-visibility vector encoded in observation-window order. */
class BinaryVisibilityState {
  public:
    /**
     * @brief Builds a validated binary state.
     * @param bits Visible-coordinate bit mask.
     * @param coordinateCount Number of valid low-order coordinates.
     */
    BinaryVisibilityState(std::uint32_t bits, std::size_t coordinateCount) : bits_(bits), coordinateCount_(coordinateCount) {
        if (coordinateCount_ == 0) {
            throw std::invalid_argument("BinaryVisibilityState requires at least one coordinate.");
        }
        if (coordinateCount_ > ObservationWindow::maxNumOffsets) {
            throw std::invalid_argument("BinaryVisibilityState supports at most 32 coordinates.");
        }
        const std::uint32_t validMask = coordinateCount_ == 32 ? std::numeric_limits<std::uint32_t>::max()
                                                               : (std::uint32_t{1} << coordinateCount_) - std::uint32_t{1};
        if ((bits_ & ~validMask) != 0) {
            throw std::invalid_argument("BinaryVisibilityState contains bits outside its coordinate domain.");
        }
    }

    /** @brief Returns the visibility mask. @return Low-order coordinate bits. */
    [[nodiscard]] std::uint32_t bits() const noexcept { return bits_; }
    /** @brief Returns the coordinate domain size. @return Number of valid bits. */
    [[nodiscard]] std::size_t coordinateCount() const noexcept { return coordinateCount_; }

    /** @brief Tests one visibility coordinate. @param coordinate Coordinate index. @return True when the sample is visible. */
    [[nodiscard]] bool isVisible(std::size_t coordinate) const {
        if (coordinate >= coordinateCount_) {
            throw std::out_of_range("BinaryVisibilityState coordinate is outside the observation window.");
        }
        return (bits_ & (std::uint32_t{1} << coordinate)) != 0;
    }

    /**
     * @brief Returns a state with additional coordinates visible.
     * @param enteringOffsetMask Coordinates entering at one anchored entry.
     * @return Updated immutable state value.
     */
    [[nodiscard]] BinaryVisibilityState withEnteringOffsets(std::uint32_t enteringOffsetMask) const {
        return BinaryVisibilityState(bits_ | enteringOffsetMask, coordinateCount_);
    }

    /** @brief Compares visibility masks and coordinate domains. @param lhs Left state. @param rhs Right state. @return True when equal. */
    friend bool operator==(const BinaryVisibilityState& lhs, const BinaryVisibilityState& rhs) = default;

  private:
    /** @brief Tag that restricts unchecked construction to trusted internal access. */
    struct UncheckedConstructionTag {};

    /**
     * @brief Builds a state after its invariants have been established internally.
     * @param bits Visible-coordinate bit mask.
     * @param coordinateCount Number of valid low-order coordinates.
     */
    BinaryVisibilityState(UncheckedConstructionTag, std::uint32_t bits, std::size_t coordinateCount) noexcept
        : bits_(bits), coordinateCount_(coordinateCount) {}

    std::uint32_t bits_ = 0;             ///< Visible-coordinate bit mask.
    std::size_t coordinateCount_ = 0;    ///< Number of meaningful low-order bits.

    friend struct detail::BinaryVisibilityStateAccess;
};

/** @brief One inclusion node and the observation coordinates that enter there. */
struct AnchoredEntryMask {
    NodeId node = InvalidNode;                 ///< Anchored-entry node.
    std::uint32_t enteringOffsetMask = 0;      ///< Coordinates entering at `node`.

    /** @brief Compares node and coordinate mask. @param lhs Left entry. @param rhs Right entry. @return True when equal. */
    friend bool operator==(const AnchoredEntryMask& lhs, const AnchoredEntryMask& rhs) = default;
};

/** @brief Window-coordinate map whose missing entries represent out-of-domain samples. */
class AnchoredEntryMap {
  public:
    /** @brief Builds a coordinate-preserving map. @param anchorPixel Anchor pixel. @param entries Optional anchored entry per coordinate. */
    AnchoredEntryMap(PixelId anchorPixel, std::vector<std::optional<NodeId>> entries) : anchorPixel_(anchorPixel), entries_(std::move(entries)) {}

    /** @brief Returns the anchor pixel. @return Row-major anchor identifier. */
    [[nodiscard]] PixelId anchorPixel() const noexcept { return anchorPixel_; }
    /** @brief Returns the coordinate count. @return Number of mapped coordinates. */
    [[nodiscard]] std::size_t size() const noexcept { return entries_.size(); }
    /** @brief Returns one coordinate entry. @param coordinate Coordinate index. @return Optional anchored-entry node. */
    [[nodiscard]] const std::optional<NodeId>& operator[](std::size_t coordinate) const noexcept { return entries_[coordinate]; }
    /** @brief Returns every coordinate entry. @return Borrowed optional-node span. */
    [[nodiscard]] std::span<const std::optional<NodeId>> entries() const noexcept { return entries_; }

  private:
    PixelId anchorPixel_ = InvalidPixel;                 ///< Row-major anchor identifier.
    std::vector<std::optional<NodeId>> entries_;         ///< Anchored entry for each observation coordinate.
};

using AnchorBranch = std::vector<NodeId>;
using AnchoredEntrySet = std::set<NodeId>;
using OrderedAnchoredEntries = std::vector<AnchoredEntryMask>;

/** @brief One anchor-specific signed change attached to an anchored entry. */
template <class Value> struct EventDelta {
    PixelId anchorPixel = InvalidPixel; ///< Anchor whose local state changed.
    NodeId anchoredEntry = InvalidNode; ///< Node where the new coordinates enter.
    Value value{};                      ///< Signed decision difference.

    /** @brief Compares every event field. @param lhs Left event. @param rhs Right event. @return True when equal. */
    friend bool operator==(const EventDelta& lhs, const EventDelta& rhs) = default;
};

/** @brief Sum of anchor-specific event deltas attached directly to one node. */
template <class Value> struct LocalAttributeIncrement {
    NodeId node = InvalidNode; ///< Dense node slot.
    Value value{};             ///< Non-aggregated signed increment.

    /** @brief Compares node and increment value. @param lhs Left increment. @param rhs Right increment. @return True when equal. */
    friend bool operator==(const LocalAttributeIncrement& lhs, const LocalAttributeIncrement& rhs) = default;
};

/** @brief Final bottom-up aggregated attribute value for one node. */
template <class Value> struct NodeAttribute {
    NodeId node = InvalidNode; ///< Dense node slot.
    Value value{};             ///< Aggregated node attribute.

    /** @brief Compares node and attribute value. @param lhs Left attribute. @param rhs Right attribute. @return True when equal. */
    friend bool operator==(const NodeAttribute& lhs, const NodeAttribute& rhs) = default;
};

/** @brief Pure local decision from a binary visibility state to one value. */
template <class Decision>
concept LocalDecision = requires(const Decision& decision, BinaryVisibilityState state) {
    typename Decision::Value;
    { decision.evaluateLocalDecision(state) } -> std::same_as<typename Decision::Value>;
};

/** @brief Additive event algebra kept separate from the local decision. */
template <class Algebra, class Value>
concept EventAlgebra = requires(const Algebra& algebra, Value& target, const Value& source) {
    { algebra.additiveIdentity() } -> std::same_as<Value>;
    { algebra.addAssign(target, source) } -> std::same_as<void>;
    { algebra.subtractAssign(target, source) } -> std::same_as<void>;
};

} // namespace mmcfilters::local_attributes
