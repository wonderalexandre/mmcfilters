#pragma once

#include "ContourSideAttributeData.hpp"
#include "../../../localEvents/EventEngine.hpp"

#include <array>
#include <cstdint>
#include <span>
#include <vector>

namespace mmcfilters::attributes::computers::detail {

namespace kernel {

/** @brief Event-engine policy that accumulates directional contour-side deltas. */
struct ContourSideCountPolicy {
    /** @brief Per-node accumulator used by the event engine. */
    using Bucket = ContourSideCounts;

    /** @brief Converts a local visibility state into side counts. @param state Five-sample local state. @return Counts represented by the state. */
    static Bucket value(uint32_t state) {
        Bucket counts;
        if ((state & uint32_t{1}) == 0) {
            return counts;
        }
        counts.north = (state & (uint32_t{1} << 1)) == 0 ? 1 : 0;
        counts.west = (state & (uint32_t{1} << 2)) == 0 ? 1 : 0;
        counts.east = (state & (uint32_t{1} << 3)) == 0 ? 1 : 0;
        counts.south = (state & (uint32_t{1} << 4)) == 0 ? 1 : 0;
        counts.exposedSides = counts.north + counts.west + counts.east + counts.south;
        counts.contourPixels = counts.exposedSides > 0 ? 1 : 0;
        return counts;
    }

    /** @brief Adds an initial local state. @param bucket Accumulator to update. @param state Initial local state. */
    void applyInitial(Bucket& bucket, uint32_t state) const { addAssign(bucket, value(state)); }
    /**
     * @brief Applies one local-state transition.
     * @param bucket Accumulator to update.
     * @param oldState State before the transition.
     * @param newState State after the transition.
     */
    void applyTransition(Bucket& bucket, uint32_t oldState, uint32_t newState) const { addDelta(bucket, value(newState), value(oldState)); }
    /** @brief Merges a child accumulator into its parent. @param parent Parent accumulator. @param child Child accumulator. */
    void merge(Bucket& parent, const Bucket& child) const { addAssign(parent, child); }

  private:
    /** @brief Adds one side-count bucket to another. @param target Accumulator to update. @param source Counts to add. */
    static void addAssign(Bucket& target, const Bucket& source) {
        target.contourPixels += source.contourPixels;
        target.exposedSides += source.exposedSides;
        target.north += source.north;
        target.west += source.west;
        target.east += source.east;
        target.south += source.south;
    }

    /**
     * @brief Adds the difference between two local side-count states.
     * @param target Accumulator to update.
     * @param current Counts after the transition.
     * @param previous Counts before the transition.
     */
    static void addDelta(Bucket& target, const Bucket& current, const Bucket& previous) {
        target.contourPixels += current.contourPixels - previous.contourPixels;
        target.exposedSides += current.exposedSides - previous.exposedSides;
        target.north += current.north - previous.north;
        target.west += current.west - previous.west;
        target.east += current.east - previous.east;
        target.south += current.south - previous.south;
    }
};

/**
 * @brief Computes per-node contour-side counts through the established event kernel.
 * @param tree Established tree topology and grid domain.
 * @return Per-node directional contour-side counts.
 */
inline std::vector<ContourSideCounts> computeContourSideCounts(const MorphologicalTree& tree) {
    constexpr std::array<local_events::WindowOffset, 5> contourWindow{{{0, 0}, {-1, 0}, {0, -1}, {0, 1}, {1, 0}}};
    return local_events::detail::kernel::computeWithPolicy(tree, std::span<const local_events::WindowOffset>(contourWindow), ContourSideCountPolicy{});
}

} // namespace kernel

/**
 * @brief Side-contour counters computed through local entry events.
 *
 * @details
 * This class is the contour-count primitive for the local event layer. It
 * evaluates a five-sample window:
 *
 * - bit 0: anchor pixel;
 * - bit 1: north neighbour;
 * - bit 2: west neighbour;
 * - bit 3: east neighbour;
 * - bit 4: south neighbour.
 *
 * For each visible anchor, it counts which side neighbours are absent from the
 * support. This yields both the number of contour pixels and the number of
 * exposed sides, including directional side counts. Out-of-domain neighbours
 * do not enter the local event state, so border sides are naturally treated as
 * exposed.
 */
class ContourSideLocalEventComputation {
  public:
    /**
     * @brief Reused event-engine window offset type.
     */
    using WindowOffset = local_events::WindowOffset;

    /** @brief Defines the `ContourSideCounts` alias used by the component. */
    using ContourSideCounts = ::mmcfilters::attributes::computers::detail::ContourSideCounts;

  public:
    /**
     * @brief Computes side-contour counters for every internal node slot.
     *
     * @param tree Tree topology used by the operation.
     * @return The computed side-contour counters for every internal node slot.
     */
    [[nodiscard]] static std::vector<ContourSideCounts> computeContourSideCounts(const MorphologicalTree& tree) {
        constexpr std::array<WindowOffset, 5> contourWindow{{{0, 0}, {-1, 0}, {0, -1}, {0, 1}, {1, 0}}};
        MMCFILTERS_CONTRACT_CHECKED_ONLY(local_events::detail::validateEventEngineInput(tree, contourWindow));
        return kernel::computeContourSideCounts(tree);
    }

    /**
     * @brief Projects side counters to scalar contour-pixel counts.
     *
     * @param counts Counters used by the operation.
     * @return The projected side counters to scalar contour-pixel counts.
     */
    [[nodiscard]] static std::vector<int> projectContourPixels(std::span<const ContourSideCounts> counts) {
        std::vector<int> projected;
        projected.reserve(counts.size());
        for (const ContourSideCounts& count : counts) {
            projected.push_back(count.contourPixels);
        }
        return projected;
    }

    /**
     * @brief Projects side counters to 4-neighbour exposed-side perimeter.
     *
     * @param counts Counters used by the operation.
     * @return The projected side counters to 4-neighbour exposed-side perimeter.
     */
    [[nodiscard]] static std::vector<int> projectExposedSides(std::span<const ContourSideCounts> counts) {
        std::vector<int> projected;
        projected.reserve(counts.size());
        for (const ContourSideCounts& count : counts) {
            projected.push_back(count.exposedSides);
        }
        return projected;
    }
};

} // namespace mmcfilters::attributes::computers::detail
