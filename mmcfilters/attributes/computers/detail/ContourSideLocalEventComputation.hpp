#pragma once

#include "ContourSideAttributeData.hpp"
#include "../../../localEvents/EventEngine.hpp"

#include <cstdint>
#include <span>
#include <vector>

namespace mmcfilters::attributes::computers::detail {

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

  private:
    /**
     * @brief Policy that accumulates contour-side counters directly.
     */
    struct ContourSideCountPolicy {
        /** @brief Defines the `Bucket` alias used by the component. */
        using Bucket = ContourSideCounts;

        /**
         * @brief Returns the contour-side contribution of one five-bit state.
         *
         * @param state State read or updated by the operation.
         * @return The contour-side contribution of one five-bit state.
         */
        static Bucket value(uint32_t state) {
            Bucket counts;
            const bool anchorVisible = (state & uint32_t{1}) != 0;
            if (!anchorVisible) {
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

        /**
         * @brief Adds the initial local-state contribution for one anchor.
         *
         * @param bucket Bucket read or updated by the operation.
         * @param state State read or updated by the operation.
         */
        void applyInitial(Bucket& bucket, uint32_t state) const { addAssign(bucket, value(state)); }

        /**
         * @brief Applies the difference between consecutive local states.
         *
         * @param bucket Bucket read or updated by the operation.
         * @param oldState Local state before the transition.
         * @param newState Local state after the transition.
         */
        void applyTransition(Bucket& bucket, uint32_t oldState, uint32_t newState) const { addDelta(bucket, value(newState), value(oldState)); }

        /**
         * @brief Accumulates child side counters into the parent support.
         *
         * @param parent Parent node used by the operation.
         * @param child Child node used by the operation.
         */
        void merge(Bucket& parent, const Bucket& child) const { addAssign(parent, child); }

      private:
        /**
         * @brief Adds a source contour-side accumulator into a destination accumulator.
         *
         * @param target Destination value or object.
         * @param source Source value or object.
         */
        static void addAssign(Bucket& target, const Bucket& source) {
            target.contourPixels += source.contourPixels;
            target.exposedSides += source.exposedSides;
            target.north += source.north;
            target.west += source.west;
            target.east += source.east;
            target.south += source.south;
        }

        /**
         * @brief Adds delta.
         *
         * @param target Destination value or object.
         * @param current Current item in the traversal or local event.
         * @param previous Previous item in the local event sequence.
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

  public:
    /**
     * @brief Computes side-contour counters for every internal node slot.
     *
     * @param tree Tree topology used by the operation.
     * @return The computed side-contour counters for every internal node slot.
     */
    [[nodiscard]] static std::vector<ContourSideCounts> computeContourSideCounts(const MorphologicalTree& tree) {
        const std::vector<WindowOffset> contourWindow = {
            {0, 0}, {-1, 0}, {0, -1}, {0, 1}, {1, 0},
        };

        return local_events::EventEngine::computeWithPolicy(tree, contourWindow, ContourSideCountPolicy{});
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
