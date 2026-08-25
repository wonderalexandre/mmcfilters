#pragma once

#include "ContourSideAttributeData.hpp"
#include "../../../localAttributes/FiniteWindowLocalAttributeComputer.hpp"

#include <array>
#include <cstdint>
#include <span>
#include <vector>

namespace mmcfilters::attributes::computers::detail {

namespace kernel {

/** @brief Pure local decision for directional contour-side contributions. */
struct ContourSideLocalDecision {
    /** @brief Additive-group value produced by the rule. */
    using Value = ContourSideCounts;

    /** @brief Converts a local visibility state into side counts. @param visibilityState Five-sample visibility state. @return Counts represented by the state.
     */
    [[nodiscard]] Value evaluateLocalDecision(local_attributes::BinaryVisibilityState visibilityState) const {
        Value counts;
        const uint32_t state = visibilityState.bits();
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

};

/** @brief Additive algebra for directional contour-side events. */
struct ContourSideEventAlgebra {
    using Value = ContourSideCounts;

    [[nodiscard]] Value additiveIdentity() const { return {}; }

    /** @brief Adds one group value. @param target Value to update. @param source Value to add. */
    void addAssign(Value& target, const Value& source) const {
        target.contourPixels += source.contourPixels;
        target.exposedSides += source.exposedSides;
        target.north += source.north;
        target.west += source.west;
        target.east += source.east;
        target.south += source.south;
    }

    /** @brief Subtracts one group value. @param target Value to update. @param source Value to subtract. */
    void subtractAssign(Value& target, const Value& source) const {
        target.contourPixels -= source.contourPixels;
        target.exposedSides -= source.exposedSides;
        target.north -= source.north;
        target.west -= source.west;
        target.east -= source.east;
        target.south -= source.south;
    }
};

/** @brief Immutable canonical five-sample contour window. */
inline const local_attributes::ObservationWindow contourObservationWindow{{{0, 0}, {-1, 0}, {0, -1}, {0, 1}, {1, 0}}};

/**
 * @brief Computes per-node contour-side counts through the established event kernel.
 * @param tree Established tree topology and grid domain.
 * @return Per-node directional contour-side counts.
 */
inline std::vector<ContourSideCounts> computeContourSideCounts(const MorphologicalTree& tree) {
    const auto nodeAttributes =
        local_attributes::detail::kernel::computeFiniteWindowLocalAttribute(tree, contourObservationWindow, ContourSideLocalDecision{},
                                                                             ContourSideEventAlgebra{});
    std::vector<ContourSideCounts> counts;
    counts.reserve(nodeAttributes.size());
    for (const auto& nodeAttribute : nodeAttributes) {
        counts.push_back(nodeAttribute.value);
    }
    return counts;
}

} // namespace kernel

/**
 * @brief Side-contour counters computed through local entry events.
 *
 * @details
 * This class is the contour-count primitive for the finite-window layer. It
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
 * do not enter the binary visibility state, so border sides are naturally treated as
 * exposed.
 */
class ContourSideFiniteWindowComputation {
  public:
    /**
     * @brief Reused observation-window offset type.
     */
    using WindowOffset = local_attributes::WindowOffset;

    /** @brief Defines the `ContourSideCounts` alias used by the component. */
    using ContourSideCounts = ::mmcfilters::attributes::computers::detail::ContourSideCounts;

  public:
    /**
     * @brief Computes side-contour counters for every internal node slot.
     *
     * @param tree Tree topology.
     * @return The computed side-contour counters for every internal node slot.
     */
    [[nodiscard]] static std::vector<ContourSideCounts> computeContourSideCounts(const MorphologicalTree& tree) {
        MMCFILTERS_CONTRACT_CHECKED_ONLY(local_attributes::detail::validateFiniteWindowLocalAttributeInput(tree));
        return kernel::computeContourSideCounts(tree);
    }

    /**
     * @brief Projects side counters to scalar contour-pixel counts.
     *
     * @param counts Counters.
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
     * @param counts Counters.
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
