#pragma once

/**
 * @file ResidualTreeCandidateContext.hpp
 * @brief Reusable candidate data shared by residual-tree construction modes.
 */

#include "../../../utils/Common.hpp"
#include "../../../utils/GenerationStampSet.hpp"

#include <cstddef>
#include <span>
#include <vector>

namespace mmcfilters::sdrt::detail {

/** @brief Owns reusable candidate preparation data shared by both residual-tree modes. */
struct ResidualTreeCandidateContext {
    /** @brief Boundary pixel marks. */
    GenerationStampSet boundaryPixelMarks;
    /** @brief Boundary smallest node marks. */
    GenerationStampSet boundarySmallestNodeMarks;
    /** @brief Support smallest node marks. */
    GenerationStampSet supportSmallestNodeMarks;
    /** @brief Dense node identifier of the boundary smallest nodes. */
    std::vector<NodeId> boundarySmallestNodes;
    /** @brief Pixel identifier of the boundary pixels. */
    std::vector<PixelId> boundaryPixels;
    /** @brief Pixel identifier of the support pixels. */
    std::span<const PixelId> supportPixels;
    /** @brief Dense node identifier of the support smallest nodes. */
    std::vector<NodeId> supportSmallestNodes;
    /** @brief Flat-zone representative marks. */
    GenerationStampSet flatZoneRepresentativeMarks;
    /** @brief Pixel identifier of the flat-zone representatives selected for merging. */
    std::vector<PixelId> flatZoneMergeRepresentatives;
    /** @brief Pixel identifier of the selected flat-zone representative. */
    PixelId flatZoneRepresentative = InvalidPixel;
    /** @brief Dense identifier of the dual extremal smallest node. */
    NodeId dualExtremalSmallestNode = InvalidNode;
    /** @brief Dense identifier of the whole support smallest node. */
    NodeId wholeSupportSmallestNode = InvalidNode;

    /**
     * @brief Constructs a `ResidualTreeCandidateContext` instance.
     *
     * @param numPixels Num pixels.
     * @param maxNodeSlots Max node slots.
     */
    ResidualTreeCandidateContext(std::size_t numPixels, std::size_t maxNodeSlots)
        : boundaryPixelMarks(numPixels), boundarySmallestNodeMarks(maxNodeSlots), supportSmallestNodeMarks(maxNodeSlots), flatZoneRepresentativeMarks(numPixels) {
        boundarySmallestNodes.reserve(32);
        boundaryPixels.reserve(64);
        supportSmallestNodes.reserve(8);
        flatZoneMergeRepresentatives.reserve(8);
    }
};

} // namespace mmcfilters::sdrt::detail
