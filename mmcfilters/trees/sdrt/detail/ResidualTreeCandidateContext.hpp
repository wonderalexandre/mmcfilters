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
    /** @brief Stores the boundary pixel marks. */
    GenerationStampSet boundaryPixelMarks;
    /** @brief Stores the boundary owner marks. */
    GenerationStampSet boundaryOwnerMarks;
    /** @brief Stores the support owner marks. */
    GenerationStampSet supportOwnerMarks;
    /** @brief Stores the boundary owners. */
    std::vector<NodeId> boundaryOwners;
    /** @brief Stores the boundary pixels. */
    std::vector<NodeId> boundaryPixels;
    /** @brief Stores the support pixels. */
    std::span<const NodeId> supportPixels;
    /** @brief Stores the support owners. */
    std::vector<NodeId> supportOwners;
    /** @brief Stores the flat-zone representative marks. */
    GenerationStampSet flatZoneRepresentativeMarks;
    /** @brief Stores the flat-zone representatives selected for merging. */
    std::vector<NodeId> flatZoneMergeRepresentatives;
    /** @brief Stores the selected flat-zone representative. */
    NodeId flatZoneRepresentative = InvalidNode;
    /** @brief Stores the dual extremal owner. */
    NodeId dualExtremalOwner = InvalidNode;
    /** @brief Stores the whole support owner. */
    NodeId wholeSupportOwner = InvalidNode;

    /**
     * @brief Constructs a `ResidualTreeCandidateContext` instance.
     *
     * @param numPixels Num pixels used by the operation.
     * @param maxNodeSlots Max node slots used by the operation.
     */
    ResidualTreeCandidateContext(std::size_t numPixels, std::size_t maxNodeSlots)
        : boundaryPixelMarks(numPixels), boundaryOwnerMarks(maxNodeSlots), supportOwnerMarks(maxNodeSlots), flatZoneRepresentativeMarks(numPixels) {
        boundaryOwners.reserve(32);
        boundaryPixels.reserve(64);
        supportOwners.reserve(8);
        flatZoneMergeRepresentatives.reserve(8);
    }
};

} // namespace mmcfilters::sdrt::detail
