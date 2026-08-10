#pragma once

/**
 * @file SaturatedResidualEligibility.hpp
 * @brief State and exact complement traversal exclusive to saturated residual trees.
 */

#include "ResidualTreeCandidateContext.hpp"
#include "SaturatedDynamicLca.hpp"
#include "../ResidualTreePolicies.hpp"
#include "../../../utils/GenerationStampSet.hpp"
#include "../../../utils/RegularGridAdjacency2D.hpp"

#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcfilters::sdrt::detail {

/** @brief Scratch state and exact complement-connectivity tests used only by saturated construction. */
template <AltitudeValue T> class SaturatedResidualEligibilityState {
  public:
    /** @brief Dynamic LCA index used by saturated fast certificates. */
    using DynamicLcaIndex = typename SaturatedLcaTypes<T>::DynamicLcaIndex;
    /** @brief Candidate support and boundary scratch shared with the engine. */
    using CandidateContextScratch = ResidualTreeCandidateContext;

    DynamicLcaIndex maxLca;                         ///< Dynamic index for the mutable max-tree.
    DynamicLcaIndex minLca;                         ///< Dynamic index for the mutable min-tree.
    GenerationStampSet supportMarks;                ///< Candidate-support pixel marks.
    GenerationStampSet visitedMarks;                ///< Complement traversal marks.
    std::vector<NodeId> frontier;                   ///< DFS stack or multi-source traversal queue.
    std::vector<int> complementLabelByPixel;        ///< Multi-source complement label by pixel.
    std::vector<int> complementLabelParent;         ///< Union-find parent by complement label.
    std::vector<std::uint8_t> complementLabelRank;  ///< Union-find rank by complement label.

    /**
     * @brief Allocates saturated-only certification state.
     * @param numPixels Size of the common pixel domain.
     * @param adjacency Established symmetric pixel adjacency.
     * @param infinityPixel Exterior seed excluded from candidate supports.
     * @param fallbackPolicy Exact complement traversal policy.
     */
    SaturatedResidualEligibilityState(std::size_t numPixels, const RegularGridAdjacency2D& adjacency,
                                      NodeId infinityPixel, SaturatedMinMaxFallbackPolicy fallbackPolicy)
        : supportMarks(numPixels), visitedMarks(numPixels), complementLabelByPixel(numPixels, -1),
          adjacency_(&adjacency), infinityPixel_(infinityPixel), fallbackPolicy_(fallbackPolicy) {
        frontier.reserve(numPixels);
    }

    /**
     * @brief Checks complement connectivity from one reference seed.
     *
     * @param context Prepared candidate context.
     * @return `true` when every complement boundary pixel is reachable from the exterior seed.
     */
    [[nodiscard]] bool exactSingleSourceComplementTraversal(const CandidateContextScratch& context) {
        visitedMarks.resetAll();
        frontier.clear();
        if (supportMarks.isMarked(static_cast<std::size_t>(infinityPixel_))) {
            return false;
        }
        if (context.boundaryPixels.empty()) {
            throw std::runtime_error("Min/max residual complement traversal requires a non-empty external boundary.");
        }
        const NodeId seed = context.boundaryPixels.front();
        frontier.push_back(seed);
        visitedMarks.mark(static_cast<std::size_t>(seed));
        std::size_t reachedBoundaryPixels = 1;
        while (!frontier.empty() && reachedBoundaryPixels < context.boundaryPixels.size()) {
            const NodeId pixel = frontier.back();
            frontier.pop_back();
            for (NodeId neighbor : adjacency_->getNeighborIndices(pixel)) {
                const std::size_t index = static_cast<std::size_t>(neighbor);
                if (supportMarks.isMarked(index) || visitedMarks.isMarked(index)) {
                    continue;
                }
                visitedMarks.mark(index);
                frontier.push_back(neighbor);
                if (context.boundaryPixelMarks.isMarked(index)) {
                    ++reachedBoundaryPixels;
                }
            }
        }
        return reachedBoundaryPixels == context.boundaryPixels.size();
    }

    /**
     * @brief Checks complement connectivity from all boundary seeds.
     *
     * @param context Prepared candidate context.
     * @return `true` when all complement boundary seeds belong to one connected component.
     */
    [[nodiscard]] bool exactMultiSourceComplementTraversal(const CandidateContextScratch& context) {
        visitedMarks.resetAll();
        frontier.clear();
        if (supportMarks.isMarked(static_cast<std::size_t>(infinityPixel_))) {
            return false;
        }
        if (context.boundaryPixels.empty()) {
            throw std::runtime_error("Min/max residual complement traversal requires a non-empty external boundary.");
        }

        const std::size_t numLabels = context.boundaryPixels.size();
        complementLabelParent.resize(numLabels);
        std::iota(complementLabelParent.begin(), complementLabelParent.end(), 0);
        complementLabelRank.assign(numLabels, std::uint8_t{0});

        for (std::size_t label = 0; label < numLabels; ++label) {
            const NodeId pixel = context.boundaryPixels[label];
            const std::size_t pixelIndex = static_cast<std::size_t>(pixel);
            visitedMarks.mark(pixelIndex);
            complementLabelByPixel[pixelIndex] = static_cast<int>(label);
            frontier.push_back(pixel);
        }

        std::size_t activeLabels = numLabels;
        const auto findLabel = [this](int label) {
            int root = label;
            while (complementLabelParent[static_cast<std::size_t>(root)] != root) {
                root = complementLabelParent[static_cast<std::size_t>(root)];
            }
            while (label != root) {
                const int parent = complementLabelParent[static_cast<std::size_t>(label)];
                complementLabelParent[static_cast<std::size_t>(label)] = root;
                label = parent;
            }
            return root;
        };
        const auto uniteLabels = [this, &findLabel, &activeLabels](int first, int second) {
            int lhs = findLabel(first);
            int rhs = findLabel(second);
            if (lhs == rhs) {
                return;
            }
            if (complementLabelRank[static_cast<std::size_t>(lhs)] < complementLabelRank[static_cast<std::size_t>(rhs)]) {
                std::swap(lhs, rhs);
            }
            complementLabelParent[static_cast<std::size_t>(rhs)] = lhs;
            if (complementLabelRank[static_cast<std::size_t>(lhs)] == complementLabelRank[static_cast<std::size_t>(rhs)]) {
                ++complementLabelRank[static_cast<std::size_t>(lhs)];
            }
            --activeLabels;
        };

        std::size_t head = 0;
        while (head < frontier.size() && activeLabels > 1) {
            const NodeId pixel = frontier[head++];
            const int pixelLabel = complementLabelByPixel[static_cast<std::size_t>(pixel)];
            for (NodeId neighbor : adjacency_->getNeighborIndices(pixel)) {
                const std::size_t neighborIndex = static_cast<std::size_t>(neighbor);
                if (supportMarks.isMarked(neighborIndex)) {
                    continue;
                }
                if (!visitedMarks.isMarked(neighborIndex)) {
                    visitedMarks.mark(neighborIndex);
                    complementLabelByPixel[neighborIndex] = pixelLabel;
                    frontier.push_back(neighbor);
                    continue;
                }
                uniteLabels(pixelLabel, complementLabelByPixel[neighborIndex]);
                if (activeLabels == 1) {
                    break;
                }
            }
        }
        return activeLabels == 1;
    }

    /**
     * @brief Runs the configured exact complement-connectivity traversal.
     *
     * @param context Prepared candidate context.
     * @return `true` when the candidate support has connected complement.
     */
    [[nodiscard]] bool exactComplementTraversal(const CandidateContextScratch& context) {
        if (fallbackPolicy_ == SaturatedMinMaxFallbackPolicy::SingleSourceDepthFirst) {
            return exactSingleSourceComplementTraversal(context);
        }
        return exactMultiSourceComplementTraversal(context);
    }

  private:
    const RegularGridAdjacency2D* adjacency_;       ///< Non-owning adjacency used by exact traversal.
    NodeId infinityPixel_;                          ///< Exterior seed pixel.
    SaturatedMinMaxFallbackPolicy fallbackPolicy_; ///< Selected exact traversal strategy.
};

/** @brief Empty eligibility state compiled into unrestricted residual construction. */
struct UnrestrictedResidualEligibilityState {
    /**
     * @brief Accepts the common construction signature without retaining saturated-only state.
     * @param numPixels Unused common pixel-domain size.
     * @param adjacency Unused common adjacency.
     * @param infinityPixel Unused exterior seed.
     * @param fallbackPolicy Unused exact traversal policy.
     */
    UnrestrictedResidualEligibilityState([[maybe_unused]] std::size_t numPixels,
                                         [[maybe_unused]] const RegularGridAdjacency2D& adjacency,
                                         [[maybe_unused]] NodeId infinityPixel,
                                         [[maybe_unused]] SaturatedMinMaxFallbackPolicy fallbackPolicy) noexcept {}
};

} // namespace mmcfilters::sdrt::detail
