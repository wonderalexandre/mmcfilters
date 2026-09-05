#pragma once

/**
 * @file SaturatedResidualEligibility.hpp
 * @brief Complete eligibility evaluation for saturated residual-tree candidates.
 */

#include "ResidualTreeCandidateContext.hpp"
#include "SaturatedDynamicLca.hpp"
#include "../ResidualTreeBuildStatistics.hpp"
#include "../ResidualTreePolicies.hpp"
#include "../../ValuedMorphologicalTree.hpp"
#include "../../../utils/GenerationStampSet.hpp"
#include "../../../utils/RegularGridAdjacency2D.hpp"

#include <cstddef>
#include <cstdint>
#include <numeric>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcfilters::sdrt::detail {

/**
 * @brief Evaluates complement connectivity for saturated residual-tree candidates.
 *
 * The evaluator owns both the fast dual-tree LCA certificate and the exact
 * complement traversal used when that certificate is inconclusive.
 */
template <AltitudeValue T> class SaturatedResidualEligibility {
  public:
    using Tree = ValuedMorphologicalTree<T>;                         ///< Mutable component-tree type.
    using DynamicLcaIndex = typename SaturatedLcaTypes<T>::DynamicLcaIndex; ///< Policy-selected LCA index.

    /**
     * @brief Creates and binds all saturated-only certification state.
     * @param numPixels Size of the common pixel domain.
     * @param adjacency Established symmetric pixel adjacency.
     * @param infinityPixel Declared infinity pixel excluded from candidate supports.
     * @param lcaPolicy Dynamic LCA strategy.
     * @param fallbackPolicy Exact complement traversal strategy.
     * @param minTree Current min-tree.
     * @param maxTree Current max-tree.
     */
    SaturatedResidualEligibility(std::size_t numPixels, const RegularGridAdjacency2D& adjacency, PixelId infinityPixel,
                                 SaturatedMinMaxLcaPolicy lcaPolicy, SaturatedMinMaxFallbackPolicy fallbackPolicy,
                                 const Tree& minTree, const Tree& maxTree)
        : supportMarks_(numPixels), visitedMarks_(numPixels), complementLabelByPixel_(numPixels, -1), adjacency_(&adjacency),
          infinityPixel_(infinityPixel), fallbackPolicy_(fallbackPolicy) {
        frontier_.reserve(numPixels);
        minLca_.bind(minTree, lcaPolicy);
        maxLca_.bind(maxTree, lcaPolicy);
    }

    /**
     * @brief Tests whether one prepared candidate is eligible for saturated contraction.
     * @param primal Component tree that owns the candidate.
     * @param dual Opposite-polarity component tree.
     * @param candidateNode Current primal-tree leaf.
     * @param primalIsMaxTree Whether `primal` is the max-tree.
     * @param containsInfinityPixel Whether the candidate support contains the infinity pixel.
     * @param context Prepared common support and boundary data.
     * @param statistics Build diagnostics updated by exact fallback traversals.
     * @return `true` when the candidate support has connected complement.
     */
    [[nodiscard]] bool isEligible(const Tree& primal, const Tree& dual, NodeId candidateNode, bool primalIsMaxTree, bool containsInfinityPixel,
                                  const ResidualTreeCandidateContext& context, ResidualTreeBuildStatistics& statistics) {
        if (containsInfinityPixel) {
            return false;
        }
        if (context.boundarySmallestNodes.empty()) {
            throw std::logic_error("Saturated residual eligibility requires prepared boundary smallest nodes.");
        }

        NodeId common = context.boundarySmallestNodes.front();
        DynamicLcaIndex& dualLca = primalIsMaxTree ? minLca_ : maxLca_;
        for (std::size_t index = 1; index < context.boundarySmallestNodes.size(); ++index) {
            common = lowestCommonAncestor(dual, dualLca, common, context.boundarySmallestNodes[index]);
        }

        const T sourceLevel = primal.nodeAltitude(candidateNode);
        const T connectingLevel = dual.nodeAltitude(common);
        const bool certifiedByDualTree = primalIsMaxTree ? connectingLevel < sourceLevel : connectingLevel > sourceLevel;
        if (certifiedByDualTree) {
            return true;
        }

        const MorphologicalTree& dualTopology = dual.topology();
        bool commonContainsSupport = false;
        for (NodeId smallestNodeId : context.supportSmallestNodes) {
            NodeId cursor = smallestNodeId;
            while (true) {
                if (cursor == common) {
                    commonContainsSupport = true;
                    break;
                }
                const NodeId parent = dualTopology.parent(cursor);
                if (parent == cursor) {
                    break;
                }
                cursor = parent;
            }
            if (commonContainsSupport) {
                break;
            }
        }
        if (!commonContainsSupport) {
            return true;
        }

        ++statistics.complementTraversalCertificates;
        supportMarks_.resetAll();
        for (PixelId pixel : context.supportPixels) {
            supportMarks_.mark(static_cast<std::size_t>(pixel));
        }
        return exactComplementTraversal(context);
    }

    /**
     * @brief Updates the opposite-tree LCA index after committing one candidate.
     * @param primalIsMaxTree Whether the committed candidate came from the max-tree.
     * @param changedRoots Roots whose opposite-tree topology changed.
     */
    void noteDualMutation(bool primalIsMaxTree, std::span<const NodeId> changedRoots) {
        DynamicLcaIndex& dualLca = primalIsMaxTree ? minLca_ : maxLca_;
        dualLca.noteMutation(changedRoots);
    }

  private:
    DynamicLcaIndex maxLca_;                         ///< Dynamic index for the mutable max-tree.
    DynamicLcaIndex minLca_;                         ///< Dynamic index for the mutable min-tree.
    GenerationStampSet supportMarks_;                ///< Candidate-support pixel marks.
    GenerationStampSet visitedMarks_;                ///< Complement traversal marks.
    std::vector<PixelId> frontier_;                  ///< DFS stack or multi-source traversal queue.
    std::vector<int> complementLabelByPixel_;        ///< Multi-source complement label by pixel.
    std::vector<int> complementLabelParent_;         ///< Union-find parent by complement label.
    std::vector<std::uint8_t> complementLabelRank_;  ///< Union-find rank by complement label.
    const RegularGridAdjacency2D* adjacency_;         ///< Non-owning adjacency used by exact traversal.
    PixelId infinityPixel_;                         ///< Declared infinity pixel.
    SaturatedMinMaxFallbackPolicy fallbackPolicy_;   ///< Selected exact traversal strategy.

    /**
     * @brief Computes an LCA by following parent links in the current tree.
     *
     * Mutation-specific exception to MorphologicalTree's static LCA cache:
     * strict min/max altitude order selects which parent link to follow, permitting
     * an O(height) query with no global
     * rebuild after a candidate edit. It also resolves dirty snapshot queries.
     * Keep this fallback confined to saturated certification.
     * @param tree Mutable component tree containing both nodes.
     * @param first First node identifier.
     * @param second Second node identifier.
     * @return Lowest common ancestor of `first` and `second` in the current tree.
     */
    [[nodiscard]] static NodeId parentClimbLowestCommonAncestor(const Tree& tree, NodeId first, NodeId second) {
        const MorphologicalTree& topology = tree.topology();
        const bool minTree = topology.kind() == MorphologicalTreeKind::MinTree;
        NodeId lhs = first;
        NodeId rhs = second;
        while (lhs != rhs) {
            const T lhsAltitude = tree.nodeAltitude(lhs);
            const T rhsAltitude = tree.nodeAltitude(rhs);
            bool climbLhs = false;
            bool climbRhs = false;
            if (lhsAltitude == rhsAltitude) {
                climbLhs = true;
                climbRhs = true;
            } else if ((minTree && lhsAltitude < rhsAltitude) || (!minTree && lhsAltitude > rhsAltitude)) {
                climbLhs = true;
            } else {
                climbRhs = true;
            }
            if (climbLhs) {
                const NodeId parent = topology.parent(lhs);
                if (parent == lhs) {
                    throw std::runtime_error("Min/max residual altitude-aligned LCA reached one root before convergence.");
                }
                lhs = parent;
            }
            if (climbRhs) {
                const NodeId parent = topology.parent(rhs);
                if (parent == rhs) {
                    throw std::runtime_error("Min/max residual altitude-aligned LCA reached one root before convergence.");
                }
                rhs = parent;
            }
            if (!topology.isAlive(lhs) || !topology.isAlive(rhs)) {
                throw std::runtime_error("Min/max residual altitude-aligned LCA reached a dead node.");
            }
        }
        return lhs;
    }

    /**
     * @brief Resolves an LCA through the configured index or by following parent links.
     * @param tree Mutable component tree containing both nodes.
     * @param index Policy-selected dynamic LCA index.
     * @param first First node identifier.
     * @param second Second node identifier.
     * @return Lowest common ancestor of `first` and `second`.
     */
    [[nodiscard]] static NodeId lowestCommonAncestor(const Tree& tree, DynamicLcaIndex& index, NodeId first, NodeId second) {
        const auto indexed = index.findLowestCommonAncestor(first, second);
        return indexed.has_value() ? *indexed : parentClimbLowestCommonAncestor(tree, first, second);
    }

    /**
     * @brief Checks complement connectivity from one boundary seed.
     * @param context Prepared support and deduplicated boundary pixels.
     * @return `true` when one traversal reaches every boundary pixel outside the support.
     */
    [[nodiscard]] bool exactSingleSourceComplementTraversal(const ResidualTreeCandidateContext& context) {
        visitedMarks_.resetAll();
        frontier_.clear();
        if (supportMarks_.isMarked(static_cast<std::size_t>(infinityPixel_))) {
            return false;
        }
        if (context.boundaryPixels.empty()) {
            throw std::runtime_error("Min/max residual complement traversal requires a non-empty external boundary.");
        }
        const PixelId seed = context.boundaryPixels.front();
        frontier_.push_back(seed);
        visitedMarks_.mark(static_cast<std::size_t>(seed));
        std::size_t reachedBoundaryPixels = 1;
        while (!frontier_.empty() && reachedBoundaryPixels < context.boundaryPixels.size()) {
            const PixelId pixel = frontier_.back();
            frontier_.pop_back();
            for (PixelId neighbor : adjacency_->getNeighborIndices(pixel)) {
                const std::size_t index = static_cast<std::size_t>(neighbor);
                if (supportMarks_.isMarked(index) || visitedMarks_.isMarked(index)) {
                    continue;
                }
                visitedMarks_.mark(index);
                frontier_.push_back(neighbor);
                if (context.boundaryPixelMarks.isMarked(index)) {
                    ++reachedBoundaryPixels;
                }
            }
        }
        return reachedBoundaryPixels == context.boundaryPixels.size();
    }

    /**
     * @brief Checks complement connectivity by merging simultaneous boundary traversals.
     * @param context Prepared support and deduplicated boundary pixels.
     * @return `true` when all boundary traversal labels merge outside the support.
     */
    [[nodiscard]] bool exactMultiSourceComplementTraversal(const ResidualTreeCandidateContext& context) {
        visitedMarks_.resetAll();
        frontier_.clear();
        if (supportMarks_.isMarked(static_cast<std::size_t>(infinityPixel_))) {
            return false;
        }
        if (context.boundaryPixels.empty()) {
            throw std::runtime_error("Min/max residual complement traversal requires a non-empty external boundary.");
        }

        const std::size_t numLabels = context.boundaryPixels.size();
        complementLabelParent_.resize(numLabels);
        std::iota(complementLabelParent_.begin(), complementLabelParent_.end(), 0);
        complementLabelRank_.assign(numLabels, std::uint8_t{0});
        for (std::size_t label = 0; label < numLabels; ++label) {
            const PixelId pixel = context.boundaryPixels[label];
            const std::size_t pixelIndex = static_cast<std::size_t>(pixel);
            visitedMarks_.mark(pixelIndex);
            complementLabelByPixel_[pixelIndex] = static_cast<int>(label);
            frontier_.push_back(pixel);
        }

        std::size_t activeLabels = numLabels;
        const auto findLabel = [this](int label) {
            int root = label;
            while (complementLabelParent_[static_cast<std::size_t>(root)] != root) {
                root = complementLabelParent_[static_cast<std::size_t>(root)];
            }
            while (label != root) {
                const int parent = complementLabelParent_[static_cast<std::size_t>(label)];
                complementLabelParent_[static_cast<std::size_t>(label)] = root;
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
            if (complementLabelRank_[static_cast<std::size_t>(lhs)] < complementLabelRank_[static_cast<std::size_t>(rhs)]) {
                std::swap(lhs, rhs);
            }
            complementLabelParent_[static_cast<std::size_t>(rhs)] = lhs;
            if (complementLabelRank_[static_cast<std::size_t>(lhs)] == complementLabelRank_[static_cast<std::size_t>(rhs)]) {
                ++complementLabelRank_[static_cast<std::size_t>(lhs)];
            }
            --activeLabels;
        };

        std::size_t head = 0;
        while (head < frontier_.size() && activeLabels > 1) {
            const PixelId pixel = frontier_[head++];
            const int pixelLabel = complementLabelByPixel_[static_cast<std::size_t>(pixel)];
            for (PixelId neighbor : adjacency_->getNeighborIndices(pixel)) {
                const std::size_t neighborIndex = static_cast<std::size_t>(neighbor);
                if (supportMarks_.isMarked(neighborIndex)) {
                    continue;
                }
                if (!visitedMarks_.isMarked(neighborIndex)) {
                    visitedMarks_.mark(neighborIndex);
                    complementLabelByPixel_[neighborIndex] = pixelLabel;
                    frontier_.push_back(neighbor);
                    continue;
                }
                uniteLabels(pixelLabel, complementLabelByPixel_[neighborIndex]);
                if (activeLabels == 1) {
                    break;
                }
            }
        }
        return activeLabels == 1;
    }

    /**
     * @brief Dispatches the configured exact complement traversal.
     * @param context Prepared support and boundary data.
     * @return `true` when the candidate complement is connected.
     */
    [[nodiscard]] bool exactComplementTraversal(const ResidualTreeCandidateContext& context) {
        if (fallbackPolicy_ == SaturatedMinMaxFallbackPolicy::SingleSourceDepthFirst) {
            return exactSingleSourceComplementTraversal(context);
        }
        return exactMultiSourceComplementTraversal(context);
    }
};

} // namespace mmcfilters::sdrt::detail
