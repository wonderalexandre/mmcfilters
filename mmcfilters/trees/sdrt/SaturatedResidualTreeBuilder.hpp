#pragma once

/**
 * @file SaturatedResidualTreeBuilder.hpp
 * @brief Builder dedicated to the saturated self-dual residual tree.
 */

#include "ResidualTreeBuildStatistics.hpp"
#include "ResidualTreePolicies.hpp"
#include "detail/MinMaxResidualTreeEngine.hpp"

#include <span>
#include <utility>

namespace mmcfilters::sdrt {

/**
 * @brief Builds the saturated residual tree from synchronized component trees.
 *
 * A current regional extremum is eligible only when its complement is
 * connected relative to the configured exterior seed. This API owns every
 * policy that participates in that certification.
 *
 * @tparam T Finite image and residual-node altitude type.
 */
template <AltitudeValue T> class SaturatedResidualTreeBuilder {
  public:
    using altitude_t = T;                                ///< Residual-node altitude type.
    using image_ptr_t = ImagePtr<altitude_t>;             ///< Shared input image type.
    using tree_t = WeightedMorphologicalTree<altitude_t>; ///< Component-tree seed type.

    /**
     * @brief Creates a reusable saturated builder from an extensible option object.
     * @param adjacency Symmetric connected adjacency shared by both seeds.
     * @param infinityPixel Row-major exterior seed.
     * @param options Saturation and ordering policies.
     */
    explicit SaturatedResidualTreeBuilder(RegularGridAdjacency2D adjacency, NodeId infinityPixel,
                                          SaturatedResidualTreeOptions options = {})
        : implementation_(std::move(adjacency), infinityPixel, options) {}

    /**
     * @brief Consumes synchronized min-tree and max-tree seeds.
     * @param image Image reconstructed by the resulting residual hierarchy.
     * @param minTree Min-tree seed consumed by the construction.
     * @param maxTree Max-tree seed consumed by the construction.
     */
    void build(const image_ptr_t& image, tree_t&& minTree, tree_t&& maxTree) {
        implementation_.build(image, std::move(minTree), std::move(maxTree));
    }

    /** @return Number of rows in the last completed build. */
    [[nodiscard]] int getRows() const { return implementation_.getRows(); }
    /** @return Number of columns in the last completed build. */
    [[nodiscard]] int getCols() const { return implementation_.getCols(); }
    /** @return Dense root node identifier. */
    [[nodiscard]] NodeId getRoot() const { return implementation_.getRoot(); }
    /** @return Configured row-major exterior seed. */
    [[nodiscard]] NodeId getInfinityPixel() const noexcept { return implementation_.getInfinityPixel(); }
    /** @return Configured grid adjacency. */
    [[nodiscard]] const RegularGridAdjacency2D& getAdjacency() const noexcept { return implementation_.getAdjacency(); }
    /** @return Configured equal-area tie policy. */
    [[nodiscard]] SdrtTiePolicy getTiePolicy() const noexcept { return implementation_.getTiePolicy(); }
    /** @return Configured dynamic LCA policy. */
    [[nodiscard]] SaturatedMinMaxLcaPolicy getLcaPolicy() const noexcept { return implementation_.getLcaPolicy(); }
    /** @return Configured exact complement fallback. */
    [[nodiscard]] SaturatedMinMaxFallbackPolicy getFallbackPolicy() const noexcept { return implementation_.getFallbackPolicy(); }
    /** @return Parent buffer indexed by residual node id. */
    [[nodiscard]] std::span<const NodeId> getNodeParent() const { return implementation_.getNodeParent(); }
    /** @return Direct residual owner indexed by source pixel id. */
    [[nodiscard]] std::span<const NodeId> getProperPartOwner() const { return implementation_.getProperPartOwner(); }
    /** @return Altitude buffer indexed by residual node id. */
    [[nodiscard]] std::span<const altitude_t> getAltitude() const { return implementation_.getAltitude(); }
    /** @return Correctness-oriented statistics from the last completed build. */
    [[nodiscard]] const ResidualTreeBuildStatistics& getStatistics() const { return implementation_.getStatistics(); }

    /**
     * @brief Transfers the completed result into validated native storage.
     * @param semantics Semantic descriptor assigned to the transferred hierarchy.
     * @return Validated native buffers and their topology proof.
     */
    [[nodiscard]] mmcfilters::detail::ValidatedNativeHierarchy<altitude_t> takeValidatedHierarchy(HierarchySemantics semantics) && {
        return std::move(implementation_).takeValidatedHierarchy(std::move(semantics));
    }

  private:
    detail::MinMaxResidualTreeEngine<altitude_t, true> implementation_; ///< Saturated synchronized min/max engine.
};

} // namespace mmcfilters::sdrt
