#pragma once

/**
 * @file SaturatedResidualTreeBuilder.hpp
 * @brief Builder dedicated to the saturated self-dual residual tree.
 */

#include "ResidualTreeBuildStatistics.hpp"
#include "ResidualTreePolicies.hpp"
#include "detail/SynchronizedResidualTreeEvolution.hpp"

#include <span>
#include <utility>

namespace mmcfilters::sdrt {

/**
 * @brief Builds the saturated residual tree from synchronized component trees.
 *
 * A current regional extremum is eligible only when its complement is
 * connected relative to the configured infinity pixel. This API owns every
 * policy that participates in that certification.
 *
 * @tparam T Finite image and residual-node altitude type.
 */
template <AltitudeValue T> class SaturatedResidualTreeBuilder {
  public:
    using Altitude = T;                                ///< Residual-node altitude type.
    using ImagePointer = ImagePtr<Altitude>;            ///< Shared input image type.
    using Tree = ValuedMorphologicalTree<Altitude>;     ///< Component-tree seed type.

    /**
     * @brief Creates a reusable saturated builder from an extensible option object.
     * @param adjacency Symmetric connected adjacency shared by both seeds.
     * @param infinityPixel Row-major infinity pixel.
     * @param options Saturation and ordering policies.
     */
    explicit SaturatedResidualTreeBuilder(RegularGridAdjacency2D adjacency, PixelId infinityPixel,
                                          SaturatedResidualTreeOptions options = {})
        : implementation_(std::move(adjacency), infinityPixel, options) {}

    /**
     * @brief Consumes synchronized min-tree and max-tree seeds.
     * @param image Image reconstructed by the resulting residual hierarchy.
     * @param minTree Min-tree seed consumed by the construction.
     * @param maxTree Max-tree seed consumed by the construction.
     */
    void build(const ImagePointer& image, Tree&& minTree, Tree&& maxTree) {
        implementation_.build(image, std::move(minTree), std::move(maxTree));
    }

    /** @return Number of rows in the last completed build. */
    [[nodiscard]] int rows() const { return implementation_.rows(); }
    /** @return Number of columns in the last completed build. */
    [[nodiscard]] int columns() const { return implementation_.columns(); }
    /** @return Dense root node identifier. */
    [[nodiscard]] NodeId root() const { return implementation_.root(); }
    /** @return Declared row-major infinity pixel. */
    [[nodiscard]] PixelId infinityPixel() const noexcept { return implementation_.infinityPixel(); }
    /** @return Configured grid adjacency. */
    [[nodiscard]] const RegularGridAdjacency2D& adjacency() const noexcept { return implementation_.adjacency(); }
    /** @return Total order used to define support minima. */
    [[nodiscard]] const SpatialOrder& spatialOrder() const noexcept { return implementation_.spatialOrder(); }
    /** @return Configured dynamic LCA policy. */
    [[nodiscard]] SaturatedMinMaxLcaPolicy lcaPolicy() const noexcept { return implementation_.lcaPolicy(); }
    /** @return Configured exact complement fallback. */
    [[nodiscard]] SaturatedMinMaxFallbackPolicy fallbackPolicy() const noexcept { return implementation_.fallbackPolicy(); }
    /** @return Parent buffer indexed by residual node id. */
    [[nodiscard]] std::span<const NodeId> parents() const { return implementation_.parents(); }
    /** @return Smallest residual node indexed by source pixel id. */
    [[nodiscard]] std::span<const NodeId> smallestNodeMap() const { return implementation_.smallestNodeMap(); }
    /** @return Altitude buffer indexed by residual node id. */
    [[nodiscard]] std::span<const Altitude> nodeAltitudes() const { return implementation_.nodeAltitudes(); }
    /** @return Correctness-oriented statistics from the last completed build. */
    [[nodiscard]] const ResidualTreeBuildStatistics& statistics() const { return implementation_.statistics(); }

    /**
     * @brief Transfers the completed result into validated native storage.
     * @param semantics Semantic descriptor assigned to the transferred hierarchy.
     * @return Validated native buffers and their topology proof.
     */
    [[nodiscard]] mmcfilters::detail::ValidatedNativeHierarchy<Altitude> takeValidatedHierarchy(MorphologicalTreeSemantics semantics) && {
        return std::move(implementation_).takeValidatedHierarchy(std::move(semantics));
    }

  private:
    detail::SynchronizedResidualTreeEvolution<Altitude, true> implementation_; ///< Saturated synchronized evolution.
};

} // namespace mmcfilters::sdrt
