#pragma once

#include "MorphologicalTree.hpp"
#include "TreeAltitudeAlgorithms.hpp"
#include "WeightedMorphologicalTree.hpp"
#include "detail/MorphologicalTreeConstructionTag.hpp"
#include "sdrt/SelfDualResidualTreeBuilder.hpp"

#include <optional>
#include <span>
#include <utility>

namespace mmcfilters {

/**
 * @brief Public construction facade for all high-level morphological trees.
 *
 * `MorphologicalTreeFactory` is the single public entry point for creating
 * tree instances from images, static Higra hierarchies, and self-dual residual
 * construction. The factory owns only construction orchestration; it does not
 * store mutable build state and every construction method returns a fully
 * initialized weighted owner.
 *
 * Design contract:
 *
 * - image-driven builders create the topology first and then infer one altitude
 *   value per internal node from the source image;
 * - Higra imports preserve the external Higra node-id domain for altitude
 *   mapping, then expose the result through the internal dense `NodeId` domain;
 * - SDRT construction delegates the iterative residual build to
     *   `sdrt::SelfDualResidualTreeBuilder<T>` and only wraps the materialized buffers
 *   into the standard weighted tree representation;
 * - callers that only need topology should use `WeightedMorphologicalTree<std::uint8_t>::topology()`
 *   instead of constructing a bare `MorphologicalTree` directly;
 * - max/min image factories derive altitude type from the image type. `uint8_t`
 *   uses counting sort; other supported altitude values use comparison sort;
 * - static Higra imports can preserve a non-canonical altitude type because
 *   the altitude array itself supplies that type.
 *
 * Keeping construction in this factory prevents `MorphologicalTree` and
 * `WeightedMorphologicalTree<std::uint8_t>` from depending on public factory logic while still
 * allowing their tag-based constructors to remain available for header-only
 * composition inside this library.
 */
class MorphologicalTreeFactory {
private:
    /**
     * @brief Creates the internal construction capability token.
     *
     * The token restricts direct calls to low-level constructors. Public users
     * create trees through this factory, which is the only class allowed to
     * instantiate the tag.
     */
    static constexpr detail::MorphologicalTreeConstructionTag tag() noexcept {
        return detail::MorphologicalTreeConstructionTag{};
    }

public:
    /**
     * @brief Builds a typed weighted max-tree from an image.
     *
     * @param img Non-null, non-empty image. Its row/column domain becomes the
     * proper-part domain of the resulting tree, and its pixel type becomes the
     * altitude type of the returned owner.
     * @param radius Radius used by the component-tree adjacency relation.
     *
     * @return A `WeightedMorphologicalTree<T>` whose topology is a max-tree
     * and whose altitude buffer is indexed by internal `NodeId`.
     *
     * @throws std::invalid_argument if the image domain is invalid.
     */
    template<AltitudeValue T>
    [[nodiscard]] static WeightedMorphologicalTree<T> createMaxTree(ImagePtr<T> img, double radius = 1.5) {
        TreeAltitudeAlgorithms::validateFiniteImageAltitudes(img, "MorphologicalTreeFactory::createMaxTree image");
        auto topology = MorphologicalTree(tag(), img, ComponentTreeKind::MAX_TREE, radius);
        return WeightedMorphologicalTree<T>(tag(), std::move(topology), std::move(img));
    }

    /**
     * @brief Builds a typed weighted min-tree from an image.
     *
     * @param img Non-null, non-empty image. Its row/column domain becomes the
     * proper-part domain of the resulting tree, and its pixel type becomes the
     * altitude type of the returned owner.
     * @param radius Radius used by the component-tree adjacency relation.
     *
     * @return A `WeightedMorphologicalTree<T>` whose topology is a min-tree
     * and whose altitude buffer is indexed by internal `NodeId`.
     *
     * @throws std::invalid_argument if the image domain is invalid.
     */
    template<AltitudeValue T>
    [[nodiscard]] static WeightedMorphologicalTree<T> createMinTree(ImagePtr<T> img, double radius = 1.5) {
        TreeAltitudeAlgorithms::validateFiniteImageAltitudes(img, "MorphologicalTreeFactory::createMinTree image");
        auto topology = MorphologicalTree(tag(), img, ComponentTreeKind::MIN_TREE, radius);
        return WeightedMorphologicalTree<T>(tag(), std::move(topology), std::move(img));
    }

    /**
     * @brief Builds a weighted tree of shapes from an 8-bit image.
     *
     * @param img Non-null, non-empty image.
     * @param interpolation Interpolation policy used by the ToS builder.
     * @param infinitySeedRow Row coordinate of the infinity seed expected by
     * the underlying ToS construction.
     * @param infinitySeedCol Column coordinate of the infinity seed expected by
     * the underlying ToS construction.
     *
     * @return A `WeightedMorphologicalTree<std::uint8_t>` whose topology is a tree of shapes
     * and whose altitude buffer is indexed by internal `NodeId`.
     *
     * @throws std::invalid_argument if the image domain or ToS parameters are
     * rejected by the underlying builder.
     */
    [[nodiscard]] static WeightedMorphologicalTree<std::uint8_t> createTreeOfShapes(
        ImageUInt8Ptr img,
        ToSInterpolation interpolation = ToSInterpolation::SelfDual,
        int infinitySeedRow = ToSDefaultInfinityRow,
        int infinitySeedCol = ToSDefaultInfinityCol) {
        auto topology = MorphologicalTree(tag(), img, interpolation, infinitySeedRow, infinitySeedCol);
        return WeightedMorphologicalTree<std::uint8_t>(tag(), std::move(topology), std::move(img));
    }

    /**
     * @brief Imports a static Higra parent/altitude hierarchy.
     *
     * Higra arrays use one compact id domain containing leaves followed by
     * internal nodes. For an image domain with `rows * cols` leaves, leaf ids
     * are `[0, rows * cols)` and internal ids are `[rows * cols, parent.size())`.
     * Every leaf must point to an internal node, every internal node must point
     * to another internal node or to itself, and exactly one internal node must
     * be self-parented as the root.
     *
     * @param higraParent Parent array in the compact Higra id domain.
     * @param higraAltitude Altitudes in the same compact Higra id domain.
     * @param rows Number of rows in the image/proper-part domain.
     * @param cols Number of columns in the image/proper-part domain.
     * @param kind Semantic tree kind represented by the imported hierarchy.
     * @param adjacency Required for max-tree and min-tree imports. Optional for
     * tree-of-shapes imports, whose topology can be interpreted without a
     * component-tree adjacency relation.
     *
     * @return A typed weighted tree whose internal altitude buffer has already
     * been remapped from Higra ids to internal `NodeId` slots.
     *
     * @throws std::invalid_argument if the hierarchy is malformed, if the
     * altitude domain does not match the parent domain, or if a max/min import
     * does not provide adjacency.
     */
    template<AltitudeValue T>
    [[nodiscard]] static WeightedMorphologicalTree<T> createFromHigraParent(std::span<const NodeId> higraParent, std::span<const T> higraAltitude, int rows, int cols, MorphologicalTreeKind kind, std::optional<AdjacencyRelation> adjacency = std::nullopt) {
        auto topology = MorphologicalTree(tag(), higraParent, rows, cols, kind, std::move(adjacency));
        return WeightedMorphologicalTree<T>(tag(), std::move(topology), detail::AltitudeInput<T>{higraAltitude, detail::AltitudeDomain::HigraNodeIds});
    }

    /**
     * @brief Builds a typed weighted self-dual residual tree from an image.
     *
     * The factory first creates the initial min-tree and max-tree using the same
     * image and radius, then passes both weighted trees to
     * `sdrt::SelfDualResidualTreeBuilder<T>`. The resulting native parent,
     * proper-part owner, root, and altitude buffers are wrapped into a
     * `WeightedMorphologicalTree<T>`. The image pixel type defines the SDRT
     * altitude type; no arbitrary altitude conversion is performed.
     *
     * @param img Non-null, non-empty image.
     * @param radius Radius used both for the initial component trees and for the
     * dual-tree adjustment adjacency.
     *
     * @return A weighted tree whose topology has type
     * `MorphologicalTreeKind::SELF_DUAL_RESIDUAL_TREE`.
     *
     * @throws std::invalid_argument if the image domain is invalid.
     * @throws std::runtime_error if the SDRT build detects an inconsistent
     * residual state.
     */
    template<AltitudeValue T>
    [[nodiscard]] static WeightedMorphologicalTree<T> createSelfDualResidualTree(ImagePtr<T> img, double radius = 1.5) {
        auto minTree = createMinTree(img, radius);
        auto maxTree = createMaxTree(img, radius);
        sdrt::SelfDualResidualTreeBuilder<T> builder(radius);
        builder.build(img, std::move(minTree), std::move(maxTree));
        auto topology = MorphologicalTree(tag(), builder.getNodeParent(), builder.getProperPartOwner(), builder.getRoot(), builder.getRows(), builder.getCols());
        return WeightedMorphologicalTree<T>(tag(), std::move(topology), detail::AltitudeInput<T>{builder.getAltitude(), detail::AltitudeDomain::InternalNodeSlots});
    }
};

} // namespace mmcfilters
