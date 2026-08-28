#pragma once

#include "MorphologicalTree.hpp"
#include "NativeHierarchy.hpp"
#include "TreeOfShapesProducer.hpp"
#include "TreeAltitudeAlgorithms.hpp"
#include "ValuedMorphologicalTree.hpp"
#include "detail/ComponentTreeProducerDetail.hpp"
#include "detail/HigraHierarchyAdapterDetail.hpp"
#include "detail/MorphologicalTreeConstructionTag.hpp"
#include "detail/NativeHierarchyValidationDetail.hpp"
#include "sdrt/SaturatedResidualTreeBuilder.hpp"
#include "sdrt/ResidualEvolution.hpp"
#include "sdrt/UnrestrictedResidualTreeBuilder.hpp"
#include "../utils/Image.hpp"
#include "../utils/Contract.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace mmcfilters {

/**
 * @brief Public construction facade for all high-level morphological trees.
 *
 * `MorphologicalTreeFactory` is the single public entry point for creating
 * tree instances from images and static Higra hierarchies. The factory owns
 * only construction orchestration; it does not
 * store mutable build state and every construction method returns a fully
 * initialized valued-tree owner.
 *
 * Design contract:
 *
 * - builders whose nodes may have an empty proper part supply
 *   explicit node altitudes together with their native topology;
 * - Higra imports preserve the external Higra node-id domain for altitude
 *   mapping, then expose the result through the internal dense `NodeId` domain;
 * - callers that only need topology should use `ValuedMorphologicalTree<std::uint8_t>::topology()`
 *   instead of constructing a bare `MorphologicalTree` directly;
 * - max/min image factories derive altitude type from the image type. `uint8_t`
 *   uses counting sort; other supported altitude values use comparison sort;
 * - static Higra imports can preserve a non-canonical altitude type because
 *   the altitude array itself supplies that type.
 *
 * Keeping construction in this factory prevents `MorphologicalTree` and
 * `ValuedMorphologicalTree<std::uint8_t>` from depending on public factory logic while still
 * allowing their tag-based constructors to remain available for header-only
 * composition inside this library.
 */
class MorphologicalTreeFactory {
  private:
    /**
     * @brief Validates image preconditions shared by component-tree factories.
     * @param img Input image.
     * @param context Factory name used in diagnostics.
     */
    template <AltitudeValue T> static void validateComponentTreeImage(const ImagePtr<T>& img, const char* context) {
        MMCFILTERS_CONTRACT_REQUIRE(img != nullptr, throw std::invalid_argument(std::string(context) + " requires a non-null image."));
        MMCFILTERS_CONTRACT_REQUIRE(img->getNumRows() > 0 && img->getNumColumns() > 0 && img->getSize() > 0,
                                    throw std::invalid_argument(std::string(context) + " requires a non-empty 2D image."));
        MMCFILTERS_CONTRACT_CHECKED_ONLY(TreeAltitudeAlgorithms::validateFiniteImageAltitudes(img, context));
    }

    /**
     * @brief Validates that a component-tree grid adjacency is connected.
     * @param adjacency Input grid adjacency.
     */
    static void validateConnectedComponentTreeAdjacency(const RegularGridAdjacency2D& adjacency) {
        const int numPixels = adjacency.getNumRows() * adjacency.getNumColumns();
        std::vector<std::uint8_t> visited(static_cast<std::size_t>(numPixels), 0);
        std::vector<PixelId> pending;
        pending.reserve(static_cast<std::size_t>(numPixels));
        pending.push_back(0);
        visited[0] = 1;
        int reached = 0;
        while (!pending.empty()) {
            const PixelId pixel = pending.back();
            pending.pop_back();
            ++reached;
            for (PixelId neighbor : adjacency.getNeighborIndices(pixel)) {
                if (visited[static_cast<std::size_t>(neighbor)] == 0) {
                    visited[static_cast<std::size_t>(neighbor)] = 1;
                    pending.push_back(neighbor);
                }
            }
        }
        if (reached != numPixels) {
            throw std::invalid_argument("Component-tree adjacency must induce one connected image-domain graph.");
        }
    }

    /**
     * @brief Preserves observable versioning from former linked constructors.
     *
     * Component-tree and Higra constructors historically called one linking
     * primitive per non-root node between two storage invalidations. Native
     * materialization links in bulk, so this policy restores the same final
     * opaque mutation token without replaying those per-node invalidations.
     */
    enum class MaterializedVersionPolicy { CanonicalNative, PreserveLinkedConstruction };

    /**
     * @brief Creates the internal construction capability token.
     *
     * The token restricts direct calls to low-level constructors. Public users
     * create trees through this factory, which is the only class allowed to
     * instantiate the tag.
     *
     * @return The created internal construction capability token.
     */
    static constexpr detail::MorphologicalTreeConstructionTag tag() noexcept { return detail::MorphologicalTreeConstructionTag{}; }

    /**
     * @brief Validates tree of shapes image.
     *
     * @param img Image.
     */
    static void validateTreeOfShapesImage(const ImageUInt8Ptr& img) {
        if (!img) {
            throw std::invalid_argument("MorphologicalTreeFactory::createTreeOfShapes requires a non-null image.");
        }
        if (img->getNumRows() <= 0 || img->getNumColumns() <= 0 || img->getSize() <= 0) {
            throw std::invalid_argument("MorphologicalTreeFactory::createTreeOfShapes requires a non-empty 2D image.");
        }
    }

    /**
     * @brief Resolves a domain-free canonical immersion against a source domain.
     *
     * A `CanonicalComplementaryGridImmersion` declares only the ordered
     * minimum/maximum pairing. Binding it to the source domain yields the
     * concrete `ComplementaryGridImmersion` that the result retains, so every
     * consumer of a retained convention observes explicit adjacencies.
     *
     * @param convention Convention whose immersion is resolved in place.
     * @param numRows Number of rows in the source domain.
     * @param numColumns Number of columns in the source domain.
     */
    static void resolveTopographicImmersion(TopographicConvention& convention, int numRows, int numColumns) {
        const auto* canonical = std::get_if<CanonicalComplementaryGridImmersion>(&convention.immersion);
        if (canonical == nullptr) {
            return;
        }
        constexpr double connectivity4Radius = 1.0;
        constexpr double connectivity8Radius = 1.5;
        const bool min4Max8 = canonical->pairing == ComplementaryPairing::Min4Max8;
        convention.immersion = ComplementaryGridImmersion{ComplementaryAdjacencies{
            RegularGridAdjacency2D(numRows, numColumns, min4Max8 ? connectivity4Radius : connectivity8Radius),
            RegularGridAdjacency2D(numRows, numColumns, min4Max8 ? connectivity8Radius : connectivity4Radius)}};
    }

    /**
     * @brief Consumes a proven owning native hierarchy.
     *
     * @param hierarchy Hierarchy data.
     * @param preservedExternalNodeIdOffset Node identifier.
     * @param versionPolicy Policy used to preserve or reset imported version metadata.
     * @return The materialized proven owning native hierarchy.
     */
    template <AltitudeValue T>
    [[nodiscard]] static ValuedMorphologicalTree<T>
    materializeNativeHierarchy(detail::ValidatedNativeHierarchy<T>&& hierarchy, std::optional<NodeId> preservedExternalNodeIdOffset = std::nullopt,
                               MaterializedVersionPolicy versionPolicy = MaterializedVersionPolicy::CanonicalNative) {
        auto storage = std::move(hierarchy).release();
        const std::size_t numNodes = storage.parent.size();
        MorphologicalTree topology(tag(), std::move(storage.parent), std::move(storage.smallestNodeMap), storage.root, storage.gridDomain2D,
                                   std::move(storage.semantics), std::move(storage.topologyProof));
        if (versionPolicy == MaterializedVersionPolicy::PreserveLinkedConstruction) {
            // One initial storage invalidation, one per non-root link, and one
            // final iterator invalidation produced `numNodes + 1`.
            topology.mutationVersion_ = numNodes + 1;
        }
        if (preservedExternalNodeIdOffset) {
            topology.preserveExternalNodeIdOffset(*preservedExternalNodeIdOffset);
        }
        return ValuedMorphologicalTree<T>(tag(), std::move(topology), std::move(storage.nodeAltitudes));
    }

  public:
    /**
     * @brief Builds a typed valuedTree max-tree from an image.
     *
     * @param img Non-null, non-empty image. Its row/column domain becomes the
     * pixel domain of the resulting tree, and its pixel type becomes the
     * altitude type of the returned owner.
     * @param radius Radius used by the component-tree adjacency relation.
     *
     * @return A `ValuedMorphologicalTree<T>` whose topology is a max-tree
     * and whose altitude buffer is indexed by internal `NodeId`.
     *
     * @throws std::invalid_argument if the image domain is invalid.
     */
    template <AltitudeValue T> [[nodiscard]] static ValuedMorphologicalTree<T> createMaxTree(ImagePtr<T> img, double radius = 1.5) {
        validateComponentTreeImage(img, "MorphologicalTreeFactory::createMaxTree");
        RegularGridAdjacency2D adjacency(img->getNumRows(), img->getNumColumns(), radius);
        MMCFILTERS_CONTRACT_REQUIRE(img->getSize() == 1 || radius >= 1.0,
                                    throw std::invalid_argument("Component-tree adjacency must induce one connected image-domain graph."));
        return materializeNativeHierarchy<T>(
            detail::kernel::ComponentTreeProducer<T>(detail::kernel::ComponentTreePolarity::MaxTree, std::move(adjacency)).build(img), std::nullopt,
            MaterializedVersionPolicy::PreserveLinkedConstruction);
    }

    /**
     * @brief Builds a max-tree with an explicit regular-grid 2D adjacency.
     *
     * The relation may be radius-, rectangle-, line-, or
     * structuring-element-based. Its grid dimensions must match `img`.
     *
     * @param img Image.
     * @param adjacency Adjacency relation.
     * @return The resulting max-tree with an explicit regular-grid 2D adjacency.
     */
    template <AltitudeValue T> [[nodiscard]] static ValuedMorphologicalTree<T> createMaxTree(ImagePtr<T> img, RegularGridAdjacency2D adjacency) {
        validateComponentTreeImage(img, "MorphologicalTreeFactory::createMaxTree");
        MMCFILTERS_CONTRACT_REQUIRE(adjacency.getNumRows() == img->getNumRows() && adjacency.getNumColumns() == img->getNumColumns(),
                                    throw std::invalid_argument("Component-tree adjacency domain must match the input image domain."));
        MMCFILTERS_CONTRACT_CHECKED_ONLY(validateConnectedComponentTreeAdjacency(adjacency));
        return materializeNativeHierarchy<T>(
            detail::kernel::ComponentTreeProducer<T>(detail::kernel::ComponentTreePolarity::MaxTree, std::move(adjacency)).build(img), std::nullopt,
            MaterializedVersionPolicy::PreserveLinkedConstruction);
    }

    /**
     * @brief Builds a typed valuedTree min-tree from an image.
     *
     * @param img Non-null, non-empty image. Its row/column domain becomes the
     * pixel domain of the resulting tree, and its pixel type becomes the
     * altitude type of the returned owner.
     * @param radius Radius used by the component-tree adjacency relation.
     *
     * @return A `ValuedMorphologicalTree<T>` whose topology is a min-tree
     * and whose altitude buffer is indexed by internal `NodeId`.
     *
     * @throws std::invalid_argument if the image domain is invalid.
     */
    template <AltitudeValue T> [[nodiscard]] static ValuedMorphologicalTree<T> createMinTree(ImagePtr<T> img, double radius = 1.5) {
        validateComponentTreeImage(img, "MorphologicalTreeFactory::createMinTree");
        RegularGridAdjacency2D adjacency(img->getNumRows(), img->getNumColumns(), radius);
        MMCFILTERS_CONTRACT_REQUIRE(img->getSize() == 1 || radius >= 1.0,
                                    throw std::invalid_argument("Component-tree adjacency must induce one connected image-domain graph."));
        return materializeNativeHierarchy<T>(
            detail::kernel::ComponentTreeProducer<T>(detail::kernel::ComponentTreePolarity::MinTree, std::move(adjacency)).build(img), std::nullopt,
            MaterializedVersionPolicy::PreserveLinkedConstruction);
    }

    /**
     * @brief Builds a min-tree with an explicit regular-grid 2D adjacency.
     *
     * @param img Image.
     * @param adjacency Adjacency relation.
     * @return The resulting min-tree with an explicit regular-grid 2D adjacency.
     */
    template <AltitudeValue T> [[nodiscard]] static ValuedMorphologicalTree<T> createMinTree(ImagePtr<T> img, RegularGridAdjacency2D adjacency) {
        validateComponentTreeImage(img, "MorphologicalTreeFactory::createMinTree");
        MMCFILTERS_CONTRACT_REQUIRE(adjacency.getNumRows() == img->getNumRows() && adjacency.getNumColumns() == img->getNumColumns(),
                                    throw std::invalid_argument("Component-tree adjacency domain must match the input image domain."));
        MMCFILTERS_CONTRACT_CHECKED_ONLY(validateConnectedComponentTreeAdjacency(adjacency));
        return materializeNativeHierarchy<T>(
            detail::kernel::ComponentTreeProducer<T>(detail::kernel::ComponentTreePolarity::MinTree, std::move(adjacency)).build(img), std::nullopt,
            MaterializedVersionPolicy::PreserveLinkedConstruction);
    }

    /**
     * @brief Builds the unrestricted residual tree.
     *
     * Current regional extrema from the synchronized max-tree and min-tree are
     * selected by the self-dual residual schedule. Both component trees use the
     * same symmetric connected adjacency.
     *
     * @param img Image processed by the operation.
     * @param adjacency Adjacency relation.
     * @param options Spatial-order options for unrestricted construction.
     * @return An unrestricted residual tree over `adjacency`.
     */
    template <AltitudeValue T>
    [[nodiscard]] static ValuedMorphologicalTree<T>
    createUnrestrictedResidualTree(ImagePtr<T> img, RegularGridAdjacency2D adjacency,
                                   sdrt::UnrestrictedResidualTreeOptions options = {}) {
        auto minTree = createMinTree(img, adjacency);
        auto maxTree = createMaxTree(img, adjacency);
        sdrt::UnrestrictedResidualTreeBuilder<T> builder(adjacency, options);
        builder.build(img, std::move(minTree), std::move(maxTree));
        return materializeNativeHierarchy<T>(
            std::move(builder).takeValidatedHierarchy(makeMorphologicalTreeSemantics(
                MorphologicalTreeKind::UnrestrictedResidualTree, SharedAdjacencyContext{std::move(adjacency)})));
    }

    /**
     * @brief Builds the unrestricted residual tree with a radius adjacency.
     *
     * @param img Image processed by the operation.
     * @param radius Radius of the regular-grid adjacency.
     * @param options Spatial-order options for unrestricted construction.
     * @return An unrestricted residual tree using the radius adjacency.
     */
    template <AltitudeValue T>
    [[nodiscard]] static ValuedMorphologicalTree<T>
    createUnrestrictedResidualTree(ImagePtr<T> img, double radius = 1.5, sdrt::UnrestrictedResidualTreeOptions options = {}) {
        TreeAltitudeAlgorithms::validateFiniteImageAltitudes(img, "MorphologicalTreeFactory::createUnrestrictedResidualTree image");
        RegularGridAdjacency2D adjacency(img->getNumRows(), img->getNumColumns(), radius);
        return createUnrestrictedResidualTree(std::move(img), std::move(adjacency), options);
    }

    /**
     * @brief Builds the saturated residual tree.
     *
     * Only current regional extrema whose complement is connected to
     * `infinityPixel` are eligible. The supplied adjacency must be symmetric,
     * connected on the image domain, and shared by both component trees.
     *
     * @param img Image processed by the operation.
     * @param adjacency Adjacency relation.
     * @param infinityPixel Row-major infinity pixel.
     * @param options Saturation and spatial-order policies.
     * @return A saturated residual tree relative to `infinityPixel`.
     */
    template <AltitudeValue T>
    [[nodiscard]] static ValuedMorphologicalTree<T>
    createSaturatedResidualTree(ImagePtr<T> img, RegularGridAdjacency2D adjacency, PixelId infinityPixel,
                                sdrt::SaturatedResidualTreeOptions options = {}) {
        auto minTree = createMinTree(img, adjacency);
        auto maxTree = createMaxTree(img, adjacency);
        sdrt::SaturatedResidualTreeBuilder<T> builder(adjacency, infinityPixel, options);
        builder.build(img, std::move(minTree), std::move(maxTree));
        return materializeNativeHierarchy<T>(
            std::move(builder).takeValidatedHierarchy(makeMorphologicalTreeSemantics(
                MorphologicalTreeKind::SaturatedResidualTree, SaturatedResidualContext{std::move(adjacency), infinityPixel})));
    }

    /**
     * @brief Builds the saturated residual tree with a radius adjacency.
     *
     * @param img Image processed by the operation.
     * @param infinityPixel Row-major infinity pixel.
     * @param radius Radius of the regular-grid adjacency.
     * @param options Saturation and spatial-order policies.
     * @return A saturated residual tree using the radius adjacency.
     */
    template <AltitudeValue T>
    [[nodiscard]] static ValuedMorphologicalTree<T>
    createSaturatedResidualTree(ImagePtr<T> img, PixelId infinityPixel, double radius = 1.5,
                                sdrt::SaturatedResidualTreeOptions options = {}) {
        TreeAltitudeAlgorithms::validateFiniteImageAltitudes(img, "MorphologicalTreeFactory::createSaturatedResidualTree image");
        RegularGridAdjacency2D adjacency(img->getNumRows(), img->getNumColumns(), radius);
        return createSaturatedResidualTree(std::move(img), std::move(adjacency), infinityPixel, options);
    }

    /**
     * @brief Builds a valued tree of shapes from an 8-bit image.
     *
     * The published altitude type must match the altitude encoding declared by
     * the convention. The default convention selects the canonical
     * minimum-4/maximum-8 complementary-grid immersion without domain padding,
     * uses infinity pixel zero, and publishes unchanged 8-bit source levels.
     * `ToSGrayLevel` altitudes are available under
     * `TopographicAltitudeEncoding::ExactDoubled`, which is the only encoding a
     * self-dual span immersion admits.
     *
     * @tparam Altitude Published node-altitude type.
     * @param img Non-null, non-empty image.
     * @param convention Complete topographic convention retained by the result.
     * @return A `ValuedMorphologicalTree<Altitude>` whose topology is a tree of
     * shapes and whose altitude buffer is indexed by internal `NodeId`.
     *
     * @throws std::invalid_argument if the image domain or ToS parameters are
     * rejected by the underlying builder, or if `Altitude` contradicts the
     * declared altitude encoding.
     */
    template <class Altitude = std::uint8_t>
    [[nodiscard]] static ValuedMorphologicalTree<Altitude> createTreeOfShapes(ImageUInt8Ptr img, TopographicConvention convention = {}) {
        validateTreeOfShapesImage(img);
        resolveTopographicImmersion(convention, img->getNumRows(), img->getNumColumns());

        TreeOfShapesProducer producer(convention);
        TreeOfShapesBuildResult<Altitude> result = producer.build<Altitude>(img);
        MorphologicalTreeSemantics semantics = makeMorphologicalTreeSemantics(MorphologicalTreeKind::TreeOfShapes, std::move(convention));
        return materializeNativeHierarchy<Altitude>(std::move(result).takeValidatedHierarchy(std::move(semantics)));
    }

    /**
     * @brief Materializes the common output contract of a native producer.
     *
     * @param hierarchy Hierarchy data.
     * @return The materialized common output contract of a native producer.
     */
    template <AltitudeValue T> [[nodiscard]] static ValuedMorphologicalTree<T> createFromNativeHierarchy(NativeHierarchyView<T> hierarchy) {
        return materializeNativeHierarchy<T>(detail::validateAndCopyNativeHierarchy<T>(std::move(hierarchy)));
    }

    /**
     * @brief Imports a valued connected-subset tree from native buffers.
     *
     * Node ids and pixel ids have independent domains. Consequently, any node,
     * including the root, may have an empty proper part as long as the complete
     * buffers encode one connected rooted hierarchy and every node has
     * non-empty support.
     *
     * @param parent Parent node indexed by internal node identifier.
     * @param smallestNodeMap Pixel-indexed inclusion-smallest node map.
     * @param nodeAltitudes Altitude data indexed by internal node identifier.
     * @param root Root node of the traversal.
     * @param semantics Hierarchy semantics validated by the operation.
     * @return The imported valued connected-subset tree.
     */
    template <AltitudeValue T>
    [[nodiscard]] static ValuedMorphologicalTree<T> createFromNativeTopology(std::span<const NodeId> parent, std::span<const NodeId> smallestNodeMap,
                                                                               std::span<const T> nodeAltitudes, NodeId root, MorphologicalTreeSemantics semantics) {
        return createFromNativeHierarchy<T>(NativeHierarchyView<T>{parent, smallestNodeMap, nodeAltitudes, root, std::nullopt, std::move(semantics)});
    }

    /**
     * @brief Imports a native hierarchy with a regular 2D pixel layout.
     *
     * The row/column shape is optional metadata: it enables reconstruction and
     * geometry-dependent algorithms but does not alter the topology contract.
     *
     * @param parent Parent node indexed by internal node identifier.
     * @param smallestNodeMap Pixel-indexed inclusion-smallest node map.
     * @param nodeAltitudes Altitude data indexed by internal node identifier.
     * @param root Root node of the traversal.
     * @param rows Number of rows in the domain.
     * @param columns Number of columns in the domain.
     * @param semantics Hierarchy semantics validated by the operation.
     * @return The imported native hierarchy with a regular 2D pixel layout.
     */
    template <AltitudeValue T>
    [[nodiscard]] static ValuedMorphologicalTree<T> createFromNativeTopology(std::span<const NodeId> parent, std::span<const NodeId> smallestNodeMap,
                                                                               std::span<const T> nodeAltitudes, NodeId root, int rows, int columns,
                                                                               MorphologicalTreeSemantics semantics) {
        return createFromNativeHierarchy<T>(
            NativeHierarchyView<T>{parent, smallestNodeMap, nodeAltitudes, root, GridDomain2D{rows, columns}, std::move(semantics)});
    }

    /**
     * @brief Imports a static Higra parent/altitude hierarchy.
     *
     * Higra arrays use one compact id domain containing leaves followed by
     * internal nodes. For an image domain with `rows * columns` leaves, leaf ids
     * are `[0, rows * columns)` and internal ids are `[rows * columns, parent.size())`.
     * Every leaf must point to an internal node, every internal node must point
     * to another internal node or to itself, and exactly one internal node must
     * be self-parented as the root.
     *
     * @param parent Parent array in the compact Higra id domain.
     * @param nodeAltitudes Altitudes in the same compact Higra id domain.
     * @param rows Number of rows in the image/pixel domain.
     * @param columns Number of columns in the image/pixel domain.
     * @param kind Semantic tree kind represented by the imported hierarchy.
     * @param adjacency Required for max-tree and min-tree imports. Optional for
     * tree-of-shapes imports, whose topology can be interpreted without a
     * component-tree adjacency relation.
     *
     * @return A typed valued tree whose internal altitude buffer has already
     * been remapped from Higra ids to internal `NodeId` slots.
     *
     * @throws std::invalid_argument if the hierarchy is malformed, if the
     * altitude domain does not match the parent domain, or if a max/min import
     * does not provide adjacency.
     */
    template <AltitudeValue T>
    [[nodiscard]] static ValuedMorphologicalTree<T> createFromHigraParent(std::span<const NodeId> parent, std::span<const T> nodeAltitudes, int rows,
                                                                            int columns, MorphologicalTreeKind kind,
                                                                            std::optional<RegularGridAdjacency2D> adjacency = std::nullopt) {
        auto imported = detail::adaptHigraHierarchy<T>(parent, nodeAltitudes, rows, columns, kind, std::move(adjacency));
        return materializeNativeHierarchy<T>(std::move(imported.hierarchy), imported.internalNodeOffset, MaterializedVersionPolicy::PreserveLinkedConstruction);
    }
};

} // namespace mmcfilters
