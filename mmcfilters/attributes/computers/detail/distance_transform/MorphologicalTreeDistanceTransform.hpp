#pragma once

#include "ExactSquaredEuclideanDistanceTransform2D.hpp"
#include "MorphologicalTreeContourScheduler.hpp"
#include "MorphologicalTreeRegionIndex.hpp"
#include "NodeDistanceFieldProvider.hpp"
#include "../../../../trees/MorphologicalTree.hpp"
#include "../../../../utils/Common.hpp"

#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

namespace mmcfilters::attributes::computers::detail::distance_transform {

/**
 * @brief Borrowed semantic view of one node's computed distance transform.
 *
 * The view is valid only during the consumer callback issued by
 * `MorphologicalTreeDistanceTransform::forEachNode`. It exposes translated
 * per-pixel distance labels without leaking the separable lower-envelope
 * workspace or hierarchical contour schedule.
 */
class NodeDistanceTransformFrame {
  public:
    NodeDistanceTransformFrame(const NodeDistanceTransformFrame&) = delete;
    NodeDistanceTransformFrame& operator=(const NodeDistanceTransformFrame&) = delete;
    NodeDistanceTransformFrame(NodeDistanceTransformFrame&&) = delete;
    NodeDistanceTransformFrame& operator=(NodeDistanceTransformFrame&&) = delete;

    /**
     * @brief Returns the node represented by this frame.
     */
    [[nodiscard]] NodeId node() const noexcept { return node_; }

    /**
     * @brief Returns the number of pixels in the active node support.
     */
    [[nodiscard]] std::size_t supportSize() const noexcept { return supportPixels_.size(); }

    /**
     * @brief Returns the minimal global-coordinate support box.
     */
    [[nodiscard]] const RegionBox2D& boundingBox() const noexcept { return box_; }

    /**
     * @brief Feeds any number of reducers through one support enumeration.
     *
     * The squared-distance label is read once per support pixel and fanned out
     * to every reducer. This is the production extension point for additional
     * attributes derived from the same node EDT.
     */
    template <NodeDistanceTransformReducer... Reducers>
        requires(sizeof...(Reducers) > 0)
    void reduceSamples(Reducers&... reducers) const {
        (reducers.beginNode(node_), ...);
        for (PixelId pixel : supportPixels_) {
            const SquaredDistance squaredDistance = transform_.establishedSquaredDistance(transformBox_.localPixel(pixel, globalColumns_));
            (reducers.consumeSample(pixel, squaredDistance), ...);
        }
        (reducers.endNode(node_), ...);
    }

    /**
     * @brief Convenience callback over the samples of one borrowed frame.
     */
    template <class Consumer> void forEachSample(Consumer&& consumer) const {
        for (PixelId pixel : supportPixels_) {
            consumer(pixel, transform_.establishedSquaredDistance(transformBox_.localPixel(pixel, globalColumns_)));
        }
    }

  private:
    friend class MorphologicalTreeDistanceTransform;

    NodeDistanceTransformFrame(NodeId node, std::span<const PixelId> supportPixels, RegionBox2D box, RegionBox2D transformBox, int globalColumns,
                               const ExactSquaredEuclideanDistanceTransform2D& transform)
        : node_(node), supportPixels_(supportPixels), box_(box), transformBox_(transformBox), globalColumns_(globalColumns), transform_(transform) {}

    NodeId node_ = InvalidNode;
    std::span<const PixelId> supportPixels_;
    RegionBox2D box_{};
    RegionBox2D transformBox_{};
    int globalColumns_ = 0;
    const ExactSquaredEuclideanDistanceTransform2D& transform_;
};

/**
 * @brief Uniform node-wise exact EDT backend for every 2D MorphologicalTree.
 *
 * Every node is evaluated from the same support-driven procedure. Proper parts
 * are indexed once, foreground 4-neighbour contours are scheduled bottom-up
 * from compact local changes, and one reusable exact separable EDT workspace is
 * initialized once per heavy path and updated by exact dirty-line transforms.
 * Tree kind, altitude, and construction adjacency never select another
 * implementation.
 */
class MorphologicalTreeDistanceTransform {
  public:
    /**
     * @brief Declares that every emitted sample is an exact squared distance.
     */
    static constexpr DistanceFieldAccuracy accuracy = DistanceFieldAccuracy::Exact;

    /**
     * @brief Computes the exact distance-transform frame of one live node.
     *
     * The hierarchy support index and shared boundary-lifetime predicate are
     * built once. Only the requested support bounding box receives a separable
     * EDT, so this path avoids transforming unrelated live nodes. The borrowed
     * frame remains valid only during the synchronous consumer callback.
     *
     * @tparam Consumer Synchronous callback accepting a `NodeDistanceTransformFrame`.
     * @param tree Established tree with a non-empty regular 2D domain.
     * @param node Live node whose field is requested.
     * @param consumer Callback invoked exactly once with the requested frame.
     */
    template <class Consumer> static void forNode(const MorphologicalTree& tree, NodeId node, Consumer&& consumer) {
        const MorphologicalTreeRegionIndex regions(tree);
        const GridDomain2D& domain = ::mmcfilters::detail::CommittedTreeAccess::gridDomain2D(tree);
        if (!tree.isAlive(node)) {
            throw std::out_of_range("Single-node distance-transform computation requires a live node id.");
        }

        const std::size_t mutationVersion = tree.getMutationVersion();
        const contours::detail::MorphologicalTreeBoundaryLifetimeIndex boundaries(tree);
        const std::span<const PixelId> supportPixels = regions.establishedSupport(node);
        const RegionBox2D box = regions.establishedBoundingBox(node);

        std::vector<PixelId> localContourPixels;
        localContourPixels.reserve(supportPixels.size());
        for (PixelId pixel : supportPixels) {
            if (boundaries.establishedIsBoundaryAt(pixel, node)) {
                localContourPixels.push_back(box.localPixel(pixel, domain.columns));
            }
        }
        if (localContourPixels.empty()) {
            throw std::logic_error("Single-node distance-transform computation found an empty foreground contour.");
        }

        ExactSquaredEuclideanDistanceTransform2D transform(box.rows(), box.columns());
        transform.computeEstablished(localContourPixels);

        const NodeDistanceTransformFrame frame(node, supportPixels, box, box, domain.columns, transform);
        consumer(frame);
        tree.requireMutationVersion(mutationVersion, "MorphologicalTreeDistanceTransform::forNode callback");
    }

    /**
     * @brief Computes one independent distance-transform frame per live node.
     *
     * @param tree Established tree with a non-empty regular 2D domain.
     * @param consumer Callback invoked synchronously with each borrowed frame.
     */
    template <class Consumer> static void forEachNode(const MorphologicalTree& tree, Consumer&& consumer) {
        MorphologicalTreeRegionIndex regions(tree);
        const GridDomain2D& domain = ::mmcfilters::detail::CommittedTreeAccess::gridDomain2D(tree);
        ExactSquaredEuclideanDistanceTransform2D transform(1, 1);
        std::vector<PixelId> localContourPixels;
        std::vector<PixelId> localAdditionPixels;
        std::vector<PixelId> localRemovalPixels;
        RegionBox2D transformBox{};

        MorphologicalTreeContourScheduler::forEachEstablishedNode(tree, regions, [&](const NodeContourFrame& contourFrame) {
            const NodeId node = contourFrame.node();
            const std::span<const PixelId> supportPixels = regions.establishedSupport(node);
            const RegionBox2D box = regions.establishedBoundingBox(node);

            if (contourFrame.startsHeavyPath()) {
                transformBox = regions.establishedBoundingBox(contourFrame.heavyPathTop());
                localContourPixels.clear();
                localContourPixels.reserve(contourFrame.pixels().size());
                for (PixelId pixel : contourFrame.pixels()) {
                    localContourPixels.push_back(transformBox.localPixel(pixel, domain.columns));
                }
                if (localContourPixels.empty()) {
                    throw std::logic_error("Distance-transform computation found a live node with an empty foreground contour.");
                }

                transform.resetDomain(transformBox.rows(), transformBox.columns());
                transform.computeEstablished(localContourPixels);
            } else {
                localAdditionPixels.clear();
                localRemovalPixels.clear();
                localAdditionPixels.reserve(contourFrame.additionsFromHeavyChild().size());
                localRemovalPixels.reserve(contourFrame.removalsFromHeavyChild().size());
                for (PixelId pixel : contourFrame.additionsFromHeavyChild()) {
                    localAdditionPixels.push_back(transformBox.localPixel(pixel, domain.columns));
                }
                for (PixelId pixel : contourFrame.removalsFromHeavyChild()) {
                    localRemovalPixels.push_back(transformBox.localPixel(pixel, domain.columns));
                }
                transform.applyEstablishedSiteDelta(localAdditionPixels, localRemovalPixels);
            }

            const NodeDistanceTransformFrame frame(node, supportPixels, box, transformBox, domain.columns, transform);
            consumer(frame);
        });
    }

    /**
     * @brief Runs multiple node reducers over one EDT and one support scan.
     *
     * Adding a reducer increases only the constant fan-out per sample; it does
     * not trigger another contour reconstruction, transform, or support walk.
     */
    template <NodeDistanceTransformReducer... Reducers>
        requires(sizeof...(Reducers) > 0)
    static void reduce(const MorphologicalTree& tree, Reducers&... reducers) {
        forEachNode(tree, [&](const NodeDistanceTransformFrame& frame) { frame.reduceSamples(reducers...); });
    }
};

/**
 * @brief Semantic name used by consumers that require exact distance fields.
 */
using ExactNodeDistanceFieldProvider = MorphologicalTreeDistanceTransform;

} // namespace mmcfilters::attributes::computers::detail::distance_transform
