#pragma once

#include "DistanceTransformProfiles.hpp"
#include "NodeDistanceField.hpp"
#include "../attributes/computers/detail/distance_transform/DistanceFieldHistogram.hpp"
#include "../attributes/computers/detail/distance_transform/MorphologicalTreeDistanceTransform.hpp"
#include "../attributes/computers/detail/distance_transform_approx/MorphologicalTreeApproximateDistanceTransform.hpp"
#include "../trees/MorphologicalTree.hpp"

#include <cstdint>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcfilters {

/**
 * @brief Public topology-only computation of node distance fields and profiles.
 *
 * Exact single-node fields use a directed provider that transforms only the
 * requested support bounding box after shared topology indexes are built.
 * Approximate fields and both profile backends retain their hierarchy-wide
 * traversal semantics.
 */
class DistanceTransformComputation {
    /** @cond */
  private:
    using InternalHistogram = attributes::computers::detail::distance_transform::DistanceFieldHistogram;
    using ExactTransform = attributes::computers::detail::distance_transform::ExactNodeDistanceFieldProvider;
    using ApproximateTransform = attributes::computers::detail::distance_transform_approx::MorphologicalTreeApproximateDistanceTransform;
    /** @endcond */

  public:
    /**
     * @brief Materializes the compact field of one live node.
     * @param tree Regular 2D morphological-tree topology.
     * @param node Live internal node to materialize.
     * @param accuracy Exact or approximate numerical backend.
     * @return Self-contained support pixels and aligned squared distances.
     *
     * Exact accuracy transforms only the selected node. Approximate accuracy
     * traverses the shared hierarchy-wide DIFT and copies only this node.
     */
    [[nodiscard]] static NodeDistanceField computeNodeField(const MorphologicalTree& tree, NodeId node,
                                                            DistanceTransformAccuracy accuracy = DistanceTransformAccuracy::Exact) {
        requireTreeAndNode(tree, node);
        requireAccuracy(accuracy);
        std::vector<PixelId> supportPixels;
        std::vector<std::int64_t> squaredDistances;
        bool found = false;
        if (accuracy == DistanceTransformAccuracy::Exact) {
            static_cast<void>(ExactTransform::forNode(tree, node, [&](const auto& frame) {
                supportPixels.reserve(frame.supportSize());
                squaredDistances.reserve(frame.supportSize());
                frame.forEachSample([&](PixelId pixel, std::int64_t squaredDistance) {
                    supportPixels.push_back(pixel);
                    squaredDistances.push_back(squaredDistance);
                });
                found = true;
            }));
        } else {
            static_cast<void>(
                ApproximateTransform::forEachNodeField(tree, [&](NodeId currentNode, std::span<const PixelId> pixels, const auto& squaredDistanceAt) {
                    if (currentNode != node) {
                        return;
                    }
                    supportPixels.assign(pixels.begin(), pixels.end());
                    squaredDistances.reserve(pixels.size());
                    for (PixelId pixel : pixels) {
                        squaredDistances.push_back(squaredDistanceAt(pixel));
                    }
                    found = true;
                }));
        }
        if (!found || supportPixels.empty() || supportPixels.size() != squaredDistances.size()) {
            throw std::logic_error("Distance-transform computation did not materialize the requested live-node field.");
        }
        return NodeDistanceField(node, tree.numRows(), tree.numColumns(), accuracy, std::move(supportPixels), std::move(squaredDistances));
    }

    /**
     * @brief Computes sorted sparse histograms for every live node.
     * @param tree Regular 2D morphological-tree topology.
     * @param accuracy Exact or approximate numerical backend.
     * @return Self-contained CSR-like all-node sparse profiles.
     *
     * The approximate backend selects a histogram-only observer, avoiding
     * scalar- and spatial-moment storage and updates.
     */
    [[nodiscard]] static DistanceTransformProfiles computeProfiles(const MorphologicalTree& tree,
                                                                   DistanceTransformAccuracy accuracy = DistanceTransformAccuracy::Exact) {
        requireTree(tree);
        requireAccuracy(accuracy);
        std::vector<DistanceTransformHistogramBin> bins;
        std::vector<DistanceTransformProfiles::Span> spans(static_cast<std::size_t>(tree.numInternalNodeSlots()));
        std::vector<std::uint8_t> alive(static_cast<std::size_t>(tree.numInternalNodeSlots()), std::uint8_t{0});
        for (NodeId node : tree.aliveNodeIds()) {
            alive[static_cast<std::size_t>(node)] = 1;
        }

        const auto appendHistogram = [&](NodeId node, const InternalHistogram& histogram) {
            DistanceTransformProfiles::Span& span = spans[static_cast<std::size_t>(node)];
            span.offset = bins.size();
            span.size = histogram.levelCount();
            for (const auto& [squaredDistance, count] : histogram.bins()) {
                bins.push_back({squaredDistance, count});
            }
        };

        if (accuracy == DistanceTransformAccuracy::Exact) {
            static_cast<void>(ExactTransform::forEachNode(tree, [&](const auto& frame) {
                InternalHistogram histogram;
                frame.forEachSample([&](PixelId, std::int64_t squaredDistance) { histogram.add<false>(squaredDistance); });
                appendHistogram(frame.node(), histogram);
            }));
        } else {
            static_cast<void>(ApproximateTransform::template forEachNodeSelectedStatistics<false, true, false>(
                tree, [&](NodeId node, const auto&, const auto&, const auto& observer, PixelId representative) {
                    appendHistogram(node, observer.histogramFor(representative));
                }));
        }
        return DistanceTransformProfiles(accuracy, std::move(bins), std::move(spans), std::move(alive));
    }

    /** @cond */
  private:
    static void requireTree(const MorphologicalTree& tree) {
        tree.requireNotEditing("DistanceTransformComputation");
        if (!tree.hasGridDomain2D() || tree.numRows() <= 0 || tree.numColumns() <= 0) {
            throw std::invalid_argument("Distance-transform computation requires a non-empty regular 2D tree domain.");
        }
    }

    static void requireTreeAndNode(const MorphologicalTree& tree, NodeId node) {
        requireTree(tree);
        if (!tree.isAlive(node)) {
            throw std::out_of_range("Distance-transform field computation requires a live internal node id.");
        }
    }

    static void requireAccuracy(DistanceTransformAccuracy accuracy) {
        if (accuracy != DistanceTransformAccuracy::Exact && accuracy != DistanceTransformAccuracy::Approximate) {
            throw std::invalid_argument("Distance-transform computation requires Exact or Approximate accuracy.");
        }
    }
    /** @endcond */
};

} // namespace mmcfilters
