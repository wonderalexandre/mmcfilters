#pragma once

#include "AdaptiveSquaredEuclideanDIFT2D.hpp"
#include "DynamicDistanceFieldStatistics2D.hpp"
#include "MorphologicalTreePropagationEdgeIndex.hpp"
#include "MorphologicalTreeTopologicalLevelIndex.hpp"
#include "../distance_transform/MorphologicalTreeRegionIndex.hpp"
#include "../../../../contours/detail/MorphologicalTreeBoundaryLifetimeIndex.hpp"
#include "../../../../trees/MorphologicalTree.hpp"
#include "../../../../trees/detail/CommittedTreeAccess.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <concepts>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace mmcfilters::attributes::computers::detail::distance_transform_approx {

/**
 * @brief Production policy for the topology-general approximate DIFT.
 */
struct ProductionApproximateDiftExecutionPolicy {
    using Queue = SquaredDistanceBucketQueue;
    inline static constexpr bool validateInternalInvariants = false;
    inline static constexpr bool assumeOpenImpliesSupport = true;
    inline static constexpr bool useLinearRemovalOffsets = false;
    inline static constexpr bool eraseInvalidatedPixels = false;
    inline static constexpr bool allowComparableOwnerFastPath = true;
};

/**
 * @brief Fully checked reference policy used by diagnostics and tests.
 */
struct ReferenceApproximateDiftExecutionPolicy {
    using Queue = SquaredDistanceBucketQueue;
    inline static constexpr bool validateInternalInvariants = true;
    inline static constexpr bool assumeOpenImpliesSupport = false;
    inline static constexpr bool useLinearRemovalOffsets = false;
    inline static constexpr bool eraseInvalidatedPixels = true;
    inline static constexpr bool allowComparableOwnerFastPath = false;
};

/**
 * @brief Checked form of the production invariants for tests and diagnostics.
 */
struct AuditedProductionApproximateDiftExecutionPolicy {
    using Queue = SquaredDistanceBucketQueue;
    inline static constexpr bool validateInternalInvariants = true;
    inline static constexpr bool assumeOpenImpliesSupport = true;
    inline static constexpr bool useLinearRemovalOffsets = false;
    inline static constexpr bool eraseInvalidatedPixels = false;
    inline static constexpr bool allowComparableOwnerFastPath = true;
};

/**
 * @brief Uniform topology-level approximate EDT for every 2D tree kind.
 *
 * Every proper-part pixel is inserted exactly once into one global DIFT. Nodes
 * are processed in pairwise-incomparable bottom-up topological levels. Multiple
 * frontier regions coexist safely because production either proves every A8
 * endpoint-owner pair comparable or activates each edge only when its owner
 * LCA is processed. Foreground A4 contour events remain exact and
 * topology-only; only the adaptive IFT propagation is approximate.
 */
template <class PropagationPolicy, class ExecutionPolicy> class BasicMorphologicalTreeApproximateDistanceTransform {
    using MomentsObserver =
        std::conditional_t<ExecutionPolicy::validateInternalInvariants, DynamicDistanceFieldMoments2D, UncheckedDynamicDistanceFieldMoments2D>;
    template <bool TrackMoments, bool TrackHistogram, bool TrackSpatialMoments>
    using SelectedStatisticsObserver =
        BasicDynamicDistanceFieldStatistics2D<TrackMoments, TrackHistogram, TrackSpatialMoments, ExecutionPolicy::validateInternalInvariants>;
    using CompleteStatisticsObserver = SelectedStatisticsObserver<true, true, true>;
    using PlateauTracker = std::conditional_t<ExecutionPolicy::validateInternalInvariants, MaximumPlateauTracker2D, UncheckedMaximumPlateauTracker2D>;

  public:
    [[nodiscard]] static std::vector<ApproxSquaredDistance> computeMaxSquaredDistance(const MorphologicalTree& tree) {
        std::vector<ApproxSquaredDistance> values(static_cast<std::size_t>(tree.numInternalNodeSlots()), ApproxSquaredDistance{0});
        forEachNodeMaximum(tree, [&values](NodeId node, ApproxSquaredDistance value) { values[static_cast<std::size_t>(node)] = value; });
        return values;
    }

    /**
     * @brief Emits every approximate per-pixel node field synchronously.
     *
     * The consumer receives `(node, supportPixels, squaredDistanceAt)`. The
     * accessor and support span are borrowed and must not escape the callback.
     */
    template <class FieldConsumer> static void forEachNodeField(const MorphologicalTree& tree, FieldConsumer&& consumeField) {
        distance_transform::MorphologicalTreeRegionIndex regions(tree);
        forEachNodeWithPolicies<NoopDistanceFieldObserver2D, NoopMaximumPixelTracker2D>(
            tree,
            [&regions, &consumeField](NodeId node, const distance_transform::DistanceFieldExtremum&, const distance_transform::DistanceFieldMaximumPlateau&,
                                      const NoopDistanceFieldObserver2D&, PixelId, const auto& transform) {
                const auto squaredDistanceAt = [&transform](PixelId pixel) { return transform.squaredDistance(pixel); };
                consumeField(node, regions.support(node), squaredDistanceAt);
            });
    }

    /**
     * @brief Emits only the selected sufficient statistics for every live node.
     *
     * This policy entry point is used by attribute and profile materializers
     * that need a subset of moments, histogram, spatial moments, maximum
     * localization, or maximum-plateau geometry. Disabled summaries allocate
     * no per-pixel storage and receive no label-update work.
     *
     * The consumer receives `(node, extremum, plateau, observer,
     * representative)`. It may access only observer summaries enabled by the
     * corresponding template flags. Borrowed references must not escape the
     * callback.
     *
     * @tparam TrackMoments Retain scalar distance moments.
     * @tparam TrackHistogram Retain the sparse squared-distance histogram.
     * @tparam TrackSpatialMoments Retain distance-weighted spatial moments.
     * @tparam TrackLocation Retain the canonical maximum pixel.
     * @tparam TrackPlateau Retain complete maximum-plateau geometry.
     * @tparam StatisticsConsumer Synchronous node consumer.
     * @param tree Established regular 2D morphological tree.
     * @param consumeStatistics Consumer invoked once per live node.
     */
    template <bool TrackMoments, bool TrackHistogram, bool TrackSpatialMoments, bool TrackLocation = false, bool TrackPlateau = false, class StatisticsConsumer>
    static void forEachNodeSelectedStatistics(const MorphologicalTree& tree, StatisticsConsumer&& consumeStatistics) {
        static_assert(TrackMoments || TrackHistogram || TrackSpatialMoments, "Selected approximate statistics require at least one field summary.");
        static_assert(!(TrackLocation && TrackPlateau), "Maximum localization and plateau tracking are mutually exclusive policies.");
        using Observer = SelectedStatisticsObserver<TrackMoments, TrackHistogram, TrackSpatialMoments>;
        using MaximumTracker =
            std::conditional_t<TrackPlateau, PlateauTracker, std::conditional_t<TrackLocation, MaximumPixelTracker2D, NoopMaximumPixelTracker2D>>;
        forEachNodeWithPolicies<Observer, MaximumTracker>(tree, [&consumeStatistics](NodeId node, const distance_transform::DistanceFieldExtremum& extremum,
                                                                                     const distance_transform::DistanceFieldMaximumPlateau& plateau,
                                                                                     const Observer& observer, PixelId representative) {
            if constexpr (ExecutionPolicy::validateInternalInvariants) {
                const std::uint64_t supportSize = observer.supportSizeFor(representative);
                if constexpr (TrackMoments) {
                    if (observer.momentsFor(representative).count() != supportSize) {
                        throw std::logic_error("Selected approximate moments do not cover the completed node support.");
                    }
                }
                if constexpr (TrackHistogram) {
                    if (observer.histogramFor(representative).count() != supportSize) {
                        throw std::logic_error("Selected approximate histogram does not cover the completed node support.");
                    }
                }
                if constexpr (TrackSpatialMoments) {
                    if (observer.spatialMomentsFor(representative).count() != supportSize) {
                        throw std::logic_error("Selected approximate spatial moments do not cover the completed node support.");
                    }
                }
            }
            consumeStatistics(node, extremum, plateau, observer, representative);
        });
    }

    /**
     * @brief Produces the approximate local maximum for each live node.
     *
     * The provider owns the adaptive DIFT traversal and emits one sufficient
     * local statistic per node. The consumer decides how that statistic is
     * materialized, so attribute storage remains outside the algorithm.
     */
    template <class MaximumConsumer>
        requires std::invocable<MaximumConsumer&, NodeId, ApproxSquaredDistance>
    static void forEachNodeMaximum(const MorphologicalTree& tree, MaximumConsumer&& consumeMaximum) {
        forEachNodeWithPolicies<NoopDistanceFieldObserver2D, NoopMaximumPixelTracker2D>(
            tree,
            [&consumeMaximum](NodeId node, const distance_transform::DistanceFieldExtremum& extremum, const distance_transform::DistanceFieldMaximumPlateau&,
                              const NoopDistanceFieldObserver2D&, PixelId) { consumeMaximum(node, extremum.squaredDistance); });
    }

    /**
     * @brief Produces the approximate Bedt maximum and its canonical attaining pixel.
     */
    template <class ExtremumConsumer>
        requires std::invocable<ExtremumConsumer&, NodeId, const distance_transform::DistanceFieldExtremum&>
    static void forEachNodeExtremum(const MorphologicalTree& tree, ExtremumConsumer&& consumeExtremum) {
        forEachNodeWithPolicies<NoopDistanceFieldObserver2D, MaximumPixelTracker2D>(
            tree,
            [&consumeExtremum](NodeId node, const distance_transform::DistanceFieldExtremum& extremum, const distance_transform::DistanceFieldMaximumPlateau&,
                               const NoopDistanceFieldObserver2D&, PixelId) { consumeExtremum(node, extremum); });
    }

    /**
     * @brief Produces the complete approximate maximum-distance plateau for each live node.
     */
    template <class PlateauConsumer>
        requires std::invocable<PlateauConsumer&, NodeId, const distance_transform::DistanceFieldMaximumPlateau&>
    static void forEachNodePlateau(const MorphologicalTree& tree, PlateauConsumer&& consumePlateau) {
        forEachNodeWithPolicies<NoopDistanceFieldObserver2D, PlateauTracker>(
            tree,
            [&consumePlateau](NodeId node, const distance_transform::DistanceFieldExtremum&, const distance_transform::DistanceFieldMaximumPlateau& plateau,
                              const NoopDistanceFieldObserver2D&, PixelId) { consumePlateau(node, plateau); });
    }

    /**
     * @brief Produces the Bedt maximum and additive moments of every approximate node field.
     */
    template <class SummaryConsumer>
        requires std::invocable<SummaryConsumer&, NodeId, ApproxSquaredDistance, const distance_transform::DistanceFieldMoments&>
    static void forEachNodeSummary(const MorphologicalTree& tree, SummaryConsumer&& consumeSummary) {
        forEachNodeWithPolicies<MomentsObserver, NoopMaximumPixelTracker2D>(
            tree, [&consumeSummary](NodeId node, const distance_transform::DistanceFieldExtremum& extremum,
                                    const distance_transform::DistanceFieldMaximumPlateau&, const MomentsObserver& observer, PixelId representative) {
                const distance_transform::DistanceFieldMoments& moments = observer.momentsFor(representative);
                if (moments.count() != observer.supportSizeFor(representative)) {
                    throw std::logic_error("Topology-level approximate DIFT completed a node with non-finite distance labels.");
                }
                consumeSummary(node, extremum.squaredDistance, moments);
            });
    }

    /**
     * @brief Produces Bedt localization together with additive approximate field moments.
     */
    template <class SummaryConsumer>
        requires std::invocable<SummaryConsumer&, NodeId, const distance_transform::DistanceFieldExtremum&, const distance_transform::DistanceFieldMoments&>
    static void forEachNodeSummaryAndExtremum(const MorphologicalTree& tree, SummaryConsumer&& consumeSummary) {
        forEachNodeWithPolicies<MomentsObserver, MaximumPixelTracker2D>(
            tree, [&consumeSummary](NodeId node, const distance_transform::DistanceFieldExtremum& extremum,
                                    const distance_transform::DistanceFieldMaximumPlateau&, const MomentsObserver& observer, PixelId representative) {
                const distance_transform::DistanceFieldMoments& moments = observer.momentsFor(representative);
                if (moments.count() != observer.supportSizeFor(representative)) {
                    throw std::logic_error("Topology-level approximate DIFT completed a localized node with non-finite distance labels.");
                }
                consumeSummary(node, extremum, moments);
            });
    }

    /**
     * @brief Produces complete Bedt plateau geometry together with additive approximate field moments.
     */
    template <class SummaryConsumer>
        requires std::invocable<SummaryConsumer&, NodeId, const distance_transform::DistanceFieldMaximumPlateau&,
                                const distance_transform::DistanceFieldMoments&>
    static void forEachNodeSummaryAndPlateau(const MorphologicalTree& tree, SummaryConsumer&& consumeSummary) {
        forEachNodeWithPolicies<MomentsObserver, PlateauTracker>(tree, [&consumeSummary](NodeId node, const distance_transform::DistanceFieldExtremum&,
                                                                                         const distance_transform::DistanceFieldMaximumPlateau& plateau,
                                                                                         const MomentsObserver& observer, PixelId representative) {
            const distance_transform::DistanceFieldMoments& moments = observer.momentsFor(representative);
            if (moments.count() != observer.supportSizeFor(representative)) {
                throw std::logic_error("Topology-level approximate DIFT completed a plateau node with non-finite distance labels.");
            }
            consumeSummary(node, plateau, moments);
        });
    }

    /**
     * @brief Produces complete moments, histogram, and spatial summaries for every node field.
     */
    template <class SummaryConsumer>
        requires std::invocable<SummaryConsumer&, NodeId, ApproxSquaredDistance, const distance_transform::DistanceFieldMoments&,
                                const distance_transform::DistanceFieldHistogram&, const distance_transform::DistanceWeightedSpatialMoments&>
    static void forEachNodeCompleteStatistics(const MorphologicalTree& tree, SummaryConsumer&& consumeSummary) {
        forEachNodeWithPolicies<CompleteStatisticsObserver, NoopMaximumPixelTracker2D>(
            tree,
            [&consumeSummary](NodeId node, const distance_transform::DistanceFieldExtremum& extremum, const distance_transform::DistanceFieldMaximumPlateau&,
                              const CompleteStatisticsObserver& observer, PixelId representative) {
                validateCompleteStatistics(observer, representative);
                consumeSummary(node, extremum.squaredDistance, observer.momentsFor(representative), observer.histogramFor(representative),
                               observer.spatialMomentsFor(representative));
            });
    }

    /**
     * @brief Produces complete field statistics and the canonical maximum location.
     */
    template <class SummaryConsumer>
        requires std::invocable<SummaryConsumer&, NodeId, const distance_transform::DistanceFieldExtremum&, const distance_transform::DistanceFieldMoments&,
                                const distance_transform::DistanceFieldHistogram&, const distance_transform::DistanceWeightedSpatialMoments&>
    static void forEachNodeCompleteStatisticsAndExtremum(const MorphologicalTree& tree, SummaryConsumer&& consumeSummary) {
        forEachNodeWithPolicies<CompleteStatisticsObserver, MaximumPixelTracker2D>(
            tree,
            [&consumeSummary](NodeId node, const distance_transform::DistanceFieldExtremum& extremum, const distance_transform::DistanceFieldMaximumPlateau&,
                              const CompleteStatisticsObserver& observer, PixelId representative) {
                validateCompleteStatistics(observer, representative);
                consumeSummary(node, extremum, observer.momentsFor(representative), observer.histogramFor(representative),
                               observer.spatialMomentsFor(representative));
            });
    }

    /**
     * @brief Produces complete field statistics and maximum-plateau geometry.
     */
    template <class SummaryConsumer>
        requires std::invocable<SummaryConsumer&, NodeId, const distance_transform::DistanceFieldMaximumPlateau&,
                                const distance_transform::DistanceFieldMoments&, const distance_transform::DistanceFieldHistogram&,
                                const distance_transform::DistanceWeightedSpatialMoments&>
    static void forEachNodeCompleteStatisticsAndPlateau(const MorphologicalTree& tree, SummaryConsumer&& consumeSummary) {
        forEachNodeWithPolicies<CompleteStatisticsObserver, PlateauTracker>(
            tree,
            [&consumeSummary](NodeId node, const distance_transform::DistanceFieldExtremum&, const distance_transform::DistanceFieldMaximumPlateau& plateau,
                              const CompleteStatisticsObserver& observer, PixelId representative) {
                validateCompleteStatistics(observer, representative);
                consumeSummary(node, plateau, observer.momentsFor(representative), observer.histogramFor(representative),
                               observer.spatialMomentsFor(representative));
            });
    }

  private:
    static void validateCompleteStatistics(const CompleteStatisticsObserver& observer, PixelId representative) {
        if constexpr (ExecutionPolicy::validateInternalInvariants) {
            const std::uint64_t supportSize = observer.supportSizeFor(representative);
            if (observer.momentsFor(representative).count() != supportSize || observer.histogramFor(representative).count() != supportSize ||
                observer.spatialMomentsFor(representative).count() != supportSize) {
                throw std::logic_error("Topology-level approximate DIFT completed a node with inconsistent complete field statistics.");
            }
        }
    }

    template <class Observer, class MaximumTracker, class NodeConsumer>
    static void forEachNodeWithPolicies(const MorphologicalTree& tree, NodeConsumer&& consumeNode) {
        const GridDomain2D domain = requireDomain(tree);
        if constexpr (ExecutionPolicy::allowComparableOwnerFastPath) {
            if (allA8EdgeOwnersComparable(tree, domain)) {
                forEachNodeWithPoliciesImpl<false, Observer, MaximumTracker>(tree, domain, consumeNode);
                return;
            }
        }
        forEachNodeWithPoliciesImpl<true, Observer, MaximumTracker>(tree, domain, consumeNode);
    }

    template <bool CheckActiveEdges, class Observer, class MaximumTracker, class NodeConsumer>
    static void forEachNodeWithPoliciesImpl(const MorphologicalTree& tree, GridDomain2D domain, NodeConsumer& consumeNode) {
        using Transform = EdtDIFT2D<PropagationPolicy, Observer, MaximumTracker, typename ExecutionPolicy::Queue, ApproxSquaredDistance, CheckActiveEdges,
                                    ExecutionPolicy::assumeOpenImpliesSupport, ExecutionPolicy::useLinearRemovalOffsets,
                                    ExecutionPolicy::validateInternalInvariants, ExecutionPolicy::eraseInvalidatedPixels>;
        const std::size_t mutationVersion = tree.getMutationVersion();
        const contours::detail::MorphologicalTreeBoundaryLifetimeIndex boundaries(tree);
        const MorphologicalTreeTopologicalLevelIndex levels(tree);
        std::unique_ptr<MorphologicalTreePropagationEdgeIndex> propagationEdges;
        if constexpr (CheckActiveEdges) {
            propagationEdges = std::make_unique<MorphologicalTreePropagationEdgeIndex>(tree);
        }

        Transform transform(domain.rows, domain.columns,
                            CheckActiveEdges ? PropagationEdgeInitialization::NoEdges : PropagationEdgeInitialization::AllDomainEdges);

        const std::size_t numSlots = static_cast<std::size_t>(tree.numInternalNodeSlots());
        std::vector<std::unique_ptr<std::vector<PixelId>>> contoursByNode(numSlots);
        std::vector<PixelId> contourPositions(static_cast<std::size_t>(tree.numPixels()), InvalidPixel);
        std::vector<std::uint8_t> additionMarks(static_cast<std::size_t>(tree.numPixels()), std::uint8_t{0});
        std::vector<PixelId> newlyOpenedPixels;
        newlyOpenedPixels.reserve(static_cast<std::size_t>(tree.numPixels()));

        for (int level = 0; level < levels.establishedNumLevels(); ++level) {
            const std::span<const NodeId> levelNodes = levels.establishedNodesAtLevel(level);
            newlyOpenedPixels.clear();

            for (NodeId node : levelNodes) {
                buildNodeContour(tree, node, boundaries, contoursByNode, contourPositions);
                transform.removeSeeds(boundaries.establishedRemovals(node));

                const std::span<const PixelId> additions = boundaries.establishedAdditions(node);
                markPixels(additions, additionMarks, tree.numPixels());
                std::size_t consumedAdditions = 0;
                for (PixelId pixel : ::mmcfilters::detail::CommittedTreeAccess::properParts(tree, node)) {
                    transform.addPixelToBinaryImage(pixel);
                    std::uint8_t& mark = additionMarks[static_cast<std::size_t>(pixel)];
                    if (mark != 0) {
                        transform.seed(pixel);
                        mark = 2;
                        ++consumedAdditions;
                    } else {
                        transform.open(pixel);
                        newlyOpenedPixels.push_back(pixel);
                    }
                }
                if constexpr (ExecutionPolicy::validateInternalInvariants) {
                    if (consumedAdditions != additions.size()) {
                        throw std::logic_error("Topology-level approximate DIFT found a contour addition outside its owner proper part.");
                    }
                }
                clearMarks(additions, additionMarks);
            }

            if constexpr (CheckActiveEdges) {
                for (NodeId node : levelNodes) {
                    for (const PropagationEdge2D& edge : propagationEdges->establishedActivations(node)) {
                        transform.activateEdge(edge.first, edge.second, edge.firstDirectionBit, edge.secondDirectionBit);
                    }
                }
            }
            for (PixelId pixel : newlyOpenedPixels) {
                transform.insertSupportNeighbours(pixel);
            }

            transform.run();
            for (NodeId node : levelNodes) {
                const auto& contour = contoursByNode[static_cast<std::size_t>(node)];
                if constexpr (ExecutionPolicy::validateInternalInvariants) {
                    if (!contour || contour->empty()) {
                        throw std::logic_error("Topology-level approximate DIFT produced an empty live-node contour.");
                    }
                }
                distance_transform::DistanceFieldExtremum extremum;
                distance_transform::DistanceFieldMaximumPlateau plateau;
                if constexpr (std::is_same_v<MaximumTracker, NoopMaximumPixelTracker2D>) {
                    extremum.squaredDistance = transform.maximumRootDistance(*contour);
                } else if constexpr (std::is_same_v<MaximumTracker, PlateauTracker>) {
                    plateau = transform.maximumRootPlateau(*contour);
                    extremum = distance_transform::DistanceFieldExtremum{plateau.squaredDistance, plateau.pixel};
                } else {
                    extremum = transform.maximumRootExtremum(*contour);
                }
                if constexpr (std::invocable<NodeConsumer&, NodeId, const distance_transform::DistanceFieldExtremum&,
                                             const distance_transform::DistanceFieldMaximumPlateau&, const Observer&, PixelId, const Transform&>) {
                    consumeNode(node, extremum, plateau, transform.observer(), contour->front(), transform);
                } else {
                    consumeNode(node, extremum, plateau, transform.observer(), contour->front());
                }
            }
            tree.requireMutationVersion(mutationVersion, "MorphologicalTreeApproximateDistanceTransform level update");
        }
    }

    /**
     * @brief Detects when the unrestricted Opt3 propagation graph is exact.
     *
     * If the inclusion-smallest owners of every A8 edge are comparable, the
     * edge is already internal as soon as both endpoints enter the support.
     * Consequently no sibling can communicate prematurely and all domain
     * edges may stay active, exactly as in the paper implementation. Trees that
     * do not satisfy this topology-only condition retain LCA-gated activation.
     */
    [[nodiscard]] static bool allA8EdgeOwnersComparable(const MorphologicalTree& tree, GridDomain2D domain) {
        const std::span<const NodeId> owners = tree.smallestNodeMap();
        constexpr std::array<std::pair<int, int>, 4> forwardOffsets{{{0, 1}, {1, -1}, {1, 0}, {1, 1}}};
        for (PixelId pixel = 0; pixel < tree.numPixels(); ++pixel) {
            const int row = pixel / domain.columns;
            const int column = pixel % domain.columns;
            const NodeId firstOwner = owners[static_cast<std::size_t>(pixel)];
            for (const auto [rowOffset, columnOffset] : forwardOffsets) {
                const int neighbourRow = row + rowOffset;
                const int neighbourColumn = column + columnOffset;
                if (neighbourRow < 0 || neighbourRow >= domain.rows || neighbourColumn < 0 || neighbourColumn >= domain.columns) {
                    continue;
                }
                const PixelId neighbour = neighbourRow * domain.columns + neighbourColumn;
                const NodeId secondOwner = owners[static_cast<std::size_t>(neighbour)];
                if (firstOwner == secondOwner) {
                    continue;
                }
                if (!::mmcfilters::detail::CommittedTreeAccess::isAncestor(tree, firstOwner, secondOwner) &&
                    !::mmcfilters::detail::CommittedTreeAccess::isAncestor(tree, secondOwner, firstOwner)) {
                    return false;
                }
            }
        }
        return true;
    }

    [[nodiscard]] static GridDomain2D requireDomain(const MorphologicalTree& tree) {
        tree.requireNotEditing("MorphologicalTreeApproximateDistanceTransform");
        if (!tree.hasGridDomain2D()) {
            throw std::invalid_argument("Approximate tree distance transform requires a regular 2D pixel domain.");
        }
        const GridDomain2D domain = *tree.gridDomain2D();
        if (domain.rows <= 0 || domain.columns <= 0 || domain.columns > std::numeric_limits<int>::max() / domain.rows ||
            domain.rows * domain.columns != tree.numPixels()) {
            throw std::invalid_argument("Approximate tree distance transform requires a consistent non-empty 2D pixel domain.");
        }
        return domain;
    }

    static void buildNodeContour(const MorphologicalTree& tree, NodeId node, const contours::detail::MorphologicalTreeBoundaryLifetimeIndex& boundaries,
                                 std::vector<std::unique_ptr<std::vector<PixelId>>>& contoursByNode, std::vector<PixelId>& positions) {
        NodeId storageChild = InvalidNode;
        std::size_t largestContour = 0;
        std::size_t combinedSize = boundaries.establishedAdditions(node).size();
        for (NodeId child : ::mmcfilters::detail::CommittedTreeAccess::children(tree, node)) {
            const auto& childContour = contoursByNode[static_cast<std::size_t>(child)];
            if constexpr (ExecutionPolicy::validateInternalInvariants) {
                if (!childContour) {
                    throw std::logic_error("Topology-level contour scheduling found an unmaterialized child contour.");
                }
            }
            combinedSize = sumSizes(combinedSize, childContour->size());
            if (storageChild == InvalidNode || childContour->size() > largestContour) {
                storageChild = child;
                largestContour = childContour->size();
            }
        }

        std::unique_ptr<std::vector<PixelId>> nodeContour;
        if (storageChild == InvalidNode) {
            nodeContour = std::make_unique<std::vector<PixelId>>();
        } else {
            nodeContour = std::move(contoursByNode[static_cast<std::size_t>(storageChild)]);
        }
        nodeContour->reserve(combinedSize);

        for (NodeId child : ::mmcfilters::detail::CommittedTreeAccess::children(tree, node)) {
            if (child == storageChild) {
                continue;
            }
            auto& childContour = contoursByNode[static_cast<std::size_t>(child)];
            for (PixelId pixel : *childContour) {
                positions[static_cast<std::size_t>(pixel)] = positionFromSize(nodeContour->size());
                nodeContour->push_back(pixel);
            }
            childContour.reset();
        }

        for (PixelId pixel : boundaries.establishedAdditions(node)) {
            if constexpr (ExecutionPolicy::validateInternalInvariants) {
                if (!tree.isPixel(pixel) || positions[static_cast<std::size_t>(pixel)] != InvalidPixel) {
                    throw std::logic_error("Topology-level contour scheduling found an invalid or duplicate local addition.");
                }
            }
            positions[static_cast<std::size_t>(pixel)] = positionFromSize(nodeContour->size());
            nodeContour->push_back(pixel);
        }

        for (PixelId pixel : boundaries.establishedRemovals(node)) {
            const PixelId position = positions[static_cast<std::size_t>(pixel)];
            if constexpr (ExecutionPolicy::validateInternalInvariants) {
                if (!tree.isPixel(pixel) || position == InvalidPixel || static_cast<std::size_t>(position) >= nodeContour->size() ||
                    (*nodeContour)[static_cast<std::size_t>(position)] != pixel) {
                    throw std::logic_error("Topology-level contour scheduling removed a pixel outside the active contour.");
                }
            }
            const PixelId movedPixel = nodeContour->back();
            (*nodeContour)[static_cast<std::size_t>(position)] = movedPixel;
            positions[static_cast<std::size_t>(movedPixel)] = position;
            nodeContour->pop_back();
            positions[static_cast<std::size_t>(pixel)] = InvalidPixel;
        }
        contoursByNode[static_cast<std::size_t>(node)] = std::move(nodeContour);
    }

    static void markPixels(std::span<const PixelId> pixels, std::vector<std::uint8_t>& marks, int numPixels) {
        for (PixelId pixel : pixels) {
            if constexpr (ExecutionPolicy::validateInternalInvariants) {
                if (pixel < 0 || pixel >= numPixels) {
                    throw std::logic_error("Topology-level approximate DIFT found an invalid contour pixel.");
                }
            }
            std::uint8_t& mark = marks[static_cast<std::size_t>(pixel)];
            if constexpr (ExecutionPolicy::validateInternalInvariants) {
                if (mark != 0) {
                    throw std::logic_error("Topology-level approximate DIFT found a duplicate contour addition.");
                }
            }
            mark = 1;
        }
    }

    static void clearMarks(std::span<const PixelId> pixels, std::vector<std::uint8_t>& marks) noexcept {
        for (PixelId pixel : pixels) {
            marks[static_cast<std::size_t>(pixel)] = 0;
        }
    }

    [[nodiscard]] static std::size_t sumSizes(std::size_t first, std::size_t second) {
        if constexpr (ExecutionPolicy::validateInternalInvariants) {
            if (second > std::numeric_limits<std::size_t>::max() - first) {
                throw std::overflow_error("Topology-level contour storage exceeds the addressable size domain.");
            }
        }
        return first + second;
    }

    [[nodiscard]] static PixelId positionFromSize(std::size_t position) {
        if constexpr (ExecutionPolicy::validateInternalInvariants) {
            if (position > static_cast<std::size_t>(std::numeric_limits<PixelId>::max())) {
                throw std::overflow_error("Topology-level contour position exceeds the PixelId domain.");
            }
        }
        return static_cast<PixelId>(position);
    }
};

/**
 * @brief The single approximate distance-transform provider used by production APIs.
 */
using MorphologicalTreeApproximateDistanceTransform =
    BasicMorphologicalTreeApproximateDistanceTransform<AdaptiveA8SquaredEuclideanPolicy, ProductionApproximateDiftExecutionPolicy>;

/**
 * @brief Fully checked reference reserved for tests and diagnostics.
 */
using ReferenceMorphologicalTreeApproximateDistanceTransform =
    BasicMorphologicalTreeApproximateDistanceTransform<AdaptiveA8SquaredEuclideanPolicy, ReferenceApproximateDiftExecutionPolicy>;

/**
 * @brief Audited production policy reserved for invariant tests.
 */
using AuditedProductionMorphologicalTreeApproximateDistanceTransform =
    BasicMorphologicalTreeApproximateDistanceTransform<AdaptiveA8SquaredEuclideanPolicy, AuditedProductionApproximateDiftExecutionPolicy>;

} // namespace mmcfilters::attributes::computers::detail::distance_transform_approx
