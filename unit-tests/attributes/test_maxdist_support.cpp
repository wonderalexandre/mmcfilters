#include "support/TestSupport.hpp"

#include "attributes/support/maxdist/component_tree_dift/Geometry.hpp"
#include "attributes/support/maxdist/component_tree_dift/PQueue.hpp"
#include "mmcfilters/attributes/computers/MaxDistExactComputer.hpp"
#include "mmcfilters/attributes/computers/MaxDistComputer.hpp"
#include "mmcfilters/attributes/computers/detail/distance_transform/ExactSquaredEuclideanDistanceTransform2D.hpp"
#include "mmcfilters/attributes/computers/detail/distance_transform/MorphologicalTreeContourScheduler.hpp"
#include "mmcfilters/attributes/computers/detail/distance_transform/MorphologicalTreeDistanceTransform.hpp"
#include "mmcfilters/attributes/computers/detail/distance_transform/MorphologicalTreeRegionIndex.hpp"
#include "mmcfilters/attributes/computers/detail/distance_transform_approx/AdaptiveSquaredEuclideanDIFT2D.hpp"
#include "mmcfilters/attributes/computers/detail/distance_transform_approx/SquaredDistanceBucketQueue.hpp"
#include "mmcfilters/attributes/computers/detail/distance_transform_approx/MorphologicalTreeApproximateDistanceTransform.hpp"
#include "mmcfilters/contours/detail/MorphologicalTreeBoundaryLifetimeIndex.hpp"
#include "attributes/support/maxdist/MaxDistOracles.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <random>
#include <type_traits>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests::maxdist_oracle::component_tree_dift;
using namespace mmcfilters::unit_tests;

using DistanceTransformFrame = attributes::computers::detail::distance_transform::NodeDistanceTransformFrame;
using ExactDistanceTransformWorkspace = attributes::computers::detail::distance_transform::ExactSquaredEuclideanDistanceTransform2D;
static_assert(!std::is_copy_constructible_v<DistanceTransformFrame>);
static_assert(!std::is_move_constructible_v<DistanceTransformFrame>);
static_assert(!std::is_copy_constructible_v<ExactDistanceTransformWorkspace>);
static_assert(!std::is_move_constructible_v<ExactDistanceTransformWorkspace>);

namespace {

using BoundaryLifetimeIndex = contours::detail::MorphologicalTreeBoundaryLifetimeIndex;
using ContourScheduler = attributes::computers::detail::distance_transform::MorphologicalTreeContourScheduler;
using RegionIndex = attributes::computers::detail::distance_transform::MorphologicalTreeRegionIndex;
using SquaredDistance = attributes::computers::detail::distance_transform::SquaredDistance;
using DistanceFieldExtremum = attributes::computers::detail::distance_transform::DistanceFieldExtremum;
using DistanceFieldMaximumPlateau = attributes::computers::detail::distance_transform::DistanceFieldMaximumPlateau;
using ExactProvider = attributes::computers::detail::distance_transform::ExactNodeDistanceFieldProvider;
using ApproximateDift = attributes::computers::detail::distance_transform_approx::EdtDIFT2D<>;
using ApproximateQueue = attributes::computers::detail::distance_transform_approx::SquaredDistanceBucketQueue;
using ApproximateTreeTransform = attributes::computers::detail::distance_transform_approx::ReferenceMorphologicalTreeApproximateDistanceTransform;
using AuditedFastApproximateTreeTransform =
    attributes::computers::detail::distance_transform_approx::AuditedProductionMorphologicalTreeApproximateDistanceTransform;
using HistogramOnlyObserver = attributes::computers::detail::distance_transform_approx::BasicDynamicDistanceFieldStatistics2D<false, true, false, false>;
using SpatialOnlyObserver = attributes::computers::detail::distance_transform_approx::BasicDynamicDistanceFieldStatistics2D<false, false, true, false>;
static_assert(!HistogramOnlyObserver::tracksMoments && HistogramOnlyObserver::tracksHistogram && !HistogramOnlyObserver::tracksSpatialMoments);
static_assert(!SpatialOnlyObserver::tracksMoments && !SpatialOnlyObserver::tracksHistogram && SpatialOnlyObserver::tracksSpatialMoments);

struct DenseMaximumReducer {
    explicit DenseMaximumReducer(int numNodeSlots) : values(static_cast<std::size_t>(numNodeSlots), SquaredDistance{-1}) {}

    void beginNode(NodeId) {
        ++beginCount;
        current = 0;
    }

    void consumeSample(PixelId, SquaredDistance squaredDistance) {
        ++sampleCount;
        current = std::max(current, squaredDistance);
    }

    void endNode(NodeId node) {
        ++endCount;
        values[static_cast<std::size_t>(node)] = current;
    }

    std::vector<SquaredDistance> values;
    SquaredDistance current = 0;
    std::size_t beginCount = 0;
    std::size_t sampleCount = 0;
    std::size_t endCount = 0;
};

struct DenseDistanceSumReducer {
    explicit DenseDistanceSumReducer(int numNodeSlots) : values(static_cast<std::size_t>(numNodeSlots), SquaredDistance{-1}) {}

    void beginNode(NodeId) {
        ++beginCount;
        current = 0;
    }

    void consumeSample(PixelId, SquaredDistance squaredDistance) {
        ++sampleCount;
        current += squaredDistance;
    }

    void endNode(NodeId node) {
        ++endCount;
        values[static_cast<std::size_t>(node)] = current;
    }

    std::vector<SquaredDistance> values;
    SquaredDistance current = 0;
    std::size_t beginCount = 0;
    std::size_t sampleCount = 0;
    std::size_t endCount = 0;
};

template <std::floating_point Real> struct ProjectedMaximumReducer {
    ProjectedMaximumReducer(std::span<Real> output, const AttributeNames& names) : output(output), names(names) {}

    void beginNode(NodeId) { current = 0; }
    void consumeSample(PixelId, SquaredDistance squaredDistance) { current = std::max(current, squaredDistance); }
    void endNode(NodeId node) { output[names.linearIndex(node, MaxSquaredDistExact)] = static_cast<Real>(current); }

    std::span<Real> output;
    const AttributeNames& names;
    SquaredDistance current = 0;
};

static_assert(attributes::computers::detail::distance_transform::NodeDistanceTransformReducer<DenseMaximumReducer>);
static_assert(attributes::computers::detail::distance_transform::NodeDistanceTransformReducer<DenseDistanceSumReducer>);
static_assert(attributes::computers::detail::distance_transform::NodeDistanceFieldProviderFor<ExactProvider, DenseMaximumReducer>);
static_assert(ExactProvider::accuracy == attributes::computers::detail::distance_transform::DistanceFieldAccuracy::Exact);

void requireRegionIndexMatchesTree(const MorphologicalTree& tree, const std::string& label) {
    RegionIndex regions(tree);
    requireEqual(regions.numIndexedPixels(), static_cast<std::size_t>(tree.numPixels()), label + " indexed pixel count");

    for (NodeId node : tree.aliveNodeIds()) {
        const std::span<const PixelId> indexedSupport = regions.support(node);
        const std::vector<PixelId> actual(indexedSupport.begin(), indexedSupport.end());
        const auto supportRange = tree.nodeSupport(node);
        const std::vector<PixelId> expected(supportRange.begin(), supportRange.end());
        requireVectorEqual(actual, expected, label + " support interval node " + std::to_string(node));

        int rowMin = std::numeric_limits<int>::max();
        int rowMax = -1;
        int columnMin = std::numeric_limits<int>::max();
        int columnMax = -1;
        const int columns = tree.gridDomain2D()->columns;
        for (PixelId pixel : expected) {
            rowMin = std::min(rowMin, pixel / columns);
            rowMax = std::max(rowMax, pixel / columns);
            columnMin = std::min(columnMin, pixel % columns);
            columnMax = std::max(columnMax, pixel % columns);
        }

        const auto& box = regions.boundingBox(node);
        requireEqual(box.rowMin, rowMin, label + " bbox row min node " + std::to_string(node));
        requireEqual(box.rowMax, rowMax, label + " bbox row max node " + std::to_string(node));
        requireEqual(box.columnMin, columnMin, label + " bbox column min node " + std::to_string(node));
        requireEqual(box.columnMax, columnMax, label + " bbox column max node " + std::to_string(node));

        const auto parentInterval = regions.supportInterval(node);
        std::vector<attributes::computers::detail::distance_transform::RegionSupportInterval> childIntervals;
        for (NodeId child : tree.children(node)) {
            const auto childInterval = regions.supportInterval(child);
            require(childInterval.begin >= parentInterval.begin && childInterval.end <= parentInterval.end,
                    label + " child support interval contained by parent");
            childIntervals.push_back(childInterval);
        }
        for (std::size_t lhs = 0; lhs < childIntervals.size(); ++lhs) {
            for (std::size_t rhs = lhs + 1; rhs < childIntervals.size(); ++rhs) {
                require(childIntervals[lhs].end <= childIntervals[rhs].begin || childIntervals[rhs].end <= childIntervals[lhs].begin,
                        label + " sibling support intervals are disjoint");
            }
        }
    }

    const auto rootInterval = regions.supportInterval(tree.root());
    requireEqual(rootInterval.begin, std::size_t{0}, label + " root support begin");
    requireEqual(rootInterval.end, static_cast<std::size_t>(tree.numPixels()), label + " root support end");
}

void requireBoundaryLifetimesAndScheduleMatchOracle(const MorphologicalTree& tree, const std::string& label) {
    BoundaryLifetimeIndex lifetimes(tree);
    RegionIndex regions(tree);
    NodeId previousScheduledNode = InvalidNode;
    std::vector<std::uint8_t> transitionMask(static_cast<std::size_t>(tree.numPixels()), std::uint8_t{0});
    ContourScheduler::forEachNode(tree, regions, [&](const auto& frame) {
        require(tree.isAncestor(frame.heavyPathTop(), frame.node()), label + " heavy-path top is an ancestor");
        if (frame.startsHeavyPath()) {
            requireEqual(frame.heavyChild(), InvalidNode, label + " heavy-path start has no heavy child");
            std::fill(transitionMask.begin(), transitionMask.end(), std::uint8_t{0});
            for (PixelId pixel : frame.pixels()) {
                transitionMask[static_cast<std::size_t>(pixel)] = 1;
            }
        } else {
            requireEqual(frame.heavyChild(), previousScheduledNode, label + " transition names the immediately preceding heavy child");
            require(previousScheduledNode != InvalidNode && tree.parent(previousScheduledNode) == frame.node(),
                    label + " heavy child is emitted immediately before its parent");
            for (PixelId pixel : frame.additionsFromHeavyChild()) {
                require(transitionMask[static_cast<std::size_t>(pixel)] == 0, label + " heavy transition adds an inactive site");
                transitionMask[static_cast<std::size_t>(pixel)] = 1;
            }
            for (PixelId pixel : frame.removalsFromHeavyChild()) {
                require(transitionMask[static_cast<std::size_t>(pixel)] != 0, label + " heavy transition removes an active site");
                transitionMask[static_cast<std::size_t>(pixel)] = 0;
            }
        }
        for (PixelId pixel = 0; pixel < tree.numPixels(); ++pixel) {
            const bool scheduledBoundary = std::find(frame.pixels().begin(), frame.pixels().end(), pixel) != frame.pixels().end();
            requireEqual(transitionMask[static_cast<std::size_t>(pixel)] != 0, scheduledBoundary, label + " heavy transition reconstructs the emitted contour");
        }

        const auto exact = maxdist_oracle::exactNodeSquaredDistanceTransform(tree, frame.node());
        std::vector<PixelId> actual(frame.pixels().begin(), frame.pixels().end());
        std::vector<PixelId> expected = exact.contourPixels;
        std::sort(actual.begin(), actual.end());
        std::sort(expected.begin(), expected.end());
        requireVectorEqual(actual, expected, label + " hierarchical contour node " + std::to_string(frame.node()));
        std::vector<std::uint8_t> expectedMask(static_cast<std::size_t>(tree.numPixels()), std::uint8_t{0});
        for (PixelId pixel : expected) {
            expectedMask[static_cast<std::size_t>(pixel)] = 1;
        }
        for (PixelId pixel = 0; pixel < tree.numPixels(); ++pixel) {
            requireEqual(lifetimes.isBoundaryAt(pixel, frame.node()), expectedMask[static_cast<std::size_t>(pixel)] != 0,
                         label + " boundary lifetime predicate node " + std::to_string(frame.node()));
        }
        previousScheduledNode = frame.node();
    });
}

void requireExactTransformMatchesBruteForce() {
    using attributes::computers::detail::distance_transform::ExactSquaredEuclideanDistanceTransform2D;
    using attributes::computers::detail::distance_transform::SquaredDistance;

    requireThrows<std::invalid_argument>([] { static_cast<void>(ExactSquaredEuclideanDistanceTransform2D(0, 2)); },
                                         "exact EDT must reject a zero-sized domain");
    requireThrows<std::invalid_argument>([] { static_cast<void>(ExactSquaredEuclideanDistanceTransform2D(-1, 2)); },
                                         "exact EDT must reject a negative domain extent");
    requireThrows<std::overflow_error>([] { static_cast<void>(ExactSquaredEuclideanDistanceTransform2D(50000, 50000)); },
                                       "exact EDT must reject a domain whose pixel count overflows");

    std::mt19937 generator(0x45445432U);
    for (int rows = 1; rows <= 8; ++rows) {
        for (int columns = 1; columns <= 9; ++columns) {
            ExactSquaredEuclideanDistanceTransform2D transform(rows, columns);
            std::uniform_int_distribution<int> siteDistribution(0, rows * columns - 1);
            for (int trial = 0; trial < 8; ++trial) {
                std::vector<PixelId> sites;
                const int numSites = 1 + trial % std::min(rows * columns, 7);
                for (int siteIndex = 0; siteIndex < numSites; ++siteIndex) {
                    sites.push_back(siteDistribution(generator));
                }
                transform.compute(sites);

                for (PixelId pixel = 0; pixel < rows * columns; ++pixel) {
                    const SquaredDistance row = pixel / columns;
                    const SquaredDistance column = pixel % columns;
                    SquaredDistance expected = std::numeric_limits<SquaredDistance>::max();
                    for (PixelId site : sites) {
                        const SquaredDistance deltaRow = row - site / columns;
                        const SquaredDistance deltaColumn = column - site % columns;
                        expected = std::min(expected, deltaRow * deltaRow + deltaColumn * deltaColumn);
                    }
                    requireEqual(transform.squaredDistance(pixel), expected, "exact separable EDT must match brute force on every pixel");
                }
            }
        }
    }

    ExactSquaredEuclideanDistanceTransform2D transform(2, 2);
    requireThrows<std::logic_error>([&] { static_cast<void>(transform.squaredDistance(0)); }, "exact EDT must reject a query before its first computation");
    const std::array<PixelId, 1> validSite{0};
    transform.compute(validSite);
    const SquaredDistance retainedDistance = transform.squaredDistance(3);
    requireThrows<std::invalid_argument>([&] { transform.compute(std::span<const PixelId>{}); }, "exact EDT must reject an empty boundary-site set");
    const std::array<PixelId, 1> invalidSite{4};
    requireThrows<std::out_of_range>([&] { transform.compute(invalidSite); }, "exact EDT must reject an out-of-domain boundary site");
    const std::array<PixelId, 1> negativeSite{-1};
    requireThrows<std::out_of_range>([&] { transform.compute(negativeSite); }, "exact EDT must reject a negative boundary site");
    requireEqual(transform.squaredDistance(3), retainedDistance, "invalid EDT input must preserve the previous valid result");
    requireThrows<std::out_of_range>([&] { static_cast<void>(transform.squaredDistance(4)); }, "exact EDT must reject an invalid query pixel");

    transform.resetDomain(1, 3);
    requireEqual(transform.rows(), 1, "reconfigured EDT row count");
    requireEqual(transform.columns(), 3, "reconfigured EDT column count");
    requireThrows<std::logic_error>([&] { static_cast<void>(transform.squaredDistance(0)); }, "reconfiguring EDT must invalidate the previous result");
    const std::array<PixelId, 1> middleSite{1};
    transform.compute(middleSite);
    requireEqual(transform.squaredDistance(0), SquaredDistance{1}, "reconfigured EDT left distance");
    requireEqual(transform.squaredDistance(1), SquaredDistance{0}, "reconfigured EDT site distance");
    requireEqual(transform.squaredDistance(2), SquaredDistance{1}, "reconfigured EDT right distance");
    requireThrows<std::invalid_argument>([&] { transform.resetDomain(0, 3); }, "invalid reconfiguration must be rejected");
    requireEqual(transform.squaredDistance(2), SquaredDistance{1}, "invalid reconfiguration must preserve the prior valid result");

    ExactSquaredEuclideanDistanceTransform2D transactional(2, 2);
    const std::array<PixelId, 1> transactionalInitialSite{0};
    transactional.compute(transactionalInitialSite);
    const auto captureTransactionalField = [&] {
        std::vector<SquaredDistance> field;
        field.reserve(4);
        for (PixelId pixel = 0; pixel < 4; ++pixel) {
            field.push_back(transactional.squaredDistance(pixel));
        }
        return field;
    };
    const auto requireTransactionalField = [&](const std::vector<SquaredDistance>& expected, const std::string& label) {
        requireVectorEqual(captureTransactionalField(), expected, label);
    };
    const std::vector<SquaredDistance> initialTransactionalField = captureTransactionalField();

    const std::array<PixelId, 1> invalidDeltaPixel{4};
    requireThrows<std::out_of_range>([&] { static_cast<void>(transactional.applySiteDelta(invalidDeltaPixel, {})); },
                                     "incremental EDT must reject an invalid addition pixel");
    requireTransactionalField(initialTransactionalField, "invalid addition preserves the incremental EDT field");
    requireThrows<std::out_of_range>([&] { static_cast<void>(transactional.applySiteDelta({}, invalidDeltaPixel)); },
                                     "incremental EDT must reject an invalid removal pixel");
    requireTransactionalField(initialTransactionalField, "invalid removal preserves the incremental EDT field");

    const std::array<PixelId, 2> duplicateBatchAdditions{1, 1};
    requireThrows<std::invalid_argument>([&] { static_cast<void>(transactional.applySiteDelta(duplicateBatchAdditions, {})); },
                                         "incremental EDT must reject a duplicate addition within one batch");
    requireTransactionalField(initialTransactionalField, "duplicate addition rejection preserves the incremental EDT field");
    const std::array<PixelId, 2> duplicateBatchRemovals{0, 0};
    requireThrows<std::invalid_argument>([&] { static_cast<void>(transactional.applySiteDelta({}, duplicateBatchRemovals)); },
                                         "incremental EDT must reject a duplicate removal within one batch");
    requireTransactionalField(initialTransactionalField, "duplicate removal rejection preserves the incremental EDT field");

    requireThrows<std::invalid_argument>([&] { static_cast<void>(transactional.applySiteDelta({}, transactionalInitialSite)); },
                                         "incremental EDT must reject removal of its last active site");
    requireTransactionalField(initialTransactionalField, "last-site removal rejection preserves the incremental EDT field");

    const std::array<PixelId, 1> transactionalReplacementSite{3};
    transactional.applySiteDelta(transactionalReplacementSite, transactionalInitialSite);
    requireTransactionalField({SquaredDistance{2}, SquaredDistance{1}, SquaredDistance{1}, SquaredDistance{0}},
                              "incremental EDT atomically replaces its last active site");

    constexpr int incrementalRows = 7;
    constexpr int incrementalColumns = 8;
    ExactSquaredEuclideanDistanceTransform2D incremental(incrementalRows, incrementalColumns);
    std::vector<std::uint8_t> activeSites(static_cast<std::size_t>(incrementalRows * incrementalColumns), std::uint8_t{0});
    const std::array<PixelId, 3> initialSites{0, 3 * incrementalColumns + 4, incrementalRows * incrementalColumns - 1};
    for (PixelId site : initialSites) {
        activeSites[static_cast<std::size_t>(site)] = 1;
    }
    incremental.compute(initialSites);

    const auto requireIncrementalFieldExact = [&] {
        for (PixelId pixel = 0; pixel < incrementalRows * incrementalColumns; ++pixel) {
            const SquaredDistance row = pixel / incrementalColumns;
            const SquaredDistance column = pixel % incrementalColumns;
            SquaredDistance expected = std::numeric_limits<SquaredDistance>::max();
            for (PixelId site = 0; site < incrementalRows * incrementalColumns; ++site) {
                if (activeSites[static_cast<std::size_t>(site)] == 0) {
                    continue;
                }
                const SquaredDistance deltaRow = row - site / incrementalColumns;
                const SquaredDistance deltaColumn = column - site % incrementalColumns;
                expected = std::min(expected, deltaRow * deltaRow + deltaColumn * deltaColumn);
            }
            requireEqual(incremental.squaredDistance(pixel), expected, "incremental dirty-line EDT must match brute force");
        }
    };
    requireIncrementalFieldExact();

    const std::array<PixelId, 2> firstAdditions{incrementalColumns + 1, 5 * incrementalColumns + 2};
    const std::array<PixelId, 1> firstRemovals{3 * incrementalColumns + 4};
    for (PixelId site : firstAdditions) {
        activeSites[static_cast<std::size_t>(site)] = 1;
    }
    for (PixelId site : firstRemovals) {
        activeSites[static_cast<std::size_t>(site)] = 0;
    }
    incremental.applySiteDelta(firstAdditions, firstRemovals);
    requireIncrementalFieldExact();

    const PixelId cancelledSite = 2 * incrementalColumns + 6;
    const std::array<PixelId, 1> cancelledAddition{cancelledSite};
    const std::array<PixelId, 1> cancelledRemoval{cancelledSite};
    incremental.applySiteDelta(cancelledAddition, cancelledRemoval);
    requireIncrementalFieldExact();

    std::mt19937 incrementalGenerator(0x49445432U);
    for (int step = 0; step < 160; ++step) {
        std::vector<PixelId> additions;
        std::vector<PixelId> removals;
        const PixelId candidate = static_cast<PixelId>(incrementalGenerator() % activeSites.size());
        const std::size_t activeCount = static_cast<std::size_t>(std::count(activeSites.begin(), activeSites.end(), std::uint8_t{1}));
        if (activeSites[static_cast<std::size_t>(candidate)] == 0) {
            additions.push_back(candidate);
            activeSites[static_cast<std::size_t>(candidate)] = 1;
        } else if (activeCount > 1) {
            removals.push_back(candidate);
            activeSites[static_cast<std::size_t>(candidate)] = 0;
        }
        static_cast<void>(incremental.applySiteDelta(additions, removals));
        requireIncrementalFieldExact();
    }

    PixelId activeSite = InvalidPixel;
    PixelId inactiveSite = InvalidPixel;
    for (PixelId pixel = 0; pixel < incrementalRows * incrementalColumns; ++pixel) {
        (activeSites[static_cast<std::size_t>(pixel)] != 0 ? activeSite : inactiveSite) = pixel;
    }
    if (inactiveSite == InvalidPixel) {
        const std::array<PixelId, 1> makeInactive{activeSite};
        static_cast<void>(incremental.applySiteDelta({}, makeInactive));
        activeSites[static_cast<std::size_t>(activeSite)] = 0;
        inactiveSite = activeSite;
        activeSite = InvalidPixel;
        for (PixelId pixel = 0; pixel < incrementalRows * incrementalColumns; ++pixel) {
            if (activeSites[static_cast<std::size_t>(pixel)] != 0) {
                activeSite = pixel;
                break;
            }
        }
    }
    const SquaredDistance retainedIncrementalDistance = incremental.squaredDistance(0);
    const std::array<PixelId, 1> duplicateAddition{activeSite};
    requireThrows<std::invalid_argument>([&] { static_cast<void>(incremental.applySiteDelta(duplicateAddition, {})); },
                                         "incremental EDT must reject an already-active addition");
    const std::array<PixelId, 1> inactiveRemoval{inactiveSite};
    requireThrows<std::invalid_argument>([&] { static_cast<void>(incremental.applySiteDelta({}, inactiveRemoval)); },
                                         "incremental EDT must reject an inactive removal");
    requireEqual(incremental.squaredDistance(0), retainedIncrementalDistance, "invalid incremental updates preserve the prior exact field");
}

void requireProductionMatchesExactOracle(const MorphologicalTree& tree, const std::string& label) {
    requireRegionIndexMatchesTree(tree, label);
    requireBoundaryLifetimesAndScheduleMatchOracle(tree, label);
    const std::vector<std::int64_t> expected = maxdist_oracle::exactMaxSquaredDistance(tree);
    auto [names, values] = AttributeComputation::computeSingleTopologyAttribute(tree, MaxSquaredDistExact);
    for (NodeId node : tree.aliveNodeIds()) {
        requireEqual(values[names.linearIndex(node, MaxSquaredDistExact)], static_cast<float>(expected[static_cast<std::size_t>(node)]),
                     label + " exact oracle node " + std::to_string(node));
    }

    const std::vector<AttributeOrGroup> summaryAttributes{MaxSquaredDistExact,
                                                          DistSquaredSumExact,
                                                          DistSquaredMeanExact,
                                                          DistRmsExact,
                                                          DistSquaredVarianceExact,
                                                          MaxDistCenterRowExact,
                                                          MaxDistCenterColumnExact,
                                                          MaxDistPlateauAreaExact,
                                                          MaxDistPlateauCentroidRowExact,
                                                          MaxDistPlateauCentroidColumnExact};
    auto [summaryNames, summaryValues] = AttributeComputation::computeTopologyAttributes<double>(tree, summaryAttributes);
    for (NodeId node : tree.aliveNodeIds()) {
        const maxdist_oracle::ExactNodeDistanceTransform exact = maxdist_oracle::exactNodeSquaredDistanceTransform(tree, node);
        long double sum = 0.0L;
        long double sumOfSquares = 0.0L;
        DistanceFieldExtremum expectedExtremum;
        DistanceFieldMaximumPlateau expectedPlateau;
        for (std::size_t sampleIndex = 0; sampleIndex < exact.squaredDistances.size(); ++sampleIndex) {
            const SquaredDistance distance = exact.squaredDistances[sampleIndex];
            const long double value = static_cast<long double>(distance);
            sum += value;
            sumOfSquares += value * value;
            attributes::computers::detail::distance_transform::updateDistanceFieldExtremum(expectedExtremum, exact.supportPixels[sampleIndex], distance);
            attributes::computers::detail::distance_transform::updateDistanceFieldMaximumPlateau(expectedPlateau, exact.supportPixels[sampleIndex], distance,
                                                                                                 tree.numColumns());
        }
        const long double count = static_cast<long double>(exact.squaredDistances.size());
        const long double mean = sum / count;
        const long double variance = std::max(0.0L, sumOfSquares / count - mean * mean);
        requireEqual(summaryValues[summaryNames.linearIndex(node, MaxSquaredDistExact)], static_cast<double>(expected[static_cast<std::size_t>(node)]),
                     label + " exact summary maximum node " + std::to_string(node));
        requireNear(summaryValues[summaryNames.linearIndex(node, DistSquaredSumExact)], static_cast<double>(sum), 1.0e-12,
                    label + " exact squared-distance sum node " + std::to_string(node));
        requireNear(summaryValues[summaryNames.linearIndex(node, DistSquaredMeanExact)], static_cast<double>(mean), 1.0e-12,
                    label + " exact squared-distance mean node " + std::to_string(node));
        requireNear(summaryValues[summaryNames.linearIndex(node, DistRmsExact)], static_cast<double>(std::sqrt(mean)), 1.0e-12,
                    label + " exact distance RMS node " + std::to_string(node));
        requireNear(summaryValues[summaryNames.linearIndex(node, DistSquaredVarianceExact)], static_cast<double>(variance), 1.0e-12,
                    label + " exact squared-distance variance node " + std::to_string(node));
        requireEqual(summaryValues[summaryNames.linearIndex(node, MaxDistCenterRowExact)], static_cast<double>(expectedExtremum.pixel / tree.numColumns()),
                     label + " exact maximum-distance center row node " + std::to_string(node));
        requireEqual(summaryValues[summaryNames.linearIndex(node, MaxDistCenterColumnExact)], static_cast<double>(expectedExtremum.pixel % tree.numColumns()),
                     label + " exact maximum-distance center column node " + std::to_string(node));
        requireEqual(summaryValues[summaryNames.linearIndex(node, MaxDistPlateauAreaExact)], static_cast<double>(expectedPlateau.count),
                     label + " exact maximum-distance plateau area node " + std::to_string(node));
        requireNear(summaryValues[summaryNames.linearIndex(node, MaxDistPlateauCentroidRowExact)], static_cast<double>(expectedPlateau.centroidRow()), 1.0e-12,
                    label + " exact maximum-distance plateau centroid row node " + std::to_string(node));
        requireNear(summaryValues[summaryNames.linearIndex(node, MaxDistPlateauCentroidColumnExact)], static_cast<double>(expectedPlateau.centroidColumn()),
                    1.0e-12, label + " exact maximum-distance plateau centroid column node " + std::to_string(node));
    }

    attributes::computers::detail::distance_transform::MorphologicalTreeDistanceTransform::forEachNode(
        tree, [&](const attributes::computers::detail::distance_transform::NodeDistanceTransformFrame& frame) {
            const maxdist_oracle::ExactNodeDistanceTransform exact = maxdist_oracle::exactNodeSquaredDistanceTransform(tree, frame.node());
            std::size_t sampleIndex = 0;
            SquaredDistance actualMaximum = 0;
            frame.forEachSample([&](PixelId pixel, SquaredDistance squaredDistance) {
                require(sampleIndex < exact.supportPixels.size(), label + " production EDT emitted too many samples");
                requireEqual(pixel, exact.supportPixels[sampleIndex], label + " production EDT support pixel");
                requireEqual(squaredDistance, exact.squaredDistances[sampleIndex], label + " production EDT per-pixel distance");
                actualMaximum = std::max(actualMaximum, squaredDistance);
                ++sampleIndex;
            });
            requireEqual(sampleIndex, exact.supportPixels.size(), label + " production EDT sample cardinality");
            requireEqual(actualMaximum, *std::max_element(exact.squaredDistances.begin(), exact.squaredDistances.end()), label + " streamed maximum");
        });
}

void requireSingleNodeExactPathMatchesOracle(const MorphologicalTree& tree, NodeId node, const std::string& label) {
    const maxdist_oracle::ExactNodeDistanceTransform exact = maxdist_oracle::exactNodeSquaredDistanceTransform(tree, node);
    std::vector<PixelId> actualPixels;
    std::vector<SquaredDistance> actualDistances;
    ExactProvider::forNode(tree, node, [&](const DistanceTransformFrame& frame) {
        requireEqual(frame.node(), node, label + " single-node frame id");
        frame.forEachSample([&](PixelId pixel, SquaredDistance squaredDistance) {
            actualPixels.push_back(pixel);
            actualDistances.push_back(squaredDistance);
        });
    });
    requireVectorEqual(actualPixels, exact.supportPixels, label + " single-node support pixels");
    requireVectorEqual(actualDistances, exact.squaredDistances, label + " single-node squared distances");
}

void requireReducerProjectionPrecisionPolicy(const MorphologicalTree& tree) {
    constexpr SquaredDistance beyondFloatIntegerPrecision = SquaredDistance{16777217};
    const AttributeNames names = AttributeNames::fromList({MaxSquaredDistExact});

    std::vector<float> floatBuffer(static_cast<std::size_t>(tree.numInternalNodeSlots()), 0.0f);
    ProjectedMaximumReducer<float> floatReducer(std::span<float>(floatBuffer), names);
    floatReducer.beginNode(tree.root());
    floatReducer.consumeSample(PixelId{0}, beyondFloatIntegerPrecision);
    floatReducer.endNode(tree.root());
    requireEqual(floatBuffer[names.linearIndex(tree.root(), MaxSquaredDistExact)], static_cast<float>(beyondFloatIntegerPrecision),
                 "float32 reducer follows the public floating-point projection policy");
    require(static_cast<SquaredDistance>(floatBuffer[names.linearIndex(tree.root(), MaxSquaredDistExact)]) != beyondFloatIntegerPrecision,
            "float32 projection documents loss above 2^24");

    std::vector<double> doubleBuffer(static_cast<std::size_t>(tree.numInternalNodeSlots()), 0.0);
    ProjectedMaximumReducer<double> doubleReducer(std::span<double>(doubleBuffer), names);
    doubleReducer.beginNode(tree.root());
    doubleReducer.consumeSample(PixelId{0}, beyondFloatIntegerPrecision);
    doubleReducer.endNode(tree.root());
    requireEqual(doubleBuffer[names.linearIndex(tree.root(), MaxSquaredDistExact)], static_cast<double>(beyondFloatIntegerPrecision),
                 "float64 reducer preserves the exact squared integer");
}

template <class Queue> void requireApproximateQueueSemantics(const std::string& label) {
    Queue queue(16, 6);
    queue.insert(PixelId{0}, 4);
    queue.insert(PixelId{1}, 2);
    queue.insert(PixelId{2}, 2);
    requireEqual(queue.popMinimumFifo(), PixelId{1}, label + " preserves FIFO ties");
    requireEqual(queue.popMinimumFifo(), PixelId{2}, label + " pops the second FIFO tie");
    queue.update(PixelId{0}, 1);
    requireEqual(queue.popMinimumFifo(), PixelId{0}, label + " applies decrease-key");
    queue.insert(PixelId{3}, 16);
    queue.insert(PixelId{4}, 16);
    queue.erase(PixelId{3});
    requireEqual(queue.popMinimumFifo(), PixelId{4}, label + " preserves the tail after erasing a FIFO head");
    queue.update(PixelId{5}, 3);
    require(queue.contains(PixelId{5}), label + " inserts an absent element during update");
    requireEqual(queue.popMinimumFifo(), PixelId{5}, label + " pops an element inserted by update");
    queue.insert(PixelId{0}, 7);
    queue.insert(PixelId{1}, 7);
    queue.erase(PixelId{1});
    requireEqual(queue.popMinimumFifo(), PixelId{0}, label + " preserves the head after erasing a FIFO tail");
    require(queue.empty(), label + " can be reused across empty-queue transitions");
}

void requireQueueConstructionContracts() {
    requireThrows<std::invalid_argument>([] { static_cast<void>(ApproximateQueue(-1, 2)); },
                                         "production queue constructor always rejects a negative bucket domain");
    requireThrows<std::invalid_argument>([] { static_cast<void>(ApproximateQueue(4, 0)); },
                                         "production queue constructor always rejects an empty element domain");
}

void requireAdaptiveApproximateProviderBoundary() {
    requireApproximateQueueSemantics<ApproximateQueue>("production adaptive DIFT queue");
    requireQueueConstructionContracts();

    ApproximateDift transform(5, 5);
    std::vector<PixelId> contour;
    for (PixelId pixel = 0; pixel < 25; ++pixel) {
        transform.addPixelToBinaryImage(pixel);
        const int row = pixel / 5;
        const int column = pixel % 5;
        if (row == 0 || row == 4 || column == 0 || column == 4) {
            contour.push_back(pixel);
            transform.seed(pixel);
        } else {
            transform.open(pixel);
        }
    }
    transform.run();
    requireEqual(transform.squaredDistance(PixelId{12}), std::int64_t{4}, "adaptive DIFT full-square center distance");
    requireEqual(transform.maximumRootDistance(contour), std::int64_t{4}, "adaptive DIFT full-square maximum");

    auto constantImage = ImageUInt8::create(5, 5);
    constantImage->fill(7);
    auto constantTree = makeComponentTree(constantImage, true);
    const auto denseResult = ApproximateTreeTransform::computeMaxSquaredDistance(*constantTree);
    std::vector<std::int64_t> consumedValues(static_cast<std::size_t>(constantTree->numInternalNodeSlots()), -1);
    ApproximateTreeTransform::forEachNodeMaximum(
        *constantTree, [&consumedValues](NodeId node, std::int64_t value) { consumedValues[static_cast<std::size_t>(node)] = value; });
    auto [names, publicValues] = AttributeComputation::computeSingleTopologyAttribute<double>(*constantTree, MaxSquaredDist);

    for (NodeId node : constantTree->aliveNodeIds()) {
        const std::int64_t expected = denseResult[static_cast<std::size_t>(node)];
        requireEqual(consumedValues[static_cast<std::size_t>(node)], expected, "adaptive provider consumer receives the standalone local decision");
        requireEqual(publicValues[names.linearIndex(node, MaxSquaredDist)], static_cast<double>(expected),
                     "public approximate attribute materializes the provider decision");
    }
    requireEqual(consumedValues[static_cast<std::size_t>(constantTree->root())], std::int64_t{4}, "adaptive approximate provider constant-tree maximum");
    std::vector<DistanceFieldExtremum> extrema(static_cast<std::size_t>(constantTree->numInternalNodeSlots()));
    ApproximateTreeTransform::forEachNodeExtremum(
        *constantTree, [&extrema](NodeId node, const DistanceFieldExtremum& extremum) { extrema[static_cast<std::size_t>(node)] = extremum; });
    requireEqual(extrema[static_cast<std::size_t>(constantTree->root())].squaredDistance, std::int64_t{4},
                 "localized approximate provider constant-tree maximum");
    requireEqual(extrema[static_cast<std::size_t>(constantTree->root())].pixel, PixelId{12}, "localized approximate provider constant-tree center");

    std::vector<DistanceFieldMaximumPlateau> plateaus(static_cast<std::size_t>(constantTree->numInternalNodeSlots()));
    ApproximateTreeTransform::forEachNodePlateau(
        *constantTree, [&plateaus](NodeId node, const DistanceFieldMaximumPlateau& plateau) { plateaus[static_cast<std::size_t>(node)] = plateau; });
    const DistanceFieldMaximumPlateau& constantPlateau = plateaus[static_cast<std::size_t>(constantTree->root())];
    requireEqual(constantPlateau.squaredDistance, std::int64_t{4}, "constant approximate plateau maximum");
    requireEqual(constantPlateau.pixel, PixelId{12}, "constant approximate plateau canonical center");
    requireEqual(constantPlateau.count, std::uint64_t{1}, "constant approximate plateau area");
    requireNear(static_cast<double>(constantPlateau.centroidRow()), 2.0, 1.0e-12, "constant approximate plateau centroid row");
    requireNear(static_cast<double>(constantPlateau.centroidColumn()), 2.0, 1.0e-12, "constant approximate plateau centroid column");

    std::vector<long double> summarySums(static_cast<std::size_t>(constantTree->numInternalNodeSlots()), -1.0L);
    ApproximateTreeTransform::forEachNodeSummaryAndExtremum(
        *constantTree, [&summarySums](NodeId node, const DistanceFieldExtremum& extremum,
                                      const attributes::computers::detail::distance_transform::DistanceFieldMoments& moments) {
            requireEqual(extremum.squaredDistance, std::int64_t{4}, "constant approximate summary Bedt maximum");
            requireEqual(extremum.pixel, PixelId{12}, "constant approximate summary canonical center");
            requireEqual(moments.count(), std::uint64_t{25}, "constant approximate summary support count");
            requireNear(static_cast<double>(moments.sum()), 12.0, 1.0e-12, "constant approximate summary squared-distance sum");
            requireNear(static_cast<double>(moments.mean()), 12.0 / 25.0, 1.0e-12, "constant approximate summary squared-distance mean");
            requireNear(static_cast<double>(moments.rms()), std::sqrt(12.0 / 25.0), 1.0e-12, "constant approximate summary distance RMS");
            requireNear(static_cast<double>(moments.populationVariance()), 0.7296, 1.0e-12, "constant approximate summary squared-distance variance");
            summarySums[static_cast<std::size_t>(node)] = moments.sum();
        });
    requireNear(static_cast<double>(summarySums[static_cast<std::size_t>(constantTree->root())]), 12.0, 1.0e-12,
                "constant approximate provider materializes the root summary");

    const std::vector<AttributeOrGroup> approximateAttributes{MaxSquaredDist,
                                                              DistSquaredSum,
                                                              DistSquaredMean,
                                                              DistRms,
                                                              DistSquaredVariance,
                                                              MaxDistCenterRow,
                                                              MaxDistCenterColumn,
                                                              MaxDistPlateauArea,
                                                              MaxDistPlateauCentroidRow,
                                                              MaxDistPlateauCentroidColumn};
    auto [summaryNames, summaryValues] = AttributeComputation::computeTopologyAttributes<double>(*constantTree, approximateAttributes);
    const NodeId root = constantTree->root();
    requireEqual(summaryValues[summaryNames.linearIndex(root, MaxSquaredDist)], 4.0, "public approximate summary maximum");
    requireEqual(summaryValues[summaryNames.linearIndex(root, DistSquaredSum)], 12.0, "public approximate squared-distance sum");
    requireNear(summaryValues[summaryNames.linearIndex(root, DistSquaredMean)], 12.0 / 25.0, 1.0e-12, "public approximate squared-distance mean");
    requireNear(summaryValues[summaryNames.linearIndex(root, DistRms)], std::sqrt(12.0 / 25.0), 1.0e-12, "public approximate distance RMS");
    requireNear(summaryValues[summaryNames.linearIndex(root, DistSquaredVariance)], 0.7296, 1.0e-12, "public approximate squared-distance variance");
    requireEqual(summaryValues[summaryNames.linearIndex(root, MaxDistCenterRow)], 2.0, "public approximate maximum-distance center row");
    requireEqual(summaryValues[summaryNames.linearIndex(root, MaxDistCenterColumn)], 2.0, "public approximate maximum-distance center column");
    requireEqual(summaryValues[summaryNames.linearIndex(root, MaxDistPlateauArea)], 1.0, "public approximate maximum-distance plateau area");
    requireEqual(summaryValues[summaryNames.linearIndex(root, MaxDistPlateauCentroidRow)], 2.0, "public approximate maximum-distance plateau centroid row");
    requireEqual(summaryValues[summaryNames.linearIndex(root, MaxDistPlateauCentroidColumn)], 2.0,
                 "public approximate maximum-distance plateau centroid column");

    const std::vector<AttributeOrGroup> centerOnlyAttributes{MaxDistCenterRow, MaxDistCenterColumn};
    auto [centerNames, centerValues] = AttributeComputation::computeTopologyAttributes<double>(*constantTree, centerOnlyAttributes);
    requireEqual(centerValues[centerNames.linearIndex(root, MaxDistCenterRow)], 2.0, "public approximate center-only row");
    requireEqual(centerValues[centerNames.linearIndex(root, MaxDistCenterColumn)], 2.0, "public approximate center-only column");

    const std::vector<AttributeOrGroup> plateauOnlyAttributes{MaxDistPlateauArea, MaxDistPlateauCentroidRow, MaxDistPlateauCentroidColumn};
    auto [plateauNames, plateauValues] = AttributeComputation::computeTopologyAttributes<double>(*constantTree, plateauOnlyAttributes);
    requireEqual(plateauValues[plateauNames.linearIndex(root, MaxDistPlateauArea)], 1.0, "public approximate plateau-only area");
    requireEqual(plateauValues[plateauNames.linearIndex(root, MaxDistPlateauCentroidRow)], 2.0, "public approximate plateau-only centroid row");
    requireEqual(plateauValues[plateauNames.linearIndex(root, MaxDistPlateauCentroidColumn)], 2.0, "public approximate plateau-only centroid column");
}

void requireApproximateSummariesMatchExactSmallGrid(const MorphologicalTree& tree, const std::string& label) {
    const std::size_t numSlots = static_cast<std::size_t>(tree.numInternalNodeSlots());
    std::vector<std::int64_t> maxima(numSlots, -1);
    std::vector<long double> sums(numSlots, -1.0L);
    std::vector<long double> means(numSlots, -1.0L);
    std::vector<long double> rmsValues(numSlots, -1.0L);
    std::vector<long double> variances(numSlots, -1.0L);
    std::vector<PixelId> centers(numSlots, InvalidPixel);
    std::vector<DistanceFieldMaximumPlateau> plateaus(numSlots);
    ApproximateTreeTransform::forEachNodeSummaryAndPlateau(tree, [&](NodeId node, const DistanceFieldMaximumPlateau& plateau,
                                                                     const attributes::computers::detail::distance_transform::DistanceFieldMoments& moments) {
        const std::size_t index = static_cast<std::size_t>(node);
        maxima[index] = plateau.squaredDistance;
        centers[index] = plateau.pixel;
        plateaus[index] = plateau;
        sums[index] = moments.sum();
        means[index] = moments.mean();
        rmsValues[index] = moments.rms();
        variances[index] = moments.populationVariance();
    });
    std::vector<std::int64_t> auditedMaxima(numSlots, -1);
    AuditedFastApproximateTreeTransform::forEachNodeMaximum(
        tree, [&auditedMaxima](NodeId node, std::int64_t maximum) { auditedMaxima[static_cast<std::size_t>(node)] = maximum; });
    for (NodeId node : tree.aliveNodeIds()) {
        requireEqual(auditedMaxima[static_cast<std::size_t>(node)], maxima[static_cast<std::size_t>(node)],
                     label + " audited fast-path maximum node " + std::to_string(node));
    }

    const std::vector<AttributeOrGroup> requested{MaxSquaredDist,
                                                  DistSquaredSum,
                                                  DistSquaredMean,
                                                  DistRms,
                                                  DistSquaredVariance,
                                                  MaxDistCenterRow,
                                                  MaxDistCenterColumn,
                                                  MaxDistPlateauArea,
                                                  MaxDistPlateauCentroidRow,
                                                  MaxDistPlateauCentroidColumn};
    auto [names, values] = AttributeComputation::computeTopologyAttributes<double>(tree, requested);
    for (NodeId node : tree.aliveNodeIds()) {
        const std::size_t index = static_cast<std::size_t>(node);
        requireEqual(values[names.linearIndex(node, MaxSquaredDist)], static_cast<double>(maxima[index]), label + " approximate public/provider maximum");
        requireNear(values[names.linearIndex(node, DistSquaredSum)], static_cast<double>(sums[index]), 1.0e-12, label + " approximate public/provider sum");
        requireNear(values[names.linearIndex(node, DistSquaredMean)], static_cast<double>(means[index]), 1.0e-12, label + " approximate public/provider mean");
        requireNear(values[names.linearIndex(node, DistRms)], static_cast<double>(rmsValues[index]), 1.0e-12, label + " approximate public/provider RMS");
        requireNear(values[names.linearIndex(node, DistSquaredVariance)], static_cast<double>(variances[index]), 1.0e-12,
                    label + " approximate public/provider variance");
        requireEqual(values[names.linearIndex(node, MaxDistCenterRow)], static_cast<double>(centers[index] / tree.numColumns()),
                     label + " approximate public/provider center row");
        requireEqual(values[names.linearIndex(node, MaxDistCenterColumn)], static_cast<double>(centers[index] % tree.numColumns()),
                     label + " approximate public/provider center column");
        requireEqual(values[names.linearIndex(node, MaxDistPlateauArea)], static_cast<double>(plateaus[index].count),
                     label + " approximate public/provider plateau area");
        requireNear(values[names.linearIndex(node, MaxDistPlateauCentroidRow)], static_cast<double>(plateaus[index].centroidRow()), 1.0e-12,
                    label + " approximate public/provider plateau centroid row");
        requireNear(values[names.linearIndex(node, MaxDistPlateauCentroidColumn)], static_cast<double>(plateaus[index].centroidColumn()), 1.0e-12,
                    label + " approximate public/provider plateau centroid column");

        const maxdist_oracle::ExactNodeDistanceTransform exact = maxdist_oracle::exactNodeSquaredDistanceTransform(tree, node);
        long double exactSum = 0.0L;
        long double exactSumOfSquares = 0.0L;
        DistanceFieldExtremum exactExtremum;
        DistanceFieldMaximumPlateau exactPlateau;
        for (std::size_t sampleIndex = 0; sampleIndex < exact.squaredDistances.size(); ++sampleIndex) {
            const SquaredDistance distance = exact.squaredDistances[sampleIndex];
            const long double value = static_cast<long double>(distance);
            exactSum += value;
            exactSumOfSquares += value * value;
            attributes::computers::detail::distance_transform::updateDistanceFieldExtremum(exactExtremum, exact.supportPixels[sampleIndex], distance);
            attributes::computers::detail::distance_transform::updateDistanceFieldMaximumPlateau(exactPlateau, exact.supportPixels[sampleIndex], distance,
                                                                                                 tree.numColumns());
        }
        const long double count = static_cast<long double>(exact.squaredDistances.size());
        const long double exactMean = exactSum / count;
        const long double exactVariance = std::max(0.0L, exactSumOfSquares / count - exactMean * exactMean);
        requireNear(static_cast<double>(sums[index]), static_cast<double>(exactSum), 1.0e-12, label + " approximate/exact small-grid sum");
        requireNear(static_cast<double>(means[index]), static_cast<double>(exactMean), 1.0e-12, label + " approximate/exact small-grid mean");
        requireNear(static_cast<double>(rmsValues[index]), static_cast<double>(std::sqrt(exactMean)), 1.0e-12, label + " approximate/exact small-grid RMS");
        requireNear(static_cast<double>(variances[index]), static_cast<double>(exactVariance), 1.0e-12, label + " approximate/exact small-grid variance");
        requireEqual(centers[index], exactExtremum.pixel, label + " approximate/exact small-grid canonical center");
        requireEqual(plateaus[index].count, exactPlateau.count, label + " approximate/exact small-grid plateau area");
        requireNear(static_cast<double>(plateaus[index].centroidRow()), static_cast<double>(exactPlateau.centroidRow()), 1.0e-12,
                    label + " approximate/exact small-grid plateau centroid row");
        requireNear(static_cast<double>(plateaus[index].centroidColumn()), static_cast<double>(exactPlateau.centroidColumn()), 1.0e-12,
                    label + " approximate/exact small-grid plateau centroid column");
    }
}

void requireCompositeReducersMatchExactOracle(const MorphologicalTree& tree, const std::string& label) {
    DenseMaximumReducer maximumReducer(tree.numInternalNodeSlots());
    DenseDistanceSumReducer sumReducer(tree.numInternalNodeSlots());
    attributes::computers::detail::distance_transform::MorphologicalTreeDistanceTransform::reduce(tree, maximumReducer, sumReducer);

    requireEqual(maximumReducer.beginCount, static_cast<std::size_t>(tree.numNodes()), label + " maximum reducer begin count");
    requireEqual(maximumReducer.endCount, static_cast<std::size_t>(tree.numNodes()), label + " maximum reducer end count");
    requireEqual(sumReducer.beginCount, static_cast<std::size_t>(tree.numNodes()), label + " sum reducer begin count");
    requireEqual(sumReducer.endCount, static_cast<std::size_t>(tree.numNodes()), label + " sum reducer end count");
    requireEqual(maximumReducer.sampleCount, sumReducer.sampleCount, label + " reducers receive the same shared samples");

    for (NodeId node : tree.aliveNodeIds()) {
        const maxdist_oracle::ExactNodeDistanceTransform exact = maxdist_oracle::exactNodeSquaredDistanceTransform(tree, node);
        const SquaredDistance expectedMaximum = *std::max_element(exact.squaredDistances.begin(), exact.squaredDistances.end());
        SquaredDistance expectedSum = 0;
        for (SquaredDistance squaredDistance : exact.squaredDistances) {
            expectedSum += squaredDistance;
        }
        requireEqual(maximumReducer.values[static_cast<std::size_t>(node)], expectedMaximum, label + " maximum reducer value");
        requireEqual(sumReducer.values[static_cast<std::size_t>(node)], expectedSum, label + " sum reducer value");
    }
    for (NodeId node = 0; node < tree.numInternalNodeSlots(); ++node) {
        if (!tree.isAlive(node)) {
            requireEqual(maximumReducer.values[static_cast<std::size_t>(node)], SquaredDistance{-1}, label + " dead-slot maximum remains untouched");
            requireEqual(sumReducer.values[static_cast<std::size_t>(node)], SquaredDistance{-1}, label + " dead-slot sum remains untouched");
        }
    }
}

void requireUniformBackendAndOracles() {
    auto constantImage = ImageUInt8::create(5, 5);
    constantImage->fill(7);
    auto constantTree = makeComponentTree(constantImage, true);
    requireReducerProjectionPrecisionPolicy(*constantTree);
    requireProductionMatchesExactOracle(*constantTree, "constant 5x5 max-tree");
    auto [constantNames, constantValues] = AttributeComputation::computeSingleTopologyAttribute(*constantTree, MaxSquaredDistExact);
    requireEqual(constantValues[constantNames.linearIndex(constantTree->root(), MaxSquaredDistExact)], 4.0f, "constant 5x5 squared MAX_SQUARED_DIST_EXACT");

    const std::vector<AttributeOrGroup> exactAttributes{MaxSquaredDistExact,
                                                        DistSquaredSumExact,
                                                        DistSquaredMeanExact,
                                                        DistRmsExact,
                                                        DistSquaredVarianceExact,
                                                        MaxDistCenterRowExact,
                                                        MaxDistCenterColumnExact,
                                                        MaxDistPlateauAreaExact,
                                                        MaxDistPlateauCentroidRowExact,
                                                        MaxDistPlateauCentroidColumnExact};
    auto [exactNames, exactValues] = AttributeComputation::computeTopologyAttributes<double>(*constantTree, exactAttributes);
    const NodeId constantRoot = constantTree->root();
    requireEqual(exactValues[exactNames.linearIndex(constantRoot, MaxSquaredDistExact)], 4.0, "constant public exact maximum");
    requireEqual(exactValues[exactNames.linearIndex(constantRoot, DistSquaredSumExact)], 12.0, "constant public exact squared-distance sum");
    requireNear(exactValues[exactNames.linearIndex(constantRoot, DistSquaredMeanExact)], 12.0 / 25.0, 1.0e-12, "constant public exact squared-distance mean");
    requireNear(exactValues[exactNames.linearIndex(constantRoot, DistRmsExact)], std::sqrt(12.0 / 25.0), 1.0e-12, "constant public exact distance RMS");
    requireNear(exactValues[exactNames.linearIndex(constantRoot, DistSquaredVarianceExact)], 0.7296, 1.0e-12,
                "constant public exact squared-distance variance");
    requireEqual(exactValues[exactNames.linearIndex(constantRoot, MaxDistCenterRowExact)], 2.0, "constant public exact maximum-distance center row");
    requireEqual(exactValues[exactNames.linearIndex(constantRoot, MaxDistCenterColumnExact)], 2.0, "constant public exact maximum-distance center column");
    requireEqual(exactValues[exactNames.linearIndex(constantRoot, MaxDistPlateauAreaExact)], 1.0, "constant public exact maximum-distance plateau area");
    requireEqual(exactValues[exactNames.linearIndex(constantRoot, MaxDistPlateauCentroidRowExact)], 2.0,
                 "constant public exact maximum-distance plateau centroid row");
    requireEqual(exactValues[exactNames.linearIndex(constantRoot, MaxDistPlateauCentroidColumnExact)], 2.0,
                 "constant public exact maximum-distance plateau centroid column");

    const std::vector<AttributeOrGroup> exactCenterOnlyAttributes{MaxDistCenterRowExact, MaxDistCenterColumnExact};
    auto [exactCenterNames, exactCenterValues] = AttributeComputation::computeTopologyAttributes<double>(*constantTree, exactCenterOnlyAttributes);
    requireEqual(exactCenterValues[exactCenterNames.linearIndex(constantRoot, MaxDistCenterRowExact)], 2.0, "constant public exact center-only row");
    requireEqual(exactCenterValues[exactCenterNames.linearIndex(constantRoot, MaxDistCenterColumnExact)], 2.0, "constant public exact center-only column");

    const std::vector<AttributeOrGroup> exactPlateauOnlyAttributes{MaxDistPlateauAreaExact, MaxDistPlateauCentroidRowExact, MaxDistPlateauCentroidColumnExact};
    auto [exactPlateauNames, exactPlateauValues] = AttributeComputation::computeTopologyAttributes<double>(*constantTree, exactPlateauOnlyAttributes);
    requireEqual(exactPlateauValues[exactPlateauNames.linearIndex(constantRoot, MaxDistPlateauAreaExact)], 1.0, "constant public exact plateau-only area");

    auto tiedImage = ImageUInt8::create(2, 2);
    tiedImage->fill(9);
    auto tiedTree = makeComponentTree(tiedImage, true);
    const std::vector<AttributeOrGroup> tiedCenters{
        MaxDistCenterRowExact,          MaxDistCenterColumnExact,          MaxDistCenterRow,   MaxDistCenterColumn,       MaxDistPlateauAreaExact,
        MaxDistPlateauCentroidRowExact, MaxDistPlateauCentroidColumnExact, MaxDistPlateauArea, MaxDistPlateauCentroidRow, MaxDistPlateauCentroidColumn};
    auto [tiedNames, tiedValues] = AttributeComputation::computeTopologyAttributes<double>(*tiedTree, tiedCenters);
    const NodeId tiedRoot = tiedTree->root();
    requireEqual(tiedValues[tiedNames.linearIndex(tiedRoot, MaxDistCenterRowExact)], 0.0, "exact tied maximum selects the smallest row-major row");
    requireEqual(tiedValues[tiedNames.linearIndex(tiedRoot, MaxDistCenterColumnExact)], 0.0, "exact tied maximum selects the smallest row-major column");
    requireEqual(tiedValues[tiedNames.linearIndex(tiedRoot, MaxDistCenterRow)], 0.0, "approximate tied maximum selects the smallest row-major row");
    requireEqual(tiedValues[tiedNames.linearIndex(tiedRoot, MaxDistCenterColumn)], 0.0, "approximate tied maximum selects the smallest row-major column");
    requireEqual(tiedValues[tiedNames.linearIndex(tiedRoot, MaxDistPlateauAreaExact)], 4.0, "exact tied maximum plateau contains every pixel");
    requireEqual(tiedValues[tiedNames.linearIndex(tiedRoot, MaxDistPlateauArea)], 4.0, "approximate tied maximum plateau contains every pixel");
    requireEqual(tiedValues[tiedNames.linearIndex(tiedRoot, MaxDistPlateauCentroidRowExact)], 0.5, "exact tied maximum plateau centroid row");
    requireEqual(tiedValues[tiedNames.linearIndex(tiedRoot, MaxDistPlateauCentroidColumnExact)], 0.5, "exact tied maximum plateau centroid column");
    requireEqual(tiedValues[tiedNames.linearIndex(tiedRoot, MaxDistPlateauCentroidRow)], 0.5, "approximate tied maximum plateau centroid row");
    requireEqual(tiedValues[tiedNames.linearIndex(tiedRoot, MaxDistPlateauCentroidColumn)], 0.5, "approximate tied maximum plateau centroid column");

    DenseMaximumReducer maximumReducer(constantTree->numInternalNodeSlots());
    DenseDistanceSumReducer sumReducer(constantTree->numInternalNodeSlots());
    attributes::computers::detail::distance_transform::MorphologicalTreeDistanceTransform::reduce(*constantTree, maximumReducer, sumReducer);
    requireEqual(maximumReducer.beginCount, std::size_t{1}, "maximum reducer begins once");
    requireEqual(maximumReducer.endCount, std::size_t{1}, "maximum reducer ends once");
    requireEqual(sumReducer.beginCount, std::size_t{1}, "sum reducer begins once");
    requireEqual(sumReducer.endCount, std::size_t{1}, "sum reducer ends once");
    requireEqual(maximumReducer.sampleCount, std::size_t{25}, "maximum reducer sample count");
    requireEqual(sumReducer.sampleCount, std::size_t{25}, "sum reducer receives the same shared samples");
    requireEqual(maximumReducer.values[static_cast<std::size_t>(constantTree->root())], SquaredDistance{4}, "shared maximum reducer value");
    requireEqual(sumReducer.values[static_cast<std::size_t>(constantTree->root())], SquaredDistance{12}, "shared sum reducer value");
    auto independentConstantTree = makeComponentTree(constantImage, true);
    RegionIndex foreignRegions(*independentConstantTree);
    requireThrows<std::invalid_argument>(
        [&] {
            static_cast<void>(ContourScheduler::forEachNode(*constantTree, foreignRegions,
                                                            [](const attributes::computers::detail::distance_transform::NodeContourFrame&) {}));
        },
        "contour scheduler must reject a region index from another tree instance");

    auto fixture = makeComponentTreeFixture();
    for (bool isMaxtree : {true, false}) {
        auto valuedTree = makeValuedComponentTree(fixture, isMaxtree);
        const MorphologicalTree& tree = valuedTree->topology();
        requireProductionMatchesExactOracle(tree, isMaxtree ? "max-tree" : "min-tree");
        requireSingleNodeExactPathMatchesOracle(tree, tree.root(), isMaxtree ? "max-tree" : "min-tree");
        NodeId lastLiveNode = tree.root();
        for (NodeId liveNode : tree.aliveNodeIds()) {
            lastLiveNode = liveNode;
        }
        if (lastLiveNode != tree.root()) {
            requireSingleNodeExactPathMatchesOracle(tree, lastLiveNode, isMaxtree ? "max-tree leaf" : "min-tree leaf");
        }
        requireApproximateSummariesMatchExactSmallGrid(tree, isMaxtree ? "max-tree" : "min-tree");

        const auto& altitude = valuedTree->nodeAltitudes();
        const std::vector<int> optimizedOracle =
            maxdist_oracle::componentTreeDiftMaxSquaredDistance(tree, std::span<const std::uint8_t>(altitude.data(), altitude.size()));
        auto [names, values] = AttributeComputation::computeSingleTopologyAttribute(tree, MaxSquaredDistExact);
        for (NodeId node : tree.aliveNodeIds()) {
            requireEqual(values[names.linearIndex(node, MaxSquaredDistExact)], static_cast<float>(optimizedOracle[static_cast<std::size_t>(node)]),
                         isMaxtree ? "max-tree optimized oracle" : "min-tree optimized oracle");
        }
    }

    auto editedTree = makeComponentTree(fixture, true);
    [[maybe_unused]] RegionIndex staleRegions(*editedTree);
    [[maybe_unused]] BoundaryLifetimeIndex staleBoundaryLifetimes(*editedTree);
    editedTree->mergeNodeIntoParent(5);
    if constexpr (contract::validationsEnabled) {
        requireThrows<std::logic_error>([&] { static_cast<void>(staleRegions.support(editedTree->root())); },
                                        "region index must reject access after a committed tree mutation");
        requireThrows<std::logic_error>([&] { static_cast<void>(staleBoundaryLifetimes.lifetime(0)); },
                                        "boundary-lifetime index must reject access after a committed tree mutation");
    }
    requireThrows<std::logic_error>(
        [&] {
            static_cast<void>(
                ContourScheduler::forEachNode(*editedTree, staleRegions, [](const attributes::computers::detail::distance_transform::NodeContourFrame&) {}));
        },
        "contour scheduler must reject a stale region index");
    RegionIndex editedRegions(*editedTree);
    requireThrows<std::out_of_range>([&] { static_cast<void>(editedRegions.support(5)); }, "region index must reject a dead dense node slot");
    requireProductionMatchesExactOracle(*editedTree, "max-tree with a dead node slot");
    requireSingleNodeExactPathMatchesOracle(*editedTree, editedTree->root(), "max-tree with a dead node slot");
    requireCompositeReducersMatchExactOracle(*editedTree, "composite reducers on a tree with a dead node slot");

    for (TestTopographicImmersion immersion :
         {TestTopographicImmersion::SelfDualSpan, TestTopographicImmersion::Min4Max8, TestTopographicImmersion::Min8Max4}) {
        auto treeOfShapes = makeTreeOfShapes(fixture, immersion);
        requireProductionMatchesExactOracle(*treeOfShapes, "tree of shapes projected-support convention");
        requireSingleNodeExactPathMatchesOracle(*treeOfShapes, treeOfShapes->root(), "tree of shapes projected-support convention");
        requireApproximateSummariesMatchExactSmallGrid(*treeOfShapes, "tree of shapes projected-support convention");
    }

    auto virtualRootImage = makeImage(2, 2, {1, 1, 0, 0});
    auto virtualRootTreeOfShapes = makeTreeOfShapes(virtualRootImage, TestTopographicImmersion::SelfDualSpan);
    requireEqual(virtualRootTreeOfShapes->properPartCardinality(virtualRootTreeOfShapes->root()), 0, "virtual-root ToS fixture has an empty root proper part");
    requireProductionMatchesExactOracle(*virtualRootTreeOfShapes, "tree of shapes with a virtual empty-proper-part root");
    requireSingleNodeExactPathMatchesOracle(*virtualRootTreeOfShapes, virtualRootTreeOfShapes->root(), "tree of shapes with a virtual empty-proper-part root");

    const std::vector<NodeId> genericParent{0, 0, 0};
    const std::vector<NodeId> genericSmallestNode{1, 1, 2, 2};
    const std::vector<std::uint8_t> genericAltitude{7, 2, 9};
    auto genericTree = MorphologicalTreeFactory::createFromNativeTopology(
        std::span<const NodeId>(genericParent), std::span<const NodeId>(genericSmallestNode), std::span<const std::uint8_t>(genericAltitude), NodeId{0}, 1, 4,
        MorphologicalTreeSemantics{MorphologicalTreeKind::Generic, NodeAltitudeOrder::Unconstrained, NoConstructionContext{}});
    requireProductionMatchesExactOracle(genericTree.topology(), "unconstrained generic tree without construction context");
    requireSingleNodeExactPathMatchesOracle(genericTree.topology(), genericTree.topology().root(), "unconstrained generic tree without construction context");
    requireApproximateSummariesMatchExactSmallGrid(genericTree.topology(), "unconstrained generic tree without construction context");
    auto transposedGenericTree = MorphologicalTreeFactory::createFromNativeTopology(
        std::span<const NodeId>(genericParent), std::span<const NodeId>(genericSmallestNode), std::span<const std::uint8_t>(genericAltitude), NodeId{0}, 4, 1,
        MorphologicalTreeSemantics{MorphologicalTreeKind::Generic, NodeAltitudeOrder::Unconstrained, NoConstructionContext{}});
    requireProductionMatchesExactOracle(transposedGenericTree.topology(), "single-column unconstrained generic tree");
    requireSingleNodeExactPathMatchesOracle(transposedGenericTree.topology(), transposedGenericTree.topology().root(),
                                            "single-column unconstrained generic tree");
    requireApproximateSummariesMatchExactSmallGrid(transposedGenericTree.topology(), "single-column unconstrained generic tree");
    const std::vector<NodeId> nonzeroRootParent{2, 2, 2};
    const std::vector<NodeId> nonzeroRootSmallestNode{0, 0, 1, 1};
    const std::vector<std::uint8_t> nonzeroRootAltitude{4, 8, 0};
    auto nonzeroRootTree = MorphologicalTreeFactory::createFromNativeTopology(
        std::span<const NodeId>(nonzeroRootParent), std::span<const NodeId>(nonzeroRootSmallestNode), std::span<const std::uint8_t>(nonzeroRootAltitude),
        NodeId{2}, 1, 4, MorphologicalTreeSemantics{MorphologicalTreeKind::Generic, NodeAltitudeOrder::Unconstrained, NoConstructionContext{}});
    requireProductionMatchesExactOracle(nonzeroRootTree.topology(), "generic tree with a nonzero root slot");
    requireSingleNodeExactPathMatchesOracle(nonzeroRootTree.topology(), nonzeroRootTree.topology().root(), "generic tree with a nonzero root slot");

    constexpr int ringRows = 9;
    constexpr int ringColumns = 10;
    const std::vector<NodeId> ringParent{0, 0};
    std::vector<NodeId> ringSmallestNode(static_cast<std::size_t>(ringRows * ringColumns), NodeId{0});
    for (int row = 1; row <= 7; ++row) {
        for (int column = 2; column <= 8; ++column) {
            const bool inHole = row >= 3 && row <= 5 && column >= 4 && column <= 6;
            if (!inHole) {
                ringSmallestNode[static_cast<std::size_t>(row * ringColumns + column)] = 1;
            }
        }
    }
    const std::vector<std::uint8_t> ringAltitude{0, 1};
    auto ringTree = MorphologicalTreeFactory::createFromNativeTopology(
        std::span<const NodeId>(ringParent), std::span<const NodeId>(ringSmallestNode), std::span<const std::uint8_t>(ringAltitude), NodeId{0}, ringRows,
        ringColumns, MorphologicalTreeSemantics{MorphologicalTreeKind::Generic, NodeAltitudeOrder::Unconstrained, NoConstructionContext{}});
    requireProductionMatchesExactOracle(ringTree.topology(), "shifted generic ring with a hole");
    requireSingleNodeExactPathMatchesOracle(ringTree.topology(), NodeId{1}, "shifted generic ring with a hole");
    requireCompositeReducersMatchExactOracle(ringTree.topology(), "composite reducers on a shifted ring");
    requireApproximateSummariesMatchExactSmallGrid(ringTree.topology(), "shifted generic ring with a hole");

    const std::vector<NodeId> emptyChainParent{0, 0, 1};
    const std::vector<NodeId> emptyChainSmallestNode{2};
    const std::vector<std::uint8_t> emptyChainAltitude{0, 0, 0};
    auto emptyProperPartChain = MorphologicalTreeFactory::createFromNativeTopology(
        std::span<const NodeId>(emptyChainParent), std::span<const NodeId>(emptyChainSmallestNode), std::span<const std::uint8_t>(emptyChainAltitude),
        NodeId{0}, 1, 1, MorphologicalTreeSemantics{MorphologicalTreeKind::Generic, NodeAltitudeOrder::Unconstrained, NoConstructionContext{}});
    requireProductionMatchesExactOracle(emptyProperPartChain.topology(), "generic unary chain with empty proper parts");
    requireSingleNodeExactPathMatchesOracle(emptyProperPartChain.topology(), NodeId{1}, "generic unary chain with empty proper parts");
    requireApproximateSummariesMatchExactSmallGrid(emptyProperPartChain.topology(), "generic unary chain with empty proper parts");
    auto unrestrictedResidual = MorphologicalTreeFactory::createUnrestrictedResidualTree(fixture);
    requireProductionMatchesExactOracle(unrestrictedResidual.topology(), "unrestricted residual tree");
    requireSingleNodeExactPathMatchesOracle(unrestrictedResidual.topology(), unrestrictedResidual.topology().root(), "unrestricted residual tree");
    requireApproximateSummariesMatchExactSmallGrid(unrestrictedResidual.topology(), "unrestricted residual tree");

    auto saturatedResidual = MorphologicalTreeFactory::createSaturatedResidualTree(fixture, PixelId{0});
    requireProductionMatchesExactOracle(saturatedResidual.topology(), "saturated residual tree");
    requireSingleNodeExactPathMatchesOracle(saturatedResidual.topology(), saturatedResidual.topology().root(), "saturated residual tree");
    requireApproximateSummariesMatchExactSmallGrid(saturatedResidual.topology(), "saturated residual tree");
}

void requireRandomizedSmallGridAgreement() {
    std::mt19937 generator(0x4D415844U);
    std::uniform_int_distribution<int> grayLevel(0, 7);

    for (int trial = 0; trial < 24; ++trial) {
        const int rows = 3 + trial % 3;
        const int columns = 3 + (trial / 3) % 3;
        auto image = ImageUInt8::create(rows, columns);
        for (PixelId pixel = 0; pixel < rows * columns; ++pixel) {
            (*image)[pixel] = static_cast<std::uint8_t>(grayLevel(generator));
        }

        for (bool isMaxtree : {true, false}) {
            auto valuedTree = makeValuedComponentTree(image, isMaxtree);
            const MorphologicalTree& tree = valuedTree->topology();
            const std::string label = (isMaxtree ? "random max-tree trial " : "random min-tree trial ") + std::to_string(trial);
            requireProductionMatchesExactOracle(tree, label);
            requireApproximateSummariesMatchExactSmallGrid(tree, label);

            const auto& altitude = valuedTree->nodeAltitudes();
            const std::vector<int> optimizedOracle =
                maxdist_oracle::componentTreeDiftMaxSquaredDistance(tree, std::span<const std::uint8_t>(altitude.data(), altitude.size()));
            const std::vector<std::int64_t> exactOracle = maxdist_oracle::exactMaxSquaredDistance(tree);
            for (NodeId node : tree.aliveNodeIds()) {
                requireEqual(optimizedOracle[static_cast<std::size_t>(node)], exactOracle[static_cast<std::size_t>(node)],
                             label + " optimized component-tree oracle");
            }
        }

        if (trial % 4 == 0) {
            for (TestTopographicImmersion immersion :
                 {TestTopographicImmersion::SelfDualSpan, TestTopographicImmersion::Min4Max8, TestTopographicImmersion::Min8Max4}) {
                auto treeOfShapes = makeTreeOfShapes(image, immersion);
                requireProductionMatchesExactOracle(*treeOfShapes, "random tree of shapes projected-support trial " + std::to_string(trial));
            }
        }
    }
}

void requireDeepRegionIndexIsIterative() {
    constexpr int numNodes = 4096;
    std::vector<NodeId> parent(static_cast<std::size_t>(numNodes));
    std::vector<std::uint8_t> altitude(static_cast<std::size_t>(numNodes), std::uint8_t{0});
    parent[0] = 0;
    for (NodeId node = 1; node < numNodes; ++node) {
        parent[static_cast<std::size_t>(node)] = node - 1;
    }
    const std::array<NodeId, 1> smallestNode{numNodes - 1};
    auto deepTree = MorphologicalTreeFactory::createFromNativeTopology(
        std::span<const NodeId>(parent), std::span<const NodeId>(smallestNode), std::span<const std::uint8_t>(altitude), NodeId{0}, 1, 1,
        MorphologicalTreeSemantics{MorphologicalTreeKind::Generic, NodeAltitudeOrder::Unconstrained, NoConstructionContext{}});

    RegionIndex regions(deepTree.topology());
    requireEqual(regions.numIndexedPixels(), std::size_t{1}, "deep iterative region index stores the pixel once");
    requireEqual(regions.supportInterval(0).size(), std::size_t{1}, "deep iterative region index root support");
    requireEqual(regions.supportInterval(numNodes - 1).size(), std::size_t{1}, "deep iterative region index leaf support");
    requireEqual(regions.boundingBox(0).area(), std::size_t{1}, "deep iterative region index root box");

    std::size_t frames = 0;
    attributes::computers::detail::distance_transform::MorphologicalTreeDistanceTransform::forEachNode(
        deepTree.topology(), [&frames](const attributes::computers::detail::distance_transform::NodeDistanceTransformFrame&) { ++frames; });
    requireEqual(frames, static_cast<std::size_t>(numNodes), "deep exact traversal visits every node iteratively");
}

void requireDeepOfflineBoundaryLcaIsIterative() {
    constexpr NodeId chainLeaf = 32768;
    constexpr NodeId incomparableSibling = chainLeaf + 1;
    constexpr int numNodes = incomparableSibling + 1;
    constexpr PixelId centerPixel = 4;

    std::vector<NodeId> parent(static_cast<std::size_t>(numNodes), NodeId{0});
    for (NodeId node = 2; node <= chainLeaf; ++node) {
        parent[static_cast<std::size_t>(node)] = node - 1;
    }
    parent[static_cast<std::size_t>(incomparableSibling)] = 0;

    std::array<NodeId, 9> smallestNode{};
    smallestNode.fill(incomparableSibling);
    smallestNode[static_cast<std::size_t>(centerPixel)] = chainLeaf;
    std::vector<std::uint8_t> altitude(static_cast<std::size_t>(numNodes), std::uint8_t{0});
    auto deepTree = MorphologicalTreeFactory::createFromNativeTopology(
        std::span<const NodeId>(parent), std::span<const NodeId>(smallestNode), std::span<const std::uint8_t>(altitude), NodeId{0}, 3, 3,
        MorphologicalTreeSemantics{MorphologicalTreeKind::Generic, NodeAltitudeOrder::Unconstrained, NoConstructionContext{}});
    const MorphologicalTree& tree = deepTree.topology();

    BoundaryLifetimeIndex lifetimes(tree);
    const auto& centerLifetime = lifetimes.lifetime(centerPixel);
    requireEqual(centerLifetime.owner, chainLeaf, "deep offline LCA preserves the center owner");
    requireEqual(centerLifetime.stopExclusive, tree.root(), "deep incomparable owners stop at their root LCA");
    requireVectorEqual(std::vector<PixelId>(lifetimes.additions(chainLeaf).begin(), lifetimes.additions(chainLeaf).end()), {centerPixel},
                       "deep offline LCA center addition event");
    requireVectorEqual(std::vector<PixelId>(lifetimes.removals(tree.root()).begin(), lifetimes.removals(tree.root()).end()), {centerPixel},
                       "deep offline LCA center removal event");

    RegionIndex regions(tree);
    std::uint64_t frames = 0;
    bool sawLeaf = false;
    bool sawMiddle = false;
    bool sawRoot = false;
    const NodeId middleNode = chainLeaf / 2;
    ContourScheduler::forEachNode(tree, regions, [&](const auto& frame) {
        ++frames;
        if (frame.node() == chainLeaf || frame.node() == middleNode) {
            requireVectorEqual(std::vector<PixelId>(frame.pixels().begin(), frame.pixels().end()), {centerPixel},
                               "deep chain contour remains the center site below the root");
            sawLeaf = sawLeaf || frame.node() == chainLeaf;
            sawMiddle = sawMiddle || frame.node() == middleNode;
        } else if (frame.node() == tree.root()) {
            std::vector<PixelId> actual(frame.pixels().begin(), frame.pixels().end());
            std::sort(actual.begin(), actual.end());
            requireVectorEqual(actual, {PixelId{0}, PixelId{1}, PixelId{2}, PixelId{3}, PixelId{5}, PixelId{6}, PixelId{7}, PixelId{8}},
                               "deep offline LCA root contour excludes its interior center");
            sawRoot = true;
        }
    });
    requireEqual(frames, static_cast<std::uint64_t>(numNodes), "deep offline LCA scheduler visits every live node iteratively");
    require(sawLeaf && sawMiddle && sawRoot, "deep offline LCA contour oracle observes the selected hierarchy frames");
}

void requireShiftedBoxRoundTrip() {
    Box2D domain(10, 20, 12, 21);
    requireEqual(domain.width(), 3, "shifted Box2D width");
    requireEqual(domain.height(), 2, "shifted Box2D height");

    for (int idx = 0; idx < domain.width() * domain.height(); ++idx) {
        const Point2D point = domain.point(idx);
        requireEqual(domain.index(point), idx, "shifted Box2D index(point(idx)) round-trip");
        requireEqual(domain.index(point.x(), point.y()), idx, "shifted Box2D index(x, y) round-trip");
    }

    requireEqual(domain.index(9, 20), -1, "shifted Box2D rejects x before top-left");
    requireEqual(domain.index(13, 21), -1, "shifted Box2D rejects x after bottom-right");
    requireEqual(domain.index(10, 19), -1, "shifted Box2D rejects y before top-left");
    requireEqual(domain.index(12, 22), -1, "shifted Box2D rejects y after bottom-right");
}

void requirePQueueLifoClearsSingletonBucket() {
    PQueue queue(16, 3);
    queue.setCost(0, 5);
    queue.insert(0);

    requireEqual(queue.popMinLIFO(), 0, "PQueue singleton LIFO pop");
    require(queue.isEmpty(), "PQueue singleton LIFO pop empties queue");

    queue.setCost(1, 5);
    queue.insert(1);
    requireEqual(queue.minElemFIFO(), 1, "PQueue insert after singleton LIFO must reset bucket head");
    requireEqual(queue.popMinFIFO(), 1, "PQueue FIFO after singleton LIFO reinsertion");
    require(queue.isEmpty(), "PQueue FIFO pop after singleton LIFO reinsertion empties queue");
}

void requirePQueueLifoMaintainsRemainingTail() {
    PQueue queue(16, 4);
    queue.setCost(0, 7);
    queue.setCost(1, 7);
    queue.insert(0);
    queue.insert(1);

    requireEqual(queue.popMinLIFO(), 1, "PQueue multi-element LIFO returns newest bucket element");
    require(!queue.isEmpty(), "PQueue multi-element LIFO leaves older bucket element queued");
    requireEqual(queue.minElemFIFO(), 0, "PQueue multi-element LIFO preserves bucket head");
    requireEqual(queue.popMinFIFO(), 0, "PQueue FIFO removes remaining element after LIFO");
    require(queue.isEmpty(), "PQueue LIFO then FIFO empties queue");

    queue.setCost(2, 7);
    queue.insert(2);
    requireEqual(queue.minElemFIFO(), 2, "PQueue insert after LIFO/FIFO sequence must reset bucket head");
    requireEqual(queue.popMinLIFO(), 2, "PQueue LIFO after LIFO/FIFO reinsertion");
    require(queue.isEmpty(), "PQueue LIFO after LIFO/FIFO reinsertion empties queue");
}

void requirePQueueMixedFifoLifoMaintainsLinks() {
    PQueue queue(16, 4);
    queue.setCost(0, 9);
    queue.setCost(1, 9);
    queue.insert(0);
    queue.insert(1);

    requireEqual(queue.popMinFIFO(), 0, "PQueue FIFO removes oldest bucket element");
    requireEqual(queue.popMinLIFO(), 1, "PQueue LIFO removes remaining element after FIFO");
    require(queue.isEmpty(), "PQueue FIFO then LIFO empties queue");

    queue.setCost(2, 9);
    queue.insert(2);
    requireEqual(queue.minElemFIFO(), 2, "PQueue insert after FIFO/LIFO sequence must reset bucket head");
}

void requirePQueueAcceptsMaxBucketCost() {
    PQueue queue(5, 2);
    queue.setCost(1, 5);
    queue.insert(1);

    requireEqual(queue.minValue(), 5, "PQueue max bucket cost min value");
    requireEqual(queue.maxValue(), 5, "PQueue max bucket cost max value");
    requireEqual(queue.minElemFIFO(), 1, "PQueue max bucket cost head");
    requireEqual(queue.popMinFIFO(), 1, "PQueue max bucket cost pop");
    require(queue.isEmpty(), "PQueue max bucket cost pop empties queue");
}

} // namespace

int main() {
    requireExactTransformMatchesBruteForce();
    requireAdaptiveApproximateProviderBoundary();
    requireUniformBackendAndOracles();
    requireRandomizedSmallGridAgreement();
    requireDeepRegionIndexIsIterative();
    requireDeepOfflineBoundaryLcaIsIterative();
    requireShiftedBoxRoundTrip();
    requirePQueueLifoClearsSingletonBucket();
    requirePQueueLifoMaintainsRemainingTail();
    requirePQueueMixedFifoLifoMaintainsLinks();
    requirePQueueAcceptsMaxBucketCost();
    return 0;
}
