#pragma once

#include "AttributeComputerDomain.hpp"
#include "AttributeComputerFamily.hpp"
#include "../detail/AttributeKernelSupport.hpp"
#include "detail/distance_transform_approx/MorphologicalTreeApproximateDistanceTransform.hpp"
#include "../../trees/MorphologicalTree.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <span>
#include <stdexcept>
#include <string_view>

/** @cond */
namespace mmcfilters::attributes::computers::detail {

inline void requireMaxDistCapabilities(const MorphologicalTree& tree) {
    if (!tree.hasGridDomain2D()) {
        throw std::invalid_argument("Approximate distance-field attributes require a regular 2D pixel domain.");
    }
}

struct ApproximateDistanceTransformAttributeRequest {
    bool maxDist = false;
    bool maxSquaredDist = false;
    bool squaredSum = false;
    bool squaredMean = false;
    bool rms = false;
    bool squaredVariance = false;
    bool centerRow = false;
    bool centerColumn = false;
    bool plateauArea = false;
    bool plateauCentroidRow = false;
    bool plateauCentroidColumn = false;
    bool distanceSum = false;
    bool distanceMean = false;
    bool distanceVariance = false;
    bool median = false;
    bool mode = false;
    bool q25 = false;
    bool q75 = false;
    bool q90 = false;
    bool entropy = false;
    bool positiveArea = false;
    bool levelCount = false;
    bool weightedCentroidRow = false;
    bool weightedCentroidColumn = false;
    bool weightedMu20 = false;
    bool weightedMu02 = false;
    bool weightedMu11 = false;
    bool weightedAxisOrientation = false;
    bool weightedEccentricity = false;

    [[nodiscard]] bool any() const noexcept { return maximum() || moments() || profile() || spatial() || geometry(); }
    [[nodiscard]] bool maximum() const noexcept { return maxDist || maxSquaredDist; }
    [[nodiscard]] bool moments() const noexcept {
        return squaredSum || squaredMean || rms || squaredVariance || distanceSum || distanceMean || distanceVariance;
    }
    [[nodiscard]] bool profile() const noexcept { return median || mode || q25 || q75 || q90 || entropy || positiveArea || levelCount; }
    [[nodiscard]] bool spatial() const noexcept {
        return weightedCentroidRow || weightedCentroidColumn || weightedMu20 || weightedMu02 || weightedMu11 || weightedAxisOrientation || weightedEccentricity;
    }
    [[nodiscard]] bool localization() const noexcept { return centerRow || centerColumn; }
    [[nodiscard]] bool plateau() const noexcept { return plateauArea || plateauCentroidRow || plateauCentroidColumn; }
    [[nodiscard]] bool geometry() const noexcept { return localization() || plateau(); }
    [[nodiscard]] bool needsCoordinates() const noexcept { return geometry() || spatial(); }

    [[nodiscard]] static ApproximateDistanceTransformAttributeRequest from(std::span<const Attribute> requestedAttributes) {
        const auto contains = [requestedAttributes](Attribute attribute) {
            return std::find(requestedAttributes.begin(), requestedAttributes.end(), attribute) != requestedAttributes.end();
        };
        return {.maxDist = contains(MaxDist),
                .maxSquaredDist = contains(MaxSquaredDist),
                .squaredSum = contains(DistSquaredSum),
                .squaredMean = contains(DistSquaredMean),
                .rms = contains(DistRms),
                .squaredVariance = contains(DistSquaredVariance),
                .centerRow = contains(MaxDistCenterRow),
                .centerColumn = contains(MaxDistCenterColumn),
                .plateauArea = contains(MaxDistPlateauArea),
                .plateauCentroidRow = contains(MaxDistPlateauCentroidRow),
                .plateauCentroidColumn = contains(MaxDistPlateauCentroidColumn),
                .distanceSum = contains(DistSum),
                .distanceMean = contains(DistMean),
                .distanceVariance = contains(DistVariance),
                .median = contains(DistMedian),
                .mode = contains(DistMode),
                .q25 = contains(DistQ25),
                .q75 = contains(DistQ75),
                .q90 = contains(DistQ90),
                .entropy = contains(DistEntropy),
                .positiveArea = contains(DistPositiveArea),
                .levelCount = contains(DistLevelCount),
                .weightedCentroidRow = contains(DistWeightedCentroidRow),
                .weightedCentroidColumn = contains(DistWeightedCentroidColumn),
                .weightedMu20 = contains(DistWeightedCentralMoment20),
                .weightedMu02 = contains(DistWeightedCentralMoment02),
                .weightedMu11 = contains(DistWeightedCentralMoment11),
                .weightedAxisOrientation = contains(DistWeightedAxisOrientation),
                .weightedEccentricity = contains(DistWeightedEccentricity)};
    }
};

template <std::floating_point Real, class SquaredDistance>
inline void materializeApproximateMaximum(const AttributeComputeContext<Real>& context, const ApproximateDistanceTransformAttributeRequest& request,
                                          NodeId node, SquaredDistance squaredDistance) {
    if (request.maxDist) {
        context.buffer[context.attrNames.linearIndex(node, MaxDist)] = static_cast<Real>(std::sqrt(static_cast<long double>(squaredDistance)));
    }
    if (request.maxSquaredDist) {
        context.buffer[context.attrNames.linearIndex(node, MaxSquaredDist)] = static_cast<Real>(squaredDistance);
    }
}

template <std::floating_point Real>
inline void materializeApproximateDistanceFieldExtremum(const AttributeComputeContext<Real>& context,
                                                        const ApproximateDistanceTransformAttributeRequest& request, NodeId node,
                                                        const distance_transform::DistanceFieldExtremum& extremum, int numColumns) {
    if (extremum.pixel == InvalidPixel) {
        throw std::logic_error("Approximate distance-field localization produced an empty live-node support.");
    }
    const auto materialize = [&context, node](Attribute attribute, auto value) {
        context.buffer[context.attrNames.linearIndex(node, attribute)] = static_cast<Real>(value);
    };
    materializeApproximateMaximum(context, request, node, extremum.squaredDistance);
    if (request.centerRow) {
        materialize(MaxDistCenterRow, extremum.pixel / numColumns);
    }
    if (request.centerColumn) {
        materialize(MaxDistCenterColumn, extremum.pixel % numColumns);
    }
}

template <std::floating_point Real>
inline void materializeApproximateDistanceFieldPlateau(const AttributeComputeContext<Real>& context,
                                                       const ApproximateDistanceTransformAttributeRequest& request, NodeId node,
                                                       const distance_transform::DistanceFieldMaximumPlateau& plateau, int numColumns) {
    if (plateau.pixel == InvalidPixel || plateau.count == 0) {
        throw std::logic_error("Approximate distance-field plateau reduction produced an empty live-node support.");
    }
    const auto materialize = [&context, node](Attribute attribute, auto value) {
        context.buffer[context.attrNames.linearIndex(node, attribute)] = static_cast<Real>(value);
    };
    materializeApproximateMaximum(context, request, node, plateau.squaredDistance);
    if (request.centerRow) {
        materialize(MaxDistCenterRow, plateau.pixel / numColumns);
    }
    if (request.centerColumn) {
        materialize(MaxDistCenterColumn, plateau.pixel % numColumns);
    }
    if (request.plateauArea) {
        materialize(MaxDistPlateauArea, plateau.count);
    }
    if (request.plateauCentroidRow) {
        materialize(MaxDistPlateauCentroidRow, plateau.centroidRow());
    }
    if (request.plateauCentroidColumn) {
        materialize(MaxDistPlateauCentroidColumn, plateau.centroidColumn());
    }
}

template <std::floating_point Real>
inline void materializeApproximateDistanceFieldMoments(const AttributeComputeContext<Real>& context,
                                                       const ApproximateDistanceTransformAttributeRequest& request, NodeId node,
                                                       const distance_transform::DistanceFieldMoments& moments) {
    const auto materialize = [&context, node](Attribute attribute, long double value) {
        context.buffer[context.attrNames.linearIndex(node, attribute)] = static_cast<Real>(value);
    };
    if (request.squaredSum) {
        materialize(DistSquaredSum, moments.sum());
    }
    if (request.squaredMean) {
        materialize(DistSquaredMean, moments.mean());
    }
    if (request.rms) {
        materialize(DistRms, moments.rms());
    }
    if (request.squaredVariance) {
        materialize(DistSquaredVariance, moments.populationVariance());
    }
    if (request.distanceSum) {
        materialize(DistSum, moments.distanceSum());
    }
    if (request.distanceMean) {
        materialize(DistMean, moments.distanceMean());
    }
    if (request.distanceVariance) {
        materialize(DistVariance, moments.distancePopulationVariance());
    }
}

template <std::floating_point Real>
inline void materializeApproximateDistanceFieldHistogram(const AttributeComputeContext<Real>& context,
                                                         const ApproximateDistanceTransformAttributeRequest& request, NodeId node,
                                                         const distance_transform::DistanceFieldHistogram& histogram) {
    const auto materialize = [&context, node](Attribute attribute, auto value) {
        context.buffer[context.attrNames.linearIndex(node, attribute)] = static_cast<Real>(value);
    };
    if (request.median) {
        materialize(DistMedian, histogram.quantile(0.5L));
    }
    if (request.mode) {
        materialize(DistMode, histogram.mode());
    }
    if (request.q25) {
        materialize(DistQ25, histogram.quantile(0.25L));
    }
    if (request.q75) {
        materialize(DistQ75, histogram.quantile(0.75L));
    }
    if (request.q90) {
        materialize(DistQ90, histogram.quantile(0.9L));
    }
    if (request.entropy) {
        materialize(DistEntropy, histogram.entropyBits());
    }
    if (request.positiveArea) {
        materialize(DistPositiveArea, histogram.positiveArea());
    }
    if (request.levelCount) {
        materialize(DistLevelCount, histogram.levelCount());
    }
}

template <std::floating_point Real>
inline void materializeApproximateDistanceFieldSpatialMoments(const AttributeComputeContext<Real>& context,
                                                              const ApproximateDistanceTransformAttributeRequest& request, NodeId node,
                                                              const distance_transform::DistanceWeightedSpatialMoments& spatialMoments) {
    const auto materialize = [&context, node](Attribute attribute, auto value) {
        context.buffer[context.attrNames.linearIndex(node, attribute)] = static_cast<Real>(value);
    };
    if (request.weightedCentroidRow) {
        materialize(DistWeightedCentroidRow, spatialMoments.centroidRow());
    }
    if (request.weightedCentroidColumn) {
        materialize(DistWeightedCentroidColumn, spatialMoments.centroidColumn());
    }
    if (request.weightedMu20) {
        materialize(DistWeightedCentralMoment20, spatialMoments.centralMoment20());
    }
    if (request.weightedMu02) {
        materialize(DistWeightedCentralMoment02, spatialMoments.centralMoment02());
    }
    if (request.weightedMu11) {
        materialize(DistWeightedCentralMoment11, spatialMoments.centralMoment11());
    }
    if (request.weightedAxisOrientation) {
        materialize(DistWeightedAxisOrientation, spatialMoments.axisOrientationDegrees());
    }
    if (request.weightedEccentricity) {
        materialize(DistWeightedEccentricity, spatialMoments.eccentricity());
    }
}

template <bool TrackMoments, bool TrackHistogram, bool TrackSpatialMoments, std::floating_point Real>
inline void computeSelectedApproximateDistanceTransformAttributes(const AttributeComputeContext<Real>& context,
                                                                  const ApproximateDistanceTransformAttributeRequest& request, int numColumns) {
    using Transform = distance_transform_approx::MorphologicalTreeApproximateDistanceTransform;
    const auto consume = [&context, request, numColumns](NodeId node, const distance_transform::DistanceFieldExtremum& extremum,
                                                         const distance_transform::DistanceFieldMaximumPlateau& plateau, const auto& observer,
                                                         PixelId representative) {
        if (request.plateau()) {
            materializeApproximateDistanceFieldPlateau(context, request, node, plateau, numColumns);
        } else if (request.localization()) {
            materializeApproximateDistanceFieldExtremum(context, request, node, extremum, numColumns);
        } else {
            materializeApproximateMaximum(context, request, node, extremum.squaredDistance);
        }
        if constexpr (TrackMoments) {
            materializeApproximateDistanceFieldMoments(context, request, node, observer.momentsFor(representative));
        }
        if constexpr (TrackHistogram) {
            materializeApproximateDistanceFieldHistogram(context, request, node, observer.histogramFor(representative));
        }
        if constexpr (TrackSpatialMoments) {
            materializeApproximateDistanceFieldSpatialMoments(context, request, node, observer.spatialMomentsFor(representative));
        }
    };

    if (request.plateau()) {
        Transform::template forEachNodeSelectedStatistics<TrackMoments, TrackHistogram, TrackSpatialMoments, false, true>(context.tree, consume);
        return;
    }
    if (request.localization()) {
        Transform::template forEachNodeSelectedStatistics<TrackMoments, TrackHistogram, TrackSpatialMoments, true, false>(context.tree, consume);
        return;
    }
    Transform::template forEachNodeSelectedStatistics<TrackMoments, TrackHistogram, TrackSpatialMoments>(context.tree, consume);
}

template <std::floating_point Real>
inline void dispatchSelectedApproximateDistanceTransformAttributes(const AttributeComputeContext<Real>& context,
                                                                   const ApproximateDistanceTransformAttributeRequest& request, int numColumns) {
    if (request.profile()) {
        if (request.spatial()) {
            if (request.moments()) {
                computeSelectedApproximateDistanceTransformAttributes<true, true, true>(context, request, numColumns);
                return;
            }
            computeSelectedApproximateDistanceTransformAttributes<false, true, true>(context, request, numColumns);
            return;
        }
        if (request.moments()) {
            computeSelectedApproximateDistanceTransformAttributes<true, true, false>(context, request, numColumns);
            return;
        }
        computeSelectedApproximateDistanceTransformAttributes<false, true, false>(context, request, numColumns);
        return;
    }
    if (request.moments()) {
        computeSelectedApproximateDistanceTransformAttributes<true, false, true>(context, request, numColumns);
        return;
    }
    computeSelectedApproximateDistanceTransformAttributes<false, false, true>(context, request, numColumns);
}

template <std::floating_point Real>
inline void computeApproximateDistanceTransformAttributes(const AttributeComputeContext<Real>& context,
                                                          const ApproximateDistanceTransformAttributeRequest& request) {
    using Transform = distance_transform_approx::MorphologicalTreeApproximateDistanceTransform;
    const int numColumns = request.needsCoordinates() ? context.tree.numColumns() : 0;
    if (request.profile() || request.spatial()) {
        dispatchSelectedApproximateDistanceTransformAttributes(context, request, numColumns);
        return;
    }
    if (request.moments()) {
        if (request.plateau()) {
            Transform::forEachNodeSummaryAndPlateau(context.tree,
                                                    [&context, request, numColumns](NodeId node, const distance_transform::DistanceFieldMaximumPlateau& plateau,
                                                                                    const distance_transform::DistanceFieldMoments& moments) {
                                                        materializeApproximateDistanceFieldPlateau(context, request, node, plateau, numColumns);
                                                        materializeApproximateDistanceFieldMoments(context, request, node, moments);
                                                    });
            return;
        }
        if (request.localization()) {
            Transform::forEachNodeSummaryAndExtremum(context.tree,
                                                     [&context, request, numColumns](NodeId node, const distance_transform::DistanceFieldExtremum& extremum,
                                                                                     const distance_transform::DistanceFieldMoments& moments) {
                                                         materializeApproximateDistanceFieldExtremum(context, request, node, extremum, numColumns);
                                                         materializeApproximateDistanceFieldMoments(context, request, node, moments);
                                                     });
            return;
        }
        Transform::forEachNodeSummary(context.tree, [&context, request](NodeId node, distance_transform_approx::ApproxSquaredDistance maximum,
                                                                        const distance_transform::DistanceFieldMoments& moments) {
            materializeApproximateMaximum(context, request, node, maximum);
            materializeApproximateDistanceFieldMoments(context, request, node, moments);
        });
        return;
    }
    if (request.plateau()) {
        Transform::forEachNodePlateau(context.tree,
                                      [&context, request, numColumns](NodeId node, const distance_transform::DistanceFieldMaximumPlateau& plateau) {
                                          materializeApproximateDistanceFieldPlateau(context, request, node, plateau, numColumns);
                                      });
        return;
    }
    if (request.localization()) {
        Transform::forEachNodeExtremum(context.tree, [&context, request, numColumns](NodeId node, const distance_transform::DistanceFieldExtremum& extremum) {
            materializeApproximateDistanceFieldExtremum(context, request, node, extremum, numColumns);
        });
        return;
    }
    Transform::forEachNodeMaximum(context.tree, [&context, request](NodeId node, distance_transform_approx::ApproxSquaredDistance value) {
        materializeApproximateMaximum(context, request, node, value);
    });
}

} // namespace mmcfilters::attributes::computers::detail
/** @endcond */

namespace mmcfilters::attributes::computers {

/** @brief Topology-only approximate distance-transform attribute computer. */
class MaxDistComputer {
  public:
    /// Family name used in dependency-plan diagnostics.
    static constexpr std::string_view familyName = "max-dist";

    /// Stable family id used by the scheduler.
    static constexpr AttributeComputerFamily family = AttributeComputerFamily::MaxDist;

    /// Approximate distance-transform attributes require topology/support, not altitudes.
    static constexpr AttributeComputerDomain domain = AttributeComputerDomain::Topology;

    /// Canonical list of scalar descriptors materialized by this computer.
    inline static constexpr std::array<Attribute, 29> producedAttributes{MaxDist,
                                                                         MaxSquaredDist,
                                                                         DistSquaredSum,
                                                                         DistSquaredMean,
                                                                         DistRms,
                                                                         DistSquaredVariance,
                                                                         MaxDistCenterRow,
                                                                         MaxDistCenterColumn,
                                                                         MaxDistPlateauArea,
                                                                         MaxDistPlateauCentroidRow,
                                                                         MaxDistPlateauCentroidColumn,
                                                                         DistSum,
                                                                         DistMean,
                                                                         DistVariance,
                                                                         DistMedian,
                                                                         DistMode,
                                                                         DistQ25,
                                                                         DistQ75,
                                                                         DistQ90,
                                                                         DistEntropy,
                                                                         DistPositiveArea,
                                                                         DistLevelCount,
                                                                         DistWeightedCentroidRow,
                                                                         DistWeightedCentroidColumn,
                                                                         DistWeightedCentralMoment20,
                                                                         DistWeightedCentralMoment02,
                                                                         DistWeightedCentralMoment11,
                                                                         DistWeightedAxisOrientation,
                                                                         DistWeightedEccentricity};

    /** @brief Validates the topology-only geometric capability contract. @param tree Tree to validate. */
    static void requireSupportedTreeKind(const MorphologicalTree& tree) { detail::requireMaxDistCapabilities(tree); }

    /** @brief Computes requested approximate distance-transform projections. @tparam Real Output scalar type. @param context Established compute context. */
    template <std::floating_point Real> static void compute(const AttributeComputeContext<Real>& context) {
        const detail::ApproximateDistanceTransformAttributeRequest request =
            detail::ApproximateDistanceTransformAttributeRequest::from(context.requestedAttributes);
        MMCFILTERS_CONTRACT_CHECKED_ONLY(requireAttributeBufferShape(context.tree, context.buffer, context.attrNames);
                                         requireRequestedAttributeColumns(context); detail::requireMaxDistCapabilities(context.tree));
        if (request.any()) {
            static_cast<void>(detail::computeApproximateDistanceTransformAttributes(context, request));
        }
    }

    /** @brief Materializes approximate unit-support values. @tparam Real Output scalar type. @param context Established unit-row context. */
    template <std::floating_point Real> static void computeUnitRows(const UnitAttributeComputeContext<Real>& context) {
        requireUnitAttributeBufferShape(context.tree, context.unitPixels, context.buffer, context.attrNames);
        const int numColumns = context.tree.numColumns();
        for (NodeId leafIndex = 0; leafIndex < static_cast<NodeId>(context.unitPixels.size()); ++leafIndex) {
            const PixelId pixel = context.unitPixels[static_cast<std::size_t>(leafIndex)];
            for (Attribute attribute : producedAttributes) {
                if (!requestsAttribute(context.requestedAttributes, attribute)) {
                    continue;
                }
                Real value = Real{0};
                if (attribute == MaxDistCenterRow) {
                    value = static_cast<Real>(pixel / numColumns);
                } else if (attribute == MaxDistCenterColumn) {
                    value = static_cast<Real>(pixel % numColumns);
                } else if (attribute == MaxDistPlateauArea) {
                    value = Real{1};
                } else if (attribute == MaxDistPlateauCentroidRow) {
                    value = static_cast<Real>(pixel / numColumns);
                } else if (attribute == MaxDistPlateauCentroidColumn) {
                    value = static_cast<Real>(pixel % numColumns);
                } else if (attribute == DistLevelCount) {
                    value = Real{1};
                } else if (attribute == DistWeightedCentroidRow) {
                    value = static_cast<Real>(pixel / numColumns);
                } else if (attribute == DistWeightedCentroidColumn) {
                    value = static_cast<Real>(pixel % numColumns);
                } else if (attribute == DistWeightedEccentricity) {
                    value = Real{1};
                }
                context.buffer[context.attrNames.linearIndex(leafIndex, attribute)] = value;
            }
        }
    }
};

} // namespace mmcfilters::attributes::computers
