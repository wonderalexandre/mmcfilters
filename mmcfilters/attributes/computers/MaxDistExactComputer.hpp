#pragma once

#include "AttributeComputerDomain.hpp"
#include "AttributeComputerFamily.hpp"
#include "../detail/AttributeKernelSupport.hpp"
#include "detail/distance_transform/DistanceFieldMoments.hpp"
#include "detail/distance_transform/DistanceFieldHistogram.hpp"
#include "detail/distance_transform/DistanceWeightedSpatialMoments.hpp"
#include "detail/distance_transform/MorphologicalTreeDistanceTransform.hpp"
#include "../../trees/MorphologicalTree.hpp"
#include "../../utils/Common.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <span>
#include <stdexcept>
#include <string_view>

/** @cond */
namespace mmcfilters::attributes::computers::detail {

/** @brief Requested scalar projections of the shared node EDT family. */
struct DistanceTransformAttributeRequest {
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

    /** @brief Tests whether at least one distance-transform projection is requested. */
    [[nodiscard]] bool any() const noexcept { return maximum() || moments() || profile() || spatial() || geometry(); }

    /** @brief Tests whether either scalar projection of the maximum is requested. */
    [[nodiscard]] bool maximum() const noexcept { return maxDist || maxSquaredDist; }

    /** @brief Tests whether the full additive distance distribution is needed. */
    [[nodiscard]] bool moments() const noexcept {
        return squaredSum || squaredMean || rms || squaredVariance || distanceSum || distanceMean || distanceVariance;
    }

    /** @brief Tests whether an ordered distance histogram is required. */
    [[nodiscard]] bool profile() const noexcept { return median || mode || q25 || q75 || q90 || entropy || positiveArea || levelCount; }

    /** @brief Tests whether distance-weighted spatial moments are required. */
    [[nodiscard]] bool spatial() const noexcept {
        return weightedCentroidRow || weightedCentroidColumn || weightedMu20 || weightedMu02 || weightedMu11 || weightedAxisOrientation || weightedEccentricity;
    }

    /** @brief Tests whether the row-major maximum location is needed. */
    [[nodiscard]] bool localization() const noexcept { return centerRow || centerColumn; }

    /** @brief Tests whether the complete maximum-distance plateau is needed. */
    [[nodiscard]] bool plateau() const noexcept { return plateauArea || plateauCentroidRow || plateauCentroidColumn; }

    /** @brief Tests whether any maximum-distance geometry is needed. */
    [[nodiscard]] bool geometry() const noexcept { return localization() || plateau(); }
    [[nodiscard]] bool needsCoordinates() const noexcept { return geometry() || spatial(); }

    /** @brief Builds the reducer request from a public attribute subset. */
    [[nodiscard]] static DistanceTransformAttributeRequest from(std::span<const Attribute> requestedAttributes) {
        const auto contains = [requestedAttributes](Attribute attribute) {
            return std::find(requestedAttributes.begin(), requestedAttributes.end(), attribute) != requestedAttributes.end();
        };
        return {.maxDist = contains(MaxDistExact),
                .maxSquaredDist = contains(MaxSquaredDistExact),
                .squaredSum = contains(DistSquaredSumExact),
                .squaredMean = contains(DistSquaredMeanExact),
                .rms = contains(DistRmsExact),
                .squaredVariance = contains(DistSquaredVarianceExact),
                .centerRow = contains(MaxDistCenterRowExact),
                .centerColumn = contains(MaxDistCenterColumnExact),
                .plateauArea = contains(MaxDistPlateauAreaExact),
                .plateauCentroidRow = contains(MaxDistPlateauCentroidRowExact),
                .plateauCentroidColumn = contains(MaxDistPlateauCentroidColumnExact),
                .distanceSum = contains(DistSumExact),
                .distanceMean = contains(DistMeanExact),
                .distanceVariance = contains(DistVarianceExact),
                .median = contains(DistMedianExact),
                .mode = contains(DistModeExact),
                .q25 = contains(DistQ25Exact),
                .q75 = contains(DistQ75Exact),
                .q90 = contains(DistQ90Exact),
                .entropy = contains(DistEntropyExact),
                .positiveArea = contains(DistPositiveAreaExact),
                .levelCount = contains(DistLevelCountExact),
                .weightedCentroidRow = contains(DistWeightedCentroidRowExact),
                .weightedCentroidColumn = contains(DistWeightedCentroidColumnExact),
                .weightedMu20 = contains(DistWeightedCentralMoment20Exact),
                .weightedMu02 = contains(DistWeightedCentralMoment02Exact),
                .weightedMu11 = contains(DistWeightedCentralMoment11Exact),
                .weightedAxisOrientation = contains(DistWeightedAxisOrientationExact),
                .weightedEccentricity = contains(DistWeightedEccentricityExact)};
    }
};

/**
 * @brief Validates the topology-only geometric contract of EDT attributes.
 *
 * Tree kind, altitude order, node altitudes, and construction adjacency are
 * deliberately irrelevant. Every established MorphologicalTree with a regular
 * 2D domain follows the same production distance-transform computation.
 */
inline void requireMaxDistExactCapabilities(const MorphologicalTree& tree) {
    if (!tree.hasGridDomain2D()) {
        throw std::invalid_argument("Exact distance-field attributes require a regular 2D pixel domain.");
    }
}

/** @brief Projects the exact node EDT stream to its maximum squared distance. */
template <std::floating_point Real> class MaxSquaredDistanceAttributeReducer {
  public:
    MaxSquaredDistanceAttributeReducer(std::span<Real> buffer, const AttributeNames& attrNames, DistanceTransformAttributeRequest request)
        : buffer_(buffer), attrNames_(attrNames), request_(request) {}

    void beginNode(NodeId) { currentMaximum_ = 0; }

    void consumeSample(PixelId, distance_transform::SquaredDistance squaredDistance) { currentMaximum_ = std::max(currentMaximum_, squaredDistance); }

    void endNode(NodeId node) {
        if (request_.maxDist) {
            materialize(node, MaxDistExact, std::sqrt(static_cast<long double>(currentMaximum_)));
        }
        if (request_.maxSquaredDist) {
            materialize(node, MaxSquaredDistExact, currentMaximum_);
        }
    }

  private:
    template <class Value> void materialize(NodeId node, Attribute attribute, Value value) {
        buffer_[attrNames_.linearIndex(node, attribute)] = static_cast<Real>(value);
    }

    std::span<Real> buffer_;
    const AttributeNames& attrNames_;
    DistanceTransformAttributeRequest request_;
    distance_transform::SquaredDistance currentMaximum_ = 0;
};

/** @brief Projects the exact node EDT stream to its maximum and canonical location. */
template <std::floating_point Real> class DistanceFieldExtremumAttributeReducer {
  public:
    DistanceFieldExtremumAttributeReducer(std::span<Real> buffer, const AttributeNames& attrNames, DistanceTransformAttributeRequest request, int numColumns)
        : buffer_(buffer), attrNames_(attrNames), request_(request), numColumns_(numColumns) {}

    void beginNode(NodeId) { extremum_ = {}; }

    void consumeSample(PixelId pixel, distance_transform::SquaredDistance squaredDistance) {
        distance_transform::updateDistanceFieldExtremum(extremum_, pixel, squaredDistance);
    }

    void endNode(NodeId node) {
        requireLocation();
        if (request_.maxDist) {
            materialize(node, MaxDistExact, std::sqrt(static_cast<long double>(extremum_.squaredDistance)));
        }
        if (request_.maxSquaredDist) {
            materialize(node, MaxSquaredDistExact, extremum_.squaredDistance);
        }
        if (request_.centerRow) {
            materialize(node, MaxDistCenterRowExact, extremum_.pixel / numColumns_);
        }
        if (request_.centerColumn) {
            materialize(node, MaxDistCenterColumnExact, extremum_.pixel % numColumns_);
        }
    }

  private:
    void requireLocation() const {
        if (extremum_.pixel == InvalidPixel) {
            throw std::logic_error("Exact distance-field localization produced an empty live-node support.");
        }
    }

    template <class Value> void materialize(NodeId node, Attribute attribute, Value value) {
        buffer_[attrNames_.linearIndex(node, attribute)] = static_cast<Real>(value);
    }

    std::span<Real> buffer_;
    const AttributeNames& attrNames_;
    DistanceTransformAttributeRequest request_;
    int numColumns_ = 0;
    distance_transform::DistanceFieldExtremum extremum_;
};

/** @brief Projects the exact node EDT stream to the complete maximum plateau. */
template <std::floating_point Real> class DistanceFieldPlateauAttributeReducer {
  public:
    DistanceFieldPlateauAttributeReducer(std::span<Real> buffer, const AttributeNames& attrNames, DistanceTransformAttributeRequest request, int numColumns)
        : buffer_(buffer), attrNames_(attrNames), request_(request), numColumns_(numColumns) {}

    void beginNode(NodeId) { plateau_ = {}; }

    void consumeSample(PixelId pixel, distance_transform::SquaredDistance squaredDistance) {
        distance_transform::updateDistanceFieldMaximumPlateau(plateau_, pixel, squaredDistance, numColumns_);
    }

    void endNode(NodeId node) {
        requirePlateau();
        if (request_.maxDist) {
            materialize(node, MaxDistExact, std::sqrt(static_cast<long double>(plateau_.squaredDistance)));
        }
        if (request_.maxSquaredDist) {
            materialize(node, MaxSquaredDistExact, plateau_.squaredDistance);
        }
        if (request_.centerRow) {
            materialize(node, MaxDistCenterRowExact, plateau_.pixel / numColumns_);
        }
        if (request_.centerColumn) {
            materialize(node, MaxDistCenterColumnExact, plateau_.pixel % numColumns_);
        }
        if (request_.plateauArea) {
            materialize(node, MaxDistPlateauAreaExact, plateau_.count);
        }
        if (request_.plateauCentroidRow) {
            materialize(node, MaxDistPlateauCentroidRowExact, plateau_.centroidRow());
        }
        if (request_.plateauCentroidColumn) {
            materialize(node, MaxDistPlateauCentroidColumnExact, plateau_.centroidColumn());
        }
    }

  private:
    void requirePlateau() const {
        if (plateau_.pixel == InvalidPixel || plateau_.count == 0) {
            throw std::logic_error("Exact distance-field plateau reduction produced an empty live-node support.");
        }
    }

    template <class Value> void materialize(NodeId node, Attribute attribute, Value value) {
        buffer_[attrNames_.linearIndex(node, attribute)] = static_cast<Real>(value);
    }

    std::span<Real> buffer_;
    const AttributeNames& attrNames_;
    DistanceTransformAttributeRequest request_;
    int numColumns_ = 0;
    distance_transform::DistanceFieldMaximumPlateau plateau_;
};

/** @brief Projects one exact node EDT stream to all requested scalar moments. */
template <std::floating_point Real, bool TrackLocation, bool TrackPlateau> class DistanceFieldSummaryAttributeReducer {
  public:
    DistanceFieldSummaryAttributeReducer(std::span<Real> buffer, const AttributeNames& attrNames, DistanceTransformAttributeRequest request, int numColumns)
        : buffer_(buffer), attrNames_(attrNames), request_(request), numColumns_(numColumns) {}

    void beginNode(NodeId) {
        if (request_.moments()) {
            moments_.clear();
        }
        if (request_.profile()) {
            histogram_.clear();
        }
        if (request_.spatial()) {
            spatialMoments_.clear();
        }
        maximum_ = 0;
        if constexpr (TrackPlateau) {
            plateau_ = {};
        } else if constexpr (TrackLocation) {
            extremum_ = {};
        }
    }

    void consumeSample(PixelId pixel, distance_transform::SquaredDistance squaredDistance) {
        if (request_.moments()) {
            moments_.add<false>(squaredDistance);
        }
        if (request_.profile()) {
            histogram_.add<false>(squaredDistance);
        }
        if (request_.spatial()) {
            spatialMoments_.add<false>(pixel, squaredDistance, numColumns_);
        }
        maximum_ = std::max(maximum_, squaredDistance);
        if constexpr (TrackPlateau) {
            distance_transform::updateDistanceFieldMaximumPlateau(plateau_, pixel, squaredDistance, numColumns_);
        } else if constexpr (TrackLocation) {
            distance_transform::updateDistanceFieldExtremum(extremum_, pixel, squaredDistance);
        }
    }

    void endNode(NodeId node) {
        const std::uint64_t sampleCount = request_.moments() ? moments_.count() : (request_.profile() ? histogram_.count() : spatialMoments_.count());
        if (sampleCount == 0) {
            throw std::logic_error("Exact distance-field reduction produced an empty live-node support.");
        }
        if (request_.maxDist) {
            materialize(node, MaxDistExact, std::sqrt(static_cast<long double>(maximum_)));
        }
        if (request_.maxSquaredDist) {
            materialize(node, MaxSquaredDistExact, static_cast<long double>(maximum_));
        }
        if (request_.squaredSum) {
            materialize(node, DistSquaredSumExact, moments_.sum());
        }
        if (request_.squaredMean) {
            materialize(node, DistSquaredMeanExact, moments_.mean());
        }
        if (request_.rms) {
            materialize(node, DistRmsExact, moments_.rms());
        }
        if (request_.squaredVariance) {
            materialize(node, DistSquaredVarianceExact, moments_.populationVariance());
        }
        if (request_.distanceSum) {
            materialize(node, DistSumExact, moments_.distanceSum());
        }
        if (request_.distanceMean) {
            materialize(node, DistMeanExact, moments_.distanceMean());
        }
        if (request_.distanceVariance) {
            materialize(node, DistVarianceExact, moments_.distancePopulationVariance());
        }
        if (request_.median) {
            materialize(node, DistMedianExact, histogram_.quantile(0.5L));
        }
        if (request_.mode) {
            materialize(node, DistModeExact, histogram_.mode());
        }
        if (request_.q25) {
            materialize(node, DistQ25Exact, histogram_.quantile(0.25L));
        }
        if (request_.q75) {
            materialize(node, DistQ75Exact, histogram_.quantile(0.75L));
        }
        if (request_.q90) {
            materialize(node, DistQ90Exact, histogram_.quantile(0.9L));
        }
        if (request_.entropy) {
            materialize(node, DistEntropyExact, histogram_.entropyBits());
        }
        if (request_.positiveArea) {
            materialize(node, DistPositiveAreaExact, histogram_.positiveArea());
        }
        if (request_.levelCount) {
            materialize(node, DistLevelCountExact, histogram_.levelCount());
        }
        if (request_.weightedCentroidRow) {
            materialize(node, DistWeightedCentroidRowExact, spatialMoments_.centroidRow());
        }
        if (request_.weightedCentroidColumn) {
            materialize(node, DistWeightedCentroidColumnExact, spatialMoments_.centroidColumn());
        }
        if (request_.weightedMu20) {
            materialize(node, DistWeightedCentralMoment20Exact, spatialMoments_.centralMoment20());
        }
        if (request_.weightedMu02) {
            materialize(node, DistWeightedCentralMoment02Exact, spatialMoments_.centralMoment02());
        }
        if (request_.weightedMu11) {
            materialize(node, DistWeightedCentralMoment11Exact, spatialMoments_.centralMoment11());
        }
        if (request_.weightedAxisOrientation) {
            materialize(node, DistWeightedAxisOrientationExact, spatialMoments_.axisOrientationDegrees());
        }
        if (request_.weightedEccentricity) {
            materialize(node, DistWeightedEccentricityExact, spatialMoments_.eccentricity());
        }
        if constexpr (TrackPlateau) {
            if (plateau_.pixel == InvalidPixel || plateau_.count == 0) {
                throw std::logic_error("Exact distance-field summary produced an empty maximum plateau.");
            }
            if (request_.centerRow) {
                materialize(node, MaxDistCenterRowExact, static_cast<long double>(plateau_.pixel / numColumns_));
            }
            if (request_.centerColumn) {
                materialize(node, MaxDistCenterColumnExact, static_cast<long double>(plateau_.pixel % numColumns_));
            }
            if (request_.plateauArea) {
                materialize(node, MaxDistPlateauAreaExact, static_cast<long double>(plateau_.count));
            }
            if (request_.plateauCentroidRow) {
                materialize(node, MaxDistPlateauCentroidRowExact, plateau_.centroidRow());
            }
            if (request_.plateauCentroidColumn) {
                materialize(node, MaxDistPlateauCentroidColumnExact, plateau_.centroidColumn());
            }
        } else if constexpr (TrackLocation) {
            if (extremum_.pixel == InvalidPixel) {
                throw std::logic_error("Exact distance-field summary produced an empty maximum location.");
            }
            if (request_.centerRow) {
                materialize(node, MaxDistCenterRowExact, static_cast<long double>(extremum_.pixel / numColumns_));
            }
            if (request_.centerColumn) {
                materialize(node, MaxDistCenterColumnExact, static_cast<long double>(extremum_.pixel % numColumns_));
            }
        }
    }

  private:
    void materialize(NodeId node, Attribute attribute, long double value) { buffer_[attrNames_.linearIndex(node, attribute)] = static_cast<Real>(value); }

    std::span<Real> buffer_;
    const AttributeNames& attrNames_;
    DistanceTransformAttributeRequest request_;
    int numColumns_ = 0;
    distance_transform::DistanceFieldMoments moments_;
    distance_transform::DistanceFieldHistogram histogram_;
    distance_transform::DistanceWeightedSpatialMoments spatialMoments_;
    distance_transform::SquaredDistance maximum_ = 0;
    distance_transform::DistanceFieldExtremum extremum_;
    distance_transform::DistanceFieldMaximumPlateau plateau_;
};

namespace kernel {

/** @brief Computes all requested projections through one shared EDT traversal. */
template <std::floating_point Real>
inline void computeDistanceTransformAttributes(const AttributeComputeContext<Real>& context, const DistanceTransformAttributeRequest& request, int numColumns) {
    if (!request.any()) {
        return;
    }

    if (request.moments() || request.profile() || request.spatial()) {
        if (request.plateau()) {
            DistanceFieldSummaryAttributeReducer<Real, false, true> summaryReducer(context.buffer, context.attrNames, request, numColumns);
            distance_transform::ExactNodeDistanceFieldProvider::reduce(context.tree, summaryReducer);
        } else if (request.localization()) {
            DistanceFieldSummaryAttributeReducer<Real, true, false> summaryReducer(context.buffer, context.attrNames, request, numColumns);
            distance_transform::ExactNodeDistanceFieldProvider::reduce(context.tree, summaryReducer);
        } else {
            DistanceFieldSummaryAttributeReducer<Real, false, false> summaryReducer(context.buffer, context.attrNames, request, numColumns);
            distance_transform::ExactNodeDistanceFieldProvider::reduce(context.tree, summaryReducer);
        }
    } else if (request.plateau()) {
        DistanceFieldPlateauAttributeReducer<Real> plateauReducer(context.buffer, context.attrNames, request, numColumns);
        distance_transform::ExactNodeDistanceFieldProvider::reduce(context.tree, plateauReducer);
    } else if (request.localization()) {
        DistanceFieldExtremumAttributeReducer<Real> extremumReducer(context.buffer, context.attrNames, request, numColumns);
        distance_transform::ExactNodeDistanceFieldProvider::reduce(context.tree, extremumReducer);
    } else if (request.maximum()) {
        MaxSquaredDistanceAttributeReducer<Real> maximumReducer(context.buffer, context.attrNames, request);
        distance_transform::ExactNodeDistanceFieldProvider::reduce(context.tree, maximumReducer);
    }
}

} // namespace kernel

} // namespace mmcfilters::attributes::computers::detail
/** @endcond */

namespace mmcfilters::attributes::computers {

/**
 * @brief Topology-only scalar computer for attributes derived from the node EDT.
 *
 * Requested exact distance-derived projections share a single transform
 * traversal and one support enumeration per node.
 */
class MaxDistExactComputer {
  public:
    /// Family name used in dependency-plan diagnostics.
    static constexpr std::string_view familyName = "max-dist-exact";

    /// Stable family id used by the scheduler.
    static constexpr AttributeComputerFamily family = AttributeComputerFamily::MaxDistExact;

    /// Distance-transform attributes require topology/support, not altitudes.
    static constexpr AttributeComputerDomain domain = AttributeComputerDomain::Topology;

    /// Canonical list of scalar descriptors materialized by this computer.
    inline static constexpr std::array<Attribute, 29> producedAttributes{MaxDistExact,
                                                                         MaxSquaredDistExact,
                                                                         DistSquaredSumExact,
                                                                         DistSquaredMeanExact,
                                                                         DistRmsExact,
                                                                         DistSquaredVarianceExact,
                                                                         MaxDistCenterRowExact,
                                                                         MaxDistCenterColumnExact,
                                                                         MaxDistPlateauAreaExact,
                                                                         MaxDistPlateauCentroidRowExact,
                                                                         MaxDistPlateauCentroidColumnExact,
                                                                         DistSumExact,
                                                                         DistMeanExact,
                                                                         DistVarianceExact,
                                                                         DistMedianExact,
                                                                         DistModeExact,
                                                                         DistQ25Exact,
                                                                         DistQ75Exact,
                                                                         DistQ90Exact,
                                                                         DistEntropyExact,
                                                                         DistPositiveAreaExact,
                                                                         DistLevelCountExact,
                                                                         DistWeightedCentroidRowExact,
                                                                         DistWeightedCentroidColumnExact,
                                                                         DistWeightedCentralMoment20Exact,
                                                                         DistWeightedCentralMoment02Exact,
                                                                         DistWeightedCentralMoment11Exact,
                                                                         DistWeightedAxisOrientationExact,
                                                                         DistWeightedEccentricityExact};

    /** @brief Validates the topology-only geometric capability contract. @param tree Tree to validate. */
    static void requireSupportedTreeKind(const MorphologicalTree& tree) { detail::requireMaxDistExactCapabilities(tree); }

    /** @brief Computes requested distance-transform projections in one pass. @tparam Real Output scalar type. @param context Established compute context. */
    template <std::floating_point Real> static void compute(const AttributeComputeContext<Real>& context) {
        const detail::DistanceTransformAttributeRequest request = detail::DistanceTransformAttributeRequest::from(context.requestedAttributes);
        MMCFILTERS_CONTRACT_CHECKED_ONLY(requireAttributeBufferShape(context.tree, context.buffer, context.attrNames);
                                         requireRequestedAttributeColumns(context); detail::requireMaxDistExactCapabilities(context.tree));
        const int numColumns = request.needsCoordinates() ? context.tree.numColumns() : 0;
        detail::kernel::computeDistanceTransformAttributes(context, request, numColumns);
    }

    /** @brief Materializes zero for one-pixel unit supports. @tparam Real Output scalar type. @param context Established unit-row context. */
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
                if (attribute == MaxDistCenterRowExact) {
                    value = static_cast<Real>(pixel / numColumns);
                } else if (attribute == MaxDistCenterColumnExact) {
                    value = static_cast<Real>(pixel % numColumns);
                } else if (attribute == MaxDistPlateauAreaExact) {
                    value = Real{1};
                } else if (attribute == MaxDistPlateauCentroidRowExact) {
                    value = static_cast<Real>(pixel / numColumns);
                } else if (attribute == MaxDistPlateauCentroidColumnExact) {
                    value = static_cast<Real>(pixel % numColumns);
                } else if (attribute == DistLevelCountExact) {
                    value = Real{1};
                } else if (attribute == DistWeightedCentroidRowExact) {
                    value = static_cast<Real>(pixel / numColumns);
                } else if (attribute == DistWeightedCentroidColumnExact) {
                    value = static_cast<Real>(pixel % numColumns);
                } else if (attribute == DistWeightedEccentricityExact) {
                    value = Real{1};
                }
                context.buffer[context.attrNames.linearIndex(leafIndex, attribute)] = value;
            }
        }
    }
};

} // namespace mmcfilters::attributes::computers
