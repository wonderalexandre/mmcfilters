#include <mmcfilters/attributes/Attributes.hpp>
#include <mmcfilters/trees/MorphologicalTreeFactory.hpp>
#include <mmcfilters/utils/Image.hpp>

#include <array>
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

template <class T, class U> void requireEqual(const T& actual, const U& expected, const std::string& message) {
    if (!(actual == expected)) {
        throw std::runtime_error(message);
    }
}

template <class T, class U, class V> void requireNear(T actual, U expected, V tolerance, const std::string& message) {
    if (std::abs(static_cast<double>(actual) - static_cast<double>(expected)) > static_cast<double>(tolerance)) {
        throw std::runtime_error(message);
    }
}

void requireHiddenDependenciesAreAbsent(const mmcfilters::AttributeNames& names, std::initializer_list<mmcfilters::Attribute> hiddenDependencies,
                                        const std::string& label) {
    for (mmcfilters::Attribute dependency : hiddenDependencies) {
        require(!names.contains(dependency), label + " must not expose hidden dependency");
    }
}

void requireFiniteLiveNodeValues(const mmcfilters::MorphologicalTree& tree, const mmcfilters::AttributeNames& names, const std::vector<float>& values,
                                 mmcfilters::Attribute attribute, const std::string& label) {
    for (mmcfilters::NodeId nodeId : tree.aliveNodeIds()) {
        const float value = values[names.linearIndex(nodeId, attribute)];
        require(std::isfinite(value), label + " must be finite for every live node");
    }
}

} // namespace

int main() {
    using namespace mmcfilters;

    std::array<std::uint8_t, 16> pixels{3, 3, 2, 2, 3, 4, 4, 2, 1, 4, 5, 2, 1, 1, 5, 0};
    auto image = ImageUInt8::fromExternal(pixels.data(), 4, 4);
    auto valuedTree = MorphologicalTreeFactory::createMaxTree(image, 1.5);

    auto valuedTreeAttributes = AttributeComputation::computeAttributes(valuedTree, std::vector<AttributeOrGroup>{Area, GrayLevelHeight, MaxDistExact});
    requireEqual(valuedTreeAttributes.first.NUM_ATTRIBUTES, 3, "valuedTree public attribute stride");
    require(valuedTreeAttributes.second.size() ==
                static_cast<std::size_t>(valuedTree.topology().numInternalNodeSlots()) * static_cast<std::size_t>(valuedTreeAttributes.first.NUM_ATTRIBUTES),
            "valuedTree public attribute buffer shape");
    requireEqual(valuedTreeAttributes.second[valuedTreeAttributes.first.linearIndex(valuedTree.topology().root(), Area)], 16.0f,
                 "root AREA through public valuedTree facade");

    auto topologyAttributes =
        AttributeComputation::computeTopologyAttributes(valuedTree.topology(), std::vector<AttributeOrGroup>{Area, BoxWidth, BalanceNode});
    requireEqual(topologyAttributes.first.NUM_ATTRIBUTES, 3, "topology public attribute stride");
    requireEqual(topologyAttributes.second[topologyAttributes.first.linearIndex(valuedTree.topology().root(), Area)], 16.0f,
                 "root AREA through explicit topology facade");

    auto topologyMaxDistExact = AttributeComputation::computeSingleTopologyAttribute(valuedTree.topology(), MaxDistExact);
    auto valuedMaxDistExact = AttributeComputation::computeSingleAttribute(valuedTree, MaxDistExact);
    require(topologyMaxDistExact.second == valuedMaxDistExact.second,
            "MAX_DIST_EXACT must use the same topology-only computation through topology and valued facades");

    auto topologyMaxDist = AttributeComputation::computeSingleTopologyAttribute(valuedTree.topology(), MaxDist);
    auto valuedMaxDist = AttributeComputation::computeSingleAttribute(valuedTree, MaxDist);
    require(topologyMaxDist.second == valuedMaxDist.second,
            "MAX_DIST must use the same topology-only computation through topology and valued facades");
    requireFiniteLiveNodeValues(valuedTree.topology(), topologyMaxDist.first, topologyMaxDist.second, MaxDist, "MAX_DIST request");

    auto bothDistanceAttributes = AttributeComputation::computeTopologyAttributes(
        valuedTree.topology(), std::vector<AttributeOrGroup>{MaxDistExact, MaxDist, MaxSquaredDistExact, MaxSquaredDist});
    requireEqual(bothDistanceAttributes.first.NUM_ATTRIBUTES, 4, "linear and squared maximum-distance public stride");
    require(bothDistanceAttributes.first.contains(MaxDistExact) && bothDistanceAttributes.first.contains(MaxDist) &&
                bothDistanceAttributes.first.contains(MaxSquaredDistExact) && bothDistanceAttributes.first.contains(MaxSquaredDist),
            "linear and squared exact/approximate maxima remain separate public attributes");
    for (NodeId node : valuedTree.topology().aliveNodeIds()) {
        const float approximate = bothDistanceAttributes.second[bothDistanceAttributes.first.linearIndex(node, MaxDist)];
        const float approximateSquared = bothDistanceAttributes.second[bothDistanceAttributes.first.linearIndex(node, MaxSquaredDist)];
        const float exact = bothDistanceAttributes.second[bothDistanceAttributes.first.linearIndex(node, MaxDistExact)];
        const float exactSquared = bothDistanceAttributes.second[bothDistanceAttributes.first.linearIndex(node, MaxSquaredDistExact)];
        requireNear(approximate * approximate, approximateSquared, 1.0e-5f, "approximate MAX_DIST squared contract");
        requireNear(exact * exact, exactSquared, 1.0e-5f, "exact MAX_DIST_EXACT squared contract");
    }

    const std::vector<AttributeOrGroup> distanceSummaryRequest{
        MaxDistExact, MaxSquaredDistExact, DistSquaredSumExact, DistSquaredMeanExact, DistRmsExact, DistSquaredVarianceExact,
        MaxDistCenterRowExact, MaxDistCenterColumnExact, MaxDist, MaxSquaredDist, DistSquaredSum, DistSquaredMean, DistRms,
        DistSquaredVariance, MaxDistCenterRow, MaxDistCenterColumn, MaxDistPlateauAreaExact,
        MaxDistPlateauCentroidRowExact, MaxDistPlateauCentroidColumnExact, MaxDistPlateauArea, MaxDistPlateauCentroidRow,
        MaxDistPlateauCentroidColumn,
        DistSum, DistMean, DistVariance, DistMedian, DistMode, DistQ25, DistQ75, DistQ90, DistEntropy,
        DistPositiveArea, DistLevelCount, DistWeightedCentroidRow, DistWeightedCentroidColumn,
        DistWeightedCentralMoment20, DistWeightedCentralMoment02, DistWeightedCentralMoment11,
        DistWeightedAxisOrientation, DistWeightedEccentricity,
        DistSumExact, DistMeanExact, DistVarianceExact, DistMedianExact, DistModeExact, DistQ25Exact,
        DistQ75Exact, DistQ90Exact, DistEntropyExact, DistPositiveAreaExact, DistLevelCountExact,
        DistWeightedCentroidRowExact, DistWeightedCentroidColumnExact, DistWeightedCentralMoment20Exact,
        DistWeightedCentralMoment02Exact, DistWeightedCentralMoment11Exact, DistWeightedAxisOrientationExact,
        DistWeightedEccentricityExact};
    auto distanceSummary = AttributeComputation::computeTopologyAttributes(valuedTree.topology(), distanceSummaryRequest);
    requireEqual(distanceSummary.first.NUM_ATTRIBUTES, 58, "exact and approximate distance-summary/geometry public stride");
    for (AttributeOrGroup item : distanceSummaryRequest) {
        const Attribute attribute = std::get<Attribute>(item);
        require(distanceSummary.first.contains(attribute), "distance-summary request exposes every requested scalar");
        requireFiniteLiveNodeValues(valuedTree.topology(), distanceSummary.first, distanceSummary.second, attribute,
                                    "distance-summary request");
    }
    auto approximateDistanceTransformGroup = AttributeComputation::computeTopologyAttributes(
        valuedTree.topology(), std::vector<AttributeOrGroup>{AttributeGroup::DistTransf});
    auto exactDistanceTransformGroup = AttributeComputation::computeTopologyAttributes(
        valuedTree.topology(), std::vector<AttributeOrGroup>{AttributeGroup::DistTransfExact});
    requireEqual(approximateDistanceTransformGroup.first.NUM_ATTRIBUTES, 29, "DIST_TRANSF public stride");
    requireEqual(exactDistanceTransformGroup.first.NUM_ATTRIBUTES, 29, "DIST_TRANSF_EXACT public stride");
    for (AttributeOrGroup item : distanceSummaryRequest) {
        const Attribute attribute = std::get<Attribute>(item);
        const bool exact = AttributeNames::toString(attribute).ends_with("_EXACT");
        require((exact ? exactDistanceTransformGroup.first : approximateDistanceTransformGroup.first).contains(attribute),
                exact ? "DIST_TRANSF_EXACT exposes every exact scalar" : "DIST_TRANSF exposes every approximate scalar");
    }

    auto valuedDistanceMean = AttributeComputation::computeSingleAttribute<double>(valuedTree, DistSquaredMeanExact);
    auto topologyDistanceMean = AttributeComputation::computeSingleTopologyAttribute<double>(valuedTree.topology(), DistSquaredMeanExact);
    require(valuedDistanceMean.second == topologyDistanceMean.second,
            "DIST_SQUARED_MEAN_EXACT must use the same topology-only computation through topology and valued facades");

    auto sampledArea = AttributeComputation::computeSampledNodeAttribute(valuedTree, Area, AltitudeDifference<std::uint8_t>{1}, 1);
    requireEqual(sampledArea.first.NUM_ATTRIBUTES, 3, "sampled public attribute stride");
    require(sampledArea.second.size() ==
                static_cast<std::size_t>(valuedTree.topology().numInternalNodeSlots()) * static_cast<std::size_t>(sampledArea.first.NUM_ATTRIBUTES),
            "sampled public attribute buffer shape");

    const NodeId root = valuedTree.topology().root();
    requireEqual(sampledArea.second[sampledArea.first.linearIndex(root, Area, 0)], 16.0f,
                 "current-node AREA sample through public valued-tree facade");

    auto meanOnly = AttributeComputation::computeAttributes(valuedTree, std::vector<AttributeOrGroup>{MeanGrayLevel});
    requireEqual(meanOnly.first.NUM_ATTRIBUTES, 1, "MeanGrayLevel request exposes one public attribute");
    require(meanOnly.first.contains(MeanGrayLevel), "MeanGrayLevel request exposes requested attribute");
    requireHiddenDependenciesAreAbsent(meanOnly.first, {Area, Volume}, "MeanGrayLevel request");
    requireFiniteLiveNodeValues(valuedTree.topology(), meanOnly.first, meanOnly.second, MeanGrayLevel, "MeanGrayLevel request");

    auto eccentricityOnly = AttributeComputation::computeAttributes(valuedTree, std::vector<AttributeOrGroup>{Eccentricity});
    requireEqual(eccentricityOnly.first.NUM_ATTRIBUTES, 1, "ECCENTRICITY valuedTree request exposes one public attribute");
    require(eccentricityOnly.first.contains(Eccentricity), "ECCENTRICITY valuedTree request exposes requested attribute");
    requireHiddenDependenciesAreAbsent(eccentricityOnly.first, {Area, CentralMoment20, CentralMoment02, CentralMoment11},
                                       "ECCENTRICITY valuedTree request");
    requireFiniteLiveNodeValues(valuedTree.topology(), eccentricityOnly.first, eccentricityOnly.second, Eccentricity, "ECCENTRICITY valuedTree request");

    auto topologyEccentricityOnly = AttributeComputation::computeTopologyAttributes(valuedTree.topology(), std::vector<AttributeOrGroup>{Eccentricity});
    requireEqual(topologyEccentricityOnly.first.NUM_ATTRIBUTES, 1, "ECCENTRICITY topology request exposes one public attribute");
    require(topologyEccentricityOnly.first.contains(Eccentricity), "ECCENTRICITY topology request exposes requested attribute");
    requireHiddenDependenciesAreAbsent(topologyEccentricityOnly.first, {Area, CentralMoment20, CentralMoment02, CentralMoment11},
                                       "ECCENTRICITY topology request");
    requireFiniteLiveNodeValues(valuedTree.topology(), topologyEccentricityOnly.first, topologyEccentricityOnly.second, Eccentricity,
                                "ECCENTRICITY topology request");

    auto huOnly = AttributeComputation::computeAttributes(valuedTree, std::vector<AttributeOrGroup>{HuMoment1});
    requireEqual(huOnly.first.NUM_ATTRIBUTES, 1, "HU_MOMENT_1 request exposes one public attribute");
    require(huOnly.first.contains(HuMoment1), "HU_MOMENT_1 request exposes requested attribute");
    requireHiddenDependenciesAreAbsent(
        huOnly.first,
        {Area, CentralMoment20, CentralMoment02, CentralMoment11, CentralMoment30, CentralMoment03, CentralMoment21, CentralMoment12},
        "HU_MOMENT_1 request");
    requireFiniteLiveNodeValues(valuedTree.topology(), huOnly.first, huOnly.second, HuMoment1, "HU_MOMENT_1 request");

    auto rectangularityOnly = AttributeComputation::computeTopologyAttributes(valuedTree.topology(), std::vector<AttributeOrGroup>{Rectangularity});
    requireEqual(rectangularityOnly.first.NUM_ATTRIBUTES, 1, "RECTANGULARITY topology request exposes one public attribute");
    require(rectangularityOnly.first.contains(Rectangularity), "RECTANGULARITY topology request exposes requested attribute");
    requireHiddenDependenciesAreAbsent(rectangularityOnly.first, {Area}, "RECTANGULARITY topology request");
    requireFiniteLiveNodeValues(valuedTree.topology(), rectangularityOnly.first, rectangularityOnly.second, Rectangularity, "RECTANGULARITY topology request");

    auto grayGroup = AttributeComputation::computeAttributes(valuedTree, std::vector<AttributeOrGroup>{AttributeGroup::GrayLevel});
    requireEqual(grayGroup.first.NUM_ATTRIBUTES, 5, "GRAY_LEVEL group public attribute stride");
    require(!grayGroup.first.contains(Area), "GRAY_LEVEL group must not expose hidden AREA dependency");

    return 0;
}
