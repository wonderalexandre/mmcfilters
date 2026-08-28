#include "support/TestSupport.hpp"
#include "attributes/support/maxdist/MaxDistOracles.hpp"

#include "mmcfilters/attributes/Attributes.hpp"
#include "mmcfilters/attributes/computers/detail/distance_transform_approx/MorphologicalTreeApproximateDistanceTransform.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <numbers>
#include <span>
#include <string>
#include <vector>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

namespace {

constexpr double Tolerance = 1.0e-12;

void requireRootStudy(const MorphologicalTree& tree, bool exact) {
    const std::vector<AttributeOrGroup> attributes = exact ? std::vector<AttributeOrGroup>{MaxDistExact,
                                                                                           MaxSquaredDistExact,
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
                                                                                           DistWeightedEccentricityExact}
                                                           : std::vector<AttributeOrGroup>{MaxDist,
                                                                                           MaxSquaredDist,
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
    const auto [names, values] = AttributeComputation::computeTopologyAttributes<double>(tree, attributes);
    const NodeId root = tree.root();
    const auto value = [&](Attribute approximate, Attribute exactAttribute) { return values[names.linearIndex(root, exact ? exactAttribute : approximate)]; };

    const double expectedEntropy = -(16.0 / 25.0) * std::log2(16.0 / 25.0) - (8.0 / 25.0) * std::log2(8.0 / 25.0) - (1.0 / 25.0) * std::log2(1.0 / 25.0);
    requireNear(value(MaxDist, MaxDistExact), 2.0, Tolerance, "maximum distance in pixels");
    requireNear(value(MaxSquaredDist, MaxSquaredDistExact), 4.0, Tolerance, "maximum squared distance in squared pixels");
    requireNear(value(DistSum, DistSumExact), 10.0, Tolerance, "distance sum");
    requireNear(value(DistMean, DistMeanExact), 10.0 / 25.0, Tolerance, "distance mean");
    requireNear(value(DistVariance, DistVarianceExact), 12.0 / 25.0 - 0.16, Tolerance, "distance variance");
    requireNear(value(DistMedian, DistMedianExact), 0.0, Tolerance, "distance median");
    requireNear(value(DistMode, DistModeExact), 0.0, Tolerance, "distance mode");
    requireNear(value(DistQ25, DistQ25Exact), 0.0, Tolerance, "distance q25");
    requireNear(value(DistQ75, DistQ75Exact), 1.0, Tolerance, "distance q75");
    requireNear(value(DistQ90, DistQ90Exact), 1.0, Tolerance, "distance q90");
    requireNear(value(DistEntropy, DistEntropyExact), expectedEntropy, Tolerance, "distance entropy");
    requireEqual(value(DistPositiveArea, DistPositiveAreaExact), 9.0, "positive-distance area");
    requireEqual(value(DistLevelCount, DistLevelCountExact), 3.0, "distance level count");
    requireNear(value(DistWeightedCentroidRow, DistWeightedCentroidRowExact), 2.0, Tolerance, "weighted centroid row");
    requireNear(value(DistWeightedCentroidColumn, DistWeightedCentroidColumnExact), 2.0, Tolerance, "weighted centroid column");
    requireNear(value(DistWeightedCentralMoment20, DistWeightedCentralMoment20Exact), 6.0, Tolerance, "weighted mu20");
    requireNear(value(DistWeightedCentralMoment02, DistWeightedCentralMoment02Exact), 6.0, Tolerance, "weighted mu02");
    requireNear(value(DistWeightedCentralMoment11, DistWeightedCentralMoment11Exact), 0.0, Tolerance, "weighted mu11");
    requireNear(value(DistWeightedAxisOrientation, DistWeightedAxisOrientationExact), 0.0, Tolerance, "weighted orientation");
    requireNear(value(DistWeightedEccentricity, DistWeightedEccentricityExact), 1.0, Tolerance, "weighted eccentricity");
}

void requireWeightedMomentAxisConvention() {
    auto image = ImageUInt8::create(3, 5);
    image->fill(7);
    const auto tree = makeComponentTree(image, true);
    const auto [names, values] = AttributeComputation::computeTopologyAttributes<double>(
        *tree, std::vector<AttributeOrGroup>{DistWeightedCentralMoment20, DistWeightedCentralMoment02, DistWeightedAxisOrientation, DistWeightedEccentricity,
                                             DistWeightedCentralMoment20Exact, DistWeightedCentralMoment02Exact, DistWeightedAxisOrientationExact,
                                             DistWeightedEccentricityExact});
    const NodeId root = tree->root();
    const auto value = [&](Attribute attribute) { return values[names.linearIndex(root, attribute)]; };

    for (Attribute attribute : {DistWeightedCentralMoment20, DistWeightedCentralMoment20Exact}) {
        requireNear(value(attribute), 2.0, Tolerance, "weighted mu20 uses the column axis");
    }
    for (Attribute attribute : {DistWeightedCentralMoment02, DistWeightedCentralMoment02Exact}) {
        requireNear(value(attribute), 0.0, Tolerance, "weighted mu02 uses the row axis");
    }
    for (Attribute attribute : {DistWeightedAxisOrientation, DistWeightedAxisOrientationExact}) {
        requireNear(value(attribute), 0.0, Tolerance, "weighted orientation uses the column axis");
    }
    for (Attribute attribute : {DistWeightedEccentricity, DistWeightedEccentricityExact}) {
        requireNear(value(attribute), 1.0e6, Tolerance, "line-degenerate weighted eccentricity is finite");
    }
}

long double lowerQuantile(std::vector<long double> values, long double probability) {
    std::sort(values.begin(), values.end());
    const std::size_t rank = probability <= 0.0L ? 0 : static_cast<std::size_t>(std::ceil(probability * static_cast<long double>(values.size()))) - 1;
    return values[std::min(rank, values.size() - 1)];
}

void requireExactStudyMatchesIndependentOracle(const MorphologicalTree& tree, const std::string& label) {
    const std::vector<AttributeOrGroup> request{MaxDistExact,
                                                MaxSquaredDistExact,
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
    const auto [names, values] = AttributeComputation::computeTopologyAttributes<double>(tree, request);

    for (NodeId node : tree.aliveNodeIds()) {
        const auto oracle = maxdist_oracle::exactNodeSquaredDistanceTransform(tree, node);
        std::vector<long double> distances;
        distances.reserve(oracle.squaredDistances.size());
        std::map<std::int64_t, std::uint64_t> histogram;
        long double distanceSum = 0.0L;
        long double squaredDistanceSum = 0.0L;
        long double rowSum = 0.0L;
        long double columnSum = 0.0L;
        long double weightedRowSum = 0.0L;
        long double weightedColumnSum = 0.0L;
        long double weightedRowSquaredSum = 0.0L;
        long double weightedColumnSquaredSum = 0.0L;
        long double weightedRowColumnSum = 0.0L;
        for (std::size_t index = 0; index < oracle.supportPixels.size(); ++index) {
            const PixelId pixel = oracle.supportPixels[index];
            const std::int64_t squaredDistance = oracle.squaredDistances[index];
            const long double distance = std::sqrt(static_cast<long double>(squaredDistance));
            const long double row = static_cast<long double>(pixel / tree.numColumns());
            const long double column = static_cast<long double>(pixel % tree.numColumns());
            distances.push_back(distance);
            ++histogram[squaredDistance];
            distanceSum += distance;
            squaredDistanceSum += squaredDistance;
            rowSum += row;
            columnSum += column;
            weightedRowSum += distance * row;
            weightedColumnSum += distance * column;
            weightedRowSquaredSum += distance * row * row;
            weightedColumnSquaredSum += distance * column * column;
            weightedRowColumnSum += distance * row * column;
        }

        const long double count = static_cast<long double>(distances.size());
        const long double mean = distanceSum / count;
        const long double variance = std::max(0.0L, squaredDistanceSum / count - mean * mean);
        std::uint64_t modeCount = 0;
        long double mode = 0.0L;
        long double entropy = 0.0L;
        std::uint64_t positiveArea = 0;
        for (const auto& [squaredDistance, binCount] : histogram) {
            if (binCount > modeCount) {
                modeCount = binCount;
                mode = std::sqrt(static_cast<long double>(squaredDistance));
            }
            if (squaredDistance > 0) {
                positiveArea += binCount;
            }
            const long double probability = static_cast<long double>(binCount) / count;
            entropy -= probability * std::log2(probability);
        }

        const long double weightedCentroidRow = distanceSum > 0.0L ? weightedRowSum / distanceSum : rowSum / count;
        const long double weightedCentroidColumn = distanceSum > 0.0L ? weightedColumnSum / distanceSum : columnSum / count;
        const long double mu20 = distanceSum > 0.0L ? std::max(0.0L, weightedColumnSquaredSum - weightedColumnSum * weightedColumnSum / distanceSum) : 0.0L;
        const long double mu02 = distanceSum > 0.0L ? std::max(0.0L, weightedRowSquaredSum - weightedRowSum * weightedRowSum / distanceSum) : 0.0L;
        const long double mu11 = distanceSum > 0.0L ? weightedRowColumnSum - weightedRowSum * weightedColumnSum / distanceSum : 0.0L;
        const long double anisotropy = std::hypot(mu20 - mu02, 2.0L * mu11);
        const long double scale = std::max({1.0L, std::abs(mu20), std::abs(mu02)});
        const long double orientation =
            anisotropy <= 64.0L * std::numeric_limits<long double>::epsilon() * scale
                ? 0.0L
                : std::fmod(std::abs(0.5L * std::atan2(2.0L * mu11, mu20 - mu02) * 180.0L / std::numbers::pi_v<long double>), 360.0L);
        const long double trace = mu20 + mu02;
        const long double discriminant = std::sqrt(std::max(0.0L, (mu20 - mu02) * (mu20 - mu02) + 4.0L * mu11 * mu11));
        const long double major = 0.5L * (trace + discriminant);
        const long double minor = std::max(0.0L, 0.5L * (trace - discriminant));
        const long double eccentricity = major <= std::numeric_limits<long double>::epsilon()
                                             ? 1.0L
                                             : (minor <= std::numeric_limits<long double>::epsilon() ? 1.0e6L : std::min(1.0e6L, major / minor));
        const auto actual = [&](Attribute attribute) { return values[names.linearIndex(node, attribute)]; };
        const std::string nodeLabel = label + " node " + std::to_string(node) + " ";
        const std::int64_t maximumSquaredDistance = histogram.rbegin()->first;
        requireNear(actual(MaxDistExact), std::sqrt(static_cast<double>(maximumSquaredDistance)), Tolerance, nodeLabel + "maximum distance");
        requireEqual(actual(MaxSquaredDistExact), static_cast<double>(maximumSquaredDistance), nodeLabel + "maximum squared distance");
        requireNear(actual(DistSumExact), static_cast<double>(distanceSum), Tolerance, nodeLabel + "distance sum");
        requireNear(actual(DistMeanExact), static_cast<double>(mean), Tolerance, nodeLabel + "distance mean");
        requireNear(actual(DistVarianceExact), static_cast<double>(variance), Tolerance, nodeLabel + "distance variance");
        requireNear(actual(DistMedianExact), static_cast<double>(lowerQuantile(distances, 0.5L)), Tolerance, nodeLabel + "median");
        requireNear(actual(DistModeExact), static_cast<double>(mode), Tolerance, nodeLabel + "mode");
        requireNear(actual(DistQ25Exact), static_cast<double>(lowerQuantile(distances, 0.25L)), Tolerance, nodeLabel + "q25");
        requireNear(actual(DistQ75Exact), static_cast<double>(lowerQuantile(distances, 0.75L)), Tolerance, nodeLabel + "q75");
        requireNear(actual(DistQ90Exact), static_cast<double>(lowerQuantile(distances, 0.9L)), Tolerance, nodeLabel + "q90");
        requireNear(actual(DistEntropyExact), static_cast<double>(entropy), Tolerance, nodeLabel + "entropy");
        requireEqual(actual(DistPositiveAreaExact), static_cast<double>(positiveArea), nodeLabel + "positive area");
        requireEqual(actual(DistLevelCountExact), static_cast<double>(histogram.size()), nodeLabel + "level count");
        requireNear(actual(DistWeightedCentroidRowExact), static_cast<double>(weightedCentroidRow), Tolerance, nodeLabel + "weighted centroid row");
        requireNear(actual(DistWeightedCentroidColumnExact), static_cast<double>(weightedCentroidColumn), Tolerance, nodeLabel + "weighted centroid column");
        requireNear(actual(DistWeightedCentralMoment20Exact), static_cast<double>(mu20), Tolerance, nodeLabel + "weighted mu20");
        requireNear(actual(DistWeightedCentralMoment02Exact), static_cast<double>(mu02), Tolerance, nodeLabel + "weighted mu02");
        requireNear(actual(DistWeightedCentralMoment11Exact), static_cast<double>(mu11), Tolerance, nodeLabel + "weighted mu11");
        requireNear(actual(DistWeightedAxisOrientationExact), static_cast<double>(orientation), Tolerance, nodeLabel + "orientation");
        requireNear(actual(DistWeightedEccentricityExact), static_cast<double>(eccentricity), Tolerance, nodeLabel + "eccentricity");
    }
}

void requireApproximateStudyMatchesInternalField(const MorphologicalTree& tree, const std::string& label) {
    const std::vector<AttributeOrGroup> request{MaxDist,
                                                MaxSquaredDist,
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
                                                DistWeightedCentralMoment11};
    const auto [names, values] = AttributeComputation::computeTopologyAttributes<double>(tree, request);

    using ApproximateTransform = attributes::computers::detail::distance_transform_approx::MorphologicalTreeApproximateDistanceTransform;
    ApproximateTransform::forEachNodeField(tree, [&](NodeId node, std::span<const PixelId> supportPixels, const auto& squaredDistanceAt) {
        std::vector<long double> distances;
        std::map<std::int64_t, std::uint64_t> histogram;
        long double distanceSum = 0.0L;
        long double squaredDistanceSum = 0.0L;
        long double weightedRowSum = 0.0L;
        long double weightedColumnSum = 0.0L;
        long double weightedRowSquaredSum = 0.0L;
        long double weightedColumnSquaredSum = 0.0L;
        long double weightedRowColumnSum = 0.0L;
        for (PixelId pixel : supportPixels) {
            const std::int64_t squaredDistance = squaredDistanceAt(pixel);
            const long double distance = std::sqrt(static_cast<long double>(squaredDistance));
            const long double row = static_cast<long double>(pixel / tree.numColumns());
            const long double column = static_cast<long double>(pixel % tree.numColumns());
            distances.push_back(distance);
            ++histogram[squaredDistance];
            distanceSum += distance;
            squaredDistanceSum += squaredDistance;
            weightedRowSum += distance * row;
            weightedColumnSum += distance * column;
            weightedRowSquaredSum += distance * row * row;
            weightedColumnSquaredSum += distance * column * column;
            weightedRowColumnSum += distance * row * column;
        }
        const long double count = static_cast<long double>(supportPixels.size());
        const long double mean = distanceSum / count;
        const long double centroidRow = distanceSum > 0.0L ? weightedRowSum / distanceSum : 0.0L;
        const long double centroidColumn = distanceSum > 0.0L ? weightedColumnSum / distanceSum : 0.0L;
        const long double mu20 = distanceSum > 0.0L ? std::max(0.0L, weightedColumnSquaredSum - weightedColumnSum * weightedColumnSum / distanceSum) : 0.0L;
        const long double mu02 = distanceSum > 0.0L ? std::max(0.0L, weightedRowSquaredSum - weightedRowSum * weightedRowSum / distanceSum) : 0.0L;
        const long double mu11 = distanceSum > 0.0L ? weightedRowColumnSum - weightedRowSum * weightedColumnSum / distanceSum : 0.0L;
        std::uint64_t modeCount = 0;
        long double mode = 0.0L;
        long double entropy = 0.0L;
        std::uint64_t positiveArea = 0;
        for (const auto& [squaredDistance, binCount] : histogram) {
            if (binCount > modeCount) {
                modeCount = binCount;
                mode = std::sqrt(static_cast<long double>(squaredDistance));
            }
            positiveArea += squaredDistance > 0 ? binCount : 0;
            const long double probability = static_cast<long double>(binCount) / count;
            entropy -= probability * std::log2(probability);
        }
        const auto actual = [&](Attribute attribute) { return values[names.linearIndex(node, attribute)]; };
        const std::string nodeLabel = label + " node " + std::to_string(node) + " ";
        const std::int64_t maximumSquaredDistance = histogram.rbegin()->first;
        requireNear(actual(MaxDist), std::sqrt(static_cast<double>(maximumSquaredDistance)), Tolerance, nodeLabel + "maximum distance");
        requireEqual(actual(MaxSquaredDist), static_cast<double>(maximumSquaredDistance), nodeLabel + "maximum squared distance");
        requireNear(actual(DistSum), static_cast<double>(distanceSum), Tolerance, nodeLabel + "distance sum");
        requireNear(actual(DistMean), static_cast<double>(mean), Tolerance, nodeLabel + "distance mean");
        requireNear(actual(DistVariance), static_cast<double>(std::max(0.0L, squaredDistanceSum / count - mean * mean)), Tolerance,
                    nodeLabel + "distance variance");
        requireNear(actual(DistMedian), static_cast<double>(lowerQuantile(distances, 0.5L)), Tolerance, nodeLabel + "median");
        requireNear(actual(DistMode), static_cast<double>(mode), Tolerance, nodeLabel + "mode");
        requireNear(actual(DistQ25), static_cast<double>(lowerQuantile(distances, 0.25L)), Tolerance, nodeLabel + "q25");
        requireNear(actual(DistQ75), static_cast<double>(lowerQuantile(distances, 0.75L)), Tolerance, nodeLabel + "q75");
        requireNear(actual(DistQ90), static_cast<double>(lowerQuantile(distances, 0.9L)), Tolerance, nodeLabel + "q90");
        requireNear(actual(DistEntropy), static_cast<double>(entropy), Tolerance, nodeLabel + "entropy");
        requireEqual(actual(DistPositiveArea), static_cast<double>(positiveArea), nodeLabel + "positive area");
        requireEqual(actual(DistLevelCount), static_cast<double>(histogram.size()), nodeLabel + "level count");
        if (distanceSum > 0.0L) {
            requireNear(actual(DistWeightedCentroidRow), static_cast<double>(centroidRow), Tolerance, nodeLabel + "weighted centroid row");
            requireNear(actual(DistWeightedCentroidColumn), static_cast<double>(centroidColumn), Tolerance, nodeLabel + "weighted centroid column");
        }
        requireNear(actual(DistWeightedCentralMoment20), static_cast<double>(mu20), Tolerance, nodeLabel + "weighted mu20");
        requireNear(actual(DistWeightedCentralMoment02), static_cast<double>(mu02), Tolerance, nodeLabel + "weighted mu02");
        requireNear(actual(DistWeightedCentralMoment11), static_cast<double>(mu11), Tolerance, nodeLabel + "weighted mu11");
    });
}

} // namespace

int main() {
    auto image = ImageUInt8::create(5, 5);
    image->fill(7);
    const auto tree = makeComponentTree(image, true);

    requireRootStudy(*tree, true);
    requireRootStudy(*tree, false);
    requireWeightedMomentAxisConvention();

    requireExactStudyMatchesIndependentOracle(*tree, "constant square");
    const auto fixtureTree = makeComponentTree(makeComponentTreeFixture(), true);
    requireExactStudyMatchesIndependentOracle(*fixtureTree, "component fixture");
    requireApproximateStudyMatchesInternalField(*fixtureTree, "component fixture approximate");

    auto concaveImage = ImageUInt8::create(7, 7);
    concaveImage->fill(0);
    for (PixelId pixel : std::vector<PixelId>{8, 9, 10, 15, 22, 23, 24, 25, 26, 17, 19, 31, 32, 33}) {
        (*concaveImage)[pixel] = 9;
    }
    const auto concaveTree = makeComponentTree(concaveImage, true);
    requireExactStudyMatchesIndependentOracle(*concaveTree, "concave support");
    requireApproximateStudyMatchesInternalField(*concaveTree, "concave support approximate");

    return 0;
}
