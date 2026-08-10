#include "ApiBenchmark.hpp"

#include "mmcfilters/filters/AttributeFilters.hpp"
#include "mmcfilters/filters/DepthStableRegionComputer.hpp"
#include "mmcfilters/filters/ExtinctionValues.hpp"
#include "mmcfilters/filters/MSERComputer.hpp"
#include "mmcfilters/filters/UltimateAttributeOpening.hpp"

#include <cstdint>
#include <vector>

namespace mmcfilters::benchmarks::api {
namespace {

struct UaoResult {
    ImageUInt8Ptr maxContrast;
    ImageInt32Ptr associated;
};

[[nodiscard]] std::uint64_t boolMaskChecksum(const std::vector<bool>& mask) noexcept {
    Fnv1a64 hash;
    hash.append(mask.size());
    for (bool value : mask) {
        hash.append(static_cast<std::uint8_t>(value));
    }
    return hash.value();
}

[[nodiscard]] std::uint64_t byteMaskChecksum(const std::vector<std::uint8_t>& mask) noexcept {
    Fnv1a64 hash;
    hash.appendVector(mask);
    return hash.value();
}

[[nodiscard]] std::uint64_t uaoChecksum(const UaoResult& result) {
    Fnv1a64 hash;
    appendImage(hash, result.maxContrast);
    appendImage(hash, result.associated);
    return hash.value();
}

} // namespace

void addFilterScenarios(Context& context, std::vector<ScenarioResult>& results) {
    results.push_back(benchmarkScenario(
        "filters", "direct_criterion", TimingScope::EstablishedInput, context.options.repetitions, context.maxTreeMetrics,
        [&] {
            AttributeFilters<std::uint8_t> filters(context.maxTree);
            return filters.filteringByDirectRule(context.keepCriterion);
        },
        [](const ImageUInt8Ptr& result) { return imageChecksum(result); }));

    results.push_back(benchmarkScenario(
        "filters", "pruning_max_attribute", TimingScope::EstablishedInput, context.options.repetitions, context.maxTreeMetrics,
        [&] {
            AttributeFilters<std::uint8_t> filters(context.maxTree);
            return filters.filteringByPruningMax(context.area.second.data(), context.areaThreshold);
        },
        [](const ImageUInt8Ptr& result) { return imageChecksum(result); }));

    results.push_back(benchmarkScenario(
        "filters", "viterbi_attribute", TimingScope::EstablishedInput, context.options.repetitions, context.maxTreeMetrics,
        [&] {
            AttributeFilters<std::uint8_t> filters(context.maxTree);
            return filters.filteringByViterbiRule(context.area.second.data(), context.areaThreshold);
        },
        [](const ImageUInt8Ptr& result) { return imageChecksum(result); }));

    if (!atLeast(context.options.profile, Profile::Core)) {
        return;
    }

    results.push_back(benchmarkScenario(
        "filters", "subtractive_criterion", TimingScope::EstablishedInput, context.options.repetitions, context.maxTreeMetrics,
        [&] {
            AttributeFilters<std::uint8_t> filters(context.maxTree);
            return filters.filteringBySubtractiveRule(context.keepCriterion);
        },
        [](const ImageUInt8Ptr& result) { return imageChecksum(result); }));

    results.push_back(benchmarkScenario(
        "filters", "subtractive_score", TimingScope::EstablishedInput, context.options.repetitions, context.maxTreeMetrics,
        [&] {
            AttributeFilters<std::uint8_t> filters(context.maxTree);
            return filters.filteringBySubtractiveScoreRule(context.scores);
        },
        [](const ImageFloatPtr& result) { return imageChecksum(result); }));

    results.push_back(benchmarkScenario(
        "filters", "pruning_min_criterion", TimingScope::EstablishedInput, context.options.repetitions, context.maxTreeMetrics,
        [&] {
            AttributeFilters<std::uint8_t> filters(context.maxTree);
            return filters.filteringByPruningMin(context.keepCriterion);
        },
        [](const ImageUInt8Ptr& result) { return imageChecksum(result); }));

    results.push_back(benchmarkScenario(
        "filters", "pruning_max_criterion", TimingScope::EstablishedInput, context.options.repetitions, context.maxTreeMetrics,
        [&] {
            AttributeFilters<std::uint8_t> filters(context.maxTree);
            return filters.filteringByPruningMax(context.keepCriterion);
        },
        [](const ImageUInt8Ptr& result) { return imageChecksum(result); }));

    results.push_back(benchmarkScenario(
        "filters", "pruning_min_attribute", TimingScope::EstablishedInput, context.options.repetitions, context.maxTreeMetrics,
        [&] {
            AttributeFilters<std::uint8_t> filters(context.maxTree);
            return filters.filteringByPruningMin(context.area.second.data(), context.areaThreshold);
        },
        [](const ImageUInt8Ptr& result) { return imageChecksum(result); }));

    results.push_back(benchmarkScenario(
        "filters", "adaptive_mser_criterion", TimingScope::EstablishedInput, context.options.repetitions, context.maxTreeMetrics,
        [&] {
            AttributeFilters<std::uint8_t> filters(context.maxTree);
            return filters.getAdaptiveCriterion(context.keepCriterion, AltitudeDiff<std::uint8_t>{4});
        },
        [](const std::vector<bool>& result) { return boolMaskChecksum(result); }));

    results.push_back(benchmarkScenario(
        "filters", "adaptive_depth_criterion", TimingScope::EstablishedInput, context.options.repetitions, context.maxTreeMetrics,
        [&] {
            AttributeFilters<std::uint8_t> filters(context.maxTree);
            return filters.getAdaptiveCriterionByDepth(context.keepCriterion, 4);
        },
        [](const std::vector<bool>& result) { return boolMaskChecksum(result); }));

    results.push_back(benchmarkScenario(
        "filters", "ultimate_attribute_opening", TimingScope::EstablishedInput, context.options.repetitions, context.maxTreeMetrics,
        [&] {
            UltimateAttributeOpening<std::uint8_t, double> uao(context.maxTree, context.area.second.data());
            uao.execute(context.areaThreshold);
            return UaoResult{uao.getMaxContrastImage(), uao.getAssociatedImage()};
        },
        [](const UaoResult& result) { return uaoChecksum(result); }));

    results.push_back(benchmarkScenario(
        "filters", "extinction_top_16", TimingScope::EstablishedInput, context.options.repetitions, context.maxTreeMetrics,
        [&] {
            ExtinctionValues<std::uint8_t, double> extinction(context.maxTree, context.area.second.data());
            return extinction.filtering(ExtinctionSelectionPolicy<double>::byTopK(16));
        },
        [](const ImageUInt8Ptr& result) { return imageChecksum(result); }));

    results.push_back(benchmarkScenario(
        "filters", "extinction_contour_rank", TimingScope::EstablishedInput, context.options.repetitions, context.maxTreeMetrics,
        [&] {
            ExtinctionValues<std::uint8_t, double> extinction(context.maxTree, context.area.second.data());
            return extinction.contourMap(ExtinctionSelectionPolicy<double>::byTopK(16), ExtinctionContourScorePolicy::RankScore);
        },
        [](const ImagePtr<double>& result) { return imageChecksum(result); }));

    results.push_back(benchmarkScenario(
        "filters", "mser", TimingScope::EstablishedInput, context.options.repetitions, context.maxTreeMetrics,
        [&] {
            MSERComputer<std::uint8_t, double> mser(context.maxTree, context.area.second.data());
            return mser.computeMSER(AltitudeDiff<std::uint8_t>{4});
        },
        [](const std::vector<std::uint8_t>& result) { return byteMaskChecksum(result); }));

    results.push_back(benchmarkScenario(
        "filters", "depth_stability", TimingScope::EstablishedInput, context.options.repetitions, context.maxTreeMetrics,
        [&] {
            DepthStableRegionComputer<double> depth(context.maxTree.topology(), context.area.second.data());
            return depth.computeByDepth(4);
        },
        [](const std::vector<std::uint8_t>& result) { return byteMaskChecksum(result); }));

    if (!atLeast(context.options.profile, Profile::Publication)) {
        return;
    }

    WorkloadMetrics saliencyMetrics = context.maxTreeMetrics;
    saliencyMetrics.edges = countUndirectedEdges(context.adjacency);
    results.push_back(benchmarkScenario(
        "filters", "extinction_ranked_formal_saliency", TimingScope::EstablishedInput, context.options.repetitions, saliencyMetrics,
        [&] {
            ExtinctionValues<std::uint8_t, double> extinction(context.maxTree, context.area.second.data());
            return extinction.computeRankedFormalSaliencyEdgeMap(context.adjacency);
        },
        [](const EdgeSaliencyMap<int>& result) { return edgeMapChecksum(result); }));
}

} // namespace mmcfilters::benchmarks::api
