#include "ApiBenchmark.hpp"

#include "mmcfilters/trees/adjust/CasfComponentTrees.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

namespace mmcfilters::benchmarks::api {
namespace {

using adjust::CasfComponentTreesAttribute;

template <AltitudeValue T> struct OwnedCasfResult {
    std::unique_ptr<adjust::CasfComponentTrees<T>> casf;
    ImagePtr<T> image;
};

struct ThresholdPlan {
    std::vector<double> thresholds;
    bool hasPruningCandidates = false;
};

struct EditSnapshot {
    std::int64_t minNodes = 0;
    std::int64_t maxNodes = 0;
    std::int64_t completeCommits = 0;
    std::int64_t incrementalCommits = 0;
};

template <AltitudeValue T> [[nodiscard]] WorkloadMetrics casfMetrics(const adjust::CasfComponentTrees<T>& casf, std::int64_t edges) {
    const MorphologicalTree& minTopology = casf.minTree().topology();
    const MorphologicalTree& maxTopology = casf.maxTree().topology();
    return WorkloadMetrics{
        .pixels = static_cast<std::int64_t>(minTopology.getNumRowsOfGridDomain2D()) * minTopology.getNumColsOfGridDomain2D(),
        .properParts = minTopology.getNumTotalProperParts(),
        .primaryNodeSlots = minTopology.getNumInternalNodeSlots(),
        .primaryLiveNodes = minTopology.getNumNodes(),
        .secondaryNodeSlots = maxTopology.getNumInternalNodeSlots(),
        .secondaryLiveNodes = maxTopology.getNumNodes(),
        .edges = edges,
    };
}

template <AltitudeValue T> [[nodiscard]] EditSnapshot editSnapshot(const adjust::CasfComponentTrees<T>& casf) {
    const MorphologicalTree& minTopology = casf.minTree().topology();
    const MorphologicalTree& maxTopology = casf.maxTree().topology();
    const TreeEditValidationStatistics& minStatistics = minTopology.getEditValidationStatistics();
    const TreeEditValidationStatistics& maxStatistics = maxTopology.getEditValidationStatistics();
    return EditSnapshot{
        .minNodes = minTopology.getNumNodes(),
        .maxNodes = maxTopology.getNumNodes(),
        .completeCommits = static_cast<std::int64_t>(minStatistics.completeValidationCommits + maxStatistics.completeValidationCommits),
        .incrementalCommits = static_cast<std::int64_t>(minStatistics.incrementalValidationCommits + maxStatistics.incrementalValidationCommits),
    };
}

[[nodiscard]] OutcomeMetrics outcomeBetween(const EditSnapshot& before, const EditSnapshot& after, const ThresholdPlan& plan, std::size_t steps,
                                            bool requireIncrementalEdit) {
    OutcomeMetrics outcome{
        .steps = static_cast<std::int64_t>(steps),
        .primaryNodesRemoved = before.minNodes - after.minNodes,
        .secondaryNodesRemoved = before.maxNodes - after.maxNodes,
        .completeValidationCommits = after.completeCommits - before.completeCommits,
        .incrementalValidationCommits = after.incrementalCommits - before.incrementalCommits,
        .lightThreshold = plan.thresholds[0],
        .mediumThreshold = plan.thresholds[1],
        .heavyThreshold = plan.thresholds[2],
    };
    if (outcome.completeValidationCommits != 0) {
        throw std::runtime_error("CASF benchmark detected a complete validation commit in the incremental hot path.");
    }
    if (requireIncrementalEdit && plan.hasPruningCandidates && outcome.incrementalValidationCommits == 0) {
        throw std::runtime_error("CASF benchmark quantiles did not exercise an incremental edit.");
    }
    return outcome;
}

template <AltitudeValue T>
void appendNonRootAttributeValues(const WeightedMorphologicalTree<T>& tree, Attribute attribute, std::vector<double>& values) {
    const auto computed = AttributeComputation::computeSingleAttribute<double>(tree, attribute);
    for (NodeId node : tree.topology().getAliveNodeIds()) {
        if (!tree.topology().isRoot(node)) {
            values.push_back(computed.second[static_cast<std::size_t>(node)]);
        }
    }
}

template <AltitudeValue T>
[[nodiscard]] ThresholdPlan thresholdsFromQuantiles(const adjust::CasfComponentTrees<T>& casf, CasfComponentTreesAttribute attribute,
                                                    const std::vector<double>& quantiles) {
    const Attribute scalarAttribute = attribute == CasfComponentTreesAttribute::AREA ? AREA
                                      : attribute == CasfComponentTreesAttribute::BOUNDING_BOX_WIDTH
                                          ? BOX_WIDTH
                                      : attribute == CasfComponentTreesAttribute::BOUNDING_BOX_HEIGHT ? BOX_HEIGHT
                                                                                                      : DIAGONAL_LENGTH;
    std::vector<double> values;
    values.reserve(static_cast<std::size_t>(casf.minTree().topology().getNumNodes() + casf.maxTree().topology().getNumNodes()));
    appendNonRootAttributeValues(casf.minTree(), scalarAttribute, values);
    appendNonRootAttributeValues(casf.maxTree(), scalarAttribute, values);
    if (values.empty()) {
        return ThresholdPlan{.thresholds = std::vector<double>(quantiles.size(), 0.0), .hasPruningCandidates = false};
    }

    std::sort(values.begin(), values.end());
    std::vector<double> thresholds;
    thresholds.reserve(quantiles.size());
    for (std::size_t index = 0; index < quantiles.size(); ++index) {
        std::size_t selected = static_cast<std::size_t>(std::floor(quantiles[index] * static_cast<double>(values.size() - 1)));
        if (!thresholds.empty() && values[selected] <= thresholds.back()) {
            const auto greater = std::upper_bound(values.begin() + static_cast<std::ptrdiff_t>(selected), values.end(), thresholds.back());
            if (greater != values.end()) {
                selected = static_cast<std::size_t>(std::distance(values.begin(), greater));
            }
        }
        thresholds.push_back(values[selected]);
    }
    return ThresholdPlan{std::move(thresholds), true};
}

template <AltitudeValue T>
[[nodiscard]] OutcomeMetrics observeSequence(const ImagePtr<T>& image, CasfComponentTreesAttribute attribute, const ThresholdPlan& plan,
                                             const std::vector<double>& prefix = {}) {
    adjust::CasfComponentTrees<T> casf(image, attribute, 1.5);
    if (!prefix.empty()) {
        static_cast<void>(casf.filter(prefix));
    }
    const EditSnapshot before = editSnapshot(casf);
    static_cast<void>(casf.filter(plan.thresholds));
    return outcomeBetween(before, editSnapshot(casf), plan, plan.thresholds.size(), true);
}

template <AltitudeValue T>
[[nodiscard]] OutcomeMetrics observeStep(const ImagePtr<T>& image, CasfComponentTreesAttribute attribute, const ThresholdPlan& plan,
                                         std::size_t stepIndex) {
    adjust::CasfComponentTrees<T> casf(image, attribute, 1.5);
    const std::vector<double> prefix(plan.thresholds.begin(), plan.thresholds.begin() + static_cast<std::ptrdiff_t>(stepIndex));
    if (!prefix.empty()) {
        static_cast<void>(casf.filter(prefix));
    }
    const EditSnapshot before = editSnapshot(casf);
    static_cast<void>(casf.filter({plan.thresholds[stepIndex]}));
    return outcomeBetween(before, editSnapshot(casf), plan, 1, false);
}

template <AltitudeValue T> void appendCasfState(Fnv1a64& hash, const adjust::CasfComponentTrees<T>& casf, bool requireIncrementalEdit) {
    hash.append(weightedTreeChecksum(casf.minTree()));
    hash.append(weightedTreeChecksum(casf.maxTree()));
    const EditSnapshot snapshot = editSnapshot(casf);
    if (snapshot.completeCommits != 0) {
        throw std::runtime_error("CASF benchmark detected a complete validation commit in the incremental hot path.");
    }
    if (requireIncrementalEdit && snapshot.incrementalCommits == 0) {
        throw std::runtime_error("CASF benchmark threshold sequence did not exercise an incremental edit.");
    }
    hash.append(snapshot.completeCommits);
    hash.append(snapshot.incrementalCommits);
}

template <AltitudeValue T> [[nodiscard]] std::uint64_t casfStateChecksum(const adjust::CasfComponentTrees<T>& casf, bool requireIncrementalEdit) {
    Fnv1a64 hash;
    appendCasfState(hash, casf, requireIncrementalEdit);
    return hash.value();
}

template <AltitudeValue T>
[[nodiscard]] std::uint64_t casfImageChecksum(const adjust::CasfComponentTrees<T>& casf, const ImagePtr<T>& image, bool requireIncrementalEdit) {
    Fnv1a64 hash;
    appendImage(hash, image);
    appendCasfState(hash, casf, requireIncrementalEdit);
    return hash.value();
}

template <AltitudeValue T>
void addTypedCasfPublicationScenarios(Context& context, std::vector<ScenarioResult>& results, const ImagePtr<T>& image, std::string typeName,
                                      CasfComponentTreesAttribute attribute) {
    adjust::CasfComponentTrees<T> baseline(image, attribute, 1.5);
    const ThresholdPlan plan = thresholdsFromQuantiles(baseline, attribute, context.options.casfQuantiles);
    const WorkloadMetrics metrics = casfMetrics(baseline, context.maxTreeMetrics.edges);
    results.push_back(benchmarkScenario(
        "casf", "construction_" + typeName, TimingScope::EndToEnd, context.options.repetitions, metrics,
        [&] { return std::make_unique<adjust::CasfComponentTrees<T>>(image, attribute, 1.5); },
        [](const std::unique_ptr<adjust::CasfComponentTrees<T>>& result) { return casfStateChecksum(*result, false); }));

    ScenarioResult sequence = benchmarkPreparedScenario(
        "casf", "incremental_sequence_" + typeName, context.options.repetitions, metrics,
        [&] { return std::make_unique<adjust::CasfComponentTrees<T>>(image, attribute, 1.5); },
        [&](std::unique_ptr<adjust::CasfComponentTrees<T>>& state) { return state->filter(plan.thresholds); },
        [&](const std::unique_ptr<adjust::CasfComponentTrees<T>>& state, const ImagePtr<T>& result) {
            return casfImageChecksum(*state, result, plan.hasPruningCandidates);
        });
    sequence.outcome = observeSequence(image, attribute, plan);
    results.push_back(std::move(sequence));
}

void addAreaStepScenario(Context& context, std::vector<ScenarioResult>& results, const WorkloadMetrics& metrics, const ThresholdPlan& plan,
                         std::string name, std::size_t stepIndex) {
    const std::vector<double> prefix(plan.thresholds.begin(), plan.thresholds.begin() + static_cast<std::ptrdiff_t>(stepIndex));
    const std::vector<double> threshold{plan.thresholds[stepIndex]};
    ScenarioResult scenario = benchmarkPreparedScenario(
        "casf", std::move(name), context.options.repetitions, metrics,
        [&] {
            auto casf = std::make_unique<adjust::CasfComponentTrees<std::uint8_t>>(context.imageUInt8, CasfComponentTreesAttribute::AREA, 1.5);
            if (!prefix.empty()) {
                static_cast<void>(casf->filter(prefix));
            }
            return casf;
        },
        [&](std::unique_ptr<adjust::CasfComponentTrees<std::uint8_t>>& state) { return state->filter(threshold); },
        [&](const std::unique_ptr<adjust::CasfComponentTrees<std::uint8_t>>& state, const ImageUInt8Ptr& result) {
            return casfImageChecksum(*state, result, false);
        });
    scenario.outcome = observeStep(context.imageUInt8, CasfComponentTreesAttribute::AREA, plan, stepIndex);
    results.push_back(std::move(scenario));
}

} // namespace

void addCasfScenarios(Context& context, std::vector<ScenarioResult>& results) {
    adjust::CasfComponentTrees<std::uint8_t> baseline(context.imageUInt8, CasfComponentTreesAttribute::AREA, 1.5);
    const ThresholdPlan areaPlan = thresholdsFromQuantiles(baseline, CasfComponentTreesAttribute::AREA, context.options.casfQuantiles);
    const WorkloadMetrics metrics = casfMetrics(baseline, context.maxTreeMetrics.edges);

    results.push_back(benchmarkScenario(
        "casf", "construction_area", TimingScope::EndToEnd, context.options.repetitions, metrics,
        [&] { return std::make_unique<adjust::CasfComponentTrees<std::uint8_t>>(context.imageUInt8, CasfComponentTreesAttribute::AREA, 1.5); },
        [](const std::unique_ptr<adjust::CasfComponentTrees<std::uint8_t>>& result) { return casfStateChecksum(*result, false); }));

    ScenarioResult areaSequence = benchmarkPreparedScenario(
        "casf", "incremental_area_sequence", context.options.repetitions, metrics,
        [&] { return std::make_unique<adjust::CasfComponentTrees<std::uint8_t>>(context.imageUInt8, CasfComponentTreesAttribute::AREA, 1.5); },
        [&](std::unique_ptr<adjust::CasfComponentTrees<std::uint8_t>>& state) { return state->filter(areaPlan.thresholds); },
        [&](const std::unique_ptr<adjust::CasfComponentTrees<std::uint8_t>>& state, const ImageUInt8Ptr& result) {
            return casfImageChecksum(*state, result, areaPlan.hasPruningCandidates);
        });
    areaSequence.outcome = observeSequence(context.imageUInt8, CasfComponentTreesAttribute::AREA, areaPlan);
    results.push_back(std::move(areaSequence));

    if (!atLeast(context.options.profile, Profile::Core)) {
        return;
    }

    ScenarioResult pipeline = benchmarkScenario(
        "casf", "pipeline_area_sequence", TimingScope::EndToEnd, context.options.repetitions, metrics,
        [&] {
            auto casf = std::make_unique<adjust::CasfComponentTrees<std::uint8_t>>(context.imageUInt8, CasfComponentTreesAttribute::AREA, 1.5);
            auto image = casf->filter(areaPlan.thresholds);
            return OwnedCasfResult<std::uint8_t>{std::move(casf), std::move(image)};
        },
        [&](const OwnedCasfResult<std::uint8_t>& result) {
            return casfImageChecksum(*result.casf, result.image, areaPlan.hasPruningCandidates);
        });
    pipeline.outcome = observeSequence(context.imageUInt8, CasfComponentTreesAttribute::AREA, areaPlan);
    results.push_back(std::move(pipeline));

    addAreaStepScenario(context, results, metrics, areaPlan, "incremental_area_light_step", 0);
    addAreaStepScenario(context, results, metrics, areaPlan, "incremental_area_medium_step_after_light", 1);
    addAreaStepScenario(context, results, metrics, areaPlan, "incremental_area_heavy_step_after_medium", 2);

    adjust::CasfComponentTrees<std::uint8_t> boundingBoxBaseline(context.imageUInt8, CasfComponentTreesAttribute::BOUNDING_BOX_DIAGONAL, 1.5);
    const ThresholdPlan boundingBoxPlan =
        thresholdsFromQuantiles(boundingBoxBaseline, CasfComponentTreesAttribute::BOUNDING_BOX_DIAGONAL, context.options.casfQuantiles);
    ScenarioResult boundingBox = benchmarkPreparedScenario(
        "casf", "incremental_bounding_box_diagonal_sequence", context.options.repetitions, metrics,
        [&] {
            return std::make_unique<adjust::CasfComponentTrees<std::uint8_t>>(context.imageUInt8, CasfComponentTreesAttribute::BOUNDING_BOX_DIAGONAL, 1.5);
        },
        [&](std::unique_ptr<adjust::CasfComponentTrees<std::uint8_t>>& state) { return state->filter(boundingBoxPlan.thresholds); },
        [&](const std::unique_ptr<adjust::CasfComponentTrees<std::uint8_t>>& state, const ImageUInt8Ptr& result) {
            return casfImageChecksum(*state, result, boundingBoxPlan.hasPruningCandidates);
        });
    boundingBox.outcome = observeSequence(context.imageUInt8, CasfComponentTreesAttribute::BOUNDING_BOX_DIAGONAL, boundingBoxPlan);
    results.push_back(std::move(boundingBox));

    if (!atLeast(context.options.profile, Profile::Publication)) {
        return;
    }

    addTypedCasfPublicationScenarios<std::int32_t>(context, results, context.imageInt32, "int32_area", CasfComponentTreesAttribute::AREA);
    addTypedCasfPublicationScenarios<float>(context, results, context.imageFloat, "float_bounding_box_diagonal",
                                            CasfComponentTreesAttribute::BOUNDING_BOX_DIAGONAL);
}

} // namespace mmcfilters::benchmarks::api
