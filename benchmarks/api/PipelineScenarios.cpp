#include "ApiBenchmark.hpp"

#include "mmcfilters/filters/AttributeFilters.hpp"

#include <algorithm>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace mmcfilters::benchmarks::api {
namespace {

template <AltitudeValue T> struct PipelineResult {
    WeightedMorphologicalTree<T> tree;
    ComputedAttributeData<double> attributes;
    ImagePtr<T> image;
};

template <AltitudeValue T> [[nodiscard]] std::uint64_t pipelineChecksum(const PipelineResult<T>& result) {
    Fnv1a64 hash;
    hash.append(weightedTreeChecksum(result.tree));
    hash.append(semanticAttributeChecksum(result.attributes));
    appendImage(hash, result.image);
    return hash.value();
}

template <AltitudeValue T>
[[nodiscard]] std::vector<bool> criterionFromAttribute(const WeightedMorphologicalTree<T>& tree,
                                                       const ComputedAttributeData<double>& attributes, Attribute attribute,
                                                       double threshold) {
    const int stride = attributes.first.NUM_ATTRIBUTES;
    const int offset = attributes.first.getIndex(attribute);
    std::vector<bool> criterion(static_cast<std::size_t>(tree.topology().getNumInternalNodeSlots()), false);
    for (NodeId node : tree.topology().getAliveNodeIds()) {
        const std::size_t row = static_cast<std::size_t>(node);
        criterion[row] = attributes.second[row * static_cast<std::size_t>(stride) + static_cast<std::size_t>(offset)] >= threshold;
    }
    criterion[static_cast<std::size_t>(tree.topology().getRoot())] = true;
    return criterion;
}

template <AltitudeValue T>
[[nodiscard]] double relativeThreshold(const WeightedMorphologicalTree<T>& tree, const ComputedAttributeData<double>& attributes,
                                       Attribute attribute, double fraction) {
    const int stride = attributes.first.NUM_ATTRIBUTES;
    const int offset = attributes.first.getIndex(attribute);
    double maximum = 0.0;
    for (NodeId node : tree.topology().getAliveNodeIds()) {
        const std::size_t row = static_cast<std::size_t>(node);
        maximum = std::max(maximum, attributes.second[row * static_cast<std::size_t>(stride) + static_cast<std::size_t>(offset)]);
    }
    return maximum * fraction;
}

template <class TreeFactory>
[[nodiscard]] PipelineResult<std::uint8_t> directPipeline(TreeFactory&& createTree, std::vector<AttributeOrGroup> request,
                                                         Attribute criterionAttribute, double thresholdFraction) {
    auto tree = createTree();
    auto attributes = AttributeComputation::computeAttributes<double>(tree, request);
    const double threshold = relativeThreshold(tree, attributes, criterionAttribute, thresholdFraction);
    auto criterion = criterionFromAttribute(tree, attributes, criterionAttribute, threshold);
    AttributeFilters<std::uint8_t> filters(tree);
    auto image = filters.filteringByDirectRule(criterion);
    return PipelineResult<std::uint8_t>{std::move(tree), std::move(attributes), std::move(image)};
}

} // namespace

void addPipelineScenarios(Context& context, std::vector<ScenarioResult>& results) {
    results.push_back(benchmarkScenario(
        "pipelines", "max_tree_area_direct", TimingScope::EndToEnd, context.options.repetitions, context.maxTreeMetrics,
        [&] {
            return directPipeline(
                [&] { return MorphologicalTreeFactory::createMaxTree(context.imageUInt8, context.adjacency); }, {AttributeOrGroup{AREA}}, AREA, 1.0 / 16.0);
        },
        [](const PipelineResult<std::uint8_t>& result) { return pipelineChecksum(result); }));

    if (!atLeast(context.options.profile, Profile::Core)) {
        return;
    }

    results.push_back(benchmarkScenario(
        "pipelines", "max_tree_max_dist_direct", TimingScope::EndToEnd, context.options.repetitions, context.maxTreeMetrics,
        [&] {
            return directPipeline(
                [&] { return MorphologicalTreeFactory::createMaxTree(context.imageUInt8, context.adjacency); }, {AttributeOrGroup{MAX_DIST}}, MAX_DIST, 0.25);
        },
        [](const PipelineResult<std::uint8_t>& result) { return pipelineChecksum(result); }));

    results.push_back(benchmarkScenario(
        "pipelines", "max_tree_shape_group_direct_area", TimingScope::EndToEnd, context.options.repetitions, context.maxTreeMetrics,
        [&] {
            return directPipeline([&] { return MorphologicalTreeFactory::createMaxTree(context.imageUInt8, context.adjacency); },
                                  {AttributeOrGroup{AttributeGroup::SHAPE}}, AREA, 1.0 / 16.0);
        },
        [](const PipelineResult<std::uint8_t>& result) { return pipelineChecksum(result); }));

    results.push_back(benchmarkScenario(
        "pipelines", "min_tree_area_pruning_min", TimingScope::EndToEnd, context.options.repetitions, metricsOf(context.minTree),
        [&] {
            auto tree = MorphologicalTreeFactory::createMinTree(context.imageUInt8, context.adjacency);
            auto attributes = AttributeComputation::computeSingleAttribute<double>(tree, AREA);
            AttributeFilters<std::uint8_t> filters(tree);
            auto image = filters.filteringByPruningMin(attributes.second.data(), context.areaThreshold);
            return PipelineResult<std::uint8_t>{std::move(tree), std::move(attributes), std::move(image)};
        },
        [](const PipelineResult<std::uint8_t>& result) { return pipelineChecksum(result); }));

    if (!atLeast(context.options.profile, Profile::Publication)) {
        return;
    }

    auto treeOfShapesBaseline = MorphologicalTreeFactory::createTreeOfShapes(context.imageUInt8, ToSInterpolation::SelfDual);
    WorkloadMetrics treeOfShapesMetrics = metricsOf(treeOfShapesBaseline);
    treeOfShapesMetrics.edges = context.maxTreeMetrics.edges;
    results.push_back(benchmarkScenario(
        "pipelines", "tree_of_shapes_boundary_group_direct_area", TimingScope::EndToEnd, context.options.repetitions, treeOfShapesMetrics,
        [&] {
            return directPipeline([&] { return MorphologicalTreeFactory::createTreeOfShapes(context.imageUInt8, ToSInterpolation::SelfDual); },
                                  {AttributeOrGroup{AREA}, AttributeOrGroup{AttributeGroup::BOUNDARY}}, AREA, 1.0 / 16.0);
        },
        [](const PipelineResult<std::uint8_t>& result) { return pipelineChecksum(result); }));
}

} // namespace mmcfilters::benchmarks::api
