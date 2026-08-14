#include "ApiBenchmark.hpp"

#include "mmcfilters/filters/AttributeFilters.hpp"

#include <algorithm>
#include <cstdint>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace mmcfilters::benchmarks::api {
namespace {

template <AltitudeValue T> struct PipelineResult {
    ValuedMorphologicalTree<T> tree;
    ComputedAttributeData<double> attributes;
    ImagePtr<T> image;
};

template <AltitudeValue T> [[nodiscard]] std::uint64_t pipelineChecksum(const PipelineResult<T>& result) {
    Fnv1a64 hash;
    hash.append(valuedTreeChecksum(result.tree));
    hash.append(semanticAttributeChecksum(result.attributes));
    appendImage(hash, result.image);
    return hash.value();
}

template <AltitudeValue T>
[[nodiscard]] std::vector<bool> preservationDecisionsFromAttribute(const ValuedMorphologicalTree<T>& tree,
                                                       const ComputedAttributeData<double>& attributes, Attribute attribute,
                                                       double threshold) {
    const int stride = attributes.first.NUM_ATTRIBUTES;
    const int offset = attributes.first.getIndex(attribute);
    std::vector<bool> preservationDecisions(static_cast<std::size_t>(tree.topology().numInternalNodeSlots()), false);
    for (NodeId node : tree.topology().aliveNodeIds()) {
        const std::size_t row = static_cast<std::size_t>(node);
        preservationDecisions[row] = attributes.second[row * static_cast<std::size_t>(stride) + static_cast<std::size_t>(offset)] >= threshold;
    }
    preservationDecisions[static_cast<std::size_t>(tree.topology().root())] = true;
    return preservationDecisions;
}

template <AltitudeValue T>
[[nodiscard]] double relativeThreshold(const ValuedMorphologicalTree<T>& tree, const ComputedAttributeData<double>& attributes,
                                       Attribute attribute, double fraction) {
    const int stride = attributes.first.NUM_ATTRIBUTES;
    const int offset = attributes.first.getIndex(attribute);
    double maximum = 0.0;
    for (NodeId node : tree.topology().aliveNodeIds()) {
        const std::size_t row = static_cast<std::size_t>(node);
        maximum = std::max(maximum, attributes.second[row * static_cast<std::size_t>(stride) + static_cast<std::size_t>(offset)]);
    }
    return maximum * fraction;
}

template <class TreeFactory>
[[nodiscard]] auto directPipeline(TreeFactory&& createTree, std::vector<AttributeOrGroup> request, Attribute preservationAttribute,
                                  double thresholdFraction) {
    auto tree = createTree();
    using Altitude = typename std::remove_cvref_t<decltype(tree)>::AltitudeType;
    auto attributes = AttributeComputation::computeAttributes<double>(tree, request);
    const double threshold = relativeThreshold(tree, attributes, preservationAttribute, thresholdFraction);
    auto preservationDecisions = preservationDecisionsFromAttribute(tree, attributes, preservationAttribute, threshold);
    DirectAttributeFilter<Altitude> filter(tree);
    auto image = filter.applyDirectAttributeFilter(NodePreservationMask(std::move(preservationDecisions)));
    return PipelineResult<Altitude>{std::move(tree), std::move(attributes), std::move(image)};
}

} // namespace

void addPipelineScenarios(Context& context, std::vector<ScenarioResult>& results) {
    results.push_back(benchmarkScenario(
        "pipelines", "max_tree_area_direct", TimingScope::EndToEnd, context.options.repetitions, context.maxTreeMetrics,
        [&] {
            return directPipeline(
                [&] { return MorphologicalTreeFactory::createMaxTree(context.imageUInt8, context.adjacency); }, {AttributeOrGroup{Area}}, Area, 1.0 / 16.0);
        },
        [](const PipelineResult<std::uint8_t>& result) { return pipelineChecksum(result); }));

    if (!atLeast(context.options.profile, Profile::Core)) {
        return;
    }

    results.push_back(benchmarkScenario(
        "pipelines", "max_tree_max_dist_direct", TimingScope::EndToEnd, context.options.repetitions, context.maxTreeMetrics,
        [&] {
            return directPipeline(
                [&] { return MorphologicalTreeFactory::createMaxTree(context.imageUInt8, context.adjacency); }, {AttributeOrGroup{MaxDist}}, MaxDist, 0.25);
        },
        [](const PipelineResult<std::uint8_t>& result) { return pipelineChecksum(result); }));

    results.push_back(benchmarkScenario(
        "pipelines", "max_tree_shape_group_direct_area", TimingScope::EndToEnd, context.options.repetitions, context.maxTreeMetrics,
        [&] {
            return directPipeline([&] { return MorphologicalTreeFactory::createMaxTree(context.imageUInt8, context.adjacency); },
                                  {AttributeOrGroup{AttributeGroup::Shape}}, Area, 1.0 / 16.0);
        },
        [](const PipelineResult<std::uint8_t>& result) { return pipelineChecksum(result); }));

    results.push_back(benchmarkScenario(
        "pipelines", "min_tree_area_pruning_min", TimingScope::EndToEnd, context.options.repetitions, metricsOf(context.minTree),
        [&] {
            auto tree = MorphologicalTreeFactory::createMinTree(context.imageUInt8, context.adjacency);
            auto attributes = AttributeComputation::computeSingleAttribute<double>(tree, Area);
            AttributeFilters<std::uint8_t> filters(tree);
            auto image = filters.filteringByPruningMin(attributes.second.data(), context.areaThreshold);
            return PipelineResult<std::uint8_t>{std::move(tree), std::move(attributes), std::move(image)};
        },
        [](const PipelineResult<std::uint8_t>& result) { return pipelineChecksum(result); }));

    if (!atLeast(context.options.profile, Profile::Publication)) {
        return;
    }

    auto treeOfShapesBaseline = MorphologicalTreeFactory::createTreeOfShapes(context.imageUInt8);
    WorkloadMetrics treeOfShapesMetrics = metricsOf(treeOfShapesBaseline);
    treeOfShapesMetrics.edges = context.maxTreeMetrics.edges;
    results.push_back(benchmarkScenario(
        "pipelines", "tree_of_shapes_boundary_group_direct_area", TimingScope::EndToEnd, context.options.repetitions, treeOfShapesMetrics,
        [&] {
            return directPipeline([&] { return MorphologicalTreeFactory::createTreeOfShapes(context.imageUInt8); },
                                  {AttributeOrGroup{Area}, AttributeOrGroup{AttributeGroup::Boundary}}, Area, 1.0 / 16.0);
        },
        [](const PipelineResult<ToSGrayLevel>& result) { return pipelineChecksum(result); }));
}

} // namespace mmcfilters::benchmarks::api
