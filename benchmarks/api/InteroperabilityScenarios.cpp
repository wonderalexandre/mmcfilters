#include "ApiBenchmark.hpp"

#include "mmcfilters/trees/saliency/HierarchySaliencyMap.hpp"

#include <span>
#include <vector>

namespace mmcfilters::benchmarks::api {
namespace {

template <typename Value> [[nodiscard]] std::uint64_t vectorChecksum(const std::vector<Value>& values) noexcept {
    Fnv1a64 hash;
    hash.appendVector(values);
    return hash.value();
}

} // namespace

void addInteroperabilityScenarios(Context& context, std::vector<ScenarioResult>& results) {
    results.push_back(benchmarkScenario(
        "interoperability", "reconstruct_max_tree", TimingScope::EstablishedInput, context.options.repetitions, context.maxTreeMetrics,
        [&] { return context.maxTree.reconstructFromNodeAltitudes(); }, [](const ImageUInt8Ptr& result) { return imageChecksum(result); }));

    results.push_back(benchmarkScenario(
        "interoperability", "export_higra", TimingScope::EstablishedInput, context.options.repetitions, context.maxTreeMetrics,
        [&] { return context.maxTree.exportHigraHierarchy(); },
        [](const std::pair<std::vector<NodeId>, std::vector<std::uint8_t>>& result) {
            Fnv1a64 hash;
            appendExportedHierarchy(hash, result);
            return hash.value();
        }));

    if (!atLeast(context.options.profile, Profile::Core)) {
        return;
    }

    const auto exported = context.maxTree.exportHigraHierarchy();
    auto importedBaseline = MorphologicalTreeFactory::createFromHigraParent<std::uint8_t>(
        std::span<const NodeId>(exported.first), std::span<const std::uint8_t>(exported.second), context.options.rows, context.options.columns,
        MorphologicalTreeKind::MaxTree, context.adjacency);
    WorkloadMetrics importedMetrics = metricsOf(importedBaseline);
    importedMetrics.edges = context.maxTreeMetrics.edges;
    results.push_back(benchmarkScenario(
        "interoperability", "import_higra", TimingScope::EstablishedInput, context.options.repetitions, importedMetrics,
        [&] {
            return MorphologicalTreeFactory::createFromHigraParent<std::uint8_t>(
                std::span<const NodeId>(exported.first), std::span<const std::uint8_t>(exported.second), context.options.rows, context.options.columns,
                MorphologicalTreeKind::MaxTree, context.adjacency);
        },
        [](const ValuedMorphologicalTree<std::uint8_t>& result) { return valuedTreeChecksum(result); }));

    results.push_back(benchmarkScenario(
        "interoperability", "project_area_to_higra", TimingScope::EstablishedInput, context.options.repetitions, context.maxTreeMetrics,
        [&] { return AttributeComputation::projectNodeValuesToExportedHigra<double>(context.maxTree, context.area.first, context.area.second); },
        [](const std::vector<double>& result) { return vectorChecksum(result); }));

    results.push_back(benchmarkScenario(
        "interoperability", "map_area_to_image", TimingScope::EstablishedInput, context.options.repetitions, context.maxTreeMetrics,
        [&] { return AttributeComputation::computeAttributeMapping<double>(context.maxTree, Area); },
        [](const ImagePtr<double>& result) { return imageChecksum(result); }));

    results.push_back(benchmarkScenario(
        "interoperability", "valuedTree_view_attributes", TimingScope::EstablishedInput, context.options.repetitions, context.maxTreeMetrics,
        [&] { return AttributeComputation::computeAttributes<double>(context.maxTree.asView(), {Area, MeanGrayLevel, MaxDist}); },
        [](const ComputedAttributeData<double>& result) { return semanticAttributeChecksum(result); }));

    WorkloadMetrics saliencyMetrics = context.maxTreeMetrics;
    saliencyMetrics.edges = context.maxTreeMetrics.edges;
    results.push_back(benchmarkScenario(
        "interoperability", "topological_saliency", TimingScope::EstablishedInput, context.options.repetitions, saliencyMetrics,
        [&] { return HierarchySaliencyMap::computeTopologicalLevelEdgeMap(context.maxTree.topology(), context.adjacency); },
        [](const EdgeSaliencyMap<int>& result) { return edgeMapChecksum(result); }));

    results.push_back(benchmarkScenario(
        "interoperability", "normalized_altitude_saliency", TimingScope::EstablishedInput, context.options.repetitions, saliencyMetrics,
        [&] { return HierarchySaliencyMap::computeNormalizedAltitudeEdgeMap(context.maxTree, context.adjacency); },
        [](const EdgeSaliencyMap<double>& result) { return edgeMapChecksum(result); }));

    if (!atLeast(context.options.profile, Profile::Publication)) {
        return;
    }

    results.push_back(benchmarkScenario(
        "interoperability", "attributes_in_higra_space", TimingScope::EstablishedInput, context.options.repetitions, importedMetrics,
        [&] { return AttributeComputation::computeAttributes<double>(importedBaseline, {Area, MeanGrayLevel, MaxDist}, NodeIdSpace::Higra); },
        [](const ComputedAttributeData<double>& result) { return semanticAttributeChecksum(result); }));
}

} // namespace mmcfilters::benchmarks::api
