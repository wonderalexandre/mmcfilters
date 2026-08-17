#include "ApiBenchmark.hpp"

#include <utility>

namespace mmcfilters::benchmarks::api {
namespace {

template <AltitudeValue T, class Builder>
void addTreeConstruction(Context& context, std::vector<ScenarioResult>& results, std::string name, Builder builder) {
    auto baseline = builder();
    WorkloadMetrics metrics = metricsOf(baseline);
    metrics.edges = context.maxTreeMetrics.edges;
    results.push_back(benchmarkScenario("construction", std::move(name), TimingScope::EndToEnd, context.options.repetitions, metrics, builder,
                                        [](const ValuedMorphologicalTree<T>& tree) { return valuedTreeChecksum(tree); }));
}

} // namespace

void addConstructionScenarios(Context& context, std::vector<ScenarioResult>& results) {
    addTreeConstruction<std::uint8_t>(context, results, "max_tree_uint8", [&] { return MorphologicalTreeFactory::createMaxTree(context.imageUInt8, context.adjacency); });
    addTreeConstruction<std::uint8_t>(context, results, "min_tree_uint8", [&] { return MorphologicalTreeFactory::createMinTree(context.imageUInt8, context.adjacency); });

    if (!atLeast(context.options.profile, Profile::Core)) {
        return;
    }

    addTreeConstruction<ToSGrayLevel>(context, results, "tree_of_shapes_self_dual",
                                     [&] { return MorphologicalTreeFactory::createTreeOfShapes<ToSGrayLevel>(context.imageUInt8, selfDualSpanConvention()); });
    addTreeConstruction<ToSGrayLevel>(context, results, "tree_of_shapes_min4c_max8c",
                                     [&] {
                                          return MorphologicalTreeFactory::createTreeOfShapes<ToSGrayLevel>(
                                              context.imageUInt8,
                                              complementaryTopographicConvention(context.options.rows, context.options.columns, true));
                                      });
    addTreeConstruction<std::uint8_t>(context, results, "self_dual_residual_unrestricted", [&] {
        return MorphologicalTreeFactory::createUnrestrictedResidualTree(
            context.imageUInt8, context.adjacency,
            sdrt::UnrestrictedResidualTreeOptions{});
    });
    addTreeConstruction<std::uint8_t>(context, results, "self_dual_residual_saturated", [&] {
        return MorphologicalTreeFactory::createSaturatedResidualTree(context.imageUInt8, context.adjacency, PixelId{0},
                                                                             sdrt::SaturatedResidualTreeOptions{
                                                                                 sdrt::RowMajorSpatialOrder{}});
    });
    addTreeConstruction<std::int32_t>(context, results, "max_tree_int32", [&] { return MorphologicalTreeFactory::createMaxTree(context.imageInt32, 1.5); });
    addTreeConstruction<float>(context, results, "max_tree_float", [&] { return MorphologicalTreeFactory::createMaxTree(context.imageFloat, 1.5); });

    if (!atLeast(context.options.profile, Profile::Publication)) {
        return;
    }

    addTreeConstruction<ToSGrayLevel>(context, results, "tree_of_shapes_min8c_max4c",
                                      [&] {
                                          return MorphologicalTreeFactory::createTreeOfShapes<ToSGrayLevel>(
                                              context.imageUInt8,
                                              complementaryTopographicConvention(context.options.rows, context.options.columns, false));
                                      });
    addTreeConstruction<ToSGrayLevel>(context, results, "tree_of_shapes_self_dual", [&] {
        return MorphologicalTreeFactory::createTreeOfShapes<ToSGrayLevel>(context.imageUInt8, selfDualSpanConvention());
    });
    addTreeConstruction<std::int32_t>(context, results, "min_tree_int32", [&] { return MorphologicalTreeFactory::createMinTree(context.imageInt32, 1.5); });
    addTreeConstruction<float>(context, results, "min_tree_float", [&] { return MorphologicalTreeFactory::createMinTree(context.imageFloat, 1.5); });
}

} // namespace mmcfilters::benchmarks::api
