#include "ApiBenchmark.hpp"

#include <array>
#include <set>
#include <span>
#include <utility>

namespace mmcfilters::benchmarks::api {
namespace {

struct AttributeBatchResult {
    std::vector<ComputedAttributeData<double>> values;
};

[[nodiscard]] std::uint64_t batchChecksum(const AttributeBatchResult& batch) {
    std::vector<const ComputedAttributeData<double>*> pointers;
    pointers.reserve(batch.values.size());
    for (const ComputedAttributeData<double>& value : batch.values) {
        pointers.push_back(&value);
    }
    return semanticAttributeChecksum<double>(std::span<const ComputedAttributeData<double>* const>(pointers.data(), pointers.size()));
}

void addScalarAttribute(Context& context, std::vector<ScenarioResult>& results, Attribute attribute) {
    const std::string name = "scalar_" + lowerName(attributes::registry::name(attribute));
    results.push_back(benchmarkScenario(
        "attributes", name, TimingScope::EstablishedInput, context.options.repetitions, context.maxTreeMetrics,
        [&] { return AttributeComputation::computeSingleAttribute<double>(context.maxTree, attribute); },
        [](const ComputedAttributeData<double>& result) { return semanticAttributeChecksum(result); }));
}

void addAttributeGroup(Context& context, std::vector<ScenarioResult>& results, AttributeGroup group, bool includeSequentialReference) {
    const std::string name = groupName(group);
    ScenarioResult grouped = benchmarkScenario(
        "attribute_groups", name, TimingScope::EstablishedInput, context.options.repetitions, context.maxTreeMetrics,
        [&] { return AttributeComputation::computeAttributes<double>(context.maxTree, {group}); },
        [](const ComputedAttributeData<double>& result) { return semanticAttributeChecksum(result); });
    const std::uint64_t groupedChecksum = grouped.checksum;
    results.push_back(std::move(grouped));

    if (!includeSequentialReference) {
        return;
    }

    const std::vector<Attribute>& groupAttributes = attributes::registry::attributeGroups().at(group);
    ScenarioResult sequential = benchmarkScenario(
        "attribute_groups", name + "_sequential_scalars", TimingScope::EstablishedInput, context.options.repetitions, context.maxTreeMetrics,
        [&] {
            AttributeBatchResult batch;
            batch.values.reserve(groupAttributes.size());
            for (Attribute attribute : groupAttributes) {
                batch.values.push_back(AttributeComputation::computeSingleAttribute<double>(context.maxTree, attribute));
            }
            return batch;
        },
        [](const AttributeBatchResult& result) { return batchChecksum(result); });
    if (sequential.checksum != groupedChecksum) {
        throw std::runtime_error("Grouped and sequential attribute computations produced different semantic results for " + name + ".");
    }
    results.push_back(std::move(sequential));
}

[[nodiscard]] std::vector<Attribute> expandBundle(const AttributeBundle& bundle) {
    std::set<Attribute> unique;
    for (const AttributeOrGroup& request : bundle.requests) {
        if (const Attribute* attribute = std::get_if<Attribute>(&request)) {
            unique.insert(*attribute);
        } else {
            const AttributeGroup group = std::get<AttributeGroup>(request);
            const std::vector<Attribute>& expansion = attributes::registry::attributeGroups().at(group);
            unique.insert(expansion.begin(), expansion.end());
        }
    }
    return {unique.begin(), unique.end()};
}

void addAttributeBundle(Context& context, std::vector<ScenarioResult>& results, const AttributeBundle& bundle) {
    ScenarioResult grouped = benchmarkScenario(
        "attribute_bundles", bundle.name, TimingScope::EstablishedInput, context.options.repetitions, context.maxTreeMetrics,
        [&] { return AttributeComputation::computeAttributes<double>(context.maxTree, bundle.requests); },
        [](const ComputedAttributeData<double>& result) { return semanticAttributeChecksum(result); });
    const std::uint64_t groupedChecksum = grouped.checksum;
    results.push_back(std::move(grouped));

    const std::vector<Attribute> scalarAttributes = expandBundle(bundle);
    ScenarioResult sequential = benchmarkScenario(
        "attribute_bundles", bundle.name + "_sequential_scalars", TimingScope::EstablishedInput, context.options.repetitions, context.maxTreeMetrics,
        [&] {
            AttributeBatchResult batch;
            batch.values.reserve(scalarAttributes.size());
            for (Attribute attribute : scalarAttributes) {
                batch.values.push_back(AttributeComputation::computeSingleAttribute<double>(context.maxTree, attribute));
            }
            return batch;
        },
        [](const AttributeBatchResult& result) { return batchChecksum(result); });
    if (sequential.checksum != groupedChecksum) {
        throw std::runtime_error("Bundled and sequential attribute computations produced different semantic results for " + bundle.name + ".");
    }
    results.push_back(std::move(sequential));
}

} // namespace

void addAttributeScenarios(Context& context, std::vector<ScenarioResult>& results) {
    std::vector<Attribute> scalarAttributes;
    if (atLeast(context.options.profile, Profile::Publication)) {
        scalarAttributes.reserve(attributes::registry::ATTRIBUTE_METADATA.size());
        for (const attributes::registry::AttributeMetadata& metadata : attributes::registry::ATTRIBUTE_METADATA) {
            scalarAttributes.push_back(metadata.attribute);
        }
    } else if (atLeast(context.options.profile, Profile::Core)) {
        scalarAttributes = {AREA,          LEVEL,           VOLUME,          BOX_WIDTH,       HU_MOMENT_1,
                            BITQUADS_AREA, HEIGHT_NODE,     MAX_DIST,        CONTOUR_PERIMETER};
    } else {
        scalarAttributes = {AREA, LEVEL, MAX_DIST};
    }
    for (Attribute attribute : scalarAttributes) {
        addScalarAttribute(context, results, attribute);
    }

    const std::array<AttributeGroup, 6> groups{AttributeGroup::GRAY_LEVEL, AttributeGroup::SHAPE,         AttributeGroup::MOMENTS,
                                                AttributeGroup::BOUNDARY,   AttributeGroup::TREE_TOPOLOGY, AttributeGroup::ALL};
    for (AttributeGroup group : groups) {
        if (!atLeast(context.options.profile, Profile::Core) && group != AttributeGroup::GRAY_LEVEL && group != AttributeGroup::TREE_TOPOLOGY) {
            continue;
        }
        const bool sequential = atLeast(context.options.profile, Profile::Core) &&
                                (group != AttributeGroup::ALL || atLeast(context.options.profile, Profile::Publication));
        addAttributeGroup(context, results, group, sequential);
    }

    for (const AttributeBundle& bundle : context.options.attributeBundles) {
        addAttributeBundle(context, results, bundle);
    }

    if (!atLeast(context.options.profile, Profile::Core)) {
        return;
    }

    const std::vector<AttributeOrGroup> heterogeneousRequest{AREA, LEVEL, BOX_WIDTH, HU_MOMENT_1, BITQUADS_AREA, HEIGHT_NODE, MAX_DIST, CONTOUR_PERIMETER};
    results.push_back(benchmarkScenario(
        "attributes", "heterogeneous_request", TimingScope::EstablishedInput, context.options.repetitions, context.maxTreeMetrics,
        [&] { return AttributeComputation::computeAttributes<double>(context.maxTree, heterogeneousRequest); },
        [](const ComputedAttributeData<double>& result) { return semanticAttributeChecksum(result); }));

    results.push_back(benchmarkScenario(
        "attributes", "area_delta_radius_2", TimingScope::EstablishedInput, context.options.repetitions, context.maxTreeMetrics,
        [&] {
            return AttributeComputation::computeSingleAttributeWithDelta<double>(context.maxTree, AREA, AltitudeDiff<std::uint8_t>{1}, 2, "last-padding");
        },
        [](const ComputedAttributeDataWithDelta<double>& result) { return deltaAttributeChecksum(result); }));
}

} // namespace mmcfilters::benchmarks::api
