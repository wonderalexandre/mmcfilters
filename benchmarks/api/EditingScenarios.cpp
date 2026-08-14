#include "ApiBenchmark.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace mmcfilters::benchmarks::api {
namespace {

struct EditState {
    ValuedMorphologicalTree<std::uint8_t> tree;
    std::vector<NodeId> candidates;
};

[[nodiscard]] EditState prepareEditState(Context& context, int maximumCandidates) {
    auto tree = MorphologicalTreeFactory::createMaxTree(context.imageUInt8, context.adjacency);
    std::vector<NodeId> candidates;
    candidates.reserve(static_cast<std::size_t>(maximumCandidates));
    const MorphologicalTree& topology = tree.topology();
    for (NodeId node : topology.aliveNodeIds()) {
        if (!topology.isRoot(node) && topology.numChildren(node) == 0) {
            candidates.push_back(node);
            if (static_cast<int>(candidates.size()) == maximumCandidates) {
                break;
            }
        }
    }
    if (candidates.empty()) {
        for (NodeId node : topology.aliveNodeIds()) {
            if (!topology.isRoot(node)) {
                candidates.push_back(node);
                break;
            }
        }
    }
    if (candidates.empty()) {
        throw std::runtime_error("Editing benchmark requires at least one non-root node.");
    }
    return EditState{std::move(tree), std::move(candidates)};
}

[[nodiscard]] std::uint64_t editedTreeChecksum(const EditState& state, int operations) {
    Fnv1a64 hash;
    hash.append(operations);
    hash.append(valuedTreeChecksum(state.tree));
    return hash.value();
}

} // namespace

void addEditingScenarios(Context& context, std::vector<ScenarioResult>& results) {
    if (context.maxTree.topology().numNodes() <= 1) {
        if (context.options.suites.contains("all")) {
            return;
        }
        throw std::runtime_error("The editing suite requires an input that produces at least one non-root node.");
    }

    results.push_back(benchmarkPreparedScenario(
        "editing", "safe_prune_leaf", context.options.repetitions, context.maxTreeMetrics, [&] { return prepareEditState(context, 1); },
        [](EditState& state) {
            state.tree.pruneNode(state.candidates.front());
            return 1;
        },
        [](const EditState& state, int operations) { return editedTreeChecksum(state, operations); }));

    results.push_back(benchmarkPreparedScenario(
        "editing", "safe_merge_leaf", context.options.repetitions, context.maxTreeMetrics, [&] { return prepareEditState(context, 1); },
        [](EditState& state) {
            state.tree.mergeNodeIntoParent(state.candidates.front());
            return 1;
        },
        [](const EditState& state, int operations) { return editedTreeChecksum(state, operations); }));

    if (!atLeast(context.options.profile, Profile::Core)) {
        return;
    }

    const int batchSize = atLeast(context.options.profile, Profile::Publication) ? 64 : 16;
    results.push_back(benchmarkPreparedScenario(
        "editing", "editor_batch_merge_commit", context.options.repetitions, context.maxTreeMetrics,
        [&] { return prepareEditState(context, batchSize); },
        [](EditState& state) {
            auto editor = state.tree.edit();
            for (NodeId candidate : state.candidates) {
                editor.mergeNodeIntoParent(candidate);
            }
            editor.commit();
            return static_cast<int>(state.candidates.size());
        },
        [](const EditState& state, int operations) { return editedTreeChecksum(state, operations); }));

    if (!atLeast(context.options.profile, Profile::Publication)) {
        return;
    }

    results.push_back(benchmarkPreparedScenario(
        "editing", "safe_merge_sequence", context.options.repetitions, context.maxTreeMetrics,
        [&] { return prepareEditState(context, batchSize); },
        [](EditState& state) {
            for (NodeId candidate : state.candidates) {
                state.tree.mergeNodeIntoParent(candidate);
            }
            return static_cast<int>(state.candidates.size());
        },
        [](const EditState& state, int operations) { return editedTreeChecksum(state, operations); }));
}

} // namespace mmcfilters::benchmarks::api
