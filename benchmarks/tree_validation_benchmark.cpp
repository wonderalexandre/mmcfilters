#include "mmcfilters/trees/MorphologicalTreeFactory.hpp"

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;
using namespace mmcfilters;

struct Fixture {
    std::vector<NodeId> parent;
    std::vector<NodeId> owner;
    std::vector<std::uint8_t> altitude;
};

struct HigraFixture {
    std::vector<NodeId> parent;
    std::vector<std::uint8_t> altitude;
    int numProperParts = 0;
};

Fixture makeFixture(int numNodes, bool branching) {
    Fixture fixture;
    fixture.parent.resize(static_cast<std::size_t>(numNodes));
    fixture.altitude.assign(static_cast<std::size_t>(numNodes), std::uint8_t{0});
    fixture.parent[0] = 0;

    std::vector<int> childCount(static_cast<std::size_t>(numNodes), 0);
    for (NodeId nodeId = 1; nodeId < numNodes; ++nodeId) {
        const NodeId parentId = branching ? (nodeId - 1) / 2 : nodeId - 1;
        fixture.parent[static_cast<std::size_t>(nodeId)] = parentId;
        ++childCount[static_cast<std::size_t>(parentId)];
    }

    for (NodeId nodeId = 0; nodeId < numNodes; ++nodeId) {
        if (childCount[static_cast<std::size_t>(nodeId)] == 0) {
            fixture.owner.push_back(nodeId);
        }
    }
    return fixture;
}

HigraFixture makeHigraFixture(const Fixture& native) {
    const NodeId numProperParts = static_cast<NodeId>(native.owner.size());
    const NodeId numNodes = static_cast<NodeId>(native.parent.size());
    const NodeId numVertices = numProperParts + numNodes;

    HigraFixture higra;
    higra.parent.resize(static_cast<std::size_t>(numVertices));
    higra.altitude.assign(static_cast<std::size_t>(numVertices), std::uint8_t{0});
    higra.numProperParts = numProperParts;

    for (NodeId properPartId = 0; properPartId < numProperParts; ++properPartId) {
        higra.parent[static_cast<std::size_t>(properPartId)] = numProperParts + native.owner[static_cast<std::size_t>(properPartId)];
    }
    for (NodeId nodeId = 0; nodeId < numNodes; ++nodeId) {
        higra.parent[static_cast<std::size_t>(numProperParts + nodeId)] = numProperParts + native.parent[static_cast<std::size_t>(nodeId)];
    }
    return higra;
}

double milliseconds(Clock::time_point start, Clock::time_point finish) { return std::chrono::duration<double, std::milli>(finish - start).count(); }

void runCase(const char* shape, int numNodes, int commits, bool branching) {
    Fixture fixture = makeFixture(numNodes, branching);

    const auto constructionStart = Clock::now();
    auto weighted = MorphologicalTreeFactory::createFromNativeHierarchy(
        NativeHierarchyView<std::uint8_t>{fixture.parent, fixture.owner, fixture.altitude, 0, std::nullopt, HierarchySemantics{}});
    const auto constructionFinish = Clock::now();

    const auto commitStart = Clock::now();
    for (int iteration = 0; iteration < commits; ++iteration) {
        auto editor = weighted.edit();
        editor.commit();
    }
    const auto commitFinish = Clock::now();

    std::cout << shape << ',' << numNodes << ',' << fixture.owner.size() << ',' << commits << ',' << milliseconds(constructionStart, constructionFinish) << ','
              << milliseconds(commitStart, commitFinish) << '\n';
}

void runHigraCase(const char* shape, int numNodes, int commits, bool branching) {
    const Fixture native = makeFixture(numNodes, branching);
    const HigraFixture fixture = makeHigraFixture(native);

    const auto constructionStart = Clock::now();
    auto weighted = MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(fixture.parent), std::span<const std::uint8_t>(fixture.altitude), 1,
                                                                    fixture.numProperParts, MorphologicalTreeKind::MAX_TREE,
                                                                    RegularGridAdjacency2D(1, fixture.numProperParts, 1.5));
    const auto constructionFinish = Clock::now();

    const auto commitStart = Clock::now();
    for (int iteration = 0; iteration < commits; ++iteration) {
        auto editor = weighted.edit();
        editor.commit();
    }
    const auto commitFinish = Clock::now();

    std::cout << shape << ',' << numNodes << ',' << fixture.numProperParts << ',' << commits << ',' << milliseconds(constructionStart, constructionFinish)
              << ',' << milliseconds(commitStart, commitFinish) << '\n';
}

} // namespace

int main(int argc, char** argv) {
    const int numNodes = argc > 1 ? std::stoi(argv[1]) : 100000;
    const int commits = argc > 2 ? std::stoi(argv[2]) : 10;
    if (numNodes < 2 || commits < 1) {
        std::cerr << "usage: tree_validation_benchmark "
                     "[num_nodes>=2] [commits>=1]\n";
        return EXIT_FAILURE;
    }

    std::cout << "shape,num_nodes,num_proper_parts,"
                 "commits,construction_ms,commits_ms\n";
    runCase("deep", numNodes, commits, false);
    runCase("branching", numNodes, commits, true);
    runHigraCase("higra-deep", numNodes, commits, false);
    runHigraCase("higra-branching", numNodes, commits, true);
    return EXIT_SUCCESS;
}
