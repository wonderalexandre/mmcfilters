#include "mmcfilters/trees/MorphologicalTreeFactory.hpp"
#include <algorithm>
#include <bit>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

using namespace mmcfilters;

// Independent paper oracle: BFS flat zones and literal levelings, then the
// support-inclusion Hasse diagram. No component-tree builder internals.
namespace {
using Mask = std::uint64_t;
struct PaperNode {
    Mask support, proper;
    double eta, delta;
    int parent;
};

// Literal Eqs. (6), (8), (16)-(17): recompute flat zones by BFS after every leveling.
// No component trees, union-find, or agenda.
template <class T>
std::vector<PaperNode> literalPaper(const std::vector<T>& input, const std::vector<std::vector<int>>& neighbors) {
    const int n = static_cast<int>(input.size());
    const Mask full = (Mask{1} << n) - 1;
    auto current = input;
    std::vector<PaperNode> nodes;
    std::size_t lastArea = 0;
    for (;;) {
        std::vector<Mask> zones;
        Mask visited = 0;
        for (int seed = 0; seed < n; ++seed) {
            if (visited & (Mask{1} << seed))
                continue;
            std::vector<int> queue{seed};
            Mask zone = Mask{1} << seed;
            visited |= zone;
            for (std::size_t i = 0; i < queue.size(); ++i) {
                for (int q : neighbors[queue[i]]) {
                    const Mask bit = Mask{1} << q;
                    if (!(visited & bit) && current[q] == current[seed]) {
                        visited |= bit;
                        zone |= bit;
                        queue.push_back(q);
                    }
                }
            }
            zones.push_back(zone);
        }
        if (zones.size() == 1) {
            nodes.push_back({full, 0, double(current.front()), double(current.front()), -1});
            break;
        }
        Mask selected = 0;
        T selectedLevel{}, target{};
        auto best = std::tuple{std::numeric_limits<int>::max(), std::numeric_limits<int>::max()};
        for (Mask zone : zones) {
            const T level = current[std::countr_zero(zone)];
            bool maximum = true, minimum = true, hasBoundary = false;
            T lower{}, upper{};
            int spatialRank = n;
            for (int p = 0; p < n; ++p) {
                if (!(zone & (Mask{1} << p)))
                    continue;
                spatialRank = std::min(spatialRank, p);
                for (int q : neighbors[p]) {
                    if (zone & (Mask{1} << q))
                        continue;
                    if (!hasBoundary) {
                        lower = upper = current[q];
                        hasBoundary = true;
                    }
                    lower = std::min(lower, current[q]);
                    upper = std::max(upper, current[q]);
                    maximum = maximum && current[q] < level;
                    minimum = minimum && current[q] > level;
                }
            }
            if (!hasBoundary || (!maximum && !minimum))
                continue;
            const auto key = std::tuple{std::popcount(zone), spatialRank};
            if (key < best) {
                best = key;
                selected = zone;
                selectedLevel = level;
                target = maximum ? upper : lower;
            }
        }
        if (!selected)
            throw std::runtime_error("literal paper: no extremum");
        const auto area = static_cast<std::size_t>(std::popcount(selected));
        if (area < lastArea)
            throw std::runtime_error("literal paper: decreasing selected area");
        lastArea = area;
        nodes.push_back({selected, 0, double(selectedLevel), double(selectedLevel) - double(target), -1});
        for (int p = 0; p < n; ++p)
            if (selected & (Mask{1} << p))
                current[p] = target;
    }
    // Definition (23): derive the inclusion Hasse diagram independently.
    const int root = static_cast<int>(nodes.size()) - 1;
    for (int i = 0; i < root; ++i) {
        int parentArea = n + 1;
        for (int j = 0; j <= root; ++j) {
            if (nodes[i].support == nodes[j].support)
                continue;
            if ((nodes[i].support & nodes[j].support) != nodes[i].support)
                continue;
            const int area = std::popcount(nodes[j].support);
            if (area < parentArea) {
                nodes[i].parent = j;
                parentArea = area;
            }
        }
        if (nodes[i].parent < 0)
            throw std::runtime_error("literal paper: missing inclusion parent");
    }
    nodes[root].parent = root;
    for (int i = 0; i <= root; ++i) {
        nodes[i].proper = nodes[i].support;
        for (int j = 0; j < root; ++j)
            if (nodes[j].parent == i)
                nodes[i].proper &= ~nodes[j].support;
        if (!nodes[i].proper)
            throw std::runtime_error("literal paper: empty proper part");
        const double expectedDelta = i == root ? nodes[i].eta : nodes[i].eta - nodes[nodes[i].parent].eta;
        if (nodes[i].delta != expectedDelta)
            throw std::runtime_error("literal paper: signed residue differs from Eq. (26)");
        for (int p = 0; p < n; ++p)
            if (nodes[i].proper & (Mask{1} << p)) {
                if (double(input[p]) != nodes[i].eta)
                    throw std::runtime_error("literal paper: Eq. (25) failed");
            }
    }
    return nodes;
}

std::size_t cases = 0, nodesChecked = 0;
template <class T> void checkCase(int rows, int columns, const std::vector<T>& values, const RegularGridAdjacency2D& adjacency) {
    const int n = rows * columns;
    std::vector<std::vector<int>> neighbors(n);
    for (int p = 0; p < n; ++p)
        for (int q : adjacency.getNeighborIndices(p))
            neighbors[p].push_back(q);
    const auto expected = literalPaper(values, neighbors);
    auto image = Image<T>::create(rows, columns);
    for (int p = 0; p < n; ++p)
        (*image)[p] = values[p];
    const auto actual = MorphologicalTreeFactory::createExtremumResidualTreeByCTs(image, adjacency);
    const auto& topology = actual.topology();
    topology.validateConnectedRootedTree();
    auto fail = [&](const char* message) {
        std::cerr << "case=" << cases << " size=" << rows << 'x' << columns << " values=";
        for (T value : values)
            std::cerr << double(value) << ',';
        std::cerr << '\n';
        throw std::runtime_error(message);
    };
    if (topology.numNodes() != static_cast<int>(expected.size()))
        fail("node count");
    std::map<Mask, NodeId> actualBySupport;
    for (NodeId id : topology.aliveNodeIds()) {
        Mask support = 0;
        for (PixelId p : topology.nodeSupport(id))
            support |= Mask{1} << p;
        if (!actualBySupport.emplace(support, id).second)
            fail("duplicate support");
    }
    for (std::size_t stepIndex = 0; stepIndex < expected.size(); ++stepIndex) {
        const auto& node = expected[stepIndex];
        const auto found = actualBySupport.find(node.support);
        if (found == actualBySupport.end())
            fail("support family");
        const NodeId id = found->second;
        // This synchronized factory numbers the root as 0 and recorded supports in chronological order from 1.
        const NodeId chronologicalId = stepIndex + 1 == expected.size() ? 0 : static_cast<NodeId>(stepIndex + 1);
        if (id != chronologicalId)
            fail("selection sequence");
        if (topology.parent(id) != actualBySupport.at(expected[node.parent].support))
            fail("inclusion parent");
        Mask proper = 0;
        for (PixelId p : topology.properPart(id))
            proper |= Mask{1} << p;
        if (proper != node.proper)
            fail("proper part");
        if (double(actual.nodeAltitude(id)) != node.eta)
            fail("eta altitude");
        if (double(actual.nodeResidue(id)) != node.delta)
            fail("delta residue");
    }
    const auto reconstructed = actual.reconstructFromNodeAltitudes();
    for (int p = 0; p < n; ++p) {
        double sum = 0;
        for (const auto& node : expected)
            if (node.support & (Mask{1} << p))
                sum += node.delta;
        if ((*reconstructed)[p] != values[p] || sum != double(values[p]))
            fail("reconstruction Eqs. (22), (28)");
    }
    ++cases;
    nodesChecked += expected.size();
}

} // namespace

int main() {
    try {
        {
            cases = nodesChecked = 0;
            std::cout << "extremum canonical" << '\n';
            for (int n : {6, 9}) {
                const int rows = n == 6 ? 2 : 3, columns = 3;
                int count = 1;
                for (int i = 0; i < n; ++i)
                    count *= 3;
                for (int variant = 0; variant < 2; ++variant) {
                    const RegularGridAdjacency2D adjacency(rows, columns, variant == 1 ? 1.5 : 1.0);
                    for (int code = 0; code < count; ++code) {
                        std::vector<std::uint8_t> values(n);
                        int c = code;
                        for (auto& value : values) {
                            value = c % 3;
                            c /= 3;
                        }
                        checkCase(rows, columns, values, adjacency);
                    }
                }
            }
            std::cout << "exhaustive ternary 2x3 and 3x3: " << cases << " cases passed\n";
            std::mt19937 random(20271004);
            for (int i = 0; i < 256; ++i) {
                const RegularGridAdjacency2D adjacency(4, 4, i % 3 == 0 ? 3.2 : i % 3 == 1 ? 1.0 : 1.5);
                std::vector<std::int32_t> values(16), inverted(16);
                std::vector<double> fractional(16);
                for (int p = 0; p < 16; ++p) {
                    values[p] = int(random() % 17) - 8;
                    inverted[p] = -values[p];
                    fractional[p] = 0.25 * values[p] - 0.5;
                }
                checkCase(4, 4, values, adjacency);
                checkCase(4, 4, inverted, adjacency);
                checkCase(4, 4, fractional, adjacency);
            }
            for (std::int32_t value : {-4, 0, 7})
                checkCase(1, 1, std::vector<std::int32_t>{value}, RegularGridAdjacency2D(1, 1, 1.0));
            std::cout << "PASS total=" << cases << " nodes=" << nodesChecked
                      << " (selection sequence, support family, inclusion parents, proper parts, eta, signed delta, two reconstructions)\n";
        }
    } catch (const std::exception& error) {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
}
