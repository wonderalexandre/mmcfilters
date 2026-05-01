#include "mmcfilters/contours/ContoursComputedIncrementally.hpp"
#include "mmcfilters/trees/MorphologicalTree.hpp"
#include "stb_image.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using namespace mmcfilters;

namespace {

ImageUInt8Ptr makeBenchmarkImage(int rows, int cols) {
    auto image = ImageUInt8::create(rows, cols);
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            const int idx = row * cols + col;
            const int radial = (row - rows / 2) * (row - rows / 2) + (col - cols / 2) * (col - cols / 2);
            (*image)[idx] = static_cast<uint8_t>((radial / 97 + row * 3 + col * 5) & 0xff);
        }
    }
    return image;
}

ImageUInt8Ptr loadGrayscaleImage(const std::filesystem::path& path) {
    int width = 0;
    int height = 0;
    int channels = 0;
    unsigned char* raw = stbi_load(path.string().c_str(), &width, &height, &channels, 1);
    if (!raw) {
        throw std::runtime_error("Failed to load image: " + path.string());
    }

    auto image = ImageUInt8::create(height, width);
    const int size = width * height;
    for (int i = 0; i < size; ++i) {
        (*image)[i] = raw[i];
    }
    stbi_image_free(raw);
    return image;
}

template <typename Fn>
auto timed(Fn&& fn) {
    const auto start = std::chrono::steady_clock::now();
    auto result = fn();
    const auto end = std::chrono::steady_clock::now();
    const auto micros = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    return std::pair{std::move(result), static_cast<double>(micros) / 1000.0};
}

template <typename Fn>
void measure(const std::string& label, int repeats, Fn&& fn) {
    std::uint64_t checksum = 0;
    const auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < repeats; ++i) {
        checksum += static_cast<std::uint64_t>(fn());
    }
    const auto end = std::chrono::steady_clock::now();
    const auto micros = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "  " << label << ": " << (static_cast<double>(micros) / repeats / 1000.0) << " ms/run"
              << " checksum=" << checksum << '\n';
}

std::size_t sumContour(const ContoursComputedIncrementally::IncrementalContours& contours, NodeId node) {
    std::size_t sum = 0;
    for (int pixel : contours.getContour(node)) {
        sum += static_cast<std::size_t>(pixel + 1);
    }
    return sum;
}

std::vector<NodeId> aliveNodesInTreeOrder(const MorphologicalTree& tree) {
    std::vector<NodeId> nodes;
    nodes.reserve(static_cast<std::size_t>(tree.getNumNodes()));
    for (NodeId node : tree.getAliveNodeIds()) {
        nodes.push_back(node);
    }
    return nodes;
}

std::vector<NodeId> shuffledAliveNodes(const MorphologicalTree& tree) {
    auto nodes = aliveNodesInTreeOrder(tree);
    std::mt19937 rng(0x5eed);
    std::shuffle(nodes.begin(), nodes.end(), rng);
    return nodes;
}

std::size_t sumContoursInOrder(const std::vector<NodeId>& nodes,
                               const ContoursComputedIncrementally::IncrementalContours& contours) {
    std::size_t sum = 0;
    for (NodeId node : nodes) {
        sum += sumContour(contours, node);
    }
    return sum;
}

std::size_t sumContoursByNode(const ContoursComputedIncrementally::IncrementalContours& contours) {
    std::size_t sum = 0;
    for (auto [node, contour] : contours.contoursByNode()) {
        (void)node;
        for (int pixel : contour) {
            sum += static_cast<std::size_t>(pixel + 1);
        }
    }
    return sum;
}

double mib(std::size_t bytes) {
    return static_cast<double>(bytes) / (1024.0 * 1024.0);
}

std::size_t storageChecksum(const ContoursComputedIncrementally::IncrementalContours::StorageStats& stats) {
    return stats.addDeltaValues + stats.removeDeltaValues + stats.cachedContourValues +
           stats.cachedContourCapacity + stats.cachedContourReadyNodes + stats.approxAllocatedBytes;
}

void printStorageStats(const std::string& label,
                       const ContoursComputedIncrementally::IncrementalContours& contours) {
    const auto stats = contours.storageStats();
    std::cout << "    " << label << ": "
              << "add_delta_values=" << stats.addDeltaValues
              << " remove_delta_values=" << stats.removeDeltaValues
              << " cached_contour_values=" << stats.cachedContourValues
              << " cached_contour_capacity=" << stats.cachedContourCapacity
              << " cached_ready_nodes=" << stats.cachedContourReadyNodes
              << " approx_allocated=" << mib(stats.approxAllocatedBytes) << " MiB\n";
}

void printBenchmarkMethods() {
    std::cout << "methods:\n"
              << "  extract only -> ContoursComputedIncrementally::extractCompactContours(tree)\n"
              << "  root via getContour(root) -> contours.getContour(root)\n"
              << "  iterate all via getContour(node) -> for node in order: contours.getContour(node)\n"
              << "  iterate all via contoursByNode() -> for (auto [node, contour] : contours.contoursByNode())\n"
              << "  iterate all random via getContour(node) -> shuffled nodes + contours.getContour(node)\n"
              << "  materializeAll + iterate -> contours.materializeAll(); then contours.getContour(node)\n";
}

void printTreeStats(const MorphologicalTree& tree, double buildMs) {
    int maxChildren = 0;
    long long totalChildren = 0;
    long long totalDirectProperParts = 0;
    int maxDirectProperParts = 0;
    int minDirectProperParts = std::numeric_limits<int>::max();

    for (NodeId node : tree.getAliveNodeIds()) {
        const int children = tree.getNumChildren(node);
        const int properParts = tree.getNumProperParts(node);
        totalChildren += children;
        totalDirectProperParts += properParts;
        maxChildren = std::max(maxChildren, children);
        maxDirectProperParts = std::max(maxDirectProperParts, properParts);
        minDirectProperParts = std::min(minDirectProperParts, properParts);
    }

    const int numNodes = tree.getNumNodes();
    if (numNodes == 0) {
        minDirectProperParts = 0;
    }

    std::cout << "  build: " << buildMs << " ms\n"
              << "  nodes: live=" << numNodes
              << " slots=" << tree.getNumInternalNodeSlots()
              << " leaves=" << tree.getNumLeafNodes()
              << " proper_parts=" << tree.getNumTotalProperParts() << '\n'
              << "  branching: avg_children=" << (numNodes > 0 ? static_cast<double>(totalChildren) / numNodes : 0.0)
              << " max_children=" << maxChildren << '\n'
              << "  direct_proper_parts: avg=" << (numNodes > 0 ? static_cast<double>(totalDirectProperParts) / numNodes : 0.0)
              << " min=" << minDirectProperParts
              << " max=" << maxDirectProperParts << '\n';
}

void printSingleRunBreakdown(const MorphologicalTree& tree,
                             const std::vector<NodeId>& nodesInTreeOrder,
                             const std::vector<NodeId>& nodesInRandomOrder) {
    std::cout << "  single-run breakdown:\n";

    auto [extractedContours, extractMs] = timed([&]() {
        return ContoursComputedIncrementally::extractCompactContours(tree);
    });
    std::cout << "    extract only: " << extractMs << " ms"
              << " checksum=" << storageChecksum(extractedContours.storageStats()) << '\n';
    printStorageStats("after extract", extractedContours);

    auto [treeOrderSum, treeOrderIterMs] = timed([&]() {
        return sumContoursInOrder(nodesInTreeOrder, extractedContours);
    });
    std::cout << "    iterate all via getContour(node) after extract: " << treeOrderIterMs << " ms"
              << " checksum=" << treeOrderSum << '\n';
    printStorageStats("after tree-order iteration", extractedContours);

    auto [byNodeContours, byNodeExtractMs] = timed([&]() {
        return ContoursComputedIncrementally::extractCompactContours(tree);
    });
    const auto byNodeExtractStats = byNodeContours.storageStats();
    auto [byNodeSum, byNodeIterMs] = timed([&]() {
        return sumContoursByNode(byNodeContours);
    });
    std::cout << "    extract only for contoursByNode(): " << byNodeExtractMs << " ms"
              << " checksum=" << storageChecksum(byNodeExtractStats) << '\n';
    std::cout << "    iterate all via contoursByNode() after extract: " << byNodeIterMs << " ms"
              << " checksum=" << byNodeSum << '\n';
    printStorageStats("after contoursByNode()", byNodeContours);

    auto [randomContours, randomExtractMs] = timed([&]() {
        return ContoursComputedIncrementally::extractCompactContours(tree);
    });
    const auto randomExtractStats = randomContours.storageStats();
    auto [randomOrderSum, randomOrderIterMs] = timed([&]() {
        return sumContoursInOrder(nodesInRandomOrder, randomContours);
    });
    std::cout << "    extract only for random order: " << randomExtractMs << " ms"
              << " checksum=" << storageChecksum(randomExtractStats) << '\n';
    std::cout << "    iterate all via getContour(node) random order after extract: " << randomOrderIterMs << " ms"
              << " checksum=" << randomOrderSum << '\n';
    printStorageStats("after random-order iteration", randomContours);

    auto [globalCacheContours, globalCacheExtractMs] = timed([&]() {
        return ContoursComputedIncrementally::extractCompactContours(tree);
    });
    const auto globalCacheExtractStats = globalCacheContours.storageStats();
    auto [globalCacheChecksum, globalCacheMs] = timed([&]() {
        globalCacheContours.materializeAll();
        return storageChecksum(globalCacheContours.storageStats());
    });
    auto [globalCacheIterSum, globalCacheIterMs] = timed([&]() {
        return sumContoursInOrder(nodesInTreeOrder, globalCacheContours);
    });
    std::cout << "    extract only for materialize all: " << globalCacheExtractMs << " ms"
              << " checksum=" << storageChecksum(globalCacheExtractStats) << '\n';
    std::cout << "    materializeAll() after extract: " << globalCacheMs << " ms"
              << " checksum=" << globalCacheChecksum << '\n';
    std::cout << "    iterate all via getContour(node) after materializeAll(): " << globalCacheIterMs << " ms"
              << " checksum=" << globalCacheIterSum << '\n';
    printStorageStats("after materialize all", globalCacheContours);
}

void runCase(const std::string& label, const MorphologicalTree& tree, double buildMs, int repeats) {
    const NodeId root = tree.getRoot();
    const auto nodesInTreeOrder = aliveNodesInTreeOrder(tree);
    const auto nodesInRandomOrder = shuffledAliveNodes(tree);

    std::cout << "\n[" << label << "]\n";
    printTreeStats(tree, buildMs);

    measure("extract only [extractCompactContours]", repeats, [&]() {
        auto contours = ContoursComputedIncrementally::extractCompactContours(tree);
        return storageChecksum(contours.storageStats());
    });

    measure("extract + root via getContour(root)", repeats, [&]() {
        auto contours = ContoursComputedIncrementally::extractCompactContours(tree);
        return sumContour(contours, root);
    });

    measure("extract + iterate all via getContour(node)", repeats, [&]() {
        auto contours = ContoursComputedIncrementally::extractCompactContours(tree);
        return sumContoursInOrder(nodesInTreeOrder, contours);
    });

    measure("extract + iterate all via contoursByNode()", repeats, [&]() {
        auto contours = ContoursComputedIncrementally::extractCompactContours(tree);
        return sumContoursByNode(contours);
    });

    measure("extract + iterate all via getContour(node) random order", repeats, [&]() {
        auto contours = ContoursComputedIncrementally::extractCompactContours(tree);
        return sumContoursInOrder(nodesInRandomOrder, contours);
    });

    measure("extract + materializeAll() + iterate via getContour(node)", repeats, [&]() {
        auto contours = ContoursComputedIncrementally::extractCompactContours(tree);
        contours.materializeAll();
        return sumContoursInOrder(nodesInTreeOrder, contours);
    });

    printSingleRunBreakdown(tree, nodesInTreeOrder, nodesInRandomOrder);
}

} // namespace

int main(int argc, char** argv) {
    ImageUInt8Ptr image;
    int repeats = 5;

    if (argc > 1 && std::filesystem::exists(argv[1])) {
        image = loadGrayscaleImage(argv[1]);
        repeats = argc > 2 ? std::stoi(argv[2]) : repeats;
    } else {
        const int rows = argc > 1 ? std::stoi(argv[1]) : 256;
        const int cols = argc > 2 ? std::stoi(argv[2]) : rows;
        repeats = argc > 3 ? std::stoi(argv[3]) : repeats;
        image = makeBenchmarkImage(rows, cols);
    }

    std::cout << "domain=" << image->getNumRows() << "x" << image->getNumCols()
              << " repeats=" << repeats << '\n';
    printBenchmarkMethods();

    auto [componentTree, componentBuildMs] = timed([&]() {
        return MorphologicalTree::createComponentTree(image, true, 1.5);
    });
    runCase("component-tree max-tree radius=1.5", componentTree, componentBuildMs, repeats);

    auto [tosTree, tosBuildMs] = timed([&]() {
        return MorphologicalTree::createTreeOfShapes(image, ToSInterpolation::SelfDual);
    });
    runCase("tree-of-shapes SelfDual", tosTree, tosBuildMs, repeats);

    return 0;
}
