/**
 * Benchmark contour-edge construction, lazy edge caching, and ordered
 * boundary tracing.
 *
 * Build with `-DMMCFILTERS_BUILD_EXAMPLES=ON` and run either
 * `./build/examples/mmcfilters_contour_trace_benchmark 512 512 5` for a
 * synthetic image or
 * `./build/examples/mmcfilters_contour_trace_benchmark path/to/image.png 5`
 * for a grayscale image file.
 */
#include "mmcfilters/contours/ContourTraceComputation.hpp"
#include "mmcfilters/trees/MorphologicalTree.hpp"
#include "mmcfilters/trees/MorphologicalTreeFactory.hpp"
#include "stb_image.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace mmcfilters;

namespace {

ImageUInt8Ptr makeBenchmarkImage(int rows, int columns) {
    auto image = ImageUInt8::create(rows, columns);
    for (int row = 0; row < rows; ++row) {
        for (int column = 0; column < columns; ++column) {
            const PixelId pixel = row * columns + column;
            const int radial = (row - rows / 2) * (row - rows / 2) + (column - columns / 2) * (column - columns / 2);
            (*image)[pixel] = static_cast<uint8_t>((radial / 113 + row * 7 + column * 3) & 0xff);
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

template <typename Fn> auto timed(Fn&& fn) {
    const auto start = std::chrono::steady_clock::now();
    auto result = fn();
    const auto end = std::chrono::steady_clock::now();
    const auto micros = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    return std::pair{std::move(result), static_cast<double>(micros) / 1000.0};
}

template <typename Fn> double timedVoid(Fn&& fn) {
    const auto start = std::chrono::steady_clock::now();
    fn();
    const auto end = std::chrono::steady_clock::now();
    const auto micros = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    return static_cast<double>(micros) / 1000.0;
}

std::vector<NodeId> aliveNodesInTreeOrder(const MorphologicalTree& tree) {
    std::vector<NodeId> nodes;
    nodes.reserve(static_cast<std::size_t>(tree.numNodes()));
    for (NodeId node : tree.aliveNodeIds()) {
        nodes.push_back(node);
    }
    return nodes;
}

std::vector<NodeId> shuffledAliveNodes(const MorphologicalTree& tree) {
    auto nodes = aliveNodesInTreeOrder(tree);
    std::mt19937 rng(0x71ace);
    std::shuffle(nodes.begin(), nodes.end(), rng);
    return nodes;
}

std::size_t edgeChecksum(const ContourTraceComputation& contourTraces, NodeId node) {
    std::size_t sum = 0;
    for (ContourEdge edge : contourTraces.edges(node)) {
        sum += static_cast<std::size_t>(ContourTraceComputation::packEdge(edge.pixel, edge.side) + 1);
    }
    return sum;
}

std::size_t boundaryChecksum(const ContourTraceComputation& contourTraces, NodeId node) {
    std::size_t sum = 0;
    for (const ContourBoundary& boundary : contourTraces.boundaries(node)) {
        sum += static_cast<std::size_t>(boundary.edgeCount + 1);
        sum += static_cast<std::size_t>(boundary.doubledSignedArea >= 0 ? boundary.doubledSignedArea : -boundary.doubledSignedArea);
        for (ContourEdge edge : contourTraces.boundaryEdges(boundary)) {
            sum += static_cast<std::size_t>(ContourTraceComputation::packEdge(edge.pixel, edge.side) + 1);
        }
    }
    return sum;
}

std::size_t edgeChecksumInOrder(const std::vector<NodeId>& nodes, const ContourTraceComputation& contourTraces) {
    std::size_t sum = 0;
    for (NodeId node : nodes) {
        sum += edgeChecksum(contourTraces, node);
    }
    return sum;
}

std::size_t boundaryChecksumInOrder(const std::vector<NodeId>& nodes, const ContourTraceComputation& contourTraces) {
    std::size_t sum = 0;
    for (NodeId node : nodes) {
        sum += boundaryChecksum(contourTraces, node);
    }
    return sum;
}

void measureCase(const std::string& label, const MorphologicalTree& tree, int repeats) {
    const NodeId root = tree.root();
    const auto treeOrder = aliveNodesInTreeOrder(tree);
    const auto randomOrder = shuffledAliveNodes(tree);

    std::cout << "\n"
              << label << '\n'
              << "  nodes: live=" << tree.numNodes() << " slots=" << tree.numInternalNodeSlots() << " proper_parts=" << tree.numPixels()
              << '\n';

    auto [contourTraces, constructionMs] = timed([&]() { return ContourTraceComputation(tree); });
    (void)contourTraces;
    std::cout << "  construction only: " << constructionMs << " ms";
    std::cout << '\n';

    auto [rootEdgeComputation, rootEdgeConstructionMs] = timed([&]() { return ContourTraceComputation(tree); });
    const double rootEdgeAccessMs = timedVoid([&]() {
        volatile std::size_t checksum = edgeChecksum(rootEdgeComputation, root);
        (void)checksum;
    });
    std::cout << "  root edges(): construction=" << rootEdgeConstructionMs << " ms edge_access=" << rootEdgeAccessMs;
    std::cout << '\n';

    auto [rootBoundaryComputation, rootBoundaryConstructionMs] = timed([&]() { return ContourTraceComputation(tree); });
    const double rootBoundaryTraceMs = timedVoid([&]() {
        volatile std::size_t checksum = boundaryChecksum(rootBoundaryComputation, root);
        (void)checksum;
    });
    std::cout << "  root boundaries(): construction=" << rootBoundaryConstructionMs << " ms trace=" << rootBoundaryTraceMs;
    std::cout << '\n';

    auto [orderedEdgeComputation, orderedEdgeConstructionMs] = timed([&]() { return ContourTraceComputation(tree); });
    const double orderedEdgeAccessMs = timedVoid([&]() {
        volatile std::size_t checksum = edgeChecksumInOrder(treeOrder, orderedEdgeComputation);
        (void)checksum;
    });
    std::cout << "  all edges() tree order: construction=" << orderedEdgeConstructionMs << " ms edge_access=" << orderedEdgeAccessMs;
    std::cout << '\n';

    auto [randomBoundaryComputation, randomBoundaryConstructionMs] = timed([&]() { return ContourTraceComputation(tree); });
    const double randomBoundaryTraceMs = timedVoid([&]() {
        volatile std::size_t checksum = boundaryChecksumInOrder(randomOrder, randomBoundaryComputation);
        (void)checksum;
    });
    std::cout << "  all boundaries() random order: construction=" << randomBoundaryConstructionMs << " ms trace=" << randomBoundaryTraceMs;
    std::cout << '\n';

    double totalConstructionMs = 0.0;
    double totalTraceMs = 0.0;
    for (int i = 0; i < repeats; ++i) {
        auto [contourTraces, runConstructionMs] = timed([&]() { return ContourTraceComputation(tree); });
        totalConstructionMs += runConstructionMs;
        totalTraceMs += timedVoid([&]() { contourTraces.traceAll(); });
    }
    std::cout << "  repeated traceAll(): construction_avg=" << totalConstructionMs / repeats << " ms trace_avg=" << totalTraceMs / repeats << " ms";
    std::cout << '\n';
}

void runBenchmark(ImageUInt8Ptr image, int repeats) {
    auto [maxTree, maxBuildMs] = timed([&]() { return MorphologicalTreeFactory::createMaxTree(image, 1.5); });
    std::cout << "max-tree build: " << maxBuildMs << " ms\n";
    measureCase("max-tree contour traces", maxTree.topology(), repeats);

    auto [minTree, minBuildMs] = timed([&]() { return MorphologicalTreeFactory::createMinTree(image, 1.5); });
    std::cout << "\nmin-tree build: " << minBuildMs << " ms\n";
    measureCase("min-tree contour traces", minTree.topology(), repeats);
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 4) {
            const int rows = std::stoi(argv[1]);
            const int columns = std::stoi(argv[2]);
            const int repeats = std::stoi(argv[3]);
            if (rows <= 0 || columns <= 0 || repeats <= 0) {
                throw std::invalid_argument("rows, columns, and repeats must be positive.");
            }
            runBenchmark(makeBenchmarkImage(rows, columns), repeats);
            return 0;
        }

        if (argc == 3) {
            const int repeats = std::stoi(argv[2]);
            if (repeats <= 0) {
                throw std::invalid_argument("repeats must be positive.");
            }
            runBenchmark(loadGrayscaleImage(argv[1]), repeats);
            return 0;
        }

        std::cerr << "usage:\n"
                  << "  " << argv[0] << " rows columns repeats\n"
                  << "  " << argv[0] << " image-path repeats\n";
        return 2;
    } catch (const std::exception& ex) {
        std::cerr << "error: " << ex.what() << '\n';
        return 1;
    }
}
