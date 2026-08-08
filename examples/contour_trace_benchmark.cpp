/**
 * Benchmark side-level contour trace extraction, lazy edge materialization, and
 * loop tracing.
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

ImageUInt8Ptr makeBenchmarkImage(int rows, int cols) {
    auto image = ImageUInt8::create(rows, cols);
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            const int idx = row * cols + col;
            const int radial = (row - rows / 2) * (row - rows / 2) + (col - cols / 2) * (col - cols / 2);
            (*image)[idx] = static_cast<uint8_t>((radial / 113 + row * 7 + col * 3) & 0xff);
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
    nodes.reserve(static_cast<std::size_t>(tree.getNumNodes()));
    for (NodeId node : tree.getAliveNodeIds()) {
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

std::size_t sumEdges(const ContourTraceComputation::IncrementalContourTraces& traces, NodeId node) {
    std::size_t sum = 0;
    for (ContourTraceEdge edge : traces.getEdges(node)) {
        sum += static_cast<std::size_t>(ContourTraceComputation::packEdge(edge.pixel, edge.side) + 1);
    }
    return sum;
}

std::size_t sumLoops(const ContourTraceComputation::IncrementalContourTraces& traces, NodeId node) {
    std::size_t sum = 0;
    for (const ContourTraceLoop& loop : traces.getLoops(node)) {
        sum += static_cast<std::size_t>(loop.edgeCount + 1);
        sum += static_cast<std::size_t>(loop.signedArea2 >= 0 ? loop.signedArea2 : -loop.signedArea2);
        for (ContourTraceEdge edge : traces.getLoopEdges(loop)) {
            sum += static_cast<std::size_t>(ContourTraceComputation::packEdge(edge.pixel, edge.side) + 1);
        }
    }
    return sum;
}

std::size_t sumEdgesInOrder(const std::vector<NodeId>& nodes, const ContourTraceComputation::IncrementalContourTraces& traces) {
    std::size_t sum = 0;
    for (NodeId node : nodes) {
        sum += sumEdges(traces, node);
    }
    return sum;
}

std::size_t sumLoopsInOrder(const std::vector<NodeId>& nodes, const ContourTraceComputation::IncrementalContourTraces& traces) {
    std::size_t sum = 0;
    for (NodeId node : nodes) {
        sum += sumLoops(traces, node);
    }
    return sum;
}

#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
double mib(std::size_t bytes) { return static_cast<double>(bytes) / (1024.0 * 1024.0); }

std::size_t storageChecksum(const ContourTraceComputation::IncrementalContourTraces::StorageStats& stats) {
    return stats.addDeltaValues + stats.removeDeltaValues + stats.cachedEdgeValues + stats.cachedLoopEdges + stats.cachedLoops + stats.cachedEdgeReadyNodes +
           stats.cachedLoopReadyNodes + stats.approxAllocatedBytes;
}
#endif

void appendStorageChecksum(const ContourTraceComputation::IncrementalContourTraces& traces) {
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
    std::cout << " checksum=" << storageChecksum(traces.storageStats());
#else
    (void)traces;
#endif
}

#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
void printStorageStats(const std::string& label, const ContourTraceComputation::IncrementalContourTraces& traces) {
    const auto stats = traces.storageStats();
    std::cout << "    " << label << ": "
              << "add_delta_values=" << stats.addDeltaValues << " remove_delta_values=" << stats.removeDeltaValues
              << " cached_edge_values=" << stats.cachedEdgeValues << " cached_loop_edges=" << stats.cachedLoopEdges << " cached_loops=" << stats.cachedLoops
              << " cached_edge_ready_nodes=" << stats.cachedEdgeReadyNodes << " cached_loop_ready_nodes=" << stats.cachedLoopReadyNodes
              << " dense_outgoing_slots=" << stats.traceDenseOutgoingSlots << " sparse_outgoing_slots=" << stats.traceSparseOutgoingSlots
              << " approx_allocated=" << mib(stats.approxAllocatedBytes) << " MiB\n";
}
#endif

void measureCase(const std::string& label, const MorphologicalTree& tree, int repeats) {
    const NodeId root = tree.getRoot();
    const auto treeOrder = aliveNodesInTreeOrder(tree);
    const auto randomOrder = shuffledAliveNodes(tree);

    std::cout << "\n"
              << label << '\n'
              << "  nodes: live=" << tree.getNumNodes() << " slots=" << tree.getNumInternalNodeSlots() << " proper_parts=" << tree.getNumTotalProperParts()
              << '\n';

    auto [extractOnly, extractMs] = timed([&]() { return ContourTraceComputation::extract(tree); });
    std::cout << "  extract only: " << extractMs << " ms";
    appendStorageChecksum(extractOnly);
    std::cout << '\n';
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
    printStorageStats("after extract", extractOnly);
#endif

    auto [rootEdgeTraces, rootEdgeExtractMs] = timed([&]() { return ContourTraceComputation::extract(tree); });
    const double rootEdgeMs = timedVoid([&]() {
        volatile std::size_t checksum = sumEdges(rootEdgeTraces, root);
        (void)checksum;
    });
    std::cout << "  root getEdges(): extract=" << rootEdgeExtractMs << " ms materialize=" << rootEdgeMs;
    appendStorageChecksum(rootEdgeTraces);
    std::cout << '\n';
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
    printStorageStats("after root edges", rootEdgeTraces);
#endif

    auto [rootLoopTraces, rootLoopExtractMs] = timed([&]() { return ContourTraceComputation::extract(tree); });
    const double rootLoopMs = timedVoid([&]() {
        volatile std::size_t checksum = sumLoops(rootLoopTraces, root);
        (void)checksum;
    });
    std::cout << "  root getLoops(): extract=" << rootLoopExtractMs << " ms trace=" << rootLoopMs;
    appendStorageChecksum(rootLoopTraces);
    std::cout << '\n';
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
    printStorageStats("after root loops", rootLoopTraces);
#endif

    auto [orderedEdgeTraces, orderedEdgeExtractMs] = timed([&]() { return ContourTraceComputation::extract(tree); });
    const double orderedEdgeMs = timedVoid([&]() {
        volatile std::size_t checksum = sumEdgesInOrder(treeOrder, orderedEdgeTraces);
        (void)checksum;
    });
    std::cout << "  all getEdges() tree order: extract=" << orderedEdgeExtractMs << " ms materialize=" << orderedEdgeMs;
    appendStorageChecksum(orderedEdgeTraces);
    std::cout << '\n';
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
    printStorageStats("after all edges", orderedEdgeTraces);
#endif

    auto [randomLoopTraces, randomLoopExtractMs] = timed([&]() { return ContourTraceComputation::extract(tree); });
    const double randomLoopMs = timedVoid([&]() {
        volatile std::size_t checksum = sumLoopsInOrder(randomOrder, randomLoopTraces);
        (void)checksum;
    });
    std::cout << "  all getLoops() random order: extract=" << randomLoopExtractMs << " ms trace=" << randomLoopMs;
    appendStorageChecksum(randomLoopTraces);
    std::cout << '\n';
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
    printStorageStats("after all loops random", randomLoopTraces);
#endif

    double totalExtract = 0.0;
    double totalMaterialize = 0.0;
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
    std::size_t repeatedChecksum = 0;
#endif
    for (int i = 0; i < repeats; ++i) {
        auto [traces, runExtractMs] = timed([&]() { return ContourTraceComputation::extract(tree); });
        totalExtract += runExtractMs;
        totalMaterialize += timedVoid([&]() { traces.materializeAll(); });
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
        repeatedChecksum += storageChecksum(traces.storageStats());
#endif
    }
    std::cout << "  repeated materializeAll(): extract_avg=" << totalExtract / repeats << " ms materialize_avg=" << totalMaterialize / repeats << " ms";
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
    std::cout << " checksum=" << repeatedChecksum;
#endif
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
            const int cols = std::stoi(argv[2]);
            const int repeats = std::stoi(argv[3]);
            if (rows <= 0 || cols <= 0 || repeats <= 0) {
                throw std::invalid_argument("rows, cols, and repeats must be positive.");
            }
            runBenchmark(makeBenchmarkImage(rows, cols), repeats);
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
                  << "  " << argv[0] << " rows cols repeats\n"
                  << "  " << argv[0] << " image-path repeats\n";
        return 2;
    } catch (const std::exception& ex) {
        std::cerr << "error: " << ex.what() << '\n';
        return 1;
    }
}
