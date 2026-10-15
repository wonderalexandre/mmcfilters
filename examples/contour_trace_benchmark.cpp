/**
 * Benchmark sequential and node-local ordered contour tracing.
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

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

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

template <typename Trace> std::uint64_t traceChecksum(const Trace& trace) {
    std::uint64_t sum = 0;
    for (const ContourBoundary& boundary : trace.boundaries()) {
        sum += static_cast<std::uint64_t>(boundary.edgeCount + 1);
        sum += static_cast<std::uint64_t>(boundary.kind == ContourBoundaryKind::External ? 1 : 2);
        sum += static_cast<std::uint64_t>(boundary.doubledSignedArea >= 0 ? boundary.doubledSignedArea : -boundary.doubledSignedArea);
        for (ContourEdge edge : trace.boundaryEdges(boundary)) {
            sum += static_cast<std::uint64_t>(ContourTraceComputation::packEdge(edge.pixel, edge.side)) + 1;
        }
    }
    return sum;
}

std::uint64_t traversalChecksum(const ContourTraceComputation& contourTraces) {
    std::uint64_t sum = 0;
    for (auto [node, trace] : contourTraces) {
        sum += static_cast<std::uint64_t>(node + 1);
        sum += traceChecksum(trace);
    }
    return sum;
}

std::uint64_t callbackChecksum(const ContourTraceComputation& contourTraces) {
    std::uint64_t sum = 0;
    contourTraces.forEachTrace([&](NodeId node, ContourTraceView trace) {
        sum += static_cast<std::uint64_t>(node + 1);
        sum += traceChecksum(trace);
    });
    return sum;
}

void measureCase(const std::string& label, const MorphologicalTree& tree, int repeats) {
    std::cout << "\n"
              << label << '\n'
              << "  nodes: live=" << tree.numNodes() << " slots=" << tree.numInternalNodeSlots() << " proper_parts=" << tree.numPixels()
              << '\n';

    auto [contourTraces, constructionMs] = timed([&]() { return ContourTraceComputation(tree); });
    std::cout << "  construction only: " << constructionMs << " ms\n";

    auto [rootTrace, rootTraceMs] = timed([&]() { return contourTraces.trace(tree.root()); });
    std::cout << "  root trace(): " << rootTraceMs << " ms checksum=" << traceChecksum(rootTrace) << '\n';

    auto [firstTraversalChecksum, firstTraversalMs] = timed([&]() { return traversalChecksum(contourTraces); });
    const std::uint64_t firstCallbackChecksum = callbackChecksum(contourTraces);
    if (firstCallbackChecksum != firstTraversalChecksum) {
        throw std::logic_error("Callback and iterator trace checksums differ.");
    }
    std::cout << "  first complete traversal: " << firstTraversalMs << " ms checksum=" << firstTraversalChecksum << '\n';

    double totalConstructionMs = 0.0;
    double totalTraversalMs = 0.0;
    double totalCallbackMs = 0.0;
    for (int repeat = 0; repeat < repeats; ++repeat) {
        auto [runTraces, runConstructionMs] = timed([&]() { return ContourTraceComputation(tree); });
        totalConstructionMs += runConstructionMs;
        totalTraversalMs += timedVoid([&]() {
            volatile std::uint64_t checksum = traversalChecksum(runTraces);
            (void)checksum;
        });
        totalCallbackMs += timedVoid([&]() {
            volatile std::uint64_t checksum = callbackChecksum(runTraces);
            (void)checksum;
        });
    }
    std::cout << "  repeated traversal: construction_avg=" << totalConstructionMs / repeats
              << " ms iterator_avg=" << totalTraversalMs / repeats
              << " ms callback_avg=" << totalCallbackMs / repeats << " ms\n";
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
