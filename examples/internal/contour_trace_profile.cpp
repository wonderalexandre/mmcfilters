/**
 * External timing harness for node-local and sequential boundary tracing.
 * Build with optimization and debug symbols to use a sampling profiler.
 * All measurements stay in this executable; the library has no profiling hooks.
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

void profileCase(const std::string& label, const MorphologicalTree& tree, int repeats) {
    double totalConstructionMs = 0.0;
    double totalRootTraceMs = 0.0;
    double totalTraversalMs = 0.0;

    for (int repeat = 0; repeat < repeats; ++repeat) {
        auto [contourTraces, constructionMs] = timed([&]() { return ContourTraceComputation(tree); });
        totalConstructionMs += constructionMs;

        totalRootTraceMs += timedVoid([&]() {
            const auto trace = contourTraces.trace(tree.root());
            volatile std::size_t edgeCount = trace.edges().size();
            (void)edgeCount;
        });

        totalTraversalMs += timedVoid([&]() {
            std::size_t edgeCount = 0;
            contourTraces.forEachTrace([&](NodeId, ContourTraceView trace) { edgeCount += trace.edges().size(); });
            volatile std::size_t consumedEdgeCount = edgeCount;
            (void)consumedEdgeCount;
        });
    }

    std::cout << "\n"
              << label << '\n'
              << "  nodes: live=" << tree.numNodes() << " slots=" << tree.numInternalNodeSlots() << " proper_parts=" << tree.numPixels()
              << '\n'
              << "  construction_avg=" << totalConstructionMs / repeats << " ms"
              << " root_trace_avg=" << totalRootTraceMs / repeats << " ms"
              << " traversal_avg=" << totalTraversalMs / repeats << " ms\n";
}

void runProfile(ImageUInt8Ptr image, int repeats) {
    auto [maxTree, maxBuildMs] = timed([&]() { return MorphologicalTreeFactory::createMaxTree(image, 1.5); });
    std::cout << "max-tree build: " << maxBuildMs << " ms\n";
    profileCase("max-tree contour trace profile", maxTree.topology(), repeats);

    auto [minTree, minBuildMs] = timed([&]() { return MorphologicalTreeFactory::createMinTree(image, 1.5); });
    std::cout << "\nmin-tree build: " << minBuildMs << " ms\n";
    profileCase("min-tree contour trace profile", minTree.topology(), repeats);
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
            runProfile(makeBenchmarkImage(rows, columns), repeats);
            return 0;
        }

        if (argc == 3) {
            const int repeats = std::stoi(argv[2]);
            if (repeats <= 0) {
                throw std::invalid_argument("repeats must be positive.");
            }
            runProfile(loadGrayscaleImage(argv[1]), repeats);
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
