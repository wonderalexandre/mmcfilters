/**
 * Benchmark compact contour indexes, incremental iteration, and ordered
 * boundary tracing.
 *
 * Build with `-DMMCFILTERS_BUILD_EXAMPLES=ON` and run either
 * `./build/examples/mmcfilters_contour_benchmark 512 512 5` for a synthetic
 * image or `./build/examples/mmcfilters_contour_benchmark path/to/image.png 5`
 * for a grayscale image file.
 *
 * Timings separate contour construction, access to one node, iterator
 * traversal, callbacks, contour-trace construction, and ordered-boundary
 * traversal. Measure process memory externally; the public API has no profiling
 * counters.
 */
#include "mmcfilters/contours/ContourComputation.hpp"
#include "mmcfilters/contours/ContourTraceComputation.hpp"
#include "mmcfilters/trees/MorphologicalTree.hpp"
#include "mmcfilters/trees/MorphologicalTreeFactory.hpp"
#include "stb_image.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <limits>
#include <span>
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
            (*image)[pixel] = static_cast<uint8_t>((radial / 97 + row * 3 + column * 5) & 0xff);
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

template <typename Fn> void measure(const std::string& label, int repeats, Fn&& fn) {
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

std::size_t contourChecksum(const ContourComputation& contours, NodeId node) {
    std::size_t sum = 0;
    for (PixelId pixel : contours.contour(node)) {
        sum += static_cast<std::size_t>(pixel + 1);
    }
    return sum;
}

std::size_t allContoursChecksum(const ContourComputation& contours) {
    std::size_t sum = 0;
    for (auto [node, pixels] : contours) {
        (void)node;
        for (PixelId pixel : pixels) {
            sum += static_cast<std::size_t>(pixel + 1);
        }
    }
    return sum;
}

std::size_t callbackContoursChecksum(const ContourComputation& contours) {
    std::size_t sum = 0;
    contours.forEachContour([&](NodeId, std::span<const PixelId> pixels) {
        for (PixelId pixel : pixels) {
            sum += static_cast<std::size_t>(pixel + 1);
        }
    });
    return sum;
}

std::uint64_t allBoundariesChecksum(const MorphologicalTree& tree, const ContourTraceComputation& contourTraces) {
    std::uint64_t sum = 0;
    for (NodeId node : tree.aliveNodeIds()) {
        sum += static_cast<std::uint64_t>(node + 1);
        for (const ContourBoundary& boundary : contourTraces.boundaries(node)) {
            sum += static_cast<std::uint64_t>(boundary.edgeCount + 1);
            sum += static_cast<std::uint64_t>(boundary.kind == ContourBoundaryKind::External ? 1 : 2);
            sum += static_cast<std::uint64_t>(boundary.doubledSignedArea >= 0 ? boundary.doubledSignedArea : -boundary.doubledSignedArea);
            for (ContourEdge edge : contourTraces.boundaryEdges(boundary)) {
                sum += static_cast<std::uint64_t>(ContourTraceComputation::packEdge(edge.pixel, edge.side)) + 1;
            }
        }
    }
    return sum;
}

void printTreeStats(const MorphologicalTree& tree, double buildMs) {
    int maxChildren = 0;
    long long totalChildren = 0;
    long long totalDirectProperParts = 0;
    int maxDirectProperParts = 0;
    int minDirectProperParts = std::numeric_limits<int>::max();

    for (NodeId node : tree.aliveNodeIds()) {
        const int children = tree.numChildren(node);
        const int properParts = tree.properPartCardinality(node);
        totalChildren += children;
        totalDirectProperParts += properParts;
        maxChildren = std::max(maxChildren, children);
        maxDirectProperParts = std::max(maxDirectProperParts, properParts);
        minDirectProperParts = std::min(minDirectProperParts, properParts);
    }

    const int numNodes = tree.numNodes();
    if (numNodes == 0) {
        minDirectProperParts = 0;
    }

    std::cout << "  build: " << buildMs << " ms\n"
              << "  nodes: live=" << numNodes << " slots=" << tree.numInternalNodeSlots() << " leaves=" << tree.numLeafNodes()
              << " proper_parts=" << tree.numPixels() << '\n'
              << "  branching: avg_children=" << (numNodes > 0 ? static_cast<double>(totalChildren) / numNodes : 0.0) << " max_children=" << maxChildren << '\n'
              << "  direct_proper_parts: avg=" << (numNodes > 0 ? static_cast<double>(totalDirectProperParts) / numNodes : 0.0)
              << " min=" << minDirectProperParts << " max=" << maxDirectProperParts << '\n';
}

template <typename TreeView> void runCase(const std::string& label, const TreeView& view, double buildMs, int repeats) {
    const MorphologicalTree& tree = view.topology();
    std::cout << "\n[" << label << "]\n";
    printTreeStats(tree, buildMs);

    auto [contours, contourConstructionMs] = timed([&]() { return ContourComputation(view); });
    std::cout << "  contour construction (includes batch LCA when needed): " << contourConstructionMs << " ms\n";
    const auto iteratorChecksum = allContoursChecksum(contours);
    if (callbackContoursChecksum(contours) != iteratorChecksum) {
        throw std::logic_error("Callback and iterator contour checksums differ.");
    }
    measure("construct + root contour()", repeats, [&]() { return contourChecksum(ContourComputation(view), tree.root()); });
    measure("construct + iterate all contours", repeats, [&]() { return allContoursChecksum(ContourComputation(view)); });
    measure("construct + callback for all contours", repeats, [&]() { return callbackContoursChecksum(ContourComputation(view)); });
    measure("root contour() on existing computation", repeats, [&]() { return contourChecksum(contours, tree.root()); });
    measure("repeat traversal on existing computation", repeats, [&]() { return allContoursChecksum(contours); });
    measure("repeat callback on existing computation", repeats, [&]() { return callbackContoursChecksum(contours); });

    auto [contourTraces, traceConstructionMs] = timed([&]() { return ContourTraceComputation(view); });
    auto [boundaryChecksum, firstBoundaryTraversalMs] = timed([&]() { return allBoundariesChecksum(tree, contourTraces); });
    std::cout << "  contour-trace construction: " << traceConstructionMs << " ms\n"
              << "  first ordered-boundary trace: " << firstBoundaryTraversalMs << " ms checksum=" << boundaryChecksum << '\n';

    measure("construct + trace all boundaries", repeats, [&]() { return allBoundariesChecksum(tree, ContourTraceComputation(view)); });
    measure("repeat boundary traversal on existing computation", repeats, [&]() { return allBoundariesChecksum(tree, contourTraces); });
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
        const int columns = argc > 2 ? std::stoi(argv[2]) : rows;
        repeats = argc > 3 ? std::stoi(argv[3]) : repeats;
        image = makeBenchmarkImage(rows, columns);
    }

    if (repeats <= 0) {
        throw std::invalid_argument("Repeat count must be positive.");
    }

    std::cout << "domain=" << image->getNumRows() << "x" << image->getNumColumns() << " repeats=" << repeats << '\n';

    auto [componentTree, componentBuildMs] = timed([&]() { return MorphologicalTreeFactory::createMaxTree(image, 1.5); });
    runCase("component-tree max-tree radius=1.5", componentTree.asView(), componentBuildMs, repeats);

    auto [minTree, minBuildMs] = timed([&]() { return MorphologicalTreeFactory::createMinTree(image, 1.5); });
    runCase("component-tree min-tree radius=1.5", minTree.asView(), minBuildMs, repeats);

    auto [tosTree, tosBuildMs] = timed([&]() { return MorphologicalTreeFactory::createTreeOfShapes<ToSGrayLevel>(image, selfDualSpanConvention()); });
    runCase("tree of shapes (self-dual span)", tosTree.asView(), tosBuildMs, repeats);

    return 0;
}
