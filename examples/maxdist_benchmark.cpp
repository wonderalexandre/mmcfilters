/**
 * Benchmark `MAX_DIST` attribute computation on max-trees and min-trees.
 *
 * Build with `-DMMCFILTERS_BUILD_EXAMPLES=ON` and run:
 * `./build/examples/mmcfilters_maxdist_benchmark --sizes 128,256,512
 * --repeats 3 --radius 1.5`.
 *
 * `MAX_DIST` requires component trees with adjacency metadata. The checksum
 * keeps benchmark runs observable without printing per-node buffers.
 */
#include "mmcfilters/attributes/Attributes.hpp"
#include "mmcfilters/trees/MorphologicalTreeFactory.hpp"
#include "mmcfilters/trees/WeightedMorphologicalTree.hpp"

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace mmcfilters;

namespace {

struct Options {
    std::vector<int> sizes{128, 256, 512};
    int repeats = 3;
    double radius = 1.5;
};

struct Measurement {
    double msPerRun = 0.0;
    std::uint64_t checksum = 0;
};

ImageUInt8Ptr makeBenchmarkImage(int rows, int cols) {
    auto image = ImageUInt8::create(rows, cols);
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            const int idx = row * cols + col;
            const int radial = (row - rows / 2) * (row - rows / 2) +
                               (col - cols / 2) * (col - cols / 2);
            const int waves = (row * 17) ^ (col * 31) ^ ((row + col) * 7);
            (*image)[idx] = static_cast<uint8_t>((radial / 113 + waves) & 0xff);
        }
    }
    return image;
}

template <typename Fn>
auto timedValue(Fn&& fn) {
    const auto start = std::chrono::steady_clock::now();
    auto value = fn();
    const auto end = std::chrono::steady_clock::now();
    const auto micros = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    return std::pair{std::move(value), static_cast<double>(micros) / 1000.0};
}

std::uint64_t checksumMaxDist(const WeightedMorphologicalTree<std::uint8_t>& tree, const ComputedAttributeData<float>& data) {
    const auto& attrNames = data.first;
    const auto& buffer = data.second;

    std::uint64_t checksum = 0;
    for (NodeId nodeId : tree.topology().getAliveNodeIds()) {
        const float value = buffer[attrNames.linearIndex(nodeId, MAX_DIST)];
        checksum += static_cast<std::uint64_t>(value * 1000.0f) + static_cast<std::uint64_t>(nodeId + 1);
    }
    return checksum;
}

Measurement measureMaxDist(const WeightedMorphologicalTree<std::uint8_t>& tree, int repeats) {
    if (repeats <= 0) {
        throw std::invalid_argument("repeats must be positive.");
    }

    std::uint64_t checksum = 0;
    const auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < repeats; ++i) {
        auto data = AttributeComputation::computeSingleAttribute(tree, MAX_DIST);
        checksum += checksumMaxDist(tree, data);
    }
    const auto end = std::chrono::steady_clock::now();
    const auto micros = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    return {static_cast<double>(micros) / static_cast<double>(repeats) / 1000.0, checksum};
}

std::vector<int> parseSizes(const std::string& value) {
    std::vector<int> sizes;
    std::stringstream stream(value);
    std::string item;
    while (std::getline(stream, item, ',')) {
        const int size = std::stoi(item);
        if (size <= 0) {
            throw std::invalid_argument("--sizes values must be positive.");
        }
        sizes.push_back(size);
    }
    if (sizes.empty()) {
        throw std::invalid_argument("--sizes requires at least one value.");
    }
    return sizes;
}

Options parseOptions(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--sizes") {
            if (++i >= argc) {
                throw std::invalid_argument("--sizes requires a comma-separated value.");
            }
            options.sizes = parseSizes(argv[i]);
        } else if (arg == "--repeats") {
            if (++i >= argc) {
                throw std::invalid_argument("--repeats requires a value.");
            }
            options.repeats = std::stoi(argv[i]);
            if (options.repeats <= 0) {
                throw std::invalid_argument("--repeats must be positive.");
            }
        } else if (arg == "--radius") {
            if (++i >= argc) {
                throw std::invalid_argument("--radius requires a value.");
            }
            options.radius = std::stod(argv[i]);
            if (options.radius <= 0.0) {
                throw std::invalid_argument("--radius must be positive.");
            }
        } else if (arg == "--help" || arg == "-h") {
            std::cout
                << "usage: mmcfilters_maxdist_benchmark [--sizes 128,256,512] [--repeats 3] [--radius 1.5]\n";
            std::exit(0);
        } else {
            throw std::invalid_argument("unknown argument: " + arg);
        }
    }
    return options;
}

void runCase(int rows, int cols, bool isMaxTree, int repeats, double radius) {
    auto image = makeBenchmarkImage(rows, cols);
    auto [tree, buildMs] = timedValue([&]() {
        return isMaxTree
            ? MorphologicalTreeFactory::createMaxTree(image, radius)
            : MorphologicalTreeFactory::createMinTree(image, radius);
    });
    const Measurement maxDist = measureMaxDist(tree, repeats);

    std::cout << std::setw(5) << rows
              << " x " << std::setw(5) << cols
              << "  " << (isMaxTree ? "max" : "min")
              << "  nodes=" << std::setw(8) << tree.topology().getNumNodes()
              << "  build_ms=" << std::setw(10) << std::fixed << std::setprecision(3) << buildMs
              << "  max_dist_ms/run=" << std::setw(10) << std::fixed << std::setprecision(3) << maxDist.msPerRun
              << "  checksum=" << maxDist.checksum
              << '\n';
}

} // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);

        std::cout << "MAX_DIST benchmark\n"
                  << "repeats=" << options.repeats
                  << " radius=" << options.radius
                  << "\n\n";

        for (int size : options.sizes) {
            runCase(size, size, true, options.repeats, options.radius);
            runCase(size, size, false, options.repeats, options.radius);
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
}
