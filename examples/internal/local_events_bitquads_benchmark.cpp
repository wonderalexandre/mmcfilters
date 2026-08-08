/**
 * Internal benchmark for local-event bitquad and contour-side computations.
 *
 * Build with `-DMMCFILTERS_BUILD_EXAMPLES=ON` and run:
 * `./build/examples/mmcfilters_local_events_bitquads_benchmark 256 256 3
 * --tos`.
 *
 * The benchmark compares local-event storage and execution paths with public
 * attribute requests. It includes `detail/` headers and is not public API.
 */
#include "mmcfilters/attributes/Attributes.hpp"
#include "mmcfilters/attributes/computers/AttributeComputerRegistry.hpp"
#include "mmcfilters/attributes/computers/BitquadAttributeComputer.hpp"
#include "mmcfilters/attributes/computers/detail/BitquadLocalEventComputation.hpp"
#include "mmcfilters/attributes/computers/detail/ContourSideLocalEventComputation.hpp"
#include "mmcfilters/trees/MorphologicalTree.hpp"
#include "mmcfilters/trees/MorphologicalTreeFactory.hpp"
#include "mmcfilters/trees/WeightedMorphologicalTree.hpp"
#include "stb_image.h"

#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cctype>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using namespace mmcfilters;
using namespace mmcfilters::attributes::computers;
using namespace mmcfilters::attributes::computers::detail;

namespace {

struct Options {
    int rows = 128;
    int cols = 128;
    int repeats = 3;
    bool includeTos = false;
    bool skipSynthetic = false;
    std::vector<ToSInterpolation> tosInterpolations;
    std::vector<std::filesystem::path> imagePaths;
};

struct Measurement {
    double msPerRun = 0.0;
    std::uint64_t checksum = 0;
};

struct TreeOfShapesPolarityStats {
    int rootNodes = 0;
    int minTreeNodes = 0;
    int maxTreeNodes = 0;
    int equalAltitudeEdges = 0;
};

std::string toString(ToSInterpolation interpolation) {
    switch (interpolation) {
    case ToSInterpolation::SelfDual:
        return "SelfDual";
    case ToSInterpolation::Min4cMax8c:
        return "Min4cMax8c";
    case ToSInterpolation::Min8cMax4c:
        return "Min8cMax4c";
    }
    throw std::invalid_argument("Unsupported ToS interpolation.");
}

ToSInterpolation parseToSInterpolation(const std::string& value) {
    if (value == "SelfDual") {
        return ToSInterpolation::SelfDual;
    }
    if (value == "Min4cMax8c") {
        return ToSInterpolation::Min4cMax8c;
    }
    if (value == "Min8cMax4c") {
        return ToSInterpolation::Min8cMax4c;
    }
    throw std::invalid_argument("Unsupported --tos-interpolation value: " + value + ". Expected SelfDual, Min4cMax8c, or Min8cMax4c.");
}

std::string lowerExtension(const std::filesystem::path& path) {
    std::string extension = path.extension().string();
    std::transform(extension.begin(), extension.end(), extension.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return extension;
}

bool hasSupportedImageExtension(const std::filesystem::path& path) {
    const std::string extension = lowerExtension(path);
    return extension == ".png" || extension == ".jpg" || extension == ".jpeg" || extension == ".bmp" || extension == ".pgm" || extension == ".ppm" ||
           extension == ".tif" || extension == ".tiff";
}

ImageUInt8Ptr loadGrayscaleImage(const std::filesystem::path& path) {
    int width = 0;
    int height = 0;
    int channels = 0;
    unsigned char* raw = stbi_load(path.string().c_str(), &width, &height, &channels, 1);
    if (raw == nullptr) {
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

std::vector<std::filesystem::path> collectImagePaths(const std::filesystem::path& directory) {
    if (!std::filesystem::is_directory(directory)) {
        throw std::invalid_argument("Image directory does not exist: " + directory.string());
    }

    std::vector<std::filesystem::path> images;
    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
        if (entry.is_regular_file() && hasSupportedImageExtension(entry.path())) {
            images.push_back(entry.path());
        }
    }
    std::sort(images.begin(), images.end());
    return images;
}

ImageUInt8Ptr makePatternImage(int rows, int cols, const std::string& pattern) {
    auto image = ImageUInt8::create(rows, cols);
    std::uint32_t rng = 0x12345678u;

    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            const int idx = row * cols + col;
            uint8_t value = 0;

            if (pattern == "flat") {
                value = 128;
            } else if (pattern == "ramp") {
                value = static_cast<uint8_t>((row * 7 + col * 11) & 0xff);
            } else if (pattern == "noise") {
                rng ^= rng << 13;
                rng ^= rng >> 17;
                rng ^= rng << 5;
                value = static_cast<uint8_t>(rng & 0xff);
            } else if (pattern == "hollow") {
                const int top = rows / 4;
                const int bottom = rows - top - 1;
                const int left = cols / 4;
                const int right = cols - left - 1;
                const bool inOuter = row >= top && row <= bottom && col >= left && col <= right;
                const bool inInner = row > top + rows / 8 && row < bottom - rows / 8 && col > left + cols / 8 && col < right - cols / 8;
                value = inOuter && !inInner ? 240 : 30;
            } else if (pattern == "diagonal") {
                value = std::abs(row - col) <= std::max(1, rows / 32) ? 220 : 40;
            } else {
                throw std::invalid_argument("Unknown pattern: " + pattern);
            }

            (*image)[idx] = value;
        }
    }

    return image;
}

template <typename Fn> auto timed(Fn&& fn) {
    const auto start = std::chrono::steady_clock::now();
    auto result = fn();
    const auto end = std::chrono::steady_clock::now();
    const auto micros = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    return std::pair{std::move(result), static_cast<double>(micros) / 1000.0};
}

template <typename Fn> Measurement measure(int repeats, Fn&& fn) {
    std::uint64_t checksum = 0;
    const auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < repeats; ++i) {
        checksum += fn();
    }
    const auto end = std::chrono::steady_clock::now();
    const auto micros = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    return {static_cast<double>(micros) / static_cast<double>(repeats) / 1000.0, checksum};
}

double mib(std::size_t bytes) { return static_cast<double>(bytes) / (1024.0 * 1024.0); }

void printMemoryEstimate(const std::string& label, std::size_t bytes, int numNodes) {
    const double bytesPerNode = numNodes > 0 ? static_cast<double>(bytes) / static_cast<double>(numNodes) : 0.0;
    std::cout << "    " << label << ": " << mib(bytes) << " MiB"
              << " (" << bytesPerNode << " B/node)\n";
}

std::uint64_t mix(std::uint64_t seed, std::uint64_t value) {
    seed ^= value + 0x9e3779b97f4a7c15ull + (seed << 6) + (seed >> 2);
    return seed;
}

std::uint64_t checksumLocalBitquads(const BitquadLocalEventComputation& computer) {
    std::uint64_t checksum = 0;
    for (const auto& histogram : computer.getBitquadStateHistograms()) {
        for (int value : histogram) {
            checksum = mix(checksum, static_cast<std::uint64_t>(value));
        }
    }
    for (const auto& family : computer.getBitquadFamilyCounts()) {
        checksum = mix(checksum, static_cast<std::uint64_t>(family.empty));
        checksum = mix(checksum, static_cast<std::uint64_t>(family.q1));
        checksum = mix(checksum, static_cast<std::uint64_t>(family.q2));
        checksum = mix(checksum, static_cast<std::uint64_t>(family.qd));
        checksum = mix(checksum, static_cast<std::uint64_t>(family.q3));
        checksum = mix(checksum, static_cast<std::uint64_t>(family.q4));
    }
    return checksum;
}

std::uint64_t checksumBitquadFamilyCounts(std::span<const BitquadFamilyCounts> bitquadFamilyCounts) {
    std::uint64_t checksum = 0;
    for (const auto& family : bitquadFamilyCounts) {
        checksum = mix(checksum, static_cast<std::uint64_t>(family.empty));
        checksum = mix(checksum, static_cast<std::uint64_t>(family.q1));
        checksum = mix(checksum, static_cast<std::uint64_t>(family.q2));
        checksum = mix(checksum, static_cast<std::uint64_t>(family.qd));
        checksum = mix(checksum, static_cast<std::uint64_t>(family.q3));
        checksum = mix(checksum, static_cast<std::uint64_t>(family.q4));
    }
    return checksum;
}

std::uint64_t checksumContourSideCounts(std::span<const ContourSideCounts> counts) {
    std::uint64_t checksum = 0;
    for (const auto& count : counts) {
        checksum = mix(checksum, static_cast<std::uint64_t>(count.contourPixels));
        checksum = mix(checksum, static_cast<std::uint64_t>(count.exposedSides));
        checksum = mix(checksum, static_cast<std::uint64_t>(count.north));
        checksum = mix(checksum, static_cast<std::uint64_t>(count.west));
        checksum = mix(checksum, static_cast<std::uint64_t>(count.east));
        checksum = mix(checksum, static_cast<std::uint64_t>(count.south));
    }
    return checksum;
}

std::uint64_t checksumFloatBuffer(std::span<const float> buffer) {
    std::uint64_t checksum = 0;
    for (float value : buffer) {
        checksum = mix(checksum, std::bit_cast<std::uint32_t>(value));
    }
    return checksum;
}

AttributeNames makeDenseAttributeNames(const std::vector<Attribute>& attributes) {
    std::unordered_map<Attribute, int> offsets;
    for (int i = 0; i < static_cast<int>(attributes.size()); ++i) {
        offsets[attributes[static_cast<std::size_t>(i)]] = i;
    }
    return AttributeNames(std::move(offsets));
}

template <AltitudeValue T> std::uint64_t computeLocalScalarChecksum(const MorphologicalTree& tree, std::span<const T> altitude) {
    const auto attributes = runtimeProducedAttributes<BitquadAttributeComputer>();
    const AttributeNames names = makeDenseAttributeNames(attributes);
    std::vector<float> buffer(static_cast<std::size_t>(tree.getNumInternalNodeSlots()) * static_cast<std::size_t>(names.NUM_ATTRIBUTES), 0.0f);
    BitquadAttributeComputer::compute(
        AltitudeAttributeComputeContext<float, T>{tree, altitude, std::span<float>(buffer), names, std::span<const Attribute>(attributes)});
    return checksumFloatBuffer(buffer);
}

std::uint64_t computeLocalScalarChecksum(const MorphologicalTree& tree) {
    const auto attributes = runtimeProducedAttributes<BitquadAttributeComputer>();
    const AttributeNames names = makeDenseAttributeNames(attributes);
    std::vector<float> buffer(static_cast<std::size_t>(tree.getNumInternalNodeSlots()) * static_cast<std::size_t>(names.NUM_ATTRIBUTES), 0.0f);
    BitquadAttributeComputer::compute(AttributeComputeContext<float>{tree, std::span<float>(buffer), names, std::span<const Attribute>(attributes)});
    return checksumFloatBuffer(buffer);
}

void printMeasurement(const std::string& label, const Measurement& measurement) {
    std::cout << "    " << std::left << std::setw(34) << label << std::right << std::setw(10) << measurement.msPerRun << " ms/run"
              << " checksum=" << measurement.checksum << '\n';
}

void printTreeStats(const MorphologicalTree& tree, double buildMs) {
    std::cout << "    build=" << buildMs << " ms"
              << " live_nodes=" << tree.getNumNodes() << " slots=" << tree.getNumInternalNodeSlots() << " proper_parts=" << tree.getNumTotalProperParts()
              << '\n';

    const int numNodes = tree.getNumInternalNodeSlots();
    printMemoryEstimate("local 16-state histograms", static_cast<std::size_t>(numNodes) * sizeof(BitquadLocalEventComputation::BitquadStateHistogram),
                        numNodes);
    printMemoryEstimate("local direct family counts", static_cast<std::size_t>(numNodes) * sizeof(BitquadFamilyCounts), numNodes);
    printMemoryEstimate("bitquad scalar buffer", static_cast<std::size_t>(numNodes) * 9 * sizeof(float), numNodes);
    printMemoryEstimate("local contour side counts", static_cast<std::size_t>(numNodes) * sizeof(ContourSideCounts), numNodes);
}

TreeOfShapesPolarityStats computeTreeOfShapesPolarityStats(const WeightedMorphologicalTree<std::uint8_t>& weightedTree) {
    const MorphologicalTree& tree = weightedTree.topology();
    TreeOfShapesPolarityStats stats;
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        if (tree.isRoot(nodeId)) {
            ++stats.rootNodes;
            continue;
        }

        const NodeId parentNodeId = tree.getNodeParent(nodeId);
        const std::uint8_t nodeAltitude = weightedTree.getAltitude(nodeId);
        const std::uint8_t parentAltitude = weightedTree.getAltitude(parentNodeId);
        if (nodeAltitude > parentAltitude) {
            ++stats.maxTreeNodes;
        } else if (nodeAltitude < parentAltitude) {
            ++stats.minTreeNodes;
        } else {
            ++stats.equalAltitudeEdges;
        }
    }
    return stats;
}

void printTreeOfShapesPolicyStats(const WeightedMorphologicalTree<std::uint8_t>& weightedTree, ToSInterpolation interpolation) {
    const MorphologicalTree& tree = weightedTree.topology();
    const TreeOfShapesPolarityStats stats = computeTreeOfShapesPolarityStats(weightedTree);
    std::cout << "    tos_interpolation=" << toString(interpolation) << " decreasing_radius=" << tree.getDecreasingGridAdjacency2D()->getRadius()
              << " increasing_radius=" << tree.getIncreasingGridAdjacency2D()->getRadius() << " root_nodes=" << stats.rootNodes
              << " min_nodes=" << stats.minTreeNodes << " max_nodes=" << stats.maxTreeNodes << " equal_altitude_edges=" << stats.equalAltitudeEdges << '\n';
}

void runComponentTreeCase(const std::string& label, const ImageUInt8Ptr& image, bool isMaxTree, double radius, int repeats) {
    auto [weightedTree, buildMs] =
        timed([&]() { return isMaxTree ? MorphologicalTreeFactory::createMaxTree(image, radius) : MorphologicalTreeFactory::createMinTree(image, radius); });
    const MorphologicalTree& tree = weightedTree.topology();

    std::cout << "\n  [" << label << "]\n";
    printTreeStats(tree, buildMs);

    printMeasurement("local histograms+families", measure(repeats, [&]() {
                         BitquadLocalEventComputation computer(tree);
                         return checksumLocalBitquads(computer);
                     }));

    printMeasurement("local delta families", measure(repeats, [&]() {
                         const auto familyCounts = BitquadLocalEventComputation::computeBitquadFamilyCounts(tree);
                         return checksumBitquadFamilyCounts(familyCounts);
                     }));

    const Measurement localScalar = measure(repeats, [&]() { return computeLocalScalarChecksum(tree); });
    printMeasurement("local delta scalar computer", localScalar);

    printMeasurement("local contour side counts", measure(repeats, [&]() {
                         const auto counts = ContourSideLocalEventComputation::computeContourSideCounts(tree);
                         return checksumContourSideCounts(counts);
                     }));
}

void runTreeOfShapesCase(const std::string& label, const ImageUInt8Ptr& image, ToSInterpolation interpolation, int repeats) {
    auto [weightedTree, buildMs] = timed([&]() { return MorphologicalTreeFactory::createTreeOfShapes(image, interpolation); });
    const MorphologicalTree& tree = weightedTree.topology();

    std::cout << "\n  [" << label << "]\n";
    printTreeStats(tree, buildMs);
    printTreeOfShapesPolicyStats(weightedTree, interpolation);

    printMeasurement("local histograms+families", measure(repeats, [&]() {
                         BitquadLocalEventComputation computer(tree);
                         return checksumLocalBitquads(computer);
                     }));

    printMeasurement("local delta families", measure(repeats, [&]() {
                         const auto familyCounts = BitquadLocalEventComputation::computeBitquadFamilyCounts(tree);
                         return checksumBitquadFamilyCounts(familyCounts);
                     }));

    printMeasurement("local delta scalar computer", measure(repeats, [&]() { return computeLocalScalarChecksum(tree, weightedTree.altitudeSpan()); }));

    printMeasurement("local contour side counts", measure(repeats, [&]() {
                         const auto counts = ContourSideLocalEventComputation::computeContourSideCounts(tree);
                         return checksumContourSideCounts(counts);
                     }));
}

void runImage(const Options& options, const std::string& label, const ImageUInt8Ptr& image) {
    std::cout << "\n" << label << " domain=" << image->getNumRows() << "x" << image->getNumCols() << " repeats=" << options.repeats << '\n';

    runComponentTreeCase("max-tree radius=1.0", image, true, 1.0, options.repeats);
    runComponentTreeCase("min-tree radius=1.0", image, false, 1.0, options.repeats);
    runComponentTreeCase("max-tree radius=1.5", image, true, 1.5, options.repeats);
    runComponentTreeCase("min-tree radius=1.5", image, false, 1.5, options.repeats);

    if (options.includeTos) {
        for (ToSInterpolation interpolation : options.tosInterpolations) {
            runTreeOfShapesCase("tree-of-shapes " + toString(interpolation), image, interpolation, options.repeats);
        }
    }
}

void runPattern(const Options& options, const std::string& pattern) {
    runImage(options, "pattern=" + pattern, makePatternImage(options.rows, options.cols, pattern));
}

void runRealImage(const Options& options, const std::filesystem::path& path) {
    ImageUInt8Ptr image;
    try {
        image = loadGrayscaleImage(path);
    } catch (const std::exception& ex) {
        std::cout << "\nimage=" << path.string() << " skipped reason=\"" << ex.what() << "\"\n";
        return;
    }
    runImage(options, "image=" + path.string(), image);
}

Options parseOptions(int argc, char** argv) {
    Options options;
    std::vector<std::string> positional;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--tos") {
            options.includeTos = true;
        } else if (arg == "--tos-interpolation") {
            if (++i >= argc) {
                throw std::invalid_argument("--tos-interpolation requires SelfDual, Min4cMax8c, or Min8cMax4c.");
            }
            options.includeTos = true;
            options.tosInterpolations.push_back(parseToSInterpolation(argv[i]));
        } else if (arg == "--skip-synthetic") {
            options.skipSynthetic = true;
        } else if (arg == "--image") {
            if (++i >= argc) {
                throw std::invalid_argument("--image requires a file path.");
            }
            const std::filesystem::path path = argv[i];
            if (!std::filesystem::is_regular_file(path)) {
                throw std::invalid_argument("Image file does not exist: " + path.string());
            }
            options.imagePaths.push_back(path);
        } else if (arg == "--image-dir") {
            if (++i >= argc) {
                throw std::invalid_argument("--image-dir requires a directory path.");
            }
            auto images = collectImagePaths(argv[i]);
            options.imagePaths.insert(options.imagePaths.end(), images.begin(), images.end());
        } else if (arg.rfind("--", 0) == 0) {
            throw std::invalid_argument("Unknown option: " + arg);
        } else {
            positional.push_back(arg);
        }
    }

    if (positional.size() > 3) {
        throw std::invalid_argument("Expected at most three positional arguments: rows cols repeats.");
    }
    if (!positional.empty()) {
        options.rows = std::stoi(positional[0]);
        options.cols = options.rows;
    }
    if (positional.size() > 1) {
        options.cols = std::stoi(positional[1]);
    }
    if (positional.size() > 2) {
        options.repeats = std::stoi(positional[2]);
    }
    if (options.rows <= 0 || options.cols <= 0 || options.repeats <= 0) {
        throw std::invalid_argument("rows, cols, and repeats must be positive.");
    }
    if (options.includeTos && options.tosInterpolations.empty()) {
        options.tosInterpolations.push_back(ToSInterpolation::SelfDual);
    }
    std::sort(options.tosInterpolations.begin(), options.tosInterpolations.end());
    options.tosInterpolations.erase(std::unique(options.tosInterpolations.begin(), options.tosInterpolations.end()), options.tosInterpolations.end());
    std::sort(options.imagePaths.begin(), options.imagePaths.end());
    options.imagePaths.erase(std::unique(options.imagePaths.begin(), options.imagePaths.end()), options.imagePaths.end());
    if (options.skipSynthetic && options.imagePaths.empty()) {
        throw std::invalid_argument("--skip-synthetic requires at least one --image or --image-dir input.");
    }
    return options;
}

} // namespace

int main(int argc, char** argv) {
    const Options options = parseOptions(argc, argv);

    std::cout << std::fixed << std::setprecision(3) << "local event bitquad benchmark\n"
              << "usage: mmcfilters_local_events_bitquads_benchmark [rows] [cols] [repeats] [--tos] "
                 "[--tos-interpolation name] [--image path] [--image-dir dir] [--skip-synthetic]\n"
              << "memory numbers are type-size estimates, not process RSS\n";

    if (!options.skipSynthetic) {
        for (const std::string& pattern : {"flat", "ramp", "noise", "hollow", "diagonal"}) {
            runPattern(options, pattern);
        }
    }

    for (const auto& imagePath : options.imagePaths) {
        runRealImage(options, imagePath);
    }

    return 0;
}
