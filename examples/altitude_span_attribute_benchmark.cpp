/**
 * Benchmark canonical weighted-tree attribute computation against
 * `WeightedTreeView<T>` altitude-span computation for equivalent altitude
 * buffers.
 *
 * Build with `-DMMCFILTERS_BUILD_EXAMPLES=ON` and run:
 * `./build/examples/mmcfilters_altitude_span_attribute_benchmark
 * --sizes 128,256,512 --repeats 3 --radius 1.5 --suite core`.
 *
 * The reported checksum verifies that each typed altitude-span path returns the
 * same node-indexed attribute values as the canonical uint8 owner path.
 */
#include "mmcfilters/attributes/Attributes.hpp"
#include "mmcfilters/trees/MorphologicalTreeFactory.hpp"
#include "mmcfilters/trees/WeightedMorphologicalTree.hpp"
#include "mmcfilters/trees/WeightedTreeView.hpp"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using namespace mmcfilters;

namespace {

struct Options {
    std::vector<int> sizes{128, 256, 512};
    int repeats = 3;
    double radius = 1.5;
    std::string suite = "core";
};

struct RequestCase {
    std::string name;
    std::vector<AttributeOrGroup> request;
    std::vector<Attribute> checksumAttributes;
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
            const int radial = (row - rows / 2) * (row - rows / 2) + (col - cols / 2) * (col - cols / 2);
            const int waves = (row * 17) ^ (col * 31) ^ ((row + col) * 7);
            (*image)[idx] = static_cast<uint8_t>((radial / 113 + waves) & 0xff);
        }
    }
    return image;
}

std::string formatBytes(std::size_t bytes) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(2) << (static_cast<double>(bytes) / 1024.0) << " KiB";
    return out.str();
}

std::uint64_t checksumComputed(const MorphologicalTree& tree, const ComputedAttributeData<float>& data, const std::vector<Attribute>& attributes) {
    std::uint64_t checksum = 1469598103934665603ull;
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        for (Attribute attribute : attributes) {
            const float value = data.second[data.first.linearIndex(nodeId, attribute)];
            std::uint64_t encodedValue = 0;
            if (std::isnan(value)) {
                encodedValue = 0x7fc00000ull;
            } else if (std::isinf(value)) {
                encodedValue = std::signbit(value) ? 0xff800000ull : 0x7f800000ull;
            } else {
                const auto scaled = static_cast<std::int64_t>(std::llround(static_cast<double>(value) * 1000.0));
                encodedValue = static_cast<std::uint64_t>(scaled);
            }
            checksum ^= encodedValue + 0x9e3779b97f4a7c15ull;
            checksum *= 1099511628211ull;
            checksum ^= static_cast<std::uint64_t>(nodeId + 1);
            checksum *= 1099511628211ull;
        }
    }
    return checksum;
}

template <class Fn> Measurement measure(int repeats, Fn&& fn) {
    if (repeats <= 0) {
        throw std::invalid_argument("repeats must be positive.");
    }

    std::uint64_t checksum = 0;
    const auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < repeats; ++i) {
        checksum += fn();
    }
    const auto end = std::chrono::steady_clock::now();
    const auto micros = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    return {static_cast<double>(micros) / static_cast<double>(repeats) / 1000.0, checksum};
}

template <class T> std::vector<T> makeEquivalentAltitude(const WeightedMorphologicalTree<std::uint8_t>& weighted) {
    std::vector<T> altitude(static_cast<std::size_t>(weighted.topology().getNumInternalNodeSlots()), T{});
    for (NodeId nodeId : weighted.topology().getAliveNodeIds()) {
        altitude[static_cast<std::size_t>(nodeId)] = static_cast<T>(weighted.getAltitude(nodeId));
    }
    return altitude;
}

template <class T> const char* typeName() {
    if constexpr (std::is_same_v<T, std::uint8_t>) {
        return "uint8_t";
    } else if constexpr (std::is_same_v<T, std::int32_t>) {
        return "int32_t";
    } else if constexpr (std::is_same_v<T, float>) {
        return "float";
    } else if constexpr (std::is_same_v<T, double>) {
        return "double";
    } else {
        return "unknown";
    }
}

template <class T>
void runAltitudeSpanType(const WeightedMorphologicalTree<std::uint8_t>& weighted, const RequestCase& requestCase, std::uint64_t baselineChecksum, int repeats) {
    const MorphologicalTree& tree = weighted.topology();
    const std::vector<T> altitude = makeEquivalentAltitude<T>(weighted);
    const WeightedTreeView<T> view(tree, std::span<const T>(altitude));
    const std::size_t altitudeBytes = altitude.size() * sizeof(T);
    const Measurement measurement = measure(repeats, [&]() {
        const auto data = AttributeComputation::computeAttributesFromAltitudeSpan(view, requestCase.request);
        return checksumComputed(tree, data, requestCase.checksumAttributes);
    });

    const bool checksumMatches = (measurement.checksum == baselineChecksum);
    std::cout << "    altitude-span " << std::setw(8) << typeName<T>() << "  altitude=" << std::setw(10) << formatBytes(altitudeBytes)
              << "  ms/run=" << std::setw(9) << std::fixed << std::setprecision(3) << measurement.msPerRun << "  checksum=" << measurement.checksum
              << "  match=" << (checksumMatches ? "yes" : "no") << '\n';
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
        } else if (arg == "--suite") {
            if (++i >= argc) {
                throw std::invalid_argument("--suite requires one of: core, topology, all.");
            }
            options.suite = argv[i];
            if (options.suite != "core" && options.suite != "topology" && options.suite != "all") {
                throw std::invalid_argument("--suite requires one of: core, topology, all.");
            }
        } else if (arg == "--help" || arg == "-h") {
            std::cout << "usage: mmcfilters_altitude_span_attribute_benchmark "
                      << "[--sizes 128,256,512] [--repeats 3] [--radius 1.5] "
                      << "[--suite core|topology|all]\n";
            std::exit(0);
        } else {
            throw std::invalid_argument("unknown argument: " + arg);
        }
    }
    return options;
}

std::vector<Attribute> expandChecksumAttributes(const std::vector<AttributeOrGroup>& request) {
    std::set<Attribute> expanded;
    for (const AttributeOrGroup& item : request) {
        if (std::holds_alternative<Attribute>(item)) {
            expanded.insert(std::get<Attribute>(item));
            continue;
        }
        const auto it = ATTRIBUTE_GROUPS.find(std::get<AttributeGroup>(item));
        if (it == ATTRIBUTE_GROUPS.end()) {
            throw std::invalid_argument("benchmark request uses an unknown AttributeGroup.");
        }
        expanded.insert(it->second.begin(), it->second.end());
    }
    return {expanded.begin(), expanded.end()};
}

RequestCase makeRequestCase(std::string name, std::vector<AttributeOrGroup> request) {
    std::vector<Attribute> checksumAttributes = expandChecksumAttributes(request);
    return {std::move(name), std::move(request), std::move(checksumAttributes)};
}

std::vector<RequestCase> coreRequestCases() {
    return {
        {
            "level_gray",
            {LEVEL, GRAY_HEIGHT},
            expandChecksumAttributes({LEVEL, GRAY_HEIGHT}),
        },
        {
            "volume_relative",
            {VOLUME, RELATIVE_VOLUME},
            expandChecksumAttributes({VOLUME, RELATIVE_VOLUME}),
        },
        {
            "level_volume_maxdist",
            {LEVEL, VOLUME, MAX_DIST},
            expandChecksumAttributes({LEVEL, VOLUME, MAX_DIST}),
        },
    };
}

std::vector<RequestCase> topologyOnlyRequestCases() {
    return {
        makeRequestCase("gray_level", {AttributeGroup::GRAY_LEVEL}),
        makeRequestCase("shape", {AttributeGroup::SHAPE}),
        makeRequestCase("moments", {AttributeGroup::MOMENTS}),
        makeRequestCase("boundary", {AttributeGroup::BOUNDARY}),
        makeRequestCase("tree_topology", {AttributeGroup::TREE_TOPOLOGY}),
        makeRequestCase("all", {AttributeGroup::ALL}),
        makeRequestCase("mixed_level_shape_boundary", {LEVEL, AttributeGroup::MOMENTS, AttributeGroup::BOUNDARY}),
    };
}

std::vector<RequestCase> requestCases(const std::string& suite) {
    if (suite == "core") {
        return coreRequestCases();
    }
    if (suite == "topology") {
        return topologyOnlyRequestCases();
    }
    std::vector<RequestCase> cases = coreRequestCases();
    auto extra = topologyOnlyRequestCases();
    cases.insert(cases.end(), std::make_move_iterator(extra.begin()), std::make_move_iterator(extra.end()));
    return cases;
}

void runCase(int rows, int cols, bool isMaxTree, const Options& options) {
    auto image = makeBenchmarkImage(rows, cols);
    const auto buildStart = std::chrono::steady_clock::now();
    auto weighted = isMaxTree ? MorphologicalTreeFactory::createMaxTree(image, options.radius) : MorphologicalTreeFactory::createMinTree(image, options.radius);
    const auto buildEnd = std::chrono::steady_clock::now();
    const auto buildMicros = std::chrono::duration_cast<std::chrono::microseconds>(buildEnd - buildStart).count();

    const MorphologicalTree& tree = weighted.topology();
    const std::size_t slots = static_cast<std::size_t>(tree.getNumInternalNodeSlots());
    std::cout << "\n"
              << rows << "x" << cols << " " << (isMaxTree ? "max" : "min") << " nodes=" << tree.getNumNodes() << " slots=" << slots
              << " build_ms=" << std::fixed << std::setprecision(3) << (static_cast<double>(buildMicros) / 1000.0) << '\n';

    for (const RequestCase& requestCase : requestCases(options.suite)) {
        const std::size_t resultBytes = slots * requestCase.checksumAttributes.size() * sizeof(float);
        std::cout << "  request=" << requestCase.name << " attrs=" << requestCase.checksumAttributes.size() << " result_buffer=" << formatBytes(resultBytes)
                  << '\n';

        const Measurement baseline = measure(options.repeats, [&]() {
            const auto data = AttributeComputation::computeAttributes(weighted, requestCase.request);
            return checksumComputed(tree, data, requestCase.checksumAttributes);
        });

        std::cout << "    canonical"
                  << "          "
                  << "  altitude=" << std::setw(10) << formatBytes(slots * sizeof(std::uint8_t)) << "  ms/run=" << std::setw(9) << std::fixed
                  << std::setprecision(3) << baseline.msPerRun << "  checksum=" << baseline.checksum << "  match=yes" << '\n';

        runAltitudeSpanType<std::uint8_t>(weighted, requestCase, baseline.checksum, options.repeats);
        runAltitudeSpanType<std::int32_t>(weighted, requestCase, baseline.checksum, options.repeats);
        runAltitudeSpanType<float>(weighted, requestCase, baseline.checksum, options.repeats);
        runAltitudeSpanType<double>(weighted, requestCase, baseline.checksum, options.repeats);
    }
}

} // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);
        std::cout << "Altitude-span attribute benchmark\n"
                  << "repeats=" << options.repeats << " radius=" << options.radius << " suite=" << options.suite << '\n';

        for (int size : options.sizes) {
            runCase(size, size, true, options);
            runCase(size, size, false, options);
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
}
