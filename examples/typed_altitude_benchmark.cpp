/**
 * Benchmark typed component-tree construction and attribute computation.
 *
 * Build with `-DMMCFILTERS_BUILD_EXAMPLES=ON` and run:
 * `./build/examples/mmcfilters_typed_altitude_benchmark --sizes 256,512
 * --repeats 3 --radius 1.5 --suite both`.
 *
 * Output is CSV comparing `uint8_t`, `int32_t`, and `float` altitude owners.
 * This demonstrates the typed C++ boundary; Python remains uint8-only.
 */
#include "mmcfilters/attributes/Attributes.hpp"
#include "mmcfilters/trees/MorphologicalTreeFactory.hpp"
#include "mmcfilters/utils/Image.hpp"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

using namespace mmcfilters;

namespace {

struct Options {
    std::vector<int> sizes{256, 512, 1024};
    int repeats = 3;
    double radius = 1.5;
    std::string suite = "both";
};

struct Measurement {
    double msPerRun = 0.0;
    std::uint64_t checksum = 0;
};

struct AttributeRequest {
    std::string name;
    std::vector<AttributeOrGroup> request;
    std::vector<Attribute> checksumAttributes;
};

template<class T>
const char* typeName() {
    if constexpr (std::is_same_v<T, std::uint8_t>) {
        return "uint8";
    } else if constexpr (std::is_same_v<T, std::int32_t>) {
        return "int32";
    } else if constexpr (std::is_same_v<T, float>) {
        return "float32";
    } else {
        return "unknown";
    }
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
                throw std::invalid_argument("--suite requires one of: light, all, both.");
            }
            options.suite = argv[i];
            if (options.suite != "light" && options.suite != "all" && options.suite != "both") {
                throw std::invalid_argument("--suite requires one of: light, all, both.");
            }
        } else if (arg == "--help" || arg == "-h") {
            std::cout
                << "usage: mmcfilters_typed_altitude_benchmark "
                << "[--sizes 256,512,1024] [--repeats 3] [--radius 1.5] "
                << "[--suite light|all|both]\n";
            std::exit(0);
        } else {
            throw std::invalid_argument("unknown argument: " + arg);
        }
    }
    return options;
}

template<class T>
ImagePtr<T> makeBenchmarkImage(int rows, int cols) {
    auto image = Image<T>::create(rows, cols);
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            const int idx = row * cols + col;
            const int radial = (row - rows / 2) * (row - rows / 2) +
                               (col - cols / 2) * (col - cols / 2);
            const int waves = (row * 17) ^ (col * 31) ^ ((row + col) * 7);
            const int value = (radial / 113 + waves) & 0xff;
            (*image)[idx] = static_cast<T>(value);
        }
    }
    return image;
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

std::vector<AttributeRequest> requestsForSuite(const std::string& suite) {
    std::vector<AttributeRequest> requests;
    if (suite == "light" || suite == "both") {
        std::vector<AttributeOrGroup> light{LEVEL, AREA, VOLUME, GRAY_HEIGHT};
        requests.push_back({"light", light, expandChecksumAttributes(light)});
    }
    if (suite == "all" || suite == "both") {
        std::vector<AttributeOrGroup> all{AttributeGroup::ALL};
        requests.push_back({"all", all, expandChecksumAttributes(all)});
    }
    return requests;
}

template<class Fn>
Measurement measure(int repeats, Fn&& fn) {
    std::uint64_t checksum = 0;
    const auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < repeats; ++i) {
        checksum += fn();
    }
    const auto end = std::chrono::steady_clock::now();
    const auto micros = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    return {static_cast<double>(micros) / static_cast<double>(repeats) / 1000.0, checksum};
}

std::uint64_t checksumTree(const MorphologicalTree& tree) {
    std::uint64_t checksum = 1469598103934665603ull;
    checksum ^= static_cast<std::uint64_t>(tree.getNumNodes());
    checksum *= 1099511628211ull;
    checksum ^= static_cast<std::uint64_t>(tree.getNumInternalNodeSlots());
    checksum *= 1099511628211ull;
    checksum ^= static_cast<std::uint64_t>(tree.getNumTotalProperParts());
    checksum *= 1099511628211ull;
    return checksum;
}

std::uint64_t checksumComputed(
    const MorphologicalTree& tree,
    const ComputedAttributeData<float>& data,
    const std::vector<Attribute>& attributes)
{
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

template<class T>
std::uint64_t checksumAltitude(const WeightedMorphologicalTree<T>& weighted) {
    std::uint64_t checksum = checksumTree(weighted.topology());
    for (NodeId nodeId : weighted.topology().getAliveNodeIds()) {
        const auto scaled = static_cast<std::int64_t>(
            std::llround(static_cast<double>(weighted.getAltitude(nodeId)) * 1000.0));
        checksum ^= static_cast<std::uint64_t>(scaled) + 0x9e3779b97f4a7c15ull;
        checksum *= 1099511628211ull;
    }
    return checksum;
}

template<class T>
auto createTree(const ImagePtr<T>& image, bool maxTree, double radius) {
    return maxTree
        ? MorphologicalTreeFactory::createMaxTree(image, radius)
        : MorphologicalTreeFactory::createMinTree(image, radius);
}

template<class T>
void runTypedCase(
    int size,
    bool maxTree,
    const Options& options,
    const std::vector<AttributeRequest>& requests)
{
    auto image = makeBenchmarkImage<T>(size, size);
    const Measurement build = measure(options.repeats, [&]() {
        auto weighted = createTree<T>(image, maxTree, options.radius);
        return checksumAltitude(weighted);
    });

    auto weighted = createTree<T>(image, maxTree, options.radius);
    const MorphologicalTree& tree = weighted.topology();
    const std::size_t slots = static_cast<std::size_t>(tree.getNumInternalNodeSlots());
    const std::size_t altitudeBytes = slots * sizeof(T);

    for (const AttributeRequest& request : requests) {
        const std::size_t resultBytes =
            slots * request.checksumAttributes.size() * sizeof(float);
        const Measurement attr = measure(options.repeats, [&]() {
            const auto data = AttributeComputation::computeAttributes(
                weighted,
                request.request);
            return checksumComputed(tree, data, request.checksumAttributes);
        });

        std::cout
            << size << ','
            << (maxTree ? "max" : "min") << ','
            << typeName<T>() << ','
            << request.name << ','
            << tree.getNumNodes() << ','
            << tree.getNumInternalNodeSlots() << ','
            << request.checksumAttributes.size() << ','
            << altitudeBytes << ','
            << resultBytes << ','
            << std::fixed << std::setprecision(3) << build.msPerRun << ','
            << std::fixed << std::setprecision(3) << attr.msPerRun << ','
            << build.checksum << ','
            << attr.checksum
            << '\n';
    }
}

void runSize(int size, const Options& options, const std::vector<AttributeRequest>& requests) {
    for (bool maxTree : {true, false}) {
        runTypedCase<std::uint8_t>(size, maxTree, options, requests);
        runTypedCase<std::int32_t>(size, maxTree, options, requests);
        runTypedCase<float>(size, maxTree, options, requests);
    }
}

} // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);
        const auto requests = requestsForSuite(options.suite);
        std::cout
            << "size,tree,dtype,request,nodes,slots,attrs,altitude_bytes,result_bytes,"
            << "build_ms,attributes_ms,build_checksum,attributes_checksum\n";
        for (int size : options.sizes) {
            runSize(size, options, requests);
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
}
