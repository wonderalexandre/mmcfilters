#include "mmcfilters/attributes/AttributeComputation.hpp"
#include "mmcfilters/filters/AttributeFilters.hpp"
#include "mmcfilters/filters/DepthStableRegionComputer.hpp"
#include "mmcfilters/filters/ExtinctionValues.hpp"
#include "mmcfilters/filters/MSERComputer.hpp"
#include "mmcfilters/filters/UltimateAttributeOpening.hpp"
#include "mmcfilters/trees/MorphologicalTreeFactory.hpp"
#include "mmcfilters/utils/Contract.hpp"
#include "mmcfilters/utils/Image.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
constexpr std::uint64_t fnvPrime = 1099511628211ULL;

struct ScenarioResult {
    std::string name;
    double medianMilliseconds = 0.0;
    std::uint64_t checksum = 0;
};

struct UaoResult {
    mmcfilters::ImageUInt8Ptr maxContrast;
    mmcfilters::ImageInt32Ptr associated;
};

mmcfilters::ImageUInt8Ptr makeInput(int rows, int columns) {
    auto image = mmcfilters::ImageUInt8::create(rows, columns);
    for (int row = 0; row < rows; ++row) {
        for (int column = 0; column < columns; ++column) {
            const int centeredRow = row - rows / 2;
            const int centeredColumn = column - columns / 2;
            const int radial = centeredRow * centeredRow + centeredColumn * centeredColumn;
            (*image)[row * columns + column] = static_cast<std::uint8_t>((radial + 17 * row + 31 * column) & 0xff);
        }
    }
    return image;
}

template <typename Value>
    requires std::is_trivially_copyable_v<Value>
void hashValue(std::uint64_t& hash, const Value& value) noexcept {
    const auto* bytes = reinterpret_cast<const unsigned char*>(&value);
    for (std::size_t index = 0; index < sizeof(Value); ++index) {
        hash ^= static_cast<std::uint64_t>(bytes[index]);
        hash *= fnvPrime;
    }
}

template <typename Value> void hashVector(std::uint64_t& hash, const std::vector<Value>& values) noexcept {
    hashValue(hash, values.size());
    for (const Value& value : values) {
        hashValue(hash, value);
    }
}

template <typename Value> void hashImage(std::uint64_t& hash, const mmcfilters::ImagePtr<Value>& image) noexcept {
    const int size = image->getSize();
    hashValue(hash, size);
    for (int index = 0; index < size; ++index) {
        hashValue(hash, (*image)[index]);
    }
}

void hashAttributeLayout(std::uint64_t& hash, const mmcfilters::AttributeNames& names) {
    std::vector<std::pair<mmcfilters::Attribute, int>> entries(names.indexMap.begin(), names.indexMap.end());
    std::sort(entries.begin(), entries.end(), [](const auto& lhs, const auto& rhs) {
        return static_cast<int>(lhs.first) < static_cast<int>(rhs.first);
    });
    hashValue(hash, names.NUM_ATTRIBUTES);
    for (const auto& [attribute, offset] : entries) {
        hashValue(hash, attribute);
        hashValue(hash, offset);
    }
}

void hashNodeAttributeSampleLayout(std::uint64_t& hash, const mmcfilters::NodeAttributeSampleLayout& names) {
    std::vector<std::pair<mmcfilters::NodeAttributeSampleKey, int>> entries(names.indexMap.begin(), names.indexMap.end());
    std::sort(entries.begin(), entries.end(), [](const auto& lhs, const auto& rhs) {
        if (lhs.first.sampleOffset != rhs.first.sampleOffset) {
            return lhs.first.sampleOffset < rhs.first.sampleOffset;
        }
        return static_cast<int>(lhs.first.attribute) < static_cast<int>(rhs.first.attribute);
    });
    hashValue(hash, names.NUM_ATTRIBUTES);
    for (const auto& [key, offset] : entries) {
        hashValue(hash, key.attribute);
        hashValue(hash, key.sampleOffset);
        hashValue(hash, offset);
    }
}

template <typename Operation, typename Checksum>
ScenarioResult benchmarkScenario(std::string name, int repetitions, Operation&& operation, Checksum&& computeChecksum) {
    auto warmup = operation();
    const std::uint64_t expectedChecksum = computeChecksum(warmup);

    std::vector<double> samples;
    samples.reserve(static_cast<std::size_t>(repetitions));
    for (int repetition = 0; repetition < repetitions; ++repetition) {
        const auto start = Clock::now();
        auto result = operation();
        const auto finish = Clock::now();
        samples.push_back(std::chrono::duration<double, std::milli>(finish - start).count());

        const std::uint64_t checksum = computeChecksum(result);
        if (checksum != expectedChecksum) {
            throw std::runtime_error(name + " produced a non-deterministic result across repetitions.");
        }
    }

    std::sort(samples.begin(), samples.end());
    return ScenarioResult{std::move(name), samples[samples.size() / 2], expectedChecksum};
}

template <typename Real> std::uint64_t attributeChecksum(const mmcfilters::ComputedAttributeData<Real>& result) {
    std::uint64_t hash = fnvOffsetBasis;
    hashAttributeLayout(hash, result.first);
    hashValue(hash, result.nodeIdSpace);
    hashVector(hash, result.second);
    return hash;
}

template <typename Real> std::uint64_t sampledNodeAttributeChecksum(const mmcfilters::SampledNodeAttributeData<Real>& result) {
    std::uint64_t hash = fnvOffsetBasis;
    hashNodeAttributeSampleLayout(hash, result.first);
    hashValue(hash, result.nodeIdSpace);
    hashVector(hash, result.second);
    return hash;
}

template <typename Value> std::uint64_t imageChecksum(const mmcfilters::ImagePtr<Value>& result) noexcept {
    std::uint64_t hash = fnvOffsetBasis;
    hashImage(hash, result);
    return hash;
}

std::uint64_t maskChecksum(const std::vector<std::uint8_t>& result) noexcept {
    std::uint64_t hash = fnvOffsetBasis;
    hashVector(hash, result);
    return hash;
}

std::uint64_t uaoChecksum(const UaoResult& result) noexcept {
    std::uint64_t hash = fnvOffsetBasis;
    hashImage(hash, result.maxContrast);
    hashImage(hash, result.associated);
    return hash;
}

template <mmcfilters::AltitudeValue T>
std::uint64_t higraChecksum(const std::pair<std::vector<mmcfilters::NodeId>, std::vector<T>>& result) noexcept {
    std::uint64_t hash = fnvOffsetBasis;
    hashVector(hash, result.first);
    hashVector(hash, result.second);
    return hash;
}

} // namespace

int main(int argc, char** argv) {
    const int rows = argc > 1 ? std::atoi(argv[1]) : 256;
    const int columns = argc > 2 ? std::atoi(argv[2]) : 256;
    const int repetitions = argc > 3 ? std::atoi(argv[3]) : 9;
    if (rows <= 0 || columns <= 0 || repetitions <= 0) {
        std::cerr << "usage: validation_sensitive_algorithms_benchmark [rows>0] [columns>0] [repetitions>0]\n";
        return EXIT_FAILURE;
    }

    try {
        // The image, topology, altitude and shared increasing attribute are the
        // established input state. Their construction is intentionally outside
        // every timed region so each scenario measures only its named operation.
        const auto image = makeInput(rows, columns);
        const auto valuedTree = mmcfilters::MorphologicalTreeFactory::createMaxTree(image, 1.5);
        const auto area = mmcfilters::AttributeComputation::computeSingleAttribute<double>(valuedTree, mmcfilters::Area);
        const double attributeThreshold = static_cast<double>(rows) * static_cast<double>(columns) / 16.0;

        std::vector<ScenarioResult> results;
        results.reserve(8);

        results.push_back(benchmarkScenario(
            "viterbi_filter", repetitions,
            [&]() {
                mmcfilters::AttributeFilters<std::uint8_t> filters(valuedTree);
                return filters.filteringByViterbiRule(area.second.data(), attributeThreshold);
            },
            [](const auto& result) { return imageChecksum(result); }));

        results.push_back(benchmarkScenario(
            "ultimate_attribute_opening", repetitions,
            [&]() {
                mmcfilters::UltimateAttributeOpening<std::uint8_t, double> uao(valuedTree, area.second.data());
                uao.execute(attributeThreshold);
                return UaoResult{uao.getMaxContrastImage(), uao.getAssociatedImage()};
            },
            [](const UaoResult& result) { return uaoChecksum(result); }));

        results.push_back(benchmarkScenario(
            "extinction_filter", repetitions,
            [&]() {
                mmcfilters::ExtinctionValues<std::uint8_t, double> extinction(valuedTree, area.second.data());
                return extinction.filtering(mmcfilters::ExtinctionSelectionPolicy<double>::byTopK(16));
            },
            [](const auto& result) { return imageChecksum(result); }));

        results.push_back(benchmarkScenario(
            "mser", repetitions,
            [&]() {
                mmcfilters::MSERComputer<std::uint8_t, double> mser(valuedTree, area.second.data());
                return mser.computeMSER(mmcfilters::AltitudeDifference<std::uint8_t>{4});
            },
            [](const std::vector<std::uint8_t>& result) { return maskChecksum(result); }));

        results.push_back(benchmarkScenario(
            "depth_stability", repetitions,
            [&]() {
                mmcfilters::DepthStableRegionComputer<double> depth(valuedTree.topology(), area.second.data());
                return depth.computeByDepth(4);
            },
            [](const std::vector<std::uint8_t>& result) { return maskChecksum(result); }));

        results.push_back(benchmarkScenario(
            "sampled_node_attribute", repetitions,
            [&]() {
                return mmcfilters::AttributeComputation::computeSampledNodeAttribute<double>(
                    valuedTree, mmcfilters::Area, mmcfilters::AltitudeDifference<std::uint8_t>{1}, 2);
            },
            [](const auto& result) { return sampledNodeAttributeChecksum(result); }));

        results.push_back(benchmarkScenario(
            "bitquad_attributes", repetitions,
            [&]() {
                return mmcfilters::AttributeComputation::computeAttributes<double>(
                    valuedTree, {mmcfilters::BitquadArea, mmcfilters::BitquadPerimeter, mmcfilters::BitquadCircularity});
            },
            [](const auto& result) { return attributeChecksum(result); }));

        results.push_back(benchmarkScenario(
            "higra_export", repetitions, [&]() { return valuedTree.exportHigraHierarchy(); },
            [](const auto& result) { return higraChecksum(result); }));

        std::cout << std::fixed << std::setprecision(3)
                  << "contract_mode=" << (mmcfilters::contract::validationsEnabled ? "CHECKED" : "UNCHECKED") << '\n'
                  << "rows=" << rows << '\n'
                  << "columns=" << columns << '\n'
                  << "repetitions=" << repetitions << '\n'
                  << "internal_node_slots=" << valuedTree.topology().numInternalNodeSlots() << '\n';
        for (const ScenarioResult& result : results) {
            std::cout << result.name << "_median_ms=" << result.medianMilliseconds << '\n'
                      << result.name << "_checksum=" << result.checksum << '\n';
        }
    } catch (const std::exception& error) {
        std::cerr << "benchmark failed: " << error.what() << '\n';
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
