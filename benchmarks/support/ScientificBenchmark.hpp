#pragma once

#include "mmcfilters/attributes/AttributeRegistry.hpp"
#include "mmcfilters/attributes/AttributeResultTypes.hpp"
#include "mmcfilters/trees/ValuedMorphologicalTree.hpp"
#include "mmcfilters/utils/Contract.hpp"
#include "mmcfilters/utils/Image.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace mmcfilters::benchmarks {

using BenchmarkClock = std::chrono::steady_clock;

enum class TimingScope { EndToEnd, EstablishedInput };

/** @brief Structural size of the input processed by a benchmark scenario. */
struct WorkloadMetrics {
    std::int64_t pixels = 0;             ///< Number of pixels in the input domain.
    std::int64_t properParts = 0;        ///< Number of proper-part samples represented by the workload.
    std::int64_t primaryNodeSlots = 0;   ///< Allocated node slots in the primary tree.
    std::int64_t primaryLiveNodes = 0;   ///< Live nodes in the primary tree.
    std::int64_t secondaryNodeSlots = 0; ///< Allocated node slots in the secondary tree, when present.
    std::int64_t secondaryLiveNodes = 0; ///< Live nodes in the secondary tree, when present.
    std::int64_t edges = 0;              ///< Undirected adjacency edges represented by the workload.
};

/** @brief Structural edits and thresholds observed during a benchmark scenario. */
struct OutcomeMetrics {
    std::int64_t steps = 0;                        ///< Number of editing or filtering steps.
    std::int64_t primaryNodesRemoved = 0;          ///< Nodes removed from the primary tree.
    std::int64_t secondaryNodesRemoved = 0;        ///< Nodes removed from the secondary tree.
    std::int64_t completeValidationCommits = 0;    ///< Commits that required complete validation.
    std::int64_t incrementalValidationCommits = 0; ///< Commits that used incremental validation.
    double lightThreshold = 0.0;                   ///< Threshold selected for the light-edit step.
    double mediumThreshold = 0.0;                  ///< Threshold selected for the medium-edit step.
    double heavyThreshold = 0.0;                   ///< Threshold selected for the heavy-edit step.
};

/** @brief Timing, checksum, and workload metadata for one benchmark scenario. */
struct ScenarioResult {
    std::string suite;                                            ///< Benchmark suite containing the scenario.
    std::string name;                                             ///< Scenario name within the suite and scope.
    TimingScope scope = TimingScope::EstablishedInput;            ///< Boundary of the timed operation.
    double medianMilliseconds = 0.0;                              ///< Median elapsed time in milliseconds.
    double madMilliseconds = 0.0;                                 ///< Median absolute deviation in milliseconds.
    double minimumMilliseconds = 0.0;                             ///< Minimum elapsed time in milliseconds.
    double maximumMilliseconds = 0.0;                             ///< Maximum elapsed time in milliseconds.
    std::uint64_t checksum = 0;                                   ///< Deterministic checksum of the scientific result.
    WorkloadMetrics metrics;                                      ///< Structural size of the measured workload.
    OutcomeMetrics outcome;                                       ///< Editing outcomes associated with the result.
    std::vector<double> samplesMilliseconds;                      ///< Individual elapsed-time samples.

    /** @brief Returns the stable suite, scope, and scenario identifier. @return Fully qualified scenario name. */
    [[nodiscard]] std::string fullName() const {
        return suite + "." + (scope == TimingScope::EndToEnd ? "end_to_end." : "established_input.") + name;
    }
};

/** @brief Incremental 64-bit FNV-1a checksum used to prevent dead-code elimination and verify result equality. */
class Fnv1a64 {
  private:
    static constexpr std::uint64_t offsetBasis = 14695981039346656037ULL; ///< FNV-1a 64-bit offset basis.
    static constexpr std::uint64_t prime = 1099511628211ULL;              ///< FNV-1a 64-bit prime.
    std::uint64_t value_ = offsetBasis;                                  ///< Current incremental hash state.

  public:
    /**
     * @brief Appends an uninterpreted byte sequence to the checksum.
     * @param bytes First byte in the sequence.
     * @param size Number of bytes to append.
     */
    void appendBytes(const unsigned char* bytes, std::size_t size) noexcept {
        for (std::size_t index = 0; index < size; ++index) {
            value_ ^= static_cast<std::uint64_t>(bytes[index]);
            value_ *= prime;
        }
    }

    /**
     * @brief Appends the object representation of one trivially copyable value.
     * @param value Value appended to the checksum.
     */
    template <typename Value>
        requires std::is_trivially_copyable_v<Value>
    void append(const Value& value) noexcept {
        appendBytes(reinterpret_cast<const unsigned char*>(&value), sizeof(Value));
    }

    /** @brief Appends a length-prefixed string. @param value String appended to the checksum. */
    void append(std::string_view value) noexcept {
        append(value.size());
        appendBytes(reinterpret_cast<const unsigned char*>(value.data()), value.size());
    }

    /** @brief Appends a length-prefixed vector. @param values Values appended in storage order. */
    template <typename Value> void appendVector(const std::vector<Value>& values) noexcept {
        append(values.size());
        for (const Value& value : values) {
            append(value);
        }
    }

    /** @brief Returns the current checksum value. @return Current 64-bit FNV-1a state. */
    [[nodiscard]] std::uint64_t value() const noexcept { return value_; }
};

template <typename Value> void appendImage(Fnv1a64& hash, const ImagePtr<Value>& image) {
    if (!image) {
        throw std::invalid_argument("Scientific benchmark checksum requires a non-null image.");
    }
    hash.append(image->getNumRows());
    hash.append(image->getNumColumns());
    hash.append(image->getSize());
    for (int index = 0; index < image->getSize(); ++index) {
        hash.append((*image)[index]);
    }
}

template <typename Value> [[nodiscard]] std::uint64_t imageChecksum(const ImagePtr<Value>& image) {
    Fnv1a64 hash;
    appendImage(hash, image);
    return hash.value();
}

[[nodiscard]] inline std::uint64_t benchmarkInputChecksum(const ImageUInt8Ptr& image) {
    if (!image) {
        throw std::invalid_argument("Scientific benchmark input checksum requires a non-null image.");
    }
    const std::string header = std::to_string(image->getNumRows()) + "x" + std::to_string(image->getNumColumns()) + ":";
    Fnv1a64 hash;
    hash.appendBytes(reinterpret_cast<const unsigned char*>(header.data()), header.size());
    for (int index = 0; index < image->getSize(); ++index) {
        hash.append((*image)[index]);
    }
    return hash.value();
}

template <typename Value>
void appendExportedHierarchy(Fnv1a64& hash, const std::pair<std::vector<NodeId>, std::vector<Value>>& hierarchy) noexcept {
    hash.appendVector(hierarchy.first);
    hash.appendVector(hierarchy.second);
}

template <AltitudeValue T> [[nodiscard]] std::uint64_t valuedTreeChecksum(const ValuedMorphologicalTree<T>& tree) {
    Fnv1a64 hash;
    hash.append(tree.topology().numRows());
    hash.append(tree.topology().numColumns());
    hash.append(tree.topology().numInternalNodeSlots());
    hash.append(tree.topology().numNodes());
    hash.append(tree.topology().root());
    hash.append(tree.topology().numPixels());
    for (PixelId pixel = 0; pixel < tree.topology().numPixels(); ++pixel) {
        hash.append(tree.topology().smallestNode(pixel));
    }
    appendExportedHierarchy(hash, tree.exportHigraHierarchy());
    return hash.value();
}

template <AltitudeValue T> [[nodiscard]] WorkloadMetrics metricsOf(const ValuedMorphologicalTree<T>& tree) {
    const MorphologicalTree& topology = tree.topology();
    return WorkloadMetrics{
        .pixels = static_cast<std::int64_t>(topology.numRows()) * topology.numColumns(),
        .properParts = topology.numPixels(),
        .primaryNodeSlots = topology.numInternalNodeSlots(),
        .primaryLiveNodes = topology.numNodes(),
    };
}

template <class Map> [[nodiscard]] std::uint64_t edgeMapChecksum(const Map& map) noexcept {
    Fnv1a64 hash;
    hash.append(map.numRows);
    hash.append(map.numColumns);
    hash.append(map.adjacencyRadius);
    hash.appendVector(map.sources);
    hash.appendVector(map.targets);
    hash.appendVector(map.values);
    return hash.value();
}

template <std::floating_point Real>
[[nodiscard]] std::uint64_t semanticAttributeChecksum(std::span<const ComputedAttributeData<Real>* const> results) {
    struct Column {
        Attribute attribute;
        const ComputedAttributeData<Real>* result;
        int offset;
    };

    if (results.empty()) {
        throw std::invalid_argument("Semantic attribute checksum requires at least one result.");
    }

    std::vector<Column> columns;
    std::size_t rowCount = 0;
    NodeIdSpace nodeIdSpace = results.front()->nodeIdSpace;
    for (const ComputedAttributeData<Real>* result : results) {
        if (result == nullptr || result->first.NUM_ATTRIBUTES <= 0 || result->second.size() % static_cast<std::size_t>(result->first.NUM_ATTRIBUTES) != 0) {
            throw std::invalid_argument("Semantic attribute checksum received an invalid attribute layout.");
        }
        const std::size_t currentRows = result->second.size() / static_cast<std::size_t>(result->first.NUM_ATTRIBUTES);
        if (rowCount == 0) {
            rowCount = currentRows;
        } else if (currentRows != rowCount) {
            throw std::invalid_argument("Semantic attribute checksum requires a common node domain.");
        }
        if (result->nodeIdSpace != nodeIdSpace) {
            throw std::invalid_argument("Semantic attribute checksum requires a common node-id space.");
        }
        for (const auto& [attribute, offset] : result->first.indexMap) {
            columns.push_back(Column{attribute, result, offset});
        }
    }

    std::sort(columns.begin(), columns.end(), [](const Column& lhs, const Column& rhs) {
        return static_cast<int>(lhs.attribute) < static_cast<int>(rhs.attribute);
    });
    for (std::size_t index = 1; index < columns.size(); ++index) {
        if (columns[index - 1].attribute == columns[index].attribute) {
            throw std::invalid_argument("Semantic attribute checksum received a duplicate scalar attribute.");
        }
    }

    Fnv1a64 hash;
    hash.append(nodeIdSpace);
    hash.append(rowCount);
    hash.append(columns.size());
    for (const Column& column : columns) {
        hash.append(column.attribute);
        const int stride = column.result->first.NUM_ATTRIBUTES;
        for (std::size_t row = 0; row < rowCount; ++row) {
            hash.append(column.result->second[row * static_cast<std::size_t>(stride) + static_cast<std::size_t>(column.offset)]);
        }
    }
    return hash.value();
}

template <std::floating_point Real> [[nodiscard]] std::uint64_t semanticAttributeChecksum(const ComputedAttributeData<Real>& result) {
    const std::array<const ComputedAttributeData<Real>*, 1> results{&result};
    return semanticAttributeChecksum<Real>(std::span<const ComputedAttributeData<Real>* const>(results));
}

template <std::floating_point Real> [[nodiscard]] std::uint64_t sampledNodeAttributeChecksum(const SampledNodeAttributeData<Real>& result) {
    struct Column {
        NodeAttributeSampleKey key;
        int offset;
    };

    if (result.first.NUM_ATTRIBUTES <= 0 || result.second.size() % static_cast<std::size_t>(result.first.NUM_ATTRIBUTES) != 0) {
        throw std::invalid_argument("Sampled node-attribute checksum received an invalid layout.");
    }
    std::vector<Column> columns;
    columns.reserve(result.first.indexMap.size());
    for (const auto& [key, offset] : result.first.indexMap) {
        columns.push_back(Column{key, offset});
    }
    std::sort(columns.begin(), columns.end(), [](const Column& lhs, const Column& rhs) {
        if (lhs.key.sampleOffset != rhs.key.sampleOffset) {
            return lhs.key.sampleOffset < rhs.key.sampleOffset;
        }
        return static_cast<int>(lhs.key.attribute) < static_cast<int>(rhs.key.attribute);
    });

    const std::size_t rowCount = result.second.size() / static_cast<std::size_t>(result.first.NUM_ATTRIBUTES);
    Fnv1a64 hash;
    hash.append(result.nodeIdSpace);
    hash.append(rowCount);
    hash.append(columns.size());
    for (const Column& column : columns) {
        hash.append(column.key.attribute);
        hash.append(column.key.sampleOffset);
        for (std::size_t row = 0; row < rowCount; ++row) {
            hash.append(result.second[row * static_cast<std::size_t>(result.first.NUM_ATTRIBUTES) + static_cast<std::size_t>(column.offset)]);
        }
    }
    return hash.value();
}

[[nodiscard]] inline double median(std::vector<double> values) {
    if (values.empty()) {
        throw std::invalid_argument("Scientific benchmark median requires at least one sample.");
    }
    std::sort(values.begin(), values.end());
    const std::size_t middle = values.size() / 2;
    if ((values.size() & 1U) != 0U) {
        return values[middle];
    }
    return (values[middle - 1] + values[middle]) / 2.0;
}

[[nodiscard]] inline double medianAbsoluteDeviation(const std::vector<double>& values, double center) {
    std::vector<double> deviations;
    deviations.reserve(values.size());
    for (double value : values) {
        deviations.push_back(std::abs(value - center));
    }
    return median(std::move(deviations));
}

template <class Operation, class Checksum>
[[nodiscard]] ScenarioResult benchmarkScenario(std::string suite, std::string name, TimingScope scope, int repetitions, WorkloadMetrics metrics,
                                               Operation&& operation, Checksum&& checksumOf) {
    if (repetitions <= 0) {
        throw std::invalid_argument("Scientific benchmark repetitions must be positive.");
    }

    auto warmup = operation();
    const std::uint64_t expectedChecksum = checksumOf(warmup);

    std::vector<double> samples;
    samples.reserve(static_cast<std::size_t>(repetitions));
    for (int repetition = 0; repetition < repetitions; ++repetition) {
        const auto start = BenchmarkClock::now();
        auto result = operation();
        const auto finish = BenchmarkClock::now();
        samples.push_back(std::chrono::duration<double, std::milli>(finish - start).count());

        const std::uint64_t checksum = checksumOf(result);
        if (checksum != expectedChecksum) {
            throw std::runtime_error(suite + "." + name + " produced a non-deterministic result.");
        }
    }

    const double center = median(samples);
    const auto [minimum, maximum] = std::minmax_element(samples.begin(), samples.end());
    ScenarioResult result{std::move(suite), std::move(name), scope, center, medianAbsoluteDeviation(samples, center), *minimum, *maximum,
                          expectedChecksum, metrics};
    result.samplesMilliseconds = std::move(samples);
    return result;
}

template <class Prepare, class Operation, class Checksum>
[[nodiscard]] ScenarioResult benchmarkPreparedScenario(std::string suite, std::string name, int repetitions, WorkloadMetrics metrics, Prepare&& prepare,
                                                       Operation&& operation, Checksum&& checksumOf) {
    if (repetitions <= 0) {
        throw std::invalid_argument("Scientific benchmark repetitions must be positive.");
    }

    auto warmState = prepare();
    auto warmResult = operation(warmState);
    const std::uint64_t expectedChecksum = checksumOf(warmState, warmResult);

    std::vector<double> samples;
    samples.reserve(static_cast<std::size_t>(repetitions));
    for (int repetition = 0; repetition < repetitions; ++repetition) {
        auto state = prepare();
        const auto start = BenchmarkClock::now();
        auto result = operation(state);
        const auto finish = BenchmarkClock::now();
        samples.push_back(std::chrono::duration<double, std::milli>(finish - start).count());

        const std::uint64_t checksum = checksumOf(state, result);
        if (checksum != expectedChecksum) {
            throw std::runtime_error(suite + "." + name + " produced a non-deterministic result.");
        }
    }

    const double center = median(samples);
    const auto [minimum, maximum] = std::minmax_element(samples.begin(), samples.end());
    ScenarioResult result{std::move(suite), std::move(name), TimingScope::EstablishedInput, center, medianAbsoluteDeviation(samples, center), *minimum, *maximum,
                          expectedChecksum, metrics};
    result.samplesMilliseconds = std::move(samples);
    return result;
}

[[nodiscard]] inline std::string_view timingScopeName(TimingScope scope) noexcept {
    return scope == TimingScope::EndToEnd ? "end_to_end" : "established_input";
}

[[nodiscard]] inline std::string escapeJson(std::string_view value) {
    std::string escaped;
    escaped.reserve(value.size());
    for (char character : value) {
        switch (character) {
        case '\\':
            escaped += "\\\\";
            break;
        case '"':
            escaped += "\\\"";
            break;
        case '\n':
            escaped += "\\n";
            break;
        case '\r':
            escaped += "\\r";
            break;
        case '\t':
            escaped += "\\t";
            break;
        default:
            escaped += character;
        }
    }
    return escaped;
}

using BenchmarkMetadata = std::vector<std::pair<std::string, std::string>>;

inline void emitKeyValue(const BenchmarkMetadata& metadata, const std::vector<ScenarioResult>& results, bool emitSamples = false,
                         std::ostream& output = std::cout) {
    for (const auto& [key, value] : metadata) {
        output << key << '=' << value << '\n';
    }
    output << "scenario_count=" << results.size() << '\n' << std::fixed << std::setprecision(6);
    for (const ScenarioResult& result : results) {
        const std::string name = result.fullName();
        output << name << "_median_ms=" << result.medianMilliseconds << '\n'
               << name << "_mad_ms=" << result.madMilliseconds << '\n'
               << name << "_minimum_ms=" << result.minimumMilliseconds << '\n'
               << name << "_maximum_ms=" << result.maximumMilliseconds << '\n'
               << name << "_checksum=" << result.checksum << '\n'
               << name << "_pixels=" << result.metrics.pixels << '\n'
               << name << "_proper_parts=" << result.metrics.properParts << '\n'
               << name << "_primary_node_slots=" << result.metrics.primaryNodeSlots << '\n'
               << name << "_primary_live_nodes=" << result.metrics.primaryLiveNodes << '\n'
               << name << "_secondary_node_slots=" << result.metrics.secondaryNodeSlots << '\n'
               << name << "_secondary_live_nodes=" << result.metrics.secondaryLiveNodes << '\n'
               << name << "_edges=" << result.metrics.edges << '\n'
               << name << "_steps=" << result.outcome.steps << '\n'
               << name << "_primary_nodes_removed=" << result.outcome.primaryNodesRemoved << '\n'
               << name << "_secondary_nodes_removed=" << result.outcome.secondaryNodesRemoved << '\n'
               << name << "_complete_validation_commits=" << result.outcome.completeValidationCommits << '\n'
               << name << "_incremental_validation_commits=" << result.outcome.incrementalValidationCommits << '\n';
        output << name << "_light_threshold=" << result.outcome.lightThreshold << '\n'
               << name << "_medium_threshold=" << result.outcome.mediumThreshold << '\n'
               << name << "_heavy_threshold=" << result.outcome.heavyThreshold << '\n';
        if (emitSamples) {
            output << name << "_samples_ms=" << std::setprecision(9);
            for (std::size_t index = 0; index < result.samplesMilliseconds.size(); ++index) {
                if (index != 0) {
                    output << ',';
                }
                output << result.samplesMilliseconds[index];
            }
            output << '\n' << std::setprecision(6);
        }
    }
}

inline void emitJsonLines(const BenchmarkMetadata& metadata, const std::vector<ScenarioResult>& results, bool emitSamples = false,
                          std::ostream& output = std::cout) {
    output << "{\"record\":\"metadata\"";
    for (const auto& [key, value] : metadata) {
        output << ",\"" << escapeJson(key) << "\":\"" << escapeJson(value) << '"';
    }
    output << ",\"scenario_count\":" << results.size() << "}\n" << std::fixed << std::setprecision(6);
    for (const ScenarioResult& result : results) {
        output << "{\"record\":\"scenario\",\"name\":\"" << escapeJson(result.fullName()) << "\",\"suite\":\"" << escapeJson(result.suite)
               << "\",\"scope\":\"" << timingScopeName(result.scope) << "\",\"median_ms\":" << result.medianMilliseconds
               << ",\"mad_ms\":" << result.madMilliseconds << ",\"minimum_ms\":" << result.minimumMilliseconds << ",\"maximum_ms\":"
               << result.maximumMilliseconds << ",\"checksum\":\"" << result.checksum << "\",\"pixels\":" << result.metrics.pixels
               << ",\"proper_parts\":" << result.metrics.properParts << ",\"primary_node_slots\":" << result.metrics.primaryNodeSlots
               << ",\"primary_live_nodes\":" << result.metrics.primaryLiveNodes << ",\"secondary_node_slots\":" << result.metrics.secondaryNodeSlots
               << ",\"secondary_live_nodes\":" << result.metrics.secondaryLiveNodes << ",\"edges\":" << result.metrics.edges
               << ",\"steps\":" << result.outcome.steps << ",\"primary_nodes_removed\":" << result.outcome.primaryNodesRemoved
               << ",\"secondary_nodes_removed\":" << result.outcome.secondaryNodesRemoved
               << ",\"complete_validation_commits\":" << result.outcome.completeValidationCommits
               << ",\"incremental_validation_commits\":" << result.outcome.incrementalValidationCommits
               << ",\"light_threshold\":" << result.outcome.lightThreshold << ",\"medium_threshold\":" << result.outcome.mediumThreshold
               << ",\"heavy_threshold\":" << result.outcome.heavyThreshold;
        if (emitSamples) {
            output << ",\"samples_ms\":[" << std::setprecision(9);
            for (std::size_t index = 0; index < result.samplesMilliseconds.size(); ++index) {
                if (index != 0) {
                    output << ',';
                }
                output << result.samplesMilliseconds[index];
            }
            output << ']' << std::setprecision(6);
        }
        output << "}\n";
    }
}

[[nodiscard]] inline std::string contractModeName() { return contract::validationsEnabled ? "CHECKED" : "UNCHECKED"; }

} // namespace mmcfilters::benchmarks
