#pragma once

#include "../support/ScientificBenchmark.hpp"

#include "mmcfilters/attributes/AttributeComputation.hpp"
#include "mmcfilters/trees/MorphologicalTreeFactory.hpp"
#include "mmcfilters/utils/RegularGridAdjacency2D.hpp"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace mmcfilters::benchmarks::api {

enum class Profile { Smoke, Core, Publication };
enum class InputPattern { Structured, Noise, Ramp, Geometric, Flat };
enum class InputSource { Synthetic, File };
enum class OutputFormat { KeyValue, JsonLines };

/** @brief Named sequence of scalar attributes and groups measured together. */
struct AttributeBundle {
    std::string name;                              ///< Stable bundle name emitted in benchmark records.
    std::vector<AttributeOrGroup> requests;        ///< Ordered scalar attributes or groups requested by the bundle.
};

/** @brief Fully resolved command-line and manifest configuration for an API benchmark run. */
struct Options {
    Profile profile = Profile::Smoke;                              ///< Size and repetition profile.
    InputPattern inputPattern = InputPattern::Structured;          ///< Synthetic input generator.
    InputSource inputSource = InputSource::Synthetic;              ///< Whether input is synthetic or file-backed.
    OutputFormat outputFormat = OutputFormat::KeyValue;            ///< Machine-readable output representation.
    bool emitSamples = false;                                      ///< Whether individual timing samples are emitted.
    std::set<std::string> suites{"all"};                            ///< Scenario suites selected for execution.
    int rows = 48;                                                  ///< Input height in pixels.
    int cols = 48;                                                  ///< Input width in pixels.
    int repetitions = 2;                                           ///< Timed repetitions per scenario.
    std::filesystem::path inputPath;                               ///< File-backed image path, when selected.
    std::filesystem::path manifestPath;                            ///< Workload manifest path, when selected.
    std::string workloadName = "command_line";                     ///< Stable workload identifier.
    std::optional<std::uint64_t> expectedInputChecksum;            ///< Optional checksum pinned by the manifest.
    std::vector<double> casfQuantiles{0.1, 0.5, 0.9};              ///< Quantiles used to derive CASF thresholds.
    std::vector<AttributeBundle> attributeBundles;                 ///< Manifest-defined scientific attribute bundles.
};

/** @brief Dimension and repetition fields explicitly supplied by a workload manifest. */
struct ManifestAppliedFields {
    bool rows = false;        ///< Whether the manifest supplied the row count.
    bool cols = false;        ///< Whether the manifest supplied the column count.
    bool repetitions = false; ///< Whether the manifest supplied the repetition count.
};

[[nodiscard]] inline std::string_view profileName(Profile profile) noexcept {
    switch (profile) {
    case Profile::Smoke:
        return "smoke";
    case Profile::Core:
        return "core";
    case Profile::Publication:
        return "publication";
    }
    return "unknown";
}

[[nodiscard]] inline std::string_view inputPatternName(InputPattern pattern) noexcept {
    switch (pattern) {
    case InputPattern::Structured:
        return "structured";
    case InputPattern::Noise:
        return "noise";
    case InputPattern::Ramp:
        return "ramp";
    case InputPattern::Geometric:
        return "geometric";
    case InputPattern::Flat:
        return "flat";
    }
    return "unknown";
}

[[nodiscard]] inline std::string inputSourceName(const Options& options) {
    if (options.inputSource == InputSource::File) {
        return "file:" + options.inputPath.string();
    }
    return "synthetic:" + std::string(inputPatternName(options.inputPattern));
}

[[nodiscard]] inline bool atLeast(Profile value, Profile minimum) noexcept {
    return static_cast<int>(value) >= static_cast<int>(minimum);
}

[[nodiscard]] inline bool includesSuite(const Options& options, std::string_view suite) {
    return options.suites.contains("all") || options.suites.contains(std::string(suite));
}

[[nodiscard]] inline std::string suiteSelectionName(const Options& options) {
    std::string result;
    for (const std::string& suite : options.suites) {
        if (!result.empty()) {
            result += ',';
        }
        result += suite;
    }
    return result;
}

[[nodiscard]] inline std::uint8_t baseInputValue(InputPattern pattern, int row, int col, int rows, int cols) noexcept {
    const std::uint32_t index = static_cast<std::uint32_t>(row * cols + col);
    switch (pattern) {
    case InputPattern::Structured: {
        const int centeredRow = row - rows / 2;
        const int centeredCol = col - cols / 2;
        const int radial = centeredRow * centeredRow + centeredCol * centeredCol;
        return static_cast<std::uint8_t>((radial + 17 * row + 31 * col) & 0xff);
    }
    case InputPattern::Noise: {
        std::uint32_t value = index + 0x9e3779b9U;
        value ^= value >> 16U;
        value *= 0x7feb352dU;
        value ^= value >> 15U;
        value *= 0x846ca68bU;
        value ^= value >> 16U;
        return static_cast<std::uint8_t>(value & 0xffU);
    }
    case InputPattern::Ramp:
        return static_cast<std::uint8_t>((3 * row + 5 * col) & 0xff);
    case InputPattern::Geometric: {
        std::uint8_t value = 24;
        const int centeredRow = row - rows / 2;
        const int centeredCol = col - cols / 2;
        const int radius = std::max(2, std::min(rows, cols) / 4);
        if (centeredRow * centeredRow + centeredCol * centeredCol <= radius * radius) {
            value = 180;
        }
        if (row >= rows / 8 && row < rows / 3 && col >= cols / 7 && col < cols / 2) {
            value = 92;
        }
        if (row >= rows / 2 && row < 7 * rows / 8 && col >= 5 * cols / 8 && col < 7 * cols / 8) {
            value = 232;
        }
        return value;
    }
    case InputPattern::Flat:
        return 127;
    }
    return 0;
}

template <AltitudeValue T>
[[nodiscard]] inline T typedInputValue(std::uint8_t base, InputPattern pattern, int row, int col, int cols) noexcept {
    if constexpr (std::is_same_v<T, std::uint8_t>) {
        return base;
    } else if constexpr (std::is_integral_v<T>) {
        const int perturbation = pattern == InputPattern::Flat ? 0 : (row * cols + col) % 13;
        return static_cast<T>(static_cast<int>(base) * 257 - 32768 + perturbation);
    } else {
        const T perturbation = pattern == InputPattern::Flat ? T{0} : static_cast<T>((row * cols + col) % 7) / static_cast<T>(1000);
        return static_cast<T>(base) / static_cast<T>(8) + perturbation;
    }
}

template <AltitudeValue T> [[nodiscard]] inline ImagePtr<T> makeInput(InputPattern pattern, int rows, int cols) {
    auto image = Image<T>::create(rows, cols);
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            const std::uint8_t base = baseInputValue(pattern, row, col, rows, cols);
            (*image)[row * cols + col] = typedInputValue<T>(base, pattern, row, col, cols);
        }
    }
    return image;
}

[[nodiscard]] ImageUInt8Ptr makeConfiguredInput(const Options& options);
[[nodiscard]] ManifestAppliedFields applyWorkloadManifest(Options& options, const std::filesystem::path& manifestPath,
                                                          std::string_view workloadName);

template <AltitudeValue T> [[nodiscard]] inline ImagePtr<T> makeTypedConfiguredInput(const Options& options, const ImageUInt8Ptr& baseImage) {
    if (options.inputSource == InputSource::Synthetic) {
        return makeInput<T>(options.inputPattern, options.rows, options.cols);
    }
    auto image = Image<T>::create(options.rows, options.cols);
    for (int index = 0; index < baseImage->getSize(); ++index) {
        (*image)[index] = static_cast<T>((*baseImage)[index]);
    }
    return image;
}

[[nodiscard]] inline std::int64_t countUndirectedEdges(const RegularGridAdjacency2D& adjacency) {
    std::int64_t edges = 0;
    const int size = adjacency.getNumRows() * adjacency.getNumCols();
    for (int pixel = 0; pixel < size; ++pixel) {
        for (int neighbor : adjacency.getNeighborIndices(pixel)) {
            if (neighbor > pixel) {
                ++edges;
            }
        }
    }
    return edges;
}

/** @brief Shared immutable inputs and established structures reused by API benchmark scenarios. */
struct Context {
    Options options;                                  ///< Resolved benchmark configuration.
    ImageUInt8Ptr imageUInt8;                         ///< Canonical 8-bit input image.
    ImageInt32Ptr imageInt32;                         ///< Signed integral form of the input.
    ImageFloatPtr imageFloat;                         ///< Floating-point form of the input.
    RegularGridAdjacency2D adjacency;                 ///< Established regular-grid adjacency.
    WeightedMorphologicalTree<std::uint8_t> maxTree; ///< Established max-tree.
    WeightedMorphologicalTree<std::uint8_t> minTree; ///< Established min-tree.
    ComputedAttributeData<double> area;               ///< Area values established on the max-tree.
    std::vector<bool> keepCriterion;                  ///< Direct-filter criterion derived from area.
    std::vector<float> scores;                        ///< Normalized area scores for ranked filters.
    double areaThreshold = 0.0;                       ///< Area threshold used by representative filters.
    WorkloadMetrics maxTreeMetrics;                   ///< Structural metrics of the established max-tree.

    /** @brief Builds all shared inputs from resolved options. @param benchmarkOptions Resolved benchmark configuration. */
    explicit Context(Options benchmarkOptions)
        : options(std::move(benchmarkOptions)), imageUInt8(makeConfiguredInput(options)), imageInt32(makeTypedConfiguredInput<std::int32_t>(options, imageUInt8)),
          imageFloat(makeTypedConfiguredInput<float>(options, imageUInt8)), adjacency(options.rows, options.cols, 1.5),
          maxTree(MorphologicalTreeFactory::createMaxTree(imageUInt8, adjacency)), minTree(MorphologicalTreeFactory::createMinTree(imageUInt8, adjacency)),
          area(AttributeComputation::computeSingleAttribute<double>(maxTree, AREA)),
          keepCriterion(static_cast<std::size_t>(maxTree.topology().getNumInternalNodeSlots()), false),
          scores(static_cast<std::size_t>(maxTree.topology().getNumInternalNodeSlots()), 0.0F),
          areaThreshold(static_cast<double>(options.rows) * static_cast<double>(options.cols) / 16.0), maxTreeMetrics(metricsOf(maxTree)) {
        maxTreeMetrics.edges = countUndirectedEdges(adjacency);
        const double rootArea = area.second[static_cast<std::size_t>(maxTree.topology().getRoot())];
        for (NodeId node : maxTree.topology().getAliveNodeIds()) {
            const double value = area.second[static_cast<std::size_t>(node)];
            keepCriterion[static_cast<std::size_t>(node)] = value >= areaThreshold;
            scores[static_cast<std::size_t>(node)] = rootArea > 0.0 ? static_cast<float>(value / rootArea) : 0.0F;
        }
        keepCriterion[static_cast<std::size_t>(maxTree.topology().getRoot())] = true;
    }
};

[[nodiscard]] inline std::string lowerName(std::string_view value) {
    std::string result(value);
    std::transform(result.begin(), result.end(), result.begin(), [](unsigned char character) { return static_cast<char>(std::tolower(character)); });
    return result;
}

[[nodiscard]] inline std::string groupName(AttributeGroup group) {
    switch (group) {
    case AttributeGroup::GRAY_LEVEL:
        return "gray_level";
    case AttributeGroup::SHAPE:
        return "shape";
    case AttributeGroup::MOMENTS:
        return "moments";
    case AttributeGroup::BOUNDARY:
        return "boundary";
    case AttributeGroup::TREE_TOPOLOGY:
        return "tree_topology";
    case AttributeGroup::ALL:
        return "all";
    }
    throw std::invalid_argument("Unknown attribute group in API benchmark.");
}

void addConstructionScenarios(Context& context, std::vector<ScenarioResult>& results);
void addAttributeScenarios(Context& context, std::vector<ScenarioResult>& results);
void addFilterScenarios(Context& context, std::vector<ScenarioResult>& results);
void addInteroperabilityScenarios(Context& context, std::vector<ScenarioResult>& results);
void addEditingScenarios(Context& context, std::vector<ScenarioResult>& results);
void addCasfScenarios(Context& context, std::vector<ScenarioResult>& results);
void addPipelineScenarios(Context& context, std::vector<ScenarioResult>& results);

} // namespace mmcfilters::benchmarks::api
