#include "ApiBenchmark.hpp"

#include <cstdlib>
#include <iostream>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace mmcfilters::benchmarks::api {
namespace {

[[nodiscard]] Profile parseProfile(std::string_view value) {
    if (value == "smoke") {
        return Profile::Smoke;
    }
    if (value == "core") {
        return Profile::Core;
    }
    if (value == "publication") {
        return Profile::Publication;
    }
    throw std::invalid_argument("Unknown profile: " + std::string(value));
}

[[nodiscard]] InputPattern parseInputPattern(std::string_view value) {
    if (value == "structured") {
        return InputPattern::Structured;
    }
    if (value == "noise") {
        return InputPattern::Noise;
    }
    if (value == "ramp") {
        return InputPattern::Ramp;
    }
    if (value == "geometric") {
        return InputPattern::Geometric;
    }
    if (value == "flat") {
        return InputPattern::Flat;
    }
    throw std::invalid_argument("Unknown input pattern: " + std::string(value));
}

[[nodiscard]] OutputFormat parseOutputFormat(std::string_view value) {
    if (value == "key-value") {
        return OutputFormat::KeyValue;
    }
    if (value == "jsonl") {
        return OutputFormat::JsonLines;
    }
    throw std::invalid_argument("Unknown output format: " + std::string(value));
}

[[nodiscard]] std::vector<double> parseCasfQuantiles(std::string_view value) {
    std::stringstream stream{std::string(value)};
    std::vector<double> quantiles;
    std::string item;
    while (std::getline(stream, item, ',')) {
        std::size_t consumed = 0;
        const double quantile = std::stod(item, &consumed);
        if (consumed != item.size() || !std::isfinite(quantile) || quantile <= 0.0 || quantile >= 1.0) {
            throw std::invalid_argument("--casf-quantiles values must be finite and strictly between zero and one.");
        }
        quantiles.push_back(quantile);
    }
    if (quantiles.size() != 3 || !std::is_sorted(quantiles.begin(), quantiles.end()) ||
        std::adjacent_find(quantiles.begin(), quantiles.end()) != quantiles.end()) {
        throw std::invalid_argument("--casf-quantiles requires three strictly increasing values.");
    }
    return quantiles;
}

[[nodiscard]] int parsePositiveInteger(std::string_view option, std::string_view value) {
    std::size_t consumed = 0;
    const int result = std::stoi(std::string(value), &consumed);
    if (consumed != value.size() || result <= 0) {
        throw std::invalid_argument(std::string(option) + " requires a positive integer.");
    }
    return result;
}

void applyProfileDefaults(Options& options, bool rowsExplicit, bool colsExplicit, bool repetitionsExplicit) {
    int rows = 48;
    int cols = 48;
    int repetitions = 2;
    if (options.profile == Profile::Core) {
        rows = 192;
        cols = 192;
        repetitions = 5;
    } else if (options.profile == Profile::Publication) {
        rows = 512;
        cols = 512;
        repetitions = 15;
    }
    if (!rowsExplicit) {
        options.rows = rows;
    }
    if (!colsExplicit) {
        options.cols = cols;
    }
    if (!repetitionsExplicit) {
        options.repetitions = repetitions;
    }
}

void addSuites(Options& options, std::string_view value, bool& suitesExplicit) {
    if (!suitesExplicit) {
        options.suites.clear();
        suitesExplicit = true;
    }
    std::stringstream stream{std::string(value)};
    std::string suite;
    while (std::getline(stream, suite, ',')) {
        if (suite.empty()) {
            throw std::invalid_argument("--suite contains an empty suite name.");
        }
        options.suites.insert(std::move(suite));
    }
}

void printHelp(std::ostream& output) {
    output << "Usage: mmcfilters_api_benchmark [options]\n"
              "  --profile smoke|core|publication\n"
              "  --suite all|construction|attributes|filters|interoperability|editing|casf|pipelines\n"
              "          May be repeated or receive a comma-separated list.\n"
              "  --input structured|noise|ramp|geometric|flat\n"
              "  --input-file PATH       Load a real grayscale image; rows and cols are verified.\n"
              "  --manifest PATH --workload NAME\n"
              "  --casf-quantiles Q1,Q2,Q3\n"
              "  --rows N --cols N --repetitions N\n"
              "  --emit-samples           Include every timed repetition in the output.\n"
              "  --format key-value|jsonl\n";
}

[[nodiscard]] Options parseOptions(int argc, char** argv) {
    Options options;
    std::optional<std::filesystem::path> manifestPath;
    std::optional<std::string> workloadName;
    for (int index = 1; index < argc; ++index) {
        const std::string_view argument = argv[index];
        if ((argument == "--manifest" || argument == "--workload") && index + 1 >= argc) {
            throw std::invalid_argument(std::string(argument) + " requires a value.");
        }
        if (argument == "--manifest") {
            manifestPath = argv[++index];
        } else if (argument == "--workload") {
            workloadName = argv[++index];
        }
    }
    if (manifestPath.has_value() != workloadName.has_value()) {
        throw std::invalid_argument("--manifest and --workload must be supplied together.");
    }
    ManifestAppliedFields manifestFields;
    if (manifestPath) {
        manifestFields = applyWorkloadManifest(options, *manifestPath, *workloadName);
    }

    bool rowsExplicit = manifestFields.rows;
    bool colsExplicit = manifestFields.cols;
    bool repetitionsExplicit = manifestFields.repetitions;
    bool suitesExplicit = false;

    for (int index = 1; index < argc; ++index) {
        const std::string_view argument = argv[index];
        const auto requireValue = [&](std::string_view option) -> std::string_view {
            if (++index >= argc) {
                throw std::invalid_argument(std::string(option) + " requires a value.");
            }
            return argv[index];
        };

        if (argument == "--help" || argument == "-h") {
            printHelp(std::cout);
            std::exit(0);
        }
        if (argument == "--profile") {
            options.profile = parseProfile(requireValue(argument));
            applyProfileDefaults(options, rowsExplicit, colsExplicit, repetitionsExplicit);
        } else if (argument == "--suite") {
            addSuites(options, requireValue(argument), suitesExplicit);
        } else if (argument == "--input") {
            options.inputPattern = parseInputPattern(requireValue(argument));
            options.inputSource = InputSource::Synthetic;
            options.inputPath.clear();
        } else if (argument == "--input-file") {
            options.inputSource = InputSource::File;
            options.inputPath = std::filesystem::absolute(std::filesystem::path(requireValue(argument))).lexically_normal();
        } else if (argument == "--manifest" || argument == "--workload") {
            static_cast<void>(requireValue(argument));
        } else if (argument == "--casf-quantiles") {
            options.casfQuantiles = parseCasfQuantiles(requireValue(argument));
        } else if (argument == "--rows") {
            options.rows = parsePositiveInteger(argument, requireValue(argument));
            rowsExplicit = true;
        } else if (argument == "--cols") {
            options.cols = parsePositiveInteger(argument, requireValue(argument));
            colsExplicit = true;
        } else if (argument == "--repetitions") {
            options.repetitions = parsePositiveInteger(argument, requireValue(argument));
            repetitionsExplicit = true;
        } else if (argument == "--format") {
            options.outputFormat = parseOutputFormat(requireValue(argument));
        } else if (argument == "--emit-samples") {
            options.emitSamples = true;
        } else {
            throw std::invalid_argument("Unknown option: " + std::string(argument));
        }
    }

    const std::set<std::string> allowedSuites{"all", "construction", "attributes", "filters", "interoperability", "editing", "casf", "pipelines"};
    for (const std::string& suite : options.suites) {
        if (!allowedSuites.contains(suite)) {
            throw std::invalid_argument("Unknown suite: " + suite);
        }
    }
    return options;
}

} // namespace

int run(int argc, char** argv) {
    Options options = parseOptions(argc, argv);
    Context context(options);
    std::vector<ScenarioResult> results;

    if (includesSuite(options, "construction")) {
        addConstructionScenarios(context, results);
    }
    if (includesSuite(options, "attributes")) {
        addAttributeScenarios(context, results);
    }
    if (includesSuite(options, "filters")) {
        addFilterScenarios(context, results);
    }
    if (includesSuite(options, "interoperability")) {
        addInteroperabilityScenarios(context, results);
    }
    if (includesSuite(options, "editing")) {
        addEditingScenarios(context, results);
    }
    if (includesSuite(options, "casf")) {
        addCasfScenarios(context, results);
    }
    if (includesSuite(options, "pipelines")) {
        addPipelineScenarios(context, results);
    }
    if (results.empty()) {
        throw std::runtime_error("The selected benchmark configuration produced no scenarios.");
    }

    std::string quantiles;
    for (double quantile : options.casfQuantiles) {
        if (!quantiles.empty()) {
            quantiles += ',';
        }
        quantiles += std::to_string(quantile);
    }
    std::string bundleNames;
    for (const AttributeBundle& bundle : options.attributeBundles) {
        if (!bundleNames.empty()) {
            bundleNames += ',';
        }
        bundleNames += bundle.name;
    }
    const BenchmarkMetadata metadata{{"benchmark", "scientific_api"},
                                     {"contract_mode", contractModeName()},
                                     {"profile", std::string(profileName(options.profile))},
                                     {"input", options.inputSource == InputSource::Synthetic ? std::string(inputPatternName(options.inputPattern))
                                                                                             : options.inputPath.filename().string()},
                                     {"input_source", inputSourceName(options)},
                                     {"input_checksum", std::to_string(benchmarkInputChecksum(context.imageUInt8))},
                                     {"workload", options.workloadName},
                                     {"manifest", options.manifestPath.empty() ? "none" : options.manifestPath.string()},
                                     {"casf_quantiles", quantiles},
                                     {"attribute_bundles", bundleNames.empty() ? "none" : bundleNames},
                                     {"suite", suiteSelectionName(options)},
                                     {"rows", std::to_string(options.rows)},
                                     {"cols", std::to_string(options.cols)},
                                     {"repetitions", std::to_string(options.repetitions)},
                                     {"samples_emitted", options.emitSamples ? "true" : "false"},
                                     {"cplusplus", std::to_string(__cplusplus)},
                                     {"compiler", __VERSION__}};
    if (options.outputFormat == OutputFormat::JsonLines) {
        emitJsonLines(metadata, results, options.emitSamples);
    } else {
        emitKeyValue(metadata, results, options.emitSamples);
    }
    return 0;
}

} // namespace mmcfilters::benchmarks::api

int main(int argc, char** argv) {
    try {
        return mmcfilters::benchmarks::api::run(argc, argv);
    } catch (const std::exception& error) {
        std::cerr << "API benchmark failed: " << error.what() << '\n';
        return 1;
    }
}
