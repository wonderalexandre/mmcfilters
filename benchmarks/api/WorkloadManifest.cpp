#include "ApiBenchmark.hpp"

#include "stb_image.h"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <fstream>
#include <memory>
#include <map>
#include <sstream>

namespace mmcfilters::benchmarks::api {
namespace {

[[nodiscard]] std::string trim(std::string_view value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string_view::npos) {
        return {};
    }
    const auto last = value.find_last_not_of(" \t\r\n");
    return std::string(value.substr(first, last - first + 1));
}

[[nodiscard]] std::string uppercase(std::string_view value) {
    std::string result(value);
    std::transform(result.begin(), result.end(), result.begin(), [](unsigned char character) { return static_cast<char>(std::toupper(character)); });
    return result;
}

[[nodiscard]] int positiveInteger(std::string_view key, std::string_view value) {
    int result = 0;
    const auto [end, error] = std::from_chars(value.data(), value.data() + value.size(), result);
    if (error != std::errc{} || end != value.data() + value.size() || result <= 0) {
        throw std::invalid_argument("Manifest field " + std::string(key) + " requires a positive integer.");
    }
    return result;
}

[[nodiscard]] std::uint64_t unsignedInteger(std::string_view key, std::string_view value) {
    std::uint64_t result = 0;
    const auto [end, error] = std::from_chars(value.data(), value.data() + value.size(), result);
    if (error != std::errc{} || end != value.data() + value.size()) {
        throw std::invalid_argument("Manifest field " + std::string(key) + " requires an unsigned decimal integer.");
    }
    return result;
}

[[nodiscard]] std::vector<std::string> commaSeparated(std::string_view value) {
    std::stringstream stream{std::string(value)};
    std::vector<std::string> result;
    std::string item;
    while (std::getline(stream, item, ',')) {
        item = trim(item);
        if (item.empty()) {
            throw std::invalid_argument("Manifest lists cannot contain empty values.");
        }
        result.push_back(std::move(item));
    }
    if (result.empty()) {
        throw std::invalid_argument("Manifest lists cannot be empty.");
    }
    return result;
}

[[nodiscard]] AttributeGroup parseGroup(std::string_view value) {
    const std::string name = uppercase(value);
    if (name == "GRAY_LEVEL") {
        return AttributeGroup::GrayLevel;
    }
    if (name == "SHAPE") {
        return AttributeGroup::Shape;
    }
    if (name == "MOMENTS") {
        return AttributeGroup::Moments;
    }
    if (name == "BOUNDARY") {
        return AttributeGroup::Boundary;
    }
    if (name == "TREE_TOPOLOGY") {
        return AttributeGroup::TreeTopology;
    }
    if (name == "ALL") {
        return AttributeGroup::All;
    }
    throw std::invalid_argument("Unknown attribute group in workload manifest: " + std::string(value));
}

[[nodiscard]] Attribute parseAttribute(std::string_view value) {
    const std::string name = uppercase(value);
    for (const attributes::registry::AttributeMetadata& metadata : attributes::registry::ATTRIBUTE_METADATA) {
        if (metadata.name == name) {
            return metadata.attribute;
        }
    }
    throw std::invalid_argument("Unknown scalar attribute in workload manifest: " + std::string(value));
}

[[nodiscard]] AttributeOrGroup parseAttributeRequest(std::string_view value) {
    constexpr std::string_view groupPrefix = "group:";
    constexpr std::string_view attributePrefix = "attribute:";
    if (value.starts_with(groupPrefix)) {
        return parseGroup(value.substr(groupPrefix.size()));
    }
    if (value.starts_with(attributePrefix)) {
        return parseAttribute(value.substr(attributePrefix.size()));
    }
    return parseAttribute(value);
}

void setProfileDefaults(Options& options, std::string_view value) {
    if (value == "smoke") {
        options.profile = Profile::Smoke;
        options.rows = 48;
        options.columns = 48;
        options.repetitions = 2;
    } else if (value == "core") {
        options.profile = Profile::Core;
        options.rows = 192;
        options.columns = 192;
        options.repetitions = 5;
    } else if (value == "publication") {
        options.profile = Profile::Publication;
        options.rows = 512;
        options.columns = 512;
        options.repetitions = 15;
    } else {
        throw std::invalid_argument("Unknown profile in workload manifest: " + std::string(value));
    }
}

void setSyntheticPattern(Options& options, std::string_view value) {
    if (value == "structured") {
        options.inputPattern = InputPattern::Structured;
    } else if (value == "noise") {
        options.inputPattern = InputPattern::Noise;
    } else if (value == "ramp") {
        options.inputPattern = InputPattern::Ramp;
    } else if (value == "geometric") {
        options.inputPattern = InputPattern::Geometric;
    } else if (value == "flat") {
        options.inputPattern = InputPattern::Flat;
    } else {
        throw std::invalid_argument("Unknown synthetic input in workload manifest: " + std::string(value));
    }
    options.inputSource = InputSource::Synthetic;
    options.inputPath.clear();
}

[[nodiscard]] std::vector<double> parseQuantiles(std::string_view value) {
    const std::vector<std::string> values = commaSeparated(value);
    if (values.size() != 3) {
        throw std::invalid_argument("casf_quantiles must contain exactly three values for light, medium, and heavy editing.");
    }
    std::vector<double> quantiles;
    quantiles.reserve(values.size());
    for (const std::string& item : values) {
        std::size_t consumed = 0;
        const double quantile = std::stod(item, &consumed);
        if (consumed != item.size() || !std::isfinite(quantile) || quantile <= 0.0 || quantile >= 1.0) {
            throw std::invalid_argument("CASF quantiles must be finite values strictly between zero and one.");
        }
        quantiles.push_back(quantile);
    }
    if (!std::is_sorted(quantiles.begin(), quantiles.end()) || std::adjacent_find(quantiles.begin(), quantiles.end()) != quantiles.end()) {
        throw std::invalid_argument("CASF quantiles must be strictly increasing.");
    }
    return quantiles;
}

[[nodiscard]] std::set<std::string> parseSuites(std::string_view value) {
    std::set<std::string> suites;
    for (std::string suite : commaSeparated(value)) {
        suites.insert(std::move(suite));
    }
    return suites;
}

} // namespace

ManifestAppliedFields applyWorkloadManifest(Options& options, const std::filesystem::path& manifestPath, std::string_view workloadName) {
    const std::filesystem::path absoluteManifest = std::filesystem::absolute(manifestPath).lexically_normal();
    std::ifstream input(absoluteManifest);
    if (!input) {
        throw std::runtime_error("Could not open scientific workload manifest: " + absoluteManifest.string());
    }

    const std::string expectedSection = "workload." + std::string(workloadName);
    std::string activeSection;
    std::map<std::string, std::string> fields;
    int version = 0;
    std::string line;
    int lineNumber = 0;
    while (std::getline(input, line)) {
        ++lineNumber;
        const std::size_t comment = line.find('#');
        if (comment != std::string::npos) {
            line.erase(comment);
        }
        line = trim(line);
        if (line.empty()) {
            continue;
        }
        if (line.front() == '[' && line.back() == ']') {
            activeSection = trim(std::string_view(line).substr(1, line.size() - 2));
            continue;
        }
        const std::size_t separator = line.find('=');
        if (separator == std::string::npos) {
            throw std::invalid_argument("Invalid workload manifest assignment at line " + std::to_string(lineNumber) + ".");
        }
        const std::string key = trim(std::string_view(line).substr(0, separator));
        const std::string value = trim(std::string_view(line).substr(separator + 1));
        if (activeSection.empty() && key == "version") {
            version = positiveInteger(key, value);
        } else if (activeSection == expectedSection) {
            if (!fields.emplace(key, value).second) {
                throw std::invalid_argument("Duplicate field " + key + " in workload " + std::string(workloadName) + ".");
            }
        }
    }
    if (version != 1) {
        throw std::invalid_argument("Scientific workload manifest must declare version=1.");
    }
    if (fields.empty()) {
        throw std::invalid_argument("Workload not found in manifest: " + std::string(workloadName));
    }

    options.manifestPath = absoluteManifest;
    options.workloadName = std::string(workloadName);
    ManifestAppliedFields applied;
    if (const auto item = fields.find("profile"); item != fields.end()) {
        setProfileDefaults(options, item->second);
    }
    if (const auto item = fields.find("rows"); item != fields.end()) {
        options.rows = positiveInteger(item->first, item->second);
        applied.rows = true;
    }
    if (const auto item = fields.find("columns"); item != fields.end()) {
        options.columns = positiveInteger(item->first, item->second);
        applied.columns = true;
    }
    if (const auto item = fields.find("repetitions"); item != fields.end()) {
        options.repetitions = positiveInteger(item->first, item->second);
        applied.repetitions = true;
    }
    if (const auto item = fields.find("input"); item != fields.end()) {
        constexpr std::string_view syntheticPrefix = "synthetic:";
        constexpr std::string_view filePrefix = "file:";
        if (item->second.starts_with(syntheticPrefix)) {
            setSyntheticPattern(options, std::string_view(item->second).substr(syntheticPrefix.size()));
        } else if (item->second.starts_with(filePrefix)) {
            options.inputSource = InputSource::File;
            std::filesystem::path path = trim(std::string_view(item->second).substr(filePrefix.size()));
            if (path.is_relative()) {
                path = absoluteManifest.parent_path() / path;
            }
            options.inputPath = std::filesystem::absolute(path).lexically_normal();
        } else {
            throw std::invalid_argument("Manifest input must start with synthetic: or file:.");
        }
    }
    if (const auto item = fields.find("suites"); item != fields.end()) {
        options.suites = parseSuites(item->second);
    }
    if (const auto item = fields.find("input_checksum"); item != fields.end()) {
        options.expectedInputChecksum = unsignedInteger(item->first, item->second);
    }
    if (const auto item = fields.find("casf_quantiles"); item != fields.end()) {
        options.casfQuantiles = parseQuantiles(item->second);
    }

    for (const auto& [key, value] : fields) {
        constexpr std::string_view prefix = "attribute_bundle.";
        if (!key.starts_with(prefix)) {
            continue;
        }
        const std::string name = key.substr(prefix.size());
        if (name.empty() || !std::all_of(name.begin(), name.end(), [](unsigned char character) { return std::isalnum(character) || character == '_'; })) {
            throw std::invalid_argument("Attribute bundle names may contain only letters, digits, and underscores.");
        }
        AttributeBundle bundle{.name = name};
        for (const std::string& request : commaSeparated(value)) {
            bundle.requests.push_back(parseAttributeRequest(request));
        }
        options.attributeBundles.push_back(std::move(bundle));
    }

    const std::set<std::string> knownFields{"profile", "input", "rows", "columns", "repetitions", "suites", "input_checksum", "casf_quantiles"};
    for (const auto& [key, value] : fields) {
        static_cast<void>(value);
        if (!knownFields.contains(key) && !key.starts_with("attribute_bundle.")) {
            throw std::invalid_argument("Unknown field " + key + " in workload " + std::string(workloadName) + ".");
        }
    }

    if (options.inputSource == InputSource::File && (!applied.rows || !applied.columns)) {
        throw std::invalid_argument("File workloads must declare rows and columns so image dimensions are never resampled implicitly.");
    }
    return applied;
}

ImageUInt8Ptr makeConfiguredInput(const Options& options) {
    ImageUInt8Ptr image;
    if (options.inputSource == InputSource::Synthetic) {
        image = makeInput<std::uint8_t>(options.inputPattern, options.rows, options.columns);
    } else {
        int width = 0;
        int height = 0;
        int channels = 0;
        using Buffer = std::unique_ptr<unsigned char, decltype(&stbi_image_free)>;
        Buffer pixels(stbi_load(options.inputPath.string().c_str(), &width, &height, &channels, 1), &stbi_image_free);
        if (!pixels || width <= 0 || height <= 0) {
            throw std::runtime_error("Could not load workload image: " + options.inputPath.string());
        }
        if (height != options.rows || width != options.columns) {
            throw std::invalid_argument("Workload image dimensions do not match the manifest: loaded " + std::to_string(height) + "x" +
                                        std::to_string(width) + ", expected " + std::to_string(options.rows) + "x" + std::to_string(options.columns) + ".");
        }
        image = ImageUInt8::create(height, width);
        std::copy_n(pixels.get(), image->getSize(), image->rawData());
    }
    const std::uint64_t checksum = benchmarkInputChecksum(image);
    if (options.expectedInputChecksum && checksum != *options.expectedInputChecksum) {
        throw std::runtime_error("Workload input checksum mismatch: expected " + std::to_string(*options.expectedInputChecksum) + ", received " +
                                 std::to_string(checksum) + ".");
    }
    return image;
}

} // namespace mmcfilters::benchmarks::api
