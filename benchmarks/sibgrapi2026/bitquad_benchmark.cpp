#include "morphotree/adjacency/adjacency8c.hpp"
#include "morphotree/attributes/bitquads/quadCountComputer.hpp"
#include "morphotree/attributes/bitquads/quadCountTreeOfShapesComputer.hpp"
#include "morphotree/core/box.hpp"
#include "morphotree/tree/mtree.hpp"
#include "morphotree/tree/treeOfShapes/kgrid.hpp"
#include "morphotree/tree/treeOfShapes/order_image.hpp"
#include "morphotree/tree/treeOfShapes/tos.hpp"

#include "mmcfilters/attributes/computers/detail/BitquadFiniteWindowComputation.hpp"
#include "mmcfilters/trees/MorphologicalTreeFactory.hpp"
#include "mmcfilters/trees/MorphologicalTreeSemantics.hpp"
#include "mmcfilters/utils/Contract.hpp"

#include "stb_image.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

namespace fs = std::filesystem;
namespace mf = mmcfilters;
namespace mt = morphotree;
namespace bitquad_detail = mmcfilters::attributes::computers::detail;

using Clock = std::chrono::steady_clock;
using Families = bitquad_detail::BitquadFamilyCounts;
using LocalCounter = bitquad_detail::BitquadFiniteWindowComputation;
using ReferenceTree = mt::MorphologicalTree<mt::uint8>;

constexpr int repetitions = 3;
volatile std::uint64_t benchmarkSink = 0;

struct Options {
    fs::path imageDirectory;
    std::string resolution;
    fs::path outputFilename;
    int firstImage = 0;
    int numberOfImages = 100;
    fs::path maximumDecisionTable;
    fs::path minimumDecisionTable;
};

struct GrayImage {
    int width = 0;
    int height = 0;
    std::vector<mt::uint8> pixels;
};

struct DifferenceSummary {
    std::uint64_t mismatchNodes = 0;
    std::array<std::uint64_t, 5> familyMismatchNodes{};
    std::uint64_t maximumAbsoluteFamilyError = 0;
};

struct TimedResult {
    double milliseconds = 0.0;
    std::uint64_t checksum = 0;
};

struct Measurement {
    std::string method;
    int run = 0;
    int orderPosition = 0;
    TimedResult result;
};

struct OutputContext {
    std::string resolution;
    int imageIndex = 0;
    std::string imageName;
    const GrayImage* image = nullptr;
    std::string hierarchy;
    std::string connectivity;
    std::uint32_t nodes = 0;
    std::uint32_t baselineNodes = 0;
    std::uint32_t comparedNodes = 0;
    double treeBuildMilliseconds = 0.0;
    double importMilliseconds = 0.0;
};

[[noreturn]] void usage(const char* program, const std::string& error = {}) {
    if (!error.empty()) {
        std::cerr << "error: " << error << "\n\n";
    }
    std::cerr << "Usage: " << program << " --image-dir DIR --resolution LABEL --output RAW.csv\n"
              << "       --dt-max-8c FILE --dt-min-8c FILE [--start N] [--count N]\n"
              << "The SIBGRAPI 2026 protocol always performs one untimed warm-up and exactly three timed repetitions.\n";
    throw std::invalid_argument("invalid command line");
}

Options parseOptions(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        auto value = [&](const char* name) -> std::string {
            if (++i >= argc) {
                usage(argv[0], std::string(name) + " requires a value");
            }
            return argv[i];
        };
        if (argument == "--image-dir") {
            options.imageDirectory = value("--image-dir");
        } else if (argument == "--resolution") {
            options.resolution = value("--resolution");
        } else if (argument == "--output") {
            options.outputFilename = value("--output");
        } else if (argument == "--start") {
            options.firstImage = std::stoi(value("--start"));
        } else if (argument == "--count") {
            options.numberOfImages = std::stoi(value("--count"));
        } else if (argument == "--dt-max-8c") {
            options.maximumDecisionTable = value("--dt-max-8c");
        } else if (argument == "--dt-min-8c") {
            options.minimumDecisionTable = value("--dt-min-8c");
        } else if (argument == "--help" || argument == "-h") {
            usage(argv[0]);
        } else {
            usage(argv[0], "unknown option: " + argument);
        }
    }

    if (options.imageDirectory.empty() || options.resolution.empty() || options.outputFilename.empty() ||
        options.maximumDecisionTable.empty() || options.minimumDecisionTable.empty()) {
        usage(argv[0], "all required options must be supplied");
    }
    if (!fs::is_directory(options.imageDirectory)) {
        usage(argv[0], "image directory does not exist: " + options.imageDirectory.string());
    }
    if (!fs::is_regular_file(options.maximumDecisionTable) || !fs::is_regular_file(options.minimumDecisionTable)) {
        usage(argv[0], "one or both decision-table files do not exist");
    }
    if (options.firstImage < 0 || options.numberOfImages < 1 || options.firstImage + options.numberOfImages > 100) {
        usage(argv[0], "require START >= 0, COUNT >= 1, and START + COUNT <= 100");
    }
    if (fs::exists(options.outputFilename)) {
        usage(argv[0], "refusing to overwrite existing output: " + options.outputFilename.string());
    }
    return options;
}

GrayImage readGrayImage(const fs::path& filename) {
    int width = 0;
    int height = 0;
    int channels = 0;
    unsigned char* raw = stbi_load(filename.string().c_str(), &width, &height, &channels, 1);
    if (raw == nullptr) {
        throw std::runtime_error("cannot load image " + filename.string() + ": " + stbi_failure_reason());
    }
    if (width <= 0 || height <= 0) {
        stbi_image_free(raw);
        throw std::runtime_error("image has an invalid domain: " + filename.string());
    }

    GrayImage image;
    image.width = width;
    image.height = height;
    image.pixels.assign(raw, raw + static_cast<std::size_t>(width) * static_cast<std::size_t>(height));
    stbi_image_free(raw);
    return image;
}

mt::Box domainOf(const GrayImage& image) {
    return mt::Box::fromSize(mt::I32Point{0, 0}, mt::UI32Point{static_cast<mt::uint32>(image.width), static_cast<mt::uint32>(image.height)});
}

mf::ValuedMorphologicalTree<std::uint8_t> importTree(const ReferenceTree& tree, int rows, int columns, mf::MorphologicalTreeKind kind) {
    const std::size_t numberOfNodes = tree.numberOfNodes();
    const std::size_t numberOfPixels = static_cast<std::size_t>(rows) * static_cast<std::size_t>(columns);
    if (tree.numberOfCNPs() != numberOfPixels) {
        throw std::runtime_error("tree proper-part domain does not match the image");
    }

    std::vector<mf::NodeId> parent(numberOfNodes);
    std::vector<mf::NodeId> owner(numberOfPixels, mf::InvalidNode);
    std::vector<std::uint8_t> altitude(numberOfNodes);
    for (std::size_t node = 0; node < numberOfNodes; ++node) {
        const auto sourceNode = tree.node(static_cast<mt::uint32>(node));
        parent[node] = sourceNode->parent() == nullptr ? static_cast<mf::NodeId>(node) : static_cast<mf::NodeId>(sourceNode->parent()->id());
        altitude[node] = static_cast<std::uint8_t>(sourceNode->level());
    }
    for (std::size_t pixel = 0; pixel < numberOfPixels; ++pixel) {
        owner[pixel] = static_cast<mf::NodeId>(tree.smallComponent(static_cast<mt::uint32>(pixel))->id());
    }

    mf::MorphologicalTreeSemantics semantics;
    if (kind == mf::MorphologicalTreeKind::MaxTree || kind == mf::MorphologicalTreeKind::MinTree) {
        semantics = mf::makeMorphologicalTreeSemantics(kind, mf::SharedAdjacencyContext{mf::RegularGridAdjacency2D(rows, columns, 1.5)});
    } else {
        semantics = mf::makeMorphologicalTreeSemantics(kind);
    }
    return mf::MorphologicalTreeFactory::createFromNativeTopology(
        std::span<const mf::NodeId>(parent), std::span<const mf::NodeId>(owner), std::span<const std::uint8_t>(altitude),
        static_cast<mf::NodeId>(tree.root()->id()), rows, columns, std::move(semantics));
}

std::array<std::int64_t, 5> familyValues(const mt::Quads& quads) {
    return {quads.q1(), quads.q2(), quads.qd(), quads.q3(), quads.q4()};
}

std::array<std::int64_t, 5> familyValues(const Families& families) {
    return {families.q1, families.q2, families.qd, families.q3, families.q4};
}

template <class First, class Second>
DifferenceSummary compareFamilies(const std::vector<First>& first, const std::vector<Second>& second,
                                  std::size_t nodesToCompare = std::numeric_limits<std::size_t>::max()) {
    const std::size_t numberOfNodes = nodesToCompare == std::numeric_limits<std::size_t>::max() ? first.size() : nodesToCompare;
    if (first.size() < numberOfNodes || second.size() < numberOfNodes ||
        (nodesToCompare == std::numeric_limits<std::size_t>::max() && first.size() != second.size())) {
        throw std::runtime_error("counter vectors do not cover the requested node domain");
    }

    DifferenceSummary difference;
    for (std::size_t node = 0; node < numberOfNodes; ++node) {
        const auto firstValues = familyValues(first[node]);
        const auto secondValues = familyValues(second[node]);
        bool nodeMismatch = false;
        for (std::size_t family = 0; family < firstValues.size(); ++family) {
            const std::uint64_t absoluteError = static_cast<std::uint64_t>(
                std::abs(static_cast<long long>(firstValues[family] - secondValues[family])));
            const bool differs = absoluteError != 0;
            difference.familyMismatchNodes[family] += differs;
            nodeMismatch = nodeMismatch || differs;
            difference.maximumAbsoluteFamilyError = std::max(difference.maximumAbsoluteFamilyError, absoluteError);
        }
        difference.mismatchNodes += nodeMismatch;
    }
    return difference;
}

template <class Value> std::uint64_t checksum(const std::vector<Value>& values) {
    std::uint64_t result = values.size();
    for (const Value& value : values) {
        const auto families = familyValues(value);
        result = result * UINT64_C(1099511628211) + static_cast<std::uint32_t>(families[0]) +
                 UINT64_C(3) * static_cast<std::uint32_t>(families[1]) + UINT64_C(5) * static_cast<std::uint32_t>(families[2]) +
                 UINT64_C(7) * static_cast<std::uint32_t>(families[3]) + UINT64_C(11) * static_cast<std::uint32_t>(families[4]);
    }
    return result;
}

template <class Computation> TimedResult measure(Computation&& computation) {
    const Clock::time_point start = Clock::now();
    auto values = computation();
    const Clock::time_point finish = Clock::now();
    TimedResult result{std::chrono::duration<double, std::milli>(finish - start).count(), checksum(values)};
    benchmarkSink ^= result.checksum;
    return result;
}

template <class First, class Second>
std::vector<Measurement> measureAlternatingPair(int orderSeed, const std::string& firstName, First&& first,
                                                const std::string& secondName, Second&& second) {
    std::vector<Measurement> measurements;
    measurements.reserve(static_cast<std::size_t>(repetitions) * 2);
    for (int run = 0; run < repetitions; ++run) {
        const bool firstRunsFirst = ((orderSeed + run) % 2) == 0;
        if (firstRunsFirst) {
            measurements.push_back({firstName, run, 0, first()});
            measurements.push_back({secondName, run, 1, second()});
        } else {
            measurements.push_back({secondName, run, 0, second()});
            measurements.push_back({firstName, run, 1, first()});
        }
    }
    return measurements;
}

void writeHeader(std::ostream& output) {
    output << "resolution,image_index,image,width,height,pixels,hierarchy,connectivity,method,run,order_position,time_ms,"
              "nodes,baseline_nodes,compared_nodes,reference_method,mismatch_nodes,q1_mismatch_nodes,q2_mismatch_nodes,"
              "qd_mismatch_nodes,q3_mismatch_nodes,q4_mismatch_nodes,max_abs_family_error,checksum,tree_build_ms,bridge_import_ms,contract_mode\n";
}

void writeMeasurements(std::ostream& output, const OutputContext& context, const std::vector<Measurement>& measurements,
                       const std::string& baselineMethod, const DifferenceSummary& baselineDifference) {
    const std::uint64_t pixels = static_cast<std::uint64_t>(context.image->width) * static_cast<std::uint64_t>(context.image->height);
    for (const Measurement& measurement : measurements) {
        const DifferenceSummary difference = measurement.method == baselineMethod ? baselineDifference : DifferenceSummary{};
        output << context.resolution << ',' << context.imageIndex << ',' << context.imageName << ',' << context.image->width << ','
               << context.image->height << ',' << pixels << ',' << context.hierarchy << ',' << context.connectivity << ',' << measurement.method << ','
               << measurement.run << ',' << measurement.orderPosition << ',' << std::fixed << std::setprecision(6)
               << measurement.result.milliseconds << ',' << context.nodes << ',' << context.baselineNodes << ',' << context.comparedNodes
               << ",proposed," << difference.mismatchNodes << ',' << difference.familyMismatchNodes[0] << ',' << difference.familyMismatchNodes[1] << ','
               << difference.familyMismatchNodes[2] << ',' << difference.familyMismatchNodes[3] << ',' << difference.familyMismatchNodes[4] << ','
               << difference.maximumAbsoluteFamilyError << ',' << measurement.result.checksum << ',' << context.treeBuildMilliseconds << ','
               << context.importMilliseconds << ',' << (mf::contract::validationsEnabled ? "CHECKED" : "UNCHECKED") << '\n';
    }
    output.flush();
}

void benchmarkComponentTree(std::ostream& output, const Options& options, int imageIndex, const std::string& imageName,
                            const GrayImage& image, bool maximumTree, const fs::path& decisionTable, int orderSeed) {
    const mt::Box domain = domainOf(image);
    const Clock::time_point buildStart = Clock::now();
    const ReferenceTree tree = maximumTree ? mt::buildMaxTree(image.pixels, std::make_shared<mt::Adjacency8C>(domain))
                                           : mt::buildMinTree(image.pixels, std::make_shared<mt::Adjacency8C>(domain));
    const Clock::time_point buildFinish = Clock::now();

    const Clock::time_point importStart = Clock::now();
    const auto imported = importTree(tree, image.height, image.width,
                                     maximumTree ? mf::MorphologicalTreeKind::MaxTree : mf::MorphologicalTreeKind::MinTree);
    const Clock::time_point importFinish = Clock::now();

    mt::CTreeQuadCountsComputer<mt::uint8> reference(domain, image.pixels, decisionTable.string());
    const std::vector<mt::Quads> referenceWarmup = reference.computeAttribute(tree);
    const std::vector<Families> proposedWarmup = LocalCounter::computeBitquadFamilyCounts(imported.topology());
    const DifferenceSummary difference = compareFamilies(referenceWarmup, proposedWarmup);
    benchmarkSink ^= checksum(referenceWarmup) ^ checksum(proposedWarmup);
    if (difference.mismatchNodes != 0) {
        throw std::runtime_error(std::string(maximumTree ? "max-tree" : "min-tree") +
                                 " specialized and proposed counts disagree on " + std::to_string(difference.mismatchNodes) + " nodes");
    }

    const std::vector<Measurement> measurements = measureAlternatingPair(
        orderSeed, "ref6", [&]() { return measure([&]() { return reference.computeAttribute(tree); }); },
        "proposed", [&]() { return measure([&]() { return LocalCounter::computeBitquadFamilyCounts(imported.topology()); }); });

    OutputContext context{options.resolution,
                          imageIndex,
                          imageName,
                          &image,
                          maximumTree ? "max_tree" : "min_tree",
                          "8",
                          static_cast<std::uint32_t>(tree.numberOfNodes()),
                          static_cast<std::uint32_t>(tree.numberOfNodes()),
                          static_cast<std::uint32_t>(tree.numberOfNodes()),
                          std::chrono::duration<double, std::milli>(buildFinish - buildStart).count(),
                          std::chrono::duration<double, std::milli>(importFinish - importStart).count()};
    writeMeasurements(output, context, measurements, "ref6", difference);
}

std::vector<mt::Quads> countWithOriginalReference5(const mt::KGrid<mt::uint8>& kgrid,
                                                   const std::vector<mt::uint32>& orderImage, const ReferenceTree& tree) {
    return mt::TreeOfShapesQuadCountsComputer<mt::uint8>(kgrid, orderImage).computeAttribute(tree);
}

void benchmarkTreeOfShapes(std::ostream& output, const Options& options, int imageIndex, const std::string& imageName,
                           const GrayImage& image, int orderSeed) {
    const mt::Box domain = domainOf(image);
    const Clock::time_point buildStart = Clock::now();
    const mt::KGrid<mt::uint8> kgrid{domain, image.pixels};
    const mt::OrderImageResult<mt::uint8> order = mt::computeOrderImage(domain, image.pixels, kgrid);
    const ReferenceTree enlarged = mt::buildEnlargedTreeOfShapes(order, kgrid);
    const ReferenceTree emerged = mt::emergeTreeOfShapes(kgrid, enlarged);
    const Clock::time_point buildFinish = Clock::now();

    const Clock::time_point importStart = Clock::now();
    const auto imported = importTree(emerged, image.height, image.width, mf::MorphologicalTreeKind::TreeOfShapes);
    const Clock::time_point importFinish = Clock::now();

    const std::vector<mt::Quads> referenceWarmup = countWithOriginalReference5(kgrid, order.orderImg, enlarged);
    const std::vector<Families> proposedWarmup = LocalCounter::computeBitquadFamilyCounts(imported.topology());
    const DifferenceSummary difference = compareFamilies(referenceWarmup, proposedWarmup, emerged.numberOfNodes());
    benchmarkSink ^= checksum(referenceWarmup) ^ checksum(proposedWarmup);

    const std::vector<Measurement> measurements = measureAlternatingPair(
        orderSeed, "ref5_original", [&]() { return measure([&]() { return countWithOriginalReference5(kgrid, order.orderImg, enlarged); }); },
        "proposed", [&]() { return measure([&]() { return LocalCounter::computeBitquadFamilyCounts(imported.topology()); }); });

    OutputContext context{options.resolution,
                          imageIndex,
                          imageName,
                          &image,
                          "tree_of_shapes",
                          "4/8",
                          static_cast<std::uint32_t>(emerged.numberOfNodes()),
                          static_cast<std::uint32_t>(enlarged.numberOfNodes()),
                          static_cast<std::uint32_t>(emerged.numberOfNodes()),
                          std::chrono::duration<double, std::milli>(buildFinish - buildStart).count(),
                          std::chrono::duration<double, std::milli>(importFinish - importStart).count()};
    writeMeasurements(output, context, measurements, "ref5_original", difference);
}

std::string imageFilename(int index) {
    std::ostringstream name;
    name << "val_" << std::setfill('0') << std::setw(3) << index << ".png";
    return name.str();
}

} // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);
        std::ofstream output(options.outputFilename, std::ios::out | std::ios::trunc);
        if (!output) {
            throw std::runtime_error("cannot open output: " + options.outputFilename.string());
        }
        writeHeader(output);

        for (int imageIndex = options.firstImage; imageIndex < options.firstImage + options.numberOfImages; ++imageIndex) {
            const std::string filename = imageFilename(imageIndex);
            const GrayImage image = readGrayImage(options.imageDirectory / filename);
            std::cerr << "begin resolution=" << options.resolution << " image=" << filename << " size=" << image.width << 'x' << image.height << '\n';

            benchmarkComponentTree(output, options, imageIndex, filename, image, true, options.maximumDecisionTable, imageIndex * 3);
            benchmarkComponentTree(output, options, imageIndex, filename, image, false, options.minimumDecisionTable, imageIndex * 3 + 1);
            benchmarkTreeOfShapes(output, options, imageIndex, filename, image, imageIndex * 3 + 2);

            std::cerr << "done resolution=" << options.resolution << " image=" << filename << '\n';
        }
        std::cerr << "repetitions=" << repetitions << " sink=" << benchmarkSink << '\n';
        return 0;
    } catch (const std::invalid_argument& error) {
        if (std::string(error.what()) != "invalid command line") {
            std::cerr << "error: " << error.what() << '\n';
        }
        return 2;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
