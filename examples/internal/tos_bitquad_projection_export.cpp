/**
 * Internal exporter for inspecting Tree-of-Shapes bitquad projections.
 *
 * Build with `-DMMCFILTERS_BUILD_EXAMPLES=ON` and run:
 * `./build/examples/mmcfilters_tos_bitquad_projection_export --synthetic
 * fixture --out-dir tmp/tos-bitquad-projections`.
 *
 * The program writes CSV files for implementation inspection. It includes
 * `detail/` headers and is not public API.
 */
#include "mmcfilters/attributes/computers/detail/BitquadAttributeData.hpp"
#include "mmcfilters/attributes/computers/detail/BitquadFiniteWindowComputation.hpp"
#include "mmcfilters/trees/MorphologicalTreeFactory.hpp"
#include "mmcfilters/trees/ValuedMorphologicalTree.hpp"
#include "mmcfilters/utils/Image.hpp"
#include "stb_image.h"

#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace mmcfilters;
using namespace mmcfilters::attributes::computers::detail;

namespace {

using FamilyIncrement = mmcfilters::attributes::computers::detail::BitquadFamilyIncrement;
using FamilyCounts = mmcfilters::attributes::computers::detail::BitquadFamilyCounts;
using StateIncrement = BitquadFiniteWindowComputation::NonemptyBitquadStateHistogramIncrement;
using StateHistogram = BitquadFiniteWindowComputation::BitquadStateHistogram;

enum class ExampleImmersion { SelfDualSpan, Min4Max8, Min8Max4 };

TopographicConvention makeConvention(ExampleImmersion immersion, int rows, int columns) {
    if (immersion == ExampleImmersion::SelfDualSpan) {
        return TopographicConvention{};
    }
    const bool minIs4 = immersion == ExampleImmersion::Min4Max8;
    return TopographicConvention{ComplementaryGridImmersion{
        ComplementaryAdjacencies{RegularGridAdjacency2D(rows, columns, minIs4 ? 1.0 : 1.5), RegularGridAdjacency2D(rows, columns, minIs4 ? 1.5 : 1.0)}}};
}

struct Options {
    std::filesystem::path imagePath;
    std::filesystem::path outDir = "tos-bitquad-projections";
    std::string synthetic = "fixture";
    ExampleImmersion immersion = ExampleImmersion::SelfDualSpan;
    int rows = 8;
    int columns = 8;
};

std::string toString(ExampleImmersion immersion) {
    switch (immersion) {
    case ExampleImmersion::SelfDualSpan:
        return "SelfDual";
    case ExampleImmersion::Min4Max8:
        return "Min4cMax8c";
    case ExampleImmersion::Min8Max4:
        return "Min8cMax4c";
    }
    throw std::invalid_argument("Unsupported Tree-of-Shapes interpolation.");
}

ExampleImmersion parseImmersion(const std::string& value) {
    if (value == "SelfDual") {
        return ExampleImmersion::SelfDualSpan;
    }
    if (value == "Min4cMax8c") {
        return ExampleImmersion::Min4Max8;
    }
    if (value == "Min8cMax4c") {
        return ExampleImmersion::Min8Max4;
    }
    throw std::invalid_argument("Expected SelfDual, Min4cMax8c, or Min8cMax4c for --interpolation.");
}

void printUsage(const char* program) {
    std::cerr << "Usage: " << program << " [--image path] [--out-dir dir]\n"
              << "       [--immersion SelfDual|Min4cMax8c|Min8cMax4c]\n"
              << "       [--synthetic fixture|ramp|checker] [--rows n] [--columns n]\n";
}

Options parseOptions(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto requireValue = [&](const char* name) -> std::string {
            if (i + 1 >= argc) {
                throw std::invalid_argument(std::string(name) + " requires a value.");
            }
            return argv[++i];
        };

        if (arg == "--image") {
            options.imagePath = requireValue("--image");
        } else if (arg == "--out-dir") {
            options.outDir = requireValue("--out-dir");
        } else if (arg == "--immersion") {
            options.immersion = parseImmersion(requireValue("--immersion"));
        } else if (arg == "--synthetic") {
            options.synthetic = requireValue("--synthetic");
        } else if (arg == "--rows") {
            options.rows = std::stoi(requireValue("--rows"));
        } else if (arg == "--columns") {
            options.columns = std::stoi(requireValue("--columns"));
        } else if (arg == "--help" || arg == "-h") {
            printUsage(argv[0]);
            std::exit(0);
        } else {
            throw std::invalid_argument("Unknown option: " + arg);
        }
    }
    return options;
}

ImageUInt8Ptr makeSyntheticImage(const Options& options) {
    if (options.synthetic == "fixture") {
        auto image = ImageUInt8::create(4, 4);
        const std::vector<std::uint8_t> values = {
            3, 3, 2, 2, 3, 4, 4, 2, 1, 4, 5, 2, 1, 1, 5, 0,
        };
        for (std::size_t i = 0; i < values.size(); ++i) {
            (*image)[static_cast<int>(i)] = values[i];
        }
        return image;
    }

    if (options.rows <= 0 || options.columns <= 0) {
        throw std::invalid_argument("Synthetic image dimensions must be positive.");
    }

    auto image = ImageUInt8::create(options.rows, options.columns);
    for (int row = 0; row < options.rows; ++row) {
        for (int column = 0; column < options.columns; ++column) {
            std::uint8_t value = 0;
            if (options.synthetic == "ramp") {
                value = static_cast<std::uint8_t>((row * 17 + column * 29) & 0xff);
            } else if (options.synthetic == "checker") {
                value = ((row + column) % 2 == 0) ? std::uint8_t{220} : std::uint8_t{35};
            } else {
                throw std::invalid_argument("Unknown synthetic image. Expected fixture, ramp, or checker.");
            }
            (*image)[row * options.columns + column] = value;
        }
    }
    return image;
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

ImageUInt8Ptr makeInputImage(const Options& options) {
    if (!options.imagePath.empty()) {
        return loadGrayscaleImage(options.imagePath);
    }
    return makeSyntheticImage(options);
}

std::ofstream openCsv(const std::filesystem::path& path) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Cannot write CSV: " + path.string());
    }
    return out;
}

void writeFamilyHeader(std::ostream& out, const std::string& prefix) {
    out << ',' << prefix << "q1" << ',' << prefix << "q2" << ',' << prefix << "qd" << ',' << prefix << "q3" << ',' << prefix << "q4";
}

template <class FamilyValue> void writeFamily(std::ostream& out, const FamilyValue& value) {
    out << ',' << value.q1 << ',' << value.q2 << ',' << value.qd << ',' << value.q3 << ',' << value.q4;
}

void writeStateHeader(std::ostream& out, const std::string& prefix, std::size_t firstState) {
    for (std::size_t state = firstState; state < 16; ++state) {
        out << ',' << prefix << "s" << state;
    }
}

void writeState(std::ostream& out, const StateHistogram& histogram) {
    for (int value : histogram.values()) {
        out << ',' << value;
    }
}

void writeState(std::ostream& out, const StateIncrement& increment) {
    for (int value : increment.bins) {
        out << ',' << value;
    }
}

void writeNodeCsv(const ValuedMorphologicalTree<ToSGrayLevel>& valuedTree, const std::filesystem::path& path,
                  const std::vector<FamilyIncrement>& familyIncrements, const std::vector<FamilyCounts>& familyCounts) {
    const MorphologicalTree& tree = valuedTree.topology();
    auto out = openCsv(path);
    out << "node,parent,altitude,direct_proper_parts";
    writeFamilyHeader(out, "increment_");
    writeFamilyHeader(out, "count_");
    out << '\n';

    for (NodeId nodeId : tree.aliveNodeIds()) {
        out << nodeId << ',' << tree.parent(nodeId) << ',' << static_cast<int>(valuedTree.nodeAltitude(nodeId)) << ',' << tree.properPartCardinality(nodeId);
        writeFamily(out, familyIncrements[static_cast<std::size_t>(nodeId)]);
        writeFamily(out, familyCounts[static_cast<std::size_t>(nodeId)]);
        out << '\n';
    }
}

void writeStateCsv(const MorphologicalTree& tree, const std::filesystem::path& path, const std::vector<StateIncrement>& stateIncrements,
                   const std::vector<StateHistogram>& stateCounts) {
    auto out = openCsv(path);
    out << "node,parent";
    writeStateHeader(out, "increment_", 1);
    writeStateHeader(out, "count_", 0);
    out << '\n';

    for (NodeId nodeId : tree.aliveNodeIds()) {
        out << nodeId << ',' << tree.parent(nodeId);
        writeState(out, stateIncrements[static_cast<std::size_t>(nodeId)]);
        writeState(out, stateCounts[static_cast<std::size_t>(nodeId)]);
        out << '\n';
    }
}

void writeProperPartCsv(const ValuedMorphologicalTree<ToSGrayLevel>& valuedTree, const std::filesystem::path& path,
                        const std::vector<FamilyIncrement>& projectedIncrements, const std::vector<FamilyCounts>& projectedCounts) {
    const MorphologicalTree& tree = valuedTree.topology();
    auto out = openCsv(path);
    out << "pixel,row,column,smallest_node_id,smallest_node_altitude";
    writeFamilyHeader(out, "smallest_node_increment_");
    writeFamilyHeader(out, "smallest_node_count_");
    out << '\n';

    for (PixelId pixel = 0; pixel < tree.numPixels(); ++pixel) {
        const auto [row, column] = ImageUtils::to2D(pixel, tree.numColumns());
        const NodeId smallestNodeId = tree.smallestNode(pixel);
        out << pixel << ',' << row << ',' << column << ',' << smallestNodeId << ',' << static_cast<int>(valuedTree.nodeAltitude(smallestNodeId));
        writeFamily(out, projectedIncrements[static_cast<std::size_t>(pixel)]);
        writeFamily(out, projectedCounts[static_cast<std::size_t>(pixel)]);
        out << '\n';
    }
}

void writeNodeSupportCsv(const MorphologicalTree& tree, const std::filesystem::path& path) {
    auto out = openCsv(path);
    out << "node,pixel,row,column,smallest_node_id\n";
    for (NodeId nodeId : tree.aliveNodeIds()) {
        for (PixelId pixel : tree.nodeSupport(nodeId)) {
            const auto [row, column] = ImageUtils::to2D(pixel, tree.numColumns());
            out << nodeId << ',' << pixel << ',' << row << ',' << column << ',' << tree.smallestNode(pixel) << '\n';
        }
    }
}

void writeMetadata(const ValuedMorphologicalTree<ToSGrayLevel>& valuedTree, const Options& options, const std::filesystem::path& path) {
    const MorphologicalTree& tree = valuedTree.topology();
    auto out = openCsv(path);
    out << "key,value\n";
    out << "rows," << tree.numRows() << '\n';
    out << "columns," << tree.numColumns() << '\n';
    out << "num_proper_parts," << tree.numPixels() << '\n';
    out << "num_internal_node_slots," << tree.numInternalNodeSlots() << '\n';
    out << "num_alive_nodes," << tree.numNodes() << '\n';
    out << "root," << tree.root() << '\n';
    out << "immersion," << toString(options.immersion) << '\n';
    out << "input," << (options.imagePath.empty() ? options.synthetic : options.imagePath.string()) << '\n';
}

} // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);
        const ImageUInt8Ptr image = makeInputImage(options);
        auto valuedTree = MorphologicalTreeFactory::createTreeOfShapes(image, makeConvention(options.immersion, image->getNumRows(), image->getNumColumns()));
        const MorphologicalTree& tree = valuedTree.topology();

        const auto stateIncrements = BitquadFiniteWindowComputation::computeNonemptyBitquadStateHistogramIncrements(tree);
        const auto nonemptyStateCounts = BitquadFiniteWindowComputation::aggregateNonemptyBitquadStateHistogramIncrements(tree, stateIncrements);
        const auto stateCounts = BitquadFiniteWindowComputation::materializeEmptyBitquadCount(tree, nonemptyStateCounts);
        const auto familyIncrements = BitquadFiniteWindowComputation::computeBitquadFamilyIncrements(tree);
        const auto familyCounts = BitquadFiniteWindowComputation::aggregateBitquadFamilyIncrements(tree, familyIncrements);
        const auto projectedIncrements = BitquadFiniteWindowComputation::projectBitquadFamilyIncrementsToProperParts(tree, familyIncrements);
        const auto projectedCounts = BitquadFiniteWindowComputation::projectBitquadFamilyCountsToProperParts(tree, familyCounts);

        std::filesystem::create_directories(options.outDir);
        writeMetadata(valuedTree, options, options.outDir / "metadata.csv");
        writeNodeCsv(valuedTree, options.outDir / "node_bitquad_families.csv", familyIncrements, familyCounts);
        writeStateCsv(tree, options.outDir / "node_bitquad_states.csv", stateIncrements, stateCounts);
        writeProperPartCsv(valuedTree, options.outDir / "proper_part_bitquad_projection.csv", projectedIncrements, projectedCounts);
        writeNodeSupportCsv(tree, options.outDir / "node_support.csv");

        std::cout << "Wrote Tree-of-Shapes bitquad projection CSV files to " << options.outDir << '\n';
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "error: " << ex.what() << '\n';
        printUsage(argv[0]);
        return 1;
    }
}
