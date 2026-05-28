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
#include "mmcfilters/attributes/computers/detail/BitquadLocalEventComputation.hpp"
#include "mmcfilters/trees/MorphologicalTreeFactory.hpp"
#include "mmcfilters/trees/WeightedMorphologicalTree.hpp"
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

using FamilyCounts = mmcfilters::attributes::computers::detail::BitquadFamilyCounts;
using StateHistogram = BitquadLocalEventComputation::BitquadStateHistogram;

struct Options {
    std::filesystem::path imagePath;
    std::filesystem::path outDir = "tos-bitquad-projections";
    std::string synthetic = "fixture";
    ToSInterpolation interpolation = ToSInterpolation::SelfDual;
    int rows = 8;
    int cols = 8;
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
    throw std::invalid_argument("Unsupported Tree-of-Shapes interpolation.");
}

ToSInterpolation parseInterpolation(const std::string& value) {
    if (value == "SelfDual") {
        return ToSInterpolation::SelfDual;
    }
    if (value == "Min4cMax8c") {
        return ToSInterpolation::Min4cMax8c;
    }
    if (value == "Min8cMax4c") {
        return ToSInterpolation::Min8cMax4c;
    }
    throw std::invalid_argument("Expected SelfDual, Min4cMax8c, or Min8cMax4c for --interpolation.");
}

void printUsage(const char* program) {
    std::cerr
        << "Usage: " << program << " [--image path] [--out-dir dir]\n"
        << "       [--interpolation SelfDual|Min4cMax8c|Min8cMax4c]\n"
        << "       [--synthetic fixture|ramp|checker] [--rows n] [--cols n]\n";
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
        } else if (arg == "--interpolation") {
            options.interpolation = parseInterpolation(requireValue("--interpolation"));
        } else if (arg == "--synthetic") {
            options.synthetic = requireValue("--synthetic");
        } else if (arg == "--rows") {
            options.rows = std::stoi(requireValue("--rows"));
        } else if (arg == "--cols") {
            options.cols = std::stoi(requireValue("--cols"));
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
            3, 3, 2, 2,
            3, 4, 4, 2,
            1, 4, 5, 2,
            1, 1, 5, 0,
        };
        for (std::size_t i = 0; i < values.size(); ++i) {
            (*image)[static_cast<int>(i)] = values[i];
        }
        return image;
    }

    if (options.rows <= 0 || options.cols <= 0) {
        throw std::invalid_argument("Synthetic image dimensions must be positive.");
    }

    auto image = ImageUInt8::create(options.rows, options.cols);
    for (int row = 0; row < options.rows; ++row) {
        for (int col = 0; col < options.cols; ++col) {
            std::uint8_t value = 0;
            if (options.synthetic == "ramp") {
                value = static_cast<std::uint8_t>((row * 17 + col * 29) & 0xff);
            } else if (options.synthetic == "checker") {
                value = ((row + col) % 2 == 0) ? std::uint8_t{220} : std::uint8_t{35};
            } else {
                throw std::invalid_argument("Unknown synthetic image. Expected fixture, ramp, or checker.");
            }
            (*image)[row * options.cols + col] = value;
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
    out << ',' << prefix << "empty"
        << ',' << prefix << "q1"
        << ',' << prefix << "q2"
        << ',' << prefix << "qd"
        << ',' << prefix << "q3"
        << ',' << prefix << "q4";
}

void writeFamily(std::ostream& out, const FamilyCounts& counts) {
    out << ',' << counts.empty
        << ',' << counts.q1
        << ',' << counts.q2
        << ',' << counts.qd
        << ',' << counts.q3
        << ',' << counts.q4;
}

void writeStateHeader(std::ostream& out, const std::string& prefix) {
    for (std::size_t state = 0; state < 16; ++state) {
        out << ',' << prefix << "s" << state;
    }
}

void writeState(std::ostream& out, const StateHistogram& histogram) {
    for (int value : histogram) {
        out << ',' << value;
    }
}

void writeNodeCsv(
    const WeightedMorphologicalTree<std::uint8_t>& weighted,
    const std::filesystem::path& path,
    const std::vector<FamilyCounts>& familyDeltas,
    const std::vector<FamilyCounts>& familyCounts) {
    const MorphologicalTree& tree = weighted.topology();
    auto out = openCsv(path);
    out << "node,parent,altitude,direct_proper_parts";
    writeFamilyHeader(out, "delta_");
    writeFamilyHeader(out, "count_");
    out << '\n';

    for (NodeId nodeId : tree.getAliveNodeIds()) {
        out << nodeId
            << ',' << tree.getNodeParent(nodeId)
            << ',' << static_cast<int>(weighted.getAltitude(nodeId))
            << ',' << tree.getNumProperParts(nodeId);
        writeFamily(out, familyDeltas[static_cast<std::size_t>(nodeId)]);
        writeFamily(out, familyCounts[static_cast<std::size_t>(nodeId)]);
        out << '\n';
    }
}

void writeStateCsv(
    const MorphologicalTree& tree,
    const std::filesystem::path& path,
    const std::vector<StateHistogram>& stateDeltas,
    const std::vector<StateHistogram>& stateCounts) {
    auto out = openCsv(path);
    out << "node,parent";
    writeStateHeader(out, "delta_");
    writeStateHeader(out, "count_");
    out << '\n';

    for (NodeId nodeId : tree.getAliveNodeIds()) {
        out << nodeId << ',' << tree.getNodeParent(nodeId);
        writeState(out, stateDeltas[static_cast<std::size_t>(nodeId)]);
        writeState(out, stateCounts[static_cast<std::size_t>(nodeId)]);
        out << '\n';
    }
}

void writeProperPartCsv(
    const WeightedMorphologicalTree<std::uint8_t>& weighted,
    const std::filesystem::path& path,
    const std::vector<FamilyCounts>& projectedDeltas,
    const std::vector<FamilyCounts>& projectedCounts) {
    const MorphologicalTree& tree = weighted.topology();
    auto out = openCsv(path);
    out << "proper_part,row,col,owner,owner_altitude";
    writeFamilyHeader(out, "owner_delta_");
    writeFamilyHeader(out, "owner_count_");
    out << '\n';

    for (NodeId properPart = 0; properPart < tree.getNumTotalProperParts(); ++properPart) {
        const auto [row, col] = ImageUtils::to2D(properPart, tree.getNumColsOfImage());
        const NodeId owner = tree.getProperPartOwner(properPart);
        out << properPart
            << ',' << row
            << ',' << col
            << ',' << owner
            << ',' << static_cast<int>(weighted.getAltitude(owner));
        writeFamily(out, projectedDeltas[static_cast<std::size_t>(properPart)]);
        writeFamily(out, projectedCounts[static_cast<std::size_t>(properPart)]);
        out << '\n';
    }
}

void writeNodeSupportCsv(const MorphologicalTree& tree, const std::filesystem::path& path) {
    auto out = openCsv(path);
    out << "node,proper_part,row,col,direct_owner\n";
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        for (NodeId properPart : tree.getConnectedComponent(nodeId)) {
            const auto [row, col] = ImageUtils::to2D(properPart, tree.getNumColsOfImage());
            out << nodeId
                << ',' << properPart
                << ',' << row
                << ',' << col
                << ',' << tree.getProperPartOwner(properPart)
                << '\n';
        }
    }
}

void writeMetadata(
    const WeightedMorphologicalTree<std::uint8_t>& weighted,
    const Options& options,
    const std::filesystem::path& path) {
    const MorphologicalTree& tree = weighted.topology();
    auto out = openCsv(path);
    out << "key,value\n";
    out << "rows," << tree.getNumRowsOfImage() << '\n';
    out << "cols," << tree.getNumColsOfImage() << '\n';
    out << "num_proper_parts," << tree.getNumTotalProperParts() << '\n';
    out << "num_internal_node_slots," << tree.getNumInternalNodeSlots() << '\n';
    out << "num_alive_nodes," << tree.getNumNodes() << '\n';
    out << "root," << tree.getRoot() << '\n';
    out << "interpolation," << toString(options.interpolation) << '\n';
    out << "input," << (options.imagePath.empty() ? options.synthetic : options.imagePath.string()) << '\n';
}

} // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);
        const ImageUInt8Ptr image = makeInputImage(options);
        auto weighted = MorphologicalTreeFactory::createTreeOfShapes(image, options.interpolation);
        const MorphologicalTree& tree = weighted.topology();

        const auto stateDeltas = BitquadLocalEventComputation::computeBitquadStateHistogramDeltas(tree);
        const auto stateCounts = BitquadLocalEventComputation::aggregateBitquadStateHistogramDeltas(tree, stateDeltas);
        const auto familyDeltas = BitquadLocalEventComputation::computeBitquadFamilyDeltas(tree);
        const auto familyCounts = BitquadLocalEventComputation::aggregateBitquadFamilyDeltas(tree, familyDeltas);
        const auto projectedDeltas = BitquadLocalEventComputation::projectBitquadFamilyCountsToProperParts(tree, familyDeltas);
        const auto projectedCounts = BitquadLocalEventComputation::projectBitquadFamilyCountsToProperParts(tree, familyCounts);

        std::filesystem::create_directories(options.outDir);
        writeMetadata(weighted, options, options.outDir / "metadata.csv");
        writeNodeCsv(weighted, options.outDir / "node_bitquad_families.csv", familyDeltas, familyCounts);
        writeStateCsv(tree, options.outDir / "node_bitquad_states.csv", stateDeltas, stateCounts);
        writeProperPartCsv(weighted, options.outDir / "proper_part_bitquad_projection.csv", projectedDeltas, projectedCounts);
        writeNodeSupportCsv(tree, options.outDir / "node_support.csv");

        std::cout << "Wrote Tree-of-Shapes bitquad projection CSV files to "
                  << options.outDir << '\n';
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "error: " << ex.what() << '\n';
        printUsage(argv[0]);
        return 1;
    }
}
