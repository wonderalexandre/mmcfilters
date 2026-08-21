/**
 * Internal validation program for Tree-of-Shapes bitquad finite-window data.
 *
 * Build with `-DMMCFILTERS_BUILD_EXAMPLES=ON` and run:
 * `./build/examples/mmcfilters_tos_bitquad_naive_validation path/to/image.pgm`.
 *
 * The program compares implementation histograms, aggregated increments, and
 * proper-part projections against a naive reference. It includes `detail/`
 * headers and is not public API.
 */
#include "mmcfilters/attributes/computers/detail/BitquadAttributeData.hpp"
#include "mmcfilters/attributes/computers/detail/BitquadFiniteWindowComputation.hpp"
#include "mmcfilters/trees/MorphologicalTreeFactory.hpp"
#include "mmcfilters/utils/Image.hpp"
#include "stb_image.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace bitquad_detail = mmcfilters::attributes::computers::detail;

namespace {

using Histogram = bitquad_detail::BitquadFiniteWindowComputation::BitquadStateHistogram;
using Families = bitquad_detail::BitquadFamilyCounts;

struct Summary {
    int mismatchedNodes = 0;
    int mismatchedEntries = 0;
    long long absDiff = 0;
};

struct FamilySummary {
    int mismatchedNodes = 0;
    int mismatchedFields = 0;
    long long absDiff = 0;
};

struct LoadedImage {
    int rows = 0;
    int columns = 0;
    std::vector<std::uint8_t> pixels;
};

LoadedImage loadImage(const std::filesystem::path& path) {
    int width = 0;
    int height = 0;
    int channels = 0;
    unsigned char* raw = stbi_load(path.string().c_str(), &width, &height, &channels, 1);
    if (raw == nullptr) {
        throw std::runtime_error("Failed to load image: " + path.string());
    }
    if (width <= 0 || height <= 0) {
        stbi_image_free(raw);
        throw std::invalid_argument("Loaded image must have positive dimensions.");
    }

    LoadedImage image;
    image.rows = height;
    image.columns = width;
    image.pixels.assign(raw, raw + static_cast<std::size_t>(height * width));
    stbi_image_free(raw);
    return image;
}

mmcfilters::ImageUInt8Ptr makeImage(int rows, int columns, std::span<const std::uint8_t> values) {
    auto image = mmcfilters::ImageUInt8::create(rows, columns);
    for (int i = 0; i < rows * columns; ++i) {
        (*image)[i] = values[static_cast<std::size_t>(i)];
    }
    return image;
}

enum class ExampleImmersion { SelfDualSpan, Min4Max8, Min8Max4 };

mmcfilters::TopographicConvention makeConvention(ExampleImmersion immersion, int rows, int columns) {
    if (immersion == ExampleImmersion::SelfDualSpan) {
        return mmcfilters::selfDualSpanConvention();
    }
    const bool minIs4 = immersion == ExampleImmersion::Min4Max8;
    return mmcfilters::TopographicConvention{
        mmcfilters::ComplementaryGridImmersion{mmcfilters::ComplementaryAdjacencies{
            mmcfilters::RegularGridAdjacency2D(rows, columns, minIs4 ? 1.0 : 1.5),
            mmcfilters::RegularGridAdjacency2D(rows, columns, minIs4 ? 1.5 : 1.0)}},
        mmcfilters::TopographicDomainExtension::ExteriorRing, mmcfilters::PixelId{0},
        mmcfilters::TopographicAltitudeEncoding::ExactDoubled};
}

const char* immersionName(ExampleImmersion immersion) {
    switch (immersion) {
    case ExampleImmersion::SelfDualSpan:
        return "SelfDual";
    case ExampleImmersion::Min4Max8:
        return "Min4cMax8c";
    case ExampleImmersion::Min8Max4:
        return "Min8cMax4c";
    }
    return "unknown";
}

std::vector<Histogram> naiveSupportBitquadHistograms(const mmcfilters::MorphologicalTree& tree) {
    const int rows = tree.numRows();
    const int columns = tree.numColumns();
    const std::array<std::pair<int, int>, 4> offsets = {{
        {0, 0},
        {0, 1},
        {1, 0},
        {1, 1},
    }};

    std::vector<Histogram> histograms(static_cast<std::size_t>(tree.numInternalNodeSlots()));

    std::vector<std::uint8_t> stateByNode(static_cast<std::size_t>(tree.numInternalNodeSlots()), 0);
    std::vector<mmcfilters::NodeId> touchedNodes;
    touchedNodes.reserve(128);

    for (int row = -1; row < rows; ++row) {
        for (int column = -1; column < columns; ++column) {
            touchedNodes.clear();
            for (std::size_t bit = 0; bit < offsets.size(); ++bit) {
                const int qRow = row + offsets[bit].first;
                const int qColumn = column + offsets[bit].second;
                if (qRow < 0 || qRow >= rows || qColumn < 0 || qColumn >= columns) {
                    continue;
                }

                const int q = mmcfilters::ImageUtils::to1D(qRow, qColumn, columns);
                for (mmcfilters::NodeId nodeId = tree.smallestNode(q); nodeId != mmcfilters::InvalidNode;
                     nodeId = tree.isRoot(nodeId) ? mmcfilters::InvalidNode : tree.parent(nodeId)) {
                    auto& state = stateByNode[static_cast<std::size_t>(nodeId)];
                    if (state == 0) {
                        touchedNodes.push_back(nodeId);
                    }
                    state |= static_cast<std::uint8_t>(std::uint8_t{1} << bit);
                }
            }

            for (mmcfilters::NodeId nodeId : touchedNodes) {
                auto& state = stateByNode[static_cast<std::size_t>(nodeId)];
                histograms[static_cast<std::size_t>(nodeId)].count(state) += 1;
                state = 0;
            }
        }
    }

    const int totalCells = (rows + 1) * (columns + 1);
    for (mmcfilters::NodeId nodeId : tree.aliveNodeIds()) {
        Histogram& histogram = histograms[static_cast<std::size_t>(nodeId)];
        int nonEmpty = 0;
        for (std::uint8_t state = 1; state < 16; ++state) {
            nonEmpty += histogram.count(state);
        }
        histogram.count(0) = totalCells - nonEmpty;
    }

    return histograms;
}

Summary compareHistograms(const mmcfilters::MorphologicalTree& tree, std::span<const Histogram> actual,
                          std::span<const Histogram> expected, const char* label) {
    Summary summary;
    int printed = 0;
    for (mmcfilters::NodeId nodeId : tree.aliveNodeIds()) {
        bool nodeDiffers = false;
        for (std::size_t state = 0; state < 16; ++state) {
            const int diff = actual[static_cast<std::size_t>(nodeId)].count(static_cast<std::uint8_t>(state)) -
                             expected[static_cast<std::size_t>(nodeId)].count(static_cast<std::uint8_t>(state));
            if (diff == 0) {
                continue;
            }
            nodeDiffers = true;
            ++summary.mismatchedEntries;
            summary.absDiff += std::abs(diff);
        }
        if (nodeDiffers) {
            ++summary.mismatchedNodes;
            if (printed < 5) {
                std::cout << "  " << label << " mismatch node=" << nodeId << " states:";
                for (std::size_t state = 0; state < 16; ++state) {
                    const int a = actual[static_cast<std::size_t>(nodeId)].count(static_cast<std::uint8_t>(state));
                    const int e = expected[static_cast<std::size_t>(nodeId)].count(static_cast<std::uint8_t>(state));
                    if (a != e) {
                        std::cout << " s" << state << "=" << a << "/" << e;
                    }
                }
                std::cout << "\n";
                ++printed;
            }
        }
    }
    return summary;
}

std::array<int, 5> familyValues(const Families& families) {
    return {{
        families.q1,
        families.q2,
        families.qd,
        families.q3,
        families.q4,
    }};
}

FamilySummary compareFamilies(const mmcfilters::MorphologicalTree& tree, std::span<const Families> actual,
                              std::span<const Families> expected, const char* label) {
    FamilySummary summary;
    int printed = 0;
    for (mmcfilters::NodeId nodeId : tree.aliveNodeIds()) {
        const auto a = familyValues(actual[static_cast<std::size_t>(nodeId)]);
        const auto e = familyValues(expected[static_cast<std::size_t>(nodeId)]);
        bool nodeDiffers = false;
        for (std::size_t i = 0; i < a.size(); ++i) {
            const int diff = a[i] - e[i];
            if (diff == 0) {
                continue;
            }
            nodeDiffers = true;
            ++summary.mismatchedFields;
            summary.absDiff += std::abs(diff);
        }
        if (nodeDiffers) {
            ++summary.mismatchedNodes;
            if (printed < 5) {
                std::cout << "  " << label << " mismatch node=" << nodeId << " actual={q1=" << a[0] << " q2=" << a[1] << " qd=" << a[2] << " q3=" << a[3]
                          << " q4=" << a[4] << "}"
                          << " expected={q1=" << e[0] << " q2=" << e[1] << " qd=" << e[2] << " q3=" << e[3] << " q4=" << e[4] << "}\n";
                ++printed;
            }
        }
    }
    return summary;
}

FamilySummary compareProjectedFamilies(const mmcfilters::MorphologicalTree& tree, std::span<const Families> projected,
                                       std::span<const Families> nodeExpected, const char* label) {
    FamilySummary summary;
    int printed = 0;
    for (mmcfilters::PixelId pixel = 0; pixel < tree.numPixels(); ++pixel) {
        const mmcfilters::NodeId smallestNodeId = tree.smallestNode(pixel);
        const auto a = familyValues(projected[static_cast<std::size_t>(pixel)]);
        const auto e = familyValues(nodeExpected[static_cast<std::size_t>(smallestNodeId)]);
        bool differs = false;
        for (std::size_t i = 0; i < a.size(); ++i) {
            const int diff = a[i] - e[i];
            if (diff == 0) {
                continue;
            }
            differs = true;
            ++summary.mismatchedFields;
            summary.absDiff += std::abs(diff);
        }
        if (differs) {
            ++summary.mismatchedNodes;
            if (printed < 5) {
                std::cout << "  " << label << " mismatch proper_part=" << pixel << " smallestNodeId=" << smallestNodeId << "\n";
                ++printed;
            }
        }
    }
    return summary;
}

void printSummary(const char* name, const Summary& summary) {
    std::cout << name << " mismatched_nodes=" << summary.mismatchedNodes << " mismatched_entries=" << summary.mismatchedEntries
              << " abs_diff=" << summary.absDiff << "\n";
}

void printSummary(const char* name, const FamilySummary& summary) {
    std::cout << name << " mismatched_items=" << summary.mismatchedNodes << " mismatched_fields=" << summary.mismatchedFields << " abs_diff=" << summary.absDiff
              << "\n";
}

bool compareOne(const std::filesystem::path& imagePath, int rows, int columns, std::span<const std::uint8_t> pixels, ExampleImmersion immersion) {
    auto valuedTree = mmcfilters::MorphologicalTreeFactory::createTreeOfShapes<mmcfilters::ToSGrayLevel>(
        makeImage(rows, columns, pixels), makeConvention(immersion, rows, columns));
    const mmcfilters::MorphologicalTree& tree = valuedTree.topology();

    const auto naiveHistograms = naiveSupportBitquadHistograms(tree);
    const auto naiveFamilies = bitquad_detail::BitquadFiniteWindowComputation::computeBitquadFamilyCounts(naiveHistograms);

    const auto mafHistograms = bitquad_detail::BitquadFiniteWindowComputation::computeBitquadStateHistograms(tree);
    const auto mafFamilies = bitquad_detail::BitquadFiniteWindowComputation::computeBitquadFamilyCounts(tree);
    const auto mafFamilyIncrements = bitquad_detail::BitquadFiniteWindowComputation::computeBitquadFamilyIncrements(tree);
    const auto mafAggregatedFamilies = bitquad_detail::BitquadFiniteWindowComputation::aggregateBitquadFamilyIncrements(tree, mafFamilyIncrements);
    const auto mafStateIncrements = bitquad_detail::BitquadFiniteWindowComputation::computeNonemptyBitquadStateHistogramIncrements(tree);
    const auto mafAggregatedNonemptyHistograms =
        bitquad_detail::BitquadFiniteWindowComputation::aggregateNonemptyBitquadStateHistogramIncrements(tree, mafStateIncrements);
    const auto mafAggregatedHistograms = bitquad_detail::BitquadFiniteWindowComputation::materializeEmptyBitquadCount(tree, mafAggregatedNonemptyHistograms);
    const auto mafProjectedFamilies = bitquad_detail::BitquadFiniteWindowComputation::projectBitquadFamilyCountsToProperParts(tree, mafFamilies);

    std::cout << "image=\"" << imagePath.string() << "\" size=" << rows << "x" << columns << " immersion=" << immersionName(immersion)
              << " nodes=" << tree.numNodes() << " node_slots=" << tree.numInternalNodeSlots() << " pixels=" << tree.numPixels() << "\n";

    const Summary stateSummary = compareHistograms(tree, mafHistograms, naiveHistograms, "state");
    const Summary aggregatedStateSummary = compareHistograms(tree, mafAggregatedHistograms, naiveHistograms, "aggregated-state-increment");
    const FamilySummary familySummary = compareFamilies(tree, mafFamilies, naiveFamilies, "family");
    const FamilySummary aggregatedFamilySummary = compareFamilies(tree, mafAggregatedFamilies, naiveFamilies, "aggregated-family-increment");
    const FamilySummary projectedSummary = compareProjectedFamilies(tree, mafProjectedFamilies, naiveFamilies, "proper-part-projection");

    printSummary("state", stateSummary);
    printSummary("aggregated-state-increment", aggregatedStateSummary);
    printSummary("family", familySummary);
    printSummary("aggregated-family-increment", aggregatedFamilySummary);
    printSummary("proper-part-projection", projectedSummary);

    return stateSummary.mismatchedNodes == 0 && aggregatedStateSummary.mismatchedNodes == 0 && familySummary.mismatchedNodes == 0 &&
           aggregatedFamilySummary.mismatchedNodes == 0 && projectedSummary.mismatchedNodes == 0;
}

} // namespace

int main(int argc, char** argv) {
    const std::filesystem::path imagePath = argc > 1 ? std::filesystem::path(argv[1]) : std::filesystem::path("dat/lena.pgm");
    if (argc > 2) {
        std::cerr << "usage: " << argv[0] << " [image]\n";
        return 64;
    }

    const LoadedImage image = loadImage(imagePath);
    bool ok = true;
    ok = compareOne(imagePath, image.rows, image.columns, image.pixels, ExampleImmersion::SelfDualSpan) && ok;
    ok = compareOne(imagePath, image.rows, image.columns, image.pixels, ExampleImmersion::Min4Max8) && ok;
    ok = compareOne(imagePath, image.rows, image.columns, image.pixels, ExampleImmersion::Min8Max4) && ok;
    return ok ? 0 : 2;
}
