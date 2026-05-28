/**
 * Internal validation program for Tree-of-Shapes bitquad local-event data.
 *
 * Build with `-DMMCFILTERS_BUILD_EXAMPLES=ON` and run:
 * `./build/examples/mmcfilters_tos_bitquad_naive_validation path/to/image.pgm`.
 *
 * The program compares implementation histograms, aggregated deltas, and
 * proper-part projections against a naive reference. It includes `detail/`
 * headers and is not public API.
 */
#include "mmcfilters/attributes/computers/detail/BitquadAttributeData.hpp"
#include "mmcfilters/attributes/computers/detail/BitquadLocalEventComputation.hpp"
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

namespace maf = mmcfilters;
namespace bitquad_detail = mmcfilters::attributes::computers::detail;

namespace {

using Histogram = bitquad_detail::BitquadLocalEventComputation::BitquadStateHistogram;
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
    int cols = 0;
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
    image.cols = width;
    image.pixels.assign(raw, raw + static_cast<std::size_t>(height * width));
    stbi_image_free(raw);
    return image;
}

maf::ImageUInt8Ptr makeImage(int rows, int cols, std::span<const std::uint8_t> values) {
    auto image = maf::ImageUInt8::create(rows, cols);
    for (int i = 0; i < rows * cols; ++i) {
        (*image)[i] = values[static_cast<std::size_t>(i)];
    }
    return image;
}

const char* interpolationName(maf::ToSInterpolation interpolation) {
    switch (interpolation) {
        case maf::ToSInterpolation::SelfDual:
            return "SelfDual";
        case maf::ToSInterpolation::Min4cMax8c:
            return "Min4cMax8c";
        case maf::ToSInterpolation::Min8cMax4c:
            return "Min8cMax4c";
    }
    return "unknown";
}

std::vector<Histogram> naiveSupportBitquadHistograms(const maf::MorphologicalTree& tree) {
    const int rows = tree.getNumRowsOfImage();
    const int cols = tree.getNumColsOfImage();
    const std::array<std::pair<int, int>, 4> offsets = {{
        {0, 0},
        {1, 0},
        {0, 1},
        {1, 1},
    }};

    std::vector<Histogram> histograms(static_cast<std::size_t>(tree.getNumInternalNodeSlots()));

    std::vector<std::uint8_t> stateByNode(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), 0);
    std::vector<maf::NodeId> touchedNodes;
    touchedNodes.reserve(128);

    for (int row = -1; row < rows; ++row) {
        for (int col = -1; col < cols; ++col) {
            touchedNodes.clear();
            for (std::size_t bit = 0; bit < offsets.size(); ++bit) {
                const int qRow = row + offsets[bit].first;
                const int qCol = col + offsets[bit].second;
                if (qRow < 0 || qRow >= rows || qCol < 0 || qCol >= cols) {
                    continue;
                }

                const int q = maf::ImageUtils::to1D(qRow, qCol, cols);
                for (maf::NodeId nodeId = tree.getProperPartOwner(q);
                     nodeId != maf::InvalidNode;
                     nodeId = tree.isRoot(nodeId) ? maf::InvalidNode : tree.getNodeParent(nodeId)) {
                    auto& state = stateByNode[static_cast<std::size_t>(nodeId)];
                    if (state == 0) {
                        touchedNodes.push_back(nodeId);
                    }
                    state |= static_cast<std::uint8_t>(std::uint8_t{1} << bit);
                }
            }

            for (maf::NodeId nodeId : touchedNodes) {
                auto& state = stateByNode[static_cast<std::size_t>(nodeId)];
                histograms[static_cast<std::size_t>(nodeId)][static_cast<std::size_t>(state)] += 1;
                state = 0;
            }
        }
    }

    const int totalCells = (rows + 1) * (cols + 1);
    for (maf::NodeId nodeId : tree.getAliveNodeIds()) {
        Histogram& histogram = histograms[static_cast<std::size_t>(nodeId)];
        int nonEmpty = 0;
        for (std::size_t state = 1; state < histogram.size(); ++state) {
            nonEmpty += histogram[state];
        }
        histogram[0] = totalCells - nonEmpty;
    }

    return histograms;
}

Summary compareHistograms(
    const maf::MorphologicalTree& tree,
    std::span<const Histogram> actual,
    std::span<const Histogram> expected,
    const char* label) {
    Summary summary;
    int printed = 0;
    for (maf::NodeId nodeId : tree.getAliveNodeIds()) {
        bool nodeDiffers = false;
        for (std::size_t state = 0; state < 16; ++state) {
            const int diff =
                actual[static_cast<std::size_t>(nodeId)][state] -
                expected[static_cast<std::size_t>(nodeId)][state];
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
                    const int a = actual[static_cast<std::size_t>(nodeId)][state];
                    const int e = expected[static_cast<std::size_t>(nodeId)][state];
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

std::array<int, 6> familyValues(const Families& families) {
    return {{
        families.empty,
        families.q1,
        families.q2,
        families.qd,
        families.q3,
        families.q4,
    }};
}

FamilySummary compareFamilies(
    const maf::MorphologicalTree& tree,
    std::span<const Families> actual,
    std::span<const Families> expected,
    const char* label) {
    FamilySummary summary;
    int printed = 0;
    for (maf::NodeId nodeId : tree.getAliveNodeIds()) {
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
                std::cout << "  " << label << " mismatch node=" << nodeId
                          << " actual={empty=" << a[0] << " q1=" << a[1]
                          << " q2=" << a[2] << " qd=" << a[3]
                          << " q3=" << a[4] << " q4=" << a[5] << "}"
                          << " expected={empty=" << e[0] << " q1=" << e[1]
                          << " q2=" << e[2] << " qd=" << e[3]
                          << " q3=" << e[4] << " q4=" << e[5] << "}\n";
                ++printed;
            }
        }
    }
    return summary;
}

FamilySummary compareProjectedFamilies(
    const maf::MorphologicalTree& tree,
    std::span<const Families> projected,
    std::span<const Families> nodeExpected,
    const char* label) {
    FamilySummary summary;
    int printed = 0;
    for (maf::NodeId properPart = 0; properPart < tree.getNumTotalProperParts(); ++properPart) {
        const maf::NodeId owner = tree.getProperPartOwner(properPart);
        const auto a = familyValues(projected[static_cast<std::size_t>(properPart)]);
        const auto e = familyValues(nodeExpected[static_cast<std::size_t>(owner)]);
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
                std::cout << "  " << label << " mismatch proper_part=" << properPart
                          << " owner=" << owner << "\n";
                ++printed;
            }
        }
    }
    return summary;
}

void printSummary(const char* name, const Summary& summary) {
    std::cout << name
              << " mismatched_nodes=" << summary.mismatchedNodes
              << " mismatched_entries=" << summary.mismatchedEntries
              << " abs_diff=" << summary.absDiff << "\n";
}

void printSummary(const char* name, const FamilySummary& summary) {
    std::cout << name
              << " mismatched_items=" << summary.mismatchedNodes
              << " mismatched_fields=" << summary.mismatchedFields
              << " abs_diff=" << summary.absDiff << "\n";
}

bool compareOne(
    const std::filesystem::path& imagePath,
    int rows,
    int cols,
    std::span<const std::uint8_t> pixels,
    maf::ToSInterpolation interpolation) {
    auto weighted = maf::MorphologicalTreeFactory::createTreeOfShapes(
        makeImage(rows, cols, pixels),
        interpolation);
    const maf::MorphologicalTree& tree = weighted.topology();

    const auto naiveHistograms = naiveSupportBitquadHistograms(tree);
    const auto naiveFamilies = bitquad_detail::BitquadLocalEventComputation::computeBitquadFamilyCounts(naiveHistograms);

    const auto mafHistograms = bitquad_detail::BitquadLocalEventComputation::computeBitquadStateHistograms(tree);
    const auto mafFamilies = bitquad_detail::BitquadLocalEventComputation::computeBitquadFamilyCounts(tree);
    const auto mafFamilyDeltas = bitquad_detail::BitquadLocalEventComputation::computeBitquadFamilyDeltas(tree);
    const auto mafAggregatedFamilies = bitquad_detail::BitquadLocalEventComputation::aggregateBitquadFamilyDeltas(tree, mafFamilyDeltas);
    const auto mafStateDeltas = bitquad_detail::BitquadLocalEventComputation::computeBitquadStateHistogramDeltas(tree);
    const auto mafAggregatedHistograms = bitquad_detail::BitquadLocalEventComputation::aggregateBitquadStateHistogramDeltas(tree, mafStateDeltas);
    const auto mafProjectedFamilies = bitquad_detail::BitquadLocalEventComputation::projectBitquadFamilyCountsToProperParts(tree, mafFamilies);

    std::cout << "image=\"" << imagePath.string() << "\" size=" << rows << "x" << cols
              << " interpolation=" << interpolationName(interpolation)
              << " nodes=" << tree.getNumNodes()
              << " node_slots=" << tree.getNumInternalNodeSlots()
              << " pixels=" << tree.getNumTotalProperParts() << "\n";

    const Summary stateSummary = compareHistograms(tree, mafHistograms, naiveHistograms, "state");
    const Summary aggregatedStateSummary = compareHistograms(tree, mafAggregatedHistograms, naiveHistograms, "aggregated-state-delta");
    const FamilySummary familySummary = compareFamilies(tree, mafFamilies, naiveFamilies, "family");
    const FamilySummary aggregatedFamilySummary = compareFamilies(tree, mafAggregatedFamilies, naiveFamilies, "aggregated-family-delta");
    const FamilySummary projectedSummary = compareProjectedFamilies(tree, mafProjectedFamilies, naiveFamilies, "proper-part-projection");

    printSummary("state", stateSummary);
    printSummary("aggregated-state-delta", aggregatedStateSummary);
    printSummary("family", familySummary);
    printSummary("aggregated-family-delta", aggregatedFamilySummary);
    printSummary("proper-part-projection", projectedSummary);

    return stateSummary.mismatchedNodes == 0 &&
           aggregatedStateSummary.mismatchedNodes == 0 &&
           familySummary.mismatchedNodes == 0 &&
           aggregatedFamilySummary.mismatchedNodes == 0 &&
           projectedSummary.mismatchedNodes == 0;
}

} // namespace

int main(int argc, char** argv) {
    const std::filesystem::path imagePath = argc > 1
        ? std::filesystem::path(argv[1])
        : std::filesystem::path("dat/lena.pgm");
    if (argc > 2) {
        std::cerr << "usage: " << argv[0] << " [image]\n";
        return 64;
    }

    const LoadedImage image = loadImage(imagePath);
    bool ok = true;
    ok = compareOne(imagePath, image.rows, image.cols, image.pixels, maf::ToSInterpolation::SelfDual) && ok;
    ok = compareOne(imagePath, image.rows, image.cols, image.pixels, maf::ToSInterpolation::Min4cMax8c) && ok;
    ok = compareOne(imagePath, image.rows, image.cols, image.pixels, maf::ToSInterpolation::Min8cMax4c) && ok;
    return ok ? 0 : 2;
}
