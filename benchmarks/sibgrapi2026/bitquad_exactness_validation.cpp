#include "mmcfilters/attributes/computers/detail/BitquadFiniteWindowComputation.hpp"
#include "mmcfilters/trees/MorphologicalTreeFactory.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string_view>
#include <vector>

namespace {

namespace detail = mmcfilters::attributes::computers::detail;

using detail::BitquadCode;
using detail::BitquadFamilyCounts;
using detail::BitquadFiniteWindowComputation;
using mmcfilters::ImageUInt8;
using mmcfilters::ImageUInt8Ptr;
using mmcfilters::MorphologicalTree;
using mmcfilters::MorphologicalTreeFactory;
using mmcfilters::NodeId;
using mmcfilters::PixelId;

struct ValidationSummary {
    std::uint64_t images = 0;
    std::uint64_t trees = 0;
    std::uint64_t nodes = 0;
    std::uint64_t mismatchedNodes = 0;
    std::uint64_t maximumAbsoluteFamilyError = 0;
};

std::array<int, 5> familyValues(const BitquadFamilyCounts& counts) {
    return {counts.q1, counts.q2, counts.qd, counts.q3, counts.q4};
}

std::vector<std::uint8_t> supportMask(const MorphologicalTree& tree, NodeId node) {
    std::vector<std::uint8_t> mask(static_cast<std::size_t>(tree.numPixels()), 0);
    for (NodeId subtreeNode : tree.subtreeNodes(node)) {
        for (PixelId pixel : tree.properPart(subtreeNode)) {
            mask[static_cast<std::size_t>(pixel)] = 1;
        }
    }
    return mask;
}

BitquadFamilyCounts exhaustiveFamilyCounts(const MorphologicalTree& tree, NodeId node) {
    const int rows = tree.numRows();
    const int columns = tree.numColumns();
    const std::vector<std::uint8_t> mask = supportMask(tree, node);
    detail::BitquadStateHistogram histogram;

    for (int row = -1; row < rows; ++row) {
        for (int column = -1; column < columns; ++column) {
            BitquadCode code = 0;
            for (int bit = 0; bit < 4; ++bit) {
                const int sampleRow = row + bit / 2;
                const int sampleColumn = column + bit % 2;
                if (sampleRow < 0 || sampleRow >= rows || sampleColumn < 0 || sampleColumn >= columns) {
                    continue;
                }
                const PixelId pixel = sampleRow * columns + sampleColumn;
                if (mask[static_cast<std::size_t>(pixel)] != 0) {
                    code = static_cast<BitquadCode>(code | (BitquadCode{1} << bit));
                }
            }
            ++histogram.count(code);
        }
    }
    return BitquadFiniteWindowComputation::projectBitquadFamilyCounts(histogram);
}

void validateTree(const MorphologicalTree& tree, std::string_view hierarchy, std::uint32_t pattern, ValidationSummary& summary) {
    const std::vector<BitquadFamilyCounts> actual = BitquadFiniteWindowComputation::computeBitquadFamilyCounts(tree);
    ++summary.trees;

    for (NodeId node : tree.aliveNodeIds()) {
        ++summary.nodes;
        const BitquadFamilyCounts expected = exhaustiveFamilyCounts(tree, node);
        const auto actualValues = familyValues(actual[static_cast<std::size_t>(node)]);
        const auto expectedValues = familyValues(expected);
        bool differs = false;
        for (std::size_t family = 0; family < actualValues.size(); ++family) {
            const std::uint64_t difference = static_cast<std::uint64_t>(
                std::abs(static_cast<long long>(actualValues[family]) - static_cast<long long>(expectedValues[family])));
            differs = differs || difference != 0;
            summary.maximumAbsoluteFamilyError = std::max(summary.maximumAbsoluteFamilyError, difference);
        }
        if (differs) {
            ++summary.mismatchedNodes;
            std::cerr << "mismatch hierarchy=" << hierarchy << " pattern=" << pattern << " node=" << node << '\n';
        }
    }
}

ImageUInt8Ptr makeBinaryImage(std::uint32_t pattern) {
    auto image = ImageUInt8::create(3, 3);
    for (PixelId pixel = 0; pixel < 9; ++pixel) {
        (*image)[pixel] = (pattern & (std::uint32_t{1} << pixel)) != 0 ? std::uint8_t{255} : std::uint8_t{0};
    }
    return image;
}

} // namespace

int main() {
    ValidationSummary summary;
    for (std::uint32_t pattern = 0; pattern < (std::uint32_t{1} << 9); ++pattern) {
        const ImageUInt8Ptr image = makeBinaryImage(pattern);
        const auto maxTree = MorphologicalTreeFactory::createMaxTree(image, 1.5);
        const auto minTree = MorphologicalTreeFactory::createMinTree(image, 1.5);
        const auto treeOfShapes = MorphologicalTreeFactory::createTreeOfShapes(image);

        validateTree(maxTree.topology(), "max_tree", pattern, summary);
        validateTree(minTree.topology(), "min_tree", pattern, summary);
        validateTree(treeOfShapes.topology(), "tree_of_shapes", pattern, summary);
        ++summary.images;
    }

    std::cout << "binary_patterns=" << summary.images << '\n'
              << "hierarchy_instances=" << summary.trees << '\n'
              << "validated_nodes=" << summary.nodes << '\n'
              << "mismatched_nodes=" << summary.mismatchedNodes << '\n'
              << "maximum_absolute_family_error=" << summary.maximumAbsoluteFamilyError << '\n';
    return summary.mismatchedNodes == 0 ? 0 : 1;
}
