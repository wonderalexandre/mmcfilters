#include "../support/TestSupport.hpp"

#include "mmcfilters/trees/adjust/CasfComponentTrees.hpp"
#include "mmcfilters/trees/adjust/DualMinMaxTreeIncrementalFilter.hpp"

#include <iostream>
#include <string>

using namespace mmcfilters;
using namespace mmcfilters::adjust;
using namespace mmcfilters::unit_tests;

namespace {

ImageUInt8Ptr makeCasfFixture() {
    return makeImage(
        6,
        6,
        {
            2, 2, 2, 2, 1, 1,
            2, 5, 5, 2, 3, 1,
            2, 5, 6, 2, 3, 1,
            2, 5, 5, 2, 3, 1,
            3, 3, 3, 3, 3, 1,
            2, 2, 2, 1, 1, 1,
        });
}

void requireValidWeightedTree(const WeightedMorphologicalTree& tree, const std::string& label) {
    tree.topology().validateConnectedRootedTree();
    tree.validateMonotoneAltitude();
    require(tree.topology().getRoot() != InvalidNode, label + " must keep a valid root");
}

void test_empty_and_zero_threshold_preserve_image() {
    const auto image = makeCasfFixture();
    CasfComponentTrees casf(image, CasfComponentTreesAttribute::AREA);

    const auto filteredEmpty = casf.filter({});
    requireVectorEqual(collectImageValues(filteredEmpty), collectImageValues(image), "empty threshold CASF image");

    const auto filteredZero = casf.filter({0.0});
    requireVectorEqual(collectImageValues(filteredZero), collectImageValues(image), "zero threshold CASF image");
    requireValidWeightedTree(casf.minTree(), "CASF min-tree after zero threshold");
    requireValidWeightedTree(casf.maxTree(), "CASF max-tree after zero threshold");
}

void test_positive_threshold_keeps_exportable_trees() {
    const auto image = makeCasfFixture();
    CasfComponentTrees casf(image, CasfComponentTreesAttribute::AREA);

    const auto filtered = casf.filter({2.0});
    requireImageShape(filtered, 6, 6);
    requireValidWeightedTree(casf.minTree(), "CASF min-tree after positive threshold");
    requireValidWeightedTree(casf.maxTree(), "CASF max-tree after positive threshold");

    const auto [minParent, minAltitude] = casf.exportMinTree();
    const auto [maxParent, maxAltitude] = casf.exportMaxTree();

    requireEqual(static_cast<int>(minParent.size()), static_cast<int>(minAltitude.size()), "exported min-tree size");
    requireEqual(static_cast<int>(maxParent.size()), static_cast<int>(maxAltitude.size()), "exported max-tree size");
    require(static_cast<int>(minParent.size()) >= image->getSize(), "exported min-tree must include all leaves");
    require(static_cast<int>(maxParent.size()) >= image->getSize(), "exported max-tree must include all leaves");
}

void test_bounding_box_attribute_path_executes() {
    const auto image = makeCasfFixture();
    CasfComponentTrees casf(image, CasfComponentTreesAttribute::BOUNDING_BOX_DIAGONAL);

    const auto filtered = casf.filter({2.0});
    requireImageShape(filtered, 6, 6);
    requireValidWeightedTree(casf.minTree(), "CASF bbox min-tree after positive threshold");
    requireValidWeightedTree(casf.maxTree(), "CASF bbox max-tree after positive threshold");
    static_cast<void>(casf.exportMinTree());
    static_cast<void>(casf.exportMaxTree());
}

} // namespace

int main() {
    try {
        test_empty_and_zero_threshold_preserve_image();
        test_positive_threshold_keeps_exportable_trees();
        test_bounding_box_attribute_path_executes();
    } catch (const std::exception& ex) {
        std::cerr << "component tree adjust test failed: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}
