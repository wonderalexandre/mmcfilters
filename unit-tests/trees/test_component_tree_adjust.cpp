#include "../support/TestSupport.hpp"

#include "mmcfilters/trees/adjust/CasfComponentTrees.hpp"
#include "mmcfilters/trees/adjust/DualMinMaxTreeIncrementalFilter.hpp"

#include <cstdint>
#include <iostream>
#include <initializer_list>
#include <string>
#include <type_traits>

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

template<class T>
ImagePtr<T> makeTypedCasfFixture(std::initializer_list<T> values) {
    requireEqual(static_cast<int>(values.size()), 36, "typed CASF fixture buffer size");
    auto image = Image<T>::create(6, 6);
    int index = 0;
    for (T value : values) {
        (*image)[index++] = value;
    }
    return image;
}

template<class T>
void requireValidWeightedTree(const WeightedMorphologicalTree<T>& tree, const std::string& label) {
    tree.topology().validateConnectedRootedTree();
    tree.validateMonotoneAltitude();
    require(tree.topology().getRoot() != InvalidNode, label + " must keep a valid root");
}

void test_empty_and_zero_threshold_preserve_image() {
    const auto image = makeCasfFixture();
    CasfComponentTrees<std::uint8_t> casf(image, CasfComponentTreesAttribute::AREA);

    const auto filteredEmpty = casf.filter({});
    requireVectorEqual(collectImageValues(filteredEmpty), collectImageValues(image), "empty threshold CASF image");

    const auto filteredZero = casf.filter({0.0});
    requireVectorEqual(collectImageValues(filteredZero), collectImageValues(image), "zero threshold CASF image");
    requireValidWeightedTree(casf.minTree(), "CASF min-tree after zero threshold");
    requireValidWeightedTree(casf.maxTree(), "CASF max-tree after zero threshold");
}

void test_positive_threshold_keeps_exportable_trees() {
    const auto image = makeCasfFixture();
    CasfComponentTrees<std::uint8_t> casf(image, CasfComponentTreesAttribute::AREA);

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
    CasfComponentTrees<std::uint8_t> casf(image, CasfComponentTreesAttribute::BOUNDING_BOX_DIAGONAL);

    const auto filtered = casf.filter({2.0});
    requireImageShape(filtered, 6, 6);
    requireValidWeightedTree(casf.minTree(), "CASF bbox min-tree after positive threshold");
    requireValidWeightedTree(casf.maxTree(), "CASF bbox max-tree after positive threshold");
    static_cast<void>(casf.exportMinTree());
    static_cast<void>(casf.exportMaxTree());
}

void test_int32_casf_keeps_typed_altitude_state() {
    const auto image = makeTypedCasfFixture<std::int32_t>({
        -20, -20, -20, -20, 400, 400,
        -20, 120, 120, -20, 300, 400,
        -20, 120, 900, -20, 300, 400,
        -20, 120, 120, -20, 300, 400,
        300, 300, 300, 300, 300, 400,
        -10, -10, -10, 400, 400, 400,
    });
    CasfComponentTrees<std::int32_t> casf(image, CasfComponentTreesAttribute::AREA);

    const auto filteredEmpty = casf.filter({});
    requireVectorEqual(collectImageValues(filteredEmpty), collectImageValues(image), "empty int32 CASF image");

    const auto filtered = casf.filter({2.0});
    requireImageShape(filtered, 6, 6);
    requireValidWeightedTree(casf.minTree(), "typed int32 CASF min-tree after positive threshold");
    requireValidWeightedTree(casf.maxTree(), "typed int32 CASF max-tree after positive threshold");

    const auto [minParent, minAltitude] = casf.exportMinTree();
    const auto [maxParent, maxAltitude] = casf.exportMaxTree();
    static_assert(std::is_same_v<decltype(minAltitude), const std::vector<std::int32_t>>);
    static_assert(std::is_same_v<decltype(maxAltitude), const std::vector<std::int32_t>>);
    requireEqual(minParent.size(), minAltitude.size(), "typed int32 exported min-tree size");
    requireEqual(maxParent.size(), maxAltitude.size(), "typed int32 exported max-tree size");
}

void test_float_casf_uses_sparse_altitude_backend() {
    const auto image = makeTypedCasfFixture<float>({
        0.20f, 0.20f, 0.20f, 0.20f, 1.10f, 1.10f,
        0.20f, 0.55f, 0.55f, 0.20f, 0.80f, 1.10f,
        0.20f, 0.55f, 1.35f, 0.20f, 0.80f, 1.10f,
        0.20f, 0.55f, 0.55f, 0.20f, 0.80f, 1.10f,
        0.80f, 0.80f, 0.80f, 0.80f, 0.80f, 1.10f,
        0.25f, 0.25f, 0.25f, 1.10f, 1.10f, 1.10f,
    });
    static_assert(!DualMinMaxTreeIncrementalFilter<float>::usesDenseLevelBackend());

    CasfComponentTrees<float> casf(image, CasfComponentTreesAttribute::BOUNDING_BOX_DIAGONAL);
    const auto filteredEmpty = casf.filter({});
    requireVectorEqual(collectImageValues(filteredEmpty), collectImageValues(image), "empty float CASF image");

    const auto filtered = casf.filter({2.5});
    requireImageShape(filtered, 6, 6);
    requireValidWeightedTree(casf.minTree(), "typed float CASF min-tree after positive threshold");
    requireValidWeightedTree(casf.maxTree(), "typed float CASF max-tree after positive threshold");

    const auto [minParent, minAltitude] = casf.exportMinTree();
    const auto [maxParent, maxAltitude] = casf.exportMaxTree();
    static_assert(std::is_same_v<decltype(minAltitude), const std::vector<float>>);
    static_assert(std::is_same_v<decltype(maxAltitude), const std::vector<float>>);
    requireEqual(minParent.size(), minAltitude.size(), "typed float exported min-tree size");
    requireEqual(maxParent.size(), maxAltitude.size(), "typed float exported max-tree size");
}

} // namespace

int main() {
    try {
        test_empty_and_zero_threshold_preserve_image();
        test_positive_threshold_keeps_exportable_trees();
        test_bounding_box_attribute_path_executes();
        test_int32_casf_keeps_typed_altitude_state();
        test_float_casf_uses_sparse_altitude_backend();
    } catch (const std::exception& ex) {
        std::cerr << "component tree adjust test failed: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}
