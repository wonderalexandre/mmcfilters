#pragma once


#include <array>
#include "../mmcfilters/attributes/ComputerMSER.hpp"
#include "../mmcfilters/attributes/AttributeComputedIncrementally.hpp"
#include "../mmcfilters/filters/UltimateAttributeOpening.hpp"
#include "../mmcfilters/utils/Common.hpp"
#include "MorphologicalTreePybind.hpp"
#include "PybindUtils.hpp"
#include <pybind11/pybind11.h>
namespace mmcfilters {


namespace py = pybind11;

/**
 * @brief Adaptação Pybind da Ultimate Attribute Opening.
 */
class UltimateAttributeOpeningPybind: public UltimateAttributeOpening{
    using FloatArray = py::array_t<float, py::array::c_style | py::array::forcecast>;

    MorphologicalTreePybindPtr treeOwner_;
    std::shared_ptr<WeightedMorphologicalTree> weightedOwner_;

public:
    using UltimateAttributeOpening::UltimateAttributeOpening;

    UltimateAttributeOpeningPybind(MorphologicalTreePybindPtr tree, FloatArray attr) : 
        UltimateAttributeOpening(
            *tree,
            [&]() {
                PybindUtils::require1DArray(attr.request(), tree->getNumInternalNodeSlots(), "attr");
                return PybindUtils::toShared_ptr(attr);
            }()),
        treeOwner_(std::move(tree)) {}

    UltimateAttributeOpeningPybind(std::shared_ptr<WeightedMorphologicalTree> weighted, FloatArray attr) :
        UltimateAttributeOpening(
            *weighted,
            [&]() {
                PybindUtils::require1DArray(attr.request(), weighted->topology().getNumInternalNodeSlots(), "attr");
                return PybindUtils::toShared_ptr(attr);
            }()),
        weightedOwner_(std::move(weighted)) {}

    py::array_t<uint8_t> getMaxContrastImage(){
        return PybindUtils::toNumpy(UltimateAttributeOpening::getMaxContrastImage());
    }       

    py::array_t<int32_t> getAssociatedImage(){
        auto imgOut = UltimateAttributeOpening::getAssociatedImage();
        return PybindUtils::toNumpy(imgOut);
    }
    py::array_t<uint8_t> getAssociatedColorImage(){
        return PybindUtils::toNumpy(UltimateAttributeOpening::getAssociatedColorImage());
    }



};

} // namespace mmcfilters
