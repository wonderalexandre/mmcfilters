#pragma once


#include <array>
#include "../mmcfilters/filters/ComputerMSER.hpp"
#include "../mmcfilters/attributes/AttributeComputation.hpp"
#include "../mmcfilters/filters/UltimateAttributeOpening.hpp"
#include "MorphologicalTreePybind.hpp"
#include "PybindUtils.hpp"
#include <pybind11/pybind11.h>
namespace mmcfilters {


namespace py = pybind11;

/**
 * @brief Pybind11 wrapper for Ultimate Attribute Opening.
 */
class UltimateAttributeOpeningPybind : public UltimateAttributeOpening<std::uint8_t> {
    using FloatArray = py::array_t<float, py::array::c_style | py::array::forcecast>;

    std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weightedOwner_;

public:
    using UltimateAttributeOpening<std::uint8_t>::UltimateAttributeOpening;

    UltimateAttributeOpeningPybind(std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weighted, FloatArray attr) :
        UltimateAttributeOpening<std::uint8_t>(
            *weighted,
            [&]() {
                PybindUtils::require1DArray(attr.request(), weighted->topology().getNumInternalNodeSlots(), "attr");
                return PybindUtils::toShared_ptr(attr);
            }()),
        weightedOwner_(std::move(weighted)) {}

    py::array_t<uint8_t> getMaxContrastImage(){
        return PybindUtils::toNumpy(UltimateAttributeOpening<std::uint8_t>::getMaxContrastImage());
    }       

    py::array_t<int32_t> getAssociatedImage(){
        auto imgOut = UltimateAttributeOpening<std::uint8_t>::getAssociatedImage();
        return PybindUtils::toNumpy(imgOut);
    }
    py::array_t<uint8_t> getAssociatedColorImage(){
        return PybindUtils::toNumpy(UltimateAttributeOpening<std::uint8_t>::getAssociatedColorImage());
    }



};

} // namespace mmcfilters
