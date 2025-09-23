#pragma once


#include <array>
#include "../mmcfilters/trees/NodeMT.hpp"
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

public:
    using UltimateAttributeOpening::UltimateAttributeOpening;

    UltimateAttributeOpeningPybind(MorphologicalTreePybindPtr tree,  py::array_t<float> &attr) : 
        UltimateAttributeOpening(tree, PybindUtils::toShared_ptr(attr) ){}

    py::array_t<uint8_t> getMaxConstrastImage(){
        return PybindUtils::toNumpy(UltimateAttributeOpening::getMaxConstrastImage());
    }       

    py::array_t<int32_t> getAssociatedImage(){
        auto imgOut = UltimateAttributeOpening::getAssociatedImage();
        return PybindUtils::toNumpyInt(imgOut->rawData(), imgOut->getSize());
    }
    py::array_t<uint8_t> getAssociatedColorImage(){
        return PybindUtils::toNumpy(UltimateAttributeOpening::getAssociatedColorImage());
    }



};

} // namespace mmcfilters
