
#ifndef ULTIMATE_ATTR_OPENING_PYBIND_H
#define ULTIMATE_ATTR_OPENING_PYBIND_H

#include <array>
#include "../include/NodeMT.hpp"
#include "../include/ComputerMSER.hpp"
#include "../include/AttributeComputedIncrementally.hpp"
#include "../include/UltimateAttributeOpening.hpp"
#include "../include/Common.hpp"
#include "../pybind/MorphologicalTreePybind.hpp"
#include "../pybind/PybindUtils.hpp"
#include <pybind11/pybind11.h>


namespace py = pybind11;

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

#endif