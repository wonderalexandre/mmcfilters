
#ifndef ULTIMATE_ATTR_OPENING_PYBIND_H
#define ULTIMATE_ATTR_OPENING_PYBIND_H

#include <array>
#include "../include/NodeCT.hpp"
#include "../include/ComputerMSER.hpp"
#include "../include/AttributeComputedIncrementally.hpp"
#include "../include/UltimateAttributeOpening.hpp"

#include "../pybind/ComponentTreePybind.hpp"
#include "../pybind/PybindUtils.hpp"
#include <pybind11/pybind11.h>


namespace py = pybind11;

class UltimateAttributeOpeningPybind: public UltimateAttributeOpening{

public:
    using UltimateAttributeOpening::UltimateAttributeOpening;

    UltimateAttributeOpeningPybind(ComponentTreePybind* tree,  std::vector<float> attrs_increasing) : UltimateAttributeOpening(tree, attrs_increasing){}

    py::array_t<int> getMaxConstrastImage(){
        return PybindUtils::toNumpy(UltimateAttributeOpening::getMaxConstrastImage(), this->tree->getNumColsOfImage() * this->tree->getNumRowsOfImage());
    }       

    py::array_t<int> getAssociatedImage(){
        return PybindUtils::toNumpy(UltimateAttributeOpening::getAssociatedImage(), this->tree->getNumColsOfImage() * this->tree->getNumRowsOfImage());
    }
    py::array_t<int> getAssociatedColorImage(){
        return PybindUtils::toNumpy(UltimateAttributeOpening::getAssociatedColorImage(), this->tree->getNumColsOfImage() * this->tree->getNumRowsOfImage() * 3);
    }



};

#endif