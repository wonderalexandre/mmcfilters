#include "../include/AttributeOpeningPrimitivesFamily.hpp"

#include "../pybind/ComponentTreePybind.hpp"
#include "../pybind/PybindUtils.hpp"

#include <vector>
#include <list>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#ifndef ATTRIBUTE_OPENING_PRIMITIVES_FAMILY_PYBIND_H
#define ATTRIBUTE_OPENING_PRIMITIVES_FAMILY_PYBIND_H

class AttributeOpeningPrimitivesFamilyPybind: public AttributeOpeningPrimitivesFamily{

    public:
    using AttributeOpeningPrimitivesFamily::AttributeOpeningPrimitivesFamily;

    AttributeOpeningPrimitivesFamilyPybind(ComponentTreePybind* tree, py::array_t<float> attr, float maxCriterion)
        : AttributeOpeningPrimitivesFamily(tree, static_cast<float*>(attr.request().ptr), maxCriterion) {}

    AttributeOpeningPrimitivesFamilyPybind(ComponentTreePybind* tree, py::array_t<float> attr, float maxCriterion, int deltaMSER)
        : AttributeOpeningPrimitivesFamily(tree, static_cast<float*>(attr.request().ptr), maxCriterion, deltaMSER) {}

    py::array_t<int> getPrimitive(float threshold){
        int* imgOut = new int[this->tree->getNumRowsOfImage() * this->tree->getNumColsOfImage()];
        AttributeFilters::filteringByPruningMin(this->tree, this->attrs_increasing, threshold, imgOut);
        return PybindUtils::toNumpy(imgOut, this->tree->getNumRowsOfImage()*this->tree->getNumColsOfImage());

    }

    py::array_t<int> getRestOfNumpyImage(){
        return PybindUtils::toNumpy(this->restOfImage, this->tree->getNumRowsOfImage()*this->tree->getNumColsOfImage());
    }

};

#endif
