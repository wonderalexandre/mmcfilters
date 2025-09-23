#pragma once

#include "../mmcfilters/filters/AttributeOpeningPrimitivesFamily.hpp"

#include "MorphologicalTreePybind.hpp"
#include "PybindUtils.hpp"
#include "../mmcfilters/utils/Common.hpp"

#include <vector>
#include <list>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace mmcfilters {

/**
 * @brief Adaptador Pybind para manipular famílias de primitivas de abertura.
 */
class AttributeOpeningPrimitivesFamilyPybind: public AttributeOpeningPrimitivesFamily{

    public:
    using AttributeOpeningPrimitivesFamily::AttributeOpeningPrimitivesFamily;

    AttributeOpeningPrimitivesFamilyPybind(MorphologicalTreePybindPtr tree, py::array_t<float>& attr, float maxCriterion)
        : AttributeOpeningPrimitivesFamily(tree, PybindUtils::toShared_ptr(attr), maxCriterion){ }
          
      

    AttributeOpeningPrimitivesFamilyPybind(MorphologicalTreePybindPtr tree, py::array_t<float>& attr, float maxCriterion, int deltaMSER)
        : AttributeOpeningPrimitivesFamily(tree, PybindUtils::toShared_ptr(attr), maxCriterion, deltaMSER){ }

    py::array_t<uint8_t> getPrimitive(float threshold){
        ImageUInt8Ptr imgOut = ImageUInt8::create(this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage());
        AttributeFilters::filteringByPruningMin(this->tree, this->attrs_increasing, threshold, imgOut);
        return PybindUtils::toNumpy(imgOut);

    }

    py::array_t<uint8_t> getRestOfNumpyImage(){
        return PybindUtils::toNumpy(this->restOfImage);
    }

};

} // namespace mmcfilters
