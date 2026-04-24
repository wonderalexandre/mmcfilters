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
    using FloatArray = py::array_t<float, py::array::c_style | py::array::forcecast>;

    MorphologicalTreePybindPtr treeOwner_;
    std::shared_ptr<WeightedMorphologicalTree> weightedOwner_;

    public:
    using AttributeOpeningPrimitivesFamily::AttributeOpeningPrimitivesFamily;

    AttributeOpeningPrimitivesFamilyPybind(MorphologicalTreePybindPtr tree, FloatArray attr, float maxCriterion)
        : AttributeOpeningPrimitivesFamily(
            *tree,
            [&]() {
                PybindUtils::require1DArray(attr.request(), tree->getNumInternalNodeSlots(), "attr");
                return PybindUtils::toShared_ptr(attr);
            }(),
            maxCriterion),
          treeOwner_(std::move(tree)) { }

    AttributeOpeningPrimitivesFamilyPybind(std::shared_ptr<WeightedMorphologicalTree> weighted, FloatArray attr, float maxCriterion)
        : AttributeOpeningPrimitivesFamily(
            *weighted,
            [&]() {
                PybindUtils::require1DArray(attr.request(), weighted->topology().getNumInternalNodeSlots(), "attr");
                return PybindUtils::toShared_ptr(attr);
            }(),
            maxCriterion),
          weightedOwner_(std::move(weighted)) { }
          
      

    AttributeOpeningPrimitivesFamilyPybind(MorphologicalTreePybindPtr tree, FloatArray attr, float maxCriterion, int deltaMSER)
        : AttributeOpeningPrimitivesFamily(
            *tree,
            [&]() {
                PybindUtils::require1DArray(attr.request(), tree->getNumInternalNodeSlots(), "attr");
                return PybindUtils::toShared_ptr(attr);
            }(),
            maxCriterion,
            deltaMSER),
          treeOwner_(std::move(tree)) { }

    AttributeOpeningPrimitivesFamilyPybind(std::shared_ptr<WeightedMorphologicalTree> weighted, FloatArray attr, float maxCriterion, int deltaMSER)
        : AttributeOpeningPrimitivesFamily(
            *weighted,
            [&]() {
                PybindUtils::require1DArray(attr.request(), weighted->topology().getNumInternalNodeSlots(), "attr");
                return PybindUtils::toShared_ptr(attr);
            }(),
            maxCriterion,
            deltaMSER),
          weightedOwner_(std::move(weighted)) { }

    py::array_t<uint8_t> getPrimitive(float threshold){
        return PybindUtils::toNumpy(this->filters_->filteringByPruningMin(this->attrs_increasing, threshold));
    }

    py::array_t<uint8_t> getRestOfNumpyImage(){
        return PybindUtils::toNumpy(this->restOfImage);
    }

};

} // namespace mmcfilters
