#pragma once

#include "../mmcfilters/trees/NodeMT.hpp"
#include "../mmcfilters/filters/AttributeFilters.hpp"
#include "../mmcfilters/utils/Common.hpp"

#include "MorphologicalTreePybind.hpp"
#include "AttributeComputedIncrementallyPybind.hpp"
#include "PybindUtils.hpp"

#include <stack>
#include <vector>
#include <limits.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace mmcfilters {

#define UNDEF -999999999999

/**
 * @brief *Wrapper* Pybind11 para expor filtragens por atributos ao Python.
 */
class AttributeFiltersPybind : public AttributeFilters{


    public:
    using AttributeFilters::AttributeFilters;

    AttributeFiltersPybind(MorphologicalTreePybindPtr tree): AttributeFilters(tree){}

    py::array_t<uint8_t> filteringByPruningMin(py::array_t<float> &attr, float threshold){

        std::shared_ptr<float[]> attribute = PybindUtils::toShared_ptr(attr);

        ImageUInt8Ptr imgOutput = ImageUInt8::create(this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage());

        AttributeFilters::filteringByPruningMin(this->tree, attribute, threshold, imgOutput);

        return PybindUtils::toNumpy(imgOutput);
    }

    py::array_t<uint8_t> filteringByPruningMax(py::array_t<float> &attr, float threshold){

        std::shared_ptr<float[]> attribute = PybindUtils::toShared_ptr(attr);

        ImageUInt8Ptr imgOutput = ImageUInt8::create(this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage());
        AttributeFilters::filteringByPruningMax(this->tree, attribute, threshold, imgOutput);

        return PybindUtils::toNumpy(imgOutput);

    }

    py::array_t<uint8_t> filteringByPruningMin(std::vector<bool>& criterion){
        ImageUInt8Ptr imgOutput = ImageUInt8::create(this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage());
        AttributeFilters::filteringByPruningMin(this->tree, criterion, imgOutput);

        return PybindUtils::toNumpy(imgOutput);
    }

    py::array_t<uint8_t> filteringByDirectRule(std::vector<bool>& criterion){
        ImageUInt8Ptr imgOutput = ImageUInt8::create(this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage());
        AttributeFilters::filteringByDirectRule(this->tree, criterion, imgOutput);

        return PybindUtils::toNumpy(imgOutput);
    }

    py::array_t<uint8_t> filteringByPruningMax(std::vector<bool>& criterion){
        ImageUInt8Ptr imgOutput = ImageUInt8::create(this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage());
        AttributeFilters::filteringByPruningMax(this->tree, criterion, imgOutput);

        return PybindUtils::toNumpy(imgOutput);

    }

    std::vector<bool> getAdaptativeCriterion(std::vector<bool>& criterion, int delta){
        return AttributeFilters::getAdaptativeCriterion(criterion, delta);
    }



    py::array_t<uint8_t> filteringBySubtractiveRule(std::vector<bool>& criterion){
        ImageUInt8Ptr imgOutput = ImageUInt8::create(this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage());
        AttributeFilters::filteringBySubtractiveRule(this->tree, criterion, imgOutput);

        return PybindUtils::toNumpy(imgOutput);

    }

    py::array_t<float> filteringBySubtractiveScoreRule(std::vector<float>& prob){
        ImageFloatPtr imgOutput = ImageFloat::create(this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage());
        AttributeFilters::filteringBySubtractiveScoreRule(this->tree, prob, imgOutput);

        return PybindUtils::toNumpy(imgOutput);

    }

    py::array_t<uint8_t> filteringByExtinctionValue(py::array_t<float>& attr, int leafToKeep){
        ExtinctionValues ev(this->tree, PybindUtils::toShared_ptr(attr));
        ImageUInt8Ptr filteredImagePtr =  ev.filtering(leafToKeep);
        return PybindUtils::toNumpy(filteredImagePtr);
    }

    py::array_t<float> saliencyMapByExtinctionValue(py::array_t<float>& attr, int leafToKeep, bool unweighted=false){
        
        ExtinctionValues ev(this->tree, PybindUtils::toShared_ptr(attr));
        auto saliencyMapPtr = ev.saliencyMap(leafToKeep, unweighted);
        return PybindUtils::toNumpy(saliencyMapPtr);
    }



};

} // namespace mmcfilters
