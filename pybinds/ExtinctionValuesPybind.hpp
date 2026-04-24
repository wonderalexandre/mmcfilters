#pragma once

#include "../mmcfilters/filters/ExtinctionValues.hpp"
#include "MorphologicalTreePybind.hpp"
#include "PybindUtils.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <algorithm> 
#include <cmath>
#include <iostream>
namespace mmcfilters {

namespace py = pybind11;

/**
 * @brief *Wrapper* Pybind11 para cálculo e uso de valores de extinção.
 */
class ExtinctionValuesPybind : public ExtinctionValues{
    using FloatArray = py::array_t<float, py::array::c_style | py::array::forcecast>;

    MorphologicalTreePybindPtr treeOwner_;
    std::shared_ptr<WeightedMorphologicalTree> weightedOwner_;

    public:
    using ExtinctionValues::ExtinctionValues;
    
    ExtinctionValuesPybind(MorphologicalTreePybindPtr tree, FloatArray attribute)
        : ExtinctionValues(
            *tree,
            [&]() {
                PybindUtils::require1DArray(attribute.request(), tree->getNumInternalNodeSlots(), "attribute");
                return PybindUtils::toShared_ptr(attribute);
            }()),
          treeOwner_(std::move(tree)) { }

    ExtinctionValuesPybind(std::shared_ptr<WeightedMorphologicalTree> weighted, FloatArray attribute)
        : ExtinctionValues(
            *weighted,
            [&]() {
                PybindUtils::require1DArray(attribute.request(), weighted->tree.getNumInternalNodeSlots(), "attribute");
                return PybindUtils::toShared_ptr(attribute);
            }()),
          weightedOwner_(std::move(weighted)) { }

    py::array_t<float> saliencyMap(int leafToKeep, bool unweighted=true) {

        auto saliencyMapPtr = ExtinctionValues::saliencyMap(leafToKeep, unweighted);
        return PybindUtils::toNumpy(saliencyMapPtr);
    }

    std::vector<py::tuple> getExtinctionValuesPy()  {
        auto &vec = ExtinctionValues::getExtinctionValues();
        std::vector<py::tuple> out;
        out.reserve(vec.size());
        for (const auto &item : vec) {
            out.push_back(py::make_tuple(
                item.leaf,
                item.cutoffNode,
                item.extinction
            ));
        }
        return out;
    }

    py::array_t<uint8_t> filtering(int leafToKeep) {

        ImageUInt8Ptr filteredImagePtr =  ExtinctionValues::filtering(leafToKeep);
        return PybindUtils::toNumpy(filteredImagePtr);
    }

};

} // namespace mmcfilters
