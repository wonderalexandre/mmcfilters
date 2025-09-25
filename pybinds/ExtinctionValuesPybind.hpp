#pragma once

#include "../mmcfilters/filters/ExtinctionValues.hpp"
#include "../mmcfilters/trees/NodeMT.hpp"


#include "MorphologicalTreePybind.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <algorithm> 
#include <cmath>
#include <iostream>
namespace mmcfilters {

namespace py = pybind11;


class ExtinctionValuesPybind;
using ExtinctionValuesPybindPtr = std::shared_ptr<ExtinctionValuesPybind>;

/**
 * @brief *Wrapper* Pybind11 para cálculo e uso de valores de extinção.
 */
class ExtinctionValuesPybind : public ExtinctionValues{

    public:
    using ExtinctionValues::ExtinctionValues;
    
    ExtinctionValuesPybind(MorphologicalTreePybindPtr tree, py::array_t<float>& attribute)
        : ExtinctionValues(tree, PybindUtils::toShared_ptr(attribute)) { }

    py::array_t<float> saliencyMap(int leafToKeep, bool unweighted=true) {

        auto saliencyMapPtr = ExtinctionValues::saliencyMap(leafToKeep, unweighted);
        return PybindUtils::toNumpy(saliencyMapPtr);
    }

    // Return a Python list of (leaf, cutoffNode, extinction) tuples for easy unpacking
    std::vector<py::tuple> getExtinctionValuesPy()  {
        auto &vec = ExtinctionValues::getExtinctionValues();
        std::vector<py::tuple> out;
        out.reserve(vec.size());
        for (const auto &item : vec) {
            out.push_back(py::make_tuple(this->tree->proxy(item.leaf), tree->proxy(item.cutoffNode), item.extinction));
        }
        return out;
    }

    py::array_t<uint8_t> filtering(int leafToKeep) {

        ImageUInt8Ptr filteredImagePtr =  ExtinctionValues::filtering(leafToKeep);
        return PybindUtils::toNumpy(filteredImagePtr);
    }

};

} // namespace mmcfilters
