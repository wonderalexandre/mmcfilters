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
 * @brief Pybind11 wrapper exposing extinction-value computation to Python.
 */
class ExtinctionValuesPybind : public ExtinctionValues<std::uint8_t> {
    using FloatArray = py::array_t<float, py::array::c_style | py::array::forcecast>;

    std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weightedOwner_;

    public:
    using ExtinctionValues<std::uint8_t>::ExtinctionValues;

    ExtinctionValuesPybind(std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weighted, FloatArray attribute)
        : ExtinctionValues<std::uint8_t>(
            *weighted,
            [&]() {
                PybindUtils::require1DArray(attribute.request(), weighted->topology().getNumInternalNodeSlots(), "attribute");
                return PybindUtils::toShared_ptr(attribute);
            }()),
          weightedOwner_(std::move(weighted)) { }

    py::array_t<float> saliencyMap(int leafToKeep, bool unweighted=true) {

        auto saliencyMapPtr = ExtinctionValues<std::uint8_t>::saliencyMap(leafToKeep, unweighted);
        return PybindUtils::toNumpy(saliencyMapPtr);
    }

    std::vector<py::tuple> getExtinctionValuesPy()  {
        auto &vec = ExtinctionValues<std::uint8_t>::getExtinctionValues();
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

        ImageUInt8Ptr filteredImagePtr =  ExtinctionValues<std::uint8_t>::filtering(leafToKeep);
        return PybindUtils::toNumpy(filteredImagePtr);
    }

};

} // namespace mmcfilters
