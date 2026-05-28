#pragma once

#include "../mmcfilters/filters/ExtinctionValues.hpp"
#include "MorphologicalTreePybind.hpp"
#include "PybindUtils.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <concepts>
#include <variant>
namespace mmcfilters {

namespace py = pybind11;

/**
 * @brief Pybind11 wrapper exposing extinction-value computation to Python.
 */
class ExtinctionValuesPybind {
    std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weightedOwner_;
    std::variant<ExtinctionValues<std::uint8_t, float>, ExtinctionValues<std::uint8_t, double>> extinction_;

    template <std::floating_point Real>
    static ExtinctionValues<std::uint8_t, Real> makeExtinction(WeightedMorphologicalTree<std::uint8_t>& weighted, py::array attribute) {
        auto typed = PybindUtils::requireNodeAttributeArray<Real>(std::move(attribute), weighted.topology(), "attribute");
        return ExtinctionValues<std::uint8_t, Real>(weighted, PybindUtils::toSharedPtr<Real>(typed));
    }

    static std::variant<ExtinctionValues<std::uint8_t, float>, ExtinctionValues<std::uint8_t, double>> makeExtinction(WeightedMorphologicalTree<std::uint8_t>& weighted, py::array attribute) {
        if (PybindUtils::parseFloatingArrayDType(attribute, "attribute") == PybindUtils::FloatingDType::Float64) {
            return makeExtinction<double>(weighted, std::move(attribute));
        }
        return makeExtinction<float>(weighted, std::move(attribute));
    }

public:
    ExtinctionValuesPybind(std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weighted, py::array attribute)
        : weightedOwner_(std::move(weighted)),
          extinction_(makeExtinction(*weightedOwner_, std::move(attribute))) { }

    py::array saliencyMap(int leafToKeep, bool unweighted=true) {
        return std::visit(
            [leafToKeep, unweighted](auto& extinction) -> py::array {
                return PybindUtils::toNumpy(extinction.saliencyMap(leafToKeep, unweighted));
            },
            extinction_);
    }

    std::vector<py::tuple> getExtinctionValuesPy()  {
        return std::visit(
            [](auto& extinction) {
                auto& vec = extinction.getExtinctionValues();
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
            },
            extinction_);
    }

    py::array_t<uint8_t> filtering(int leafToKeep) {
        return std::visit(
            [leafToKeep](auto& extinction) {
                return PybindUtils::toNumpy(extinction.filtering(leafToKeep));
            },
            extinction_);
    }

};

} // namespace mmcfilters
