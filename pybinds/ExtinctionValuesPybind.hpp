#pragma once

#include "../mmcfilters/filters/ExtinctionValues.hpp"
#include "PybindConversions.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <concepts>
#include <type_traits>
#include <variant>
namespace mmcfilters {

namespace py = pybind11;

/**
 * @brief Python-facing extinction selection policy converted to the native Real type on dispatch.
 */
struct ExtinctionSelectionPolicyPybind {
    /** @brief Enumerates the supported kind values. */
    enum class Kind { TopK, MinimumExtinction };

    /** @brief Stores the kind. */
    Kind kind = Kind::TopK;
    /** @brief Stores the extrema to keep. */
    int extremaToKeep = 0;
    /** @brief Stores the threshold. */
    double threshold = 0.0;

    /**
     * @brief Creates a policy that retains the top-k ranked extrema.
     *
     * @param extremaToKeep Number of highest-ranked extrema to retain.
     * @return Configured top-k selection policy.
     */
    static ExtinctionSelectionPolicyPybind byTopK(int extremaToKeep) {
        ExtinctionSelectionPolicyPybind policy;
        policy.kind = Kind::TopK;
        policy.extremaToKeep = extremaToKeep;
        return policy;
    }

    /**
     * @brief Creates a policy that retains extrema above a minimum extinction value.
     *
     * @param threshold Threshold applied by the operation.
     * @return Configured minimum-extinction selection policy.
     */
    static ExtinctionSelectionPolicyPybind byThreshold(double threshold) {
        ExtinctionSelectionPolicyPybind policy;
        policy.kind = Kind::MinimumExtinction;
        policy.threshold = threshold;
        return policy;
    }

    /**
     * @brief Converts to native.
     *
     * @return Converted to native.
     */
    template <std::floating_point Real> [[nodiscard]] ExtinctionSelectionPolicy<Real> toNative() const {
        switch (kind) {
        case Kind::TopK:
            return ExtinctionSelectionPolicy<Real>::byTopK(extremaToKeep);
        case Kind::MinimumExtinction:
            return ExtinctionSelectionPolicy<Real>::byThreshold(static_cast<Real>(threshold));
        }
        ExtinctionSelectionPolicy<Real> invalid;
        return invalid;
    }
};

/**
 * @brief Pybind11 wrapper exposing extinction-value computation to Python.
 */
class ExtinctionValuesPybind {
    /** @brief References the weighted owner used by the component. */
    std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weightedOwner_;
    /** @brief Stores the extinction. */
    std::variant<ExtinctionValues<std::uint8_t, float>, ExtinctionValues<std::uint8_t, double>> extinction_;

    /**
     * @brief Creates an extinction-value computer for the runtime attribute type.
     *
     * @param weighted Weighted tree used by the operation.
     * @param attribute Attribute requested by the operation.
     * @return Extinction-value computer variant for the array element type.
     */
    template <std::floating_point Real>
    static ExtinctionValues<std::uint8_t, Real> makeExtinction(WeightedMorphologicalTree<std::uint8_t>& weighted, py::array attribute) {
        auto typed = pybind_utils::requireNodeAttributeArray<Real>(std::move(attribute), weighted.topology(), "attribute");
        return ExtinctionValues<std::uint8_t, Real>(weighted, pybind_utils::toSharedPtr<Real>(typed));
    }


    /**
     * @brief Creates an extinction-value computer for the runtime attribute type.
     *
     * @param weighted Weighted tree used by the operation.
     * @param attribute Attribute requested by the operation.
     * @return Extinction-value computer variant for the array element type.
     */
    static std::variant<ExtinctionValues<std::uint8_t, float>, ExtinctionValues<std::uint8_t, double>>
    makeExtinction(WeightedMorphologicalTree<std::uint8_t>& weighted, py::array attribute) {
        if (pybind_utils::parseFloatingArrayDType(attribute, "attribute") == pybind_utils::FloatingDType::Float64) {
            return makeExtinction<double>(weighted, std::move(attribute));
        }
        return makeExtinction<float>(weighted, std::move(attribute));
    }

  public:
    /**
     * @brief Constructs `ExtinctionValuesPybind` from the supplied inputs.
     *
     * @param weighted Weighted tree used by the operation.
     * @param attribute Attribute requested by the operation.
     */
    ExtinctionValuesPybind(std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weighted, py::array attribute)
        : weightedOwner_(std::move(weighted)), extinction_(makeExtinction(*weightedOwner_, std::move(attribute))) {}

    /**
     * @brief Returns map.
     *
     * @param selection Policy used to select extrema from the ranked candidates.
     * @param scorePolicy Policy used to score and rank candidate extrema.
     * @return Map.
     */
    py::array contourMap(const ExtinctionSelectionPolicyPybind& selection, ExtinctionContourScorePolicy scorePolicy) {
        return std::visit(
            [&selection, scorePolicy](auto& extinction) -> py::array {
                using Extinction = std::decay_t<decltype(extinction)>;
                using Real = typename Extinction::value_type;
                return pybind_utils::toNumpy(extinction.contourMap(selection.toNative<Real>(), scorePolicy));
            },
            extinction_);
    }

    /**
     * @brief Returns regional extrema py.
     *
     * @return Regional extrema py.
     */
    std::vector<py::tuple> getRegionalExtremaPy() {
        return std::visit(
            [](auto& extinction) {
                const auto& vec = extinction.getRegionalExtrema();
                std::vector<py::tuple> out;
                out.reserve(vec.size());
                for (const auto& item : vec) {
                    out.push_back(py::make_tuple(item.leaf, item.cutoffNode, item.extinction));
                }
                return out;
            },
            extinction_);
    }


    /**
     * @brief Applies filtering.
     *
     * @param selection Policy used to select extrema from the ranked candidates.
     * @return Image or array produced by the operation.
     */
    py::array_t<uint8_t> filtering(const ExtinctionSelectionPolicyPybind& selection) {
        return std::visit(
            [&selection](auto& extinction) {
                using Extinction = std::decay_t<decltype(extinction)>;
                using Real = typename Extinction::value_type;
                return pybind_utils::toNumpy(extinction.filtering(selection.toNative<Real>()));
            },
            extinction_);
    }
};

} // namespace mmcfilters
