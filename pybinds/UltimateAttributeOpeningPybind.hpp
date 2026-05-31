#pragma once
#include <concepts>
#include <variant>
#include "../mmcfilters/filters/UltimateAttributeOpening.hpp"
#include "MorphologicalTreePybind.hpp"
#include "PybindUtils.hpp"
#include <pybind11/pybind11.h>
namespace mmcfilters {


namespace py = pybind11;

/**
 * @brief Pybind11 wrapper for Ultimate Attribute Opening.
 */
class UltimateAttributeOpeningPybind {
    std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weightedOwner_;
    std::variant<UltimateAttributeOpening<std::uint8_t, float>, UltimateAttributeOpening<std::uint8_t, double>> uao_;

    template <std::floating_point Real>
    static UltimateAttributeOpening<std::uint8_t, Real> makeUao(WeightedMorphologicalTree<std::uint8_t>& weighted, py::array attr) {
        auto typed = PybindUtils::requireNodeAttributeArray<Real>(std::move(attr), weighted.topology());
        return UltimateAttributeOpening<std::uint8_t, Real>(weighted, PybindUtils::toSharedPtr<Real>(typed));
    }

    static std::variant<UltimateAttributeOpening<std::uint8_t, float>, UltimateAttributeOpening<std::uint8_t, double>> makeUao(WeightedMorphologicalTree<std::uint8_t>& weighted, py::array attr) {
        if (PybindUtils::parseFloatingArrayDType(attr, "attr") == PybindUtils::FloatingDType::Float64) {
            return makeUao<double>(weighted, std::move(attr));
        }
        return makeUao<float>(weighted, std::move(attr));
    }

public:
    UltimateAttributeOpeningPybind(std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weighted, py::array attr)
        : weightedOwner_(std::move(weighted)),
          uao_(makeUao(*weightedOwner_, std::move(attr))) {}

    void execute(double maxCriterion) {
        std::visit(
            [maxCriterion](auto& uao) {
                using Uao = std::decay_t<decltype(uao)>;
                using Real = typename Uao::attribute_value_type;
                uao.execute(static_cast<Real>(maxCriterion));
            },
            uao_);
    }

    void executeWithMSER(double maxCriterion, int deltaMSER) {
        std::visit(
            [maxCriterion, deltaMSER](auto& uao) {
                using Uao = std::decay_t<decltype(uao)>;
                using Real = typename Uao::attribute_value_type;
                uao.executeWithMSER(static_cast<Real>(maxCriterion), deltaMSER);
            },
            uao_);
    }

    void executeWithDepthStability(double maxCriterion, int depthDelta) {
        std::visit(
            [maxCriterion, depthDelta](auto& uao) {
                using Uao = std::decay_t<decltype(uao)>;
                using Real = typename Uao::attribute_value_type;
                uao.executeWithDepthStability(static_cast<Real>(maxCriterion), depthDelta);
            },
            uao_);
    }

    py::array_t<uint8_t> getMaxContrastImage(){
        return std::visit(
            [](auto& uao) {
                return PybindUtils::toNumpy(uao.getMaxContrastImage());
            },
            uao_);
    }       

    py::array_t<int32_t> getAssociatedImage(){
        return std::visit(
            [](auto& uao) {
                return PybindUtils::toNumpy(uao.getAssociatedImage());
            },
            uao_);
    }
    py::array_t<uint8_t> getAssociatedColorImage(){
        return std::visit(
            [](auto& uao) {
                return PybindUtils::toNumpy(uao.getAssociatedColorImage());
            },
            uao_);
    }



};

} // namespace mmcfilters
