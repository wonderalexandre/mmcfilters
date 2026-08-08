#pragma once
#include <concepts>
#include <variant>
#include "../mmcfilters/filters/UltimateAttributeOpening.hpp"
#include "PybindConversions.hpp"
#include <pybind11/pybind11.h>
namespace mmcfilters {

namespace py = pybind11;

/**
 * @brief Pybind11 wrapper for Ultimate Attribute Opening.
 */
class UltimateAttributeOpeningPybind {
    /** @brief References the weighted owner used by the component. */
    std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weightedOwner_;
    /** @brief Stores the uao. */
    std::variant<UltimateAttributeOpening<std::uint8_t, float>, UltimateAttributeOpening<std::uint8_t, double>> uao_;

    /**
     * @brief Creates an ultimate-attribute-opening computer for the runtime attribute type.
     *
     * @param weighted Weighted tree used by the operation.
     * @param attr Attribute requested by the operation.
     * @return Ultimate-attribute-opening computer variant for the array element type.
     */
    template <std::floating_point Real>
    static UltimateAttributeOpening<std::uint8_t, Real> makeUao(WeightedMorphologicalTree<std::uint8_t>& weighted, py::array attr) {
        auto typed = pybind_utils::requireNodeAttributeArray<Real>(std::move(attr), weighted.topology());
        return UltimateAttributeOpening<std::uint8_t, Real>(weighted, pybind_utils::toSharedPtr<Real>(typed));
    }

    /**
     * @brief Creates an ultimate-attribute-opening computer for the runtime attribute type.
     *
     * @param weighted Weighted tree used by the operation.
     * @param attr Attribute requested by the operation.
     * @return Ultimate-attribute-opening computer variant for the array element type.
     */
    static std::variant<UltimateAttributeOpening<std::uint8_t, float>, UltimateAttributeOpening<std::uint8_t, double>>
    makeUao(WeightedMorphologicalTree<std::uint8_t>& weighted, py::array attr) {
        if (pybind_utils::parseFloatingArrayDType(attr, "attr") == pybind_utils::FloatingDType::Float64) {
            return makeUao<double>(weighted, std::move(attr));
        }
        return makeUao<float>(weighted, std::move(attr));
    }

  public:
    /**
     * @brief Constructs `UltimateAttributeOpeningPybind` from the supplied inputs.
     *
     * @param weighted Weighted tree used by the operation.
     * @param attr Attribute requested by the operation.
     */
    UltimateAttributeOpeningPybind(std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weighted, py::array attr)
        : weightedOwner_(std::move(weighted)), uao_(makeUao(*weightedOwner_, std::move(attr))) {}

    /**
     * @brief Executes the ultimate attribute opening with the supplied criterion.
     *
     * @param maxCriterion Largest attribute criterion evaluated by the filter.
     */
    void execute(double maxCriterion) {
        std::visit(
            [maxCriterion](auto& uao) {
                using Uao = std::decay_t<decltype(uao)>;
                using Real = typename Uao::attribute_value_type;
                uao.execute(static_cast<Real>(maxCriterion));
            },
            uao_);
    }

    /**
     * @brief Executes with mser.
     *
     * @param maxCriterion Largest attribute criterion evaluated by the filter.
     * @param deltaMSER Altitude delta used by the MSER stability criterion.
     */
    void executeWithMSER(double maxCriterion, int deltaMSER) {
        std::visit(
            [maxCriterion, deltaMSER](auto& uao) {
                using Uao = std::decay_t<decltype(uao)>;
                using Real = typename Uao::attribute_value_type;
                uao.executeWithMSER(static_cast<Real>(maxCriterion), deltaMSER);
            },
            uao_);
    }

    /**
     * @brief Executes with depth stability.
     *
     * @param maxCriterion Largest attribute criterion evaluated by the filter.
     * @param depthDelta Topological-depth radius of the stability window.
     */
    void executeWithDepthStability(double maxCriterion, int depthDelta) {
        std::visit(
            [maxCriterion, depthDelta](auto& uao) {
                using Uao = std::decay_t<decltype(uao)>;
                using Real = typename Uao::attribute_value_type;
                uao.executeWithDepthStability(static_cast<Real>(maxCriterion), depthDelta);
            },
            uao_);
    }

    /**
     * @brief Returns max contrast image.
     *
     * @return Max contrast image.
     */
    py::array_t<uint8_t> getMaxContrastImage() {
        return std::visit([](auto& uao) { return pybind_utils::toNumpy(uao.getMaxContrastImage()); }, uao_);
    }

    /**
     * @brief Returns associated image.
     *
     * @return Associated image.
     */
    py::array_t<int32_t> getAssociatedImage() {
        return std::visit([](auto& uao) { return pybind_utils::toNumpy(uao.getAssociatedImage()); }, uao_);
    }
    /**
     * @brief Returns associated color image.
     *
     * @return Associated color image.
     */
    py::array_t<uint8_t> getAssociatedColorImage() {
        return std::visit([](auto& uao) { return pybind_utils::toNumpy(uao.getAssociatedColorImage()); }, uao_);
    }
};

} // namespace mmcfilters
