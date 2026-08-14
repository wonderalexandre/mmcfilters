#pragma once
#include <concepts>
#include <variant>
#include "../mmcfilters/filters/UltimateAttributeOpening.hpp"
#include "PybindConversions.hpp"
#include "PythonValuedMorphologicalTree.hpp"
#include <pybind11/pybind11.h>
namespace mmcfilters {

namespace py = pybind11;

/**
 * @brief Pybind11 wrapper for Ultimate Attribute Opening.
 */
class UltimateAttributeOpeningPybind {
    /** @brief References the valued-tree owner used by the component. */
    std::shared_ptr<pybindings::PythonValuedMorphologicalTree> valuedTreeOwner_;
    /** @brief Keeps the NumPy attribute buffer alive for the native borrowed pointer. */
    py::array attributeOwner_;
    /** @brief Defines all native UAO evaluators supported at runtime. */
    using UaoStorage =
        std::variant<UltimateAttributeOpening<std::uint8_t, float>, UltimateAttributeOpening<std::uint8_t, double>,
                     UltimateAttributeOpening<ToSGrayLevel, float>, UltimateAttributeOpening<ToSGrayLevel, double>>;
    UaoStorage uao_; ///< Concrete ultimate-attribute-opening implementation.

    /**
     * @brief Creates an ultimate-attribute-opening computer for the runtime attribute type.
     *
     * @param valuedTree Valued tree.
     * @param attr Attribute requested by the operation.
     * @return Ultimate-attribute-opening computer variant for the array element type.
     */
    template <std::floating_point Real, AltitudeValue Altitude>
    static UltimateAttributeOpening<Altitude, Real> makeUao(ValuedMorphologicalTree<Altitude>& valuedTree, py::array attr) {
        auto typed = pybind_utils::requireNodeAttributeArray<Real>(std::move(attr), valuedTree.topology());
        return UltimateAttributeOpening<Altitude, Real>(valuedTree, static_cast<const Real*>(typed.request().ptr));
    }

    /**
     * @brief Creates an ultimate-attribute-opening computer for the runtime attribute type.
     *
     * @param valuedTree Valued tree.
     * @param attr Attribute requested by the operation.
     * @return Ultimate-attribute-opening computer variant for the array element type.
     */
    static UaoStorage makeUao(pybindings::PythonValuedMorphologicalTree& valuedTree, py::array attr) {
        if (pybind_utils::parseFloatingArrayDType(attr, "attr") == pybind_utils::FloatingDType::Float64) {
            return valuedTree.visit([&](auto& concreteTree) -> UaoStorage { return makeUao<double>(*concreteTree, std::move(attr)); });
        }
        return valuedTree.visit([&](auto& concreteTree) -> UaoStorage { return makeUao<float>(*concreteTree, std::move(attr)); });
    }

  public:
    /**
     * @brief Constructs `UltimateAttributeOpeningPybind` from the supplied inputs.
     *
     * @param valuedTree Valued tree.
     * @param attr Attribute requested by the operation.
     */
    UltimateAttributeOpeningPybind(std::shared_ptr<pybindings::PythonValuedMorphologicalTree> valuedTree, py::array attr)
        : valuedTreeOwner_(std::move(valuedTree)), attributeOwner_(std::move(attr)), uao_(makeUao(*valuedTreeOwner_, attributeOwner_)) {}

    /**
     * @brief Executes the ultimate attribute opening up to an attribute threshold.
     *
     * @param maximumAttributeThreshold Largest increasing-attribute value evaluated by the filter.
     */
    void execute(double maximumAttributeThreshold) {
        std::visit(
            [maximumAttributeThreshold](auto& uao) {
                using Uao = std::decay_t<decltype(uao)>;
                using Real = typename Uao::attribute_value_type;
                uao.execute(static_cast<Real>(maximumAttributeThreshold));
            },
            uao_);
    }

    /**
     * @brief Executes UAO with an altitude-stability primitive-selection mask.
     *
     * @param maximumAttributeThreshold Largest increasing-attribute value evaluated by the filter.
     * @param altitudeWindowRadius Radius of the altitude-based stability window.
     */
    void executeWithMSER(double maximumAttributeThreshold, int altitudeWindowRadius) {
        std::visit(
            [maximumAttributeThreshold, altitudeWindowRadius](auto& uao) {
                using Uao = std::decay_t<decltype(uao)>;
                using Real = typename Uao::attribute_value_type;
                uao.executeWithMSER(static_cast<Real>(maximumAttributeThreshold), altitudeWindowRadius);
            },
            uao_);
    }

    /**
     * @brief Executes UAO with a depth-stability primitive-selection mask.
     *
     * @param maximumAttributeThreshold Largest increasing-attribute value evaluated by the filter.
     * @param depthWindowRadius Topological-depth radius of the stability window.
     */
    void executeWithDepthStability(double maximumAttributeThreshold, int depthWindowRadius) {
        std::visit(
            [maximumAttributeThreshold, depthWindowRadius](auto& uao) {
                using Uao = std::decay_t<decltype(uao)>;
                using Real = typename Uao::attribute_value_type;
                uao.executeWithDepthStability(static_cast<Real>(maximumAttributeThreshold), depthWindowRadius);
            },
            uao_);
    }

    /**
     * @brief Returns max contrast image.
     *
     * @return Max contrast image.
     */
    py::array getMaxContrastImage() {
        return std::visit([](auto& uao) -> py::array { return pybind_utils::toNumpy(uao.getMaxContrastImage()); }, uao_);
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
