#pragma once

#include "../mmcfilters/filters/AttributeFilters.hpp"

#include "ExtinctionValuesPybind.hpp"
#include "PybindConversions.hpp"
#include "PythonValuedMorphologicalTree.hpp"

#include <concepts>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace mmcfilters {

/**
 * @brief Pybind11 wrapper exposing attribute filtering operators to Python.
 */
class AttributeFiltersPybind {
    /** @brief References the valued-tree owner used by the component. */
    std::shared_ptr<pybindings::PythonValuedMorphologicalTree> valuedTreeOwner_;
    /** @brief Runtime storage for the supported altitude specializations. */
    using FilterStorage =
        std::variant<std::unique_ptr<AttributeFilters<std::uint8_t>>, std::unique_ptr<AttributeFilters<ToSGrayLevel>>>;
    FilterStorage filters_; ///< Concrete filter implementation.

    /** @brief Validates and retains a Python tree owner. @param valuedTree Candidate owner. @return Valid non-null owner. */
    static std::shared_ptr<pybindings::PythonValuedMorphologicalTree>
    requireOwner(std::shared_ptr<pybindings::PythonValuedMorphologicalTree> valuedTree) {
        if (!valuedTree) {
            throw std::invalid_argument("AttributeFilters requires a non-null ValuedMorphologicalTree.");
        }
        return valuedTree;
    }

    /** @brief Creates the concrete native filter. @param valuedTree Runtime-valued tree. @return Filter variant. */
    static FilterStorage makeFilters(pybindings::PythonValuedMorphologicalTree& valuedTree) {
        return valuedTree.visit([](auto& concreteTree) -> FilterStorage {
            using Tree = typename std::remove_cvref_t<decltype(concreteTree)>::element_type;
            using Altitude = typename Tree::AltitudeType;
            return std::make_unique<AttributeFilters<Altitude>>(*concreteTree);
        });
    }

    /** @brief Invokes a callable on the active filter. @tparam Function Callable type. @param function Callable to invoke. @return Callable result. */
    template <class Function> decltype(auto) withFilters(Function&& function) {
        return std::visit([&](auto& filters) -> decltype(auto) { return std::forward<Function>(function)(*filters); }, filters_);
    }

    /**
     * @brief Returns the borrowed tree topology.
     *
     * @return The borrowed tree topology.
     */
    const MorphologicalTree& topology() const noexcept { return valuedTreeOwner_->topology(); }


    /**
     * @brief Applies by pruning min typed.
     *
     * @param attr Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @return Image or array produced by the operation.
     */
    template <std::floating_point Real> py::array filteringByPruningMinTyped(py::array attr, Real threshold) {
        auto typed = pybind_utils::requireNodeAttributeArray<Real>(std::move(attr), topology());
        const auto* values = static_cast<const Real*>(typed.request().ptr);
        return withFilters([&](auto& filters) -> py::array { return pybind_utils::toNumpy(filters.filteringByPruningMin(values, threshold)); });
    }

    /**
     * @brief Applies by pruning max typed.
     *
     * @param attr Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @return Image or array produced by the operation.
     */
    template <std::floating_point Real> py::array filteringByPruningMaxTyped(py::array attr, Real threshold) {
        auto typed = pybind_utils::requireNodeAttributeArray<Real>(std::move(attr), topology());
        const auto* values = static_cast<const Real*>(typed.request().ptr);
        return withFilters([&](auto& filters) -> py::array { return pybind_utils::toNumpy(filters.filteringByPruningMax(values, threshold)); });
    }

    /**
     * @brief Applies by viterbi rule typed.
     *
     * @param attr Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @return Image or array produced by the operation.
     */
    template <std::floating_point Real> py::array filteringByViterbiRuleTyped(py::array attr, Real threshold) {
        auto typed = pybind_utils::requireNodeAttributeArray<Real>(std::move(attr), topology());
        const py::buffer_info buffer = typed.request();
        return withFilters([&](auto& filters) -> py::array {
            return pybind_utils::toNumpy(filters.filteringByViterbiRule(static_cast<const Real*>(buffer.ptr), threshold));
        });
    }

    /**
     * @brief Applies by extinction value typed.
     *
     * @param attr Attribute requested by the operation.
     * @param selection Policy used to select extrema from the ranked candidates.
     * @return Image or array produced by the operation.
     */
    template <std::floating_point Real> py::array filteringByExtinctionValueTyped(py::array attr, const ExtinctionSelectionPolicyPybind& selection) {
        auto typed = pybind_utils::requireNodeAttributeArray<Real>(std::move(attr), topology());
        const auto* values = static_cast<const Real*>(typed.request().ptr);
        return valuedTreeOwner_->visit([&](const auto& concreteTree) -> py::array {
            using Tree = typename std::remove_cvref_t<decltype(concreteTree)>::element_type;
            using Altitude = typename Tree::AltitudeType;
            ExtinctionValues<Altitude, Real> extinction(*concreteTree, values);
            return pybind_utils::toNumpy(extinction.filtering(selection.toNative<Real>()));
        });
    }

    /**
     * @brief Returns map by extinction value typed.
     *
     * @param attr Attribute requested by the operation.
     * @param selection Policy used to select extrema from the ranked candidates.
     * @param scorePolicy Policy used to score and rank candidate extrema.
     * @return Map by extinction value typed.
     */
    template <std::floating_point Real>
    py::array contourMapByExtinctionValueTyped(py::array attr, const ExtinctionSelectionPolicyPybind& selection, ExtinctionContourScorePolicy scorePolicy) {
        auto typed = pybind_utils::requireNodeAttributeArray<Real>(std::move(attr), topology());
        const auto* values = static_cast<const Real*>(typed.request().ptr);
        return valuedTreeOwner_->visit([&](const auto& concreteTree) -> py::array {
            using Tree = typename std::remove_cvref_t<decltype(concreteTree)>::element_type;
            using Altitude = typename Tree::AltitudeType;
            ExtinctionValues<Altitude, Real> extinction(*concreteTree, values);
            return pybind_utils::toNumpy(extinction.contourMap(selection.toNative<Real>(), scorePolicy));
        });
    }

    /**
     * @brief Validates a node-preservation mask against the active tree.
     *
     * @param nodePreservationMask Per-node decisions where true preserves a node.
     * @param tree Tree topology.
     * @param argumentName Argument name included in validation error messages.
     */
    static void requireNodePreservationMask(const NodePreservationMask& nodePreservationMask, const MorphologicalTree& tree,
                                            std::string_view argumentName = "node_preservation_mask") {
        pybind_utils::requireVectorSize(nodePreservationMask.decisions(), static_cast<std::size_t>(tree.numInternalNodeSlots()), argumentName);
    }

  public:
    /** @brief Returns the retained valued-tree owner. @return Shared owner of the filtered tree. */
    [[nodiscard]] std::shared_ptr<pybindings::PythonValuedMorphologicalTree> treeOwner() const noexcept { return valuedTreeOwner_; }

    /**
     * @brief Constructs `AttributeFiltersPybind` from the supplied inputs.
     *
     * @param valuedTree Valued tree.
     */
    explicit AttributeFiltersPybind(std::shared_ptr<pybindings::PythonValuedMorphologicalTree> valuedTree)
        : valuedTreeOwner_(requireOwner(std::move(valuedTree))), filters_(makeFilters(*valuedTreeOwner_)) {}

    /**
     * @brief Applies by pruning min.
     *
     * @param attr Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @return Image or array produced by the operation.
     */
    py::array filteringByPruningMin(py::array attr, double threshold) {
        if (pybind_utils::parseFloatingArrayDType(attr, "attr") == pybind_utils::FloatingDType::Float64) {
            return filteringByPruningMinTyped<double>(std::move(attr), threshold);
        }
        return filteringByPruningMinTyped<float>(std::move(attr), static_cast<float>(threshold));
    }

    /**
     * @brief Applies by pruning max.
     *
     * @param attr Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @return Image or array produced by the operation.
     */
    py::array filteringByPruningMax(py::array attr, double threshold) {
        if (pybind_utils::parseFloatingArrayDType(attr, "attr") == pybind_utils::FloatingDType::Float64) {
            return filteringByPruningMaxTyped<double>(std::move(attr), threshold);
        }
        return filteringByPruningMaxTyped<float>(std::move(attr), static_cast<float>(threshold));
    }

    /**
     * @brief Applies by viterbi rule.
     *
     * @param attr Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @return Image or array produced by the operation.
     */
    py::array filteringByViterbiRule(py::array attr, double threshold) {
        if (pybind_utils::parseFloatingArrayDType(attr, "attr") == pybind_utils::FloatingDType::Float64) {
            return filteringByViterbiRuleTyped<double>(std::move(attr), threshold);
        }
        return filteringByViterbiRuleTyped<float>(std::move(attr), static_cast<float>(threshold));
    }

    /**
     * @brief Applies by pruning min.
     *
     * @param nodePreservationMask Per-node decisions where true preserves a node.
     * @return Image or array produced by the operation.
     */
    py::array filteringByPruningMin(const NodePreservationMask& nodePreservationMask) {
        requireNodePreservationMask(nodePreservationMask, topology());
        return withFilters([&](auto& filters) -> py::array { return pybind_utils::toNumpy(filters.filteringByPruningMin(nodePreservationMask)); });
    }

    /**
     * @brief Applies by pruning max.
     *
     * @param nodePreservationMask Per-node decisions where true preserves a node.
     * @return Image or array produced by the operation.
     */
    py::array filteringByPruningMax(const NodePreservationMask& nodePreservationMask) {
        requireNodePreservationMask(nodePreservationMask, topology());
        return withFilters([&](auto& filters) -> py::array { return pybind_utils::toNumpy(filters.filteringByPruningMax(nodePreservationMask)); });
    }

    /**
     * @brief Applies by extinction value.
     *
     * @param attr Attribute requested by the operation.
     * @param selection Policy used to select extrema from the ranked candidates.
     * @return Image or array produced by the operation.
     */
    py::array filteringByExtinctionValue(py::array attr, const ExtinctionSelectionPolicyPybind& selection) {
        if (pybind_utils::parseFloatingArrayDType(attr, "attr") == pybind_utils::FloatingDType::Float64) {
            return filteringByExtinctionValueTyped<double>(std::move(attr), selection);
        }
        return filteringByExtinctionValueTyped<float>(std::move(attr), selection);
    }

    /**
     * @brief Returns map by extinction value.
     *
     * @param attr Attribute requested by the operation.
     * @param selection Policy used to select extrema from the ranked candidates.
     * @param scorePolicy Policy used to score and rank candidate extrema.
     * @return Map by extinction value.
     */
    py::array contourMapByExtinctionValue(py::array attr, const ExtinctionSelectionPolicyPybind& selection, ExtinctionContourScorePolicy scorePolicy) {
        if (pybind_utils::parseFloatingArrayDType(attr, "attr") == pybind_utils::FloatingDType::Float64) {
            return contourMapByExtinctionValueTyped<double>(std::move(attr), selection, scorePolicy);
        }
        return contourMapByExtinctionValueTyped<float>(std::move(attr), selection, scorePolicy);
    }
};

} // namespace mmcfilters
