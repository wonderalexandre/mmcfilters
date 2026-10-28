#pragma once

#include "../mmcfilters/filters/AttributeFilters.hpp"

#include "ExtinctionValuesPybind.hpp"
#include "PybindConversions.hpp"
#include "PythonValuedMorphologicalTree.hpp"

#include <concepts>
#include <cstddef>
#include <memory>
#include <span>
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
    std::size_t topologyMutationVersion_; ///< Topology version shared by every filtering rule.
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

    /** @brief Rejects every rule after a topology edit. @param context Public operation name. */
    void requireStableTree(const char* context) const { topology().requireMutationVersion(topologyMutationVersion_, context); }

    /**
     * @brief Applies soft subtractive reconstruction for one score dtype.
     * @tparam Real Floating-point score and output type.
     * @param scores Dense one-dimensional NumPy preservation scores.
     * @return Reconstructed NumPy image with dtype Real.
     */
    template <std::floating_point Real> py::array filteringBySoftSubtractiveRuleTyped(py::array scores) {
        auto typed = pybind_utils::requireNodeAttributeArray<Real>(std::move(scores), topology(), "node_preservation_scores");
        const std::span<const Real> values(static_cast<const Real*>(typed.request().ptr), static_cast<std::size_t>(topology().numInternalNodeSlots()));
        return valuedTreeOwner_->visit([&](const auto& native) -> py::array {
            using Tree = typename std::remove_cvref_t<decltype(native)>::element_type;
            using Altitude = typename Tree::AltitudeType;
            return pybind_utils::toNumpy(SoftSubtractiveAttributeFilter<Altitude, Real>(*native).applySoftSubtractiveAttributeFilter(values));
        });
    }

    /**
     * @brief Applies by pruning min typed.
     *
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @return Image or array produced by the operation.
     */
    template <std::floating_point Real> py::array filteringByPruningMinTyped(py::array attribute, Real threshold) {
        auto typed = pybind_utils::requireNodeAttributeArray<Real>(std::move(attribute), topology(), "attribute");
        const auto* values = static_cast<const Real*>(typed.request().ptr);
        return withFilters([&](auto& filters) -> py::array { return pybind_utils::toNumpy(filters.filteringByPruningMin(values, threshold)); });
    }

    /**
     * @brief Applies by pruning max typed.
     *
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @return Image or array produced by the operation.
     */
    template <std::floating_point Real> py::array filteringByPruningMaxTyped(py::array attribute, Real threshold) {
        auto typed = pybind_utils::requireNodeAttributeArray<Real>(std::move(attribute), topology(), "attribute");
        const auto* values = static_cast<const Real*>(typed.request().ptr);
        return withFilters([&](auto& filters) -> py::array { return pybind_utils::toNumpy(filters.filteringByPruningMax(values, threshold)); });
    }

    /**
     * @brief Applies by viterbi rule typed.
     *
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @return Image or array produced by the operation.
     */
    template <std::floating_point Real> py::array filteringByViterbiRuleTyped(py::array attribute, Real threshold) {
        auto typed = pybind_utils::requireNodeAttributeArray<Real>(std::move(attribute), topology(), "attribute");
        const py::buffer_info buffer = typed.request();
        return withFilters([&](auto& filters) -> py::array {
            return pybind_utils::toNumpy(filters.filteringByViterbiRule(static_cast<const Real*>(buffer.ptr), threshold));
        });
    }

    /**
     * @brief Applies by extinction value typed.
     *
     * @param attribute Attribute requested by the operation.
     * @param selection Policy used to select extrema from the ranked candidates.
     * @return Image or array produced by the operation.
     */
    template <std::floating_point Real> py::array filteringByExtinctionValueTyped(py::array attribute, const ExtinctionSelectionPolicyPybind& selection) {
        auto typed = pybind_utils::requireNodeAttributeArray<Real>(std::move(attribute), topology(), "attribute");
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
     * @param attribute Attribute requested by the operation.
     * @param selection Policy used to select extrema from the ranked candidates.
     * @param scorePolicy Policy used to score and rank candidate extrema.
     * @return Map by extinction value typed.
     */
    template <std::floating_point Real>
    py::array contourMapByExtinctionValueTyped(py::array attribute, const ExtinctionSelectionPolicyPybind& selection, ExtinctionContourScorePolicy scorePolicy) {
        auto typed = pybind_utils::requireNodeAttributeArray<Real>(std::move(attribute), topology(), "attribute");
        const auto* values = static_cast<const Real*>(typed.request().ptr);
        return valuedTreeOwner_->visit([&](const auto& concreteTree) -> py::array {
            using Tree = typename std::remove_cvref_t<decltype(concreteTree)>::element_type;
            using Altitude = typename Tree::AltitudeType;
            ExtinctionValues<Altitude, Real> extinction(*concreteTree, values);
            return pybind_utils::toNumpy(extinction.contourMap(selection.toNative<Real>(), scorePolicy));
        });
    }

  public:
    /**
     * @brief Constructs `AttributeFiltersPybind` from the supplied inputs.
     *
     * @param valuedTree Valued tree.
     */
    explicit AttributeFiltersPybind(std::shared_ptr<pybindings::PythonValuedMorphologicalTree> valuedTree)
        : valuedTreeOwner_(requireOwner(std::move(valuedTree))), topologyMutationVersion_(topology().getMutationVersion()),
          filters_(makeFilters(*valuedTreeOwner_)) {}

    /**
     * @brief Applies direct reconstruction with an implicit root reference.
     * @param nodePreservationMask Boolean NumPy preservation decisions.
     * @return Reconstructed image in the tree altitude dtype.
     */
    py::array filteringByDirectRule(py::array nodePreservationMask) {
        requireStableTree("AttributeFilters.filtering_by_direct_rule");
        const auto decisions = pybind_utils::requireNodePreservationArray(nodePreservationMask, topology());
        return valuedTreeOwner_->visit([&](const auto& native) -> py::array {
            using Tree = typename std::remove_cvref_t<decltype(native)>::element_type;
            using Altitude = typename Tree::AltitudeType;
            return pybind_utils::toNumpy(DirectAttributeFilter<Altitude>(*native).applyDirectAttributeFilter(decisions));
        });
    }

    /**
     * @brief Applies subtractive reconstruction from Boolean decisions.
     * @param nodePreservationMask Boolean NumPy preservation decisions.
     * @return Reconstructed image in the signed altitude-difference dtype.
     */
    py::array filteringBySubtractiveRule(py::array nodePreservationMask) {
        requireStableTree("AttributeFilters.filtering_by_subtractive_rule");
        const auto decisions = pybind_utils::requireNodePreservationArray(nodePreservationMask, topology());
        return valuedTreeOwner_->visit([&](const auto& native) -> py::array {
            using Tree = typename std::remove_cvref_t<decltype(native)>::element_type;
            using OutputValue = AltitudeDifference<typename Tree::AltitudeType>;
            const auto contributions = detail::attribute_filtering::computeGatedNodeContributions(native->asView(), decisions);
            auto pixels = detail::tree_altitude::reconstructNodeContributionValuesEstablished(native->topology(), std::span<const OutputValue>(contributions));
            return pybind_utils::toNumpyOwned2D(std::move(pixels), topology().numRows(), topology().numColumns());
        });
    }

    /**
     * @brief Applies soft subtractive reconstruction from scores in [0, 1].
     * @param nodePreservationScores Dense float32 or float64 preservation scores.
     * @return Reconstructed image with the input floating-point dtype.
     */
    py::array filteringBySoftSubtractiveRule(py::array nodePreservationScores) {
        requireStableTree("AttributeFilters.filtering_by_soft_subtractive_rule");
        if (pybind_utils::parseFloatingArrayDType(nodePreservationScores, "node_preservation_scores") == pybind_utils::FloatingDType::Float64) {
            return filteringBySoftSubtractiveRuleTyped<double>(std::move(nodePreservationScores));
        }
        return filteringBySoftSubtractiveRuleTyped<float>(std::move(nodePreservationScores));
    }

    /**
     * @brief Applies by pruning min.
     *
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @return Image or array produced by the operation.
     */
    py::array filteringByPruningMin(py::array attribute, double threshold) {
        requireStableTree("AttributeFilters.filtering_by_pruning_min_rule");
        if (pybind_utils::parseFloatingArrayDType(attribute, "attribute") == pybind_utils::FloatingDType::Float64) {
            return filteringByPruningMinTyped<double>(std::move(attribute), threshold);
        }
        return filteringByPruningMinTyped<float>(std::move(attribute), static_cast<float>(threshold));
    }

    /**
     * @brief Applies by pruning max.
     *
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @return Image or array produced by the operation.
     */
    py::array filteringByPruningMax(py::array attribute, double threshold) {
        requireStableTree("AttributeFilters.filtering_by_pruning_max_rule");
        if (pybind_utils::parseFloatingArrayDType(attribute, "attribute") == pybind_utils::FloatingDType::Float64) {
            return filteringByPruningMaxTyped<double>(std::move(attribute), threshold);
        }
        return filteringByPruningMaxTyped<float>(std::move(attribute), static_cast<float>(threshold));
    }

    /**
     * @brief Applies by viterbi rule.
     *
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @return Image or array produced by the operation.
     */
    py::array filteringByViterbiRule(py::array attribute, double threshold) {
        requireStableTree("AttributeFilters.filtering_by_viterbi_rule");
        if (pybind_utils::parseFloatingArrayDType(attribute, "attribute") == pybind_utils::FloatingDType::Float64) {
            return filteringByViterbiRuleTyped<double>(std::move(attribute), threshold);
        }
        return filteringByViterbiRuleTyped<float>(std::move(attribute), static_cast<float>(threshold));
    }

    /**
     * @brief Applies by pruning min.
     *
     * @param nodePreservationMask Per-node decisions where true preserves a node.
     * @return Image or array produced by the operation.
     */
    py::array filteringByPruningMin(py::array nodePreservationMask) {
        requireStableTree("AttributeFilters.filtering_by_pruning_min_rule");
        const auto decisions = pybind_utils::requireNodePreservationArray(nodePreservationMask, topology());
        return withFilters([&](auto& filters) -> py::array { return pybind_utils::toNumpy(filters.filteringByPruningMin(decisions)); });
    }

    /**
     * @brief Applies by pruning max.
     *
     * @param nodePreservationMask Per-node decisions where true preserves a node.
     * @return Image or array produced by the operation.
     */
    py::array filteringByPruningMax(py::array nodePreservationMask) {
        requireStableTree("AttributeFilters.filtering_by_pruning_max_rule");
        const auto decisions = pybind_utils::requireNodePreservationArray(nodePreservationMask, topology());
        return withFilters([&](auto& filters) -> py::array { return pybind_utils::toNumpy(filters.filteringByPruningMax(decisions)); });
    }

    /**
     * @brief Applies by extinction value.
     *
     * @param attribute Attribute requested by the operation.
     * @param selection Policy used to select extrema from the ranked candidates.
     * @return Image or array produced by the operation.
     */
    py::array filteringByExtinctionValue(py::array attribute, const ExtinctionSelectionPolicyPybind& selection) {
        requireStableTree("AttributeFilters.filtering_by_extinction_rule");
        if (pybind_utils::parseFloatingArrayDType(attribute, "attribute") == pybind_utils::FloatingDType::Float64) {
            return filteringByExtinctionValueTyped<double>(std::move(attribute), selection);
        }
        return filteringByExtinctionValueTyped<float>(std::move(attribute), selection);
    }

    /**
     * @brief Returns map by extinction value.
     *
     * @param attribute Attribute requested by the operation.
     * @param selection Policy used to select extrema from the ranked candidates.
     * @param scorePolicy Policy used to score and rank candidate extrema.
     * @return Map by extinction value.
     */
    py::array contourMapByExtinctionValue(py::array attribute, const ExtinctionSelectionPolicyPybind& selection, ExtinctionContourScorePolicy scorePolicy) {
        requireStableTree("AttributeFilters.contour_map_by_extinction");
        if (pybind_utils::parseFloatingArrayDType(attribute, "attribute") == pybind_utils::FloatingDType::Float64) {
            return contourMapByExtinctionValueTyped<double>(std::move(attribute), selection, scorePolicy);
        }
        return contourMapByExtinctionValueTyped<float>(std::move(attribute), selection, scorePolicy);
    }
};

} // namespace mmcfilters
