#pragma once

#include "../mmcfilters/filters/AttributeFilters.hpp"

#include "ExtinctionValuesPybind.hpp"
#include "PybindConversions.hpp"

#include <concepts>
#include <string>
#include <string_view>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace mmcfilters {

/**
 * @brief Pybind11 wrapper exposing attribute filtering operators to Python.
 */
class AttributeFiltersPybind : public AttributeFilters<std::uint8_t> {
    /** @brief References the weighted owner used by the component. */
    std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weightedOwner_;

    /**
     * @brief Returns the borrowed tree topology.
     *
     * @return The borrowed tree topology.
     */
    const MorphologicalTree& topology() const noexcept { return this->tree; }

    /**
     * @brief Applies by pruning min typed.
     *
     * @param attr Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @return Image or array produced by the operation.
     */
    template <std::floating_point Real> py::array_t<uint8_t> filteringByPruningMinTyped(py::array attr, Real threshold) {
        auto typed = pybind_utils::requireNodeAttributeArray<Real>(std::move(attr), topology());
        return pybind_utils::toNumpy(AttributeFilters<std::uint8_t>::filteringByPruningMin(pybind_utils::toSharedPtr<Real>(typed), threshold));
    }

    /**
     * @brief Applies by pruning max typed.
     *
     * @param attr Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @return Image or array produced by the operation.
     */
    template <std::floating_point Real> py::array_t<uint8_t> filteringByPruningMaxTyped(py::array attr, Real threshold) {
        auto typed = pybind_utils::requireNodeAttributeArray<Real>(std::move(attr), topology());
        return pybind_utils::toNumpy(AttributeFilters<std::uint8_t>::filteringByPruningMax(pybind_utils::toSharedPtr<Real>(typed), threshold));
    }

    /**
     * @brief Applies by viterbi rule typed.
     *
     * @param attr Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @return Image or array produced by the operation.
     */
    template <std::floating_point Real> py::array_t<uint8_t> filteringByViterbiRuleTyped(py::array attr, Real threshold) {
        auto typed = pybind_utils::requireNodeAttributeArray<Real>(std::move(attr), topology());
        const py::buffer_info buffer = typed.request();
        return pybind_utils::toNumpy(AttributeFilters<std::uint8_t>::filteringByViterbiRule(static_cast<const Real*>(buffer.ptr), threshold));
    }

    /**
     * @brief Applies by extinction value typed.
     *
     * @param attr Attribute requested by the operation.
     * @param selection Policy used to select extrema from the ranked candidates.
     * @return Image or array produced by the operation.
     */
    template <std::floating_point Real> py::array_t<uint8_t> filteringByExtinctionValueTyped(py::array attr, const ExtinctionSelectionPolicyPybind& selection) {
        auto typed = pybind_utils::requireNodeAttributeArray<Real>(std::move(attr), topology());
        ExtinctionValues<std::uint8_t, Real> ev(*this->weightedOwner_, pybind_utils::toSharedPtr<Real>(typed));
        return pybind_utils::toNumpy(ev.filtering(selection.toNative<Real>()));
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
        ExtinctionValues<std::uint8_t, Real> ev(*this->weightedOwner_, pybind_utils::toSharedPtr<Real>(typed));
        return pybind_utils::toNumpy(ev.contourMap(selection.toNative<Real>(), scorePolicy));
    }

    /**
     * @brief Validates node criterion.
     *
     * @param criterion Per-node selection criterion.
     * @param tree Tree topology used by the operation.
     * @param argumentName Argument name included in validation error messages.
     */
    static void requireNodeCriterion(const std::vector<bool>& criterion, const MorphologicalTree& tree, std::string_view argumentName = "criterion") {
        pybind_utils::requireVectorSize(criterion, static_cast<std::size_t>(tree.getNumInternalNodeSlots()), argumentName);
    }

    /**
     * @brief Validates node scores.
     *
     * @param prob Per-node probability values used by the decision rule.
     * @param tree Tree topology used by the operation.
     * @param argumentName Argument name included in validation error messages.
     */
    static void requireNodeScores(const std::vector<float>& prob, const MorphologicalTree& tree, std::string_view argumentName = "prob") {
        pybind_utils::requireVectorSize(prob, static_cast<std::size_t>(tree.getNumInternalNodeSlots()), argumentName);
    }

  public:
    using AttributeFilters<std::uint8_t>::AttributeFilters;

    /**
     * @brief Constructs `AttributeFiltersPybind` from the supplied inputs.
     *
     * @param weighted Weighted tree used by the operation.
     */
    explicit AttributeFiltersPybind(std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weighted)
        : AttributeFilters<std::uint8_t>(*weighted), weightedOwner_(std::move(weighted)) {}

    /**
     * @brief Applies by pruning min.
     *
     * @param attr Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @return Image or array produced by the operation.
     */
    py::array_t<uint8_t> filteringByPruningMin(py::array attr, double threshold) {
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
    py::array_t<uint8_t> filteringByPruningMax(py::array attr, double threshold) {
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
    py::array_t<uint8_t> filteringByViterbiRule(py::array attr, double threshold) {
        if (pybind_utils::parseFloatingArrayDType(attr, "attr") == pybind_utils::FloatingDType::Float64) {
            return filteringByViterbiRuleTyped<double>(std::move(attr), threshold);
        }
        return filteringByViterbiRuleTyped<float>(std::move(attr), static_cast<float>(threshold));
    }

    /**
     * @brief Applies by pruning min.
     *
     * @param criterion Per-node selection criterion.
     * @return Image or array produced by the operation.
     */
    py::array_t<uint8_t> filteringByPruningMin(std::vector<bool>& criterion) {
        requireNodeCriterion(criterion, topology());
        return pybind_utils::toNumpy(AttributeFilters<std::uint8_t>::filteringByPruningMin(criterion));
    }

    /**
     * @brief Applies by direct rule.
     *
     * @param criterion Per-node selection criterion.
     * @return Image or array produced by the operation.
     */
    py::array_t<uint8_t> filteringByDirectRule(std::vector<bool>& criterion) {
        requireNodeCriterion(criterion, topology());
        return pybind_utils::toNumpy(AttributeFilters<std::uint8_t>::filteringByDirectRule(criterion));
    }

    /**
     * @brief Applies by pruning max.
     *
     * @param criterion Per-node selection criterion.
     * @return Image or array produced by the operation.
     */
    py::array_t<uint8_t> filteringByPruningMax(std::vector<bool>& criterion) {
        requireNodeCriterion(criterion, topology());
        return pybind_utils::toNumpy(AttributeFilters<std::uint8_t>::filteringByPruningMax(criterion));
    }

    /**
     * @brief Returns adaptive criterion.
     *
     * @param criterion Per-node selection criterion.
     * @param delta Delta offset used by the operation.
     * @return Adaptive criterion.
     */
    std::vector<bool> getAdaptiveCriterion(std::vector<bool>& criterion, int delta) {
        requireNodeCriterion(criterion, topology());
        return AttributeFilters<std::uint8_t>::getAdaptiveCriterion(criterion, delta);
    }

    /**
     * @brief Returns adaptive criterion by depth.
     *
     * @param criterion Per-node selection criterion.
     * @param depthDelta Topological-depth radius of the stability window.
     * @return Adaptive criterion by depth.
     */
    std::vector<bool> getAdaptiveCriterionByDepth(std::vector<bool>& criterion, int depthDelta) {
        requireNodeCriterion(criterion, topology());
        return AttributeFilters<std::uint8_t>::getAdaptiveCriterionByDepth(criterion, depthDelta);
    }

    /**
     * @brief Applies by subtractive rule.
     *
     * @param criterion Per-node selection criterion.
     * @return Image or array produced by the operation.
     */
    py::array_t<uint8_t> filteringBySubtractiveRule(std::vector<bool>& criterion) {
        requireNodeCriterion(criterion, topology());
        return pybind_utils::toNumpy(AttributeFilters<std::uint8_t>::filteringBySubtractiveRule(criterion));
    }

    /**
     * @brief Applies by subtractive score rule.
     *
     * @param prob Per-node probability values used by the decision rule.
     * @return Image or array produced by the operation.
     */
    py::array_t<float> filteringBySubtractiveScoreRule(std::vector<float>& prob) {
        requireNodeScores(prob, topology());
        return pybind_utils::toNumpy(AttributeFilters<std::uint8_t>::filteringBySubtractiveScoreRule(prob));
    }

    /**
     * @brief Applies by extinction value.
     *
     * @param attr Attribute requested by the operation.
     * @param selection Policy used to select extrema from the ranked candidates.
     * @return Image or array produced by the operation.
     */
    py::array_t<uint8_t> filteringByExtinctionValue(py::array attr, const ExtinctionSelectionPolicyPybind& selection) {
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
