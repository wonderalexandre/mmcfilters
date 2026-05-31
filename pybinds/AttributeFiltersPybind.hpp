#pragma once

#include "../mmcfilters/filters/AttributeFilters.hpp"
#include "../mmcfilters/filters/ExtinctionValues.hpp"

#include "MorphologicalTreePybind.hpp"
#include "PybindUtils.hpp"

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
    std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weightedOwner_;

    const MorphologicalTree& topology() const noexcept {
        return this->tree;
    }

    template <std::floating_point Real>
    py::array_t<uint8_t> filteringByPruningMinTyped(py::array attr, Real threshold) {
        auto typed = PybindUtils::requireNodeAttributeArray<Real>(std::move(attr), topology());
        return PybindUtils::toNumpy(AttributeFilters<std::uint8_t>::filteringByPruningMin(PybindUtils::toSharedPtr<Real>(typed), threshold));
    }

    template <std::floating_point Real>
    py::array_t<uint8_t> filteringByPruningMaxTyped(py::array attr, Real threshold) {
        auto typed = PybindUtils::requireNodeAttributeArray<Real>(std::move(attr), topology());
        return PybindUtils::toNumpy(AttributeFilters<std::uint8_t>::filteringByPruningMax(PybindUtils::toSharedPtr<Real>(typed), threshold));
    }

    template <std::floating_point Real>
    py::array_t<uint8_t> filteringByViterbiRuleTyped(py::array attr, Real threshold) {
        auto typed = PybindUtils::requireNodeAttributeArray<Real>(std::move(attr), topology());
        const py::buffer_info buffer = typed.request();
        return PybindUtils::toNumpy(AttributeFilters<std::uint8_t>::filteringByViterbiRule(
            static_cast<const Real*>(buffer.ptr),
            threshold));
    }

    template <std::floating_point Real>
    py::array_t<uint8_t> filteringByExtinctionValueTyped(py::array attr, int leafToKeep) {
        auto typed = PybindUtils::requireNodeAttributeArray<Real>(std::move(attr), topology());
        ExtinctionValues<std::uint8_t, Real> ev(*this->weightedOwner_, PybindUtils::toSharedPtr<Real>(typed));
        return PybindUtils::toNumpy(ev.filtering(leafToKeep));
    }

    template <std::floating_point Real>
    py::array saliencyMapByExtinctionValueTyped(py::array attr, int leafToKeep, bool unweighted) {
        auto typed = PybindUtils::requireNodeAttributeArray<Real>(std::move(attr), topology());
        ExtinctionValues<std::uint8_t, Real> ev(*this->weightedOwner_, PybindUtils::toSharedPtr<Real>(typed));
        return PybindUtils::toNumpy(ev.saliencyMap(leafToKeep, unweighted));
    }

    static void requireNodeCriterion(const std::vector<bool>& criterion, const MorphologicalTree& tree, std::string_view argumentName = "criterion") {
        PybindUtils::requireVectorSize(criterion, static_cast<std::size_t>(tree.getNumInternalNodeSlots()), argumentName);
    }

    static void requireNodeScores(const std::vector<float>& prob, const MorphologicalTree& tree, std::string_view argumentName = "prob") {
        PybindUtils::requireVectorSize(prob, static_cast<std::size_t>(tree.getNumInternalNodeSlots()), argumentName);
    }

    public:
    using AttributeFilters<std::uint8_t>::AttributeFilters;

    explicit AttributeFiltersPybind(std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> weighted)
        : AttributeFilters<std::uint8_t>(*weighted), weightedOwner_(std::move(weighted)) {}

    py::array_t<uint8_t> filteringByPruningMin(py::array attr, double threshold){
        if (PybindUtils::parseFloatingArrayDType(attr, "attr") == PybindUtils::FloatingDType::Float64) {
            return filteringByPruningMinTyped<double>(std::move(attr), threshold);
        }
        return filteringByPruningMinTyped<float>(std::move(attr), static_cast<float>(threshold));
    }

    py::array_t<uint8_t> filteringByPruningMax(py::array attr, double threshold){
        if (PybindUtils::parseFloatingArrayDType(attr, "attr") == PybindUtils::FloatingDType::Float64) {
            return filteringByPruningMaxTyped<double>(std::move(attr), threshold);
        }
        return filteringByPruningMaxTyped<float>(std::move(attr), static_cast<float>(threshold));
    }

    py::array_t<uint8_t> filteringByViterbiRule(py::array attr, double threshold){
        if (PybindUtils::parseFloatingArrayDType(attr, "attr") == PybindUtils::FloatingDType::Float64) {
            return filteringByViterbiRuleTyped<double>(std::move(attr), threshold);
        }
        return filteringByViterbiRuleTyped<float>(std::move(attr), static_cast<float>(threshold));
    }

    py::array_t<uint8_t> filteringByPruningMin(std::vector<bool>& criterion){
        requireNodeCriterion(criterion, topology());
        return PybindUtils::toNumpy(AttributeFilters<std::uint8_t>::filteringByPruningMin(criterion));
    }

    py::array_t<uint8_t> filteringByDirectRule(std::vector<bool>& criterion){
        requireNodeCriterion(criterion, topology());
        return PybindUtils::toNumpy(AttributeFilters<std::uint8_t>::filteringByDirectRule(criterion));
    }

    py::array_t<uint8_t> filteringByPruningMax(std::vector<bool>& criterion){
        requireNodeCriterion(criterion, topology());
        return PybindUtils::toNumpy(AttributeFilters<std::uint8_t>::filteringByPruningMax(criterion));

    }

    std::vector<bool> getAdaptiveCriterion(std::vector<bool>& criterion, int delta){
        requireNodeCriterion(criterion, topology());
        return AttributeFilters<std::uint8_t>::getAdaptiveCriterion(criterion, delta);
    }

    std::vector<bool> getAdaptiveCriterionByDepth(std::vector<bool>& criterion, int depthDelta){
        requireNodeCriterion(criterion, topology());
        return AttributeFilters<std::uint8_t>::getAdaptiveCriterionByDepth(criterion, depthDelta);
    }



    py::array_t<uint8_t> filteringBySubtractiveRule(std::vector<bool>& criterion){
        requireNodeCriterion(criterion, topology());
        return PybindUtils::toNumpy(AttributeFilters<std::uint8_t>::filteringBySubtractiveRule(criterion));

    }

    py::array_t<float> filteringBySubtractiveScoreRule(std::vector<float>& prob){
        requireNodeScores(prob, topology());
        return PybindUtils::toNumpy(AttributeFilters<std::uint8_t>::filteringBySubtractiveScoreRule(prob));

    }

    py::array_t<uint8_t> filteringByExtinctionValue(py::array attr, int leafToKeep){
        if (PybindUtils::parseFloatingArrayDType(attr, "attr") == PybindUtils::FloatingDType::Float64) {
            return filteringByExtinctionValueTyped<double>(std::move(attr), leafToKeep);
        }
        return filteringByExtinctionValueTyped<float>(std::move(attr), leafToKeep);
    }

    py::array saliencyMapByExtinctionValue(py::array attr, int leafToKeep, bool unweighted=false){
        if (PybindUtils::parseFloatingArrayDType(attr, "attr") == PybindUtils::FloatingDType::Float64) {
            return saliencyMapByExtinctionValueTyped<double>(std::move(attr), leafToKeep, unweighted);
        }
        return saliencyMapByExtinctionValueTyped<float>(std::move(attr), leafToKeep, unweighted);
    }



};

} // namespace mmcfilters
