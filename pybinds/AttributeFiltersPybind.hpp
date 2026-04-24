#pragma once

#include "../mmcfilters/filters/AttributeFilters.hpp"
#include "../mmcfilters/utils/Common.hpp"

#include "MorphologicalTreePybind.hpp"
#include "AttributeComputedIncrementallyPybind.hpp"
#include "PybindUtils.hpp"

#include <stack>
#include <vector>
#include <limits.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace mmcfilters {

#define UNDEF -999999999999

/**
 * @brief *Wrapper* Pybind11 para expor filtragens por atributos ao Python.
 */
class AttributeFiltersPybind : public AttributeFilters{
    using FloatArray = py::array_t<float, py::array::c_style | py::array::forcecast>;

    MorphologicalTreePybindPtr treeOwner_;
    std::shared_ptr<WeightedMorphologicalTree> weightedOwner_;

    const MorphologicalTree& topology() const noexcept {
        return this->tree;
    }

    static void requireNodeAttributeArray(const FloatArray& attr, const MorphologicalTree& tree, std::string_view argumentName = "attr") {
        PybindUtils::require1DArray(attr.request(), tree.getNumInternalNodeSlots(), argumentName);
    }

    static void requireNodeCriterion(const std::vector<bool>& criterion, const MorphologicalTree& tree, std::string_view argumentName = "criterion") {
        PybindUtils::requireVectorSize(criterion, static_cast<std::size_t>(tree.getNumInternalNodeSlots()), argumentName);
    }

    static void requireNodeScores(const std::vector<float>& prob, const MorphologicalTree& tree, std::string_view argumentName = "prob") {
        PybindUtils::requireVectorSize(prob, static_cast<std::size_t>(tree.getNumInternalNodeSlots()), argumentName);
    }

    public:
    using AttributeFilters::AttributeFilters;

    explicit AttributeFiltersPybind(MorphologicalTreePybindPtr tree)
        : AttributeFilters(*tree), treeOwner_(std::move(tree)) {}

    explicit AttributeFiltersPybind(std::shared_ptr<WeightedMorphologicalTree> weighted)
        : AttributeFilters(*weighted), weightedOwner_(std::move(weighted)) {}

    py::array_t<uint8_t> filteringByPruningMin(FloatArray attr, float threshold){
        requireNodeAttributeArray(attr, topology());

        std::shared_ptr<float[]> attribute = PybindUtils::toShared_ptr(attr);
        return PybindUtils::toNumpy(AttributeFilters::filteringByPruningMin(attribute, threshold));
    }

    py::array_t<uint8_t> filteringByPruningMax(FloatArray attr, float threshold){
        requireNodeAttributeArray(attr, topology());

        std::shared_ptr<float[]> attribute = PybindUtils::toShared_ptr(attr);
        return PybindUtils::toNumpy(AttributeFilters::filteringByPruningMax(attribute, threshold));

    }

    py::array_t<uint8_t> filteringByPruningMin(std::vector<bool>& criterion){
        requireNodeCriterion(criterion, topology());
        return PybindUtils::toNumpy(AttributeFilters::filteringByPruningMin(criterion));
    }

    py::array_t<uint8_t> filteringByDirectRule(std::vector<bool>& criterion){
        requireNodeCriterion(criterion, topology());
        return PybindUtils::toNumpy(AttributeFilters::filteringByDirectRule(criterion));
    }

    py::array_t<uint8_t> filteringByPruningMax(std::vector<bool>& criterion){
        requireNodeCriterion(criterion, topology());
        return PybindUtils::toNumpy(AttributeFilters::filteringByPruningMax(criterion));

    }

    std::vector<bool> getAdaptiveCriterion(std::vector<bool>& criterion, int delta){
        requireNodeCriterion(criterion, topology());
        return AttributeFilters::getAdaptiveCriterion(criterion, delta);
    }



    py::array_t<uint8_t> filteringBySubtractiveRule(std::vector<bool>& criterion){
        requireNodeCriterion(criterion, topology());
        return PybindUtils::toNumpy(AttributeFilters::filteringBySubtractiveRule(criterion));

    }

    py::array_t<float> filteringBySubtractiveScoreRule(std::vector<float>& prob){
        requireNodeScores(prob, topology());
        return PybindUtils::toNumpy(AttributeFilters::filteringBySubtractiveScoreRule(prob));

    }

    py::array_t<uint8_t> filteringByExtinctionValue(FloatArray attr, int leafToKeep){
        requireNodeAttributeArray(attr, topology());
        if (this->weightedOwner_ != nullptr) {
            ExtinctionValues ev(*this->weightedOwner_, PybindUtils::toShared_ptr(attr));
            return PybindUtils::toNumpy(ev.filtering(leafToKeep));
        }
        ExtinctionValues ev(this->tree, PybindUtils::toShared_ptr(attr));
        return PybindUtils::toNumpy(ev.filtering(leafToKeep));
    }

    py::array_t<float> saliencyMapByExtinctionValue(FloatArray attr, int leafToKeep, bool unweighted=false){
        requireNodeAttributeArray(attr, topology());
        if (this->weightedOwner_ != nullptr) {
            ExtinctionValues ev(*this->weightedOwner_, PybindUtils::toShared_ptr(attr));
            return PybindUtils::toNumpy(ev.saliencyMap(leafToKeep, unweighted));
        }
        ExtinctionValues ev(this->tree, PybindUtils::toShared_ptr(attr));
        return PybindUtils::toNumpy(ev.saliencyMap(leafToKeep, unweighted));
    }



};

} // namespace mmcfilters
