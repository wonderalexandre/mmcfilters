#pragma once

#include "../mmcfilters/attributes/AttributeComputedIncrementally.hpp"
#include "../mmcfilters/trees/WeightedMorphologicalTree.hpp"
#include "MorphologicalTreePybind.hpp"
#include "PybindUtils.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <algorithm> 
#include <cmath>
#include <iostream>
#include <numeric>
namespace mmcfilters {

class AttributeComputedIncrementallyPybind;
using AttributeComputedIncrementallyPybindPtr = std::shared_ptr<AttributeComputedIncrementallyPybind>;

/**
 * @brief Pybind wrapper exposing incremental attribute utilities to Python.
 */
class AttributeComputedIncrementallyPybind : public AttributeComputedIncrementally{
    static MorphologicalTree& topologyOf(MorphologicalTreePybind& tree) {
        return tree;
    }

    static const MorphologicalTree& topologyOf(const WeightedMorphologicalTree& tree) {
        return tree.topology();
    }

    template <class TreePtr>
    static int outputSize(const TreePtr& tree, NodeIdSpace outputSpace) {
        return topologyOf(*tree).getNodeIdSpaceSize(outputSpace);
    }

    template <class TreePtr>
    static py::array_t<float> computeSingleAttributeImpl(TreePtr tree, Attribute attribute, NodeIdSpace outputSpace) {
        auto [attributeNames, buffer] = AttributeComputedIncrementally::computeSingleAttribute(*tree, attribute, {}, outputSpace);
        return PybindUtils::toNumpyOwned(std::move(buffer), outputSize(tree, outputSpace));
    }

    template <class TreePtr>
    static std::pair<py::dict, py::array_t<float>> computeSingleAttributeWithDeltaImpl(TreePtr tree, Attribute attribute, int delta, std::string padding, NodeIdSpace outputSpace) {
        auto [attributeNames, buffer] = AttributeComputedIncrementally::computeSingleAttributeWithDelta(*tree, attribute, delta, padding, /*deps*/{}, outputSpace);

        const int numAttribute = attributeNames.NUM_ATTRIBUTES;
        const int n = outputSize(tree, outputSpace);

        std::vector<std::string> keys;
        std::vector<int> values;

        for (const auto& pair : attributeNames.indexMap) {
            const AttributeKey& attrKey = pair.first;
            int offset = pair.second;
            keys.push_back(AttributeNamesWithDelta::toString(attrKey.attr, attrKey.delta));
            values.push_back(offset);
        }

        std::vector<size_t> indices(values.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), [&values](size_t i1, size_t i2) { return values[i1] < values[i2]; });

        py::dict dict;
        for (size_t i = 0; i < indices.size(); ++i) {
            dict[py::str(keys[indices[i]])] = values[indices[i]];
        }

        return std::make_pair(dict, PybindUtils::toNumpyOwned2D(std::move(buffer), n, numAttribute));
    }

    template <class TreePtr>
    static py::array_t<float> computeAttributeMappingImpl(TreePtr tree, Attribute attribute) {
        auto imgFloatPtr = AttributeComputedIncrementally::computeAttributeMapping(*tree, attribute);
        return PybindUtils::toNumpy(imgFloatPtr);
    }

    template <class TreePtr>
    static std::pair<py::dict, py::array_t<float>> computeAttributesFromListImpl(TreePtr tree, const std::vector<AttributeOrGroup>& attributes, NodeIdSpace outputSpace) {
        auto [attributeNames, buffer] = AttributeComputedIncrementally::computeAttributes(*tree, attributes, {}, outputSpace);

        const int numAttribute = attributeNames.NUM_ATTRIBUTES;
        const int n = outputSize(tree, outputSpace);

        std::vector<std::string> keys;
        std::vector<int> values;
        for (const auto& pair : attributeNames.indexMap) {
            Attribute attribute = pair.first;
            int offset = pair.second;

            keys.push_back(attributeNames.toString(attribute));
            values.push_back(offset);
        }

        std::vector<size_t> indices(values.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), [&values](size_t i1, size_t i2) { return values[i1] < values[i2]; });

        py::dict dict;
        for (size_t i = 0; i < indices.size(); ++i) {
            dict[py::str(keys[indices[i]])] = values[indices[i]];
        }

        return std::make_pair(dict, PybindUtils::toNumpyOwned2D(std::move(buffer), n, numAttribute));
    }

    public:
    using AttributeComputedIncrementally::AttributeComputedIncrementally;

	

	static std::string describeAttribute(Attribute attribute) {
		return AttributeNames::describe(attribute);
	}

	static py::array_t<float> computeSingleAttribute(MorphologicalTreePybindPtr tree, Attribute attribute, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE){
		return computeSingleAttributeImpl(std::move(tree), attribute, outputSpace);
	}

    static py::array_t<float> computeSingleAttribute(std::shared_ptr<WeightedMorphologicalTree> tree, Attribute attribute, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE){
        return computeSingleAttributeImpl(std::move(tree), attribute, outputSpace);
    }

	static std::pair<py::dict, py::array_t<float>> computeSingleAttributeWithDelta(MorphologicalTreePybindPtr tree, Attribute attribute, int delta, std::string padding = "last-padding", NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE) {
		return computeSingleAttributeWithDeltaImpl(std::move(tree), attribute, delta, std::move(padding), outputSpace);
	}

    static std::pair<py::dict, py::array_t<float>> computeSingleAttributeWithDelta(std::shared_ptr<WeightedMorphologicalTree> tree, Attribute attribute, int delta, std::string padding = "last-padding", NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE) {
        return computeSingleAttributeWithDeltaImpl(std::move(tree), attribute, delta, std::move(padding), outputSpace);
    }

	static py::array_t<float> computeAttributeMapping(MorphologicalTreePybindPtr tree, Attribute attribute) {
		return computeAttributeMappingImpl(std::move(tree), attribute);
	} 

    static py::array_t<float> computeAttributeMapping(std::shared_ptr<WeightedMorphologicalTree> tree, Attribute attribute) {
        return computeAttributeMappingImpl(std::move(tree), attribute);
    }

	static std::pair<py::dict, py::array_t<float>> computeAttributesFromList(MorphologicalTreePybindPtr tree, const std::vector<AttributeOrGroup>& attributes, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE) {
		return computeAttributesFromListImpl(std::move(tree), attributes, outputSpace);
	}

    static std::pair<py::dict, py::array_t<float>> computeAttributesFromList(std::shared_ptr<WeightedMorphologicalTree> tree, const std::vector<AttributeOrGroup>& attributes, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE) {
        return computeAttributesFromListImpl(std::move(tree), attributes, outputSpace);
    }

};

} // namespace mmcfilters
