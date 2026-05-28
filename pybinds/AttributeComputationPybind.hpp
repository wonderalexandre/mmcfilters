#pragma once

#include "../mmcfilters/attributes/AttributeComputation.hpp"
#include "../mmcfilters/trees/WeightedMorphologicalTree.hpp"
#include "MorphologicalTreePybind.hpp"
#include "PybindUtils.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <algorithm>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

namespace mmcfilters {

/**
 * @brief Pybind wrapper exposing attribute computation utilities to Python.
 */
class AttributeComputationPybind {
    static MorphologicalTree& topologyOf(MorphologicalTreePybind& tree) {
        return tree;
    }

    static const MorphologicalTree& topologyOf(const WeightedMorphologicalTree<std::uint8_t>& tree) {
        return tree.topology();
    }

    template <class TreePtr>
    static int outputSize(const TreePtr& tree, NodeIdSpace outputSpace) {
        return topologyOf(*tree).getNodeIdSpaceSize(outputSpace);
    }

    static py::dict makeAttributeLayoutDict(const AttributeNames& attributeNames) {
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
        return dict;
    }

    static py::dict makeDeltaAttributeLayoutDict(const AttributeNamesWithDelta& attributeNames) {
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
        return dict;
    }

    template <std::floating_point Real, class TreePtr>
    static py::array_t<Real> computeSingleAttributeTyped(TreePtr tree, Attribute attribute, NodeIdSpace outputSpace) {
        auto [attributeNames, buffer] = AttributeComputation::computeSingleAttribute<Real>(
            *tree,
            attribute,
            outputSpace);
        return PybindUtils::toNumpyOwned(std::move(buffer), outputSize(tree, outputSpace));
    }

    template <class TreePtr>
    static py::array computeSingleAttributeImpl(TreePtr tree, Attribute attribute, NodeIdSpace outputSpace, py::object dtype) {
        if (PybindUtils::parseFloatingDType(std::move(dtype)) == PybindUtils::FloatingDType::Float64) {
            return computeSingleAttributeTyped<double>(std::move(tree), attribute, outputSpace);
        }
        return computeSingleAttributeTyped<float>(std::move(tree), attribute, outputSpace);
    }

    template <std::floating_point Real, class TreePtr>
    static std::pair<py::dict, py::array> computeSingleAttributeWithDeltaTyped(TreePtr tree, Attribute attribute, int delta, std::string padding, NodeIdSpace outputSpace) {
        auto [attributeNames, buffer] = AttributeComputation::computeSingleAttributeWithDelta<Real>(
            *tree,
            attribute,
            AltitudeDiff<std::uint8_t>{1},
            delta,
            std::move(padding),
            outputSpace);

        const int numAttribute = attributeNames.NUM_ATTRIBUTES;
        const int n = outputSize(tree, outputSpace);
        return std::make_pair(
            makeDeltaAttributeLayoutDict(attributeNames),
            PybindUtils::toNumpyOwned2D(std::move(buffer), n, numAttribute));
    }

    template <class TreePtr>
    static std::pair<py::dict, py::array> computeSingleAttributeWithDeltaImpl(TreePtr tree, Attribute attribute, int delta, std::string padding, NodeIdSpace outputSpace, py::object dtype) {
        if (PybindUtils::parseFloatingDType(std::move(dtype)) == PybindUtils::FloatingDType::Float64) {
            return computeSingleAttributeWithDeltaTyped<double>(std::move(tree), attribute, delta, std::move(padding), outputSpace);
        }
        return computeSingleAttributeWithDeltaTyped<float>(std::move(tree), attribute, delta, std::move(padding), outputSpace);
    }

    template <std::floating_point Real, class TreePtr>
    static py::array_t<Real> computeAttributeMappingTyped(TreePtr tree, Attribute attribute) {
        auto imagePtr = AttributeComputation::computeAttributeMapping<Real>(*tree, attribute);
        return PybindUtils::toNumpy(imagePtr);
    }

    template <class TreePtr>
    static py::array computeAttributeMappingImpl(TreePtr tree, Attribute attribute, py::object dtype) {
        if (PybindUtils::parseFloatingDType(std::move(dtype)) == PybindUtils::FloatingDType::Float64) {
            return computeAttributeMappingTyped<double>(std::move(tree), attribute);
        }
        return computeAttributeMappingTyped<float>(std::move(tree), attribute);
    }

    template <std::floating_point Real, class TreePtr>
    static std::pair<py::dict, py::array> computeAttributesFromListTyped(TreePtr tree, const std::vector<AttributeOrGroup>& attributes, NodeIdSpace outputSpace) {
        auto [attributeNames, buffer] = AttributeComputation::computeAttributes<Real>(
            *tree,
            attributes,
            outputSpace);

        const int numAttribute = attributeNames.NUM_ATTRIBUTES;
        const int n = outputSize(tree, outputSpace);
        return std::make_pair(
            makeAttributeLayoutDict(attributeNames),
            PybindUtils::toNumpyOwned2D(std::move(buffer), n, numAttribute));
    }

    template <class TreePtr>
    static std::pair<py::dict, py::array> computeAttributesFromListImpl(TreePtr tree, const std::vector<AttributeOrGroup>& attributes, NodeIdSpace outputSpace, py::object dtype) {
        if (PybindUtils::parseFloatingDType(std::move(dtype)) == PybindUtils::FloatingDType::Float64) {
            return computeAttributesFromListTyped<double>(std::move(tree), attributes, outputSpace);
        }
        return computeAttributesFromListTyped<float>(std::move(tree), attributes, outputSpace);
    }

    template <std::floating_point Real, class TreePtr>
    static py::array_t<Real> computeSingleTopologyAttributeTyped(TreePtr tree, Attribute attribute, NodeIdSpace outputSpace) {
        auto [attributeNames, buffer] = AttributeComputation::computeSingleTopologyAttribute<Real>(
            *tree,
            attribute,
            outputSpace);
        return PybindUtils::toNumpyOwned(std::move(buffer), outputSize(tree, outputSpace));
    }

    template <class TreePtr>
    static py::array computeSingleTopologyAttributeImpl(TreePtr tree, Attribute attribute, NodeIdSpace outputSpace, py::object dtype) {
        if (PybindUtils::parseFloatingDType(std::move(dtype)) == PybindUtils::FloatingDType::Float64) {
            return computeSingleTopologyAttributeTyped<double>(std::move(tree), attribute, outputSpace);
        }
        return computeSingleTopologyAttributeTyped<float>(std::move(tree), attribute, outputSpace);
    }

    template <std::floating_point Real, class TreePtr>
    static std::pair<py::dict, py::array> computeTopologyAttributesFromListTyped(TreePtr tree, const std::vector<AttributeOrGroup>& attributes, NodeIdSpace outputSpace) {
        auto [attributeNames, buffer] = AttributeComputation::computeTopologyAttributes<Real>(
            *tree,
            attributes,
            outputSpace);

        const int numAttribute = attributeNames.NUM_ATTRIBUTES;
        const int n = outputSize(tree, outputSpace);

        return std::make_pair(
            makeAttributeLayoutDict(attributeNames),
            PybindUtils::toNumpyOwned2D(std::move(buffer), n, numAttribute));
    }

    template <class TreePtr>
    static std::pair<py::dict, py::array> computeTopologyAttributesFromListImpl(TreePtr tree, const std::vector<AttributeOrGroup>& attributes, NodeIdSpace outputSpace, py::object dtype) {
        if (PybindUtils::parseFloatingDType(std::move(dtype)) == PybindUtils::FloatingDType::Float64) {
            return computeTopologyAttributesFromListTyped<double>(std::move(tree), attributes, outputSpace);
        }
        return computeTopologyAttributesFromListTyped<float>(std::move(tree), attributes, outputSpace);
    }

public:
    static std::string describeAttribute(Attribute attribute) {
        return AttributeNames::describe(attribute);
    }

    static py::array computeSingleAttribute(std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> tree, Attribute attribute, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE, py::object dtype = py::none()){
        return computeSingleAttributeImpl(std::move(tree), attribute, outputSpace, std::move(dtype));
    }

    static std::pair<py::dict, py::array> computeSingleAttributeWithDelta(std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> tree, Attribute attribute, int delta, std::string padding = "last-padding", NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE, py::object dtype = py::none()) {
        return computeSingleAttributeWithDeltaImpl(std::move(tree), attribute, delta, std::move(padding), outputSpace, std::move(dtype));
    }

    static py::array computeAttributeMapping(std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> tree, Attribute attribute, py::object dtype = py::none()) {
        return computeAttributeMappingImpl(std::move(tree), attribute, std::move(dtype));
    }

    static std::pair<py::dict, py::array> computeAttributesFromList(std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> tree, const std::vector<AttributeOrGroup>& attributes, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE, py::object dtype = py::none()) {
        return computeAttributesFromListImpl(std::move(tree), attributes, outputSpace, std::move(dtype));
    }

    static py::array computeSingleTopologyAttribute(MorphologicalTreePybindPtr tree, Attribute attribute, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE, py::object dtype = py::none()) {
        return computeSingleTopologyAttributeImpl(std::move(tree), attribute, outputSpace, std::move(dtype));
    }

    static py::array computeSingleTopologyAttribute(std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> tree, Attribute attribute, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE, py::object dtype = py::none()) {
        return computeSingleTopologyAttributeImpl(std::move(tree), attribute, outputSpace, std::move(dtype));
    }

    static std::pair<py::dict, py::array> computeTopologyAttributesFromList(MorphologicalTreePybindPtr tree, const std::vector<AttributeOrGroup>& attributes, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE, py::object dtype = py::none()) {
        return computeTopologyAttributesFromListImpl(std::move(tree), attributes, outputSpace, std::move(dtype));
    }

    static std::pair<py::dict, py::array> computeTopologyAttributesFromList(std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> tree, const std::vector<AttributeOrGroup>& attributes, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE, py::object dtype = py::none()) {
        return computeTopologyAttributesFromListImpl(std::move(tree), attributes, outputSpace, std::move(dtype));
    }

};

} // namespace mmcfilters
