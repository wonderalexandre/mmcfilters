#pragma once

#include "../mmcfilters/attributes/AttributeComputation.hpp"
#include "../mmcfilters/trees/WeightedMorphologicalTree.hpp"
#include "PybindConversions.hpp"

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

namespace mmcfilters::pybindings::attribute_computation {

namespace py = pybind11;

/**
 * @brief Conversion and dispatch functions used by the Python attribute API.
 *
 * @param tree Tree topology used by the operation.
 * @return Reference to the resulting object.
 */
inline const MorphologicalTree& topologyOf(const WeightedMorphologicalTree<std::uint8_t>& tree) { return tree.topology(); }

/**
 * @brief Returns the number of values required by the selected output node space.
 *
 * @param tree Tree topology used by the operation.
 * @param outputSpace Node-id domain used to index the output.
 * @return Required number of output values.
 */
template <class TreePtr> inline int outputSize(const TreePtr& tree, NodeIdSpace outputSpace) { return topologyOf(*tree).getNodeIdSpaceSize(outputSpace); }

/**
 * @brief Creates attribute layout dict.
 *
 * @param attributeNames Attribute information represented by `attributeNames`.
 * @return Created attribute layout dict.
 */
inline py::dict makeAttributeLayoutDict(const AttributeNames& attributeNames) {
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

/**
 * @brief Creates delta attribute layout dict.
 *
 * @param attributeNames Attribute information represented by `attributeNames`.
 * @return Created delta attribute layout dict.
 */
inline py::dict makeDeltaAttributeLayoutDict(const AttributeNamesWithDelta& attributeNames) {
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

/**
 * @brief Computes single attribute typed.
 *
 * @param tree Tree topology used by the operation.
 * @param attribute Attribute requested by the operation.
 * @param outputSpace Node-id domain used to index the output.
 * @return Computed single attribute typed.
 */
template <std::floating_point Real, class TreePtr>
inline py::array_t<Real> computeSingleAttributeTyped(TreePtr tree, Attribute attribute, NodeIdSpace outputSpace) {
    auto [attributeNames, buffer] = AttributeComputation::computeSingleAttribute<Real>(*tree, attribute, outputSpace);
    return pybind_utils::toNumpyOwned(std::move(buffer), outputSize(tree, outputSpace));
}

/**
 * @brief Computes single attribute impl.
 *
 * @param tree Tree topology used by the operation.
 * @param attribute Attribute requested by the operation.
 * @param outputSpace Node-id domain used to index the output.
 * @param dtype Requested NumPy floating-point type.
 * @return Computed single attribute impl.
 */
template <class TreePtr> inline py::array computeSingleAttributeImpl(TreePtr tree, Attribute attribute, NodeIdSpace outputSpace, py::object dtype) {
    if (pybind_utils::parseFloatingDType(std::move(dtype)) == pybind_utils::FloatingDType::Float64) {
        return computeSingleAttributeTyped<double>(std::move(tree), attribute, outputSpace);
    }
    return computeSingleAttributeTyped<float>(std::move(tree), attribute, outputSpace);
}

/**
 * @brief Computes single attribute with delta typed.
 *
 * @param tree Tree topology used by the operation.
 * @param attribute Attribute requested by the operation.
 * @param delta Delta offset used by the operation.
 * @param padding Padding strategy used by the operation.
 * @param outputSpace Node-id domain used to index the output.
 * @return Computed single attribute with delta typed.
 */
template <std::floating_point Real, class TreePtr>
inline std::pair<py::dict, py::array> computeSingleAttributeWithDeltaTyped(TreePtr tree, Attribute attribute, int delta, std::string padding,
                                                                           NodeIdSpace outputSpace) {
    auto [attributeNames, buffer] =
        AttributeComputation::computeSingleAttributeWithDelta<Real>(*tree, attribute, AltitudeDiff<std::uint8_t>{1}, delta, std::move(padding), outputSpace);

    const int numAttribute = attributeNames.NUM_ATTRIBUTES;
    const int n = outputSize(tree, outputSpace);
    return std::make_pair(makeDeltaAttributeLayoutDict(attributeNames), pybind_utils::toNumpyOwned2D(std::move(buffer), n, numAttribute));
}

/**
 * @brief Computes single attribute with delta impl.
 *
 * @param tree Tree topology used by the operation.
 * @param attribute Attribute requested by the operation.
 * @param delta Delta offset used by the operation.
 * @param padding Padding strategy used by the operation.
 * @param outputSpace Node-id domain used to index the output.
 * @param dtype Requested NumPy floating-point type.
 * @return Computed single attribute with delta impl.
 */
template <class TreePtr>
inline std::pair<py::dict, py::array> computeSingleAttributeWithDeltaImpl(TreePtr tree, Attribute attribute, int delta, std::string padding,
                                                                          NodeIdSpace outputSpace, py::object dtype) {
    if (pybind_utils::parseFloatingDType(std::move(dtype)) == pybind_utils::FloatingDType::Float64) {
        return computeSingleAttributeWithDeltaTyped<double>(std::move(tree), attribute, delta, std::move(padding), outputSpace);
    }
    return computeSingleAttributeWithDeltaTyped<float>(std::move(tree), attribute, delta, std::move(padding), outputSpace);
}

/**
 * @brief Computes attribute mapping typed.
 *
 * @param tree Tree topology used by the operation.
 * @param attribute Attribute requested by the operation.
 * @return Computed attribute mapping typed.
 */
template <std::floating_point Real, class TreePtr> inline py::array_t<Real> computeAttributeMappingTyped(TreePtr tree, Attribute attribute) {
    auto imagePtr = AttributeComputation::computeAttributeMapping<Real>(*tree, attribute);
    return pybind_utils::toNumpy(imagePtr);
}

/**
 * @brief Computes attribute mapping impl.
 *
 * @param tree Tree topology used by the operation.
 * @param attribute Attribute requested by the operation.
 * @param dtype Requested NumPy floating-point type.
 * @return Computed attribute mapping impl.
 */
template <class TreePtr> inline py::array computeAttributeMappingImpl(TreePtr tree, Attribute attribute, py::object dtype) {
    if (pybind_utils::parseFloatingDType(std::move(dtype)) == pybind_utils::FloatingDType::Float64) {
        return computeAttributeMappingTyped<double>(std::move(tree), attribute);
    }
    return computeAttributeMappingTyped<float>(std::move(tree), attribute);
}

/**
 * @brief Computes attributes from list typed.
 *
 * @param tree Tree topology used by the operation.
 * @param attributes Attributes requested by the operation.
 * @param outputSpace Node-id domain used to index the output.
 * @return Computed attributes from list typed.
 */
template <std::floating_point Real, class TreePtr>
inline std::pair<py::dict, py::array> computeAttributesFromListTyped(TreePtr tree, const std::vector<AttributeOrGroup>& attributes, NodeIdSpace outputSpace) {
    auto [attributeNames, buffer] = AttributeComputation::computeAttributes<Real>(*tree, attributes, outputSpace);

    const int numAttribute = attributeNames.NUM_ATTRIBUTES;
    const int n = outputSize(tree, outputSpace);
    return std::make_pair(makeAttributeLayoutDict(attributeNames), pybind_utils::toNumpyOwned2D(std::move(buffer), n, numAttribute));
}

/**
 * @brief Computes attributes from list impl.
 *
 * @param tree Tree topology used by the operation.
 * @param attributes Attributes requested by the operation.
 * @param outputSpace Node-id domain used to index the output.
 * @param dtype Requested NumPy floating-point type.
 * @return Computed attributes from list impl.
 */
template <class TreePtr>
inline std::pair<py::dict, py::array> computeAttributesFromListImpl(TreePtr tree, const std::vector<AttributeOrGroup>& attributes, NodeIdSpace outputSpace,
                                                                    py::object dtype) {
    if (pybind_utils::parseFloatingDType(std::move(dtype)) == pybind_utils::FloatingDType::Float64) {
        return computeAttributesFromListTyped<double>(std::move(tree), attributes, outputSpace);
    }
    return computeAttributesFromListTyped<float>(std::move(tree), attributes, outputSpace);
}

/**
 * @brief Computes single topology attribute typed.
 *
 * @param tree Tree topology used by the operation.
 * @param attribute Attribute requested by the operation.
 * @param outputSpace Node-id domain used to index the output.
 * @return Computed single topology attribute typed.
 */
template <std::floating_point Real, class TreePtr>
inline py::array_t<Real> computeSingleTopologyAttributeTyped(TreePtr tree, Attribute attribute, NodeIdSpace outputSpace) {
    auto [attributeNames, buffer] = AttributeComputation::computeSingleTopologyAttribute<Real>(*tree, attribute, outputSpace);
    return pybind_utils::toNumpyOwned(std::move(buffer), outputSize(tree, outputSpace));
}

/**
 * @brief Computes single topology attribute impl.
 *
 * @param tree Tree topology used by the operation.
 * @param attribute Attribute requested by the operation.
 * @param outputSpace Node-id domain used to index the output.
 * @param dtype Requested NumPy floating-point type.
 * @return Computed single topology attribute impl.
 */
template <class TreePtr> inline py::array computeSingleTopologyAttributeImpl(TreePtr tree, Attribute attribute, NodeIdSpace outputSpace, py::object dtype) {
    if (pybind_utils::parseFloatingDType(std::move(dtype)) == pybind_utils::FloatingDType::Float64) {
        return computeSingleTopologyAttributeTyped<double>(std::move(tree), attribute, outputSpace);
    }
    return computeSingleTopologyAttributeTyped<float>(std::move(tree), attribute, outputSpace);
}

/**
 * @brief Computes topology attributes from list typed.
 *
 * @param tree Tree topology used by the operation.
 * @param attributes Attributes requested by the operation.
 * @param outputSpace Node-id domain used to index the output.
 * @return Computed topology attributes from list typed.
 */
template <std::floating_point Real, class TreePtr>
inline std::pair<py::dict, py::array> computeTopologyAttributesFromListTyped(TreePtr tree, const std::vector<AttributeOrGroup>& attributes,
                                                                             NodeIdSpace outputSpace) {
    auto [attributeNames, buffer] = AttributeComputation::computeTopologyAttributes<Real>(*tree, attributes, outputSpace);

    const int numAttribute = attributeNames.NUM_ATTRIBUTES;
    const int n = outputSize(tree, outputSpace);

    return std::make_pair(makeAttributeLayoutDict(attributeNames), pybind_utils::toNumpyOwned2D(std::move(buffer), n, numAttribute));
}

/**
 * @brief Computes topology attributes from list impl.
 *
 * @param tree Tree topology used by the operation.
 * @param attributes Attributes requested by the operation.
 * @param outputSpace Node-id domain used to index the output.
 * @param dtype Requested NumPy floating-point type.
 * @return Computed topology attributes from list impl.
 */
template <class TreePtr>
inline std::pair<py::dict, py::array> computeTopologyAttributesFromListImpl(TreePtr tree, const std::vector<AttributeOrGroup>& attributes,
                                                                            NodeIdSpace outputSpace, py::object dtype) {
    if (pybind_utils::parseFloatingDType(std::move(dtype)) == pybind_utils::FloatingDType::Float64) {
        return computeTopologyAttributesFromListTyped<double>(std::move(tree), attributes, outputSpace);
    }
    return computeTopologyAttributesFromListTyped<float>(std::move(tree), attributes, outputSpace);
}

/**
 * @brief Describes attribute.
 *
 * @param attribute Attribute requested by the operation.
 * @return Human-readable attribute description.
 */
inline std::string describeAttribute(Attribute attribute) { return AttributeNames::describe(attribute); }

/**
 * @brief Computes single attribute.
 *
 * @param tree Tree topology used by the operation.
 * @param attribute Attribute requested by the operation.
 * @param outputSpace Node-id domain used to index the output.
 * @param dtype Requested NumPy floating-point type.
 * @return Computed single attribute.
 */
inline py::array computeSingleAttribute(std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> tree, Attribute attribute,
                                        NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE, py::object dtype = py::none()) {
    return computeSingleAttributeImpl(std::move(tree), attribute, outputSpace, std::move(dtype));
}

/**
 * @brief Computes single attribute with delta.
 *
 * @param tree Tree topology used by the operation.
 * @param attribute Attribute requested by the operation.
 * @param delta Delta offset used by the operation.
 * @param padding Padding strategy used by the operation.
 * @param outputSpace Node-id domain used to index the output.
 * @param dtype Requested NumPy floating-point type.
 * @return Computed single attribute with delta.
 */
inline std::pair<py::dict, py::array> computeSingleAttributeWithDelta(std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> tree, Attribute attribute,
                                                                      int delta, std::string padding = "last-padding",
                                                                      NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE,
                                                                      py::object dtype = py::none()) {
    return computeSingleAttributeWithDeltaImpl(std::move(tree), attribute, delta, std::move(padding), outputSpace, std::move(dtype));
}

/**
 * @brief Computes attribute mapping.
 *
 * @param tree Tree topology used by the operation.
 * @param attribute Attribute requested by the operation.
 * @param dtype Requested NumPy floating-point type.
 * @return Computed attribute mapping.
 */
inline py::array computeAttributeMapping(std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> tree, Attribute attribute, py::object dtype = py::none()) {
    return computeAttributeMappingImpl(std::move(tree), attribute, std::move(dtype));
}

/**
 * @brief Computes attributes from list.
 *
 * @param tree Tree topology used by the operation.
 * @param attributes Attributes requested by the operation.
 * @param outputSpace Node-id domain used to index the output.
 * @param dtype Requested NumPy floating-point type.
 * @return Computed attributes from list.
 */
inline std::pair<py::dict, py::array> computeAttributesFromList(std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> tree,
                                                                const std::vector<AttributeOrGroup>& attributes,
                                                                NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE, py::object dtype = py::none()) {
    return computeAttributesFromListImpl(std::move(tree), attributes, outputSpace, std::move(dtype));
}

/**
 * @brief Computes single topology attribute.
 *
 * @param tree Tree topology used by the operation.
 * @param attribute Attribute requested by the operation.
 * @param outputSpace Node-id domain used to index the output.
 * @param dtype Requested NumPy floating-point type.
 * @return Computed single topology attribute.
 */
inline py::array computeSingleTopologyAttribute(std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> tree, Attribute attribute,
                                                NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE, py::object dtype = py::none()) {
    return computeSingleTopologyAttributeImpl(std::move(tree), attribute, outputSpace, std::move(dtype));
}

/**
 * @brief Computes topology attributes from list.
 *
 * @param tree Tree topology used by the operation.
 * @param attributes Attributes requested by the operation.
 * @param outputSpace Node-id domain used to index the output.
 * @param dtype Requested NumPy floating-point type.
 * @return Computed topology attributes from list.
 */
inline std::pair<py::dict, py::array> computeTopologyAttributesFromList(std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> tree,
                                                                        const std::vector<AttributeOrGroup>& attributes,
                                                                        NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE,
                                                                        py::object dtype = py::none()) {
    return computeTopologyAttributesFromListImpl(std::move(tree), attributes, outputSpace, std::move(dtype));
}

} // namespace mmcfilters::pybindings::attribute_computation
