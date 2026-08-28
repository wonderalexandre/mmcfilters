#pragma once

#include "PythonValuedMorphologicalTree.hpp"
#include "../mmcfilters/attributes/AttributeComputation.hpp"
#include "../mmcfilters/attributes/AttributeRegistry.hpp"
#include "../mmcfilters/trees/ValuedMorphologicalTree.hpp"
#include "PybindConversions.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <algorithm>
#include <cctype>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <numeric>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace mmcfilters::pybindings::attribute_computation {

namespace py = pybind11;

/// @cond INTERNAL
namespace detail {

/** @brief Stable symbolic name of every public attribute group. @return Ordered name/group pairs. */
inline const std::vector<std::pair<std::string, AttributeGroup>>& attributeGroupNames() {
    static const std::vector<std::pair<std::string, AttributeGroup>> names = {
        {"ALL", AttributeGroup::All},           {"GRAY_LEVEL", AttributeGroup::GrayLevel},
        {"SHAPE", AttributeGroup::Shape},       {"MOMENTS", AttributeGroup::Moments},
        {"BOUNDARY", AttributeGroup::Boundary}, {"TREE_TOPOLOGY", AttributeGroup::TreeTopology},
        {"DIST_TRANSF", AttributeGroup::DistTransf}, {"DIST_TRANSF_EXACT", AttributeGroup::DistTransfExact},
    };
    return names;
}

/** @brief Case-insensitive uppercase form used to suggest near matches. @param text Symbolic name. @return Uppercase form. */
inline std::string upperCased(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(), [](unsigned char character) { return static_cast<char>(std::toupper(character)); });
    return text;
}

/**
 * @brief Reports an unknown symbolic name together with the closest known names.
 *
 * @param name Rejected symbolic name.
 * @param allowGroups Whether attribute-group names were also admissible.
 */
[[noreturn]] inline void reportUnknownAttributeName(const std::string& name, bool allowGroups) {
    const std::string wanted = upperCased(name);
    std::vector<std::string> suggestions;
    for (const attributes::registry::AttributeMetadata& item : attributes::registry::ATTRIBUTE_METADATA) {
        const std::string candidate(item.name);
        if (upperCased(candidate) == wanted || candidate.find(wanted) != std::string::npos || wanted.find(candidate) != std::string::npos) {
            suggestions.push_back(candidate);
        }
    }
    if (allowGroups) {
        for (const auto& [candidate, group] : attributeGroupNames()) {
            static_cast<void>(group);
            if (upperCased(candidate) == wanted) {
                suggestions.push_back(candidate);
            }
        }
    }

    std::string message = "Unknown attribute name '" + name + "'.";
    if (!suggestions.empty()) {
        message += " Did you mean ";
        for (std::size_t index = 0; index < suggestions.size() && index < 5; ++index) {
            message += (index == 0 ? "'" : ", '") + suggestions[index] + "'";
        }
        message += "?";
    } else {
        message += " Names are the stable symbolic names, upper case, such as 'AREA' or 'BOUNDING_BOX_HEIGHT'.";
    }
    throw py::value_error(message);
}

} // namespace detail
/// @endcond

/**
 * @brief Resolves one attribute given either an `Attribute` value or its symbolic name.
 *
 * Accepting the name keeps call sites short without introducing a second
 * vocabulary: the accepted strings are exactly the stable symbolic names that
 * the returned attribute layouts already use as keys.
 *
 * @param attribute `Attribute` value or `str`.
 * @return The resolved attribute.
 */
[[nodiscard]] inline Attribute resolveAttribute(const py::object& attribute) {
    if (py::isinstance<py::str>(attribute)) {
        const auto name = attribute.cast<std::string>();
        if (const std::optional<Attribute> parsed = attributes::registry::parse(name); parsed.has_value()) {
            return *parsed;
        }
        detail::reportUnknownAttributeName(name, false);
    }
    try {
        return attribute.cast<Attribute>();
    } catch (const py::cast_error&) {
        throw py::type_error("Expected an Attribute value or its symbolic name as str.");
    }
}

/**
 * @brief Resolves a sequence of attributes and groups given values or symbolic names.
 *
 * @param attributes Iterable of `Attribute`, `Attribute.Group`, or `str`.
 * @return The resolved request sequence.
 */
[[nodiscard]] inline std::vector<AttributeOrGroup> resolveAttributeOrGroupList(const py::object& attributes) {
    if (py::isinstance<py::str>(attributes)) {
        throw py::type_error("Expected a sequence of attributes; pass a list such as ['AREA'] rather than a bare str.");
    }
    std::vector<AttributeOrGroup> resolved;
    for (const py::handle item : attributes) {
        const auto entry = py::reinterpret_borrow<py::object>(item);
        if (py::isinstance<py::str>(entry)) {
            const auto name = entry.cast<std::string>();
            if (const std::optional<Attribute> parsed = attributes::registry::parse(name); parsed.has_value()) {
                resolved.emplace_back(*parsed);
                continue;
            }
            bool matchedGroup = false;
            for (const auto& [candidate, group] : detail::attributeGroupNames()) {
                if (candidate == name) {
                    resolved.emplace_back(group);
                    matchedGroup = true;
                    break;
                }
            }
            if (matchedGroup) {
                continue;
            }
            detail::reportUnknownAttributeName(name, true);
        }
        try {
            resolved.push_back(entry.cast<AttributeOrGroup>());
        } catch (const py::cast_error&) {
            throw py::type_error("Expected Attribute, Attribute.Group, or str entries in the attribute sequence.");
        }
    }
    return resolved;
}

/**
 * @brief Conversion and dispatch functions used by the Python attribute API.
 *
 * @param tree Tree topology.
 * @return Mutable reference to the updated object.
 */
template <AltitudeValue T> inline const MorphologicalTree& topologyOf(const ValuedMorphologicalTree<T>& tree) { return tree.topology(); }

inline const MorphologicalTree& topologyOf(const PythonValuedMorphologicalTree& tree) { return tree.topology(); }

/**
 * @brief Returns the number of values required by the selected output node space.
 *
 * @param tree Tree topology.
 * @param outputSpace Node-id domain used to index the output.
 * @return Required number of output values.
 */
template <class TreePtr> inline int outputSize(const TreePtr& tree, NodeIdSpace outputSpace) { return topologyOf(*tree).getNodeIdSpaceSize(outputSpace); }

/**
 * @brief Creates attribute layout dict.
 *
 * @param attributeNames Attribute information.
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
 * @brief Creates a sampled node-attribute layout dictionary.
 *
 * @param attributeNames Attribute information.
 * @return Created sampled node-attribute layout dictionary.
 */
inline py::dict makeNodeAttributeSampleLayoutDict(const NodeAttributeSampleLayout& attributeNames) {
    std::vector<std::string> keys;
    std::vector<int> values;

    for (const auto& pair : attributeNames.indexMap) {
        const NodeAttributeSampleKey& sampleKey = pair.first;
        int offset = pair.second;
        keys.push_back(NodeAttributeSampleLayout::toString(sampleKey.attribute, sampleKey.sampleOffset));
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
 * @param tree Tree topology.
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
 * @param tree Tree topology.
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

/** @brief Computes sampled node-attribute data for one concrete altitude type. */
template <std::floating_point Real, class TreePtr>
inline std::pair<py::dict, py::array>
computeSampledNodeAttributeTyped(TreePtr tree, Attribute attribute, std::int64_t altitudeStep, int samplingRadius,
                                 NodeAttributeSamplingPolicy samplingPolicy, MissingNodeAttributeSamplePolicy missingSamplePolicy,
                                 NodeIdSpace outputSpace) {
    using Altitude = typename std::remove_cvref_t<decltype(*tree)>::AltitudeType;
    auto [attributeNames, buffer] = AttributeComputation::computeSampledNodeAttribute<Real>(
        *tree, attribute, static_cast<AltitudeDifference<Altitude>>(altitudeStep), samplingRadius, samplingPolicy, missingSamplePolicy, outputSpace);

    const int numAttribute = attributeNames.NUM_ATTRIBUTES;
    const int n = outputSize(tree, outputSpace);
    return std::make_pair(makeNodeAttributeSampleLayoutDict(attributeNames), pybind_utils::toNumpyOwned2D(std::move(buffer), n, numAttribute));
}

/** @brief Dispatches sampled node-attribute result precision. */
template <class TreePtr>
inline std::pair<py::dict, py::array>
computeSampledNodeAttributeImpl(TreePtr tree, Attribute attribute, std::int64_t altitudeStep, int samplingRadius,
                                NodeAttributeSamplingPolicy samplingPolicy, MissingNodeAttributeSamplePolicy missingSamplePolicy,
                                NodeIdSpace outputSpace, py::object dtype) {
    if (pybind_utils::parseFloatingDType(std::move(dtype)) == pybind_utils::FloatingDType::Float64) {
        return computeSampledNodeAttributeTyped<double>(std::move(tree), attribute, altitudeStep, samplingRadius, samplingPolicy, missingSamplePolicy,
                                                        outputSpace);
    }
    return computeSampledNodeAttributeTyped<float>(std::move(tree), attribute, altitudeStep, samplingRadius, samplingPolicy, missingSamplePolicy,
                                                   outputSpace);
}

/**
 * @brief Computes attribute mapping typed.
 *
 * @param tree Tree topology.
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
 * @param tree Tree topology.
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
 * @param tree Tree topology.
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
 * @param tree Tree topology.
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
 * @param tree Tree topology.
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
 * @param tree Tree topology.
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
 * @param tree Tree topology.
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
 * @param tree Tree topology.
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
 * @param tree Tree topology.
 * @param attribute Attribute requested by the operation.
 * @param outputSpace Node-id domain used to index the output.
 * @param dtype Requested NumPy floating-point type.
 * @return Computed single attribute.
 */
inline py::array computeSingleAttribute(std::shared_ptr<PythonValuedMorphologicalTree> tree, Attribute attribute,
                                        NodeIdSpace outputSpace = NodeIdSpace::MorphologicalTree, py::object dtype = py::none()) {
    return tree->visit([&](const auto& concreteTree) { return computeSingleAttributeImpl(concreteTree, attribute, outputSpace, dtype); });
}

/** @brief Computes altitude-based samples of one node attribute. */
inline std::pair<py::dict, py::array>
computeSampledNodeAttribute(std::shared_ptr<PythonValuedMorphologicalTree> tree, Attribute attribute, std::int64_t altitudeStep,
                            int samplingRadius,
                            NodeAttributeSamplingPolicy samplingPolicy = NodeAttributeSamplingPolicy::LargestSupportDescendant,
                            MissingNodeAttributeSamplePolicy missingSamplePolicy = MissingNodeAttributeSamplePolicy::RepeatNearest,
                            NodeIdSpace outputSpace = NodeIdSpace::MorphologicalTree, py::object dtype = py::none()) {
    return tree->visit([&](const auto& concreteTree) {
        return computeSampledNodeAttributeImpl(concreteTree, attribute, altitudeStep, samplingRadius, samplingPolicy, missingSamplePolicy, outputSpace,
                                               dtype);
    });
}

/**
 * @brief Computes attribute mapping.
 *
 * @param tree Tree topology.
 * @param attribute Attribute requested by the operation.
 * @param dtype Requested NumPy floating-point type.
 * @return Computed attribute mapping.
 */
inline py::array computeAttributeMapping(std::shared_ptr<PythonValuedMorphologicalTree> tree, Attribute attribute, py::object dtype = py::none()) {
    return tree->visit([&](const auto& concreteTree) { return computeAttributeMappingImpl(concreteTree, attribute, dtype); });
}

/**
 * @brief Computes attributes from list.
 *
 * @param tree Tree topology.
 * @param attributes Attributes requested by the operation.
 * @param outputSpace Node-id domain used to index the output.
 * @param dtype Requested NumPy floating-point type.
 * @return Computed attributes from list.
 */
inline std::pair<py::dict, py::array> computeAttributesFromList(std::shared_ptr<PythonValuedMorphologicalTree> tree,
                                                                const std::vector<AttributeOrGroup>& attributes,
                                                                NodeIdSpace outputSpace = NodeIdSpace::MorphologicalTree, py::object dtype = py::none()) {
    return tree->visit([&](const auto& concreteTree) { return computeAttributesFromListImpl(concreteTree, attributes, outputSpace, dtype); });
}

/**
 * @brief Computes single topology attribute.
 *
 * @param tree Tree topology.
 * @param attribute Attribute requested by the operation.
 * @param outputSpace Node-id domain used to index the output.
 * @param dtype Requested NumPy floating-point type.
 * @return Computed single topology attribute.
 */
inline py::array computeSingleTopologyAttribute(std::shared_ptr<PythonValuedMorphologicalTree> tree, Attribute attribute,
                                                NodeIdSpace outputSpace = NodeIdSpace::MorphologicalTree, py::object dtype = py::none()) {
    return tree->visit([&](const auto& concreteTree) { return computeSingleTopologyAttributeImpl(concreteTree, attribute, outputSpace, dtype); });
}

/**
 * @brief Computes topology attributes from list.
 *
 * @param tree Tree topology.
 * @param attributes Attributes requested by the operation.
 * @param outputSpace Node-id domain used to index the output.
 * @param dtype Requested NumPy floating-point type.
 * @return Computed topology attributes from list.
 */
inline std::pair<py::dict, py::array> computeTopologyAttributesFromList(std::shared_ptr<PythonValuedMorphologicalTree> tree,
                                                                        const std::vector<AttributeOrGroup>& attributes,
                                                                        NodeIdSpace outputSpace = NodeIdSpace::MorphologicalTree,
                                                                        py::object dtype = py::none()) {
    return tree->visit([&](const auto& concreteTree) { return computeTopologyAttributesFromListImpl(concreteTree, attributes, outputSpace, dtype); });
}

/**
 * @brief Materializes the node-indexed buffer of an attribute requested by value or name.
 *
 * Filtering operators consume a dense attribute buffer. Accepting the attribute
 * itself removes the separate computation step from short call sites without
 * hiding a choice: the declared capability requirements state whether the
 * attribute reads node altitude, which selects the valued or the topology entry
 * point deterministically.
 *
 * @param tree Valued tree owning the nodes.
 * @param attributeRequest `Attribute` value or its symbolic name.
 * @return Dense `np.float64` attribute buffer indexed by internal node slot.
 */
[[nodiscard]] inline py::array attributeBufferFor(std::shared_ptr<PythonValuedMorphologicalTree> tree, const py::object& attributeRequest) {
    const Attribute attribute = resolveAttribute(attributeRequest);
    const py::object doubleDtype = py::dtype::of<double>();
    if (attributes::registry::capabilityRequirements(attribute).altitude) {
        return computeSingleAttribute(std::move(tree), attribute, NodeIdSpace::MorphologicalTree, doubleDtype);
    }
    return computeSingleTopologyAttribute(std::move(tree), attribute, NodeIdSpace::MorphologicalTree, doubleDtype);
}

} // namespace mmcfilters::pybindings::attribute_computation
