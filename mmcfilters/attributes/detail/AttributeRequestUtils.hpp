#pragma once

#include "../AttributeNames.hpp"
#include "../../trees/MorphologicalTree.hpp"

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

namespace mmcfilters::detail {

/**
 * @brief Expands a public attribute group into its scalar members.
 *
 * @param group Public attribute group.
 * @return Values produced by the operation.
 */
inline std::vector<Attribute> attributesOf(AttributeGroup group) {
    const auto it = ATTRIBUTE_GROUPS.find(group);
    if (it != ATTRIBUTE_GROUPS.end()) {
        return it->second;
    }
    throw std::runtime_error("Unknown AttributeGroup.");
}

/**
 * @brief Expands one scalar-or-group request into scalar attributes.
 *
 * @param item Attribute or group request item.
 * @return Values produced by the operation.
 */
inline std::vector<Attribute> expandAttributeOrGroup(const AttributeOrGroup& item) {
    if (std::holds_alternative<Attribute>(item)) {
        return {std::get<Attribute>(item)};
    }
    return attributesOf(std::get<AttributeGroup>(item));
}

/**
 * @brief Returns whether every scalar attribute implied by a request is accepted.
 *
 * @param attributes Attributes requested by the operation.
 * @param acceptsAttribute Attribute information represented by `acceptsAttribute`.
 * @return Whether every scalar attribute implied by a request is accepted.
 */
template <class Predicate> inline bool requestContainsOnlyAttributes(const std::vector<AttributeOrGroup>& attributes, Predicate&& acceptsAttribute) {
    for (const AttributeOrGroup& item : attributes) {
        for (const Attribute attribute : expandAttributeOrGroup(item)) {
            if (!acceptsAttribute(attribute)) {
                return false;
            }
        }
    }
    return true;
}

/**
 * @brief Expands a request into a sorted unique scalar attribute list.
 *
 * @param attributes Attributes requested by the operation.
 * @param acceptsAttribute Attribute information represented by `acceptsAttribute`.
 * @param unsupportedAttributeMessage Attribute information represented by `unsupportedAttributeMessage`.
 * @param unsupportedGroupMessage Diagnostic message for unsupported groups.
 * @return Values produced by the operation.
 */
template <class Predicate>
inline std::vector<Attribute> expandUniqueAttributeRequest(const std::vector<AttributeOrGroup>& attributes, Predicate&& acceptsAttribute,
                                                           const char* unsupportedAttributeMessage, const char* unsupportedGroupMessage) {
    std::set<Attribute> expanded;
    for (const AttributeOrGroup& item : attributes) {
        const bool isScalar = std::holds_alternative<Attribute>(item);
        for (const Attribute attribute : expandAttributeOrGroup(item)) {
            if (!acceptsAttribute(attribute)) {
                throw std::invalid_argument(isScalar ? unsupportedAttributeMessage : unsupportedGroupMessage);
            }
            expanded.insert(attribute);
        }
    }
    return {expanded.begin(), expanded.end()};
}

/**
 * @brief Tests whether a scalar attribute list contains `target`.
 *
 * @param attributes Attributes requested by the operation.
 * @param target Destination object or value.
 * @return True if a scalar attribute list contains target; otherwise false.
 */
inline bool containsAttribute(std::span<const Attribute> attributes, Attribute target) {
    return std::find(attributes.begin(), attributes.end(), target) != attributes.end();
}

/**
 * @brief Builds a dense node-major attribute buffer initialized to zero.
 *
 * @param tree Tree topology used by the operation.
 * @param attrNames Layout that maps attributes to buffer columns.
 * @return The resulting dense node-major attribute buffer initialized to zero.
 */
template <std::floating_point Real = float> inline std::vector<Real> makeAttributeValueBuffer(const MorphologicalTree& tree, const AttributeNames& attrNames) {
    return std::vector<Real>(static_cast<std::size_t>(tree.getNumInternalNodeSlots()) * static_cast<std::size_t>(attrNames.NUM_ATTRIBUTES), Real{0});
}

/**
 * @brief Builds a dense explicit-offset layout for a sorted attribute set.
 *
 * @param attributes Attributes requested by the operation.
 * @return The resulting dense explicit-offset layout for a sorted attribute set.
 */
inline AttributeNames makeAttributeNamesFromSet(const std::set<Attribute>& attributes) {
    return AttributeNames::fromList({attributes.begin(), attributes.end()});
}

/**
 * @brief Copies one scalar attribute between two flat node-major layouts.
 *
 * @param tree Tree topology used by the operation.
 * @param sourceNames Input represented by `sourceNames`.
 * @param sourceBuffer Input represented by `sourceBuffer`.
 * @param targetNames Destination represented by `targetNames`.
 * @param targetBuffer Destination represented by `targetBuffer`.
 * @param attribute Attribute requested by the operation.
 */
template <std::floating_point Real>
inline void copyAttributeValuesBetweenLayouts(const MorphologicalTree& tree, const AttributeNames& sourceNames, std::span<const Real> sourceBuffer,
                                              const AttributeNames& targetNames, std::span<Real> targetBuffer, Attribute attribute) {
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        targetBuffer[targetNames.linearIndex(nodeId, attribute)] = sourceBuffer[sourceNames.linearIndex(nodeId, attribute)];
    }
}

} // namespace mmcfilters::detail
