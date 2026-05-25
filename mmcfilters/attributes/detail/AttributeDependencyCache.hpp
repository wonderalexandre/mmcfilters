#pragma once

#include "../AttributeComputer.hpp"
#include "../AttributeNames.hpp"
#include "../AttributeResultTypes.hpp"

#include <algorithm>
#include <deque>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace mmcfilters::detail {

/**
 * @brief Non-owning view over one computed attribute result.
 *
 * @details
 * `ComputedAttributeView` is the internal cache/dependency handle used by
 * topology-backed attribute computation. It does not own either the layout or
 * the value buffer; both must outlive the computation that consumes the view.
 *
 * `nodeIdSpace` records whether the referenced rows are indexed by the dense
 * internal `MorphologicalTree` ids or by another exported/public node-id
 * convention. Dependency reuse is only valid in the internal node-id space.
 */
struct ComputedAttributeView {
    const AttributeNames* first = nullptr;
    const float* second = nullptr;
    NodeIdSpace nodeIdSpace = NodeIdSpace::MORPHOLOGICAL_TREE;

    /**
     * @brief Returns whether the view points to a valid layout/buffer pair.
     */
    [[nodiscard]] bool isValid() const noexcept {
        return first != nullptr && second != nullptr;
    }

    /**
     * @brief Returns the layout referenced by this view.
     */
    [[nodiscard]] const AttributeNames& attributeNames() const noexcept { return *first; }

    /**
     * @brief Returns the raw value buffer referenced by this view.
     */
    [[nodiscard]] const float* values() const noexcept { return second; }

    /**
     * @brief Returns the view form consumed directly by `AttributeComputer`.
     */
    [[nodiscard]] DependencySource dependencySource() const noexcept {
        return {first, second};
    }
};

/**
 * @brief Cache of already computed scalar attributes keyed by `Attribute`.
 *
 * The cache is a detail-only mechanism used to share intermediate buffers
 * materialized by `AttributePipeline` and the topology backend. Stored values
 * must remain in `NodeIdSpace::MORPHOLOGICAL_TREE` to be reusable as internal
 * dependencies.
 */
using DependencyMap = std::unordered_map<Attribute, ComputedAttributeView>;

using OwnedComputedResults = std::deque<ComputedAttributeData>;

/**
 * @brief Creates a dependency-cache view over an owning public result.
 */
[[nodiscard]] inline ComputedAttributeView makeComputedAttributeView(const ComputedAttributeData& computed) noexcept {
    return {&computed.first, computed.second.data(), computed.nodeIdSpace};
}

/**
 * @brief Checks whether a layout contains all scalar attributes in `attrs`.
 */
inline bool attributeSetContainsAll(const AttributeNames* names, const std::vector<Attribute>& attrs) {
    if (names == nullptr) {
        return false;
    }
    return std::all_of(attrs.begin(), attrs.end(), [&](Attribute attr) {
        return names->contains(attr);
    });
}

/**
 * @brief Checks whether a cached result can be reused as an internal
 * dependency.
 */
inline bool isReusableDependencyData(const ComputedAttributeView& computed, const std::vector<Attribute>& attrs) {
    return computed.nodeIdSpace == NodeIdSpace::MORPHOLOGICAL_TREE &&
           computed.isValid() &&
           attributeSetContainsAll(computed.first, attrs);
}

/**
 * @brief Registers a computed result under every scalar attribute it contains.
 */
inline void registerComputedAttributes(DependencyMap& available, const ComputedAttributeData& computed) {
    const ComputedAttributeView view = makeComputedAttributeView(computed);
    for (const auto& [attr, _] : computed.first.indexMap) {
        available[attr] = view;
    }
}

/**
 * @brief Moves one owned result into the local arena and registers its views.
 */
inline void stashComputedAttributes(OwnedComputedResults& ownedResults, DependencyMap& available, ComputedAttributeData&& computed) {
    ownedResults.emplace_back(std::move(computed));
    registerComputedAttributes(available, ownedResults.back());
}

/**
 * @brief Copies scalar attributes from one computed result into another
 * layout/buffer.
 */
inline void copyAttributesIntoBuffer(const MorphologicalTree& tree, const ComputedAttributeView& source, const std::vector<Attribute>& attrs, const AttributeNames& targetNames, float* targetBuffer) {
    if (!isReusableDependencyData(source, attrs)) {
        throw std::runtime_error("Attribute dependency buffer is not available in MorphologicalTree node-id space.");
    }

    const int numNodes = tree.getNumInternalNodeSlots();
    for (int nodeIndex = 0; nodeIndex < numNodes; ++nodeIndex) {
        for (const Attribute attr : attrs) {
            targetBuffer[targetNames.linearIndex(nodeIndex, attr)] =
                source.second[source.first->linearIndex(nodeIndex, attr)];
        }
    }
}

} // namespace mmcfilters::detail
