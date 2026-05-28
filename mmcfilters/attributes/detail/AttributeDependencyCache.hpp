#pragma once

#include "../AttributeNames.hpp"
#include "../AttributeResultTypes.hpp"
#include "AttributeKernelSupport.hpp"

#include <algorithm>
#include <concepts>
#include <deque>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace mmcfilters::detail {

/**
 * @brief Non-owning view over one computed attribute result.
 *
 * @details
 * `ComputedAttributeViewT<Real>` is the internal cache/dependency handle used by
 * topology-backed attribute computation. It does not own either the layout or
 * the value buffer; both must outlive the computation that consumes the view.
 *
 * `nodeIdSpace` records whether the referenced rows are indexed by the dense
 * internal `MorphologicalTree` ids or by another exported/public node-id
 * convention. Dependency reuse is only valid in the internal node-id space.
 */
template <std::floating_point Real = float>
struct ComputedAttributeViewT {
    const AttributeNames* first = nullptr;
    const Real* second = nullptr;
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
    [[nodiscard]] const Real* values() const noexcept { return second; }

    /**
     * @brief Returns the dependency-source view consumed by typed kernels.
     */
    [[nodiscard]] DependencySourceT<Real> dependencySource() const noexcept {
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
template <std::floating_point Real>
using DependencyMapT = std::unordered_map<Attribute, ComputedAttributeViewT<Real>>;

template <std::floating_point Real>
using OwnedComputedResultsT = std::deque<ComputedAttributeData<Real>>;

/**
 * @brief Creates a dependency-cache view over an owning public result.
 */
template <std::floating_point Real>
[[nodiscard]] inline ComputedAttributeViewT<Real> makeComputedAttributeView(const ComputedAttributeData<Real>& computed) noexcept {
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
template <std::floating_point Real>
inline bool isReusableDependencyData(const ComputedAttributeViewT<Real>& computed, const std::vector<Attribute>& attrs) {
    return computed.nodeIdSpace == NodeIdSpace::MORPHOLOGICAL_TREE &&
           computed.isValid() &&
           attributeSetContainsAll(computed.first, attrs);
}

/**
 * @brief Registers a computed result under every scalar attribute it contains.
 */
template <std::floating_point Real>
inline void registerComputedAttributes(DependencyMapT<Real>& available, const ComputedAttributeData<Real>& computed) {
    const ComputedAttributeViewT<Real> view = makeComputedAttributeView(computed);
    for (const auto& [attr, _] : computed.first.indexMap) {
        available[attr] = view;
    }
}

/**
 * @brief Moves one owned result into the local arena and registers its views.
 */
template <std::floating_point Real>
inline void stashComputedAttributes(OwnedComputedResultsT<Real>& ownedResults, DependencyMapT<Real>& available, ComputedAttributeData<Real>&& computed) {
    ownedResults.emplace_back(std::move(computed));
    registerComputedAttributes(available, ownedResults.back());
}

/**
 * @brief Copies scalar attributes from one computed result into another
 * layout/buffer.
 */
template <std::floating_point Real>
inline void copyAttributesIntoBuffer(const MorphologicalTree& tree, const ComputedAttributeViewT<Real>& source, const std::vector<Attribute>& attrs, const AttributeNames& targetNames, Real* targetBuffer) {
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
