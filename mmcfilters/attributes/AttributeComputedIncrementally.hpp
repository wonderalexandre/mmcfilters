#pragma once

#include "../utils/Common.hpp"
#include "../attributes/AttributeComputer.hpp"
#include "../attributes/AttributeNames.hpp"
#include "../trees/MorphologicalTree.hpp"
#include "../trees/TreeAltitudeOps.hpp"


namespace mmcfilters {

/**
 * @brief High-level overview of the incremental attribute design.
 *
 * @details
 * The attribute subsystem is organised around four core concepts:
 *
 * 1. `AttributeNames` / `AttributeNamesWithDelta`
 *    These classes describe how attributes are laid out inside flat buffers.
 *    They do not compute anything by themselves; they only answer questions
 *    such as "where is the value of `AREA` for node `u` stored?".
 *
 * 2. `AttributeComputer`
 *    This is the abstract interface implemented by each concrete attribute
 *    family. A computer knows:
 *    - which attributes it can produce together;
 *    - which upstream attributes it requires as dependencies;
 *    - how to fill a caller-provided flat buffer according to a given
 *      `AttributeNames` layout.
 *
 * 3. `AttributeFactory`
 *    The factory maps a public request (`Attribute` or `AttributeGroup`) to
 *    the concrete `AttributeComputer` responsible for that request. Several
 *    scalar attributes may map to the same computer when they share a common
 *    traversal or intermediate state.
 *
 * 4. `AttributeComputedIncrementally`
 *    This class is the orchestrator of the pipeline. It expands groups,
 *    resolves dependencies, orders computers according to their dependency
 *    graph, materialises intermediate results when needed, reuses compatible
 *    cached results, and finally projects the result to the requested
 *    `NodeIdSpace`.
 *
 * The intended execution model is therefore:
 * - the caller requests attributes through `AttributeComputedIncrementally`;
 * - `AttributeFactory` selects the relevant concrete computers;
 * - each `AttributeComputer` fills flat buffers described by `AttributeNames`;
 * - the orchestrator caches and reuses non-owning views of those buffers
 *   through `DependencyMap`;
 * - the final owner may then be transferred to another boundary, such as a
 *   NumPy array in the pybind layer.
 *
 * A crucial design rule is that computation always happens first in the dense
 * internal `MorphologicalTree` node-id space. Any projection to Higra or other
 * public node-id conventions happens only at the boundary of the public API.
 */
struct ComputedAttributeView {
    const AttributeNames* first = nullptr;
    const float* second = nullptr;
    NodeIdSpace nodeIdSpace = NodeIdSpace::MORPHOLOGICAL_TREE;

    /**
     * @brief Returns whether the view points to a valid layout/buffer pair.
     */
    bool isValid() const noexcept {
        return first != nullptr && second != nullptr;
    }

    /**
     * @brief Returns the layout referenced by this view.
     */
    const AttributeNames& attributeNames() const noexcept { return *first; }

    /**
     * @brief Returns the raw value buffer referenced by this view.
     */
    const float* values() const noexcept { return second; }

    /**
     * @brief Returns a dependency-compatible non-owning source.
     */
    DependencySource dependencySource() const noexcept {
        return {first, second};
    }
};

/**
 * @brief Owning result for one computed scalar attribute layout and buffer.
 *
 * @details
 * A computed attribute result is defined by three pieces of information:
 * - the attribute layout used to interpret the flat buffer;
 * - the owned buffer storing the per-node values;
 * - the `NodeIdSpace` in which the buffer indices are expressed.
 *
 * The result is move-only on purpose: the pipeline now treats one computed
 * buffer as having a single owner. Reuse inside dependency resolution happens
 * through `ComputedAttributeView`, which is a lightweight non-owning handle.
 *
 * The public fields remain named `first` and `second` so existing direct uses
 * of `.first`/`.second` continue to read naturally, while structured binding
 * is preserved through the tuple-like helpers defined below.
 */
struct ComputedAttributeData {
    AttributeNames first;
    std::vector<float> second;
    NodeIdSpace nodeIdSpace = NodeIdSpace::MORPHOLOGICAL_TREE;

    ComputedAttributeData(AttributeNames attrNames,std::vector<float> buffer,NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE)
        : first(std::move(attrNames)), second(std::move(buffer)), nodeIdSpace(outputSpace) {}

    ComputedAttributeData(const ComputedAttributeData&) = delete;
    ComputedAttributeData& operator=(const ComputedAttributeData&) = delete;
    ComputedAttributeData(ComputedAttributeData&&) noexcept = default;
    ComputedAttributeData& operator=(ComputedAttributeData&&) = delete;

    AttributeNames& attributeNames() noexcept { return this->first; }
    const AttributeNames& attributeNames() const noexcept { return this->first; }
    std::vector<float>& values() noexcept { return this->second; }
    const std::vector<float>& values() const noexcept { return this->second; }

    /**
     * @brief Returns a non-owning dependency view over this result.
     */
    ComputedAttributeView view() const noexcept { return {&this->first, this->second.data(), this->nodeIdSpace}; }
};

/**
 * @brief Owning result for one delta-augmented attribute layout and buffer.
 *
 * @details
 * This is the delta-aware counterpart of `ComputedAttributeData`. It is used
 * when one logical attribute is sampled at several ancestor/descendant offsets
 * around each node.
 */
struct ComputedAttributeDataWithDelta {
    AttributeNamesWithDelta first;
    std::vector<float> second;
    NodeIdSpace nodeIdSpace = NodeIdSpace::MORPHOLOGICAL_TREE;

    ComputedAttributeDataWithDelta(AttributeNamesWithDelta attrNames, std::vector<float> buffer, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE)
        : first(std::move(attrNames)), second(std::move(buffer)), nodeIdSpace(outputSpace) {}

    ComputedAttributeDataWithDelta(const ComputedAttributeDataWithDelta&) = delete;
    ComputedAttributeDataWithDelta& operator=(const ComputedAttributeDataWithDelta&) = delete;
    ComputedAttributeDataWithDelta(ComputedAttributeDataWithDelta&&) noexcept = default;
    ComputedAttributeDataWithDelta& operator=(ComputedAttributeDataWithDelta&&) = delete;

    AttributeNamesWithDelta& attributeNames() noexcept { return this->first; }
    const AttributeNamesWithDelta& attributeNames() const noexcept { return this->first; }
    std::vector<float>& values() noexcept { return this->second; }
    const std::vector<float>& values() const noexcept { return this->second; }
};

/**
 * @brief Cache of already computed scalar attributes keyed by `Attribute`.
 *
 * @details
 * The map stores non-owning views into previously computed attribute buffers.
 * It is used internally to avoid recomputing dependencies and may also be
 * provided by callers to seed the pipeline with precomputed results. Stored
 * values must remain in `NodeIdSpace::MORPHOLOGICAL_TREE` to be reusable as
 * internal dependencies, and the caller is responsible for keeping the owned
 * result alive while the computation runs.
 *
 * The intended seeding pattern is:
 * `auto computed = computeSingleAttribute(...);`
 * `deps[AREA] = computed.view();`
 */
using DependencyMap = std::unordered_map<Attribute, ComputedAttributeView>;

/**
 * @brief Utility entry points for incremental attribute computation on `MorphologicalTree`.
 *
 * @details
 * This class is the orchestrator of the incremental attribute framework. It
 * provides:
 * - a generic post-order traversal skeleton used by most computers;
 * - dependency resolution over `AttributeComputer` objects;
 * - cache-aware materialisation of one attribute, one group, or several
 *   heterogeneous requests at once;
 * - projection from the internal `MorphologicalTree` node-id space to other
 *   public node-id spaces such as the static Higra convention;
 * - convenience helpers that project node attributes back to the image domain.
 *
 * The canonical execution space is always the tree's dense internal node-id
 * space. Projection to Higra or other spaces only happens at the boundary of
 * the public API.
 */
class AttributeComputedIncrementally {
public:

    /**
     * @brief Applies a generic post-order computation rooted at `rootNodeId`.
     *
     * The caller provides the three phases of the incremental protocol:
     * preprocess the current node, merge each child into the parent, and
     * finalise the current node after all child merges.
     */
    template<class TreeLike, class PreProcessing, class MergeProcessing, class PostProcessing>
    static void traversePostOrder(TreeLike& tree, NodeId rootNodeId, PreProcessing&& preProcessing, MergeProcessing&& mergeProcessing, PostProcessing&& postProcessing) {
        preProcessing(rootNodeId);
        for (NodeId childNodeId : tree.getChildren(rootNodeId)) {
            AttributeComputedIncrementally::traversePostOrder(tree, childNodeId, preProcessing, mergeProcessing, postProcessing);
            mergeProcessing(rootNodeId, childNodeId);
        }
        postProcessing(rootNodeId);
    }

    /**
     * @brief Computes exactly the natural output of a concrete
     * `AttributeComputer`.
     *
     * @details
     * This is the lowest-level public entry point that still goes through the
     * dependency resolver. It is mainly useful when the caller already knows
     * which concrete computer they want to run.
     */
    static ComputedAttributeData computeAttributesByComputer(MorphologicalTree& tree, const AltitudeBuffer* altitude, const AttributeComputer& comp, const DependencyMap& available = {}, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE);
    static ComputedAttributeData computeAttributesByComputer(MorphologicalTree& tree, const AttributeComputer& comp, const DependencyMap& available = {}, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE);
    static ComputedAttributeData computeAttributesByComputer(WeightedMorphologicalTree& tree, const AttributeComputer& comp, const DependencyMap& available = {}, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE);

    /**
     * @brief Computes a single scalar attribute or a full attribute group.
     *
     * @details
     * If the request is a scalar attribute, the returned layout contains only
     * that attribute. If the request is a group, the returned layout contains
     * the full group expansion.
     */
    static ComputedAttributeData computeSingleAttribute(MorphologicalTree& tree, const AltitudeBuffer* altitude, AttributeOrGroup attr, const DependencyMap& availableDeps = {}, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE);
    static ComputedAttributeData computeSingleAttribute(MorphologicalTree& tree, AttributeOrGroup attr, const DependencyMap& availableDeps = {}, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE);
    static ComputedAttributeData computeSingleAttribute(WeightedMorphologicalTree& tree, AttributeOrGroup attr, const DependencyMap& availableDeps = {}, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE);

    /**
     * @brief Computes a delta-augmented version of one scalar attribute.
     *
     * @param delta Maximum ancestor/descendant offset to materialise.
     * @param padding Strategy used when the requested offset leaves the valid
     * ancestor/descendant chain.
     */
    static ComputedAttributeDataWithDelta computeSingleAttributeWithDelta(MorphologicalTree& tree, const AltitudeBuffer* altitude, Attribute attribute, int delta, std::string padding="last-padding", const DependencyMap& availableDeps={}, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE);
    static ComputedAttributeDataWithDelta computeSingleAttributeWithDelta(MorphologicalTree& tree, Attribute attribute, int delta, std::string padding="last-padding", const DependencyMap& availableDeps={}, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE);
    static ComputedAttributeDataWithDelta computeSingleAttributeWithDelta(WeightedMorphologicalTree& tree, Attribute attribute, int delta, std::string padding="last-padding", const DependencyMap& availableDeps={}, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE);

    /**
     * @brief Computes a heterogeneous set of scalar attributes and attribute
     * groups in one coordinated run.
     *
     * @details
     * The method expands groups, derives an order of concrete computers that
     * respects declared dependencies, reuses cached dependencies when possible,
     * and returns a single layout containing exactly the requested scalar
     * attributes.
     */
    static ComputedAttributeData computeAttributes(MorphologicalTree& tree, const AltitudeBuffer* altitude, const std::vector<AttributeOrGroup>& attributes,const DependencyMap& providedDependencies={}, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE);
    static ComputedAttributeData computeAttributes(MorphologicalTree& tree, const std::vector<AttributeOrGroup>& attributes,const DependencyMap& providedDependencies={}, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE);
    static ComputedAttributeData computeAttributes(WeightedMorphologicalTree& tree, const std::vector<AttributeOrGroup>& attributes,const DependencyMap& providedDependencies={}, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE);

    /**
     * @brief Projects a node attribute to a proper-part image in the original domain.
     */
    static ImageFloatPtr computeAttributeMapping(MorphologicalTree& tree, const AltitudeBuffer* altitude, Attribute attribute);
    static ImageFloatPtr computeAttributeMapping(MorphologicalTree& tree, Attribute attribute);
    static ImageFloatPtr computeAttributeMapping(WeightedMorphologicalTree& tree, Attribute attribute);

};

} // namespace mmcfilters

namespace std {

template <>
struct tuple_size<mmcfilters::ComputedAttributeData> : integral_constant<std::size_t, 2> {};

template <>
struct tuple_element<0, mmcfilters::ComputedAttributeData> {
    using type = mmcfilters::AttributeNames;
};

template <>
struct tuple_element<1, mmcfilters::ComputedAttributeData> {
    using type = std::vector<float>;
};

template <>
struct tuple_size<mmcfilters::ComputedAttributeDataWithDelta> : integral_constant<std::size_t, 2> {};

template <>
struct tuple_element<0, mmcfilters::ComputedAttributeDataWithDelta> {
    using type = mmcfilters::AttributeNamesWithDelta;
};

template <>
struct tuple_element<1, mmcfilters::ComputedAttributeDataWithDelta> {
    using type = std::vector<float>;
};

} // namespace std

namespace mmcfilters {

/**
 * @brief `std::get` compatibility overloads preserving tuple-like access to
 * `ComputedAttributeData`.
 */
template <std::size_t I>
decltype(auto) get(ComputedAttributeData& computed) noexcept {
    if constexpr (I == 0) {
        return (computed.first);
    } else {
        return (computed.second);
    }
}

template <std::size_t I>
decltype(auto) get(const ComputedAttributeData& computed) noexcept {
    if constexpr (I == 0) {
        return (computed.first);
    } else {
        return (computed.second);
    }
}

template <std::size_t I>
decltype(auto) get(ComputedAttributeData&& computed) noexcept {
    if constexpr (I == 0) {
        return std::move(computed.first);
    } else {
        return std::move(computed.second);
    }
}

template <std::size_t I>
decltype(auto) get(ComputedAttributeDataWithDelta& computed) noexcept {
    if constexpr (I == 0) {
        return (computed.first);
    } else {
        return (computed.second);
    }
}

template <std::size_t I>
decltype(auto) get(const ComputedAttributeDataWithDelta& computed) noexcept {
    if constexpr (I == 0) {
        return (computed.first);
    } else {
        return (computed.second);
    }
}

template <std::size_t I>
decltype(auto) get(ComputedAttributeDataWithDelta&& computed) noexcept {
    if constexpr (I == 0) {
        return std::move(computed.first);
    } else {
        return std::move(computed.second);
    }
}

} // namespace mmcfilters



#include "../attributes/AttributeFactory.hpp"
#include "../attributes/detail/AttributeComputedIncrementallyDetail.hpp"



namespace mmcfilters {

/**
 * @brief Executes one concrete computer with dependency resolution and optional
 * result projection.
 */
inline ComputedAttributeData AttributeComputedIncrementally::computeAttributesByComputer(MorphologicalTree& tree, const AltitudeBuffer* altitude, const AttributeComputer& comp, const DependencyMap& availableDeps, NodeIdSpace outputSpace) {
    DependencyMap available = availableDeps;
    detail::OwnedComputedResults ownedResults;

    for (const AttributeOrGroup& dep : comp.requiredAttributes()) {
        bool needsCompute = false;

        if (std::holds_alternative<Attribute>(dep)) {
            const Attribute attr = std::get<Attribute>(dep);
            needsCompute = !available.count(attr) || !detail::isReusableDependencyData(available.at(attr), {attr});
        } else {
            const auto& groupAttrs = detail::attributesOf(std::get<AttributeGroup>(dep));
            const Attribute representativeAttr = groupAttrs.front();
            needsCompute = !available.count(representativeAttr) ||
                           !detail::isReusableDependencyData(available.at(representativeAttr), groupAttrs);
        }

        if (needsCompute) {
            auto computed = computeSingleAttribute(tree, altitude, dep, available);
            detail::stashComputedAttributes(ownedResults, available, std::move(computed));
        }
    }

    auto dependencySources = detail::resolveDependencySources(tree, altitude, comp, available, ownedResults);
    const auto& computedAttrs = comp.attributes();

    std::unordered_map<Attribute, int> attrOffsets;
    for (int i = 0; i < static_cast<int>(computedAttrs.size()); ++i) {
        attrOffsets[computedAttrs[i]] = i;
    }
    AttributeNames attrNames(std::move(attrOffsets));

    int n = tree.getNumInternalNodeSlots();
    std::vector<float> buffer(static_cast<size_t>(n) * attrNames.NUM_ATTRIBUTES, 0.0f);

    comp.compute(tree, altitude, buffer, attrNames, dependencySources);

    return detail::projectComputedDataToNodeIdSpace(tree, {std::move(attrNames), std::move(buffer), NodeIdSpace::MORPHOLOGICAL_TREE}, outputSpace);
}

inline ComputedAttributeData AttributeComputedIncrementally::computeAttributesByComputer(MorphologicalTree& tree, const AttributeComputer& comp, const DependencyMap& availableDeps, NodeIdSpace outputSpace) {
    return computeAttributesByComputer(tree, nullptr, comp, availableDeps, outputSpace);
}

inline ComputedAttributeData AttributeComputedIncrementally::computeAttributesByComputer(WeightedMorphologicalTree& tree, const AttributeComputer& comp, const DependencyMap& availableDeps, NodeIdSpace outputSpace) {
    return computeAttributesByComputer(tree.tree, &tree.altitude, comp, availableDeps, outputSpace);
}

/**
 * @brief Computes one scalar attribute or one public attribute group.
 *
 * @details
 * The function is recursive through the dependency resolver: if the target
 * computer requires other attributes, they are materialised first and exposed
 * as dependency views. The returned buffer contains only the requested scalar
 * attributes.
 */
inline ComputedAttributeData AttributeComputedIncrementally::computeSingleAttribute(MorphologicalTree& tree, const AltitudeBuffer* altitude, AttributeOrGroup attrOrGroup, const DependencyMap& availableDeps, NodeIdSpace outputSpace) {
    const AttributeComputer& comp = AttributeFactory::create(attrOrGroup);
    DependencyMap available = availableDeps;
    detail::OwnedComputedResults ownedResults;

    std::vector<Attribute> requestedAttrs;
    if (std::holds_alternative<Attribute>(attrOrGroup)) {
        requestedAttrs.push_back(std::get<Attribute>(attrOrGroup));
    } else {
        requestedAttrs = detail::attributesOf(std::get<AttributeGroup>(attrOrGroup));
    }

    for (const AttributeOrGroup& dep : comp.requiredAttributes()) {
        bool needsCompute = false;

        if (std::holds_alternative<Attribute>(dep)) {
            const Attribute attr = std::get<Attribute>(dep);
            needsCompute = !available.count(attr) || !detail::isReusableDependencyData(available.at(attr), {attr});
        } else {
            const auto& groupAttrs = detail::attributesOf(std::get<AttributeGroup>(dep));
            const Attribute representativeAttr = groupAttrs.front();
            needsCompute = !available.count(representativeAttr) ||
                           !detail::isReusableDependencyData(available.at(representativeAttr), groupAttrs);
        }

        if (needsCompute) {
            auto computed = computeSingleAttribute(tree, altitude, dep, available);
            detail::stashComputedAttributes(ownedResults, available, std::move(computed));
        }
    }

    auto dependencySources = detail::resolveDependencySources(tree, altitude, comp, available, ownedResults);

    std::unordered_map<Attribute, int> attrOffsets;
    for (int i = 0; i < static_cast<int>(requestedAttrs.size()); ++i) {
        attrOffsets[requestedAttrs[i]] = i;
    }
    AttributeNames attrNames(std::move(attrOffsets));

    int n = tree.getNumInternalNodeSlots();
    std::vector<float> buffer(static_cast<size_t>(n) * attrNames.NUM_ATTRIBUTES, 0.0f);

    comp.compute(tree, altitude, buffer, attrNames, requestedAttrs, dependencySources);

    return detail::projectComputedDataToNodeIdSpace(tree, {std::move(attrNames), std::move(buffer), NodeIdSpace::MORPHOLOGICAL_TREE}, outputSpace);
}

inline ComputedAttributeData AttributeComputedIncrementally::computeSingleAttribute(MorphologicalTree& tree, AttributeOrGroup attrOrGroup, const DependencyMap& availableDeps, NodeIdSpace outputSpace) {
    return computeSingleAttribute(tree, nullptr, attrOrGroup, availableDeps, outputSpace);
}

inline ComputedAttributeData AttributeComputedIncrementally::computeSingleAttribute(WeightedMorphologicalTree& tree, AttributeOrGroup attrOrGroup, const DependencyMap& availableDeps, NodeIdSpace outputSpace) {
    return computeSingleAttribute(tree.tree, &tree.altitude, attrOrGroup, availableDeps, outputSpace);
}

/**
 * @brief Computes one scalar attribute augmented with ancestor/descendant
 * delta samples.
 *
 * @details
 * The base attribute is first computed in the canonical internal node-id
 * space. The method then builds a delta-aware layout storing:
 * - the node value at delta `0`;
 * - values sampled at ancestors for negative deltas;
 * - values sampled at descendants for positive deltas.
 *
 * Missing samples are handled according to the selected padding strategy.
 */
inline ComputedAttributeDataWithDelta AttributeComputedIncrementally::computeSingleAttributeWithDelta(MorphologicalTree& tree, const AltitudeBuffer* altitude, Attribute attribute, int delta, std::string padding, const DependencyMap& availableDeps, NodeIdSpace outputSpace) {
    if (padding != "last-padding" &&
        padding != "nan-padding" &&
        padding != "null-padding" &&
        padding != "zero-padding") {
        throw std::invalid_argument("Unknown padding strategy.");
    }

    auto [attributeNamesBase, attrsBase] = AttributeComputedIncrementally::computeSingleAttribute(tree, altitude, attribute, availableDeps);

    int n = tree.getNumInternalNodeSlots();
    std::vector<Attribute> attrVec = {attribute};
    AttributeNamesWithDelta attributeNamesDelta(AttributeNamesWithDelta::create(delta, attrVec));

    std::vector<float> attrsDelta(static_cast<size_t>(n) * attributeNamesDelta.NUM_ATTRIBUTES);
    std::fill_n(
        attrsDelta.data(),
        static_cast<size_t>(n) * attributeNamesDelta.NUM_ATTRIBUTES,
        std::numeric_limits<float>::quiet_NaN()
    );

    for (NodeId nodeId : tree.getAliveNodeIds()) {
        const NodeId nodeIndex = nodeId;
        int outIdx = attributeNamesDelta.linearIndex(nodeIndex, attribute, 0);
        int baseIdx = attributeNamesBase.linearIndex(nodeIndex, attribute);
        attrsDelta[outIdx] = attrsBase[baseIdx];
    }

    for (int d = 1; d <= delta; ++d) {
        auto [ascendants, descendants] = tree_altitude_ops::computeAscendantsAndDescendants(tree, altitude, d);

        for (NodeId nodeId : tree.getAliveNodeIds()) {
            const NodeId nodeIndex = nodeId;

            int ascIndex = (ascendants[nodeIndex] != InvalidNode ? ascendants[nodeIndex] : nodeIndex);
            if (ascIndex != nodeIndex) {
                int outIdxAsc = attributeNamesDelta.linearIndex(nodeIndex, attribute, -d);
                int baseIdxAsc = attributeNamesBase.linearIndex(ascIndex, attribute);
                attrsDelta[outIdxAsc] = attrsBase[baseIdxAsc];
            }

            int descIndex = (descendants[nodeIndex] != InvalidNode ? descendants[nodeIndex] : nodeIndex);
            if (descIndex != nodeIndex) {
                int outIdxDesc = attributeNamesDelta.linearIndex(nodeIndex, attribute, +d);
                int baseIdxDesc = attributeNamesBase.linearIndex(descIndex, attribute);
                attrsDelta[outIdxDesc] = attrsBase[baseIdxDesc];
            }
        }
    }

    if (padding == "last-padding" || padding == "nan-padding") {
        for (NodeId nodeId : tree.getAliveNodeIds()) {
            const NodeId nodeIndex = nodeId;

            for (int d = 1; d <= delta; ++d) {
                int outIdx = attributeNamesDelta.linearIndex(nodeIndex, attribute, -d);
                int refIdx = attributeNamesDelta.linearIndex(nodeIndex, attribute, -(d - 1));

                if (std::isnan(attrsDelta[outIdx])) {
                    if (padding == "last-padding") {
                        attrsDelta[outIdx] = attrsDelta[refIdx];
                    } else {
                        attrsDelta[outIdx] = std::numeric_limits<float>::quiet_NaN();
                    }
                }
            }

            for (int d = 1; d <= delta; ++d) {
                int outIdx = attributeNamesDelta.linearIndex(nodeIndex, attribute, +d);
                int refIdx = attributeNamesDelta.linearIndex(nodeIndex, attribute, +(d - 1));

                if (tree.isLeaf(nodeIndex) || std::isnan(attrsDelta[outIdx])) {
                    if (padding == "last-padding") {
                        attrsDelta[outIdx] = attrsDelta[refIdx];
                    } else {
                        attrsDelta[outIdx] = std::numeric_limits<float>::quiet_NaN();
                    }
                }
            }
        }
    } else if (padding == "zero-padding") {
        for (NodeId nodeId : tree.getAliveNodeIds()) {
            const NodeId nodeIndex = nodeId;
            for (int d = 1; d <= delta; ++d) {
                const int ascIdx = attributeNamesDelta.linearIndex(nodeIndex, attribute, -d);
                const int descIdx = attributeNamesDelta.linearIndex(nodeIndex, attribute, +d);
                if (std::isnan(attrsDelta[ascIdx])) {
                    attrsDelta[ascIdx] = 0.0f;
                }
                if (std::isnan(attrsDelta[descIdx])) {
                    attrsDelta[descIdx] = 0.0f;
                }
            }
        }
    }

    return detail::projectComputedDataToNodeIdSpace(tree, {std::move(attributeNamesDelta), std::move(attrsDelta), NodeIdSpace::MORPHOLOGICAL_TREE}, outputSpace);
}

inline ComputedAttributeDataWithDelta AttributeComputedIncrementally::computeSingleAttributeWithDelta(MorphologicalTree& tree, Attribute attribute, int delta, std::string padding, const DependencyMap& availableDeps, NodeIdSpace outputSpace) {
    return computeSingleAttributeWithDelta(tree, nullptr, attribute, delta, std::move(padding), availableDeps, outputSpace);
}

inline ComputedAttributeDataWithDelta AttributeComputedIncrementally::computeSingleAttributeWithDelta(WeightedMorphologicalTree& tree, Attribute attribute, int delta, std::string padding, const DependencyMap& availableDeps, NodeIdSpace outputSpace) {
    return computeSingleAttributeWithDelta(tree.tree, &tree.altitude, attribute, delta, std::move(padding), availableDeps, outputSpace);
}

/**
 * @brief Computes a heterogeneous set of scalar attributes and groups in one
 * coordinated execution.
 *
 * @details
 * The method expands the user request into scalar attributes, derives a
 * dependency-respecting execution order over concrete computers, reuses cached
 * dependencies when safe, and finally assembles one return buffer containing
 * exactly the requested scalar attributes.
 */
inline ComputedAttributeData AttributeComputedIncrementally::computeAttributes(MorphologicalTree& tree, const AltitudeBuffer* altitude, const std::vector<AttributeOrGroup>& attributes, const DependencyMap& providedDependencies, NodeIdSpace outputSpace) {
    DependencyMap available = providedDependencies;
    detail::OwnedComputedResults ownedResults;

    std::set<Attribute> uniqueExpandedAttrs;
    std::unordered_map<std::type_index, std::vector<Attribute>> attributesPerComputer;

    for (const auto& item : attributes) {
        std::vector<Attribute> attrs;
        if (std::holds_alternative<Attribute>(item)) {
            attrs.push_back(std::get<Attribute>(item));
        } else {
            auto groupAttrs = detail::attributesOf(std::get<AttributeGroup>(item));
            attrs.insert(attrs.end(), groupAttrs.begin(), groupAttrs.end());
        }

        for (const auto& attr : attrs) {
            uniqueExpandedAttrs.insert(attr);
            const AttributeComputer& comp = AttributeFactory::create(attr);
            std::type_index id(typeid(comp));
            attributesPerComputer[id].push_back(attr);
        }
    }

    std::unordered_map<Attribute, int> attrOffsets;
    int offset = 0;
    for (const auto& attr : uniqueExpandedAttrs) {
        attrOffsets[attr] = offset++;
    }

    AttributeNames attrNames(std::move(attrOffsets));
    int n = tree.getNumInternalNodeSlots();
    std::vector<float> buffer(static_cast<size_t>(n) * attrNames.NUM_ATTRIBUTES, 0.0f);

    auto orderedComputers = detail::getOrderedComputers(attributes);

    for (const AttributeComputer* comp : orderedComputers) {
        const auto& ref = *comp;
        std::type_index id(typeid(ref));

        if (!attributesPerComputer.count(id)) {
            auto computed = computeAttributesByComputer(tree, altitude, ref, available);
            detail::stashComputedAttributes(ownedResults, available, std::move(computed));
            continue;
        }

        std::vector<Attribute> userRequestedAttrs  = attributesPerComputer.at(id);
        bool alreadyAvailable = std::all_of(userRequestedAttrs.begin(), userRequestedAttrs.end(), [&](const Attribute& a) {
            return available.count(a) && detail::isReusableDependencyData(available.at(a), {a});
        });
        if (alreadyAvailable) {
            continue;
        }

        auto depsForThis = detail::resolveDependencySources(tree, altitude, ref, available, ownedResults);
        ref.compute(tree, altitude, buffer, attrNames, userRequestedAttrs, depsForThis);

        for (const auto& attr : userRequestedAttrs) {
            available[attr] = ComputedAttributeView{&attrNames, buffer.data(), NodeIdSpace::MORPHOLOGICAL_TREE};
        }
    }

    for (const auto& attr : uniqueExpandedAttrs) {
        const auto it = available.find(attr);
        if (it == available.end()) {
            throw std::runtime_error("Requested attribute was not materialised.");
        }
        detail::copyAttributesIntoBuffer(tree, it->second, {attr}, attrNames, buffer.data());
    }

    return detail::projectComputedDataToNodeIdSpace(tree, {std::move(attrNames), std::move(buffer), NodeIdSpace::MORPHOLOGICAL_TREE}, outputSpace);
}

inline ComputedAttributeData AttributeComputedIncrementally::computeAttributes(MorphologicalTree& tree, const std::vector<AttributeOrGroup>& attributes, const DependencyMap& providedDependencies, NodeIdSpace outputSpace) {
    return computeAttributes(tree, nullptr, attributes, providedDependencies, outputSpace);
}

inline ComputedAttributeData AttributeComputedIncrementally::computeAttributes(WeightedMorphologicalTree& tree, const std::vector<AttributeOrGroup>& attributes, const DependencyMap& providedDependencies, NodeIdSpace outputSpace) {
    return computeAttributes(tree.tree, &tree.altitude, attributes, providedDependencies, outputSpace);
}

/**
 * @brief Projects a scalar node attribute back to the image domain.
 *
 * @details
 * Each proper part receives the value of the smallest component that contains
 * it, producing an image whose pixels are indexed in the original spatial
 * domain rather than in the tree node domain.
 */
inline ImageFloatPtr AttributeComputedIncrementally::computeAttributeMapping(MorphologicalTree& tree, const AltitudeBuffer* altitude, Attribute attribute) {
    auto [attrNames, buffer] = AttributeComputedIncrementally::computeSingleAttribute(tree, altitude, attribute);
    ImageFloatPtr imgPtr = std::make_shared<ImageFloat>(tree.getNumRowsOfImage(), tree.getNumColsOfImage());
    float* img = imgPtr->rawData();
    for (int p = 0; p < imgPtr->getSize(); ++p) {
        const NodeId nodeId = tree.getSmallestComponent(p);
        const NodeId index = nodeId;
        img[p] = buffer[attrNames.linearIndex(index, attribute)];
    }
    return imgPtr;
}

inline ImageFloatPtr AttributeComputedIncrementally::computeAttributeMapping(MorphologicalTree& tree, Attribute attribute) {
    return computeAttributeMapping(tree, nullptr, attribute);
}

inline ImageFloatPtr AttributeComputedIncrementally::computeAttributeMapping(WeightedMorphologicalTree& tree, Attribute attribute) {
    return computeAttributeMapping(tree.tree, &tree.altitude, attribute);
}

} // namespace mmcfilters
