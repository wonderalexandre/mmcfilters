#pragma once

#include <deque>

namespace mmcfilters::detail {

using OwnedComputedResults = std::deque<ComputedAttributeData>;

/**
 * @brief Expands a public attribute group into its scalar members.
 */
inline std::vector<Attribute> attributesOf(AttributeGroup group) {
    const auto it = ATTRIBUTE_GROUPS.find(group);
    if (it != ATTRIBUTE_GROUPS.end()) {
        return it->second;
    }
    throw std::runtime_error("Unknown AttributeGroup.");
}

/**
 * @brief Checks whether a layout contains all scalar attributes in `attrs`.
 */
inline bool attributeSetContainsAll(const AttributeNames* names, const std::vector<Attribute>& attrs) {
    if (names == nullptr) {
        return false;
    }
    return std::all_of(attrs.begin(), attrs.end(), [&](Attribute attr) {
        return names->indexMap.count(attr) > 0;
    });
}

/**
 * @brief Checks whether a cached result can be reused as an internal
 * dependency.
 *
 * @details
 * Internal reuse is only safe when the cached data:
 * - exists;
 * - exposes the required scalar attributes;
 * - is expressed in the tree's dense internal node-id space.
 */
inline bool isReusableDependencyData(
    const ComputedAttributeView& computed,
    const std::vector<Attribute>& attrs
) {
    return computed.nodeIdSpace == NodeIdSpace::MORPHOLOGICAL_TREE &&
           computed.isValid() &&
           attributeSetContainsAll(computed.first, attrs);
}

/**
 * @brief Registers a computed result under every scalar attribute it contains.
 */
inline void registerComputedAttributes(
    DependencyMap& available,
    const ComputedAttributeData& computed
) {
    const ComputedAttributeView view = computed.view();
    for (const auto& [attr, _] : computed.first.indexMap) {
        available[attr] = view;
    }
}

/**
 * @brief Moves one owned result into the local arena and registers its views.
 */
inline void stashComputedAttributes(
    OwnedComputedResults& ownedResults,
    DependencyMap& available,
    ComputedAttributeData&& computed
) {
    ownedResults.emplace_back(std::move(computed));
    registerComputedAttributes(available, ownedResults.back());
}

/**
 * @brief Copies scalar attributes from one computed result into another
 * layout/buffer.
 *
 * @details
 * This helper is used when the caller requested a heterogeneous attribute set
 * and some requested values were already available in the dependency cache.
 * The values must still be copied into the final return buffer so that the
 * result is complete.
 */
inline void copyAttributesIntoBuffer(
    const MorphologicalTree& tree,
    const ComputedAttributeView& source,
    const std::vector<Attribute>& attrs,
    const AttributeNames& targetNames,
    float* targetBuffer
) {
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

/**
 * @brief Resolves the dependency views required by one concrete computer.
 *
 * @details
 * Missing or non-reusable dependencies are materialised recursively on demand,
 * then re-registered in the dependency cache so downstream computers can reuse
 * them without recomputation.
 */
inline std::vector<DependencySource> resolveDependencySources(
    MorphologicalTree& tree,
    const AltitudeBuffer* altitude,
    const AttributeComputer& comp,
    DependencyMap& available,
    OwnedComputedResults& ownedResults
) {
    std::vector<DependencySource> dependencySources;

    for (const AttributeOrGroup& dep : comp.requiredAttributes()) {
        if (std::holds_alternative<Attribute>(dep)) {
            const Attribute attr = std::get<Attribute>(dep);
            if (!available.count(attr) || !isReusableDependencyData(available.at(attr), {attr})) {
                auto computed = AttributeComputedIncrementally::computeSingleAttribute(tree, altitude, dep, available);
                stashComputedAttributes(ownedResults, available, std::move(computed));
            }
            dependencySources.push_back(available.at(attr).dependencySource());
            continue;
        }

        const auto groupAttrs = attributesOf(std::get<AttributeGroup>(dep));
        const Attribute representativeAttr = groupAttrs.front();
        const bool canReuseExistingGroup =
            available.count(representativeAttr) &&
            isReusableDependencyData(available.at(representativeAttr), groupAttrs);

        if (!canReuseExistingGroup) {
            auto computed = AttributeComputedIncrementally::computeSingleAttribute(tree, altitude, dep, available);
            stashComputedAttributes(ownedResults, available, std::move(computed));
        }

        dependencySources.push_back(available.at(representativeAttr).dependencySource());
    }

    return dependencySources;
}

/**
 * @brief Builds a topological execution order for the concrete computers
 * needed by a heterogeneous request.
 *
 * @details
 * The algorithm first discovers all scalar attributes implied by the request
 * and their transitive dependencies, grouping them by concrete computer type.
 * It then performs a depth-first traversal on the induced dependency graph of
 * computer types to obtain an order in which every computer appears after the
 * computers that provide its prerequisites.
 *
 * The returned collection stores non-owning pointers to stateless singleton
 * computers owned by `AttributeFactory`.
 */
inline std::vector<const AttributeComputer*> getOrderedComputers(const std::vector<AttributeOrGroup>& attrOrGroups) {
    using ACptr = const AttributeComputer*;
    using TIndex = std::type_index;

    std::map<TIndex, ACptr> computerMap;
    std::map<TIndex, std::set<TIndex>> dependencyGraph;
    std::set<Attribute> visitedAttrs;

    std::function<void(Attribute)> collect = [&](Attribute attr) {
        if (visitedAttrs.count(attr)) {
            return;
        }
        visitedAttrs.insert(attr);

        const AttributeComputer& comp = AttributeFactory::create(attr);
        TIndex id(typeid(comp));
        computerMap[id] = &comp;
        dependencyGraph[id];

        for (const auto& depOrGroup : comp.requiredAttributes()) {
            std::vector<Attribute> deps;
            if (std::holds_alternative<Attribute>(depOrGroup)) {
                deps.push_back(std::get<Attribute>(depOrGroup));
            } else {
                auto groupAttrs = attributesOf(std::get<AttributeGroup>(depOrGroup));
                deps.insert(deps.end(), groupAttrs.begin(), groupAttrs.end());
            }

            for (const auto& depAttr : deps) {
                collect(depAttr);
                const AttributeComputer& depComp = AttributeFactory::create(depAttr);
                std::type_index depId(typeid(depComp));
                dependencyGraph[id].insert(depId);
            }
        }
    };

    for (const auto& item : attrOrGroups) {
        std::vector<Attribute> attrs;
        if (std::holds_alternative<Attribute>(item)) {
            attrs.push_back(std::get<Attribute>(item));
        } else {
            auto groupAttrs = attributesOf(std::get<AttributeGroup>(item));
            attrs.insert(attrs.end(), groupAttrs.begin(), groupAttrs.end());
        }

        for (const auto& attr : attrs) {
            collect(attr);
        }
    }

    std::vector<ACptr> ordered;
    std::set<TIndex> visited;

    std::function<void(TIndex)> visit = [&](TIndex id) {
        if (visited.count(id)) {
            return;
        }
        visited.insert(id);
        for (const auto& depId : dependencyGraph[id]) {
            visit(depId);
        }
        ordered.push_back(computerMap.at(id));
    };

    for (const auto& [id, _] : computerMap) {
        visit(id);
    }

    return ordered;
}

/**
 * @brief Projects a scalar result from the internal node-id space to another
 * public node-id space.
 *
 * @details
 * The internal pipeline always computes in `MORPHOLOGICAL_TREE` space. This
 * projection step is only applied at the API boundary, for example to expose
 * results in the static Higra convention. Nodes that do not exist in the
 * target space are left as `NaN`.
 */
inline ComputedAttributeData projectComputedDataToNodeIdSpace(const MorphologicalTree& tree, ComputedAttributeData computed, NodeIdSpace outputSpace) {
    if (outputSpace == NodeIdSpace::MORPHOLOGICAL_TREE) {
        computed.nodeIdSpace = NodeIdSpace::MORPHOLOGICAL_TREE;
        return computed;
    }

    const int numAttributes = computed.first.NUM_ATTRIBUTES;
    const int outputSize = tree.getNodeIdSpaceSize(outputSpace);
    std::vector<float> projected(static_cast<size_t>(outputSize) * numAttributes, std::numeric_limits<float>::quiet_NaN());

    for (NodeId nodeId : tree.getAliveNodeIds()) {
        const NodeId outputNodeId = tree.getHigraNodeId(nodeId);
        if (outputNodeId == InvalidNode) {
            continue;
        }
        for (const auto& [attribute, _] : computed.first.indexMap) {
            projected[computed.first.linearIndex(outputNodeId, attribute)] =
                computed.second[computed.first.linearIndex(nodeId, attribute)];
        }
    }

    return {std::move(computed.first), std::move(projected), outputSpace};
}

/**
 * @brief Delta-aware counterpart of the scalar node-id-space projection.
 */
inline ComputedAttributeDataWithDelta projectComputedDataToNodeIdSpace(const MorphologicalTree& tree, ComputedAttributeDataWithDelta computed, NodeIdSpace outputSpace) {
    if (outputSpace == NodeIdSpace::MORPHOLOGICAL_TREE) {
        computed.nodeIdSpace = NodeIdSpace::MORPHOLOGICAL_TREE;
        return computed;
    }

    const int numAttributes = computed.first.NUM_ATTRIBUTES;
    const int outputSize = tree.getNodeIdSpaceSize(outputSpace);
    std::vector<float> projected(static_cast<size_t>(outputSize) * numAttributes, std::numeric_limits<float>::quiet_NaN());

    for (NodeId nodeId : tree.getAliveNodeIds()) {
        const NodeId outputNodeId = tree.getHigraNodeId(nodeId);
        if (outputNodeId == InvalidNode) {
            continue;
        }
        for (const auto& [attributeKey, _] : computed.first.indexMap) {
            projected[computed.first.linearIndex(outputNodeId, attributeKey.attr, attributeKey.delta)] =
                computed.second[computed.first.linearIndex(nodeId, attributeKey.attr, attributeKey.delta)];
        }
    }

    return {std::move(computed.first), std::move(projected), outputSpace};
}

} // namespace mmcfilters::detail
