#pragma once

#include "../MorphologicalTree.hpp"

namespace mmcfilters::detail::tree_altitude_ops {

inline std::vector<int32_t> computeAreasIncrementally(const MorphologicalTree& tree) {
    std::vector<int32_t> area(static_cast<size_t>(tree.getNumInternalNodeSlots()), 0);
    for (NodeId nodeId : tree.getPostOrderNodes()) {
        area[static_cast<size_t>(nodeId)] += static_cast<int32_t>(tree.getNumProperParts(nodeId));
        const NodeId parentNodeId = tree.getNodeParent(nodeId);
        if (parentNodeId != InvalidNode && parentNodeId != nodeId) {
            area[static_cast<size_t>(parentNodeId)] += area[static_cast<size_t>(nodeId)];
        }
    }
    return area;
}

inline void maxAreaDescendants(
    const MorphologicalTree& tree,
    const std::vector<int32_t>& areaByNode,
    std::vector<NodeId>& descendants,
    NodeId ascendantNodeId,
    NodeId candidateNodeId) {
    if (!tree.isNode(ascendantNodeId) || !tree.isNode(candidateNodeId)) {
        return;
    }
    NodeId& currentNodeId = descendants[static_cast<size_t>(ascendantNodeId)];
    if (ascendantNodeId == candidateNodeId) {
        return;
    }
    if (currentNodeId == InvalidNode ||
        areaByNode[static_cast<size_t>(candidateNodeId)] > areaByNode[static_cast<size_t>(currentNodeId)] ||
        (areaByNode[static_cast<size_t>(candidateNodeId)] == areaByNode[static_cast<size_t>(currentNodeId)] &&
         candidateNodeId < currentNodeId)) {
        currentNodeId = candidateNodeId;
    }
}

} // namespace mmcfilters::detail::tree_altitude_ops
