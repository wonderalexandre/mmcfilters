#pragma once

#include "MorphologicalTree.hpp"

namespace mmcfilters::tree_altitude_ops {

inline const AltitudeBuffer& requireAltitudeBuffer(const AltitudeBuffer* altitude) {
    if (altitude == nullptr) {
        throw std::logic_error("This operation requires an explicit altitude buffer. Use WeightedMorphologicalTree or provide an explicit altitude buffer.");
    }
    return *altitude;
}

inline void validateAltitudeBufferShape(const MorphologicalTree& tree, std::span<const AltitudeType> altitude) {
    if (altitude.size() != static_cast<size_t>(tree.getNumInternalNodeSlots())) {
        throw std::runtime_error("Altitude buffer size must match the dense internal-node domain.");
    }
}

inline void validateAltitudeBufferShape(const MorphologicalTree& tree, const AltitudeBuffer* altitude) {
    validateAltitudeBufferShape(tree, std::span<const AltitudeType>(requireAltitudeBuffer(altitude)));
}

inline AltitudeType getAltitude(std::span<const AltitudeType> altitude, NodeId nodeId) noexcept {
    return altitude[static_cast<size_t>(nodeId)];
}

inline AltitudeType getAltitude(const AltitudeBuffer* altitude, NodeId nodeId) {
    return getAltitude(std::span<const AltitudeType>(requireAltitudeBuffer(altitude)), nodeId);
}

inline AltitudeDiffType getNodeResidue(const MorphologicalTree& tree, std::span<const AltitudeType> altitude, NodeId nodeId) noexcept {
    const NodeId parentNodeId = tree.getNodeParent(nodeId);
    if (parentNodeId == InvalidNode || parentNodeId == nodeId) {
        return getAltitude(altitude, nodeId);
    }
    return getAltitude(altitude, nodeId) - getAltitude(altitude, parentNodeId);
}

inline AltitudeDiffType getNodeResidue(const MorphologicalTree& tree, const AltitudeBuffer* altitude, NodeId nodeId) {
    return getNodeResidue(tree, std::span<const AltitudeType>(requireAltitudeBuffer(altitude)), nodeId);
}

inline NodeId getNodeAscendant(const MorphologicalTree& tree, std::span<const AltitudeType> altitude, NodeId nodeId, int delta, bool useLevel) noexcept {
    NodeId currentNodeId = nodeId;
    if (useLevel) {
        for (int i = 0; i <= delta; i++) {
            if (tree.isMaxtree()) {
                if (getAltitude(altitude, nodeId) >= getAltitude(altitude, currentNodeId) + delta) {
                    return currentNodeId;
                }
            } else {
                if (getAltitude(altitude, nodeId) <= getAltitude(altitude, currentNodeId) - delta) {
                    return currentNodeId;
                }
            }
            if (tree.isRoot(currentNodeId)) {
                return currentNodeId;
            }
            currentNodeId = tree.getNodeParent(currentNodeId);
        }
    } else {
        int step = 0;
        while (step++ < delta) {
            if (tree.isRoot(currentNodeId)) {
                return InvalidNode;
            }
            currentNodeId = tree.getNodeParent(currentNodeId);
        }
    }
    return currentNodeId;
}

inline NodeId getNodeAscendant(const MorphologicalTree& tree, const AltitudeBuffer* altitude, NodeId nodeId, int delta, bool useLevel) {
    return getNodeAscendant(tree, std::span<const AltitudeType>(requireAltitudeBuffer(altitude)), nodeId, delta, useLevel);
}

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
    if (currentNodeId == InvalidNode ||
        areaByNode[static_cast<size_t>(candidateNodeId)] > areaByNode[static_cast<size_t>(currentNodeId)] ||
        (areaByNode[static_cast<size_t>(candidateNodeId)] == areaByNode[static_cast<size_t>(currentNodeId)] &&
         candidateNodeId < currentNodeId)) {
        currentNodeId = candidateNodeId;
    }
}

inline std::pair<std::vector<NodeId>, std::vector<NodeId>> computeAscendantsAndDescendants(
    const MorphologicalTree& tree,
    const AltitudeBuffer* altitude,
    int delta,
    bool useLevel = false) {
    if (!useLevel) {
        return tree.computeAscendantsAndDescendants(delta);
    }

    const AltitudeBuffer& altitudeBuffer = requireAltitudeBuffer(altitude);
    std::vector<NodeId> ascendants(static_cast<size_t>(tree.getNumInternalNodeSlots()), InvalidNode);
    std::vector<NodeId> descendants(static_cast<size_t>(tree.getNumInternalNodeSlots()), InvalidNode);
    const std::vector<int32_t> areaByNode = computeAreasIncrementally(tree);

    for (NodeId nodeId : tree.getAliveNodeIds()) {
        const NodeId ascendantNodeId = getNodeAscendant(tree, altitudeBuffer, nodeId, delta, true);
        if (ascendantNodeId == InvalidNode) {
            continue;
        }
        maxAreaDescendants(tree, areaByNode, descendants, ascendantNodeId, nodeId);
        if (descendants[static_cast<size_t>(ascendantNodeId)] != InvalidNode) {
            ascendants[static_cast<size_t>(nodeId)] = ascendantNodeId;
        }
    }

    return {std::move(ascendants), std::move(descendants)};
}

inline ImageUInt8Ptr reconstructImage(const MorphologicalTree& tree, std::span<const AltitudeType> altitude) {
    validateAltitudeBufferShape(tree, altitude);
    ImageUInt8Ptr image = ImageUInt8::create(tree.getNumRowsOfImage(), tree.getNumColsOfImage());
    auto imgBuffer = image->rawData();
    for (int pixelId = 0; pixelId < tree.getNumTotalProperParts(); ++pixelId) {
        const NodeId nodeId = tree.getSmallestComponent(pixelId);
        imgBuffer[static_cast<size_t>(pixelId)] = static_cast<uint8_t>(getAltitude(altitude, nodeId));
    }
    return image;
}

inline ImageUInt8Ptr reconstructImage(const MorphologicalTree& tree, const AltitudeBuffer* altitude) {
    return reconstructImage(tree, std::span<const AltitudeType>(requireAltitudeBuffer(altitude)));
}

inline std::pair<std::vector<NodeId>, std::vector<AltitudeType>> exportHigraHierarchy(
    const MorphologicalTree& tree,
    std::span<const AltitudeType> altitude) {
    validateAltitudeBufferShape(tree, altitude);

    if (tree.getRoot() == InvalidNode || !tree.isAlive(tree.getRoot())) {
        throw std::runtime_error("Cannot export a tree without a valid rooted component.");
    }

    std::vector<NodeId> exportedNodes;
    exportedNodes.reserve(static_cast<size_t>(tree.getNumNodes()));
    for (NodeId nodeId : tree.getNodeSubtree(tree.getRoot())) {
        exportedNodes.push_back(nodeId);
    }

    if (static_cast<int>(exportedNodes.size()) != tree.getNumNodes()) {
        throw std::runtime_error("Cannot export a forest or a tree with detached alive nodes to Higra's static representation.");
    }

    const NodeId numLeaves = tree.getNumTotalProperParts();
    const NodeId numAliveNodes = static_cast<NodeId>(exportedNodes.size());
    const NodeId numVertices = numLeaves + numAliveNodes;

    std::vector<NodeId> parent(static_cast<size_t>(numVertices), InvalidNode);
    std::vector<AltitudeType> exportedAltitude(static_cast<size_t>(numVertices), AltitudeType{});
    std::vector<NodeId> oldToNew(static_cast<size_t>(tree.getNumInternalNodeSlots()), InvalidNode);
    auto sortedNodes = exportedNodes;

    bool sortAscendingAltitude = true;
    for (NodeId nodeId : sortedNodes) {
        if (tree.isRoot(nodeId)) {
            continue;
        }

        const NodeId parentNodeId = tree.getNodeParent(nodeId);
        if (parentNodeId == InvalidNode || !tree.isAlive(parentNodeId)) {
            throw std::runtime_error("Cannot export a node whose parent is not part of the rooted alive component.");
        }
        if (getAltitude(altitude, nodeId) > getAltitude(altitude, parentNodeId)) {
            sortAscendingAltitude = false;
        }
    }

    std::stable_sort(
        sortedNodes.begin(),
        sortedNodes.end(),
        [&](NodeId lhs, NodeId rhs) {
            const AltitudeType altL = getAltitude(altitude, lhs);
            const AltitudeType altR = getAltitude(altitude, rhs);
            if (altL != altR) {
                return sortAscendingAltitude ? altL < altR : altL > altR;
            }
            return lhs < rhs;
        });

    for (NodeId i = 0; i < numAliveNodes; ++i) {
        const NodeId oldNodeId = sortedNodes[static_cast<size_t>(i)];
        const NodeId newNodeId = numLeaves + i;
        oldToNew[static_cast<size_t>(oldNodeId)] = newNodeId;
        exportedAltitude[static_cast<size_t>(newNodeId)] = getAltitude(altitude, oldNodeId);
    }

    for (NodeId properPart = 0; properPart < numLeaves; ++properPart) {
        const NodeId ownerNodeId = tree.getSmallestComponent(properPart);
        if (ownerNodeId == InvalidNode || !tree.isAlive(ownerNodeId)) {
            throw std::runtime_error("Each proper part must belong to one alive node when exporting to Higra.");
        }
        parent[static_cast<size_t>(properPart)] = oldToNew[static_cast<size_t>(ownerNodeId)];
        exportedAltitude[static_cast<size_t>(properPart)] = getAltitude(altitude, ownerNodeId);
    }

    for (NodeId oldNodeId : sortedNodes) {
        const NodeId newNodeId = oldToNew[static_cast<size_t>(oldNodeId)];
        const NodeId oldParentNodeId = tree.getNodeParent(oldNodeId);
        parent[static_cast<size_t>(newNodeId)] =
            oldParentNodeId == oldNodeId ? newNodeId : oldToNew[static_cast<size_t>(oldParentNodeId)];
    }

    return {std::move(parent), std::move(exportedAltitude)};
}

inline std::pair<std::vector<NodeId>, std::vector<AltitudeType>> exportHigraHierarchy(
    const MorphologicalTree& tree,
    const AltitudeBuffer* altitude) {
    return exportHigraHierarchy(tree, std::span<const AltitudeType>(requireAltitudeBuffer(altitude)));
}

inline void validateMonotoneAltitude(const MorphologicalTree& tree, std::span<const AltitudeType> altitude) {
    validateAltitudeBufferShape(tree, altitude);
    if (tree.getTreeType() == MorphologicalTree::TREE_OF_SHAPES) {
        return;
    }

    const bool increasingTowardLeaves = tree.isMaxtree();
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        if (tree.isRoot(nodeId)) {
            continue;
        }

        const NodeId parentNodeId = tree.getNodeParent(nodeId);
        if (parentNodeId == InvalidNode || !tree.isAlive(parentNodeId)) {
            throw std::runtime_error("Monotonic validation requires every alive non-root node to have an alive parent.");
        }

        if (increasingTowardLeaves) {
            if (getAltitude(altitude, parentNodeId) > getAltitude(altitude, nodeId)) {
                throw std::runtime_error("Max-tree altitude buffer must be non-decreasing from parent to child.");
            }
        } else if (getAltitude(altitude, parentNodeId) < getAltitude(altitude, nodeId)) {
            throw std::runtime_error("Min-tree altitude buffer must be non-increasing from parent to child.");
        }
    }
}

inline void validateMonotoneAltitude(const MorphologicalTree& tree, const AltitudeBuffer* altitude) {
    validateMonotoneAltitude(tree, std::span<const AltitudeType>(requireAltitudeBuffer(altitude)));
}

} // namespace mmcfilters::tree_altitude_ops
