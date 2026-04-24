#pragma once

#include "../utils/AdjacencyRelation.hpp"
#include "../utils/Common.hpp"
#include "BuilderMorphologicalTreeByUnionFind.hpp"

#include <algorithm>
#include <cstddef>
#include <deque>
#include <optional>
#include <queue>
#include <span>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>

namespace mmcfilters {

/**
 * @brief Selects the node-id domain used by attribute buffers exposed to callers.
 */
enum class NodeIdSpace {
    MORPHOLOGICAL_TREE,
    HIGRA
};

/**
 * @brief Selects the interpolation mode used to build a tree of shapes.
 */
enum class ToSInterpolation {
    SelfDual,
    Min4cMax8c
};

class MorphologicalTree;
class MorphologicalTreePybind;
class TreeEditor;
class WeightedMorphologicalTree;

/**
 * @brief Mutable topology and proper-part ownership for component-tree-like hierarchies.
 *
 * The class owns only structural state: dense internal node ids, parent/child
 * links, and the smallest component that owns each proper part.
 */
class MorphologicalTree {
    friend class MorphologicalTreePybind;
    friend class TreeEditor;
    friend class WeightedMorphologicalTree;

public:
    static constexpr int MAX_TREE = 0;
    static constexpr int MIN_TREE = 1;
    static constexpr int TREE_OF_SHAPES = 2;

private:
    NodeId rootNodeId_ = InvalidNode;
    int treeType_ = MAX_TREE; // 0=max-tree, 1=min-tree, 2=tree of shapes
    int numRows_ = 0;
    int numCols_ = 0;
    int numNodes_ = 0;
    std::optional<AdjacencyRelation> adj_;

    std::vector<NodeId> properPartOwner_;
    std::vector<NodeId> nodeParent_;
    std::vector<std::vector<NodeId>> children_;
    std::vector<uint8_t> alive_;
    std::vector<NodeId> freeNodeIds_;
    std::vector<std::vector<NodeId>> properPartsByNode_;

    std::vector<NodeId> higraNodeIdBySlot_;
    std::vector<NodeId> slotByHigraNodeId_;

    mutable bool traversalTimesValid_ = false;
    mutable std::vector<int> preOrderTime_;
    mutable std::vector<int> postOrderTime_;
    mutable size_t nodeStructureVersion_ = 0;
    mutable size_t topologyVersion_ = 0;
    mutable size_t properPartVersion_ = 0;

    void invalidateTraversalTimes() const noexcept {
        traversalTimesValid_ = false;
    }

    void invalidateHigraNodeIdMapping() {
        higraNodeIdBySlot_.clear();
        slotByHigraNodeId_.clear();
    }

    void resetIteratorVersions() const noexcept {
        ++nodeStructureVersion_;
        ++topologyVersion_;
        ++properPartVersion_;
        invalidateTraversalTimes();
    }

    void markTopologyChanged() {
        resetIteratorVersions();
        invalidateHigraNodeIdMapping();
    }

    void markProperPartChanged() {
        ++properPartVersion_;
        invalidateHigraNodeIdMapping();
    }

    void ensureNodeId(NodeId nodeId, const char* context) const {
        if (!isNode(nodeId)) {
            throw std::invalid_argument(context);
        }
    }

    void ensureAliveNodeId(NodeId nodeId, const char* context) const {
        if (!isAlive(nodeId)) {
            throw std::invalid_argument(context);
        }
    }

    void resetEmptyStorage(size_t numProperParts) {
        rootNodeId_ = InvalidNode;
        numNodes_ = 0;
        properPartOwner_.assign(numProperParts, InvalidNode);
        nodeParent_.clear();
        children_.clear();
        alive_.clear();
        freeNodeIds_.clear();
        properPartsByNode_.clear();
        invalidateHigraNodeIdMapping();
        resetIteratorVersions();
    }

    NodeId allocateNodeSlot() {
        NodeId nodeId = InvalidNode;
        if (!freeNodeIds_.empty()) {
            nodeId = freeNodeIds_.back();
            freeNodeIds_.pop_back();
            const auto index = static_cast<size_t>(nodeId);
            alive_[index] = true;
            nodeParent_[index] = nodeId;
            children_[index].clear();
            properPartsByNode_[index].clear();
        } else {
            nodeId = static_cast<NodeId>(nodeParent_.size());
            nodeParent_.push_back(nodeId);
            children_.emplace_back();
            alive_.push_back(true);
            properPartsByNode_.emplace_back();
        }
        ++numNodes_;
        return nodeId;
    }

    NodeId createDetachedNode() {
        const NodeId nodeId = allocateNodeSlot();
        if (rootNodeId_ == InvalidNode) {
            rootNodeId_ = nodeId;
        }
        markTopologyChanged();
        return nodeId;
    }

    NodeId createNode(NodeId parentNodeId) {
        const NodeId nodeId = allocateNodeSlot();
        if (parentNodeId == InvalidNode) {
            rootNodeId_ = nodeId;
            nodeParent_[static_cast<size_t>(nodeId)] = nodeId;
        } else {
            children_[static_cast<size_t>(parentNodeId)].push_back(nodeId);
            nodeParent_[static_cast<size_t>(nodeId)] = parentNodeId;
        }
        resetIteratorVersions();
        return nodeId;
    }

    void releaseNode(NodeId nodeId) {
        ensureAliveNodeId(nodeId, "releaseNode requires a live node.");
        const auto index = static_cast<size_t>(nodeId);
        if (!children_[index].empty() || !properPartsByNode_[index].empty()) {
            throw std::logic_error("releaseNode requires an empty detached node.");
        }
        if (nodeId == rootNodeId_) {
            rootNodeId_ = InvalidNode;
        }
        alive_[index] = false;
        nodeParent_[index] = InvalidNode;
        freeNodeIds_.push_back(nodeId);
        --numNodes_;
        markTopologyChanged();
    }

    void removeChildLink(NodeId parentNodeId, NodeId childNodeId) {
        if (!isNode(parentNodeId)) {
            return;
        }
        auto& siblings = children_[static_cast<size_t>(parentNodeId)];
        siblings.erase(std::remove(siblings.begin(), siblings.end(), childNodeId), siblings.end());
    }

    void linkChild(NodeId parentNodeId, NodeId childNodeId) {
        ensureAliveNodeId(parentNodeId, "linkChild requires a live parent node.");
        ensureAliveNodeId(childNodeId, "linkChild requires a live child node.");
        if (parentNodeId == childNodeId) {
            throw std::invalid_argument("A node cannot be linked as its own child.");
        }
        const NodeId oldParentId = nodeParent_[static_cast<size_t>(childNodeId)];
        if (oldParentId != InvalidNode && oldParentId != childNodeId) {
            removeChildLink(oldParentId, childNodeId);
        }
        nodeParent_[static_cast<size_t>(childNodeId)] = parentNodeId;
        children_[static_cast<size_t>(parentNodeId)].push_back(childNodeId);
        markTopologyChanged();
    }

    void rebuildProperPartsFromOwnership() {
        properPartsByNode_.assign(nodeParent_.size(), {});
        for (NodeId properPart = 0; properPart < static_cast<NodeId>(properPartOwner_.size()); ++properPart) {
            const NodeId ownerNodeId = properPartOwner_[static_cast<size_t>(properPart)];
            if (ownerNodeId != InvalidNode && isAlive(ownerNodeId)) {
                properPartsByNode_[static_cast<size_t>(ownerNodeId)].push_back(properPart);
            }
        }
        ++properPartVersion_;
    }

    int expectedProperPartCountFromDomain() const noexcept {
        if (numRows_ <= 0 || numCols_ <= 0) {
            return -1;
        }
        return numRows_ * numCols_;
    }

    void validateProperPartDomain(NodeId numProperParts) const {
        if (numProperParts < 0) {
            throw std::invalid_argument("The proper-part domain size must be non-negative.");
        }
        const int imageDomainSize = expectedProperPartCountFromDomain();
        if (imageDomainSize >= 0 && numProperParts != imageDomainSize) {
            throw std::invalid_argument("The proper-part domain must match the configured image domain.");
        }
        if (adj_.has_value()) {
            const int adjacencyDomainSize = adj_->getNumRows() * adj_->getNumCols();
            if (numProperParts != adjacencyDomainSize) {
                throw std::invalid_argument("The proper-part domain must match the adjacency relation domain.");
            }
        }
    }

    void validateParentArrayInput(std::span<const NodeId> parent, NodeId numProperParts) const {
        validateProperPartDomain(numProperParts);
        if (parent.size() <= static_cast<size_t>(numProperParts)) {
            throw std::invalid_argument("The compact parent array must contain at least one internal node.");
        }

        const NodeId numNodeSlots = static_cast<NodeId>(parent.size()) - numProperParts;
        std::vector<uint8_t> nodeAlive(static_cast<size_t>(numNodeSlots), false);
        for (NodeId nodeId = 0; nodeId < numNodeSlots; ++nodeId) {
            if (parent[static_cast<size_t>(numProperParts + nodeId)] != InvalidNode) {
                nodeAlive[static_cast<size_t>(nodeId)] = true;
            }
        }

        for (NodeId properPart = 0; properPart < numProperParts; ++properPart) {
            const NodeId ownerNodeId = parent[static_cast<size_t>(properPart)];
            if (ownerNodeId < 0 || ownerNodeId >= numNodeSlots || !nodeAlive[static_cast<size_t>(ownerNodeId)]) {
                throw std::invalid_argument("Each proper part must point to a live internal node.");
            }
        }

        NodeId rootNodeId = InvalidNode;
        for (NodeId nodeId = 0; nodeId < numNodeSlots; ++nodeId) {
            if (!nodeAlive[static_cast<size_t>(nodeId)]) {
                continue;
            }
            const NodeId parentNodeId = parent[static_cast<size_t>(numProperParts + nodeId)];
            if (parentNodeId < 0 || parentNodeId >= numNodeSlots || !nodeAlive[static_cast<size_t>(parentNodeId)]) {
                throw std::invalid_argument("Each live internal node must point to a live internal parent.");
            }
            if (parentNodeId == nodeId) {
                if (rootNodeId != InvalidNode) {
                    throw std::invalid_argument("The compact parent array must contain exactly one root.");
                }
                rootNodeId = nodeId;
            }
        }
        if (rootNodeId == InvalidNode) {
            throw std::invalid_argument("The compact parent array must contain exactly one root.");
        }

        for (NodeId nodeId = 0; nodeId < numNodeSlots; ++nodeId) {
            if (!nodeAlive[static_cast<size_t>(nodeId)]) {
                continue;
            }
            std::vector<uint8_t> seen(static_cast<size_t>(numNodeSlots), false);
            NodeId currentNodeId = nodeId;
            while (true) {
                if (currentNodeId < 0 || currentNodeId >= numNodeSlots || !nodeAlive[static_cast<size_t>(currentNodeId)]) {
                    throw std::invalid_argument("Parent chains must remain inside the live node domain.");
                }
                if (seen[static_cast<size_t>(currentNodeId)]) {
                    throw std::invalid_argument("Parent chains must be acyclic.");
                }
                seen[static_cast<size_t>(currentNodeId)] = true;
                const NodeId parentNodeId = parent[static_cast<size_t>(numProperParts + currentNodeId)];
                if (parentNodeId == currentNodeId) {
                    if (currentNodeId != rootNodeId) {
                        throw std::invalid_argument("Every live node must reach the unique root.");
                    }
                    break;
                }
                currentNodeId = parentNodeId;
            }
        }
    }

    void validateHigraTopologyInput(std::span<const NodeId> parent, NodeId numProperParts) const {
        validateProperPartDomain(numProperParts);
        if (parent.size() <= static_cast<size_t>(numProperParts)) {
            throw std::invalid_argument("The Higra parent array must contain leaves and internal nodes.");
        }

        const NodeId numHigraNodes = static_cast<NodeId>(parent.size());
        for (NodeId properPart = 0; properPart < numProperParts; ++properPart) {
            const NodeId parentNodeId = parent[static_cast<size_t>(properPart)];
            if (parentNodeId < numProperParts || parentNodeId >= numHigraNodes) {
                throw std::invalid_argument("Higra leaves must point to the internal-node domain.");
            }
        }

        NodeId rootHigraId = InvalidNode;
        for (NodeId higraNodeId = numProperParts; higraNodeId < numHigraNodes; ++higraNodeId) {
            const NodeId parentNodeId = parent[static_cast<size_t>(higraNodeId)];
            if (parentNodeId < numProperParts || parentNodeId >= numHigraNodes) {
                throw std::invalid_argument("Higra internal nodes must point to internal-node parents.");
            }
            if (parentNodeId == higraNodeId) {
                if (rootHigraId != InvalidNode) {
                    throw std::invalid_argument("The Higra parent array must contain exactly one root.");
                }
                rootHigraId = higraNodeId;
            }
        }
        if (rootHigraId == InvalidNode) {
            throw std::invalid_argument("The Higra parent array must contain exactly one root.");
        }

        for (NodeId higraNodeId = numProperParts; higraNodeId < numHigraNodes; ++higraNodeId) {
            std::vector<uint8_t> seen(static_cast<size_t>(numHigraNodes), false);
            NodeId currentNodeId = higraNodeId;
            while (true) {
                if (currentNodeId < numProperParts || currentNodeId >= numHigraNodes) {
                    throw std::invalid_argument("Higra parent chains must remain inside the internal-node domain.");
                }
                if (seen[static_cast<size_t>(currentNodeId)]) {
                    throw std::invalid_argument("Higra parent chains must be acyclic.");
                }
                seen[static_cast<size_t>(currentNodeId)] = true;
                const NodeId parentNodeId = parent[static_cast<size_t>(currentNodeId)];
                if (parentNodeId == currentNodeId) {
                    if (currentNodeId != rootHigraId) {
                        throw std::invalid_argument("Every Higra internal node must reach the unique root.");
                    }
                    break;
                }
                currentNodeId = parentNodeId;
            }
        }
    }

    void resetFromParentArray(std::span<const NodeId> parent, NodeId numProperParts) {
        validateParentArrayInput(parent, numProperParts);
        const NodeId numNodeSlots = static_cast<NodeId>(parent.size()) - numProperParts;

        resetEmptyStorage(static_cast<size_t>(numProperParts));
        nodeParent_.assign(static_cast<size_t>(numNodeSlots), InvalidNode);
        children_.assign(static_cast<size_t>(numNodeSlots), {});
        alive_.assign(static_cast<size_t>(numNodeSlots), false);
        properPartsByNode_.assign(static_cast<size_t>(numNodeSlots), {});
        freeNodeIds_.clear();

        for (NodeId nodeId = 0; nodeId < numNodeSlots; ++nodeId) {
            const NodeId parentNodeId = parent[static_cast<size_t>(numProperParts + nodeId)];
            if (parentNodeId == InvalidNode) {
                freeNodeIds_.push_back(nodeId);
                continue;
            }
            alive_[static_cast<size_t>(nodeId)] = true;
            nodeParent_[static_cast<size_t>(nodeId)] = parentNodeId;
            ++numNodes_;
            if (parentNodeId == nodeId) {
                rootNodeId_ = nodeId;
            } else {
                children_[static_cast<size_t>(parentNodeId)].push_back(nodeId);
            }
        }

        for (NodeId properPart = 0; properPart < numProperParts; ++properPart) {
            properPartOwner_[static_cast<size_t>(properPart)] = parent[static_cast<size_t>(properPart)];
        }
        rebuildProperPartsFromOwnership();
        invalidateHigraNodeIdMapping();
        resetIteratorVersions();
    }

    void resetFromHigraTopology(std::span<const NodeId> parent, NodeId numProperParts) {
        validateHigraTopologyInput(parent, numProperParts);
        const NodeId numHigraNodes = static_cast<NodeId>(parent.size());
        const NodeId numNodeSlots = numHigraNodes - numProperParts;

        resetEmptyStorage(static_cast<size_t>(numProperParts));
        nodeParent_.assign(static_cast<size_t>(numNodeSlots), InvalidNode);
        children_.assign(static_cast<size_t>(numNodeSlots), {});
        alive_.assign(static_cast<size_t>(numNodeSlots), true);
        properPartsByNode_.assign(static_cast<size_t>(numNodeSlots), {});
        freeNodeIds_.clear();
        numNodes_ = numNodeSlots;

        higraNodeIdBySlot_.assign(static_cast<size_t>(numNodeSlots), InvalidNode);
        slotByHigraNodeId_.assign(static_cast<size_t>(numHigraNodes), InvalidNode);
        for (NodeId higraNodeId = numProperParts; higraNodeId < numHigraNodes; ++higraNodeId) {
            const NodeId slotId = higraNodeId - numProperParts;
            higraNodeIdBySlot_[static_cast<size_t>(slotId)] = higraNodeId;
            slotByHigraNodeId_[static_cast<size_t>(higraNodeId)] = slotId;
        }

        for (NodeId properPart = 0; properPart < numProperParts; ++properPart) {
            properPartOwner_[static_cast<size_t>(properPart)] =
                slotByHigraNodeId_[static_cast<size_t>(parent[static_cast<size_t>(properPart)])];
        }

        for (NodeId higraNodeId = numProperParts; higraNodeId < numHigraNodes; ++higraNodeId) {
            const NodeId slotId = higraNodeId - numProperParts;
            const NodeId parentHigraId = parent[static_cast<size_t>(higraNodeId)];
            if (parentHigraId == higraNodeId) {
                nodeParent_[static_cast<size_t>(slotId)] = slotId;
                rootNodeId_ = slotId;
            } else {
                const NodeId parentSlotId = slotByHigraNodeId_[static_cast<size_t>(parentHigraId)];
                nodeParent_[static_cast<size_t>(slotId)] = parentSlotId;
                children_[static_cast<size_t>(parentSlotId)].push_back(slotId);
            }
        }

        rebuildProperPartsFromOwnership();
        resetIteratorVersions();
    }

    void build(ImageUInt8Ptr img, const IMorphologicalTreeBuilder& builder) {
        if (!img) {
            throw std::invalid_argument("MorphologicalTree construction requires a non-null image.");
        }
        numRows_ = img->getNumRows();
        numCols_ = img->getNumCols();
        resetEmptyStorage(static_cast<size_t>(img->getSize()));

        auto [parentPixels, orderedPixels, expectedNumNodes] = builder.createTreeByUnionFind(img);
        if (parentPixels.size() != static_cast<size_t>(img->getSize()) ||
            orderedPixels.size() != static_cast<size_t>(img->getSize())) {
            throw std::runtime_error("Union-find builder returned an invalid image-domain hierarchy.");
        }

        auto* pixel = img->rawData();
        for (int properPart : orderedPixels) {
            if (properPart < 0 || properPart >= img->getSize()) {
                throw std::runtime_error("Union-find builder returned an invalid pixel id.");
            }

            const int parentProperPart = parentPixels[static_cast<size_t>(properPart)];
            if (parentProperPart < 0 || parentProperPart >= img->getSize()) {
                throw std::runtime_error("Union-find builder returned an invalid parent pixel id.");
            }

            NodeId ownerNodeId = InvalidNode;
            if (parentProperPart == properPart) {
                ownerNodeId = createNode(InvalidNode);
            } else if (pixel[parentProperPart] != pixel[properPart]) {
                const NodeId parentNodeId = properPartOwner_[static_cast<size_t>(parentProperPart)];
                if (parentNodeId == InvalidNode) {
                    throw std::runtime_error("Union-find builder produced a child before its parent component.");
                }
                ownerNodeId = createNode(parentNodeId);
            } else {
                ownerNodeId = properPartOwner_[static_cast<size_t>(parentProperPart)];
                if (ownerNodeId == InvalidNode) {
                    throw std::runtime_error("Union-find builder produced a flat zone before its representative component.");
                }
            }
            properPartOwner_[static_cast<size_t>(properPart)] = ownerNodeId;
        }

        rebuildProperPartsFromOwnership();
        if (expectedNumNodes != numNodes_) {
            throw std::runtime_error("Union-find builder node count does not match the imported topology.");
        }
        validateConnectedRootedTree();
    }

    void detachNode(NodeId nodeId) {
        ensureAliveNodeId(nodeId, "detachNode requires a live node.");
        const NodeId parentNodeId = nodeParent_[static_cast<size_t>(nodeId)];
        if (parentNodeId != InvalidNode && parentNodeId != nodeId) {
            removeChildLink(parentNodeId, nodeId);
        }
        nodeParent_[static_cast<size_t>(nodeId)] = nodeId;
        markTopologyChanged();
    }

    void moveNode(NodeId nodeId, NodeId newParentId) {
        ensureAliveNodeId(nodeId, "moveNode requires a live node.");
        ensureAliveNodeId(newParentId, "moveNode requires a live parent node.");
        if (nodeId == rootNodeId_) {
            throw std::invalid_argument("moveNode cannot move the connected root.");
        }
        if (nodeId == newParentId || isAncestor(nodeId, newParentId)) {
            throw std::invalid_argument("moveNode would create a cycle.");
        }
        linkChild(newParentId, nodeId);
    }

    void attachNode(NodeId parentNodeId, NodeId detachedNodeId) {
        ensureAliveNodeId(parentNodeId, "attachNode requires a live parent node.");
        ensureAliveNodeId(detachedNodeId, "attachNode requires a live detached node.");
        if (nodeParent_[static_cast<size_t>(detachedNodeId)] != detachedNodeId) {
            throw std::invalid_argument("attachNode expects a detached self-parented node.");
        }
        if (detachedNodeId == parentNodeId || isAncestor(detachedNodeId, parentNodeId)) {
            throw std::invalid_argument("attachNode would create a cycle.");
        }
        linkChild(parentNodeId, detachedNodeId);
    }

    void moveChildren(NodeId targetNodeId, NodeId sourceNodeId) {
        ensureAliveNodeId(targetNodeId, "moveChildren requires a live target node.");
        ensureAliveNodeId(sourceNodeId, "moveChildren requires a live source node.");
        if (targetNodeId == sourceNodeId || isAncestor(sourceNodeId, targetNodeId)) {
            throw std::invalid_argument("moveChildren would create a cycle.");
        }

        auto movingChildren = std::move(children_[static_cast<size_t>(sourceNodeId)]);
        children_[static_cast<size_t>(sourceNodeId)].clear();
        auto& targetChildren = children_[static_cast<size_t>(targetNodeId)];
        for (NodeId childNodeId : movingChildren) {
            nodeParent_[static_cast<size_t>(childNodeId)] = targetNodeId;
            targetChildren.push_back(childNodeId);
        }
        markTopologyChanged();
    }

    void moveProperPart(NodeId targetNodeId, NodeId sourceNodeId, NodeId properPartId) {
        ensureAliveNodeId(targetNodeId, "moveProperPart requires a live target node.");
        ensureAliveNodeId(sourceNodeId, "moveProperPart requires a live source node.");
        if (!isProperPart(properPartId)) {
            throw std::invalid_argument("moveProperPart requires a valid proper-part id.");
        }
        if (properPartOwner_[static_cast<size_t>(properPartId)] != sourceNodeId) {
            throw std::invalid_argument("moveProperPart requires a direct proper part of the source node.");
        }
        properPartOwner_[static_cast<size_t>(properPartId)] = targetNodeId;
        rebuildProperPartsFromOwnership();
        markProperPartChanged();
    }

    void moveProperParts(NodeId targetNodeId, NodeId sourceNodeId) {
        ensureAliveNodeId(targetNodeId, "moveProperParts requires a live target node.");
        ensureAliveNodeId(sourceNodeId, "moveProperParts requires a live source node.");
        for (NodeId& ownerNodeId : properPartOwner_) {
            if (ownerNodeId == sourceNodeId) {
                ownerNodeId = targetNodeId;
            }
        }
        rebuildProperPartsFromOwnership();
        markProperPartChanged();
    }

    void setRoot(NodeId nodeId) {
        ensureAliveNodeId(nodeId, "setRoot requires a live node.");
        detachNode(nodeId);
        rootNodeId_ = nodeId;
        nodeParent_[static_cast<size_t>(nodeId)] = nodeId;
        markTopologyChanged();
    }

    void rebuildTraversalTimes() const {
        preOrderTime_.assign(nodeParent_.size(), -1);
        postOrderTime_.assign(nodeParent_.size(), -1);
        if (!isAlive(rootNodeId_)) {
            traversalTimesValid_ = true;
            return;
        }

        int timestamp = 0;
        std::vector<std::pair<NodeId, size_t>> stack;
        stack.emplace_back(rootNodeId_, 0);
        preOrderTime_[static_cast<size_t>(rootNodeId_)] = timestamp++;

        while (!stack.empty()) {
            auto& [nodeId, nextChildIndex] = stack.back();
            const auto& children = children_[static_cast<size_t>(nodeId)];
            if (nextChildIndex < children.size()) {
                const NodeId childNodeId = children[nextChildIndex++];
                if (!isAlive(childNodeId)) {
                    continue;
                }
                preOrderTime_[static_cast<size_t>(childNodeId)] = timestamp++;
                stack.emplace_back(childNodeId, 0);
            } else {
                postOrderTime_[static_cast<size_t>(nodeId)] = timestamp++;
                stack.pop_back();
            }
        }
        traversalTimesValid_ = true;
    }

    std::vector<int32_t> computeAreasIncrementally() const {
        std::vector<int32_t> area(nodeParent_.size(), 0);
        for (NodeId nodeId : getPostOrderNodes()) {
            area[static_cast<size_t>(nodeId)] += static_cast<int32_t>(getNumProperParts(nodeId));
            const NodeId parentNodeId = getNodeParent(nodeId);
            if (parentNodeId != InvalidNode && parentNodeId != nodeId) {
                area[static_cast<size_t>(parentNodeId)] += area[static_cast<size_t>(nodeId)];
            }
        }
        return area;
    }

    void maxAreaDescendants(
        const std::vector<int32_t>& areaByNode,
        std::vector<NodeId>& descendants,
        NodeId ascendantNodeId,
        NodeId candidateNodeId) const {
        if (!isNode(ascendantNodeId) || !isNode(candidateNodeId)) {
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

public:
    MorphologicalTree() = default;

    MorphologicalTree(ImageUInt8Ptr img, ToSInterpolation interpolation = ToSInterpolation::SelfDual) {
        treeType_ = TREE_OF_SHAPES;
        adj_ = std::nullopt;
        BuilderTreeOfShape builder(interpolation == ToSInterpolation::Min4cMax8c);
        build(img, builder);
    }

    explicit MorphologicalTree(ImageUInt8Ptr img, bool isMaxtree, double radius = 1.5) {
        treeType_ = isMaxtree ? MAX_TREE : MIN_TREE;
        numRows_ = img ? img->getNumRows() : 0;
        numCols_ = img ? img->getNumCols() : 0;
        adj_.emplace(numRows_, numCols_, radius);
        BuilderComponentTree builder(&*adj_, isMaxtree);
        build(img, builder);
    }

    MorphologicalTree(
        std::span<const NodeId> parent,
        NodeId numProperParts,
        int rows,
        int cols,
        bool isMaxtree,
        std::optional<AdjacencyRelation> adjacency = std::nullopt)
        : treeType_(isMaxtree ? MAX_TREE : MIN_TREE),
          numRows_(rows),
          numCols_(cols),
          adj_(std::move(adjacency)) {
        resetFromParentArray(parent, numProperParts);
    }

    MorphologicalTree(
        std::span<const NodeId> parent,
        int rows,
        int cols,
        bool isMaxtree,
        std::optional<AdjacencyRelation> adjacency = std::nullopt)
        : MorphologicalTree(parent, static_cast<NodeId>(rows * cols), rows, cols, isMaxtree, std::move(adjacency)) {}

    template <class ParentRange>
        requires requires(const ParentRange& range) {
            range.data();
            range.size();
        }
    MorphologicalTree(
        const ParentRange& parent,
        NodeId numProperParts,
        int rows,
        int cols,
        bool isMaxtree,
        std::optional<AdjacencyRelation> adjacency = std::nullopt)
        : MorphologicalTree(
              std::span<const NodeId>(parent.data(), parent.size()),
              numProperParts,
              rows,
              cols,
              isMaxtree,
              std::move(adjacency)) {}

    template <class ParentRange>
        requires requires(const ParentRange& range) {
            range.data();
            range.size();
        }
    MorphologicalTree(
        const ParentRange& parent,
        int rows,
        int cols,
        bool isMaxtree,
        std::optional<AdjacencyRelation> adjacency = std::nullopt)
        : MorphologicalTree(
              std::span<const NodeId>(parent.data(), parent.size()),
              rows,
              cols,
              isMaxtree,
              std::move(adjacency)) {}

    MorphologicalTree(const MorphologicalTree&) = delete;
    MorphologicalTree& operator=(const MorphologicalTree&) = delete;
    MorphologicalTree(MorphologicalTree&&) noexcept = default;
    MorphologicalTree& operator=(MorphologicalTree&&) noexcept = default;

    static MorphologicalTree create(int rows, int cols, bool isMaxtree, std::optional<AdjacencyRelation> adjacency = std::nullopt) {
        MorphologicalTree tree;
        tree.treeType_ = isMaxtree ? MAX_TREE : MIN_TREE;
        tree.numRows_ = rows;
        tree.numCols_ = cols;
        tree.adj_ = adjacency ? std::move(adjacency) : std::optional<AdjacencyRelation>(std::in_place, rows, cols, 1.5);
        tree.resetEmptyStorage(static_cast<size_t>(rows * cols));
        return tree;
    }

    void reset(std::span<const NodeId> parent, NodeId numProperParts) {
        resetFromParentArray(parent, numProperParts);
    }

    void reset(std::span<const NodeId> parent) {
        if (properPartOwner_.empty()) {
            throw std::logic_error("reset(parent) requires a known proper-part domain.");
        }
        resetFromParentArray(parent, static_cast<NodeId>(properPartOwner_.size()));
    }

    template <class ParentRange>
        requires requires(const ParentRange& range) {
            range.data();
            range.size();
        }
    void reset(const ParentRange& parent, NodeId numProperParts) {
        reset(std::span<const NodeId>(parent.data(), parent.size()), numProperParts);
    }

    template <class ParentRange>
        requires requires(const ParentRange& range) {
            range.data();
            range.size();
        }
    void reset(const ParentRange& parent) {
        reset(std::span<const NodeId>(parent.data(), parent.size()));
    }

    int getNumInternalNodeSlots() const noexcept { return static_cast<int>(nodeParent_.size()); }
    int getNumTotalProperParts() const noexcept { return static_cast<int>(properPartOwner_.size()); }
    int getNumNodes() const noexcept { return numNodes_; }
    NodeId getRoot() const noexcept { return rootNodeId_; }
    int getTreeType() const noexcept { return treeType_; }
    bool isMaxtree() const noexcept { return treeType_ == MAX_TREE; }
    int getNumRowsOfImage() const noexcept { return numRows_; }
    int getNumColsOfImage() const noexcept { return numCols_; }
    bool hasAdjacencyRelation() const noexcept { return adj_.has_value(); }

    AdjacencyRelation* getAdjacencyRelation() noexcept {
        return adj_ ? &*adj_ : nullptr;
    }

    const AdjacencyRelation* getAdjacencyRelation() const noexcept {
        return adj_ ? &*adj_ : nullptr;
    }

    int getNumFreeNodeSlots() const noexcept { return static_cast<int>(freeNodeIds_.size()); }

    bool hasHigraNodeIdMapping() const noexcept {
        return !higraNodeIdBySlot_.empty() && !slotByHigraNodeId_.empty();
    }

    int getNumHigraNodes() const noexcept {
        return hasHigraNodeIdMapping() ? static_cast<int>(slotByHigraNodeId_.size()) : 0;
    }

    int getNodeIdSpaceSize(NodeIdSpace nodeIdSpace) const noexcept {
        switch (nodeIdSpace) {
            case NodeIdSpace::MORPHOLOGICAL_TREE:
                return getNumInternalNodeSlots();
            case NodeIdSpace::HIGRA:
                return getNumHigraNodes();
        }
        return 0;
    }

    NodeId getHigraNodeId(NodeId nodeId) const noexcept {
        if (!hasHigraNodeIdMapping() || !isNode(nodeId)) {
            return InvalidNode;
        }
        return higraNodeIdBySlot_[static_cast<size_t>(nodeId)];
    }

    NodeId getNodeIdFromHigra(NodeId higraNodeId) const noexcept {
        if (!hasHigraNodeIdMapping() ||
            higraNodeId < 0 ||
            higraNodeId >= static_cast<NodeId>(slotByHigraNodeId_.size())) {
            return InvalidNode;
        }
        return slotByHigraNodeId_[static_cast<size_t>(higraNodeId)];
    }

    bool isNode(NodeId nodeId) const noexcept {
        return nodeId >= 0 && nodeId < static_cast<NodeId>(nodeParent_.size());
    }

    bool isProperPart(NodeId properPartId) const noexcept {
        return properPartId >= 0 && properPartId < static_cast<NodeId>(properPartOwner_.size());
    }

    bool isAlive(NodeId nodeId) const noexcept {
        return isNode(nodeId) && alive_[static_cast<size_t>(nodeId)] != 0;
    }

    bool isRoot(NodeId nodeId) const noexcept {
        return nodeId == rootNodeId_ && isAlive(nodeId) && getNodeParent(nodeId) == nodeId;
    }

    NodeId getNodeParent(NodeId nodeId) const noexcept {
        return isNode(nodeId) ? nodeParent_[static_cast<size_t>(nodeId)] : InvalidNode;
    }

    int getNumChildren(NodeId nodeId) const noexcept {
        return isAlive(nodeId) ? static_cast<int>(children_[static_cast<size_t>(nodeId)].size()) : 0;
    }

    const std::vector<NodeId>& getChildren(NodeId nodeId) const {
        static const std::vector<NodeId> empty;
        return isAlive(nodeId) ? children_[static_cast<size_t>(nodeId)] : empty;
    }

    NodeId getFirstChild(NodeId nodeId) const noexcept {
        if (!isAlive(nodeId) || children_[static_cast<size_t>(nodeId)].empty()) {
            return InvalidNode;
        }
        return children_[static_cast<size_t>(nodeId)].front();
    }

    NodeId getNextSibling(NodeId nodeId) const noexcept {
        if (!isAlive(nodeId)) {
            return InvalidNode;
        }
        const NodeId parentNodeId = getNodeParent(nodeId);
        if (!isAlive(parentNodeId) || parentNodeId == nodeId) {
            return InvalidNode;
        }
        const auto& siblings = children_[static_cast<size_t>(parentNodeId)];
        auto it = std::find(siblings.begin(), siblings.end(), nodeId);
        if (it == siblings.end() || ++it == siblings.end()) {
            return InvalidNode;
        }
        return *it;
    }

    bool hasChild(NodeId parentNodeId, NodeId childNodeId) const {
        if (!isAlive(parentNodeId)) {
            return false;
        }
        const auto& children = this->children_[static_cast<size_t>(parentNodeId)];
        return std::find(children.begin(), children.end(), childNodeId) != children.end();
    }

    bool isLeaf(NodeId nodeId) const noexcept {
        return isAlive(nodeId) && children_[static_cast<size_t>(nodeId)].empty();
    }

    int getNumLeafNodes() const {
        return static_cast<int>(getLeaves().size());
    }

    std::vector<NodeId> getLeaves() const {
        std::vector<NodeId> leaves;
        for (NodeId nodeId = 0; nodeId < static_cast<NodeId>(nodeParent_.size()); ++nodeId) {
            if (isLeaf(nodeId)) {
                leaves.push_back(nodeId);
            }
        }
        return leaves;
    }

    std::vector<NodeId> getLeafNodeIds() const {
        return getLeaves();
    }

    NodeId getSmallestComponent(NodeId properPartId) const noexcept {
        return isProperPart(properPartId) ? properPartOwner_[static_cast<size_t>(properPartId)] : InvalidNode;
    }

    int getNumProperParts(NodeId nodeId) const noexcept {
        return isAlive(nodeId) ? static_cast<int>(properPartsByNode_[static_cast<size_t>(nodeId)].size()) : 0;
    }

    const std::vector<NodeId>& getProperParts(NodeId nodeId) const {
        static const std::vector<NodeId> empty;
        return isAlive(nodeId) ? properPartsByNode_[static_cast<size_t>(nodeId)] : empty;
    }

    std::vector<NodeId> getAliveNodeIds() const {
        std::vector<NodeId> ids;
        ids.reserve(static_cast<size_t>(numNodes_));
        for (NodeId nodeId = 0; nodeId < static_cast<NodeId>(nodeParent_.size()); ++nodeId) {
            if (isAlive(nodeId)) {
                ids.push_back(nodeId);
            }
        }
        return ids;
    }

    int getNodeNumSiblings(NodeId nodeId) const noexcept {
        if (!isAlive(nodeId) || isRoot(nodeId)) {
            return 0;
        }
        const NodeId parentNodeId = getNodeParent(nodeId);
        return std::max(0, getNumChildren(parentNodeId) - 1);
    }

    std::vector<NodeId> getNodeSubtree(NodeId rootNodeId) const {
        if (!isAlive(rootNodeId)) {
            return {};
        }
        std::vector<NodeId> traversal;
        std::vector<NodeId> stack{rootNodeId};
        while (!stack.empty()) {
            const NodeId nodeId = stack.back();
            stack.pop_back();
            if (!isAlive(nodeId)) {
                continue;
            }
            traversal.push_back(nodeId);
            const auto& children = children_[static_cast<size_t>(nodeId)];
            for (auto it = children.rbegin(); it != children.rend(); ++it) {
                stack.push_back(*it);
            }
        }
        return traversal;
    }

    std::vector<NodeId> getDescendants(NodeId nodeId) const {
        auto subtree = getNodeSubtree(nodeId);
        if (!subtree.empty()) {
            subtree.erase(subtree.begin());
        }
        return subtree;
    }

    int getNodeNumDescendants(NodeId nodeId) const {
        return static_cast<int>(getDescendants(nodeId).size());
    }

    std::vector<NodeId> getPostOrderNodes() const {
        return getPostOrderNodes(rootNodeId_);
    }

    std::vector<NodeId> getPostOrderNodes(NodeId rootNodeId) const {
        if (!isAlive(rootNodeId)) {
            return {};
        }

        std::vector<NodeId> traversal;
        std::vector<std::pair<NodeId, bool>> stack;
        stack.emplace_back(rootNodeId, false);
        while (!stack.empty()) {
            auto [nodeId, expanded] = stack.back();
            stack.pop_back();
            if (!isAlive(nodeId)) {
                continue;
            }
            if (expanded) {
                traversal.push_back(nodeId);
                continue;
            }
            stack.emplace_back(nodeId, true);
            const auto& children = children_[static_cast<size_t>(nodeId)];
            for (auto it = children.rbegin(); it != children.rend(); ++it) {
                stack.emplace_back(*it, false);
            }
        }
        return traversal;
    }

    std::vector<NodeId> getIteratorBreadthFirstTraversal() const {
        return getIteratorBreadthFirstTraversal(rootNodeId_);
    }

    std::vector<NodeId> getIteratorBreadthFirstTraversal(NodeId rootNodeId) const {
        if (!isAlive(rootNodeId)) {
            return {};
        }
        std::vector<NodeId> traversal;
        std::queue<NodeId> queue;
        queue.push(rootNodeId);
        while (!queue.empty()) {
            const NodeId nodeId = queue.front();
            queue.pop();
            if (!isAlive(nodeId)) {
                continue;
            }
            traversal.push_back(nodeId);
            for (NodeId childNodeId : children_[static_cast<size_t>(nodeId)]) {
                queue.push(childNodeId);
            }
        }
        return traversal;
    }

    std::vector<NodeId> getPathToRootNodes(NodeId nodeId) const {
        if (!isAlive(nodeId)) {
            return {};
        }
        std::vector<NodeId> path;
        std::unordered_set<NodeId> seen;
        NodeId currentNodeId = nodeId;
        while (currentNodeId != InvalidNode && isAlive(currentNodeId)) {
            if (!seen.insert(currentNodeId).second) {
                return {};
            }
            path.push_back(currentNodeId);
            const NodeId parentNodeId = getNodeParent(currentNodeId);
            if (parentNodeId == currentNodeId) {
                break;
            }
            currentNodeId = parentNodeId;
        }
        return path;
    }

    std::vector<NodeId> getPathBetweenNodes(NodeId sourceNodeId, NodeId targetNodeId) const {
        if (!isAlive(sourceNodeId) || !isAlive(targetNodeId)) {
            return {};
        }
        if (sourceNodeId == targetNodeId) {
            return {sourceNodeId};
        }
        const NodeId lcaNodeId = getLowestCommonAncestor(sourceNodeId, targetNodeId);
        if (lcaNodeId == InvalidNode) {
            return {};
        }

        std::vector<NodeId> path;
        NodeId currentNodeId = sourceNodeId;
        while (currentNodeId != InvalidNode) {
            path.push_back(currentNodeId);
            if (currentNodeId == lcaNodeId) {
                break;
            }
            const NodeId parentNodeId = getNodeParent(currentNodeId);
            if (parentNodeId == currentNodeId) {
                break;
            }
            currentNodeId = parentNodeId;
        }
        if (path.empty() || path.back() != lcaNodeId) {
            return {};
        }

        std::vector<NodeId> targetTail;
        currentNodeId = targetNodeId;
        while (currentNodeId != InvalidNode && currentNodeId != lcaNodeId) {
            targetTail.push_back(currentNodeId);
            const NodeId parentNodeId = getNodeParent(currentNodeId);
            if (parentNodeId == currentNodeId) {
                return {};
            }
            currentNodeId = parentNodeId;
        }
        std::reverse(targetTail.begin(), targetTail.end());
        path.insert(path.end(), targetTail.begin(), targetTail.end());
        return path;
    }

    bool isAncestor(NodeId ancestorNodeId, NodeId nodeId) const {
        if (!isAlive(ancestorNodeId) || !isAlive(nodeId)) {
            return false;
        }
        NodeId currentNodeId = nodeId;
        std::unordered_set<NodeId> seen;
        while (currentNodeId != InvalidNode && isAlive(currentNodeId)) {
            if (currentNodeId == ancestorNodeId) {
                return true;
            }
            if (!seen.insert(currentNodeId).second) {
                return false;
            }
            const NodeId parentNodeId = getNodeParent(currentNodeId);
            if (parentNodeId == currentNodeId) {
                return false;
            }
            currentNodeId = parentNodeId;
        }
        return false;
    }

    bool isDescendant(NodeId nodeId, NodeId ancestorNodeId) const {
        return isAncestor(ancestorNodeId, nodeId);
    }

    bool isStrictAncestor(NodeId ancestorNodeId, NodeId nodeId) const {
        return ancestorNodeId != nodeId && isAncestor(ancestorNodeId, nodeId);
    }

    bool isStrictDescendant(NodeId nodeId, NodeId ancestorNodeId) const {
        return nodeId != ancestorNodeId && isDescendant(nodeId, ancestorNodeId);
    }

    bool isComparable(NodeId lhs, NodeId rhs) const {
        return isAncestor(lhs, rhs) || isAncestor(rhs, lhs);
    }

    bool isStrictComparable(NodeId lhs, NodeId rhs) const {
        return lhs != rhs && isComparable(lhs, rhs);
    }

    NodeId getLowestCommonAncestor(NodeId lhs, NodeId rhs) const {
        if (!isAlive(lhs) || !isAlive(rhs)) {
            return InvalidNode;
        }
        std::unordered_set<NodeId> lhsAncestors;
        NodeId currentNodeId = lhs;
        while (currentNodeId != InvalidNode && isAlive(currentNodeId)) {
            lhsAncestors.insert(currentNodeId);
            const NodeId parentNodeId = getNodeParent(currentNodeId);
            if (parentNodeId == currentNodeId) {
                break;
            }
            currentNodeId = parentNodeId;
        }

        currentNodeId = rhs;
        std::unordered_set<NodeId> seen;
        while (currentNodeId != InvalidNode && isAlive(currentNodeId)) {
            if (lhsAncestors.find(currentNodeId) != lhsAncestors.end()) {
                return currentNodeId;
            }
            if (!seen.insert(currentNodeId).second) {
                return InvalidNode;
            }
            const NodeId parentNodeId = getNodeParent(currentNodeId);
            if (parentNodeId == currentNodeId) {
                break;
            }
            currentNodeId = parentNodeId;
        }
        return InvalidNode;
    }

    int getNodeTimePreOrder(NodeId nodeId) const {
        if (!traversalTimesValid_) {
            rebuildTraversalTimes();
        }
        return isNode(nodeId) ? preOrderTime_[static_cast<size_t>(nodeId)] : -1;
    }

    int getNodeTimePostOrder(NodeId nodeId) const {
        if (!traversalTimesValid_) {
            rebuildTraversalTimes();
        }
        return isNode(nodeId) ? postOrderTime_[static_cast<size_t>(nodeId)] : -1;
    }

    bool hasDetachedAliveNodes() const noexcept {
        for (NodeId nodeId = 0; nodeId < static_cast<NodeId>(nodeParent_.size()); ++nodeId) {
            if (isAlive(nodeId) && nodeId != rootNodeId_ && nodeParent_[static_cast<size_t>(nodeId)] == nodeId) {
                return true;
            }
        }
        return false;
    }

    void validateConnectedRootedTree() const {
        if (!isAlive(rootNodeId_)) {
            throw std::runtime_error("A connected tree requires a live root.");
        }
        if (getNodeParent(rootNodeId_) != rootNodeId_) {
            throw std::runtime_error("The connected root must be self-parented.");
        }

        int selfParentCount = 0;
        std::vector<int> childReferenceCount(nodeParent_.size(), 0);
        for (NodeId nodeId = 0; nodeId < static_cast<NodeId>(nodeParent_.size()); ++nodeId) {
            if (!isAlive(nodeId)) {
                continue;
            }
            const NodeId parentNodeId = getNodeParent(nodeId);
            if (parentNodeId == nodeId) {
                ++selfParentCount;
                if (nodeId != rootNodeId_) {
                    throw std::runtime_error("Connected validation found a detached self-parented node.");
                }
            } else if (!isAlive(parentNodeId)) {
                throw std::runtime_error("Connected validation found a node with an invalid parent.");
            }

            std::unordered_set<NodeId> seen;
            NodeId currentNodeId = nodeId;
            while (currentNodeId != InvalidNode) {
                if (!isAlive(currentNodeId) || !seen.insert(currentNodeId).second) {
                    throw std::runtime_error("Connected validation found an invalid parent chain.");
                }
                const NodeId parentId = getNodeParent(currentNodeId);
                if (parentId == currentNodeId) {
                    if (currentNodeId != rootNodeId_) {
                        throw std::runtime_error("Every live node must reach the connected root.");
                    }
                    break;
                }
                currentNodeId = parentId;
            }

            for (NodeId childNodeId : children_[static_cast<size_t>(nodeId)]) {
                if (!isAlive(childNodeId) || getNodeParent(childNodeId) != nodeId) {
                    throw std::runtime_error("Connected validation found inconsistent child links.");
                }
                ++childReferenceCount[static_cast<size_t>(childNodeId)];
            }
        }

        if (selfParentCount != 1) {
            throw std::runtime_error("Connected validation requires exactly one self-parented root.");
        }
        for (NodeId nodeId = 0; nodeId < static_cast<NodeId>(nodeParent_.size()); ++nodeId) {
            if (!isAlive(nodeId)) {
                continue;
            }
            const int expectedReferences = nodeId == rootNodeId_ ? 0 : 1;
            if (childReferenceCount[static_cast<size_t>(nodeId)] != expectedReferences) {
                throw std::runtime_error("Connected validation found inconsistent child reference counts.");
            }
        }

        for (NodeId ownerNodeId : properPartOwner_) {
            if (!isAlive(ownerNodeId)) {
                throw std::runtime_error("Every proper part must be owned by a live node.");
            }
        }
    }

    void pruneNode(NodeId nodeId) {
        ensureAliveNodeId(nodeId, "pruneNode requires a live node.");
        if (isRoot(nodeId)) {
            throw std::invalid_argument("pruneNode cannot prune the root.");
        }
        const NodeId parentNodeId = getNodeParent(nodeId);
        if (!isAlive(parentNodeId) || parentNodeId == nodeId) {
            throw std::invalid_argument("pruneNode requires an attached non-root node.");
        }

        const auto subtree = getPostOrderNodes(nodeId);
        std::unordered_set<NodeId> subtreeNodes(subtree.begin(), subtree.end());
        for (NodeId& ownerNodeId : properPartOwner_) {
            if (subtreeNodes.find(ownerNodeId) != subtreeNodes.end()) {
                ownerNodeId = parentNodeId;
            }
        }
        rebuildProperPartsFromOwnership();

        for (NodeId subtreeNodeId : subtree) {
            const NodeId currentParentId = getNodeParent(subtreeNodeId);
            if (currentParentId != InvalidNode && currentParentId != subtreeNodeId) {
                removeChildLink(currentParentId, subtreeNodeId);
            }
            children_[static_cast<size_t>(subtreeNodeId)].clear();
            nodeParent_[static_cast<size_t>(subtreeNodeId)] = subtreeNodeId;
            releaseNode(subtreeNodeId);
        }
        markProperPartChanged();
        markTopologyChanged();
    }

    void mergeNodeIntoParent(NodeId nodeId) {
        ensureAliveNodeId(nodeId, "mergeNodeIntoParent requires a live node.");
        if (isRoot(nodeId)) {
            throw std::invalid_argument("mergeNodeIntoParent cannot merge the root.");
        }
        const NodeId parentNodeId = getNodeParent(nodeId);
        if (!isAlive(parentNodeId) || parentNodeId == nodeId) {
            throw std::invalid_argument("mergeNodeIntoParent requires an attached non-root node.");
        }

        for (NodeId& ownerNodeId : properPartOwner_) {
            if (ownerNodeId == nodeId) {
                ownerNodeId = parentNodeId;
            }
        }
        rebuildProperPartsFromOwnership();

        auto movingChildren = std::move(children_[static_cast<size_t>(nodeId)]);
        children_[static_cast<size_t>(nodeId)].clear();

        auto& siblings = children_[static_cast<size_t>(parentNodeId)];
        auto position = std::find(siblings.begin(), siblings.end(), nodeId);
        const auto insertOffset = static_cast<std::ptrdiff_t>(std::distance(siblings.begin(), position));
        if (position != siblings.end()) {
            siblings.erase(position);
        }
        auto insertPosition = siblings.begin() + std::min<std::ptrdiff_t>(insertOffset, static_cast<std::ptrdiff_t>(siblings.size()));
        for (NodeId childNodeId : movingChildren) {
            nodeParent_[static_cast<size_t>(childNodeId)] = parentNodeId;
            insertPosition = siblings.insert(insertPosition, childNodeId);
            ++insertPosition;
        }

        nodeParent_[static_cast<size_t>(nodeId)] = nodeId;
        releaseNode(nodeId);
        markProperPartChanged();
        markTopologyChanged();
    }

    std::pair<std::vector<NodeId>, std::vector<NodeId>> computeAscendantsAndDescendants(int delta) const {
        if (delta < 0) {
            throw std::invalid_argument("delta must be non-negative.");
        }

        std::vector<NodeId> ascendants(nodeParent_.size(), InvalidNode);
        std::vector<NodeId> descendants(nodeParent_.size(), InvalidNode);
        const std::vector<int32_t> areaByNode = computeAreasIncrementally();

        for (NodeId nodeId : getAliveNodeIds()) {
            NodeId currentNodeId = nodeId;
            bool valid = true;
            for (int step = 0; step < delta; ++step) {
                if (isRoot(currentNodeId)) {
                    valid = false;
                    break;
                }
                currentNodeId = getNodeParent(currentNodeId);
                if (!isAlive(currentNodeId)) {
                    valid = false;
                    break;
                }
            }
            if (!valid) {
                continue;
            }
            ascendants[static_cast<size_t>(nodeId)] = currentNodeId;
            maxAreaDescendants(areaByNode, descendants, currentNodeId, nodeId);
        }

        return {std::move(ascendants), std::move(descendants)};
    }
};

} // namespace mmcfilters
