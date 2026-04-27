#pragma once

#include "MorphologicalTree.hpp"
#include "TreeEditor.hpp"
#include "detail/TreeAltitudeOpsDetail.hpp"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace mmcfilters {

class AttributeComputer;
class AttributeComputedIncrementally;
class AttributeFilters;
class AttributeOpeningPrimitivesFamily;
class ComputerMSER;
class ExtinctionValues;
class UltimateAttributeOpening;
class WeightedTreeEditor;

/**
 * @brief Wrapper pairing `MorphologicalTree` topology with an external altitude buffer.
 *
 * `WeightedMorphologicalTree` owns the mutable topology and node-altitude state
 * internally, indexed by dense internal `NodeId`. Read-only topology access and
 * weighted convenience methods are exposed without leaking structural mutation.
 * Staged structural edits must go through `edit()`, which returns a
 * `WeightedTreeEditor`.
 *
 * `createFromHigraParent()` imports a static Higra hierarchy and preserves its
 * original node-id domain until the topology is edited. `exportHigraHierarchy()`
 * creates a new compact Higra domain for the current live rooted tree.
 */
class WeightedMorphologicalTree {
    friend class WeightedTreeEditor;
    friend class AttributeComputer;
    friend class AttributeComputedIncrementally;
    friend class AttributeFilters;
    friend class AttributeOpeningPrimitivesFamily;
    friend class ComputerMSER;
    friend class ExtinctionValues;
    friend class UltimateAttributeOpening;

private:
    MorphologicalTree tree_;
    AltitudeBuffer altitude_;

    void configureEmptyTopology(int rows, int cols, int treeType, std::optional<AdjacencyRelation> adjacency, NodeId numProperParts) {
        tree_.treeType_ = treeType;
        tree_.numRows_ = rows;
        tree_.numCols_ = cols;
        tree_.adj_ = std::move(adjacency);
        tree_.tosAdjacencyPolicy_ = std::nullopt;
        tree_.initializeEmptyStorage(static_cast<size_t>(numProperParts));
        altitude_.clear();
    }

    void assignAltitudeFromDirectProperParts(const ImageUInt8Ptr& img) {
        altitude_.assign(static_cast<size_t>(tree_.getNumInternalNodeSlots()), AltitudeType{});
        for (NodeId nodeId : tree_.getAliveNodeIds()) {
            const auto properParts = tree_.getProperParts(nodeId);
            auto it = properParts.begin();
            if (it == properParts.end()) {
                throw std::runtime_error("Cannot infer node altitude from a topology node without direct proper parts.");
            }
            altitude_[static_cast<size_t>(nodeId)] = static_cast<AltitudeType>((*img)[*it]);
        }
    }

    void assignInternalAltitude(std::span<const AltitudeType> altitudeValues) {
        if (altitudeValues.size() != static_cast<size_t>(tree_.getNumInternalNodeSlots())) {
            throw std::invalid_argument("Internal altitude buffer size must match the dense internal-node domain.");
        }
        altitude_.assign(altitudeValues.begin(), altitudeValues.end());
    }

    void importAltitudeFromHigra(std::span<const AltitudeType> higraAltitude) {
        if (static_cast<size_t>(tree_.getNumHigraNodes()) != higraAltitude.size()) {
            throw std::invalid_argument("Higra altitude buffer size must match the preserved imported Higra hierarchy.");
        }
        altitude_.assign(static_cast<size_t>(tree_.getNumInternalNodeSlots()), AltitudeType{});
        for (NodeId slotId = 0; slotId < tree_.getNumInternalNodeSlots(); ++slotId) {
            const NodeId higraNodeId = tree_.getHigraNodeId(slotId);
            altitude_[static_cast<size_t>(slotId)] = static_cast<AltitudeType>(higraAltitude[static_cast<size_t>(higraNodeId)]);
        }
    }

private:
    WeightedMorphologicalTree() = default;

    WeightedMorphologicalTree(MorphologicalTree&& topology, AltitudeBuffer altitudeBuffer) : tree_(std::move(topology)), altitude_(std::move(altitudeBuffer)) {
        validateAltitudeBufferShape();
    }

    WeightedMorphologicalTree(
        ImageUInt8Ptr img,
        ToSInterpolation interpolation = ToSInterpolation::SelfDual,
        int infinitySeedRow = ToSDefaultInfinityRow,
        int infinitySeedCol = ToSDefaultInfinityCol) {
        MorphologicalTree::requireNonEmptyImageDomain(img, "WeightedMorphologicalTree::createTreeOfShapes");
        configureEmptyTopology(
            img->getNumRows(),
            img->getNumCols(),
            MorphologicalTree::TREE_OF_SHAPES,
            std::nullopt,
            static_cast<NodeId>(img->getSize()));
        tree_.tosAdjacencyPolicy_ = MorphologicalTree::treeOfShapesAdjacencyPolicy(interpolation, tree_.numRows_, tree_.numCols_);
        BuilderTreeOfShape builderUF(interpolation, infinitySeedRow, infinitySeedCol);
        tree_.build(img, builderUF);
        assignAltitudeFromDirectProperParts(img);
        validateAltitudeBufferShape();
    }

    explicit WeightedMorphologicalTree(ImageUInt8Ptr img, bool isMaxtree, double radius = 1.5) {
        MorphologicalTree::requireNonEmptyImageDomain(img, "WeightedMorphologicalTree::createComponentTree");
        configureEmptyTopology(
            img->getNumRows(),
            img->getNumCols(),
            isMaxtree ? MorphologicalTree::MAX_TREE : MorphologicalTree::MIN_TREE,
            std::optional<AdjacencyRelation>(std::in_place, img->getNumRows(), img->getNumCols(), radius),
            static_cast<NodeId>(img->getSize()));
        BuilderComponentTree builderUF(&*tree_.adj_, isMaxtree);
        tree_.build(img, builderUF);
        assignAltitudeFromDirectProperParts(img);
        validateAltitudeBufferShape();
    }

public:
    WeightedMorphologicalTree(const WeightedMorphologicalTree&) = delete;
    WeightedMorphologicalTree& operator=(const WeightedMorphologicalTree&) = delete;
    WeightedMorphologicalTree(WeightedMorphologicalTree&&) noexcept = default;
    WeightedMorphologicalTree& operator=(WeightedMorphologicalTree&&) noexcept = default;

    static WeightedMorphologicalTree createComponentTree(ImageUInt8Ptr img, bool isMaxtree, double radius = 1.5) {
        return WeightedMorphologicalTree(img, isMaxtree, radius);
    }

    static WeightedMorphologicalTree createTreeOfShapes(
        ImageUInt8Ptr img,
        ToSInterpolation interpolation = ToSInterpolation::SelfDual,
        int infinitySeedRow = ToSDefaultInfinityRow,
        int infinitySeedCol = ToSDefaultInfinityCol) {
        return WeightedMorphologicalTree(img, interpolation, infinitySeedRow, infinitySeedCol);
    }

    /**
     * @brief Imports topology and altitude from a static Higra parent/altitude representation.
     */
    static WeightedMorphologicalTree createFromHigraParent(
        std::span<const NodeId> higraParent,
        std::span<const AltitudeType> higraAltitude,
        int rows,
        int cols,
        int treeType,
        std::optional<AdjacencyRelation> adjacency = std::nullopt) {
        WeightedMorphologicalTree weighted;
        weighted.tree_ = MorphologicalTree::createFromHigraParent(
            higraParent,
            rows,
            cols,
            treeType,
            std::move(adjacency));
        weighted.importAltitudeFromHigra(higraAltitude);
        return weighted;
    }

    /**
     * @brief Returns the explicit altitude buffer required by static weighted operations.
     */
    static const AltitudeBuffer& requireAltitudeBuffer(const AltitudeBuffer* altitude) {
        if (altitude == nullptr) {
            throw std::logic_error("This operation requires an explicit altitude buffer. Use WeightedMorphologicalTree or provide an explicit altitude buffer.");
        }
        return *altitude;
    }

    /**
     * @brief Validates that an altitude buffer covers the dense internal-node domain.
     */
    static void validateAltitudeBufferShape(const MorphologicalTree& tree, std::span<const AltitudeType> altitude) {
        if (altitude.size() != static_cast<size_t>(tree.getNumInternalNodeSlots())) {
            throw std::runtime_error("Altitude buffer size must match the dense internal-node domain.");
        }
    }

    static void validateAltitudeBufferShape(const MorphologicalTree& tree, const AltitudeBuffer* altitude) {
        validateAltitudeBufferShape(tree, std::span<const AltitudeType>(requireAltitudeBuffer(altitude)));
    }

    /**
     * @brief Reads one node altitude from an explicit altitude buffer.
     */
    static AltitudeType getAltitude(std::span<const AltitudeType> altitude, NodeId nodeId) {
        if (nodeId < 0 || static_cast<size_t>(nodeId) >= altitude.size()) {
            throw std::invalid_argument("Altitude access requires a valid internal NodeId.");
        }
        return altitude[static_cast<size_t>(nodeId)];
    }

    static AltitudeType getAltitude(const AltitudeBuffer* altitude, NodeId nodeId) {
        return getAltitude(std::span<const AltitudeType>(requireAltitudeBuffer(altitude)), nodeId);
    }

    /**
     * @brief Computes the altitude difference between one node and its parent.
     */
    static AltitudeDiffType getNodeResidue(const MorphologicalTree& tree, std::span<const AltitudeType> altitude, NodeId nodeId) {
        validateAltitudeBufferShape(tree, altitude);
        if (!tree.isAlive(nodeId)) {
            throw std::invalid_argument("Node residue requires a live internal NodeId.");
        }
        const NodeId parentNodeId = tree.getNodeParent(nodeId);
        if (parentNodeId == InvalidNode || parentNodeId == nodeId) {
            return getAltitude(altitude, nodeId);
        }
        return getAltitude(altitude, nodeId) - getAltitude(altitude, parentNodeId);
    }

    static AltitudeDiffType getNodeResidue(const MorphologicalTree& tree, const AltitudeBuffer* altitude, NodeId nodeId) {
        return getNodeResidue(tree, std::span<const AltitudeType>(requireAltitudeBuffer(altitude)), nodeId);
    }

    /**
     * @brief Finds the first ascendant at distance `delta` in altitude space.
     */
    static NodeId getNodeAscendant(const MorphologicalTree& tree, std::span<const AltitudeType> altitude, NodeId nodeId, int delta) {
        validateAltitudeBufferShape(tree, altitude);
        if (!tree.isAlive(nodeId)) {
            throw std::invalid_argument("Node ascendant search requires a live internal NodeId.");
        }
        NodeId currentNodeId = nodeId;
        while (true) {
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
    }

    /**
     * @brief Computes delta ascendants and the largest-area descendant attached to each ascendant.
     */
    static std::pair<std::vector<NodeId>, std::vector<NodeId>> computeAscendantsAndDescendants(
        const MorphologicalTree& tree,
        const AltitudeBuffer* altitude,
        int delta) {
        std::vector<NodeId> ascendants(static_cast<size_t>(tree.getNumInternalNodeSlots()), InvalidNode);
        std::vector<NodeId> descendants(static_cast<size_t>(tree.getNumInternalNodeSlots()), InvalidNode);
        validateAltitudeBufferShape(tree, altitude);
        const std::vector<int32_t> areaByNode = detail::tree_altitude_ops::computeAreasIncrementally(tree);
        const AltitudeBuffer& altitudeBuffer = requireAltitudeBuffer(altitude);

        for (NodeId nodeId : tree.getAliveNodeIds()) {
            const NodeId ascendantNodeId = getNodeAscendant(tree, std::span<const AltitudeType>(altitudeBuffer), nodeId, delta);
            if (ascendantNodeId == InvalidNode) {
                continue;
            }
            detail::tree_altitude_ops::maxAreaDescendants(tree, areaByNode, descendants, ascendantNodeId, nodeId);
            if (descendants[static_cast<size_t>(ascendantNodeId)] != InvalidNode) {
                ascendants[static_cast<size_t>(nodeId)] = ascendantNodeId;
            }
        }

        return {std::move(ascendants), std::move(descendants)};
    }

    /**
     * @brief Reconstructs an image from topology ownership and explicit node altitudes.
     */
    static ImageUInt8Ptr reconstructImage(const MorphologicalTree& tree, std::span<const AltitudeType> altitude) {
        validateAltitudeBufferShape(tree, altitude);
        ImageUInt8Ptr image = ImageUInt8::create(tree.getNumRowsOfImage(), tree.getNumColsOfImage());
        auto imgBuffer = image->rawData();
        for (int pixelId = 0; pixelId < tree.getNumTotalProperParts(); ++pixelId) {
            const NodeId nodeId = tree.getSmallestComponent(pixelId);
            imgBuffer[static_cast<size_t>(pixelId)] = static_cast<uint8_t>(getAltitude(altitude, nodeId));
        }
        return image;
    }

    static ImageUInt8Ptr reconstructImage(const MorphologicalTree& tree, const AltitudeBuffer* altitude) {
        return reconstructImage(tree, std::span<const AltitudeType>(requireAltitudeBuffer(altitude)));
    }

    /**
     * @brief Exports a live rooted topology and explicit altitudes to a compact parent/altitude representation.
     */
    static std::pair<std::vector<NodeId>, std::vector<AltitudeType>> exportHigraHierarchy(
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
            throw std::runtime_error("Cannot export a forest or a tree with detached alive nodes to a compact Higra representation.");
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
                return tree.getNodeTimePostOrder(lhs) < tree.getNodeTimePostOrder(rhs);
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
                throw std::runtime_error("Each proper part must belong to one alive node when exporting a compact Higra hierarchy.");
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

    static std::pair<std::vector<NodeId>, std::vector<AltitudeType>> exportHigraHierarchy(
        const MorphologicalTree& tree,
        const AltitudeBuffer* altitude) {
        return exportHigraHierarchy(tree, std::span<const AltitudeType>(requireAltitudeBuffer(altitude)));
    }

    /**
     * @brief Validates altitude monotonicity for max-trees and min-trees.
     */
    static void validateMonotoneAltitude(const MorphologicalTree& tree, std::span<const AltitudeType> altitude) {
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

    static void validateMonotoneAltitude(const MorphologicalTree& tree, const AltitudeBuffer* altitude) {
        validateMonotoneAltitude(tree, std::span<const AltitudeType>(requireAltitudeBuffer(altitude)));
    }

    void validateAltitudeBufferShape() const {
        validateAltitudeBufferShape(tree_, altitude_);
    }

    /**
     * @brief Returns read-only access to the owned topology.
     *
     * This accessor intentionally does not provide a mutable `MorphologicalTree`
     * handle. Use `edit()` for staged topology edits or the safe public
     * mutators for complete local changes.
     */
    const MorphologicalTree& topology() const noexcept {
        return tree_;
    }

    const AltitudeBuffer& getAltitudeBuffer() const noexcept {
        return altitude_;
    }

    void setAltitudeBuffer(AltitudeBuffer altitudeBuffer) {
        if (altitudeBuffer.size() != static_cast<size_t>(tree_.getNumInternalNodeSlots())) {
            throw std::runtime_error("Altitude buffer size must match the dense internal-node domain.");
        }
        altitude_ = std::move(altitudeBuffer);
    }

    void validateMonotoneAltitude() const {
        validateMonotoneAltitude(tree_, altitude_);
    }

    AltitudeType getAltitude(NodeId nodeId) const {
        if (!tree_.isAlive(nodeId) || static_cast<size_t>(nodeId) >= altitude_.size()) {
            throw std::invalid_argument("WeightedMorphologicalTree::getAltitude requires a live internal NodeId.");
        }
        return altitude_[static_cast<size_t>(nodeId)];
    }

    void setAltitude(NodeId nodeId, AltitudeType value) {
        if (!tree_.isAlive(nodeId) || static_cast<size_t>(nodeId) >= altitude_.size()) {
            throw std::invalid_argument("WeightedMorphologicalTree::setAltitude requires a live internal NodeId.");
        }
        altitude_[static_cast<size_t>(nodeId)] = value;
    }

    /**
     * @brief Prunes a complete subtree through the owned topology.
     *
     * The altitude buffer is kept as the canonical dense weighted state; dead
     * slots may retain their old values until a compact export is requested.
     */
    void pruneNode(NodeId nodeId) {
        tree_.pruneNode(nodeId);
    }

    /**
     * @brief Merges one node into its parent through the owned topology.
     *
     * This is the weighted counterpart of the safe public topology mutator.
     */
    void mergeNodeIntoParent(NodeId nodeId) {
        tree_.mergeNodeIntoParent(nodeId);
    }

    AltitudeDiffType getNodeResidue(NodeId nodeId) const {
        return getNodeResidue(tree_, altitude_, nodeId);
    }

    ImageUInt8Ptr reconstructionImage() const {
        return reconstructImage(tree_, altitude_);
    }

    /**
     * @brief Exports the current live rooted tree to a new compact Higra parent/altitude representation.
     */
    std::pair<std::vector<NodeId>, std::vector<AltitudeType>> exportHigraHierarchy() const {
        return exportHigraHierarchy(tree_, altitude_);
    }

    /**
     * @brief Opens the only public entrypoint for staged weighted edits.
     */
    WeightedTreeEditor edit();
};

/**
 * @brief Edit-session facade for `WeightedMorphologicalTree`.
 *
 * The weighted editor reuses the structural `TreeEditor` for topology while
 * keeping the external altitude buffer as the canonical weighted state.
 * `commit()` first validates the topology and then validates monotone altitude.
 */
class WeightedTreeEditor {
    friend class WeightedMorphologicalTree;

private:
    WeightedMorphologicalTree& weighted_;
    TreeEditor editor_;

    explicit WeightedTreeEditor(WeightedMorphologicalTree& weighted)
        : weighted_(weighted), editor_(weighted.tree_.edit()) {}

public:
    NodeId createDetachedNode(AltitudeType altitude = AltitudeType{}) {
        const NodeId nodeId = editor_.createDetachedNode();
        weighted_.altitude_.resize(static_cast<size_t>(weighted_.tree_.getNumInternalNodeSlots()), AltitudeType{});
        weighted_.altitude_[static_cast<size_t>(nodeId)] = altitude;
        return nodeId;
    }

    void setNodeAltitude(NodeId nodeId, AltitudeType altitude) {
        if (!weighted_.tree_.isAlive(nodeId)) {
            throw std::invalid_argument("WeightedTreeEditor::setNodeAltitude requires a live node.");
        }
        weighted_.altitude_[static_cast<size_t>(nodeId)] = altitude;
    }

    void detach(NodeId nodeId) {
        editor_.detach(nodeId);
    }

    void reparent(NodeId nodeId, NodeId newParentId) {
        editor_.reparent(nodeId, newParentId);
    }

    void attach(NodeId parentId, NodeId detachedNodeId) {
        editor_.attach(parentId, detachedNodeId);
    }

    void moveChildren(NodeId parentId, NodeId sourceId) {
        editor_.moveChildren(parentId, sourceId);
    }

    void moveProperPart(NodeId targetNodeId, NodeId sourceNodeId, NodeId properPartId) {
        editor_.moveProperPart(targetNodeId, sourceNodeId, properPartId);
    }

    void moveProperParts(NodeId targetNodeId, NodeId sourceNodeId) {
        editor_.moveProperParts(targetNodeId, sourceNodeId);
    }

    void removeChild(NodeId parentNodeId, NodeId childId, bool releaseNodeFlag) {
        editor_.removeChild(parentNodeId, childId, releaseNodeFlag);
    }

    void releaseNode(NodeId nodeId) {
        editor_.releaseNode(nodeId);
    }

    void setRoot(NodeId nodeId) {
        editor_.setRoot(nodeId);
    }

    bool hasDetachedAliveNodes() const noexcept {
        return editor_.hasDetachedAliveNodes();
    }

    void commit() {
        editor_.commit();
        weighted_.validateMonotoneAltitude();
    }

    void commitUnchecked() noexcept {
        editor_.commitUnchecked();
    }
};

inline WeightedTreeEditor WeightedMorphologicalTree::edit() {
    return WeightedTreeEditor(*this);
}

} // namespace mmcfilters
