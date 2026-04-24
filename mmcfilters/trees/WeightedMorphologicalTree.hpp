#pragma once

#include "MorphologicalTree.hpp"
#include "TreeAltitudeOps.hpp"
#include "TreeEditor.hpp"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace mmcfilters {

class WeightedTreeEditor;

/**
 * @brief Wrapper pairing `MorphologicalTree` topology with an external altitude buffer.
 *
 * `WeightedMorphologicalTree` keeps the mutable topology in `tree` and the
 * node-altitude state in `altitude`, indexed by dense internal `NodeId`. The
 * wrapper offers the weighted conveniences that are intentionally kept out of
 * `MorphologicalTree`.
 *
 * `createFromHigraParent()` imports a static Higra hierarchy and preserves its
 * original node-id domain until the topology is edited. `exportHigraHierarchy()`
 * creates a new compact Higra domain for the current live rooted tree.
 */
class WeightedMorphologicalTree {
    friend class WeightedTreeEditor;

    void configureEmptyTopology(int rows, int cols, int treeType, std::optional<AdjacencyRelation> adjacency, NodeId numProperParts) {
        tree.treeType_ = treeType;
        tree.numRows_ = rows;
        tree.numCols_ = cols;
        tree.adj_ = std::move(adjacency);
        tree.initializeEmptyStorage(static_cast<size_t>(numProperParts));
        altitude.clear();
    }

    void assignAltitudeFromDirectProperParts(const ImageUInt8Ptr& img) {
        altitude.assign(static_cast<size_t>(tree.getNumInternalNodeSlots()), AltitudeType{});
        for (NodeId nodeId : tree.getAliveNodeIds()) {
            const auto properParts = tree.getProperParts(nodeId);
            auto it = properParts.begin();
            if (it == properParts.end()) {
                throw std::runtime_error("Cannot infer node altitude from a topology node without direct proper parts.");
            }
            altitude[static_cast<size_t>(nodeId)] = static_cast<AltitudeType>((*img)[*it]);
        }
    }

    void assignInternalAltitude(std::span<const AltitudeType> altitudeValues) {
        if (altitudeValues.size() != static_cast<size_t>(tree.getNumInternalNodeSlots())) {
            throw std::invalid_argument("Internal altitude buffer size must match the dense internal-node domain.");
        }
        altitude.assign(altitudeValues.begin(), altitudeValues.end());
    }

    void importAltitudeFromHigra(std::span<const AltitudeType> higraAltitude) {
        if (static_cast<size_t>(tree.getNumHigraNodes()) != higraAltitude.size()) {
            throw std::invalid_argument("Higra altitude buffer size must match the preserved imported Higra hierarchy.");
        }
        altitude.assign(static_cast<size_t>(tree.getNumInternalNodeSlots()), AltitudeType{});
        for (NodeId slotId = 0; slotId < tree.getNumInternalNodeSlots(); ++slotId) {
            const NodeId higraNodeId = tree.getHigraNodeId(slotId);
            altitude[static_cast<size_t>(slotId)] = static_cast<AltitudeType>(higraAltitude[static_cast<size_t>(higraNodeId)]);
        }
    }

public:
    MorphologicalTree tree;
    AltitudeBuffer altitude;

private:
    WeightedMorphologicalTree() = default;

    WeightedMorphologicalTree(MorphologicalTree&& topology, AltitudeBuffer altitudeBuffer) : tree(std::move(topology)), altitude(std::move(altitudeBuffer)) {
        validateAltitudeBufferShape();
    }

    WeightedMorphologicalTree(ImageUInt8Ptr img, ToSInterpolation interpolation = ToSInterpolation::SelfDual) {
        configureEmptyTopology(
            img->getNumRows(),
            img->getNumCols(),
            MorphologicalTree::TREE_OF_SHAPES,
            std::nullopt,
            static_cast<NodeId>(img->getSize()));
        BuilderTreeOfShape builderUF(interpolation == ToSInterpolation::Min4cMax8c);
        tree.build(img, builderUF);
        assignAltitudeFromDirectProperParts(img);
        validateAltitudeBufferShape();
    }

    explicit WeightedMorphologicalTree(ImageUInt8Ptr img, bool isMaxtree, double radius = 1.5) {
        configureEmptyTopology(
            img->getNumRows(),
            img->getNumCols(),
            isMaxtree ? MorphologicalTree::MAX_TREE : MorphologicalTree::MIN_TREE,
            std::optional<AdjacencyRelation>(std::in_place, img->getNumRows(), img->getNumCols(), radius),
            static_cast<NodeId>(img->getSize()));
        BuilderComponentTree builderUF(&*tree.adj_, isMaxtree);
        tree.build(img, builderUF);
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

    static WeightedMorphologicalTree createTreeOfShapes(ImageUInt8Ptr img, ToSInterpolation interpolation = ToSInterpolation::SelfDual) {
        return WeightedMorphologicalTree(img, interpolation);
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
        weighted.tree = MorphologicalTree::createFromHigraParent(
            higraParent,
            rows,
            cols,
            treeType,
            std::move(adjacency));
        weighted.importAltitudeFromHigra(higraAltitude);
        return weighted;
    }

    void validateAltitudeBufferShape() const {
        tree_altitude_ops::validateAltitudeBufferShape(tree, altitude);
    }

    const AltitudeBuffer& getAltitudeBuffer() const noexcept {
        return altitude;
    }

    void setAltitudeBuffer(AltitudeBuffer altitudeBuffer) {
        if (altitudeBuffer.size() != static_cast<size_t>(tree.getNumInternalNodeSlots())) {
            throw std::runtime_error("Altitude buffer size must match the dense internal-node domain.");
        }
        altitude = std::move(altitudeBuffer);
    }

    void validateMonotoneAltitude() const {
        tree_altitude_ops::validateMonotoneAltitude(tree, altitude);
    }

    AltitudeType getAltitude(NodeId nodeId) const {
        if (!tree.isNode(nodeId) || static_cast<size_t>(nodeId) >= altitude.size()) {
            throw std::invalid_argument("WeightedMorphologicalTree::getAltitude requires a valid internal NodeId.");
        }
        return altitude[static_cast<size_t>(nodeId)];
    }

    void setAltitude(NodeId nodeId, AltitudeType value) {
        if (!tree.isNode(nodeId) || static_cast<size_t>(nodeId) >= altitude.size()) {
            throw std::invalid_argument("WeightedMorphologicalTree::setAltitude requires a valid internal NodeId.");
        }
        altitude[static_cast<size_t>(nodeId)] = value;
    }

    void pruneNode(NodeId nodeId) {
        tree.pruneNode(nodeId);
    }

    void mergeNodeIntoParent(NodeId nodeId) {
        tree.mergeNodeIntoParent(nodeId);
    }

    AltitudeDiffType getNodeResidue(NodeId nodeId) const noexcept {
        return tree_altitude_ops::getNodeResidue(tree, altitude, nodeId);
    }

    ImageUInt8Ptr reconstructionImage() const {
        return tree_altitude_ops::reconstructImage(tree, altitude);
    }

    /**
     * @brief Exports the current live rooted tree to a new compact Higra parent/altitude representation.
     */
    std::pair<std::vector<NodeId>, std::vector<AltitudeType>> exportHigraHierarchy() const {
        return tree_altitude_ops::exportHigraHierarchy(tree, altitude);
    }

    WeightedTreeEditor edit();
};

/**
 * @brief Edit-session facade for `WeightedMorphologicalTree`.
 *
 * The weighted editor reuses the structural `TreeEditor` for topology while
 * keeping the external altitude buffer as the canonical weighted state.
 */
class WeightedTreeEditor {
private:
    WeightedMorphologicalTree& weighted_;
    TreeEditor editor_;

public:
    explicit WeightedTreeEditor(WeightedMorphologicalTree& weighted)
        : weighted_(weighted), editor_(weighted.tree) {}

    NodeId createDetachedNode(AltitudeType altitude = AltitudeType{}) {
        const NodeId nodeId = editor_.createDetachedNode();
        weighted_.altitude.resize(static_cast<size_t>(weighted_.tree.getNumInternalNodeSlots()), AltitudeType{});
        weighted_.altitude[static_cast<size_t>(nodeId)] = altitude;
        return nodeId;
    }

    void setNodeAltitude(NodeId nodeId, AltitudeType altitude) {
        if (!weighted_.tree.isAlive(nodeId)) {
            throw std::invalid_argument("WeightedTreeEditor::setNodeAltitude requires a live node.");
        }
        weighted_.altitude[static_cast<size_t>(nodeId)] = altitude;
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
};

inline WeightedTreeEditor WeightedMorphologicalTree::edit() {
    return WeightedTreeEditor(*this);
}

} // namespace mmcfilters
