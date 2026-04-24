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
 * @brief Transitional Higra-style wrapper pairing topology with an external altitude buffer.
 *
 * `WeightedMorphologicalTree` keeps the mutable topology in `tree` and the
 * node-altitude state in `altitude`, indexed by dense internal `NodeId`. The
 * wrapper offers the weighted conveniences that are intentionally kept out of
 * `MorphologicalTree`.
 */
class WeightedMorphologicalTree {
    friend class WeightedTreeEditor;

    void configureEmptyTopology(int rows, int cols, int treeType, std::optional<AdjacencyRelation> adjacency, NodeId numProperParts) {
        tree.treeType_ = treeType;
        tree.numRows_ = rows;
        tree.numCols_ = cols;
        tree.adj_ = std::move(adjacency);
        tree.resetEmptyStorage(static_cast<size_t>(numProperParts));
        altitude.clear();
    }

    void configureImportedTopology(int rows, int cols, bool isMaxtree, std::optional<AdjacencyRelation> adjacency) {
        tree.treeType_ = isMaxtree ? MorphologicalTree::MAX_TREE : MorphologicalTree::MIN_TREE;
        tree.numRows_ = rows;
        tree.numCols_ = cols;
        tree.adj_ = std::move(adjacency);
        tree.numNodes_ = 0;
    }

    void assignAltitudeFromDirectProperParts(const ImageUInt8Ptr& img) {
        altitude.assign(static_cast<size_t>(tree.getNumInternalNodeSlots()), AltitudeType{});
        for (NodeId nodeId : tree.getAliveNodeIds()) {
            const auto& properParts = tree.getProperParts(nodeId);
            if (properParts.empty()) {
                throw std::runtime_error("Cannot infer node altitude from a topology node without direct proper parts.");
            }
            altitude[static_cast<size_t>(nodeId)] = static_cast<AltitudeType>((*img)[properParts.front()]);
        }
    }

    template <typename AltitudeValue>
        requires std::is_arithmetic_v<AltitudeValue>
    void importAltitudeFromHigra(std::span<const AltitudeValue> higraAltitude) {
        if (static_cast<size_t>(tree.getNumHigraNodes()) != higraAltitude.size()) {
            throw std::invalid_argument("Higra altitude buffer size must match the imported hierarchy.");
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

    WeightedMorphologicalTree() = default;

    WeightedMorphologicalTree(MorphologicalTree&& topology, AltitudeBuffer altitudeBuffer)
        : tree(std::move(topology)), altitude(std::move(altitudeBuffer)) {
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

    WeightedMorphologicalTree(
        std::span<const NodeId> parent,
        NodeId numProperParts,
        int rows,
        int cols,
        bool isMaxtree,
        std::optional<AdjacencyRelation> adjacency = std::nullopt)
        : tree(parent, numProperParts, rows, cols, isMaxtree, std::move(adjacency)),
          altitude(static_cast<size_t>(tree.getNumInternalNodeSlots()), AltitudeType{}) {}

    WeightedMorphologicalTree(
        std::span<const NodeId> parent,
        int rows,
        int cols,
        bool isMaxtree,
        std::optional<AdjacencyRelation> adjacency = std::nullopt)
        : WeightedMorphologicalTree(parent, rows * cols, rows, cols, isMaxtree, std::move(adjacency)) {}

    template <typename AltitudeValue>
        requires std::is_arithmetic_v<AltitudeValue>
    WeightedMorphologicalTree(
        std::span<const NodeId> higraParent,
        std::span<const AltitudeValue> higraAltitude,
        NodeId numProperParts,
        int rows,
        int cols,
        bool isMaxtree,
        std::optional<AdjacencyRelation> adjacency = std::nullopt) {
        configureImportedTopology(rows, cols, isMaxtree, std::move(adjacency));
        tree.resetFromHigraTopology(higraParent, numProperParts);
        importAltitudeFromHigra(higraAltitude);
    }

    template <typename AltitudeValue>
        requires std::is_arithmetic_v<AltitudeValue>
    WeightedMorphologicalTree(
        std::span<const NodeId> higraParent,
        std::span<const AltitudeValue> higraAltitude,
        int rows,
        int cols,
        bool isMaxtree,
        std::optional<AdjacencyRelation> adjacency = std::nullopt)
        : WeightedMorphologicalTree(higraParent, higraAltitude, rows * cols, rows, cols, isMaxtree, std::move(adjacency)) {}

    template <class ParentRange>
    WeightedMorphologicalTree(
        const ParentRange& parent,
        NodeId numProperParts,
        int rows,
        int cols,
        bool isMaxtree,
        std::optional<AdjacencyRelation> adjacency = std::nullopt)
        : WeightedMorphologicalTree(
            std::span<const NodeId>(parent.data(), parent.size()),
            numProperParts,
            rows,
            cols,
            isMaxtree,
            std::move(adjacency)) {}

    template <class ParentRange>
    WeightedMorphologicalTree(
        const ParentRange& parent,
        int rows,
        int cols,
        bool isMaxtree,
        std::optional<AdjacencyRelation> adjacency = std::nullopt)
        : WeightedMorphologicalTree(
            std::span<const NodeId>(parent.data(), parent.size()),
            rows,
            cols,
            isMaxtree,
            std::move(adjacency)) {}

    template <class ParentRange, class AltitudeRange>
        requires std::is_arithmetic_v<typename AltitudeRange::value_type>
    WeightedMorphologicalTree(
        const ParentRange& higraParent,
        const AltitudeRange& higraAltitude,
        NodeId numProperParts,
        int rows,
        int cols,
        bool isMaxtree,
        std::optional<AdjacencyRelation> adjacency = std::nullopt)
        : WeightedMorphologicalTree(
            std::span<const NodeId>(higraParent.data(), higraParent.size()),
            std::span<const typename AltitudeRange::value_type>(higraAltitude.data(), higraAltitude.size()),
            numProperParts,
            rows,
            cols,
            isMaxtree,
            std::move(adjacency)) {}

    template <class ParentRange, class AltitudeRange>
        requires std::is_arithmetic_v<typename AltitudeRange::value_type>
    WeightedMorphologicalTree(
        const ParentRange& higraParent,
        const AltitudeRange& higraAltitude,
        int rows,
        int cols,
        bool isMaxtree,
        std::optional<AdjacencyRelation> adjacency = std::nullopt)
        : WeightedMorphologicalTree(
            std::span<const NodeId>(higraParent.data(), higraParent.size()),
            std::span<const typename AltitudeRange::value_type>(higraAltitude.data(), higraAltitude.size()),
            rows,
            cols,
            isMaxtree,
            std::move(adjacency)) {}

    static WeightedMorphologicalTree create(int rows, int cols, bool isMaxtree, std::optional<AdjacencyRelation> adjacency = std::nullopt) {
        WeightedMorphologicalTree weighted;
        weighted.configureEmptyTopology(
            rows,
            cols,
            isMaxtree ? MorphologicalTree::MAX_TREE : MorphologicalTree::MIN_TREE,
            adjacency ? std::move(adjacency) : std::optional<AdjacencyRelation>(std::in_place, rows, cols, 1.5),
            static_cast<NodeId>(rows * cols));
        return weighted;
    }

    template <typename AltitudeValue>
        requires std::is_arithmetic_v<AltitudeValue>
    static WeightedMorphologicalTree createFromHigra(
        std::span<const NodeId> higraParent,
        std::span<const AltitudeValue> higraAltitude,
        NodeId numProperParts,
        int rows,
        int cols,
        bool isMaxtree,
        std::optional<AdjacencyRelation> adjacency = std::nullopt) {
        return WeightedMorphologicalTree(
            higraParent,
            higraAltitude,
            numProperParts,
            rows,
            cols,
            isMaxtree,
            std::move(adjacency));
    }

    template <class ParentRange, class AltitudeRange>
        requires std::is_arithmetic_v<typename AltitudeRange::value_type>
    static WeightedMorphologicalTree createFromHigra(
        const ParentRange& higraParent,
        const AltitudeRange& higraAltitude,
        NodeId numProperParts,
        int rows,
        int cols,
        bool isMaxtree,
        std::optional<AdjacencyRelation> adjacency = std::nullopt) {
        return createFromHigra(
            std::span<const NodeId>(higraParent.data(), higraParent.size()),
            std::span<const typename AltitudeRange::value_type>(higraAltitude.data(), higraAltitude.size()),
            numProperParts,
            rows,
            cols,
            isMaxtree,
            std::move(adjacency));
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

    std::pair<std::vector<NodeId>, std::vector<AltitudeType>> exportHigraHierarchy() const {
        return tree_altitude_ops::exportHigraHierarchy(tree, altitude);
    }

    template <typename AltitudeValue>
        requires std::is_arithmetic_v<AltitudeValue>
    void resetFromHigra(std::span<const NodeId> parent, std::span<const AltitudeValue> higraAltitude, NodeId numProperParts) {
        tree.resetFromHigraTopology(parent, numProperParts);
        importAltitudeFromHigra(higraAltitude);
    }

    template <typename AltitudeValue>
        requires std::is_arithmetic_v<AltitudeValue>
    void resetFromHigra(std::span<const NodeId> parent, std::span<const AltitudeValue> higraAltitude) {
        if (tree.getNumTotalProperParts() <= 0) {
            throw std::logic_error("resetFromHigra(parent, altitude) requires a known proper-part domain.");
        }
        resetFromHigra(parent, higraAltitude, static_cast<NodeId>(tree.getNumTotalProperParts()));
    }

    template <class ParentRange, class AltitudeRange>
        requires std::is_arithmetic_v<typename AltitudeRange::value_type>
    void resetFromHigra(const ParentRange& parent, const AltitudeRange& higraAltitude, NodeId numProperParts) {
        resetFromHigra(
            std::span<const NodeId>(parent.data(), parent.size()),
            std::span<const typename AltitudeRange::value_type>(higraAltitude.data(), higraAltitude.size()),
            numProperParts);
    }

    template <class ParentRange, class AltitudeRange>
        requires std::is_arithmetic_v<typename AltitudeRange::value_type>
    void resetFromHigra(const ParentRange& parent, const AltitudeRange& higraAltitude) {
        resetFromHigra(
            std::span<const NodeId>(parent.data(), parent.size()),
            std::span<const typename AltitudeRange::value_type>(higraAltitude.data(), higraAltitude.size()));
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
