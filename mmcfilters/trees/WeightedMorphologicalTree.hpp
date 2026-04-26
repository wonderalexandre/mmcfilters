#pragma once

#include "MorphologicalTree.hpp"
#include "TreeAltitudeOps.hpp"
#include "TreeEditor.hpp"

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

    void validateAltitudeBufferShape() const {
        tree_altitude_ops::validateAltitudeBufferShape(tree_, altitude_);
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
        tree_altitude_ops::validateMonotoneAltitude(tree_, altitude_);
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
        return tree_altitude_ops::getNodeResidue(tree_, altitude_, nodeId);
    }

    ImageUInt8Ptr reconstructionImage() const {
        return tree_altitude_ops::reconstructImage(tree_, altitude_);
    }

    /**
     * @brief Exports the current live rooted tree to a new compact Higra parent/altitude representation.
     */
    std::pair<std::vector<NodeId>, std::vector<AltitudeType>> exportHigraHierarchy() const {
        return tree_altitude_ops::exportHigraHierarchy(tree_, altitude_);
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
