#pragma once

#include "MorphologicalTree.hpp"
#include "TreeAltitudeAlgorithms.hpp"
#include "TreeEditor.hpp"
#include "WeightedTreeView.hpp"
#include "../utils/Image.hpp"

#include <memory>
#include <span>
#include <stdexcept>
#include <optional>
#include <utility>
#include <vector>

namespace mmcfilters {

template <AltitudeValue T> class WeightedTreeEditor;

template <AltitudeValue T> class WeightedMorphologicalTree;

namespace detail {

template <AltitudeValue T> [[nodiscard]] WeightedTreeEditor<T> beginEstablishedWeightedEdit(WeightedMorphologicalTree<T>& tree);

} // namespace detail

/**
 * @brief Wrapper pairing `MorphologicalTree` topology with an external altitude buffer.
 *
 * `WeightedMorphologicalTree` owns the mutable topology and node-altitude state
 * internally, indexed by dense internal `NodeId`. Read-only topology access and
 * weighted convenience methods are exposed without leaking structural mutation.
 * Staged structural edits must go through `edit()`, which returns a
 * `WeightedTreeEditor`.
 *
 * Static Higra imports preserve their original node-id domain until the
 * topology is edited. `exportHigraHierarchy()` creates a new compact Higra
 * domain for the current live rooted tree.
 */
template <AltitudeValue T> class WeightedMorphologicalTree {
    friend class WeightedTreeEditor<T>;

  public:
    /// Altitude scalar type stored by this weighted tree.
    using altitude_type = T;

    /// Dense altitude buffer indexed by internal `NodeId`.
    using altitude_buffer = std::vector<T>;

  private:
    /** @brief References the tree used by the component. */
    MorphologicalTree tree_;
    /** @brief Stores the altitude. */
    altitude_buffer altitude_;

    /**
     * @brief Copies an altitude buffer already indexed by internal node slots.
     *
     * @param altitudeValues Altitude or level represented by `altitudeValues`.
     */
    void assignInternalAltitude(std::span<const T> altitudeValues) {
        if (altitudeValues.size() != static_cast<size_t>(tree_.getNumInternalNodeSlots())) {
            throw std::invalid_argument("Internal altitude buffer size must match the dense internal-node domain.");
        }
        TreeAltitudeAlgorithms::validateFiniteAltitudeValues(altitudeValues, "WeightedMorphologicalTree internal altitude input");
        altitude_.assign(altitudeValues.begin(), altitudeValues.end());
    }

    /**
     * @brief Takes ownership of an altitude buffer already indexed by internal
     * node slots.
     *
     * @param altitudeValues Altitude or level represented by `altitudeValues`.
     */
    void assignInternalAltitude(altitude_buffer&& altitudeValues) {
        if (altitudeValues.size() != static_cast<size_t>(tree_.getNumInternalNodeSlots())) {
            throw std::invalid_argument("Internal altitude buffer size must match the dense internal-node domain.");
        }
        TreeAltitudeAlgorithms::validateFiniteAltitudeValues(std::span<const T>(altitudeValues), "WeightedMorphologicalTree internal altitude input");
        altitude_ = std::move(altitudeValues);
    }

    /**
     * @brief Returns whether this tree kind has no global monotone altitude order.
     *
     * @param tree Tree topology used by the operation.
     * @return Whether this tree kind has no global monotone altitude order.
     */
    static bool skipsMonotoneValidation(const MorphologicalTree& tree) noexcept { return tree.getAltitudeOrder() == AltitudeOrder::UNCONSTRAINED; }

    /**
     * @brief Checks the local parent/children altitude constraints for one edit.
     *
     * This keeps `setAltitude()` proportional to the node degree instead of
     * scanning the full tree.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @param value Value used by the operation.
     * @return True when the documented condition holds; otherwise false.
     */
    [[nodiscard]] bool validateLocalMonotoneAltitudeUpdate(NodeId nodeId, T value) const {
        if (skipsMonotoneValidation(tree_)) {
            return false;
        }

        const bool increasingFromRoot = tree_.getAltitudeOrder() == AltitudeOrder::INCREASING_FROM_ROOT;
        if (!tree_.isRoot(nodeId)) {
            const NodeId parentNodeId = tree_.getNodeParent(nodeId);
            if (parentNodeId == InvalidNode || !tree_.isAlive(parentNodeId)) {
                throw std::runtime_error("Monotone altitude update requires an alive non-root node to have an alive parent.");
            }
            const T parentAltitude = altitude_[static_cast<size_t>(parentNodeId)];
            if (increasingFromRoot && parentAltitude >= value) {
                throw std::runtime_error("Hierarchy altitude update must remain strictly increasing from parent to child.");
            }
            if (!increasingFromRoot && parentAltitude <= value) {
                throw std::runtime_error("Hierarchy altitude update must remain strictly decreasing from parent to child.");
            }
        }

        for (NodeId childId : tree_.getChildren(nodeId)) {
            const T childAltitude = altitude_[static_cast<size_t>(childId)];
            if (increasingFromRoot && value >= childAltitude) {
                throw std::runtime_error("Hierarchy altitude update must remain strictly increasing from parent to child.");
            }
            if (!increasingFromRoot && value <= childAltitude) {
                throw std::runtime_error("Hierarchy altitude update must remain strictly decreasing from parent to child.");
            }
        }
        return true;
    }

  private:
    /**
     * @brief Disables default construction.
     */
    WeightedMorphologicalTree() = delete;

  public:
    /**
     * @brief Disables copy construction.
     */
    WeightedMorphologicalTree(const WeightedMorphologicalTree&) = delete;
    /**
     * @brief Disables copy assignment.
     */
    WeightedMorphologicalTree& operator=(const WeightedMorphologicalTree&) = delete;

    /**
     * @brief Transfers committed topology and altitude ownership.
     *
     * The topology move rejects an active weighted editor before any state is
     * transferred.
     *
     * @param other Object to compare with or transfer from.
     */
    WeightedMorphologicalTree(WeightedMorphologicalTree&& other) : tree_(std::move(other.tree_)), altitude_(std::move(other.altitude_)) {}

    /**
     * @brief Move-assigns committed topology and altitude ownership.
     *
     * Both topology owners must be outside edit sessions. This keeps every
     * active `WeightedTreeEditor` bound to the storage it opened.
     *
     * @param other Object to compare with or transfer from.
     * @return Reference to the resulting object.
     */
    WeightedMorphologicalTree& operator=(WeightedMorphologicalTree&& other) {
        if (this != &other) {
            tree_ = std::move(other.tree_);
            altitude_ = std::move(other.altitude_);
        } else {
            tree_.requireNotEditing("WeightedMorphologicalTree self move assignment");
        }
        return *this;
    }

    /**
     * @brief Consumes a producer-owned altitude buffer in the internal node-id
     * domain.
     *
     * @param topology Tree topology used by the operation.
     * @param altitude Altitude data indexed by node identifier.
     */
    WeightedMorphologicalTree(detail::MorphologicalTreeConstructionTag, MorphologicalTree&& topology, altitude_buffer&& altitude) : tree_(std::move(topology)) {
        assignInternalAltitude(std::move(altitude));
        validateAltitudeBufferShape();
        validateMonotoneAltitude();
    }

    /**
     * @brief Checks that the altitude buffer covers the dense internal-node domain.
     */
    void validateAltitudeBufferShape() const { TreeAltitudeAlgorithms::validateAltitudeBufferShape(tree_, altitudeSpan()); }

    /**
     * @brief Returns read-only access to the owned topology.
     *
     * This accessor intentionally does not provide a mutable `MorphologicalTree`
     * handle. Use `edit()` for staged topology edits or the safe public
     * mutators for complete local changes.
     *
     * @return Read-only access to the owned topology.
     */
    [[nodiscard]] const MorphologicalTree& topology() const noexcept { return tree_; }

    /**
     * @brief Returns the dense altitude buffer indexed by internal `NodeId`.
     *
     * @return The dense altitude buffer indexed by internal NodeId.
     */
    [[nodiscard]] const altitude_buffer& getAltitudeBuffer() const noexcept { return altitude_; }

    /**
     * @brief Returns a read-only span over the dense altitude buffer.
     *
     * @return A read-only span over the dense altitude buffer.
     */
    [[nodiscard]] AltitudeSpan<T> altitudeSpan() const noexcept { return std::span<const T>(altitude_); }

    /**
     * @brief Creates a non-owning weighted view over this owner.
     *
     * @return The created non-owning weighted view over this owner.
     */
    [[nodiscard]] WeightedTreeView<T> asView() const { return WeightedTreeView<T>(tree_, altitudeSpan()); }

    /**
     * @brief Replaces the owned altitude buffer after full validation.
     *
     * Shape and finite-value checks are followed by validation of the
     * hierarchy's declared global altitude order.
     *
     * @param altitudeBuffer Altitude or level represented by `altitudeBuffer`.
     */
    void setAltitudeBuffer(altitude_buffer altitudeBuffer) {
        tree_.requireNotEditing("WeightedMorphologicalTree::setAltitudeBuffer");
        if (altitudeBuffer.size() != static_cast<size_t>(tree_.getNumInternalNodeSlots())) {
            throw std::runtime_error("Altitude buffer size must match the dense internal-node domain.");
        }
        TreeAltitudeAlgorithms::validateFiniteAltitudeValues(std::span<const T>(altitudeBuffer), "WeightedMorphologicalTree::setAltitudeBuffer");
        TreeAltitudeAlgorithms::validateMonotoneAltitude(tree_, std::span<const T>(altitudeBuffer));
        altitude_ = std::move(altitudeBuffer);
    }

    /**
     * @brief Validates the current altitude buffer against the topology order.
     */
    void validateMonotoneAltitude() const { TreeAltitudeAlgorithms::validateMonotoneAltitude(tree_, altitudeSpan()); }

    /**
     * @brief Returns one live node altitude from the dense buffer.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @return One live node altitude from the dense buffer.
     */
    [[nodiscard]] T getAltitude(NodeId nodeId) const {
        if (!tree_.isAlive(nodeId) || static_cast<size_t>(nodeId) >= altitude_.size()) {
            throw std::invalid_argument("WeightedMorphologicalTree::getAltitude requires a live internal NodeId.");
        }
        return altitude_[static_cast<size_t>(nodeId)];
    }

    /**
     * @brief Updates one live node altitude with local monotonicity validation.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @param value Value used by the operation.
     */
    void setAltitude(NodeId nodeId, T value) {
        tree_.requireNotEditing("WeightedMorphologicalTree::setAltitude");
        if (!tree_.isAlive(nodeId) || static_cast<size_t>(nodeId) >= altitude_.size()) {
            throw std::invalid_argument("WeightedMorphologicalTree::setAltitude requires a live internal NodeId.");
        }
        TreeAltitudeAlgorithms::validateFiniteAltitudeValue(value, static_cast<std::size_t>(nodeId), "WeightedMorphologicalTree::setAltitude");
        static_cast<void>(validateLocalMonotoneAltitudeUpdate(nodeId, value));
        altitude_[static_cast<size_t>(nodeId)] = value;
    }

    /**
     * @brief Prunes a complete subtree through the owned topology.
     *
     * The altitude buffer is kept as the canonical dense weighted state; dead
     * slots may retain their old values until a compact export is requested.
     *
     * @param nodeId Identifier of the node used by the operation.
     */
    void pruneNode(NodeId nodeId) {
        tree_.requireNotEditing("WeightedMorphologicalTree::pruneNode");
        tree_.pruneNode(nodeId);
    }

    /**
     * @brief Merges one node into its parent through the owned topology.
     *
     * This is the weighted counterpart of the safe public topology mutator. Dead
     * slots keep stale altitude values until a compact export is requested.
     *
     * @param nodeId Identifier of the node used by the operation.
     */
    void mergeNodeIntoParent(NodeId nodeId) {
        tree_.requireNotEditing("WeightedMorphologicalTree::mergeNodeIntoParent");
        tree_.mergeNodeIntoParent(nodeId);
    }

    /**
     * @brief Returns the altitude difference between a node and its parent.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @return The altitude difference between a node and its parent.
     */
    [[nodiscard]] AltitudeDiff<T> getNodeResidue(NodeId nodeId) const { return TreeAltitudeAlgorithms::getNodeResidue(tree_, altitudeSpan(), nodeId); }

    /**
     * @brief Reconstructs an image by assigning each proper part its owner altitude.
     *
     * @return The reconstructed image by assigning each proper part its owner altitude.
     */
    [[nodiscard]] ImagePtr<T> reconstructionImage() const {
        return TreeAltitudeAlgorithms::reconstructImage(tree_, altitudeSpan(), "WeightedMorphologicalTree::reconstructImage");
    }

    /**
     * @brief Exports the current live rooted tree to a new compact Higra parent/altitude representation.
     *
     * @details
     * Attribute buffers are projected through
     * `AttributeComputation::projectNodeValuesToExportedHigra()` so
     * weighted-tree export remains limited to topology and altitudes.
     *
     * @return The exported current live rooted tree to a new compact Higra parent/altitude representation.
     */
    [[nodiscard]] std::pair<std::vector<NodeId>, std::vector<T>> exportHigraHierarchy() const {
        return TreeAltitudeAlgorithms::exportHigraHierarchy(tree_, altitudeSpan());
    }

    /**
     * @brief Opens the only public entrypoint for staged weighted edits.
     *
     * @return The opened only public entrypoint for staged weighted edits.
     */
    [[nodiscard]] WeightedTreeEditor<T> edit() { return WeightedTreeEditor<T>(*this); }
};

/**
 * @brief Edit-session facade for `WeightedMorphologicalTree`.
 *
 * The weighted editor reuses the structural `TreeEditor` for topology while
 * keeping the external altitude buffer as the canonical weighted state.
 * `commit()` first validates the topology and then validates monotone altitude.
 */
template <AltitudeValue T> class WeightedTreeEditor {
    friend class WeightedMorphologicalTree<T>;
    /** @brief Grants `detail::beginEstablishedWeightedEdit` access to the enclosing type. */
    friend WeightedTreeEditor<T> detail::beginEstablishedWeightedEdit<T>(WeightedMorphologicalTree<T>& tree);

  private:
    /** @brief Records altitude changes for rollback of a weighted-tree edit. */
    struct AltitudeRollbackJournal {
        /** @brief Stores the original size. */
        std::size_t originalSize = 0;
        /** @brief Stores the captured. */
        TreeEditor::DeltaNodeSet captured;
        /** @brief Stores the values. */
        std::vector<std::pair<NodeId, T>> values;
    };

    /** @brief Stores the weighted. */
    WeightedMorphologicalTree<T>& weighted_;
    /** @brief Stores the editor. */
    TreeEditor editor_;
    /** @brief Stores the original altitude size. */
    std::size_t originalAltitudeSize_ = 0;
    /** @brief Stores the altitude rollback journal. */
    std::unique_ptr<AltitudeRollbackJournal> altitudeRollbackJournal_;
    /** @brief Stores the edit revision. */
    std::size_t editRevision_ = 0;
    /** @brief Stores the proven revision. */
    std::optional<std::size_t> provenRevision_;

    /**
     * @brief Opens the underlying topology edit session.
     *
     * @param weighted Weighted tree used by the operation.
     * @param transactional Whether the edit records rollback information.
     * @param invariantsEstablishedByConstruction Whether construction already established every edit invariant.
     */
    explicit WeightedTreeEditor(WeightedMorphologicalTree<T>& weighted, [[maybe_unused]] bool transactional = false,
                                bool invariantsEstablishedByConstruction = false)
        : weighted_(weighted), editor_(weighted.tree_, invariantsEstablishedByConstruction), originalAltitudeSize_(weighted.altitude_.size()) {}

    /**
     * @brief Ensures altitude rollback journal.
     */
    void ensureAltitudeRollbackJournal() {
        if (!editor_.canRollback()) {
            return;
        }
        editor_.ensureRollbackJournal();
        if (!altitudeRollbackJournal_) {
            altitudeRollbackJournal_ = std::make_unique<AltitudeRollbackJournal>();
            altitudeRollbackJournal_->originalSize = originalAltitudeSize_;
        }
    }

    /**
     * @brief Captures altitude for rollback.
     *
     * @param nodeId Identifier of the node used by the operation.
     */
    void captureAltitudeForRollback(NodeId nodeId) {
        ensureAltitudeRollbackJournal();
        if (!altitudeRollbackJournal_) {
            return;
        }
        if (nodeId < 0 || static_cast<std::size_t>(nodeId) >= originalAltitudeSize_) {
            return;
        }
        if (altitudeRollbackJournal_->captured.contains(nodeId)) {
            return;
        }
        altitudeRollbackJournal_->values.reserve(altitudeRollbackJournal_->values.size() + 1);
        if (!altitudeRollbackJournal_->captured.insert(nodeId)) {
            return;
        }
        altitudeRollbackJournal_->values.emplace_back(nodeId, weighted_.altitude_[static_cast<std::size_t>(nodeId)]);
    }

    /**
     * @brief Restores altitude journal.
     */
    void restoreAltitudeJournal() noexcept {
        if (!altitudeRollbackJournal_ || !editor_.canRollback()) {
            return;
        }
        weighted_.altitude_.resize(altitudeRollbackJournal_->originalSize);
        for (const auto& [nodeId, altitude] : altitudeRollbackJournal_->values) {
            weighted_.altitude_[static_cast<std::size_t>(nodeId)] = altitude;
        }
    }

    /**
     * @brief Records mutation.
     */
    void recordMutation() noexcept {
        if (editor_.invariantsEstablishedByConstruction_) {
            return;
        }
        ++editRevision_;
        provenRevision_.reset();
    }

  public:
    /**
     * @brief Disables copy construction.
     */
    WeightedTreeEditor(const WeightedTreeEditor&) = delete;
    /**
     * @brief Disables copy assignment.
     */
    WeightedTreeEditor& operator=(const WeightedTreeEditor&) = delete;

    /**
     * @brief Transfers the active weighted edit session and rollback journal.
     *
     * @param other Object to compare with or transfer from.
     */
    WeightedTreeEditor(WeightedTreeEditor&& other) noexcept
        : weighted_(other.weighted_), editor_(std::move(other.editor_)), originalAltitudeSize_(other.originalAltitudeSize_),
          altitudeRollbackJournal_(std::move(other.altitudeRollbackJournal_)), editRevision_(other.editRevision_), provenRevision_(other.provenRevision_) {
        other.altitudeRollbackJournal_.reset();
        other.provenRevision_.reset();
    }

    /**
     * @brief Disables move assignment.
     */
    WeightedTreeEditor& operator=(WeightedTreeEditor&&) = delete;

    /**
     * @brief Restores the altitude journal and closes the weighted edit session.
     */
    ~WeightedTreeEditor() { restoreAltitudeJournal(); }

    /**
     * @brief Tests whether topology and altitude can be rolled back.
     *
     * @return True if topology and altitude can be rolled back; otherwise false.
     */
    [[nodiscard]] bool canRollback() const noexcept { return editor_.canRollback(); }

    /**
     * @brief Restores topology and altitude captured by the delta journals.
     */
    void rollback() {
        if (!canRollback()) {
            throw std::logic_error("WeightedTreeEditor::rollback is unavailable for the internal journal-free editor.");
        }
        restoreAltitudeJournal();
        editor_.rollback();
    }

    /**
     * @brief Creates a detached topology node and initializes its altitude.
     *
     * The altitude buffer is reserved before mutating the topology so allocation
     * failure cannot leave a new node without a matching altitude slot.
     *
     * @param altitude Altitude data indexed by node identifier.
     * @return The created detached topology node and initializes its altitude.
     */
    [[nodiscard("Discarding a detached node id makes the staged edit impossible to complete safely")]] NodeId createDetachedNode(T altitude = T{}) {
        TreeAltitudeAlgorithms::validateFiniteAltitudeValue(altitude, 0, "WeightedTreeEditor::createDetachedNode");
        const std::size_t requiredAltitudeSize =
            static_cast<std::size_t>(weighted_.tree_.getNumInternalNodeSlots()) + (weighted_.tree_.getNumFreeNodeSlots() == 0 ? 1u : 0u);
        weighted_.altitude_.reserve(requiredAltitudeSize);
        ensureAltitudeRollbackJournal();
        const NodeId nodeId = editor_.createDetachedNode();
        weighted_.altitude_.resize(static_cast<size_t>(weighted_.tree_.getNumInternalNodeSlots()), T{});
        captureAltitudeForRollback(nodeId);
        weighted_.altitude_[static_cast<size_t>(nodeId)] = altitude;
        recordMutation();
        return nodeId;
    }

    /**
     * @brief Sets a live node altitude during a staged topology edit.
     *
     * Monotone order is intentionally checked at weighted commit time because
     * intermediate staged topologies may not yet have final parent/child
     * relations.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @param altitude Altitude data indexed by node identifier.
     */
    void setNodeAltitude(NodeId nodeId, T altitude) {
        if (!weighted_.tree_.isAlive(nodeId)) {
            throw std::invalid_argument("WeightedTreeEditor::setNodeAltitude requires a live node.");
        }
        TreeAltitudeAlgorithms::validateFiniteAltitudeValue(altitude, static_cast<std::size_t>(nodeId), "WeightedTreeEditor::setNodeAltitude");
        captureAltitudeForRollback(nodeId);
        weighted_.altitude_[static_cast<size_t>(nodeId)] = altitude;
        editor_.touch(nodeId);
        recordMutation();
    }

    /**
     * @brief Detaches one non-root node through the structural editor.
     *
     * @param nodeId Identifier of the node used by the operation.
     */
    void detach(NodeId nodeId) {
        editor_.detach(nodeId);
        recordMutation();
    }

    /**
     * @brief Reparents one node through the structural editor.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @param newParentId Parent-node value represented by `newParentId`.
     */
    void reparent(NodeId nodeId, NodeId newParentId) {
        editor_.reparent(nodeId, newParentId);
        recordMutation();
    }

    /**
     * @brief Attaches one detached node through the structural editor.
     *
     * @param parentId Identifier of the parent node.
     * @param detachedNodeId Node identifier represented by `detachedNodeId`.
     */
    void attach(NodeId parentId, NodeId detachedNodeId) {
        editor_.attach(parentId, detachedNodeId);
        recordMutation();
    }

    /**
     * @brief Moves all direct children from `sourceId` under `parentId`.
     *
     * @param parentId Identifier of the parent node.
     * @param sourceId Input represented by `sourceId`.
     */
    void moveChildren(NodeId parentId, NodeId sourceId) {
        editor_.moveChildren(parentId, sourceId);
        recordMutation();
    }

    /**
     * @brief Moves one direct proper part between nodes.
     *
     * @param targetNodeId Node identifier represented by `targetNodeId`.
     * @param sourceNodeId Node identifier represented by `sourceNodeId`.
     * @param properPartId Proper-part identifier used by the operation.
     */
    void moveProperPart(NodeId targetNodeId, NodeId sourceNodeId, NodeId properPartId) {
        editor_.moveProperPart(targetNodeId, sourceNodeId, properPartId);
        recordMutation();
    }

    /**
     * @brief Moves all direct proper parts between nodes.
     *
     * @param targetNodeId Node identifier represented by `targetNodeId`.
     * @param sourceNodeId Node identifier represented by `sourceNodeId`.
     */
    void moveProperParts(NodeId targetNodeId, NodeId sourceNodeId) {
        editor_.moveProperParts(targetNodeId, sourceNodeId);
        recordMutation();
    }

    /**
     * @brief Detaches a direct child and optionally releases an empty node slot.
     *
     * @param parentNodeId Identifier of the parent node.
     * @param childId Identifier of the child node.
     * @param releaseNodeFlag Flag controlling release node flag.
     */
    void removeChild(NodeId parentNodeId, NodeId childId, bool releaseNodeFlag) {
        editor_.removeChild(parentNodeId, childId, releaseNodeFlag);
        recordMutation();
    }

    /**
     * @brief Releases an empty detached node slot.
     *
     * @param nodeId Identifier of the node used by the operation.
     */
    void releaseNode(NodeId nodeId) {
        editor_.releaseNode(nodeId);
        recordMutation();
    }

    /**
     * @brief Promotes one node to the topology root.
     *
     * @param nodeId Identifier of the node used by the operation.
     */
    void setRoot(NodeId nodeId) {
        editor_.setRoot(nodeId);
        recordMutation();
    }

    /**
     * @brief Applies the topology prune helper inside the weighted edit session.
     *
     * @param nodeId Identifier of the node used by the operation.
     */
    void pruneNode(NodeId nodeId) {
        editor_.pruneNode(nodeId);
        recordMutation();
    }

    /**
     * @brief Applies the topology merge helper inside the weighted edit session.
     *
     * @param nodeId Identifier of the node used by the operation.
     */
    void mergeNodeIntoParent(NodeId nodeId) {
        editor_.mergeNodeIntoParent(nodeId);
        recordMutation();
    }

    /**
     * @brief Returns whether the structural edit still has detached live nodes.
     *
     * @return Whether the structural edit still has detached live nodes.
     */
    [[nodiscard]] bool hasDetachedAliveNodes() const noexcept { return editor_.hasDetachedAliveNodes(); }

    /**
     * @brief Produces a generic move-only proof for the current weighted edit
     * revision.
     *
     * Topology is checked on the mutation delta whenever every primitive is
     * supported. Altitude order is then checked only around touched nodes. A
     * topology fallback also selects complete altitude validation.
     *
     * @return The produced generic move-only proof for the current weighted edit revision.
     */
    [[nodiscard]] TreeEditor::IncrementalProof proveIncremental() {
        bool strictAltitudeExcludesCycles = !WeightedMorphologicalTree<T>::skipsMonotoneValidation(weighted_.tree_);
        if (!editor_.invariantsEstablishedByConstruction_) {
            for (NodeId node : editor_.touchedNodes_.entries()) {
                if (!weighted_.tree_.isAlive(node)) {
                    continue;
                }
                strictAltitudeExcludesCycles =
                    weighted_.validateLocalMonotoneAltitudeUpdate(node, weighted_.altitude_[static_cast<std::size_t>(node)]) && strictAltitudeExcludesCycles;
            }
        }

        auto proof = strictAltitudeExcludesCycles ? editor_.proveIncrementalWithStrictAltitudeAcyclicity() : editor_.proveIncremental();
#ifndef NDEBUG
        // The assertion-enabled oracle runs exactly once, including when the
        // topology proof already fell back to complete validation.
        weighted_.validateMonotoneAltitude();
#else
        if (proof.usedCompleteValidation()) {
            weighted_.validateMonotoneAltitude();
        }
#endif
        provenRevision_ = editRevision_;
        return proof;
    }

    /**
     * @brief Commits the exact weighted edit revision represented by `proof`.
     *
     * @param proof Validation proof consumed by the operation.
     */
    void commit(TreeEditor::IncrementalProof&& proof) {
        if (!provenRevision_ || *provenRevision_ != editRevision_) {
            throw std::logic_error("Incremental weighted proof is stale or belongs to another edit revision.");
        }
        editor_.commit(std::move(proof));
        provenRevision_.reset();
        altitudeRollbackJournal_.reset();
    }

    /**
     * @brief Validates only the staged topology.
     *
     * Altitude monotonicity is checked by `validateAndCommit()` / `commit()`.
     *
     * @return Validation result for only the staged topology.
     */
    [[nodiscard]] TreeValidationResult validate() const noexcept { return editor_.validate(); }

    /**
     * @brief Validates topology, validates altitude order, then closes the edit.
     *
     * @return Validation result for topology, validates altitude order, then closes the edit.
     */
    [[nodiscard("Inspect the validation result or use commit() for exception-based failure handling")]] TreeValidationResult validateAndCommit() noexcept {
        TreeValidationResult result = editor_.validate();
        if (!result.ok) {
            return result;
        }
        try {
            weighted_.validateMonotoneAltitude();
        } catch (const std::exception& ex) {
            return {false, ex.what()};
        } catch (...) {
            return {false, "WeightedTreeEditor monotone-altitude validation failed with an unknown error."};
        }
        editor_.finishCommit(TreeEditValidationMode::Complete);
        altitudeRollbackJournal_.reset();
        provenRevision_.reset();
        return result;
    }

    /**
     * @brief Exception-based wrapper around `validateAndCommit()`.
     */
    void commit() {
        TreeValidationResult result = validateAndCommit();
        if (!result.ok) {
            throw std::runtime_error(result.message);
        }
    }
};

namespace detail {

/**
 * @brief Internal entrypoint for algorithms that establish every edit
 * invariant during their existing local update.
 *
 * Release builds avoid duplicating that work. Assertion-enabled builds still
 * execute the complete topology and altitude oracle when the proof is issued.
 *
 * @param tree Tree topology used by the operation.
 * @return Editor for an invariant-preserving weighted edit.
 */
template <AltitudeValue T> [[nodiscard]] WeightedTreeEditor<T> beginEstablishedWeightedEdit(WeightedMorphologicalTree<T>& tree) {
    return WeightedTreeEditor<T>(tree, false, true);
}

} // namespace detail

} // namespace mmcfilters
