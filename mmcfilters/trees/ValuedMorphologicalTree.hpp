#pragma once

#include "MorphologicalTree.hpp"
#include "TreeAltitudeAlgorithms.hpp"
#include "TreeEditor.hpp"
#include "ValuedMorphologicalTreeView.hpp"
#include "../utils/Contract.hpp"
#include "../utils/Image.hpp"

#include <memory>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <optional>
#include <utility>
#include <vector>

namespace mmcfilters {

template <AltitudeValue T> class ValuedMorphologicalTreeEditor;

template <AltitudeValue T> class ValuedMorphologicalTree;

namespace detail {

template <AltitudeValue T> [[nodiscard]] ValuedMorphologicalTreeEditor<T> beginEstablishedValuedEdit(ValuedMorphologicalTree<T>& tree);

} // namespace detail

/**
 * @brief Wrapper pairing `MorphologicalTree` topology with an external altitude buffer.
 *
 * `ValuedMorphologicalTree` owns the mutable topology and node-altitude state
 * internally, indexed by dense internal `NodeId`. Read-only topology access and
 * valued-tree convenience methods are exposed without leaking structural mutation.
 * Staged structural edits must go through `edit()`, which returns a
 * `ValuedMorphologicalTreeEditor`.
 *
 * Static Higra imports preserve their original node-id domain until the
 * topology is edited. `exportHigraHierarchy()` creates a new compact Higra
 * domain for the current live rooted tree.
 */
template <AltitudeValue T> class ValuedMorphologicalTree {
    friend class ValuedMorphologicalTreeEditor<T>;

  public:
    /// Altitude scalar type stored by this valued morphological tree.
    using AltitudeType = T;

    /// Dense altitude buffer indexed by internal `NodeId`.
    using NodeAltitudeBuffer = mmcfilters::NodeAltitudeBuffer<T>;

  private:
    /** @brief References the tree used by the component. */
    MorphologicalTree tree_;
    /** @brief Altitude. */
    NodeAltitudeBuffer nodeAltitudes_;

    /**
     * @brief Copies an altitude buffer already indexed by internal node slots.
     *
     * @param altitudeValues Altitude or level.
     */
    void assignInternalNodeAltitudes(std::span<const T> altitudeValues) {
        if (altitudeValues.size() != static_cast<size_t>(tree_.numInternalNodeSlots())) {
            throw std::invalid_argument("Internal altitude buffer size must match the dense internal-node domain.");
        }
        TreeAltitudeAlgorithms::validateFiniteAltitudeValues(altitudeValues, "ValuedMorphologicalTree internal altitude input");
        nodeAltitudes_.assign(altitudeValues.begin(), altitudeValues.end());
    }

    /**
     * @brief Takes ownership of an altitude buffer already indexed by internal
     * node slots.
     *
     * @param altitudeValues Altitude or level.
     */
    void assignInternalNodeAltitudes(NodeAltitudeBuffer&& altitudeValues) {
        if (altitudeValues.size() != static_cast<size_t>(tree_.numInternalNodeSlots())) {
            throw std::invalid_argument("Internal altitude buffer size must match the dense internal-node domain.");
        }
        TreeAltitudeAlgorithms::validateFiniteAltitudeValues(std::span<const T>(altitudeValues), "ValuedMorphologicalTree internal altitude input");
        nodeAltitudes_ = std::move(altitudeValues);
    }

    /**
     * @brief Returns whether this tree kind has no global monotone altitude order.
     *
     * @param tree Tree topology.
     * @return Whether this tree kind has no global monotone altitude order.
     */
    static bool skipsMonotoneValidation(const MorphologicalTree& tree) noexcept { return tree.nodeAltitudeOrder() == NodeAltitudeOrder::Unconstrained; }

    /**
     * @brief Checks the local parent/children altitude constraints for one edit.
     *
     * This keeps `setNodeAltitude()` proportional to the node degree instead of
     * scanning the full tree.
     *
     * @param nodeId Dense internal node identifier.
     * @param value Value.
     * @return True when the documented condition holds; otherwise false.
     */
    [[nodiscard]] bool validateLocalMonotoneNodeAltitudeUpdate(NodeId nodeId, T value) const {
        if (skipsMonotoneValidation(tree_)) {
            return false;
        }

        const bool increasingFromRoot = tree_.nodeAltitudeOrder() == NodeAltitudeOrder::Increasing;
        if (!tree_.isRoot(nodeId)) {
            const NodeId parentNodeId = tree_.parent(nodeId);
            if (parentNodeId == InvalidNode || !tree_.isAlive(parentNodeId)) {
                throw std::runtime_error("Monotone altitude update requires an alive non-root node to have an alive parent.");
            }
            const T parentAltitude = nodeAltitudes_[static_cast<size_t>(parentNodeId)];
            if (increasingFromRoot && parentAltitude >= value) {
                throw std::runtime_error("Hierarchy altitude update must remain strictly increasing from parent to child.");
            }
            if (!increasingFromRoot && parentAltitude <= value) {
                throw std::runtime_error("Hierarchy altitude update must remain strictly decreasing from parent to child.");
            }
        }

        for (NodeId childId : tree_.children(nodeId)) {
            const T childAltitude = nodeAltitudes_[static_cast<size_t>(childId)];
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
    ValuedMorphologicalTree() = delete;

  public:
    /**
     * @brief Disables copy construction.
     */
    ValuedMorphologicalTree(const ValuedMorphologicalTree&) = delete;
    /**
     * @brief Disables copy assignment.
     */
    ValuedMorphologicalTree& operator=(const ValuedMorphologicalTree&) = delete;

    /**
     * @brief Transfers committed topology and altitude ownership.
     *
     * The topology move rejects an active valued-tree editor before any state is
     * transferred.
     *
     * @param other Object to compare with or transfer from.
     */
    ValuedMorphologicalTree(ValuedMorphologicalTree&& other) : tree_(std::move(other.tree_)), nodeAltitudes_(std::move(other.nodeAltitudes_)) {}

    /**
     * @brief Move-assigns committed topology and altitude ownership.
     *
     * Both topology owners must be outside edit sessions. This keeps every
     * active `ValuedMorphologicalTreeEditor` bound to the storage it opened.
     *
     * @param other Object to compare with or transfer from.
     * @return Mutable reference to the updated object.
     */
    ValuedMorphologicalTree& operator=(ValuedMorphologicalTree&& other) {
        if (this != &other) {
            tree_ = std::move(other.tree_);
            nodeAltitudes_ = std::move(other.nodeAltitudes_);
        } else {
            tree_.requireNotEditing("ValuedMorphologicalTree self move assignment");
        }
        return *this;
    }

    /**
     * @brief Consumes a producer-owned altitude buffer in the internal node-id
     * domain.
     *
     * @param topology Tree topology.
     * @param altitude Altitude data indexed by node identifier.
     */
    ValuedMorphologicalTree(detail::MorphologicalTreeConstructionTag, MorphologicalTree&& topology, NodeAltitudeBuffer&& altitude) : tree_(std::move(topology)) {
        assignInternalNodeAltitudes(std::move(altitude));
        validateNodeAltitudeBufferShape();
        validateMonotoneNodeAltitudes();
    }

    /**
     * @brief Checks that the altitude buffer covers the dense internal-node domain.
     */
    void validateNodeAltitudeBufferShape() const { TreeAltitudeAlgorithms::validateNodeAltitudeBufferShape(tree_, nodeAltitudeSpan()); }

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
    [[nodiscard]] const NodeAltitudeBuffer& nodeAltitudes() const noexcept { return nodeAltitudes_; }

    /**
     * @brief Returns a read-only span over the dense altitude buffer.
     *
     * @return A read-only span over the dense altitude buffer.
     */
    [[nodiscard]] NodeAltitudeSpan<T> nodeAltitudeSpan() const noexcept { return std::span<const T>(nodeAltitudes_); }

    /**
     * @brief Creates a non-owning valued-tree view over this owner.
     *
     * @return The created non-owning valued-tree view over this owner.
     */
    [[nodiscard]] ValuedMorphologicalTreeView<T> asView() const { return ValuedMorphologicalTreeView<T>(tree_, nodeAltitudeSpan()); }

    /**
     * @brief Replaces the owned altitude buffer after full validation.
     *
     * Shape and finite-value checks are followed by validation of the
     * hierarchy's declared global altitude order.
     *
     * @param altitudeBuffer Altitude or level.
     */
    void setNodeAltitudes(NodeAltitudeBuffer altitudeBuffer) {
        tree_.requireNotEditing("ValuedMorphologicalTree::setNodeAltitudes");
        if (altitudeBuffer.size() != static_cast<size_t>(tree_.numInternalNodeSlots())) {
            throw std::runtime_error("Altitude buffer size must match the dense internal-node domain.");
        }
        TreeAltitudeAlgorithms::validateFiniteAltitudeValues(std::span<const T>(altitudeBuffer), "ValuedMorphologicalTree::setNodeAltitudes");
        TreeAltitudeAlgorithms::validateMonotoneNodeAltitudes(tree_, std::span<const T>(altitudeBuffer));
        nodeAltitudes_ = std::move(altitudeBuffer);
    }

    /**
     * @brief Validates the current altitude buffer against the topology order.
     */
    void validateMonotoneNodeAltitudes() const { TreeAltitudeAlgorithms::validateMonotoneNodeAltitudes(tree_, nodeAltitudeSpan()); }

    /**
     * @brief Returns one live node altitude from the dense buffer.
     *
     * @param nodeId Dense internal node identifier.
     * @return One live node altitude from the dense buffer.
     */
    [[nodiscard]] T nodeAltitude(NodeId nodeId) const {
        MMCFILTERS_CONTRACT_REQUIRE(tree_.isAlive(nodeId) && static_cast<size_t>(nodeId) < nodeAltitudes_.size(),
                                    throw std::invalid_argument("ValuedMorphologicalTree::nodeAltitude requires a live internal NodeId."));
        return nodeAltitudes_[static_cast<size_t>(nodeId)];
    }

    /**
     * @brief Updates one live node altitude with local monotonicity validation.
     *
     * @param nodeId Dense internal node identifier.
     * @param value Value.
     */
    void setNodeAltitude(NodeId nodeId, T value) {
        tree_.requireNotEditing("ValuedMorphologicalTree::setNodeAltitude");
        if (!tree_.isAlive(nodeId) || static_cast<size_t>(nodeId) >= nodeAltitudes_.size()) {
            throw std::invalid_argument("ValuedMorphologicalTree::setNodeAltitude requires a live internal NodeId.");
        }
        TreeAltitudeAlgorithms::validateFiniteAltitudeValue(value, static_cast<std::size_t>(nodeId), "ValuedMorphologicalTree::setNodeAltitude");
        static_cast<void>(validateLocalMonotoneNodeAltitudeUpdate(nodeId, value));
        nodeAltitudes_[static_cast<size_t>(nodeId)] = value;
    }

    /**
     * @brief Prunes a complete subtree through the owned topology.
     *
     * The altitude buffer is kept as the canonical dense valued-tree state; dead
     * slots may retain their old values until a compact export is requested.
     *
     * @param nodeId Dense internal node identifier.
     */
    void pruneNode(NodeId nodeId) {
        tree_.requireNotEditing("ValuedMorphologicalTree::pruneNode");
        tree_.pruneNode(nodeId);
    }

    /**
     * @brief Merges one node into its parent through the owned topology.
     *
     * This is the valued-tree counterpart of the safe public topology mutator. Dead
     * slots keep stale altitude values until a compact export is requested.
     *
     * @param nodeId Dense internal node identifier.
     */
    void mergeNodeIntoParent(NodeId nodeId) {
        tree_.requireNotEditing("ValuedMorphologicalTree::mergeNodeIntoParent");
        tree_.mergeNodeIntoParent(nodeId);
    }

    /**
     * @brief Returns the altitude difference between a node and its parent.
     *
     * The project uses a fixed zero reconstruction baseline, so the root
     * residue is equal to the root altitude.
     *
     * @param nodeId Dense internal node identifier.
     * @return The altitude difference between a node and its parent.
     */
    [[nodiscard]] AltitudeDifference<T> nodeResidue(NodeId nodeId) const { return TreeAltitudeAlgorithms::nodeResidue(tree_, nodeAltitudeSpan(), nodeId); }

    /**
     * @brief Reconstructs an image by assigning each proper part its smallest-node altitude.
     *
     * @return The reconstructed image by assigning each proper part its smallest-node altitude.
     */
    [[nodiscard]] ImagePtr<T> reconstructFromNodeAltitudes() const {
        return TreeAltitudeAlgorithms::reconstructFromNodeAltitudes(tree_, nodeAltitudeSpan(), "ValuedMorphologicalTree::reconstructFromNodeAltitudes");
    }

    /**
     * @brief Reconstructs from arbitrary dense node contributions using the fixed zero baseline.
     * @tparam Contribution Arithmetic contribution and output value type.
     * @param nodeContributions Dense contribution buffer indexed by the internal node-slot domain.
     * @return Image obtained by accumulating contributions from the root to each smallest node.
     */
    template <class Contribution>
        requires(std::is_arithmetic_v<Contribution> && !std::is_same_v<std::remove_cv_t<Contribution>, bool>)
    [[nodiscard]] ImagePtr<Contribution> reconstructFromNodeContributions(std::span<const Contribution> nodeContributions) const {
        return TreeAltitudeAlgorithms::reconstructFromNodeContributions(tree_, nodeContributions,
                                                                         "ValuedMorphologicalTree::reconstructFromNodeContributions");
    }

    /**
     * @brief Exports the current live rooted tree to a new compact Higra parent/altitude representation.
     *
     * @details
     * Attribute buffers are projected through
     * `AttributeComputation::projectNodeValuesToExportedHigra()` so
     * valued-tree export remains limited to topology and altitudes.
     *
     * @return The exported current live rooted tree to a new compact Higra parent/altitude representation.
     */
    [[nodiscard]] std::pair<std::vector<NodeId>, std::vector<T>> exportHigraHierarchy() const {
        return TreeAltitudeAlgorithms::exportHigraHierarchy(tree_, nodeAltitudeSpan());
    }

    /**
     * @brief Opens the only public entrypoint for staged valued-tree edits.
     *
     * @return The opened only public entrypoint for staged valued-tree edits.
     */
    [[nodiscard]] ValuedMorphologicalTreeEditor<T> edit() { return ValuedMorphologicalTreeEditor<T>(*this); }
};

/**
 * @brief Edit-session facade for `ValuedMorphologicalTree`.
 *
 * The valued-tree editor reuses the structural `TreeEditor` for topology while
 * keeping the external altitude buffer as the canonical valued-tree state.
 * `commit()` first validates the topology and then validates monotone altitude.
 */
template <AltitudeValue T> class ValuedMorphologicalTreeEditor {
    friend class ValuedMorphologicalTree<T>;
    /** @brief Grants `detail::beginEstablishedValuedEdit` access to the enclosing type. */
    friend ValuedMorphologicalTreeEditor<T> detail::beginEstablishedValuedEdit<T>(ValuedMorphologicalTree<T>& tree);

  private:
    /** @brief Records altitude changes for rollback of a valued-tree edit. */
    struct AltitudeRollbackJournal {
        /** @brief Original size. */
        std::size_t originalSize = 0;
        /** @brief Captured. */
        TreeEditor::DeltaNodeSet captured;
        /** @brief Dense node identifier of the values. */
        std::vector<std::pair<NodeId, T>> values;
    };

    /** @brief Valued morphological tree. */
    ValuedMorphologicalTree<T>& valuedTree_;
    /** @brief Editor. */
    TreeEditor editor_;
    /** @brief Original altitude size. */
    std::size_t originalAltitudeSize_ = 0;
    /** @brief Altitude rollback journal. */
    std::unique_ptr<AltitudeRollbackJournal> altitudeRollbackJournal_;
    /** @brief Edit revision. */
    std::size_t editRevision_ = 0;
    /** @brief Proven revision. */
    std::optional<std::size_t> provenRevision_;

    /**
     * @brief Opens the underlying topology edit session.
     *
     * @param valuedTree Valued tree.
     * @param transactional Whether the edit records rollback information.
     * @param invariantsEstablishedByConstruction Whether construction already established every edit invariant.
     */
    explicit ValuedMorphologicalTreeEditor(ValuedMorphologicalTree<T>& valuedTree, [[maybe_unused]] bool transactional = false,
                                bool invariantsEstablishedByConstruction = false)
        : valuedTree_(valuedTree), editor_(valuedTree.tree_, invariantsEstablishedByConstruction), originalAltitudeSize_(valuedTree.nodeAltitudes_.size()) {}

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
     * @param nodeId Dense internal node identifier.
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
        altitudeRollbackJournal_->values.emplace_back(nodeId, valuedTree_.nodeAltitudes_[static_cast<std::size_t>(nodeId)]);
    }

    /**
     * @brief Restores altitude journal.
     */
    void restoreAltitudeJournal() noexcept {
        if (!altitudeRollbackJournal_ || !editor_.canRollback()) {
            return;
        }
        valuedTree_.nodeAltitudes_.resize(altitudeRollbackJournal_->originalSize);
        for (const auto& [nodeId, altitude] : altitudeRollbackJournal_->values) {
            valuedTree_.nodeAltitudes_[static_cast<std::size_t>(nodeId)] = altitude;
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
    ValuedMorphologicalTreeEditor(const ValuedMorphologicalTreeEditor&) = delete;
    /**
     * @brief Disables copy assignment.
     */
    ValuedMorphologicalTreeEditor& operator=(const ValuedMorphologicalTreeEditor&) = delete;

    /**
     * @brief Transfers the active valued-tree edit session and rollback journal.
     *
     * @param other Object to compare with or transfer from.
     */
    ValuedMorphologicalTreeEditor(ValuedMorphologicalTreeEditor&& other) noexcept
        : valuedTree_(other.valuedTree_), editor_(std::move(other.editor_)), originalAltitudeSize_(other.originalAltitudeSize_),
          altitudeRollbackJournal_(std::move(other.altitudeRollbackJournal_)), editRevision_(other.editRevision_), provenRevision_(other.provenRevision_) {
        other.altitudeRollbackJournal_.reset();
        other.provenRevision_.reset();
    }

    /**
     * @brief Disables move assignment.
     */
    ValuedMorphologicalTreeEditor& operator=(ValuedMorphologicalTreeEditor&&) = delete;

    /**
     * @brief Restores the altitude journal and closes the valued-tree edit session.
     */
    ~ValuedMorphologicalTreeEditor() { restoreAltitudeJournal(); }

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
            throw std::logic_error("ValuedMorphologicalTreeEditor::rollback is unavailable for the internal journal-free editor.");
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
        TreeAltitudeAlgorithms::validateFiniteAltitudeValue(altitude, 0, "ValuedMorphologicalTreeEditor::createDetachedNode");
        const std::size_t requiredAltitudeSize =
            static_cast<std::size_t>(valuedTree_.tree_.numInternalNodeSlots()) + (valuedTree_.tree_.getNumFreeNodeSlots() == 0 ? 1u : 0u);
        valuedTree_.nodeAltitudes_.reserve(requiredAltitudeSize);
        ensureAltitudeRollbackJournal();
        const NodeId nodeId = editor_.createDetachedNode();
        valuedTree_.nodeAltitudes_.resize(static_cast<size_t>(valuedTree_.tree_.numInternalNodeSlots()), T{});
        captureAltitudeForRollback(nodeId);
        valuedTree_.nodeAltitudes_[static_cast<size_t>(nodeId)] = altitude;
        recordMutation();
        return nodeId;
    }

    /**
     * @brief Sets a live node altitude during a staged topology edit.
     *
     * Monotone order is intentionally checked at valued-tree commit time because
     * intermediate staged topologies may not yet have final parent/child
     * relations.
     *
     * @param nodeId Dense internal node identifier.
     * @param altitude Altitude data indexed by node identifier.
     */
    void setNodeAltitude(NodeId nodeId, T altitude) {
        if (!valuedTree_.tree_.isAlive(nodeId)) {
            throw std::invalid_argument("ValuedMorphologicalTreeEditor::setNodeAltitude requires a live node.");
        }
        TreeAltitudeAlgorithms::validateFiniteAltitudeValue(altitude, static_cast<std::size_t>(nodeId), "ValuedMorphologicalTreeEditor::setNodeAltitude");
        captureAltitudeForRollback(nodeId);
        valuedTree_.nodeAltitudes_[static_cast<size_t>(nodeId)] = altitude;
        editor_.touch(nodeId);
        recordMutation();
    }

    /**
     * @brief Detaches one non-root node through the structural editor.
     *
     * @param nodeId Dense internal node identifier.
     */
    void detach(NodeId nodeId) {
        editor_.detach(nodeId);
        recordMutation();
    }

    /**
     * @brief Reparents one node through the structural editor.
     *
     * @param nodeId Dense internal node identifier.
     * @param newParentId Parent-node value.
     */
    void reparent(NodeId nodeId, NodeId newParentId) {
        editor_.reparent(nodeId, newParentId);
        recordMutation();
    }

    /**
     * @brief Attaches one detached node through the structural editor.
     *
     * @param parentId Identifier of the parent node.
     * @param detachedNodeId Node identifier.
     */
    void attach(NodeId parentId, NodeId detachedNodeId) {
        editor_.attach(parentId, detachedNodeId);
        recordMutation();
    }

    /**
     * @brief Moves all direct children from `sourceId` under `parentId`.
     *
     * @param parentId Identifier of the parent node.
     * @param sourceId Input.
     */
    void moveChildren(NodeId parentId, NodeId sourceId) {
        editor_.moveChildren(parentId, sourceId);
        recordMutation();
    }

    /**
     * @brief Moves one direct proper part between nodes.
     *
     * @param targetNodeId Node identifier.
     * @param sourceNodeId Node identifier.
     * @param pixel Proper-part identifier.
     */
    void movePixelToProperPart(NodeId targetNodeId, NodeId sourceNodeId, PixelId pixel) {
        editor_.movePixelToProperPart(targetNodeId, sourceNodeId, pixel);
        recordMutation();
    }

    /**
     * @brief Moves all direct proper parts between nodes.
     *
     * @param targetNodeId Node identifier.
     * @param sourceNodeId Node identifier.
     */
    void mergeProperParts(NodeId targetNodeId, NodeId sourceNodeId) {
        editor_.mergeProperParts(targetNodeId, sourceNodeId);
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
     * @param nodeId Dense internal node identifier.
     */
    void releaseNode(NodeId nodeId) {
        editor_.releaseNode(nodeId);
        recordMutation();
    }

    /**
     * @brief Promotes one node to the topology root.
     *
     * @param nodeId Dense internal node identifier.
     */
    void setRoot(NodeId nodeId) {
        editor_.setRoot(nodeId);
        recordMutation();
    }

    /**
     * @brief Applies the topology prune helper inside the valued-tree edit session.
     *
     * @param nodeId Dense internal node identifier.
     */
    void pruneNode(NodeId nodeId) {
        editor_.pruneNode(nodeId);
        recordMutation();
    }

    /**
     * @brief Applies the topology merge helper inside the valued-tree edit session.
     *
     * @param nodeId Dense internal node identifier.
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
     * @brief Produces a generic move-only proof for the current valued-tree edit
     * revision.
     *
     * Topology is checked on the mutation delta whenever every primitive is
     * supported. Altitude order is then checked only around touched nodes. A
     * topology fallback also selects complete altitude validation.
     *
     * @return The produced generic move-only proof for the current valued-tree edit revision.
     */
    [[nodiscard]] TreeEditor::IncrementalProof proveIncremental() {
        bool strictAltitudeExcludesCycles = !ValuedMorphologicalTree<T>::skipsMonotoneValidation(valuedTree_.tree_);
        if (!editor_.invariantsEstablishedByConstruction_) {
            for (NodeId node : editor_.touchedNodes_.entries()) {
                if (!valuedTree_.tree_.isAlive(node)) {
                    continue;
                }
                strictAltitudeExcludesCycles =
                    valuedTree_.validateLocalMonotoneNodeAltitudeUpdate(node, valuedTree_.nodeAltitudes_[static_cast<std::size_t>(node)]) && strictAltitudeExcludesCycles;
            }
        }

        auto proof = strictAltitudeExcludesCycles ? editor_.proveIncrementalWithStrictAltitudeAcyclicity() : editor_.proveIncremental();
#ifndef NDEBUG
        // The assertion-enabled oracle runs exactly once, including when the
        // topology proof already fell back to complete validation.
        valuedTree_.validateMonotoneNodeAltitudes();
#else
        if (proof.usedCompleteValidation()) {
            valuedTree_.validateMonotoneNodeAltitudes();
        }
#endif
        provenRevision_ = editRevision_;
        return proof;
    }

    /**
     * @brief Commits the exact valued-tree edit revision represented by `proof`.
     *
     * @param proof Validation proof consumed by the operation.
     */
    void commit(TreeEditor::IncrementalProof&& proof) {
        if (!provenRevision_ || *provenRevision_ != editRevision_) {
            throw std::logic_error("Incremental valued-tree proof is stale or belongs to another edit revision.");
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
            valuedTree_.validateMonotoneNodeAltitudes();
        } catch (const std::exception& ex) {
            return {false, ex.what()};
        } catch (...) {
            return {false, "ValuedMorphologicalTreeEditor monotone-altitude validation failed with an unknown error."};
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
 * @param tree Tree topology.
 * @return Editor for an invariant-preserving valued-tree edit.
 */
template <AltitudeValue T> [[nodiscard]] ValuedMorphologicalTreeEditor<T> beginEstablishedValuedEdit(ValuedMorphologicalTree<T>& tree) {
    return ValuedMorphologicalTreeEditor<T>(tree, false, true);
}

} // namespace detail

} // namespace mmcfilters
