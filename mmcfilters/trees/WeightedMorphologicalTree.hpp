#pragma once

#include "MorphologicalTree.hpp"
#include "TreeAltitudeAlgorithms.hpp"
#include "TreeEditor.hpp"
#include "WeightedTreeView.hpp"

#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcfilters {

template<AltitudeValue T>
class WeightedTreeEditor;

template<AltitudeValue T>
class WeightedMorphologicalTree;


namespace detail {

/**
 * @brief Declares how constructor-provided altitude values are indexed.
 *
 * `WeightedMorphologicalTree` stores altitudes in the dense internal `NodeId`
 * domain. Some builders already produce values in that domain, while static
 * imports from Higra provide values in the original compact Higra domain. This
 * enum makes the required remapping explicit at the constructor boundary.
 */
enum class AltitudeDomain {
    /**
     * @brief Values are indexed directly by internal `NodeId` slots.
     */
    InternalNodeSlots,

    /**
     * @brief Values are indexed by preserved Higra node ids and must be remapped.
     */
    HigraNodeIds
};

/**
 * @brief Lightweight altitude view passed by factory-owned construction paths.
 *
 * The view is consumed synchronously by the tagged constructor and copied into
 * the weighted tree's owned altitude buffer. It does not extend the lifetime of
 * the referenced values.
 */
template<AltitudeValue T>
struct AltitudeInput {
    std::span<const T> values;
    AltitudeDomain domain;
};

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
template<AltitudeValue T>
class WeightedMorphologicalTree {
    friend class WeightedTreeEditor<T>;

public:
    /// Altitude scalar type stored by this weighted tree.
    using altitude_type = T;

    /// Dense altitude buffer indexed by internal `NodeId`.
    using altitude_buffer = std::vector<T>;

private:
    MorphologicalTree tree_;
    altitude_buffer altitude_;

    /**
     * @brief Initializes node altitudes from one direct proper part per node.
     *
     * Image-built component trees create each internal node at an image level,
     * so any direct proper part owned by that node carries the node altitude.
     */
    void assignAltitudeFromDirectProperParts(const ImagePtr<T>& img) {
        altitude_.assign(static_cast<size_t>(tree_.getNumInternalNodeSlots()), T{});
        for (NodeId nodeId : tree_.getAliveNodeIds()) {
            const auto properParts = tree_.getProperParts(nodeId);
            auto it = properParts.begin();
            if (it == properParts.end()) {
                throw std::runtime_error("Cannot infer node altitude from a topology node without direct proper parts.");
            }
            altitude_[static_cast<size_t>(nodeId)] = static_cast<T>((*img)[*it]);
        }
    }

    /**
     * @brief Copies an altitude buffer already indexed by internal node slots.
     */
    void assignInternalAltitude(std::span<const T> altitudeValues) {
        if (altitudeValues.size() != static_cast<size_t>(tree_.getNumInternalNodeSlots())) {
            throw std::invalid_argument("Internal altitude buffer size must match the dense internal-node domain.");
        }
        TreeAltitudeAlgorithms::validateFiniteAltitudeValues(altitudeValues, "WeightedMorphologicalTree internal altitude input");
        altitude_.assign(altitudeValues.begin(), altitudeValues.end());
    }

    /**
     * @brief Remaps imported Higra-node altitudes to internal node slots.
     */
    void importAltitudeFromHigra(std::span<const T> higraAltitude) {
        if (static_cast<size_t>(tree_.getNumHigraNodes()) != higraAltitude.size()) {
            throw std::invalid_argument("Higra altitude buffer size must match the preserved imported Higra hierarchy.");
        }
        TreeAltitudeAlgorithms::validateFiniteAltitudeValues(higraAltitude, "WeightedMorphologicalTree Higra altitude input");
        altitude_.assign(static_cast<size_t>(tree_.getNumInternalNodeSlots()), T{});
        for (NodeId slotId = 0; slotId < tree_.getNumInternalNodeSlots(); ++slotId) {
            const NodeId higraNodeId = tree_.getHigraNodeId(slotId);
            altitude_[static_cast<size_t>(slotId)] = static_cast<T>(higraAltitude[static_cast<size_t>(higraNodeId)]);
        }
    }

    /**
     * @brief Returns whether this tree kind has no global monotone altitude order.
     */
    static bool skipsMonotoneValidation(const MorphologicalTree& tree) noexcept {
        const MorphologicalTreeKind treeType = tree.getTreeType();
        return treeType == MorphologicalTreeKind::TREE_OF_SHAPES || treeType == MorphologicalTreeKind::SELF_DUAL_RESIDUAL_TREE;
    }

    /**
     * @brief Checks the local parent/children altitude constraints for one edit.
     *
     * This keeps `setAltitude()` proportional to the node degree instead of
     * scanning the full tree.
     */
    void validateLocalMonotoneAltitudeUpdate(NodeId nodeId, T value) const {
        if (skipsMonotoneValidation(tree_)) {
            return;
        }

        const bool increasingTowardLeaves = tree_.getTreeType() == MorphologicalTreeKind::MAX_TREE;
        if (!tree_.isRoot(nodeId)) {
            const NodeId parentNodeId = tree_.getNodeParent(nodeId);
            if (parentNodeId == InvalidNode || !tree_.isAlive(parentNodeId)) {
                throw std::runtime_error("Monotone altitude update requires an alive non-root node to have an alive parent.");
            }
            const T parentAltitude = altitude_[static_cast<size_t>(parentNodeId)];
            if (increasingTowardLeaves && parentAltitude > value) {
                throw std::runtime_error("Max-tree altitude update must keep altitude non-decreasing from parent to child.");
            }
            if (!increasingTowardLeaves && parentAltitude < value) {
                throw std::runtime_error("Min-tree altitude update must keep altitude non-increasing from parent to child.");
            }
        }

        for (NodeId childId : tree_.getChildren(nodeId)) {
            const T childAltitude = altitude_[static_cast<size_t>(childId)];
            if (increasingTowardLeaves && value > childAltitude) {
                throw std::runtime_error("Max-tree altitude update must keep altitude non-decreasing from parent to child.");
            }
            if (!increasingTowardLeaves && value < childAltitude) {
                throw std::runtime_error("Min-tree altitude update must keep altitude non-increasing from parent to child.");
            }
        }
    }

private:
    WeightedMorphologicalTree() = delete;

public:
    WeightedMorphologicalTree(const WeightedMorphologicalTree&) = delete;
    WeightedMorphologicalTree& operator=(const WeightedMorphologicalTree&) = delete;

    /// Transfers topology and altitude ownership from another weighted tree.
    WeightedMorphologicalTree(WeightedMorphologicalTree&&) noexcept = default;

    /// Move-assigns topology and altitude ownership from another weighted tree.
    WeightedMorphologicalTree& operator=(WeightedMorphologicalTree&&) noexcept = default;

    /**
     * @brief Creates a weighted tree by inferring node altitudes from an image.
     *
     * The topology must already be built over the same image domain. For each
     * live internal node, the constructor reads one direct proper part owned by
     * that node and copies the corresponding image gray level into the dense
     * altitude buffer.
     *
     * This constructor is tag-protected because image-to-altitude inference is a
     * construction detail owned by the factory.
     */
    WeightedMorphologicalTree(detail::MorphologicalTreeConstructionTag, MorphologicalTree&& topology, ImagePtr<T> img) : tree_(std::move(topology)) {
        if (!img) {
            throw std::invalid_argument("WeightedMorphologicalTree construction requires a non-null image.");
        }
        if (img->getNumRows() <= 0 || img->getNumCols() <= 0 || img->getSize() <= 0) {
            throw std::invalid_argument("WeightedMorphologicalTree construction requires a non-empty 2D image.");
        }

        if (img->getNumRows() != tree_.getNumRowsOfImage() || img->getNumCols() != tree_.getNumColsOfImage() || img->getSize() != tree_.getNumTotalProperParts()) {
            throw std::invalid_argument("WeightedMorphologicalTree construction requires an image domain matching the topology.");
        }

        assignAltitudeFromDirectProperParts(img);
        validateAltitudeBufferShape();
    }

    /**
     * @brief Creates a weighted tree from an explicit altitude input view.
     *
     * The `AltitudeInput::domain` field determines whether values are copied
     * directly as internal-node altitudes or remapped from a preserved Higra
     * hierarchy. This is used by Higra import and SDRT materialization paths.
     *
     * This constructor is tag-protected so callers cannot accidentally bypass
     * the public factory contracts for altitude shape and domain.
     */
    WeightedMorphologicalTree(detail::MorphologicalTreeConstructionTag, MorphologicalTree&& topology, detail::AltitudeInput<T> altitude) : tree_(std::move(topology)) {
        switch (altitude.domain) {
            case detail::AltitudeDomain::InternalNodeSlots:
                assignInternalAltitude(altitude.values);
                break;
            case detail::AltitudeDomain::HigraNodeIds:
                importAltitudeFromHigra(altitude.values);
                break;
        }
        validateAltitudeBufferShape();
    }

    /**
     * @brief Checks that the altitude buffer covers the dense internal-node domain.
     */
    void validateAltitudeBufferShape() const {
        TreeAltitudeAlgorithms::validateAltitudeBufferShape(tree_, altitudeSpan());
    }

    /**
     * @brief Returns read-only access to the owned topology.
     *
     * This accessor intentionally does not provide a mutable `MorphologicalTree`
     * handle. Use `edit()` for staged topology edits or the safe public
     * mutators for complete local changes.
     */
    [[nodiscard]] const MorphologicalTree& topology() const noexcept {
        return tree_;
    }

    /**
     * @brief Returns the dense altitude buffer indexed by internal `NodeId`.
     */
    [[nodiscard]] const altitude_buffer& getAltitudeBuffer() const noexcept {
        return altitude_;
    }

    /**
     * @brief Returns a read-only span over the dense altitude buffer.
     */
    [[nodiscard]] AltitudeSpan<T> altitudeSpan() const noexcept {
        return std::span<const T>(altitude_);
    }

    /**
     * @brief Creates a non-owning weighted view over this owner.
     */
    [[nodiscard]] WeightedTreeView<T> asView() const {
        return WeightedTreeView<T>(tree_, altitudeSpan());
    }

    /**
     * @brief Replaces the owned altitude buffer after full validation.
     *
     * Shape and finite-value checks are followed by global monotone validation
     * for max/min component trees.
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
     * @brief Replaces the owned altitude buffer without checking tree-order monotonicity.
     *
     * This still validates the committed-edit boundary, dense buffer shape, and
     * finite floating-point values. It exists for callers that already maintain
     * the max-tree/min-tree altitude invariant by construction.
     */
    void setAltitudeBufferUnchecked(altitude_buffer altitudeBuffer) {
        tree_.requireNotEditing("WeightedMorphologicalTree::setAltitudeBufferUnchecked");
        if (altitudeBuffer.size() != static_cast<size_t>(tree_.getNumInternalNodeSlots())) {
            throw std::runtime_error("Altitude buffer size must match the dense internal-node domain.");
        }
        TreeAltitudeAlgorithms::validateFiniteAltitudeValues(std::span<const T>(altitudeBuffer), "WeightedMorphologicalTree::setAltitudeBufferUnchecked");
        altitude_ = std::move(altitudeBuffer);
    }

    /**
     * @brief Validates the current altitude buffer against the topology order.
     */
    void validateMonotoneAltitude() const {
        TreeAltitudeAlgorithms::validateMonotoneAltitude(tree_, altitudeSpan());
    }

    /**
     * @brief Returns one live node altitude from the dense buffer.
     */
    [[nodiscard]] T getAltitude(NodeId nodeId) const {
        if (!tree_.isAlive(nodeId) || static_cast<size_t>(nodeId) >= altitude_.size()) {
            throw std::invalid_argument("WeightedMorphologicalTree::getAltitude requires a live internal NodeId.");
        }
        return altitude_[static_cast<size_t>(nodeId)];
    }

    /**
     * @brief Updates one live node altitude with local monotonicity validation.
     */
    void setAltitude(NodeId nodeId, T value) {
        tree_.requireNotEditing("WeightedMorphologicalTree::setAltitude");
        if (!tree_.isAlive(nodeId) || static_cast<size_t>(nodeId) >= altitude_.size()) {
            throw std::invalid_argument("WeightedMorphologicalTree::setAltitude requires a live internal NodeId.");
        }
        TreeAltitudeAlgorithms::validateFiniteAltitudeValue(value, static_cast<std::size_t>(nodeId), "WeightedMorphologicalTree::setAltitude");
        validateLocalMonotoneAltitudeUpdate(nodeId, value);
        altitude_[static_cast<size_t>(nodeId)] = value;
    }

    /**
     * @brief Updates one node altitude without checking tree-order monotonicity.
     *
     * This still validates the committed-edit boundary, live-node access, and
     * finite floating-point values. Staged topology edits should continue to use
     * `WeightedTreeEditor::setNodeAltitude`.
     */
    void setAltitudeUnchecked(NodeId nodeId, T value) {
        tree_.requireNotEditing("WeightedMorphologicalTree::setAltitudeUnchecked");
        if (!tree_.isAlive(nodeId) || static_cast<size_t>(nodeId) >= altitude_.size()) {
            throw std::invalid_argument("WeightedMorphologicalTree::setAltitudeUnchecked requires a live internal NodeId.");
        }
        TreeAltitudeAlgorithms::validateFiniteAltitudeValue(value, static_cast<std::size_t>(nodeId), "WeightedMorphologicalTree::setAltitudeUnchecked");
        altitude_[static_cast<size_t>(nodeId)] = value;
    }

    /**
     * @brief Prunes a complete subtree through the owned topology.
     *
     * The altitude buffer is kept as the canonical dense weighted state; dead
     * slots may retain their old values until a compact export is requested.
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
     */
    void mergeNodeIntoParent(NodeId nodeId) {
        tree_.requireNotEditing("WeightedMorphologicalTree::mergeNodeIntoParent");
        tree_.mergeNodeIntoParent(nodeId);
    }

    /**
     * @brief Returns the altitude difference between a node and its parent.
     */
    [[nodiscard]] AltitudeDiff<T> getNodeResidue(NodeId nodeId) const {
        return TreeAltitudeAlgorithms::getNodeResidue(tree_, altitudeSpan(), nodeId);
    }

    /**
     * @brief Reconstructs an image by assigning each proper part its owner altitude.
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
     */
    [[nodiscard]] std::pair<std::vector<NodeId>, std::vector<T>> exportHigraHierarchy() const {
        return TreeAltitudeAlgorithms::exportHigraHierarchy(tree_, altitudeSpan());
    }

    /**
     * @brief Opens the only public entrypoint for staged weighted edits.
     */
    [[nodiscard("Discarding the editor leaves the weighted tree edit session open")]] WeightedTreeEditor<T> edit() {
        return WeightedTreeEditor<T>(*this);
    }
};

/**
 * @brief Edit-session facade for `WeightedMorphologicalTree`.
 *
 * The weighted editor reuses the structural `TreeEditor` for topology while
 * keeping the external altitude buffer as the canonical weighted state.
 * `commit()` first validates the topology and then validates monotone altitude.
 */
template<AltitudeValue T>
class WeightedTreeEditor {
    friend class WeightedMorphologicalTree<T>;

private:
    WeightedMorphologicalTree<T>& weighted_;
    TreeEditor editor_;

    /**
     * @brief Opens the underlying topology edit session.
     */
    explicit WeightedTreeEditor(WeightedMorphologicalTree<T>& weighted) : weighted_(weighted), editor_(weighted.tree_.edit()) {}

public:

    /**
     * @brief Creates a detached topology node and initializes its altitude.
     *
     * The altitude buffer is reserved before mutating the topology so allocation
     * failure cannot leave a new node without a matching altitude slot.
     */
    [[nodiscard("Discarding a detached node id makes the staged edit impossible to complete safely")]] NodeId createDetachedNode(T altitude = T{}) {
        TreeAltitudeAlgorithms::validateFiniteAltitudeValue(altitude, 0, "WeightedTreeEditor::createDetachedNode");
        const std::size_t requiredAltitudeSize = static_cast<std::size_t>(weighted_.tree_.getNumInternalNodeSlots()) + (weighted_.tree_.getNumFreeNodeSlots() == 0 ? 1u : 0u);
        weighted_.altitude_.reserve(requiredAltitudeSize);
        const NodeId nodeId = editor_.createDetachedNode();
        weighted_.altitude_.resize(static_cast<size_t>(weighted_.tree_.getNumInternalNodeSlots()), T{});
        weighted_.altitude_[static_cast<size_t>(nodeId)] = altitude;
        return nodeId;
    }

    /**
     * @brief Sets a live node altitude during a staged topology edit.
     *
     * Monotone order is intentionally checked at weighted commit time because
     * intermediate staged topologies may not yet have final parent/child
     * relations.
     */
    void setNodeAltitude(NodeId nodeId, T altitude) {
        if (!weighted_.tree_.isAlive(nodeId)) {
            throw std::invalid_argument("WeightedTreeEditor::setNodeAltitude requires a live node.");
        }
        TreeAltitudeAlgorithms::validateFiniteAltitudeValue(altitude, static_cast<std::size_t>(nodeId), "WeightedTreeEditor::setNodeAltitude");
        weighted_.altitude_[static_cast<size_t>(nodeId)] = altitude;
    }

    /**
     * @brief Detaches one non-root node through the structural editor.
     */
    void detach(NodeId nodeId) {
        editor_.detach(nodeId);
    }

    /**
     * @brief Reparents one node through the structural editor.
     */
    void reparent(NodeId nodeId, NodeId newParentId) {
        editor_.reparent(nodeId, newParentId);
    }

    /**
     * @brief Attaches one detached node through the structural editor.
     */
    void attach(NodeId parentId, NodeId detachedNodeId) {
        editor_.attach(parentId, detachedNodeId);
    }

    /**
     * @brief Moves all direct children from `sourceId` under `parentId`.
     */
    void moveChildren(NodeId parentId, NodeId sourceId) {
        editor_.moveChildren(parentId, sourceId);
    }

    /**
     * @brief Moves one direct proper part between nodes.
     */
    void moveProperPart(NodeId targetNodeId, NodeId sourceNodeId, NodeId properPartId) {
        editor_.moveProperPart(targetNodeId, sourceNodeId, properPartId);
    }

    /**
     * @brief Moves all direct proper parts between nodes.
     */
    void moveProperParts(NodeId targetNodeId, NodeId sourceNodeId) {
        editor_.moveProperParts(targetNodeId, sourceNodeId);
    }

    /**
     * @brief Detaches a direct child and optionally releases an empty node slot.
     */
    void removeChild(NodeId parentNodeId, NodeId childId, bool releaseNodeFlag) {
        editor_.removeChild(parentNodeId, childId, releaseNodeFlag);
    }

    /**
     * @brief Releases an empty detached node slot.
     */
    void releaseNode(NodeId nodeId) {
        editor_.releaseNode(nodeId);
    }

    /**
     * @brief Promotes one node to the topology root.
     */
    void setRoot(NodeId nodeId) {
        editor_.setRoot(nodeId);
    }

    /**
     * @brief Applies the topology prune helper inside the weighted edit session.
     */
    void pruneNode(NodeId nodeId) {
        editor_.pruneNode(nodeId);
    }

    /**
     * @brief Applies the topology merge helper inside the weighted edit session.
     */
    void mergeNodeIntoParent(NodeId nodeId) {
        editor_.mergeNodeIntoParent(nodeId);
    }

    /**
     * @brief Returns whether the structural edit still has detached live nodes.
     */
    [[nodiscard]] bool hasDetachedAliveNodes() const noexcept {
        return editor_.hasDetachedAliveNodes();
    }

    /**
     * @brief Validates only the staged topology.
     *
     * Altitude monotonicity is checked by `validateAndCommit()` / `commit()`.
     */
    [[nodiscard]] TreeValidationResult validate() const noexcept {
        return editor_.validate();
    }

    /**
     * @brief Validates topology, validates altitude order, then closes the edit.
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
        editor_.commitUnchecked();
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

    /**
     * @brief Closes the weighted edit without topology or altitude validation.
     */
    void commitUnchecked() noexcept {
        editor_.commitUnchecked();
    }
};

} // namespace mmcfilters
