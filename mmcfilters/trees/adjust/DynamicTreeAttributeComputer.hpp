#pragma once

#include "../WeightedMorphologicalTree.hpp"
#include "../../utils/Image.hpp"

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

namespace mmcfilters::adjust {

/**
 * @brief Selects the scalar value extracted from a dynamic bounding-box summary.
 *
 * All three measures are increasing attributes when applied to a component
 * support: enlarging the component cannot decrease its width, height, or
 * diagonal length.
 */
enum class BoundingBoxMeasure {
    /// Width of the subtree support bounding box.
    WIDTH,
    /// Height of the subtree support bounding box.
    HEIGHT,
    /// Euclidean diagonal length of the subtree support bounding box.
    DIAGONAL_LENGTH
};

/**
 * @brief Common protocol for attributes maintained during local tree adjustment.
 *
 * The dual min/max-tree adjuster performs structural edits that are local to
 * the affected component. Recomputing every node attribute after each edit
 * would change the practical cost profile of the connected alternating
 * sequential filter, so attributes used by the CASF follow the same incremental
 * contract as the Higra implementation:
 *
 * - `preProcessing` initializes the direct contribution of one node;
 * - `mergeProcessing` accumulates an already-current child contribution into
 *   its parent;
 * - `postProcessing` materializes the final scalar value;
 * - mutation hooks update any auxiliary local state when proper parts move or
 *   when a node slot is released.
 *
 * Buffers are indexed in this project's dense internal `NodeId` space, not in
 * Higra's global leaf+node id space. This is the main representation
 * difference relative to the source algorithm and keeps the implementation
 * compatible with `WeightedMorphologicalTree<T>`.
 */
template <AltitudeValue T> class DynamicTreeAttributeComputer {
  public:
    /// Dense per-node attribute buffer used by dynamic adjustment computers.
    using buffer_type = std::vector<double>;

    /**
     * @brief Destroys a dynamic attribute computer through the protocol base.
     */
    virtual ~DynamicTreeAttributeComputer() = default;

    /**
     * @brief Resizes an attribute buffer to the full internal node-id space.
     *
     * Released slots remain addressable while the dynamic adjustment is running,
     * so buffers are sized by the number of allocated internal node slots rather
     * than by the number of currently alive nodes.
     *
     * @param tree Weighted tree whose internal node-slot domain defines the
     * buffer size.
     * @param buffer Output buffer resized to the full internal node-slot count.
     */
    virtual void resize(const WeightedMorphologicalTree<T>& tree, buffer_type& buffer) const {
        buffer.resize(static_cast<size_t>(tree.topology().getNumInternalNodeSlots()), 0.0);
    }

    /**
     * @brief Initializes the direct contribution of one node before child merges.
     *
     * @param nodeId Live dense internal node id being recomputed.
     * @param tree Current mutable weighted tree state.
     * @param buffer Dense per-node output buffer.
     */
    virtual void preProcessing(NodeId nodeId, const WeightedMorphologicalTree<T>& tree, buffer_type& buffer) const = 0;

    /**
     * @brief Accumulates an already-current child contribution into its parent.
     *
     * @param parentId Live dense internal parent node id.
     * @param childId Live dense internal child node id whose buffer entry is
     * already current.
     * @param tree Current mutable weighted tree state.
     * @param buffer Dense per-node output buffer.
     */
    virtual void mergeProcessing(NodeId parentId, NodeId childId, const WeightedMorphologicalTree<T>& tree, buffer_type& buffer) const = 0;

    /**
     * @brief Materializes the final scalar value for one node after all child merges.
     *
     * @param nodeId Live dense internal node id being finalized.
     * @param tree Current mutable weighted tree state.
     * @param buffer Dense per-node output buffer.
     */
    virtual void postProcessing(NodeId nodeId, const WeightedMorphologicalTree<T>& tree, buffer_type& buffer) const = 0;

    /**
     * @brief Incremental hook called after all direct proper parts move from one node to another.
     * @details Default implementation is a no-op so concrete computers override
     * only the structural events they need.
     *
     * The parameters are source node id, destination node id, and current tree
     * state. Implementations that maintain direct-proper-part caches can use
     * this event to move or invalidate cached summaries without a full rebuild.
     */
    virtual void onMoveProperParts(NodeId, NodeId, const WeightedMorphologicalTree<T>&) const {}

    /**
     * @brief Incremental hook called after one proper part moves between nodes.
     * @details Default implementation is a no-op so concrete computers override
     * only the structural events they need.
     *
     * The parameters are proper-part id, source node id, destination node id,
     * and current tree state.
     */
    virtual void onMoveProperPart(NodeId, NodeId, NodeId, const WeightedMorphologicalTree<T>&) const {}

    /**
     * @brief Incremental hook called when one node slot is released from the live tree.
     * @details Default implementation is a no-op so concrete computers override
     * only the structural events they need.
     *
     * Implementations should clear or invalidate any auxiliary state indexed by
     * the released dense internal node id.
     */
    virtual void onNodeRemoved(NodeId, const WeightedMorphologicalTree<T>&) const {}

    /**
     * @brief Computes the attribute for the full current tree in post-order.
     *
     * This is used to bootstrap the buffers before the first adjustment step
     * and as a reference-compatible full computation. Adjustment steps later
     * refresh only the nodes marked by local structural changes.
     *
     * @param tree Current mutable weighted tree state.
     * @param buffer Dense output buffer. It is resized before computation.
     */
    void computeAttribute(const WeightedMorphologicalTree<T>& tree, buffer_type& buffer) const {
        resize(tree, buffer);
        const MorphologicalTree& topology = tree.topology();
        for (NodeId nodeId : topology.getPostOrderNodes()) {
            preProcessing(nodeId, tree, buffer);
            for (NodeId childId : topology.getChildren(nodeId)) {
                mergeProcessing(nodeId, childId, tree, buffer);
            }
            postProcessing(nodeId, tree, buffer);
        }
    }

    /**
     * @brief Recomputes one node assuming all direct children are already up to date.
     *
     * @param tree Current mutable weighted tree state.
     * @param nodeId Live dense internal node id to refresh.
     * @param buffer Dense output buffer already sized for the tree.
     */
    void computeAttributeOnNode(const WeightedMorphologicalTree<T>& tree, NodeId nodeId, buffer_type& buffer) const {
        preProcessing(nodeId, tree, buffer);
        for (NodeId childId : tree.topology().getChildren(nodeId)) {
            mergeProcessing(nodeId, childId, tree, buffer);
        }
        postProcessing(nodeId, tree, buffer);
    }
};

/**
 * @brief Incremental area attribute.
 *
 * Area is increasing and is defined as the number of proper parts in the full
 * support of a node. The local contribution is the number of direct proper
 * parts; child values are then accumulated bottom-up. The single-node refresh
 * contract is straightforward: once child areas are current, recomputing a
 * parent requires only its direct proper-part count and the cached values of
 * its direct children.
 */
template <AltitudeValue T> class DynamicAreaAttributeComputer : public DynamicTreeAttributeComputer<T> {
  private:
    /** @brief Defines the `base_t` alias used by the component. */
    using base_t = DynamicTreeAttributeComputer<T>;

  public:
    /// Dense per-node area buffer inherited from the dynamic attribute protocol.
    using buffer_type = typename base_t::buffer_type;

  public:
    /**
     * @brief Initializes one node area from its direct proper-part count.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @param tree Tree topology used by the operation.
     * @param buffer Buffer read or written by the operation.
     */
    void preProcessing(NodeId nodeId, const WeightedMorphologicalTree<T>& tree, buffer_type& buffer) const override {
        buffer[static_cast<size_t>(nodeId)] = static_cast<double>(tree.topology().getNumProperParts(nodeId));
    }

    /**
     * @brief Adds an already-current child area to its parent.
     *
     * @param parentId Identifier of the parent node.
     * @param childId Identifier of the child node.
     * @param buffer Buffer read or written by the operation.
     */
    void mergeProcessing(NodeId parentId, NodeId childId, const WeightedMorphologicalTree<T>&, buffer_type& buffer) const override {
        buffer[static_cast<size_t>(parentId)] += buffer[static_cast<size_t>(childId)];
    }

    /**
     * @brief Area has no finalization step beyond child accumulation.
     */
    void postProcessing(NodeId, const WeightedMorphologicalTree<T>&, buffer_type&) const override {}
};

/**
 * @brief Incremental bounding-box scalar attribute.
 *
 * The computer keeps two levels of cached state:
 *
 * - `local_`: bounding box of the node's direct proper parts, including enough
 *   extremum multiplicities to update most single-pixel removals without a
 *   scan;
 * - `subtree_`: bounding box after child aggregation for the current
 *   `pre/merge/post` pass.
 *
 * Single proper-part moves update the receiver eagerly. The donor is kept exact
 * unless the moved pixel exhausts one extremum, in which case its local summary
 * is marked dirty and rebuilt lazily on the next access. This preserves the
 * local-update structure of the Higra version while adapting coordinates to
 * this project's row-major image domain.
 */
template <AltitudeValue T> class DynamicBoundingBoxAttributeComputer : public DynamicTreeAttributeComputer<T> {
  private:
    /** @brief Defines the `base_t` alias used by the component. */
    using base_t = DynamicTreeAttributeComputer<T>;

  public:
    /// Dense per-node bounding-box attribute buffer inherited from the dynamic protocol.
    using buffer_type = typename base_t::buffer_type;

  private:
    /**
     * @brief Cached subtree box produced during the current bottom-up reduction.
     */
    struct BoxState {
        /** @brief Stores the xmin. */
        int xmin = 0;
        /** @brief Stores the xmax. */
        int xmax = -1;
        /** @brief Stores the ymin. */
        int ymin = 0;
        /** @brief Stores the ymax. */
        int ymax = -1;
        /** @brief Indicates whether the old snapshot was empty. */
        bool empty = true;
    };

    /**
     * @brief Cached box of direct proper parts plus metadata for local updates.
     *
     * Extremum multiplicities allow most single-pixel removals to be handled
     * without scanning the node's whole proper-part list. If the moved pixel was
     * the last occurrence of one extremum, the summary is marked dirty and is
     * rebuilt lazily before the next bottom-up reduction.
     */
    struct LocalBoxState {
        /** @brief Stores the xmin. */
        int xmin = 0;
        /** @brief Stores the xmax. */
        int xmax = -1;
        /** @brief Stores the ymin. */
        int ymin = 0;
        /** @brief Stores the ymax. */
        int ymax = -1;
        /** @brief Stores the xmin count. */
        int xminCount = 0;
        /** @brief Stores the xmax count. */
        int xmaxCount = 0;
        /** @brief Stores the ymin count. */
        int yminCount = 0;
        /** @brief Stores the ymax count. */
        int ymaxCount = 0;
        /** @brief Stores the proper part count. */
        int properPartCount = 0;
        /** @brief Indicates whether the new snapshot is empty. */
        bool empty = true;
        /** @brief Indicates whether the node requires attribute recomputation. */
        bool dirty = true;
    };

    /** @brief Stores the measure. */
    BoundingBoxMeasure measure_ = BoundingBoxMeasure::DIAGONAL_LENGTH;
    /** @brief Stores the local. */
    mutable std::vector<LocalBoxState> local_;
    /** @brief Stores the subtree. */
    mutable std::vector<BoxState> subtree_;

    /**
     * @brief Restores an empty local box and clears all extremum multiplicities.
     *
     * @param local Local state accumulated by the operation.
     */
    static void resetLocalBox(LocalBoxState& local) {
        local.xmin = 0;
        local.xmax = -1;
        local.ymin = 0;
        local.ymax = -1;
        local.xminCount = 0;
        local.xmaxCount = 0;
        local.yminCount = 0;
        local.ymaxCount = 0;
        local.properPartCount = 0;
        local.empty = true;
    }

    /**
     * @brief Resets one node-local cache entry and marks it synchronized.
     *
     * @param nodeId Identifier of the node used by the operation.
     */
    void resetLocalSummary(NodeId nodeId) const {
        auto& local = local_[static_cast<size_t>(nodeId)];
        resetLocalBox(local);
        local.dirty = false;
    }

    /**
     * @brief Clears the subtree cache associated with one node.
     *
     * @param nodeId Identifier of the node used by the operation.
     */
    void resetSubtreeSummary(NodeId nodeId) const {
        auto& subtree = subtree_[static_cast<size_t>(nodeId)];
        subtree.xmin = 0;
        subtree.xmax = -1;
        subtree.ymin = 0;
        subtree.ymax = -1;
        subtree.empty = true;
    }

    /**
     * @brief Enlarges a local box with one proper part and updates extremum counts.
     *
     * @param local Local state accumulated by the operation.
     * @param pixelId Pixel identifier used by the operation.
     * @param numCols Number of columns in the domain.
     */
    void expandLocalBoxWithPixel(LocalBoxState& local, NodeId pixelId, int numCols) const {
        const auto [y, x] = ImageUtils::to2D(pixelId, numCols);
        if (local.empty) {
            local.xmin = x;
            local.xmax = x;
            local.ymin = y;
            local.ymax = y;
            local.xminCount = 1;
            local.xmaxCount = 1;
            local.yminCount = 1;
            local.ymaxCount = 1;
            local.empty = false;
            return;
        }

        if (x < local.xmin) {
            local.xmin = x;
            local.xminCount = 1;
        } else if (x == local.xmin) {
            ++local.xminCount;
        }

        if (x > local.xmax) {
            local.xmax = x;
            local.xmaxCount = 1;
        } else if (x == local.xmax) {
            ++local.xmaxCount;
        }

        if (y < local.ymin) {
            local.ymin = y;
            local.yminCount = 1;
        } else if (y == local.ymin) {
            ++local.yminCount;
        }

        if (y > local.ymax) {
            local.ymax = y;
            local.ymaxCount = 1;
        } else if (y == local.ymax) {
            ++local.ymaxCount;
        }
    }

    /**
     * @brief Rebuilds the local summary from the current proper parts of the node.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @param tree Tree topology used by the operation.
     */
    void rebuildLocalBox(NodeId nodeId, const WeightedMorphologicalTree<T>& tree) const {
        auto& local = local_[static_cast<size_t>(nodeId)];
        resetLocalBox(local);
        const int numCols = tree.topology().getNumColsOfGridDomain2D();
        for (NodeId pixelId : tree.topology().getProperParts(nodeId)) {
            expandLocalBoxWithPixel(local, pixelId, numCols);
        }
        local.properPartCount = tree.topology().getNumProperParts(nodeId);
        local.dirty = false;
    }

    /**
     * @brief Ensures the node-local summary matches the current topology.
     * @details Dirty summaries are rebuilt lazily; clean summaries are also
     * rebuilt if their cached direct proper-part count no longer matches the
     * tree, which covers bulk moves.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @param tree Tree topology used by the operation.
     */
    void ensureLocalSummary(NodeId nodeId, const WeightedMorphologicalTree<T>& tree) const {
        auto& local = local_[static_cast<size_t>(nodeId)];
        const int properPartCount = tree.topology().getNumProperParts(nodeId);
        if (!local.dirty && local.properPartCount == properPartCount) {
            return;
        }
        rebuildLocalBox(nodeId, tree);
    }

    /**
     * @brief Initializes the subtree box from the current direct proper-part box.
     *
     * @param nodeId Identifier of the node used by the operation.
     */
    void copyLocalToSubtree(NodeId nodeId) const {
        const auto& local = local_[static_cast<size_t>(nodeId)];
        auto& subtree = subtree_[static_cast<size_t>(nodeId)];
        subtree.xmin = local.xmin;
        subtree.xmax = local.xmax;
        subtree.ymin = local.ymin;
        subtree.ymax = local.ymax;
        subtree.empty = local.empty;
    }

    /**
     * @brief Enlarges one subtree box with another already-current subtree box.
     *
     * @param target Destination object or value.
     * @param source Source object or value.
     */
    static void mergeSubtreeStates(BoxState& target, const BoxState& source) {
        if (source.empty) {
            return;
        }
        if (target.empty) {
            target = source;
            return;
        }
        target.xmin = std::min(target.xmin, source.xmin);
        target.xmax = std::max(target.xmax, source.xmax);
        target.ymin = std::min(target.ymin, source.ymin);
        target.ymax = std::max(target.ymax, source.ymax);
        target.empty = false;
    }

    /**
     * @brief Merges two direct-proper-part local boxes, including extremum counts.
     *
     * @param target Destination object or value.
     * @param source Source object or value.
     */
    static void mergeLocalBoxes(LocalBoxState& target, const LocalBoxState& source) {
        if (source.empty) {
            return;
        }
        if (target.empty) {
            target = source;
            return;
        }

        if (source.xmin < target.xmin) {
            target.xmin = source.xmin;
            target.xminCount = source.xminCount;
        } else if (source.xmin == target.xmin) {
            target.xminCount += source.xminCount;
        }

        if (source.xmax > target.xmax) {
            target.xmax = source.xmax;
            target.xmaxCount = source.xmaxCount;
        } else if (source.xmax == target.xmax) {
            target.xmaxCount += source.xmaxCount;
        }

        if (source.ymin < target.ymin) {
            target.ymin = source.ymin;
            target.yminCount = source.yminCount;
        } else if (source.ymin == target.ymin) {
            target.yminCount += source.yminCount;
        }

        if (source.ymax > target.ymax) {
            target.ymax = source.ymax;
            target.ymaxCount = source.ymaxCount;
        } else if (source.ymax == target.ymax) {
            target.ymaxCount += source.ymaxCount;
        }

        target.properPartCount += source.properPartCount;
        target.empty = false;
        target.dirty = false;
    }

  public:
    /**
     * @brief Creates a bounding-box computer returning the requested scalar measure.
     *
     * @param measure Scalar projection materialized from the accumulated
     * bounding-box state.
     */
    explicit DynamicBoundingBoxAttributeComputer(BoundingBoxMeasure measure = BoundingBoxMeasure::DIAGONAL_LENGTH) : measure_(measure) {}

    /**
     * @brief Resizes public and auxiliary buffers to the current tree slot space.
     *
     * Auxiliary direct/local and subtree summaries are indexed by the same
     * dense internal node-id space as the public output buffer.
     *
     * @param tree Tree topology used by the operation.
     * @param buffer Buffer read or written by the operation.
     */
    void resize(const WeightedMorphologicalTree<T>& tree, buffer_type& buffer) const override {
        base_t::resize(tree, buffer);
        const size_t size = static_cast<size_t>(tree.topology().getNumInternalNodeSlots());
        local_.resize(size);
        subtree_.resize(size);
    }

    /**
     * @brief Initializes the subtree box of one node from its direct proper parts.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @param tree Tree topology used by the operation.
     */
    void preProcessing(NodeId nodeId, const WeightedMorphologicalTree<T>& tree, buffer_type&) const override {
        ensureLocalSummary(nodeId, tree);
        copyLocalToSubtree(nodeId);
    }

    /**
     * @brief Accumulates a child subtree box into its parent subtree box.
     *
     * @param parentId Identifier of the parent node.
     * @param childId Identifier of the child node.
     */
    void mergeProcessing(NodeId parentId, NodeId childId, const WeightedMorphologicalTree<T>&, buffer_type&) const override {
        mergeSubtreeStates(subtree_[static_cast<size_t>(parentId)], subtree_[static_cast<size_t>(childId)]);
    }

    /**
     * @brief Converts the accumulated subtree box into the configured scalar measure.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @param buffer Buffer read or written by the operation.
     */
    void postProcessing(NodeId nodeId, const WeightedMorphologicalTree<T>&, buffer_type& buffer) const override {
        const auto& subtree = subtree_[static_cast<size_t>(nodeId)];
        if (subtree.empty) {
            buffer[static_cast<size_t>(nodeId)] = 0.0;
            return;
        }

        const double width = static_cast<double>(subtree.xmax - subtree.xmin + 1);
        const double height = static_cast<double>(subtree.ymax - subtree.ymin + 1);
        switch (measure_) {
        case BoundingBoxMeasure::WIDTH:
            buffer[static_cast<size_t>(nodeId)] = width;
            break;
        case BoundingBoxMeasure::HEIGHT:
            buffer[static_cast<size_t>(nodeId)] = height;
            break;
        case BoundingBoxMeasure::DIAGONAL_LENGTH:
            buffer[static_cast<size_t>(nodeId)] = std::sqrt(width * width + height * height);
            break;
        }
    }

    /**
     * @brief Updates local boxes after all proper parts move from `sourceId` to `targetId`.
     *
     * @param targetId Destination represented by `targetId`.
     * @param sourceId Input represented by `sourceId`.
     * @param tree Tree topology used by the operation.
     */
    void onMoveProperParts(NodeId targetId, NodeId sourceId, const WeightedMorphologicalTree<T>& tree) const override {
        ensureLocalSummary(targetId, tree);
        ensureLocalSummary(sourceId, tree);
        mergeLocalBoxes(local_[static_cast<size_t>(targetId)], local_[static_cast<size_t>(sourceId)]);
        resetLocalSummary(sourceId);
    }

    /**
     * @brief Updates local boxes after one proper part moves between nodes.
     * @details The target box is enlarged eagerly. The source box remains exact
     * unless the removed pixel exhausts one extremum; in that case the source
     * summary is marked dirty and rebuilt lazily.
     *
     * @param targetId Destination represented by `targetId`.
     * @param sourceId Input represented by `sourceId`.
     * @param pixelId Pixel identifier used by the operation.
     * @param tree Tree topology used by the operation.
     */
    void onMoveProperPart(NodeId targetId, NodeId sourceId, NodeId pixelId, const WeightedMorphologicalTree<T>& tree) const override {
        ensureLocalSummary(targetId, tree);
        if (sourceId != InvalidNode) {
            ensureLocalSummary(sourceId, tree);
        }

        const int numCols = tree.topology().getNumColsOfGridDomain2D();
        expandLocalBoxWithPixel(local_[static_cast<size_t>(targetId)], pixelId, numCols);
        local_[static_cast<size_t>(targetId)].properPartCount += 1;
        local_[static_cast<size_t>(targetId)].dirty = false;

        if (sourceId == InvalidNode) {
            return;
        }

        auto& source = local_[static_cast<size_t>(sourceId)];
        if (source.properPartCount <= 1) {
            resetLocalSummary(sourceId);
            return;
        }

        const auto [y, x] = ImageUtils::to2D(pixelId, numCols);
        --source.properPartCount;

        bool exhaustsXmin = false;
        bool exhaustsXmax = false;
        bool exhaustsYmin = false;
        bool exhaustsYmax = false;

        if (!source.dirty) {
            if (x == source.xmin) {
                if (source.xminCount <= 1) {
                    exhaustsXmin = true;
                } else {
                    --source.xminCount;
                }
            }
            if (x == source.xmax) {
                if (source.xmaxCount <= 1) {
                    exhaustsXmax = true;
                } else {
                    --source.xmaxCount;
                }
            }
            if (y == source.ymin) {
                if (source.yminCount <= 1) {
                    exhaustsYmin = true;
                } else {
                    --source.yminCount;
                }
            }
            if (y == source.ymax) {
                if (source.ymaxCount <= 1) {
                    exhaustsYmax = true;
                } else {
                    --source.ymaxCount;
                }
            }
        }

        source.dirty = source.dirty || exhaustsXmin || exhaustsXmax || exhaustsYmin || exhaustsYmax;
    }

    /**
     * @brief Clears auxiliary summaries associated with a released node slot.
     *
     * @param nodeId Identifier of the node used by the operation.
     */
    void onNodeRemoved(NodeId nodeId, const WeightedMorphologicalTree<T>&) const override {
        resetLocalSummary(nodeId);
        resetSubtreeSummary(nodeId);
    }
};

} // namespace mmcfilters::adjust
