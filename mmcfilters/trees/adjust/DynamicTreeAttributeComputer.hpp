#pragma once

#include "../WeightedMorphologicalTree.hpp"

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

namespace mmcfilters::adjust {

/**
 * @brief Selects the scalar value extracted from a dynamic bounding-box summary.
 */
enum class BoundingBoxMeasure {
    WIDTH,
    HEIGHT,
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
 * compatible with `WeightedMorphologicalTree`.
 */
class DynamicTreeAttributeComputer {
public:
    using buffer_type = std::vector<double>;

    virtual ~DynamicTreeAttributeComputer() = default;

    virtual void resize(const WeightedMorphologicalTree& tree, buffer_type& buffer) const {
        buffer.resize(static_cast<size_t>(tree.topology().getNumInternalNodeSlots()), 0.0);
    }

    virtual void preProcessing(NodeId nodeId, const WeightedMorphologicalTree& tree, buffer_type& buffer) const = 0;
    virtual void mergeProcessing(NodeId parentId, NodeId childId, const WeightedMorphologicalTree& tree, buffer_type& buffer) const = 0;
    virtual void postProcessing(NodeId nodeId, const WeightedMorphologicalTree& tree, buffer_type& buffer) const = 0;

    virtual void onMoveProperParts(NodeId, NodeId, const WeightedMorphologicalTree&) const {}
    virtual void onMoveProperPart(NodeId, NodeId, NodeId, const WeightedMorphologicalTree&) const {}
    virtual void onNodeRemoved(NodeId, const WeightedMorphologicalTree&) const {}

    /**
     * @brief Computes the attribute for the full current tree in post-order.
     *
     * This is used to bootstrap the buffers before the first adjustment step
     * and as a reference-compatible full computation. Adjustment steps later
     * refresh only the nodes marked by local structural changes.
     */
    void computeAttribute(const WeightedMorphologicalTree& tree, buffer_type& buffer) const {
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
     */
    void computeAttributeOnNode(const WeightedMorphologicalTree& tree, NodeId nodeId, buffer_type& buffer) const {
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
 * parts; child values are then accumulated bottom-up.
 */
class DynamicAreaAttributeComputer : public DynamicTreeAttributeComputer {
public:
    void preProcessing(NodeId nodeId, const WeightedMorphologicalTree& tree, buffer_type& buffer) const override {
        buffer[static_cast<size_t>(nodeId)] = static_cast<double>(tree.topology().getNumProperParts(nodeId));
    }

    void mergeProcessing(NodeId parentId, NodeId childId, const WeightedMorphologicalTree&, buffer_type& buffer) const override {
        buffer[static_cast<size_t>(parentId)] += buffer[static_cast<size_t>(childId)];
    }

    void postProcessing(NodeId, const WeightedMorphologicalTree&, buffer_type&) const override {}
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
class DynamicBoundingBoxAttributeComputer : public DynamicTreeAttributeComputer {
private:
    struct BoxState {
        int xmin = 0;
        int xmax = -1;
        int ymin = 0;
        int ymax = -1;
        bool empty = true;
    };

    struct LocalBoxState {
        int xmin = 0;
        int xmax = -1;
        int ymin = 0;
        int ymax = -1;
        int xminCount = 0;
        int xmaxCount = 0;
        int yminCount = 0;
        int ymaxCount = 0;
        int properPartCount = 0;
        bool empty = true;
        bool dirty = true;
    };

    BoundingBoxMeasure measure_ = BoundingBoxMeasure::DIAGONAL_LENGTH;
    mutable std::vector<LocalBoxState> local_;
    mutable std::vector<BoxState> subtree_;

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

    void resetLocalSummary(NodeId nodeId) const {
        auto& local = local_[static_cast<size_t>(nodeId)];
        resetLocalBox(local);
        local.dirty = false;
    }

    void resetSubtreeSummary(NodeId nodeId) const {
        auto& subtree = subtree_[static_cast<size_t>(nodeId)];
        subtree.xmin = 0;
        subtree.xmax = -1;
        subtree.ymin = 0;
        subtree.ymax = -1;
        subtree.empty = true;
    }

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

    void rebuildLocalBox(NodeId nodeId, const WeightedMorphologicalTree& tree) const {
        auto& local = local_[static_cast<size_t>(nodeId)];
        resetLocalBox(local);
        const int numCols = tree.topology().getNumColsOfImage();
        for (NodeId pixelId : tree.topology().getProperParts(nodeId)) {
            expandLocalBoxWithPixel(local, pixelId, numCols);
        }
        local.properPartCount = tree.topology().getNumProperParts(nodeId);
        local.dirty = false;
    }

    void ensureLocalSummary(NodeId nodeId, const WeightedMorphologicalTree& tree) const {
        auto& local = local_[static_cast<size_t>(nodeId)];
        const int properPartCount = tree.topology().getNumProperParts(nodeId);
        if (!local.dirty && local.properPartCount == properPartCount) {
            return;
        }
        rebuildLocalBox(nodeId, tree);
    }

    void copyLocalToSubtree(NodeId nodeId) const {
        const auto& local = local_[static_cast<size_t>(nodeId)];
        auto& subtree = subtree_[static_cast<size_t>(nodeId)];
        subtree.xmin = local.xmin;
        subtree.xmax = local.xmax;
        subtree.ymin = local.ymin;
        subtree.ymax = local.ymax;
        subtree.empty = local.empty;
    }

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
    explicit DynamicBoundingBoxAttributeComputer(BoundingBoxMeasure measure = BoundingBoxMeasure::DIAGONAL_LENGTH)
        : measure_(measure) {}

    void resize(const WeightedMorphologicalTree& tree, buffer_type& buffer) const override {
        DynamicTreeAttributeComputer::resize(tree, buffer);
        const size_t size = static_cast<size_t>(tree.topology().getNumInternalNodeSlots());
        local_.resize(size);
        subtree_.resize(size);
    }

    void preProcessing(NodeId nodeId, const WeightedMorphologicalTree& tree, buffer_type&) const override {
        ensureLocalSummary(nodeId, tree);
        copyLocalToSubtree(nodeId);
    }

    void mergeProcessing(NodeId parentId, NodeId childId, const WeightedMorphologicalTree&, buffer_type&) const override {
        mergeSubtreeStates(subtree_[static_cast<size_t>(parentId)], subtree_[static_cast<size_t>(childId)]);
    }

    void postProcessing(NodeId nodeId, const WeightedMorphologicalTree&, buffer_type& buffer) const override {
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

    void onMoveProperParts(NodeId targetId, NodeId sourceId, const WeightedMorphologicalTree& tree) const override {
        ensureLocalSummary(targetId, tree);
        ensureLocalSummary(sourceId, tree);
        mergeLocalBoxes(local_[static_cast<size_t>(targetId)], local_[static_cast<size_t>(sourceId)]);
        resetLocalSummary(sourceId);
    }

    void onMoveProperPart(NodeId targetId, NodeId sourceId, NodeId pixelId, const WeightedMorphologicalTree& tree) const override {
        ensureLocalSummary(targetId, tree);
        if (sourceId != InvalidNode) {
            ensureLocalSummary(sourceId, tree);
        }

        const int numCols = tree.topology().getNumColsOfImage();
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

    void onNodeRemoved(NodeId nodeId, const WeightedMorphologicalTree&) const override {
        resetLocalSummary(nodeId);
        resetSubtreeSummary(nodeId);
    }
};

} // namespace mmcfilters::adjust
