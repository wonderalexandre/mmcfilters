#pragma once

#include "../../utils/Common.hpp"

namespace mmcfilters::detail {

/**
 * @brief Runs a depth-first post-order traversal over a tree-like topology.
 *
 * @details
 * `TreeLike` only needs to expose `getChildren(NodeId)`. The helper is kept in
 * the tree detail layer because attributes, filters, and contour materializers
 * all need the same child-before-parent aggregation pattern. Callers provide
 * three hooks: one before visiting children, one after each child is completed,
 * and one after all children have been merged into the current node.
 */
template<class TreeLike, class PreProcessing, class MergeProcessing, class PostProcessing>
inline void traversePostOrder(
    TreeLike& tree,
    NodeId rootNodeId,
    PreProcessing&& preProcessing,
    MergeProcessing&& mergeProcessing,
    PostProcessing&& postProcessing)
{
    preProcessing(rootNodeId);
    for (NodeId childNodeId : tree.getChildren(rootNodeId)) {
        traversePostOrder(
            tree,
            childNodeId,
            preProcessing,
            mergeProcessing,
            postProcessing);
        mergeProcessing(rootNodeId, childNodeId);
    }
    postProcessing(rootNodeId);
}

} // namespace mmcfilters::detail
