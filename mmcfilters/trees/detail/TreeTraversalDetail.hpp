#pragma once

#include "../../utils/Common.hpp"
#include "../MorphologicalTree.hpp"
#include "CommittedTreeAccess.hpp"

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
 *
 * @param tree Tree topology used by the operation.
 * @param rootNodeId Identifier of the traversal root.
 * @param preProcessing Callback invoked before visiting a node children.
 * @param mergeProcessing Callback invoked after completing one child.
 * @param postProcessing Callback invoked after all children are merged.
 */
template <class TreeLike, class PreProcessing, class MergeProcessing, class PostProcessing>
inline void traversePostOrder(TreeLike& tree, NodeId rootNodeId, PreProcessing&& preProcessing, MergeProcessing&& mergeProcessing,
                              PostProcessing&& postProcessing) {
    preProcessing(rootNodeId);
    for (NodeId childNodeId : tree.getChildren(rootNodeId)) {
        traversePostOrder(tree, childNodeId, preProcessing, mergeProcessing, postProcessing);
        mergeProcessing(rootNodeId, childNodeId);
    }
    postProcessing(rootNodeId);
}

/**
 * @brief Post-order traversal for kernels whose tree domain is already established.
 *
 * Unlike the public/generic traversal, this path does not revalidate each
 * child access or iterator version.  A boundary caller must establish a
 * committed topology and a live root before entering it.
 */
namespace kernel {

/**
 * @brief Runs validation-free post-order traversal over an established tree.
 * @param tree Established tree topology.
 * @param rootNodeId Established live traversal root.
 * @param preProcessing Callback invoked before visiting a node's children.
 * @param mergeProcessing Callback invoked after completing one child.
 * @param postProcessing Callback invoked after all children are merged.
 */
template <class PreProcessing, class MergeProcessing, class PostProcessing>
inline void traversePostOrder(const MorphologicalTree& tree, NodeId rootNodeId, PreProcessing&& preProcessing, MergeProcessing&& mergeProcessing,
                              PostProcessing&& postProcessing) {
    preProcessing(rootNodeId);
    for (NodeId childNodeId : CommittedTreeAccess::children(tree, rootNodeId)) {
        traversePostOrder(tree, childNodeId, preProcessing, mergeProcessing, postProcessing);
        mergeProcessing(rootNodeId, childNodeId);
    }
    postProcessing(rootNodeId);
}

} // namespace kernel

} // namespace mmcfilters::detail
