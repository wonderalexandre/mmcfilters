#pragma once

#include "../MorphologicalTree.hpp"
#include "../../utils/Common.hpp"

namespace mmcfilters::detail {

/**
 * @brief Returns the first ancestor where two proper-part samples share support.
 *
 * @details
 * The entry node is the structural point where `sampleProperPart` becomes
 * visible while expanding the support rooted at `anchorProperPart` upward:
 *
 * - if the anchor owner already contains the sample owner, the entry is the
 *   anchor owner;
 * - if the sample owner is an ancestor of the anchor owner, the entry is the
 *   sample owner;
 * - otherwise the two owners are incomparable and the event enters at their
 *   lowest common ancestor.
 *
 * Invalid proper parts or proper parts without an owner return `InvalidNode`.
 *
 * @param tree Tree topology used by the operation.
 * @param anchorProperPart Proper-part data represented by `anchorProperPart`.
 * @param sampleProperPart Proper-part data represented by `sampleProperPart`.
 * @return The first ancestor where two proper-part samples share support.
 */
inline NodeId properPartEntryNode(const MorphologicalTree& tree, int anchorProperPart, int sampleProperPart) {
    if (anchorProperPart < 0 || anchorProperPart >= tree.getNumTotalProperParts() || sampleProperPart < 0 ||
        sampleProperPart >= tree.getNumTotalProperParts()) {
        return InvalidNode;
    }

    const NodeId anchorOwner = tree.getProperPartOwner(anchorProperPart);
    const NodeId sampleOwner = tree.getProperPartOwner(sampleProperPart);
    if (anchorOwner == InvalidNode || sampleOwner == InvalidNode) {
        return InvalidNode;
    }

    if (tree.isAncestor(anchorOwner, sampleOwner)) {
        return anchorOwner;
    }
    if (tree.isAncestor(sampleOwner, anchorOwner)) {
        return sampleOwner;
    }
    return tree.getLowestCommonAncestor(anchorOwner, sampleOwner);
}

} // namespace mmcfilters::detail
