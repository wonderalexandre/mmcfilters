#pragma once

#include "../MorphologicalTree.hpp"
#include "CommittedTreeAccess.hpp"
#include "../../utils/Common.hpp"
#include "../../utils/Contract.hpp"

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
namespace kernel {

/**
 * @brief Finds the entry node for two established proper-part samples.
 * @param tree Established tree topology and ownership domain.
 * @param anchorProperPart Established anchor proper-part id.
 * @param sampleProperPart Established sample proper-part id.
 * @return First node whose support contains both samples.
 */
inline NodeId properPartEntryNode(const MorphologicalTree& tree, int anchorProperPart, int sampleProperPart) {
    const NodeId anchorOwner = CommittedTreeAccess::properPartOwner(tree, anchorProperPart);
    const NodeId sampleOwner = CommittedTreeAccess::properPartOwner(tree, sampleProperPart);

    if (CommittedTreeAccess::isAncestor(tree, anchorOwner, sampleOwner)) {
        return anchorOwner;
    }
    if (CommittedTreeAccess::isAncestor(tree, sampleOwner, anchorOwner)) {
        return sampleOwner;
    }
    return CommittedTreeAccess::lowestCommonAncestor(tree, anchorOwner, sampleOwner);
}

} // namespace kernel

inline NodeId properPartEntryNode(const MorphologicalTree& tree, int anchorProperPart, int sampleProperPart) {
    MMCFILTERS_CONTRACT_CHECKED_ONLY(if (anchorProperPart < 0 || anchorProperPart >= tree.getNumTotalProperParts() || sampleProperPart < 0 ||
                                         sampleProperPart >= tree.getNumTotalProperParts()) { return InvalidNode; });
    return kernel::properPartEntryNode(tree, anchorProperPart, sampleProperPart);
}

} // namespace mmcfilters::detail
