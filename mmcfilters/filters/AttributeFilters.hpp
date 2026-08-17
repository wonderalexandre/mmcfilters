#pragma once

#include "../utils/Image.hpp"
#include "../utils/Common.hpp"
#include "../trees/TreeAltitudeAlgorithms.hpp"
#include "../trees/ValuedMorphologicalTree.hpp"
#include "../trees/ValuedMorphologicalTreeView.hpp"
#include "../trees/detail/CommittedTreeAccess.hpp"
#include "../trees/detail/TreeTraversalDetail.hpp"
#include "../utils/CommittedImageAccess.hpp"
#include "../utils/Contract.hpp"
#include "AttributeReconstructionFilters.hpp"
#include "detail/ViterbiDecision.hpp"

#include <cmath>
#include <concepts>
#include <memory>
#include <stack>
#include <stdexcept>
#include <string>
#include <vector>

namespace mmcfilters {

/**
 * @brief Family of attribute-based image filtering operators on morphological trees.
 *
 * `AttributeFilters` groups pruning, Viterbi, and extinction-related
 * operations. Direct and subtractive reconstruction are exposed separately by
 * `DirectAttributeFilter`, `SubtractiveAttributeFilter`, and
 * `SoftSubtractiveAttributeFilter`. The operators consume preservation masks or
 * attribute buffers defined on dense node ids and reconstruct proper-part
 * images as their output.
 *
 * Contract and node-domain assumptions:
 *
 * - every preservation mask, score, and attribute buffer is indexed by the dense
 *   internal `NodeId` slot domain, not by exported Higra ids;
 * - only alive nodes are reconstructed, but buffers must still cover every
 *   internal slot so callers can reuse the tree-wide attribute layout;
 * - output pixels are written through direct smallest-node mapping. A node
 *   represents the support of its full subtree, while `properPart(node)`
 *   contains only the pixels owned directly by that node;
 * - propagation rules are evaluated in explicit root-to-leaf order because the
 *   dense `NodeId` slot order is not a topological order after Higra import or
 *   structural edits;
 * - public object methods snapshot the topology mutation version at
 *   construction and reject use after the tree structure changes.
 *
 * Reconstruction rule summary:
 *
 * - Viterbi rule: a dynamic program chooses a connected preserved set from
 *   threshold-derived preserve/remove costs, then reconstructs it with the
 *   direct rule;
 * - pruning-min rule: rejected branches are cut immediately and their whole
 *   subtree is painted at the current accepted ancestor level;
 * - pruning-max rule: only subtrees whose nodes are all rejected are cut, and
 *   such subtrees are painted at the rejected subtree root level.
 */
template <AltitudeValue T> class AttributeFilters {
  protected:
    /// @cond INTERNAL
    /** @brief Defines the valuedTree altitude-view type used by the filter. */
    using AltitudeView = ValuedMorphologicalTreeView<T>;

    /** @brief Non-owning valued-tree view supplied at construction. */
    AltitudeView view_;
    /** @brief References the valued-tree owner when one was supplied. */
    const ValuedMorphologicalTree<T>* valuedTree_ = nullptr;
    /** @brief References the tree topology processed by the filter. */
    const MorphologicalTree& tree;
    /** @brief Topology mutation version captured at construction used to detect stale derived state. */
    std::size_t treeMutationVersion_ = 0;

    /**
     * @brief Returns the active valued-tree view.
     *
     * @return The active valued-tree view.
     */
    AltitudeView view() const noexcept { return view_; }

    /**
     * @brief Validates stable tree.
     *
     * @param context Operation name used in diagnostics.
     */
    void requireStableTree(const char* context) const { MMCFILTERS_CONTRACT_CHECKED_ONLY(tree.requireMutationVersion(treeMutationVersion_, context)); }

    /**
     * @brief Validates attribute pointer.
     *
     * @param attribute Attribute requested by the operation.
     * @param context Operation name used in diagnostics.
     */
    template <std::floating_point Real> static void requireAttributePointer(const Real* attribute, const char* context) {
        MMCFILTERS_CONTRACT_REQUIRE(attribute != nullptr,
                                    throw std::invalid_argument(std::string(context) + " requires a non-null attribute buffer."));
    }

    /**
     * @brief Validates a node-preservation mask against the dense node domain.
     *
     * @param tree Tree topology.
     * @param nodePreservationMask Per-node decisions where true preserves a node.
     * @param context Operation name used in diagnostics.
     */
    static void requireNodePreservationMaskSize(const MorphologicalTree& tree, const NodePreservationMask& nodePreservationMask,
                                                const char* context) {
        MMCFILTERS_CONTRACT_REQUIRE(nodePreservationMask.size() == static_cast<std::size_t>(tree.numInternalNodeSlots()),
                                    throw std::invalid_argument(std::string(context) +
                                                                " node-preservation-mask size must match the internal node slot count."));
    }

    /**
     * @brief Validates output image.
     *
     * @param tree Tree topology.
     * @param image Image.
     * @param context Operation name used in diagnostics.
     */
    template <typename TImagePtr> static void requireOutputImage(const MorphologicalTree& tree, const TImagePtr& image, const char* context) {
        MMCFILTERS_CONTRACT_REQUIRE(image != nullptr, throw std::invalid_argument(std::string(context) + " requires a non-null output image."));
        MMCFILTERS_CONTRACT_REQUIRE(image->getNumRows() == tree.numRows() && image->getNumColumns() == tree.numColumns(),
                                    throw std::invalid_argument(std::string(context) + " output image shape must match the tree image domain."));
    }

    /**
     * @brief Returns the altitude of a tree node.
     *
     * @param view Valued-tree view.
     * @param nodeId Dense internal node identifier.
     * @return Altitude associated with the node.
     */
    static T altitudeOf(const AltitudeView& view, NodeId nodeId) noexcept { return view.nodeAltitudes()[static_cast<std::size_t>(nodeId)]; }

    /**
     * @brief Writes one node's direct proper parts only.
     *
     * Reconstruction is done by writing every alive node's direct proper parts
     * exactly once. Descendant supports are intentionally not touched here; use
     * `writeSubtreeProperParts` only when a pruning rule collapses a whole
     * branch to a single gray level.
     *
     * @param tree Tree topology.
     * @param nodeId Dense internal node identifier.
     * @param output Destination buffer receiving the result.
     * @param value Value.
     */
    template <typename TValue> static void writeProperParts(const MorphologicalTree& tree, NodeId nodeId, TValue* output, TValue value) {
        for (PixelId pixel : detail::CommittedTreeAccess::properParts(tree, nodeId)) {
            output[pixel] = value;
        }
    }

    /**
     * @brief Paints all direct proper parts owned by nodes in one subtree.
     *
     * This is the physical image-domain effect of pruning a branch: every node
     * below the cut, including the cut node itself, receives one replacement
     * altitude chosen by the specific pruning convention.
     *
     * @param tree Tree topology.
     * @param nodeId Dense internal node identifier.
     * @param output Destination buffer receiving the result.
     * @param value Value.
     */
    template <typename TValue> static void writeSubtreeProperParts(const MorphologicalTree& tree, NodeId nodeId, TValue* output, TValue value) {
        for (NodeId subtreeNodeId : detail::CommittedTreeAccess::subtree(tree, nodeId)) {
            writeProperParts(tree, subtreeNodeId, output, value);
        }
    }

    /**
     * @brief Pruning-min reconstruction from an explicit node-preservation mask.
     *
     * The convention used here is ancestor-level pruning: once a child branch is
     * rejected, the whole branch is painted at the current accepted node altitude
     * and no deeper descendants are evaluated. This matches the attribute
     * overload where `nodePreservationMask[node] == (attribute[node] > threshold)`.
     *
     * @param view Tree view.
     * @param nodePreservationMask Per-node decisions where true preserves a node.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    static void filteringByPruningMinMaskImpl(AltitudeView view, const NodePreservationMask& nodePreservationMask, ImagePtr<T> imgOutputPtr) {
        const char* context = "AttributeFilters::filteringByPruningMin";
        view.requireTopologyUnchanged(context);
        const MorphologicalTree& tree = view.topology();
        requireNodePreservationMaskSize(tree, nodePreservationMask, context);
        requireOutputImage(tree, imgOutputPtr, context);
        std::stack<NodeId> stack;
        stack.push(tree.root());
        auto imgOutput = imgOutputPtr->rawData();

        while (!stack.empty()) {
            const NodeId nodeId = stack.top();
            stack.pop();
            writeProperParts(tree, nodeId, imgOutput, altitudeOf(view, nodeId));
            for (NodeId childNodeId : detail::CommittedTreeAccess::children(tree, nodeId)) {
                if (nodePreservationMask[static_cast<std::size_t>(childNodeId)]) {
                    stack.push(childNodeId);
                } else {
                    writeSubtreeProperParts(tree, childNodeId, imgOutput, altitudeOf(view, nodeId));
                }
            }
        }
    }

    /**
     * @brief Pruning-max reconstruction from an explicit node-preservation mask.
     *
     * This rule first computes a bottom-up rejected-subtree marker. A node is
     * collapsible only when the node itself is rejected and every descendant
     * branch is also collapsible. Mixed subtrees remain traversable so accepted
     * descendants can preserve their own levels.
     *
     * @param view Tree view.
     * @param nodePreservationMask Per-node decisions where true preserves a node.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    static void filteringByPruningMaxMaskImpl(AltitudeView view, const NodePreservationMask& nodePreservationMask, ImagePtr<T> imgOutputPtr) {
        const char* context = "AttributeFilters::filteringByPruningMax";
        view.requireTopologyUnchanged(context);
        const MorphologicalTree& tree = view.topology();
        requireNodePreservationMaskSize(tree, nodePreservationMask, context);
        requireOutputImage(tree, imgOutputPtr, context);
        // This internal marker means "this whole subtree can be collapsed",
        // which is the opposite of the caller's preservation decision at the leaves
        // before descendant information is merged.
        std::vector<uint8_t> collapsibleRejectedSubtree(tree.numInternalNodeSlots(), false);
        detail::traversePostOrder(
            tree, tree.root(),
            [&collapsibleRejectedSubtree, &nodePreservationMask](NodeId nodeId) -> void {
                collapsibleRejectedSubtree[nodeId] = !nodePreservationMask[static_cast<std::size_t>(nodeId)];
            },
            [&collapsibleRejectedSubtree](NodeId parentNodeId, NodeId childNodeId) -> void {
                collapsibleRejectedSubtree[parentNodeId] =
                    (collapsibleRejectedSubtree[parentNodeId] & collapsibleRejectedSubtree[childNodeId]);
            },
            [](NodeId) -> void {});

        auto imgOutput = imgOutputPtr->rawData();
        std::stack<NodeId> stack;
        stack.push(tree.root());
        while (!stack.empty()) {
            const NodeId nodeId = stack.top();
            stack.pop();
            writeProperParts(tree, nodeId, imgOutput, altitudeOf(view, nodeId));
            for (NodeId childNodeId : detail::CommittedTreeAccess::children(tree, nodeId)) {
                if (!collapsibleRejectedSubtree[childNodeId]) {
                    stack.push(childNodeId);
                } else {
                    writeSubtreeProperParts(tree, childNodeId, imgOutput, altitudeOf(view, childNodeId));
                }
            }
        }
    }

    /**
     * @brief Pruning-min reconstruction from a node attribute threshold.
     *
     * Nodes with `attribute > threshold` stay traversable. Rejected child
     * branches are reconstructed at the current node altitude, which makes this
     * overload equivalent to passing `attribute[node] > threshold` to the
     * node-preservation-mask overload.
     *
     * @param view Tree view.
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    template <std::floating_point Real>
    static void filteringByPruningMinAttributeImpl(AltitudeView view, const Real* attribute, Real threshold, ImagePtr<T> imgOutputPtr) {
        const char* context = "AttributeFilters::filteringByPruningMin";
        view.requireTopologyUnchanged(context);
        const MorphologicalTree& tree = view.topology();
        requireAttributePointer(attribute, context);
        requireOutputImage(tree, imgOutputPtr, context);
        auto imgOutput = imgOutputPtr->rawData();
        std::stack<NodeId> stack;
        stack.push(tree.root());
        while (!stack.empty()) {
            const NodeId nodeId = stack.top();
            stack.pop();
            writeProperParts(tree, nodeId, imgOutput, altitudeOf(view, nodeId));
            for (NodeId childNodeId : detail::CommittedTreeAccess::children(tree, nodeId)) {
                if (attribute[childNodeId] > threshold) {
                    stack.push(childNodeId);
                } else {
                    writeSubtreeProperParts(tree, childNodeId, imgOutput, altitudeOf(view, nodeId));
                }
            }
        }
    }

    /**
     * @brief Pruning-max reconstruction from a node attribute threshold.
     *
     * The bottom-up marker is true only for fully rejected subtrees. When such a
     * branch is found during the top-down reconstruction pass, it is painted at
     * the rejected branch root altitude rather than at the accepted ancestor.
     *
     * @param view Tree view.
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    template <std::floating_point Real>
    static void filteringByPruningMaxAttributeImpl(AltitudeView view, const Real* attribute, Real threshold, ImagePtr<T> imgOutputPtr) {
        const char* context = "AttributeFilters::filteringByPruningMax";
        view.requireTopologyUnchanged(context);
        const MorphologicalTree& tree = view.topology();
        requireAttributePointer(attribute, context);
        requireOutputImage(tree, imgOutputPtr, context);
        // This internal marker records collapsible rejected subtrees, not a
        // direct keep decision. The post-order merge implements the universal
        // quantifier over descendants.
        std::vector<uint8_t> collapsibleRejectedSubtree(tree.numInternalNodeSlots(), false);
        detail::traversePostOrder(
            tree, tree.root(),
            [&collapsibleRejectedSubtree, attribute, threshold](NodeId nodeId) -> void {
                if (attribute[nodeId] <= threshold) {
                    collapsibleRejectedSubtree[nodeId] = true;
                }
            },
            [&collapsibleRejectedSubtree](NodeId parentNodeId, NodeId childNodeId) -> void {
                collapsibleRejectedSubtree[parentNodeId] =
                    (collapsibleRejectedSubtree[parentNodeId] & collapsibleRejectedSubtree[childNodeId]);
            },
            [](NodeId) -> void {});

        auto imgOutput = imgOutputPtr->rawData();
        std::stack<NodeId> stack;
        stack.push(tree.root());
        while (!stack.empty()) {
            const NodeId nodeId = stack.top();
            stack.pop();
            writeProperParts(tree, nodeId, imgOutput, altitudeOf(view, nodeId));
            for (NodeId childNodeId : detail::CommittedTreeAccess::children(tree, nodeId)) {
                if (!collapsibleRejectedSubtree[childNodeId]) {
                    stack.push(childNodeId);
                } else {
                    writeSubtreeProperParts(tree, childNodeId, imgOutput, altitudeOf(view, childNodeId));
                }
            }
        }
    }

    /// @endcond

  public:
    /**
     * @brief Creates filtering operators over a non-owning valued tree view.
     *
     * The view must remain valid for the lifetime of the filter object. Public
     * methods reject use after the underlying topology mutates.
     *
     * @param view Tree view.
     */
    explicit AttributeFilters(AltitudeView view) : view_{view}, tree{view_.topology()}, treeMutationVersion_{tree.getMutationVersion()} {
        MMCFILTERS_CONTRACT_CHECKED_ONLY(view_.requireTopologyUnchanged("AttributeFilters"));
    }

    /**
     * @brief Creates filtering operators over an owned valued tree.
     *
     * This overload keeps a pointer to the owner so methods that require
     * tree-owned altitude state, such as MSER-assisted preservation selection, can be
     * used.
     *
     * @param valuedTree Valued tree.
     */
    explicit AttributeFilters(const ValuedMorphologicalTree<T>& valuedTree) : AttributeFilters(valuedTree.asView()) { valuedTree_ = &valuedTree; }

    /**
     * @brief Destroys the attribute-filter facade.
     */
    ~AttributeFilters() = default;

    /**
     * @brief Applies pruning-min filtering from an attribute buffer.
     *
     * Nodes with attribute values above `threshold` remain traversable; rejected
     * subtrees are reconstructed at the ancestor level selected by the pruning-min
     * rule.
     *
     * @param attr Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @return Image produced by the operation.
     */
    template <std::floating_point Real> [[nodiscard]] ImagePtr<T> filteringByPruningMin(const std::shared_ptr<Real[]>& attr, Real threshold) {
        return filteringByPruningMin(attr.get(), threshold);
    }

    /**
     * @brief Applies pruning-min filtering from a raw internal-node attribute buffer.
     *
     * @param attr Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @return Image produced by the operation.
     */
    template <std::floating_point Real> [[nodiscard]] ImagePtr<T> filteringByPruningMin(const Real* attr, Real threshold) {
        ImagePtr<T> imgOutput = Image<T>::create(this->tree.numRows(), this->tree.numColumns());
        filteringByPruningMinAttributeImpl(view(), attr, threshold, imgOutput);
        return imgOutput;
    }

    /**
     * @brief Applies pruning-max filtering from an attribute buffer.
     *
     * Nodes with attribute values above `threshold` are kept, and fully rejected
     * subtrees are reconstructed at their own subtree levels.
     *
     * @param attr Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @return Image produced by the operation.
     */
    template <std::floating_point Real> [[nodiscard]] ImagePtr<T> filteringByPruningMax(const std::shared_ptr<Real[]>& attr, Real threshold) {
        return filteringByPruningMax(attr.get(), threshold);
    }

    /**
     * @brief Applies pruning-max filtering from a raw internal-node attribute buffer.
     *
     * @param attr Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @return Image produced by the operation.
     */
    template <std::floating_point Real> [[nodiscard]] ImagePtr<T> filteringByPruningMax(const Real* attr, Real threshold) {
        ImagePtr<T> imgOutput = Image<T>::create(this->tree.numRows(), this->tree.numColumns());
        filteringByPruningMaxAttributeImpl(view(), attr, threshold, imgOutput);
        return imgOutput;
    }

    /**
     * @brief Applies Salembier-style Viterbi filtering from a raw attribute buffer.
     *
     * The internal `detail` implementation turns `attr` and `threshold` into
     * preserve/remove costs and solves the optimal connected keep mask on the
     * tree. The root is always preserved. Once the Viterbi path removes a node,
     * every descendant is removed as well; preserved descendants can only appear
     * below preserved ancestors. The final image is reconstructed by the direct
     * rule using that connected node-preservation mask.
     *
     * @param attr Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @return Image produced by the operation.
     */
    template <std::floating_point Real> [[nodiscard]] ImagePtr<T> filteringByViterbiRule(const Real* attr, Real threshold) {
        requireStableTree("AttributeFilters::filteringByViterbiRule");
        auto costs = detail::makeThresholdViterbiCosts(tree, attr, threshold);
        NodePreservationMask nodePreservationMask(detail::computeViterbiPreservationDecisions(tree, costs));
        return applyDirectAttributeFilter(view(), nodePreservationMask);
    }

    /**
     * @brief Applies pruning-min filtering from a dense node-preservation mask.
     *
     * @param nodePreservationMask Per-node decisions where true preserves a node.
     * @return Image produced by the operation.
     */
    [[nodiscard]] ImagePtr<T> filteringByPruningMin(const NodePreservationMask& nodePreservationMask) {
        ImagePtr<T> imgOutput = Image<T>::create(this->tree.numRows(), this->tree.numColumns());
        filteringByPruningMinMaskImpl(view(), nodePreservationMask, imgOutput);
        return imgOutput;
    }

    /**
     * @brief Applies pruning-max filtering from a dense node-preservation mask.
     *
     * @param nodePreservationMask Per-node decisions where true preserves a node.
     * @return Image produced by the operation.
     */
    [[nodiscard]] ImagePtr<T> filteringByPruningMax(const NodePreservationMask& nodePreservationMask) {
        ImagePtr<T> imgOutput = Image<T>::create(this->tree.numRows(), this->tree.numColumns());
        filteringByPruningMaxMaskImpl(view(), nodePreservationMask, imgOutput);
        return imgOutput;
    }

    /**
     * @brief Writes pruning-min filtering from a node-preservation mask into an output image.
     *
     * @param valuedTree Valued tree.
     * @param nodePreservationMask Per-node decisions where true preserves a node.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    static void filteringByPruningMin(const ValuedMorphologicalTreeView<T>& valuedTree, const NodePreservationMask& nodePreservationMask,
                                      ImagePtr<T> imgOutputPtr) {
        filteringByPruningMinMaskImpl(valuedTree, nodePreservationMask, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-min filtering from a valued-tree owner into an output image.
     *
     * @param valuedTree Valued tree.
     * @param nodePreservationMask Per-node decisions where true preserves a node.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    static void filteringByPruningMin(const ValuedMorphologicalTree<T>& valuedTree, const NodePreservationMask& nodePreservationMask,
                                      ImagePtr<T> imgOutputPtr) {
        filteringByPruningMin(valuedTree.asView(), nodePreservationMask, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-max filtering from a node-preservation mask into an output image.
     *
     * @param valuedTree Valued tree.
     * @param nodePreservationMask Per-node decisions where true preserves a node.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    static void filteringByPruningMax(const ValuedMorphologicalTreeView<T>& valuedTree, const NodePreservationMask& nodePreservationMask,
                                      ImagePtr<T> imgOutputPtr) {
        filteringByPruningMaxMaskImpl(valuedTree, nodePreservationMask, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-max filtering from a valued-tree owner into an output image.
     *
     * @param valuedTree Valued tree.
     * @param nodePreservationMask Per-node decisions where true preserves a node.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    static void filteringByPruningMax(const ValuedMorphologicalTree<T>& valuedTree, const NodePreservationMask& nodePreservationMask,
                                      ImagePtr<T> imgOutputPtr) {
        filteringByPruningMax(valuedTree.asView(), nodePreservationMask, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-min filtering from an owned attribute buffer into an output image.
     *
     * @param valuedTree Valued tree.
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    template <std::floating_point Real>
    static void filteringByPruningMin(const ValuedMorphologicalTreeView<T>& valuedTree, const std::shared_ptr<Real[]>& attribute, Real threshold, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMin(valuedTree, attribute.get(), threshold, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-min filtering from a valued-tree owner and owned attribute buffer.
     *
     * @param valuedTree Valued tree.
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    template <std::floating_point Real>
    static void filteringByPruningMin(const ValuedMorphologicalTree<T>& valuedTree, const std::shared_ptr<Real[]>& attribute, Real threshold,
                                      ImagePtr<T> imgOutputPtr) {
        filteringByPruningMin(valuedTree.asView(), attribute.get(), threshold, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-min filtering from a raw attribute buffer into an output image.
     *
     * @param valuedTree Valued tree.
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    template <std::floating_point Real>
    static void filteringByPruningMin(const ValuedMorphologicalTreeView<T>& valuedTree, const Real* attribute, Real threshold, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMinAttributeImpl(valuedTree, attribute, threshold, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-min filtering from a valued-tree owner and raw attribute buffer.
     *
     * @param valuedTree Valued tree.
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    template <std::floating_point Real>
    static void filteringByPruningMin(const ValuedMorphologicalTree<T>& valuedTree, const Real* attribute, Real threshold, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMin(valuedTree.asView(), attribute, threshold, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-max filtering from an owned attribute buffer into an output image.
     *
     * @param valuedTree Valued tree.
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    template <std::floating_point Real>
    static void filteringByPruningMax(const ValuedMorphologicalTreeView<T>& valuedTree, const std::shared_ptr<Real[]>& attribute, Real threshold, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMax(valuedTree, attribute.get(), threshold, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-max filtering from a valued-tree owner and owned attribute buffer.
     *
     * @param valuedTree Valued tree.
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    template <std::floating_point Real>
    static void filteringByPruningMax(const ValuedMorphologicalTree<T>& valuedTree, const std::shared_ptr<Real[]>& attribute, Real threshold,
                                      ImagePtr<T> imgOutputPtr) {
        filteringByPruningMax(valuedTree.asView(), attribute.get(), threshold, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-max filtering from a raw attribute buffer into an output image.
     *
     * @param valuedTree Valued tree.
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    template <std::floating_point Real>
    static void filteringByPruningMax(const ValuedMorphologicalTreeView<T>& valuedTree, const Real* attribute, Real threshold, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMaxAttributeImpl(valuedTree, attribute, threshold, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-max filtering from a valued-tree owner and raw attribute buffer.
     *
     * @param valuedTree Valued tree.
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    template <std::floating_point Real>
    static void filteringByPruningMax(const ValuedMorphologicalTree<T>& valuedTree, const Real* attribute, Real threshold, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMax(valuedTree.asView(), attribute, threshold, imgOutputPtr);
    }

};

} // namespace mmcfilters
