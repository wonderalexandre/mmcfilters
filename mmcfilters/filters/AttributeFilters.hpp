#pragma once

#include "../utils/Image.hpp"
#include "../utils/Common.hpp"
#include "../trees/TreeAltitudeAlgorithms.hpp"
#include "../trees/WeightedMorphologicalTree.hpp"
#include "../trees/WeightedTreeView.hpp"
#include "../trees/detail/CommittedTreeAccess.hpp"
#include "../trees/detail/TreeTraversalDetail.hpp"
#include "../utils/CommittedImageAccess.hpp"
#include "../utils/Contract.hpp"
#include "DepthStableRegionComputer.hpp"
#include "MSERComputer.hpp"
#include "detail/VariationMeasure.hpp"
#include "detail/ViterbiDecision.hpp"

#include <cmath>
#include <concepts>
#include <memory>
#include <stack>
#include <stdexcept>
#include <string>
#include <vector>

namespace mmcfilters {

namespace detail::kernel {

/**
 * @brief Direct reconstruction after topology, altitude, criterion, and output domains were established.
 * @param view Established weighted-tree view.
 * @param criterion Established keep/remove decision indexed by node id.
 * @param output Established output image with the tree grid dimensions.
 */
template <AltitudeValue T>
inline void filterDirect(WeightedTreeView<T> view, const std::vector<bool>& criterion, ImagePtr<T> output) {
    const MorphologicalTree& tree = view.topology();
    const std::span<const T> altitude = view.altitude();
    std::vector<T> mappedLevel(static_cast<std::size_t>(tree.getNumInternalNodeSlots()));
    const NodeId root = tree.getRoot();
    mappedLevel[static_cast<std::size_t>(root)] = altitude[static_cast<std::size_t>(root)];

    std::stack<NodeId> stack;
    stack.push(root);
    while (!stack.empty()) {
        const NodeId node = stack.top();
        stack.pop();
        for (NodeId child : CommittedTreeAccess::children(tree, node)) {
            mappedLevel[static_cast<std::size_t>(child)] =
                criterion[static_cast<std::size_t>(child)] ? altitude[static_cast<std::size_t>(child)] : mappedLevel[static_cast<std::size_t>(node)];
            stack.push(child);
        }
    }

    T* pixels = output->rawData();
    for (NodeId node = 0; node < tree.getNumInternalNodeSlots(); ++node) {
        if (CommittedTreeAccess::isAlive(tree, node)) {
            for (int pixel : CommittedTreeAccess::properParts(tree, node)) {
                pixels[pixel] = mappedLevel[static_cast<std::size_t>(node)];
            }
        }
    }
}

} // namespace detail::kernel

/**
 * @brief Family of attribute-based image filtering operators on morphological trees.
 *
 * `AttributeFilters` groups the classical reconstruction rules used in this
 * project: direct filtering, subtractive filtering, pruning-based filtering,
 * and a few residual/score-based variants. The operators consume criteria or
 * attribute buffers defined on dense node ids and reconstruct proper-part
 * images as their output.
 *
 * Contract and node-domain assumptions:
 *
 * - every criterion, score, and attribute buffer is indexed by the dense
 *   internal `NodeId` slot domain, not by exported Higra ids;
 * - only alive nodes are reconstructed, but buffers must still cover every
 *   internal slot so callers can reuse the tree-wide attribute layout;
 * - output pixels are written through direct proper-part ownership. A node
 *   represents the support of its full subtree, while `getProperParts(node)`
 *   contains only the pixels owned directly by that node;
 * - propagation rules are evaluated in explicit root-to-leaf order because the
 *   dense `NodeId` slot order is not a topological order after Higra import or
 *   structural edits;
 * - public object methods snapshot the topology mutation version at
 *   construction and reject use after the tree structure changes.
 *
 * Reconstruction rule summary:
 *
 * - direct rule: a kept node maps to its own altitude; a rejected node inherits
 *   the already-filtered level of its parent;
 * - subtractive rule: a kept node adds its altitude residue to the filtered
 *   parent level; a rejected node inherits that parent level;
 * - subtractive score rule: each residue is multiplied by a per-node score and
 *   accumulated in floating point;
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
    /** @brief Defines the weighted altitude-view type used by the filter. */
    using AltitudeView = WeightedTreeView<T>;

    /** @brief Stores the non-owning weighted-tree view supplied at construction. */
    AltitudeView view_;
    /** @brief References the weighted-tree owner when one was supplied. */
    const WeightedMorphologicalTree<T>* weighted_ = nullptr;
    /** @brief References the tree topology processed by the filter. */
    const MorphologicalTree& tree;
    /** @brief Stores the topology mutation version captured at construction. */
    std::size_t treeMutationVersion_ = 0;

    /**
     * @brief Returns the active weighted-tree view.
     *
     * @return The active weighted-tree view.
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
     * @brief Validates criterion size.
     *
     * @param tree Tree topology used by the operation.
     * @param criterion Per-node selection criterion.
     * @param context Operation name used in diagnostics.
     */
    static void requireCriterionSize(const MorphologicalTree& tree, const std::vector<bool>& criterion, const char* context) {
        MMCFILTERS_CONTRACT_REQUIRE(criterion.size() == static_cast<std::size_t>(tree.getNumInternalNodeSlots()),
                                    throw std::invalid_argument(std::string(context) + " criterion size must match the internal node slot count."));
    }

    /**
     * @brief Validates score size.
     *
     * @param tree Tree topology used by the operation.
     * @param scores Per-node scores used by the operation.
     * @param context Operation name used in diagnostics.
     */
    static void requireScoreSize(const MorphologicalTree& tree, const std::vector<float>& scores, const char* context) {
        MMCFILTERS_CONTRACT_REQUIRE(scores.size() == static_cast<std::size_t>(tree.getNumInternalNodeSlots()),
                                    throw std::invalid_argument(std::string(context) + " score size must match the internal node slot count."));
    }

    /**
     * @brief Validates output image.
     *
     * @param tree Tree topology used by the operation.
     * @param image Image used by the operation.
     * @param context Operation name used in diagnostics.
     */
    template <typename TImagePtr> static void requireOutputImage(const MorphologicalTree& tree, const TImagePtr& image, const char* context) {
        MMCFILTERS_CONTRACT_REQUIRE(image != nullptr, throw std::invalid_argument(std::string(context) + " requires a non-null output image."));
        MMCFILTERS_CONTRACT_REQUIRE(image->getNumRows() == tree.getNumRowsOfGridDomain2D() && image->getNumCols() == tree.getNumColsOfGridDomain2D(),
                                    throw std::invalid_argument(std::string(context) + " output image shape must match the tree image domain."));
    }

    /**
     * @brief Returns the altitude of a tree node.
     *
     * @param view Weighted-tree view used by the operation.
     * @param nodeId Identifier of the node used by the operation.
     * @return Altitude associated with the node.
     */
    static T altitudeOf(const AltitudeView& view, NodeId nodeId) noexcept { return view.altitude()[static_cast<std::size_t>(nodeId)]; }

    /**
     * @brief Computes the residue between the original and filtered values.
     *
     * @param view Weighted-tree view used by the operation.
     * @param nodeId Identifier of the node used by the operation.
     * @return Altitude residue associated with the node.
     */
    static AltitudeDiff<T> residueOf(const AltitudeView& view, NodeId nodeId) noexcept {
        const MorphologicalTree& tree = view.topology();
        const std::span<const T> altitude = view.altitude();
        const NodeId parentNodeId = detail::CommittedTreeAccess::nodeParent(tree, nodeId);
        if (parentNodeId == InvalidNode || parentNodeId == nodeId) {
            return static_cast<AltitudeDiff<T>>(altitude[static_cast<std::size_t>(nodeId)]);
        }
        return static_cast<AltitudeDiff<T>>(altitude[static_cast<std::size_t>(nodeId)]) -
               static_cast<AltitudeDiff<T>>(altitude[static_cast<std::size_t>(parentNodeId)]);
    }

    /**
     * @brief Writes one node's direct proper parts only.
     *
     * Reconstruction is done by writing every alive node's direct proper parts
     * exactly once. Descendant supports are intentionally not touched here; use
     * `writeSubtreeProperParts` only when a pruning rule collapses a whole
     * branch to a single gray level.
     *
     * @param tree Tree topology used by the operation.
     * @param nodeId Identifier of the node used by the operation.
     * @param output Destination buffer receiving the result.
     * @param value Value used by the operation.
     */
    template <typename TValue> static void writeProperParts(const MorphologicalTree& tree, NodeId nodeId, TValue* output, TValue value) {
        for (int pixel : detail::CommittedTreeAccess::properParts(tree, nodeId)) {
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
     * @param tree Tree topology used by the operation.
     * @param nodeId Identifier of the node used by the operation.
     * @param output Destination buffer receiving the result.
     * @param value Value used by the operation.
     */
    template <typename TValue> static void writeSubtreeProperParts(const MorphologicalTree& tree, NodeId nodeId, TValue* output, TValue value) {
        for (NodeId subtreeNodeId : detail::CommittedTreeAccess::subtree(tree, nodeId)) {
            writeProperParts(tree, subtreeNodeId, output, value);
        }
    }

    /**
     * @brief Accumulates scored altitude residues from the root to each node.
     *
     * With all scores equal to one, the accumulated level is the original
     * altitude for max-trees and min-trees. Scores in `[0, 1]` behave like a
     * soft subtractive filter, but the implementation accepts any finite float
     * value supplied by the caller.
     *
     * @param view Tree view used by the operation.
     * @param prob Per-node probability values.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    static void filteringBySubtractiveScoreRuleImpl(AltitudeView view, std::vector<float>& prob, ImageFloatPtr imgOutputPtr) {
        view.requireTopologyUnchanged("AttributeFilters::filteringBySubtractiveScoreRule");
        const MorphologicalTree& tree = view.topology();
        requireScoreSize(tree, prob, "AttributeFilters::filteringBySubtractiveScoreRule");
        requireOutputImage(tree, imgOutputPtr, "AttributeFilters::filteringBySubtractiveScoreRule");
        std::unique_ptr<float[]> mapLevel(new float[tree.getNumInternalNodeSlots()]);

        const NodeId rootNodeId = tree.getRoot();
        mapLevel[rootNodeId] = static_cast<float>(altitudeOf(view, rootNodeId));

        // Parent filtered levels must be available before children are visited.
        // `getAliveNodeIds()` is slot-order, not tree-order, so use an explicit
        // top-down traversal that is valid for image-built and Higra-imported trees.
        std::stack<NodeId> stack;
        stack.push(rootNodeId);
        while (!stack.empty()) {
            const NodeId nodeId = stack.top();
            stack.pop();
            for (NodeId childNodeId : detail::CommittedTreeAccess::children(tree, nodeId)) {
                mapLevel[childNodeId] = mapLevel[nodeId] + (static_cast<float>(residueOf(view, childNodeId)) * prob[childNodeId]);
                stack.push(childNodeId);
            }
        }

        auto imgOutput = imgOutputPtr->rawData();
        for (NodeId nodeId : tree.getAliveNodeIds()) {
            writeProperParts(tree, nodeId, imgOutput, mapLevel[nodeId]);
        }
    }

    /**
     * @brief Binary subtractive reconstruction from a dense keep criterion.
     *
     * `mapLevel[node]` stores the filtered altitude assigned to a node after
     * applying the rule along the root-to-node path. Using `AltitudeDiff<T>`
     * avoids unsigned wraparound while residues from min-trees are negative.
     *
     * @param view Tree view used by the operation.
     * @param criterion Per-node selection criterion.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    static void filteringBySubtractiveRuleImpl(AltitudeView view, std::vector<bool>& criterion, ImagePtr<T> imgOutputPtr) {
        const char* context = "AttributeFilters::filteringBySubtractiveRule";
        view.requireTopologyUnchanged(context);
        const MorphologicalTree& tree = view.topology();
        requireCriterionSize(tree, criterion, context);
        requireOutputImage(tree, imgOutputPtr, context);
        std::vector<AltitudeDiff<T>> mapLevel(static_cast<std::size_t>(tree.getNumInternalNodeSlots()));
        const NodeId rootNodeId = tree.getRoot();
        mapLevel[rootNodeId] = static_cast<AltitudeDiff<T>>(altitudeOf(view, rootNodeId));

        // The parent level in `mapLevel` is the already-filtered level, not the
        // original altitude. This makes traversal order part of the algorithm.
        std::stack<NodeId> stack;
        stack.push(rootNodeId);
        while (!stack.empty()) {
            const NodeId nodeId = stack.top();
            stack.pop();
            for (NodeId childNodeId : detail::CommittedTreeAccess::children(tree, nodeId)) {
                mapLevel[childNodeId] =
                    criterion[childNodeId] ? static_cast<AltitudeDiff<T>>(mapLevel[nodeId] + residueOf(view, childNodeId)) : mapLevel[nodeId];
                stack.push(childNodeId);
            }
        }

        auto imgOutput = imgOutputPtr->rawData();
        for (NodeId nodeId : tree.getAliveNodeIds()) {
            writeProperParts(tree, nodeId, imgOutput, static_cast<T>(mapLevel[nodeId]));
        }
    }

    /**
     * @brief Pruning-min reconstruction from an explicit keep criterion.
     *
     * The convention used here is ancestor-level pruning: once a child branch is
     * rejected, the whole branch is painted at the current accepted node level
     * and no deeper descendants are evaluated. This matches the attribute
     * overload where `criterion[node] == (attribute[node] > threshold)`.
     *
     * @param view Tree view used by the operation.
     * @param criterion Per-node selection criterion.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    static void filteringByPruningMinCriterionImpl(AltitudeView view, std::vector<bool>& criterion, ImagePtr<T> imgOutputPtr) {
        const char* context = "AttributeFilters::filteringByPruningMin";
        view.requireTopologyUnchanged(context);
        const MorphologicalTree& tree = view.topology();
        requireCriterionSize(tree, criterion, context);
        requireOutputImage(tree, imgOutputPtr, context);
        std::stack<NodeId> stack;
        stack.push(tree.getRoot());
        auto imgOutput = imgOutputPtr->rawData();

        while (!stack.empty()) {
            const NodeId nodeId = stack.top();
            stack.pop();
            writeProperParts(tree, nodeId, imgOutput, altitudeOf(view, nodeId));
            for (NodeId childNodeId : detail::CommittedTreeAccess::children(tree, nodeId)) {
                if (criterion[childNodeId]) {
                    stack.push(childNodeId);
                } else {
                    writeSubtreeProperParts(tree, childNodeId, imgOutput, altitudeOf(view, nodeId));
                }
            }
        }
    }

    /**
     * @brief Pruning-max reconstruction from an explicit keep criterion.
     *
     * This rule first computes a bottom-up rejected-subtree marker. A node is
     * collapsible only when the node itself is rejected and every descendant
     * branch is also collapsible. Mixed subtrees remain traversable so accepted
     * descendants can preserve their own levels.
     *
     * @param view Tree view used by the operation.
     * @param keepCriterion Per-node mask of nodes eligible to be kept.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    static void filteringByPruningMaxCriterionImpl(AltitudeView view, std::vector<bool>& keepCriterion, ImagePtr<T> imgOutputPtr) {
        const char* context = "AttributeFilters::filteringByPruningMax";
        view.requireTopologyUnchanged(context);
        const MorphologicalTree& tree = view.topology();
        requireCriterionSize(tree, keepCriterion, context);
        requireOutputImage(tree, imgOutputPtr, context);
        // Internal `criterion` means "this whole subtree can be collapsed",
        // which is the opposite of the caller's keep criterion at the leaves
        // before descendant information is merged.
        std::vector<uint8_t> criterion(tree.getNumInternalNodeSlots(), false);
        detail::traversePostOrder(
            tree, tree.getRoot(), [&criterion, &keepCriterion](NodeId nodeId) -> void { criterion[nodeId] = !keepCriterion[nodeId]; },
            [&criterion](NodeId parentNodeId, NodeId childNodeId) -> void { criterion[parentNodeId] = (criterion[parentNodeId] & criterion[childNodeId]); },
            [](NodeId) -> void {});

        auto imgOutput = imgOutputPtr->rawData();
        std::stack<NodeId> stack;
        stack.push(tree.getRoot());
        while (!stack.empty()) {
            const NodeId nodeId = stack.top();
            stack.pop();
            writeProperParts(tree, nodeId, imgOutput, altitudeOf(view, nodeId));
            for (NodeId childNodeId : detail::CommittedTreeAccess::children(tree, nodeId)) {
                if (!criterion[childNodeId]) {
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
     * criterion overload.
     *
     * @param view Tree view used by the operation.
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
        stack.push(tree.getRoot());
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
     * @param view Tree view used by the operation.
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
        // Internal `criterion` records collapsible rejected subtrees, not a
        // direct keep decision. The post-order merge implements the universal
        // quantifier over descendants.
        std::vector<uint8_t> criterion(tree.getNumInternalNodeSlots(), false);
        detail::traversePostOrder(
            tree, tree.getRoot(),
            [&criterion, attribute, threshold](NodeId nodeId) -> void {
                if (attribute[nodeId] <= threshold) {
                    criterion[nodeId] = true;
                }
            },
            [&criterion](NodeId parentNodeId, NodeId childNodeId) -> void { criterion[parentNodeId] = (criterion[parentNodeId] & criterion[childNodeId]); },
            [](NodeId) -> void {});

        auto imgOutput = imgOutputPtr->rawData();
        std::stack<NodeId> stack;
        stack.push(tree.getRoot());
        while (!stack.empty()) {
            const NodeId nodeId = stack.top();
            stack.pop();
            writeProperParts(tree, nodeId, imgOutput, altitudeOf(view, nodeId));
            for (NodeId childNodeId : detail::CommittedTreeAccess::children(tree, nodeId)) {
                if (!criterion[childNodeId]) {
                    stack.push(childNodeId);
                } else {
                    writeSubtreeProperParts(tree, childNodeId, imgOutput, altitudeOf(view, childNodeId));
                }
            }
        }
    }

    /**
     * @brief Derives an adaptive filtering criterion from attribute stability.
     *
     * @param tree Tree topology used by the operation.
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @param stabilityComputer Evaluator used to measure stability over the selected depth window.
     * @param context Operation name used in diagnostics.
     * @return Values produced by the operation.
     */
    template <std::floating_point Real, class StabilityComputer>
    static std::vector<bool> adaptiveCriterionFromAttributeStability(const MorphologicalTree& tree, const Real* attribute, Real threshold,
                                                                     StabilityComputer& stabilityComputer, const char* context) {
        requireAttributePointer(attribute, context);
        const std::vector<Real>& variation = stabilityComputer.getVariations();
        std::vector<bool> isPruned(tree.getNumInternalNodeSlots(), false);
        for (NodeId nodeId : tree.getAliveNodeIds()) {
            if (attribute[nodeId] < threshold) {
                // Lower finite variation values are more stable. If the center
                // has no complete stability window, keep the
                // historical fallback and prune the rejected node itself.
                if (!detail::isFiniteVariation(variation[static_cast<std::size_t>(nodeId)])) {
                    isPruned[nodeId] = true;
                } else {
                    isPruned[stabilityComputer.nodeWithMinimumVariationInWindow(nodeId)] = true;
                }
            }
        }
        return isPruned;
    }

    /**
     * @brief Derives an adaptive filtering criterion from mask stability.
     *
     * @param tree Tree topology used by the operation.
     * @param criterion Per-node selection criterion.
     * @param stabilityComputer Evaluator used to measure stability over the selected depth window.
     * @param context Operation name used in diagnostics.
     * @return Values produced by the operation.
     */
    template <std::floating_point Real, class StabilityComputer>
    static std::vector<bool> adaptiveCriterionFromMaskStability(const MorphologicalTree& tree, std::vector<bool>& criterion,
                                                                StabilityComputer& stabilityComputer, const char* context) {
        requireCriterionSize(tree, criterion, context);
        const std::vector<Real>& variation = stabilityComputer.getVariations();
        std::vector<bool> isPruned(tree.getNumInternalNodeSlots(), false);
        for (NodeId nodeId : tree.getAliveNodeIds()) {
            if (!criterion[nodeId]) {
                // `criterion == false` marks an attribute rejection. The variation
                // comparison moves the actual pruning decision to the locally
                // smallest finite variation on the ancestor/descendant window.
                if (!detail::isFiniteVariation(variation[static_cast<std::size_t>(nodeId)])) {
                    isPruned[nodeId] = true;
                } else {
                    isPruned[stabilityComputer.nodeWithMinimumVariationInWindow(nodeId)] = true;
                }
            }
        }
        return isPruned;
    }

    /**
     * @brief Returns adaptive criterion impl.
     *
     * @param weighted Weighted tree used by the operation.
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @param delta Delta offset used by the operation.
     * @return Adaptive criterion impl.
     */
    template <std::floating_point Real>
    static std::vector<bool> getAdaptiveCriterionImpl(const WeightedMorphologicalTree<T>& weighted, const Real* attribute, Real threshold,
                                                      AltitudeDiff<T> delta) {
        const MorphologicalTree& tree = weighted.topology();
        MSERComputer<T, Real> mser(weighted);
        (void)mser.computeMSER(delta);
        return adaptiveCriterionFromAttributeStability<Real>(tree, attribute, threshold, mser, "AttributeFilters::getAdaptiveCriterion");
    }

    /**
     * @brief Returns adaptive criterion impl.
     *
     * @param weighted Weighted tree used by the operation.
     * @param criterion Per-node selection criterion.
     * @param delta Delta offset used by the operation.
     * @return Adaptive criterion impl.
     */
    static std::vector<bool> getAdaptiveCriterionImpl(const WeightedMorphologicalTree<T>& weighted, std::vector<bool>& criterion, AltitudeDiff<T> delta) {
        const MorphologicalTree& tree = weighted.topology();
        MSERComputer<T> mser(weighted);
        (void)mser.computeMSER(delta);
        return adaptiveCriterionFromMaskStability<float>(tree, criterion, mser, "AttributeFilters::getAdaptiveCriterion");
    }

    /**
     * @brief Returns adaptive criterion by depth impl.
     *
     * @param tree Tree topology used by the operation.
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @param depthDelta Topological-depth radius of the stability window.
     * @return Adaptive criterion by depth impl.
     */
    template <std::floating_point Real>
    static std::vector<bool> getAdaptiveCriterionByDepthImpl(const MorphologicalTree& tree, const Real* attribute, Real threshold, int depthDelta) {
        DepthStableRegionComputer<Real> stabilityComputer(tree);
        (void)stabilityComputer.computeByDepth(depthDelta);
        return adaptiveCriterionFromAttributeStability<Real>(tree, attribute, threshold, stabilityComputer, "AttributeFilters::getAdaptiveCriterionByDepth");
    }

    /**
     * @brief Returns adaptive criterion by depth impl.
     *
     * @param tree Tree topology used by the operation.
     * @param criterion Per-node selection criterion.
     * @param depthDelta Topological-depth radius of the stability window.
     * @return Adaptive criterion by depth impl.
     */
    static std::vector<bool> getAdaptiveCriterionByDepthImpl(const MorphologicalTree& tree, std::vector<bool>& criterion, int depthDelta) {
        DepthStableRegionComputer<float> stabilityComputer(tree);
        (void)stabilityComputer.computeByDepth(depthDelta);
        return adaptiveCriterionFromMaskStability<float>(tree, criterion, stabilityComputer, "AttributeFilters::getAdaptiveCriterionByDepth");
    }
    /// @endcond

  public:
    /**
     * @brief Creates filtering operators over a non-owning weighted tree view.
     *
     * The view must remain valid for the lifetime of the filter object. Public
     * methods reject use after the underlying topology mutates.
     *
     * @param view Tree view used by the operation.
     */
    explicit AttributeFilters(AltitudeView view) : view_{view}, tree{view_.topology()}, treeMutationVersion_{tree.getMutationVersion()} {
        MMCFILTERS_CONTRACT_CHECKED_ONLY(view_.requireTopologyUnchanged("AttributeFilters"));
    }

    /**
     * @brief Creates filtering operators over an owned weighted tree.
     *
     * This overload keeps a pointer to the owner so methods that require
     * tree-owned altitude state, such as MSER-assisted adaptive criteria, can be
     * used.
     *
     * @param weighted Weighted tree used by the operation.
     */
    explicit AttributeFilters(const WeightedMorphologicalTree<T>& weighted) : AttributeFilters(weighted.asView()) { weighted_ = &weighted; }

    /**
     * @brief Destroys the attribute-filter facade.
     */
    ~AttributeFilters() = default;

    /**
     * @brief Builds an MSER-assisted pruning criterion from an existing keep/reject mask.
     *
     * @param criterion Dense internal-node mask. Nodes with `false` are candidates
     * for adaptive pruning.
     * @param delta Altitude distance used by the MSER stability computation.
     * @return Dense internal-node pruning mask.
     */
    [[nodiscard]] std::vector<bool> getAdaptiveCriterion(std::vector<bool>& criterion, AltitudeDiff<T> delta) {
        requireStableTree("AttributeFilters::getAdaptiveCriterion");
        if (weighted_ == nullptr) {
            throw std::logic_error(
                "AttributeFilters::getAdaptiveCriterion requires a WeightedMorphologicalTree owner because MSER uses the tree-owned altitude.");
        }
        return AttributeFilters::getAdaptiveCriterionImpl(*weighted_, criterion, delta);
    }

    /**
     * @brief Builds a depth-stability pruning criterion from an existing mask.
     *
     * `depthDelta` is a number of tree edges. The stability computation does not
     * read altitude, so it is suitable for tree-of-shapes and self-dual residual
     * trees where max/min altitude polarity is not defined.
     *
     * @param criterion Per-node selection criterion.
     * @param depthDelta Topological depth offset used by the operation.
     * @return The resulting depth-stability pruning criterion from an existing mask.
     */
    [[nodiscard]] std::vector<bool> getAdaptiveCriterionByDepth(std::vector<bool>& criterion, int depthDelta) {
        requireStableTree("AttributeFilters::getAdaptiveCriterionByDepth");
        return AttributeFilters::getAdaptiveCriterionByDepthImpl(this->tree, criterion, depthDelta);
    }

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
        ImagePtr<T> imgOutput = Image<T>::create(this->tree.getNumRowsOfGridDomain2D(), this->tree.getNumColsOfGridDomain2D());
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
        ImagePtr<T> imgOutput = Image<T>::create(this->tree.getNumRowsOfGridDomain2D(), this->tree.getNumColsOfGridDomain2D());
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
     * rule using that connected keep criterion.
     *
     * @param attr Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @return Image produced by the operation.
     */
    template <std::floating_point Real> [[nodiscard]] ImagePtr<T> filteringByViterbiRule(const Real* attr, Real threshold) {
        requireStableTree("AttributeFilters::filteringByViterbiRule");
        auto costs = detail::makeThresholdViterbiCosts(tree, attr, threshold);
        std::vector<bool> criterion = detail::computeViterbiKeepCriterion(tree, costs);

        ImagePtr<T> imgOutput = Image<T>::create(this->tree.getNumRowsOfGridDomain2D(), this->tree.getNumColsOfGridDomain2D());
        detail::kernel::filterDirect(view(), criterion, imgOutput);
        return imgOutput;
    }

    /**
     * @brief Applies pruning-min filtering from a dense internal-node criterion.
     *
     * @param criterion Per-node selection criterion.
     * @return Image produced by the operation.
     */
    [[nodiscard]] ImagePtr<T> filteringByPruningMin(std::vector<bool>& criterion) {
        ImagePtr<T> imgOutput = Image<T>::create(this->tree.getNumRowsOfGridDomain2D(), this->tree.getNumColsOfGridDomain2D());
        filteringByPruningMinCriterionImpl(view(), criterion, imgOutput);
        return imgOutput;
    }

    /**
     * @brief Applies pruning-max filtering from a dense internal-node criterion.
     *
     * @param criterion Per-node selection criterion.
     * @return Image produced by the operation.
     */
    [[nodiscard]] ImagePtr<T> filteringByPruningMax(std::vector<bool>& criterion) {
        ImagePtr<T> imgOutput = Image<T>::create(this->tree.getNumRowsOfGridDomain2D(), this->tree.getNumColsOfGridDomain2D());
        filteringByPruningMaxCriterionImpl(view(), criterion, imgOutput);
        return imgOutput;
    }

    /**
     * @brief Applies the direct filtering rule from a dense internal-node criterion.
     *
     * Accepted nodes use their own altitude; rejected nodes inherit the filtered
     * level propagated from their parent.
     *
     * @param criterion Per-node selection criterion.
     * @return Image produced by the operation.
     */
    [[nodiscard]] ImagePtr<T> filteringByDirectRule(std::vector<bool>& criterion) {
        requireStableTree("AttributeFilters::filteringByDirectRule");
        requireCriterionSize(tree, criterion, "AttributeFilters::filteringByDirectRule");
        ImagePtr<T> imgOutput = detail::CommittedImageAccess::create<T>(this->tree.getNumRowsOfGridDomain2D(), this->tree.getNumColsOfGridDomain2D());
        detail::kernel::filterDirect(view(), criterion, imgOutput);
        return imgOutput;
    }

    /**
     * @brief Applies the subtractive filtering rule from a dense internal-node criterion.
     *
     * Accepted nodes add their residue to the propagated reconstruction level;
     * rejected nodes inherit the parent level.
     *
     * @param criterion Per-node selection criterion.
     * @return Image produced by the operation.
     */
    [[nodiscard]] ImagePtr<T> filteringBySubtractiveRule(std::vector<bool>& criterion) {
        ImagePtr<T> imgOutput = Image<T>::create(this->tree.getNumRowsOfGridDomain2D(), this->tree.getNumColsOfGridDomain2D());
        filteringBySubtractiveRuleImpl(view(), criterion, imgOutput);
        return imgOutput;
    }

    /**
     * @brief Applies the subtractive score rule from dense per-node scores.
     *
     * Each node residue is weighted by `prob[nodeId]` before accumulation, and the
     * result is reconstructed as a floating-point image.
     *
     * @param prob Per-node probability values.
     * @return Image produced by the operation.
     */
    [[nodiscard]] ImageFloatPtr filteringBySubtractiveScoreRule(std::vector<float>& prob) {
        ImageFloatPtr imgOutput = ImageFloat::create(this->tree.getNumRowsOfGridDomain2D(), this->tree.getNumColsOfGridDomain2D());
        filteringBySubtractiveScoreRuleImpl(view(), prob, imgOutput);
        return imgOutput;
    }

    /**
     * @brief Writes subtractive-score filtering into a caller-owned output image.
     *
     * @param weighted Weighted tree used by the operation.
     * @param prob Per-node probability values.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    static void filteringBySubtractiveScoreRule(const WeightedTreeView<T>& weighted, std::vector<float>& prob, ImageFloatPtr imgOutputPtr) {
        filteringBySubtractiveScoreRuleImpl(weighted, prob, imgOutputPtr);
    }

    /**
     * @brief Writes subtractive-score filtering from a weighted owner into an output image.
     *
     * @param weighted Weighted tree used by the operation.
     * @param prob Per-node probability values.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    static void filteringBySubtractiveScoreRule(const WeightedMorphologicalTree<T>& weighted, std::vector<float>& prob, ImageFloatPtr imgOutputPtr) {
        filteringBySubtractiveScoreRule(weighted.asView(), prob, imgOutputPtr);
    }

    /**
     * @brief Writes subtractive-rule filtering into a caller-owned output image.
     *
     * @param weighted Weighted tree used by the operation.
     * @param criterion Per-node selection criterion.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    static void filteringBySubtractiveRule(const WeightedTreeView<T>& weighted, std::vector<bool>& criterion, ImagePtr<T> imgOutputPtr) {
        filteringBySubtractiveRuleImpl(weighted, criterion, imgOutputPtr);
    }

    /**
     * @brief Writes subtractive-rule filtering from a weighted owner into an output image.
     *
     * @param weighted Weighted tree used by the operation.
     * @param criterion Per-node selection criterion.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    static void filteringBySubtractiveRule(const WeightedMorphologicalTree<T>& weighted, std::vector<bool>& criterion, ImagePtr<T> imgOutputPtr) {
        filteringBySubtractiveRule(weighted.asView(), criterion, imgOutputPtr);
    }

    /**
     * @brief Writes direct-rule filtering into a caller-owned output image.
     *
     * @param weighted Weighted tree used by the operation.
     * @param criterion Per-node selection criterion.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    static void filteringByDirectRule(const WeightedTreeView<T>& weighted, std::vector<bool>& criterion, ImagePtr<T> imgOutputPtr) {
        MMCFILTERS_CONTRACT_CHECKED_ONLY(weighted.requireTopologyUnchanged("AttributeFilters::filteringByDirectRule"));
        requireCriterionSize(weighted.topology(), criterion, "AttributeFilters::filteringByDirectRule");
        requireOutputImage(weighted.topology(), imgOutputPtr, "AttributeFilters::filteringByDirectRule");
        detail::kernel::filterDirect(weighted, criterion, imgOutputPtr);
    }

    /**
     * @brief Writes direct-rule filtering from a weighted owner into an output image.
     *
     * @param weighted Weighted tree used by the operation.
     * @param criterion Per-node selection criterion.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    static void filteringByDirectRule(const WeightedMorphologicalTree<T>& weighted, std::vector<bool>& criterion, ImagePtr<T> imgOutputPtr) {
        const AltitudeView view = weighted.asView();
        MMCFILTERS_CONTRACT_CHECKED_ONLY(view.requireTopologyUnchanged("AttributeFilters::filteringByDirectRule"));
        requireCriterionSize(view.topology(), criterion, "AttributeFilters::filteringByDirectRule");
        requireOutputImage(view.topology(), imgOutputPtr, "AttributeFilters::filteringByDirectRule");
        detail::kernel::filterDirect(view, criterion, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-min filtering from a criterion into an output image.
     *
     * @param weighted Weighted tree used by the operation.
     * @param criterion Per-node selection criterion.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    static void filteringByPruningMin(const WeightedTreeView<T>& weighted, std::vector<bool>& criterion, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMinCriterionImpl(weighted, criterion, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-min filtering from a weighted owner into an output image.
     *
     * @param weighted Weighted tree used by the operation.
     * @param criterion Per-node selection criterion.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    static void filteringByPruningMin(const WeightedMorphologicalTree<T>& weighted, std::vector<bool>& criterion, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMin(weighted.asView(), criterion, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-max filtering from a criterion into an output image.
     *
     * @param weighted Weighted tree used by the operation.
     * @param criterion Per-node selection criterion.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    static void filteringByPruningMax(const WeightedTreeView<T>& weighted, std::vector<bool>& criterion, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMaxCriterionImpl(weighted, criterion, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-max filtering from a weighted owner into an output image.
     *
     * @param weighted Weighted tree used by the operation.
     * @param criterion Per-node selection criterion.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    static void filteringByPruningMax(const WeightedMorphologicalTree<T>& weighted, std::vector<bool>& criterion, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMax(weighted.asView(), criterion, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-min filtering from an owned attribute buffer into an output image.
     *
     * @param weighted Weighted tree used by the operation.
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    template <std::floating_point Real>
    static void filteringByPruningMin(const WeightedTreeView<T>& weighted, const std::shared_ptr<Real[]>& attribute, Real threshold, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMin(weighted, attribute.get(), threshold, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-min filtering from a weighted owner and owned attribute buffer.
     *
     * @param weighted Weighted tree used by the operation.
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    template <std::floating_point Real>
    static void filteringByPruningMin(const WeightedMorphologicalTree<T>& weighted, const std::shared_ptr<Real[]>& attribute, Real threshold,
                                      ImagePtr<T> imgOutputPtr) {
        filteringByPruningMin(weighted.asView(), attribute.get(), threshold, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-min filtering from a raw attribute buffer into an output image.
     *
     * @param weighted Weighted tree used by the operation.
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    template <std::floating_point Real>
    static void filteringByPruningMin(const WeightedTreeView<T>& weighted, const Real* attribute, Real threshold, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMinAttributeImpl(weighted, attribute, threshold, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-min filtering from a weighted owner and raw attribute buffer.
     *
     * @param weighted Weighted tree used by the operation.
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    template <std::floating_point Real>
    static void filteringByPruningMin(const WeightedMorphologicalTree<T>& weighted, const Real* attribute, Real threshold, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMin(weighted.asView(), attribute, threshold, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-max filtering from an owned attribute buffer into an output image.
     *
     * @param weighted Weighted tree used by the operation.
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    template <std::floating_point Real>
    static void filteringByPruningMax(const WeightedTreeView<T>& weighted, const std::shared_ptr<Real[]>& attribute, Real threshold, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMax(weighted, attribute.get(), threshold, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-max filtering from a weighted owner and owned attribute buffer.
     *
     * @param weighted Weighted tree used by the operation.
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    template <std::floating_point Real>
    static void filteringByPruningMax(const WeightedMorphologicalTree<T>& weighted, const std::shared_ptr<Real[]>& attribute, Real threshold,
                                      ImagePtr<T> imgOutputPtr) {
        filteringByPruningMax(weighted.asView(), attribute.get(), threshold, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-max filtering from a raw attribute buffer into an output image.
     *
     * @param weighted Weighted tree used by the operation.
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    template <std::floating_point Real>
    static void filteringByPruningMax(const WeightedTreeView<T>& weighted, const Real* attribute, Real threshold, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMaxAttributeImpl(weighted, attribute, threshold, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-max filtering from a weighted owner and raw attribute buffer.
     *
     * @param weighted Weighted tree used by the operation.
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @param imgOutputPtr Output image receiving the reconstruction.
     */
    template <std::floating_point Real>
    static void filteringByPruningMax(const WeightedMorphologicalTree<T>& weighted, const Real* attribute, Real threshold, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMax(weighted.asView(), attribute, threshold, imgOutputPtr);
    }

    /**
     * @brief Builds an MSER-assisted pruning criterion from an attribute threshold.
     *
     * @param weighted Weighted tree used by the operation.
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @param delta Delta offset or radius used by the operation.
     * @return The resulting MSER-assisted pruning criterion from an attribute threshold.
     */
    template <std::floating_point Real>
    [[nodiscard]] static std::vector<bool> getAdaptiveCriterion(const WeightedMorphologicalTree<T>& weighted, const std::shared_ptr<Real[]>& attribute,
                                                                Real threshold, AltitudeDiff<T> delta) {
        return getAdaptiveCriterionImpl(weighted, attribute.get(), threshold, delta);
    }

    /**
     * @brief Builds an MSER-assisted pruning criterion from a raw attribute buffer.
     *
     * @param weighted Weighted tree used by the operation.
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @param delta Delta offset or radius used by the operation.
     * @return The resulting MSER-assisted pruning criterion from a raw attribute buffer.
     */
    template <std::floating_point Real>
    [[nodiscard]] static std::vector<bool> getAdaptiveCriterion(const WeightedMorphologicalTree<T>& weighted, const Real* attribute, Real threshold,
                                                                AltitudeDiff<T> delta) {
        return getAdaptiveCriterionImpl(weighted, attribute, threshold, delta);
    }

    /**
     * @brief Builds an MSER-assisted pruning criterion from an existing criterion mask.
     *
     * @param weighted Weighted tree used by the operation.
     * @param criterion Per-node selection criterion.
     * @param delta Delta offset or radius used by the operation.
     * @return The resulting MSER-assisted pruning criterion from an existing criterion mask.
     */
    [[nodiscard]] static std::vector<bool> getAdaptiveCriterion(const WeightedMorphologicalTree<T>& weighted, std::vector<bool>& criterion,
                                                                AltitudeDiff<T> delta) {
        return getAdaptiveCriterionImpl(weighted, criterion, delta);
    }

    /**
     * @brief Builds a depth-stability pruning criterion from an attribute threshold.
     *
     * @param weighted Weighted tree used by the operation.
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @param depthDelta Topological depth offset used by the operation.
     * @return The resulting depth-stability pruning criterion from an attribute threshold.
     */
    template <std::floating_point Real>
    [[nodiscard]] static std::vector<bool> getAdaptiveCriterionByDepth(const WeightedMorphologicalTree<T>& weighted, const std::shared_ptr<Real[]>& attribute,
                                                                       Real threshold, int depthDelta) {
        return getAdaptiveCriterionByDepthImpl(weighted.topology(), attribute.get(), threshold, depthDelta);
    }

    /**
     * @brief Builds a depth-stability pruning criterion from a raw attribute buffer.
     *
     * @param weighted Weighted tree used by the operation.
     * @param attribute Attribute requested by the operation.
     * @param threshold Threshold applied by the operation.
     * @param depthDelta Topological depth offset used by the operation.
     * @return The resulting depth-stability pruning criterion from a raw attribute buffer.
     */
    template <std::floating_point Real>
    [[nodiscard]] static std::vector<bool> getAdaptiveCriterionByDepth(const WeightedMorphologicalTree<T>& weighted, const Real* attribute, Real threshold,
                                                                       int depthDelta) {
        return getAdaptiveCriterionByDepthImpl(weighted.topology(), attribute, threshold, depthDelta);
    }

    /**
     * @brief Builds a depth-stability pruning criterion from an existing criterion mask.
     *
     * @param weighted Weighted tree used by the operation.
     * @param criterion Per-node selection criterion.
     * @param depthDelta Topological depth offset used by the operation.
     * @return The resulting depth-stability pruning criterion from an existing criterion mask.
     */
    [[nodiscard]] static std::vector<bool> getAdaptiveCriterionByDepth(const WeightedMorphologicalTree<T>& weighted, std::vector<bool>& criterion,
                                                                       int depthDelta) {
        return getAdaptiveCriterionByDepthImpl(weighted.topology(), criterion, depthDelta);
    }
};

} // namespace mmcfilters
