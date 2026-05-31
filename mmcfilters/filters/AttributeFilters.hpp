#pragma once

#include "../utils/Image.hpp"
#include "../utils/Common.hpp"
#include "../trees/TreeAltitudeAlgorithms.hpp"
#include "../trees/WeightedMorphologicalTree.hpp"
#include "../trees/WeightedTreeView.hpp"
#include "../trees/detail/TreeTraversalDetail.hpp"
#include "DepthStableRegionComputer.hpp"
#include "MSERComputer.hpp"
#include "detail/VariationMeasure.hpp"
#include "detail/ViterbiDecision.hpp"

#include <cassert>
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
template<AltitudeValue T>
class AttributeFilters {
protected:
    /// @cond INTERNAL
    using AltitudeView = WeightedTreeView<T>;

    AltitudeView view_;
    const WeightedMorphologicalTree<T>* weighted_ = nullptr;
    const MorphologicalTree& tree;
    std::size_t treeMutationVersion_ = 0;

    AltitudeView view() const {
        return weighted_ != nullptr ? weighted_->asView() : view_;
    }

    void requireStableTree(const char* context) const {
        tree.requireMutationVersion(treeMutationVersion_, context);
    }

    template <std::floating_point Real>
    static void requireAttributePointer(const Real* attribute, const char* context) {
        if (attribute == nullptr) {
            throw std::invalid_argument(std::string(context) + " requires a non-null attribute buffer.");
        }
    }

    static void requireCriterionSize(const MorphologicalTree& tree, const std::vector<bool>& criterion, const char* context) {
        if (criterion.size() != static_cast<std::size_t>(tree.getNumInternalNodeSlots())) {
            throw std::invalid_argument(std::string(context) + " criterion size must match the internal node slot count.");
        }
    }

    static void requireScoreSize(const MorphologicalTree& tree, const std::vector<float>& scores, const char* context) {
        if (scores.size() != static_cast<std::size_t>(tree.getNumInternalNodeSlots())) {
            throw std::invalid_argument(std::string(context) + " score size must match the internal node slot count.");
        }
    }

    template <typename TImagePtr>
    static void requireOutputImage(const MorphologicalTree& tree, const TImagePtr& image, const char* context) {
        if (!image) {
            throw std::invalid_argument(std::string(context) + " requires a non-null output image.");
        }
        if (image->getNumRows() != tree.getNumRowsOfImage() || image->getNumCols() != tree.getNumColsOfImage()) {
            throw std::invalid_argument(std::string(context) + " output image shape must match the tree image domain.");
        }
    }

    static T altitudeOf(const AltitudeView& view, NodeId nodeId) {
        return view.getAltitude(nodeId);
    }

    static AltitudeDiff<T> residueOf(const AltitudeView& view, NodeId nodeId) {
        return view.getNodeResidue(nodeId);
    }

    /**
     * @brief Writes one node's direct proper parts only.
     *
     * Reconstruction is done by writing every alive node's direct proper parts
     * exactly once. Descendant supports are intentionally not touched here; use
     * `writeSubtreeProperParts` only when a pruning rule collapses a whole
     * branch to a single gray level.
     */
    template <typename TValue>
    static void writeProperParts(const MorphologicalTree& tree, NodeId nodeId, TValue* output, TValue value) {
        for (int pixel : tree.getProperParts(nodeId)) {
            output[pixel] = value;
        }
    }

    /**
     * @brief Paints all direct proper parts owned by nodes in one subtree.
     *
     * This is the physical image-domain effect of pruning a branch: every node
     * below the cut, including the cut node itself, receives one replacement
     * altitude chosen by the specific pruning convention.
     */
    template <typename TValue>
    static void writeSubtreeProperParts(const MorphologicalTree& tree, NodeId nodeId, TValue* output, TValue value) {
        for (NodeId subtreeNodeId : tree.getNodeSubtree(nodeId)) {
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
            for (NodeId childNodeId : tree.getChildren(nodeId)) {
                mapLevel[childNodeId] =
                    mapLevel[nodeId] + (static_cast<float>(residueOf(view, childNodeId)) * prob[childNodeId]);
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
            for (NodeId childNodeId : tree.getChildren(nodeId)) {
                mapLevel[childNodeId] = criterion[childNodeId]
                    ? static_cast<AltitudeDiff<T>>(mapLevel[nodeId] + residueOf(view, childNodeId))
                    : mapLevel[nodeId];
                stack.push(childNodeId);
            }
        }

        auto imgOutput = imgOutputPtr->rawData();
        for (NodeId nodeId : tree.getAliveNodeIds()) {
            writeProperParts(tree, nodeId, imgOutput, static_cast<T>(mapLevel[nodeId]));
        }
    }

    /**
     * @brief Direct reconstruction from a dense keep criterion.
     *
     * A kept node selects its own altitude. A rejected node does not remove its
     * descendants; it only forwards the filtered parent level, so accepted
     * descendants can still reintroduce their own levels.
     */
    static void filteringByDirectRuleImpl(AltitudeView view, std::vector<bool>& criterion, ImagePtr<T> imgOutputPtr) {
        const char* context = "AttributeFilters::filteringByDirectRule";
        view.requireTopologyUnchanged(context);
        const MorphologicalTree& tree = view.topology();
        requireCriterionSize(tree, criterion, context);
        requireOutputImage(tree, imgOutputPtr, context);
        std::vector<T> mapLevel(static_cast<std::size_t>(tree.getNumInternalNodeSlots()));
        const NodeId rootNodeId = tree.getRoot();
        mapLevel[rootNodeId] = altitudeOf(view, rootNodeId);

        // Direct filtering is path-local. Traverse root-to-leaf so every child
        // sees the filtered level chosen for its parent.
        std::stack<NodeId> stack;
        stack.push(rootNodeId);
        while (!stack.empty()) {
            const NodeId nodeId = stack.top();
            stack.pop();
            for (NodeId childNodeId : tree.getChildren(nodeId)) {
                mapLevel[childNodeId] = criterion[childNodeId] ? altitudeOf(view, childNodeId) : mapLevel[nodeId];
                stack.push(childNodeId);
            }
        }

        auto imgOutput = imgOutputPtr->rawData();
        for (NodeId nodeId : tree.getAliveNodeIds()) {
            writeProperParts(tree, nodeId, imgOutput, mapLevel[nodeId]);
        }
    }

    /**
     * @brief Pruning-min reconstruction from an explicit keep criterion.
     *
     * The convention used here is ancestor-level pruning: once a child branch is
     * rejected, the whole branch is painted at the current accepted node level
     * and no deeper descendants are evaluated. This matches the attribute
     * overload where `criterion[node] == (attribute[node] > threshold)`.
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
            for (NodeId childNodeId : tree.getChildren(nodeId)) {
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
            tree,
            tree.getRoot(),
            [&criterion, &keepCriterion](NodeId nodeId) -> void {
                criterion[nodeId] = !keepCriterion[nodeId];
            },
            [&criterion](NodeId parentNodeId, NodeId childNodeId) -> void {
                criterion[parentNodeId] = (criterion[parentNodeId] & criterion[childNodeId]);
            },
            [](NodeId) -> void {});

        auto imgOutput = imgOutputPtr->rawData();
        std::stack<NodeId> stack;
        stack.push(tree.getRoot());
        while (!stack.empty()) {
            const NodeId nodeId = stack.top();
            stack.pop();
            writeProperParts(tree, nodeId, imgOutput, altitudeOf(view, nodeId));
            for (NodeId childNodeId : tree.getChildren(nodeId)) {
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
            for (NodeId childNodeId : tree.getChildren(nodeId)) {
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
            tree,
            tree.getRoot(),
            [&criterion, attribute, threshold](NodeId nodeId) -> void {
                if (attribute[nodeId] <= threshold) {
                    criterion[nodeId] = true;
                }
            },
            [&criterion](NodeId parentNodeId, NodeId childNodeId) -> void {
                criterion[parentNodeId] = (criterion[parentNodeId] & criterion[childNodeId]);
            },
            [](NodeId) -> void {});

        auto imgOutput = imgOutputPtr->rawData();
        std::stack<NodeId> stack;
        stack.push(tree.getRoot());
        while (!stack.empty()) {
            const NodeId nodeId = stack.top();
            stack.pop();
            writeProperParts(tree, nodeId, imgOutput, altitudeOf(view, nodeId));
            for (NodeId childNodeId : tree.getChildren(nodeId)) {
                if (!criterion[childNodeId]) {
                    stack.push(childNodeId);
                } else {
                    writeSubtreeProperParts(tree, childNodeId, imgOutput, altitudeOf(view, childNodeId));
                }
            }
        }
    }

    template <std::floating_point Real, class StabilityComputer>
    static std::vector<bool> adaptiveCriterionFromAttributeStability(
        const MorphologicalTree& tree,
        const Real* attribute,
        Real threshold,
        StabilityComputer& stabilityComputer,
        const char* context) {
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

    template <std::floating_point Real, class StabilityComputer>
    static std::vector<bool> adaptiveCriterionFromMaskStability(
        const MorphologicalTree& tree,
        std::vector<bool>& criterion,
        StabilityComputer& stabilityComputer,
        const char* context) {
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

    template <std::floating_point Real>
    static std::vector<bool> getAdaptiveCriterionImpl(const WeightedMorphologicalTree<T>& weighted, const Real* attribute, Real threshold, AltitudeDiff<T> delta) {
        const MorphologicalTree& tree = weighted.topology();
        MSERComputer<T, Real> mser(weighted);
        (void)mser.computeMSER(delta);
        return adaptiveCriterionFromAttributeStability<Real>(
            tree,
            attribute,
            threshold,
            mser,
            "AttributeFilters::getAdaptiveCriterion");
    }

    static std::vector<bool> getAdaptiveCriterionImpl(const WeightedMorphologicalTree<T>& weighted, std::vector<bool>& criterion, AltitudeDiff<T> delta) {
        const MorphologicalTree& tree = weighted.topology();
        MSERComputer<T> mser(weighted);
        (void)mser.computeMSER(delta);
        return adaptiveCriterionFromMaskStability<float>(
            tree,
            criterion,
            mser,
            "AttributeFilters::getAdaptiveCriterion");
    }

    template <std::floating_point Real>
    static std::vector<bool> getAdaptiveCriterionByDepthImpl(const MorphologicalTree& tree, const Real* attribute, Real threshold, int depthDelta) {
        DepthStableRegionComputer<Real> stabilityComputer(tree);
        (void)stabilityComputer.computeByDepth(depthDelta);
        return adaptiveCriterionFromAttributeStability<Real>(
            tree,
            attribute,
            threshold,
            stabilityComputer,
            "AttributeFilters::getAdaptiveCriterionByDepth");
    }

    static std::vector<bool> getAdaptiveCriterionByDepthImpl(const MorphologicalTree& tree, std::vector<bool>& criterion, int depthDelta) {
        DepthStableRegionComputer<float> stabilityComputer(tree);
        (void)stabilityComputer.computeByDepth(depthDelta);
        return adaptiveCriterionFromMaskStability<float>(
            tree,
            criterion,
            stabilityComputer,
            "AttributeFilters::getAdaptiveCriterionByDepth");
    }
    /// @endcond

public:
    /**
     * @brief Creates filtering operators over a non-owning weighted tree view.
     *
     * The view must remain valid for the lifetime of the filter object. Public
     * methods reject use after the underlying topology mutates.
     */
    explicit AttributeFilters(AltitudeView view)
        : view_{view},
          tree{view_.topology()},
          treeMutationVersion_{tree.getMutationVersion()} {
        view_.requireTopologyUnchanged("AttributeFilters");
    }

    /**
     * @brief Creates filtering operators over an owned weighted tree.
     *
     * This overload keeps a pointer to the owner so methods that require
     * tree-owned altitude state, such as MSER-assisted adaptive criteria, can be
     * used.
     */
    explicit AttributeFilters(const WeightedMorphologicalTree<T>& weighted)
        : AttributeFilters(weighted.asView()) {
        weighted_ = &weighted;
    }

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
            throw std::logic_error("AttributeFilters::getAdaptiveCriterion requires a WeightedMorphologicalTree owner because MSER uses the tree-owned altitude.");
        }
        return AttributeFilters::getAdaptiveCriterionImpl(*weighted_, criterion, delta);
    }

    /**
     * @brief Builds a depth-stability pruning criterion from an existing mask.
     *
     * `depthDelta` is a number of tree edges. The stability computation does not
     * read altitude, so it is suitable for tree-of-shapes and self-dual residual
     * trees where max/min altitude polarity is not defined.
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
     */
    template <std::floating_point Real>
    [[nodiscard]] ImagePtr<T> filteringByPruningMin(const std::shared_ptr<Real[]>& attr, Real threshold) {
        return filteringByPruningMin(attr.get(), threshold);
    }

    /**
     * @brief Applies pruning-min filtering from a raw internal-node attribute buffer.
     */
    template <std::floating_point Real>
    [[nodiscard]] ImagePtr<T> filteringByPruningMin(const Real* attr, Real threshold) {
        requireStableTree("AttributeFilters::filteringByPruningMin");
        assert(attr != nullptr);
        ImagePtr<T> imgOutput = Image<T>::create(this->tree.getNumRowsOfImage(), this->tree.getNumColsOfImage());
        filteringByPruningMinAttributeImpl(view(), attr, threshold, imgOutput);
        return imgOutput;
    }

    /**
     * @brief Applies pruning-max filtering from an attribute buffer.
     *
     * Nodes with attribute values above `threshold` are kept, and fully rejected
     * subtrees are reconstructed at their own subtree levels.
     */
    template <std::floating_point Real>
    [[nodiscard]] ImagePtr<T> filteringByPruningMax(const std::shared_ptr<Real[]>& attr, Real threshold) {
        return filteringByPruningMax(attr.get(), threshold);
    }

    /**
     * @brief Applies pruning-max filtering from a raw internal-node attribute buffer.
     */
    template <std::floating_point Real>
    [[nodiscard]] ImagePtr<T> filteringByPruningMax(const Real* attr, Real threshold) {
        requireStableTree("AttributeFilters::filteringByPruningMax");
        assert(attr != nullptr);
        ImagePtr<T> imgOutput = Image<T>::create(this->tree.getNumRowsOfImage(), this->tree.getNumColsOfImage());
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
     */
    template <std::floating_point Real>
    [[nodiscard]] ImagePtr<T> filteringByViterbiRule(const Real* attr, Real threshold) {
        requireStableTree("AttributeFilters::filteringByViterbiRule");
        auto costs = detail::makeThresholdViterbiCosts(tree, attr, threshold);
        std::vector<bool> criterion = detail::computeViterbiKeepCriterion(tree, costs);

        ImagePtr<T> imgOutput = Image<T>::create(this->tree.getNumRowsOfImage(), this->tree.getNumColsOfImage());
        filteringByDirectRuleImpl(view(), criterion, imgOutput);
        return imgOutput;
    }

    /**
     * @brief Applies pruning-min filtering from a dense internal-node criterion.
     */
    [[nodiscard]] ImagePtr<T> filteringByPruningMin(std::vector<bool>& criterion) {
        requireStableTree("AttributeFilters::filteringByPruningMin");
        ImagePtr<T> imgOutput = Image<T>::create(this->tree.getNumRowsOfImage(), this->tree.getNumColsOfImage());
        filteringByPruningMinCriterionImpl(view(), criterion, imgOutput);
        return imgOutput;
    }

    /**
     * @brief Applies pruning-max filtering from a dense internal-node criterion.
     */
    [[nodiscard]] ImagePtr<T> filteringByPruningMax(std::vector<bool>& criterion) {
        requireStableTree("AttributeFilters::filteringByPruningMax");
        ImagePtr<T> imgOutput = Image<T>::create(this->tree.getNumRowsOfImage(), this->tree.getNumColsOfImage());
        filteringByPruningMaxCriterionImpl(view(), criterion, imgOutput);
        return imgOutput;
    }

    /**
     * @brief Applies the direct filtering rule from a dense internal-node criterion.
     *
     * Accepted nodes use their own altitude; rejected nodes inherit the filtered
     * level propagated from their parent.
     */
    [[nodiscard]] ImagePtr<T> filteringByDirectRule(std::vector<bool>& criterion) {
        requireStableTree("AttributeFilters::filteringByDirectRule");
        ImagePtr<T> imgOutput = Image<T>::create(this->tree.getNumRowsOfImage(), this->tree.getNumColsOfImage());
        filteringByDirectRuleImpl(view(), criterion, imgOutput);
        return imgOutput;
    }

    /**
     * @brief Applies the subtractive filtering rule from a dense internal-node criterion.
     *
     * Accepted nodes add their residue to the propagated reconstruction level;
     * rejected nodes inherit the parent level.
     */
    [[nodiscard]] ImagePtr<T> filteringBySubtractiveRule(std::vector<bool>& criterion) {
        requireStableTree("AttributeFilters::filteringBySubtractiveRule");
        ImagePtr<T> imgOutput = Image<T>::create(this->tree.getNumRowsOfImage(), this->tree.getNumColsOfImage());
        filteringBySubtractiveRuleImpl(view(), criterion, imgOutput);
        return imgOutput;
    }

    /**
     * @brief Applies the subtractive score rule from dense per-node scores.
     *
     * Each node residue is weighted by `prob[nodeId]` before accumulation, and the
     * result is reconstructed as a floating-point image.
     */
    [[nodiscard]] ImageFloatPtr filteringBySubtractiveScoreRule(std::vector<float>& prob) {
        requireStableTree("AttributeFilters::filteringBySubtractiveScoreRule");
        ImageFloatPtr imgOutput = ImageFloat::create(this->tree.getNumRowsOfImage(), this->tree.getNumColsOfImage());
        filteringBySubtractiveScoreRuleImpl(view(), prob, imgOutput);
        return imgOutput;
    }

    /**
     * @brief Writes subtractive-score filtering into a caller-owned output image.
     */
    static void filteringBySubtractiveScoreRule(const WeightedTreeView<T>& weighted, std::vector<float>& prob, ImageFloatPtr imgOutputPtr) {
        filteringBySubtractiveScoreRuleImpl(weighted, prob, imgOutputPtr);
    }

    /**
     * @brief Writes subtractive-score filtering from a weighted owner into an output image.
     */
    static void filteringBySubtractiveScoreRule(const WeightedMorphologicalTree<T>& weighted, std::vector<float>& prob, ImageFloatPtr imgOutputPtr) {
        filteringBySubtractiveScoreRule(weighted.asView(), prob, imgOutputPtr);
    }

    /**
     * @brief Writes subtractive-rule filtering into a caller-owned output image.
     */
    static void filteringBySubtractiveRule(const WeightedTreeView<T>& weighted, std::vector<bool>& criterion, ImagePtr<T> imgOutputPtr) {
        filteringBySubtractiveRuleImpl(weighted, criterion, imgOutputPtr);
    }

    /**
     * @brief Writes subtractive-rule filtering from a weighted owner into an output image.
     */
    static void filteringBySubtractiveRule(const WeightedMorphologicalTree<T>& weighted, std::vector<bool>& criterion, ImagePtr<T> imgOutputPtr) {
        filteringBySubtractiveRule(weighted.asView(), criterion, imgOutputPtr);
    }

    /**
     * @brief Writes direct-rule filtering into a caller-owned output image.
     */
    static void filteringByDirectRule(const WeightedTreeView<T>& weighted, std::vector<bool>& criterion, ImagePtr<T> imgOutputPtr) {
        filteringByDirectRuleImpl(weighted, criterion, imgOutputPtr);
    }

    /**
     * @brief Writes direct-rule filtering from a weighted owner into an output image.
     */
    static void filteringByDirectRule(const WeightedMorphologicalTree<T>& weighted, std::vector<bool>& criterion, ImagePtr<T> imgOutputPtr) {
        filteringByDirectRule(weighted.asView(), criterion, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-min filtering from a criterion into an output image.
     */
    static void filteringByPruningMin(const WeightedTreeView<T>& weighted, std::vector<bool>& criterion, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMinCriterionImpl(weighted, criterion, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-min filtering from a weighted owner into an output image.
     */
    static void filteringByPruningMin(const WeightedMorphologicalTree<T>& weighted, std::vector<bool>& criterion, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMin(weighted.asView(), criterion, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-max filtering from a criterion into an output image.
     */
    static void filteringByPruningMax(const WeightedTreeView<T>& weighted, std::vector<bool>& criterion, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMaxCriterionImpl(weighted, criterion, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-max filtering from a weighted owner into an output image.
     */
    static void filteringByPruningMax(const WeightedMorphologicalTree<T>& weighted, std::vector<bool>& criterion, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMax(weighted.asView(), criterion, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-min filtering from an owned attribute buffer into an output image.
     */
    template <std::floating_point Real>
    static void filteringByPruningMin(const WeightedTreeView<T>& weighted, const std::shared_ptr<Real[]>& attribute, Real threshold, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMin(weighted, attribute.get(), threshold, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-min filtering from a weighted owner and owned attribute buffer.
     */
    template <std::floating_point Real>
    static void filteringByPruningMin(const WeightedMorphologicalTree<T>& weighted, const std::shared_ptr<Real[]>& attribute, Real threshold, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMin(weighted.asView(), attribute.get(), threshold, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-min filtering from a raw attribute buffer into an output image.
     */
    template <std::floating_point Real>
    static void filteringByPruningMin(const WeightedTreeView<T>& weighted, const Real* attribute, Real threshold, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMinAttributeImpl(weighted, attribute, threshold, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-min filtering from a weighted owner and raw attribute buffer.
     */
    template <std::floating_point Real>
    static void filteringByPruningMin(const WeightedMorphologicalTree<T>& weighted, const Real* attribute, Real threshold, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMin(weighted.asView(), attribute, threshold, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-max filtering from an owned attribute buffer into an output image.
     */
    template <std::floating_point Real>
    static void filteringByPruningMax(const WeightedTreeView<T>& weighted, const std::shared_ptr<Real[]>& attribute, Real threshold, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMax(weighted, attribute.get(), threshold, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-max filtering from a weighted owner and owned attribute buffer.
     */
    template <std::floating_point Real>
    static void filteringByPruningMax(const WeightedMorphologicalTree<T>& weighted, const std::shared_ptr<Real[]>& attribute, Real threshold, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMax(weighted.asView(), attribute.get(), threshold, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-max filtering from a raw attribute buffer into an output image.
     */
    template <std::floating_point Real>
    static void filteringByPruningMax(const WeightedTreeView<T>& weighted, const Real* attribute, Real threshold, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMaxAttributeImpl(weighted, attribute, threshold, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-max filtering from a weighted owner and raw attribute buffer.
     */
    template <std::floating_point Real>
    static void filteringByPruningMax(const WeightedMorphologicalTree<T>& weighted, const Real* attribute, Real threshold, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMax(weighted.asView(), attribute, threshold, imgOutputPtr);
    }

    /**
     * @brief Builds an MSER-assisted pruning criterion from an attribute threshold.
     */
    template <std::floating_point Real>
    [[nodiscard]] static std::vector<bool> getAdaptiveCriterion(const WeightedMorphologicalTree<T>& weighted, const std::shared_ptr<Real[]>& attribute, Real threshold, AltitudeDiff<T> delta) {
        return getAdaptiveCriterionImpl(weighted, attribute.get(), threshold, delta);
    }

    /**
     * @brief Builds an MSER-assisted pruning criterion from a raw attribute buffer.
     */
    template <std::floating_point Real>
    [[nodiscard]] static std::vector<bool> getAdaptiveCriterion(const WeightedMorphologicalTree<T>& weighted, const Real* attribute, Real threshold, AltitudeDiff<T> delta) {
        return getAdaptiveCriterionImpl(weighted, attribute, threshold, delta);
    }

    /**
     * @brief Builds an MSER-assisted pruning criterion from an existing criterion mask.
     */
    [[nodiscard]] static std::vector<bool> getAdaptiveCriterion(const WeightedMorphologicalTree<T>& weighted, std::vector<bool>& criterion, AltitudeDiff<T> delta) {
        return getAdaptiveCriterionImpl(weighted, criterion, delta);
    }

    /**
     * @brief Builds a depth-stability pruning criterion from an attribute threshold.
     */
    template <std::floating_point Real>
    [[nodiscard]] static std::vector<bool> getAdaptiveCriterionByDepth(const WeightedMorphologicalTree<T>& weighted, const std::shared_ptr<Real[]>& attribute, Real threshold, int depthDelta) {
        return getAdaptiveCriterionByDepthImpl(weighted.topology(), attribute.get(), threshold, depthDelta);
    }

    /**
     * @brief Builds a depth-stability pruning criterion from a raw attribute buffer.
     */
    template <std::floating_point Real>
    [[nodiscard]] static std::vector<bool> getAdaptiveCriterionByDepth(const WeightedMorphologicalTree<T>& weighted, const Real* attribute, Real threshold, int depthDelta) {
        return getAdaptiveCriterionByDepthImpl(weighted.topology(), attribute, threshold, depthDelta);
    }

    /**
     * @brief Builds a depth-stability pruning criterion from an existing criterion mask.
     */
    [[nodiscard]] static std::vector<bool> getAdaptiveCriterionByDepth(const WeightedMorphologicalTree<T>& weighted, std::vector<bool>& criterion, int depthDelta) {
        return getAdaptiveCriterionByDepthImpl(weighted.topology(), criterion, depthDelta);
    }
};


} // namespace mmcfilters
