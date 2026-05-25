#pragma once

#include "../utils/Image.hpp"
#include "../utils/Common.hpp"
#include "../trees/TreeAltitudeAlgorithms.hpp"
#include "../trees/WeightedMorphologicalTree.hpp"
#include "../trees/WeightedTreeView.hpp"
#include "../trees/detail/TreeTraversalDetail.hpp"
#include "ComputerMSER.hpp"

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

    static void requireAttributePointer(const float* attribute, const char* context) {
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

    template <typename TValue>
    static void writeProperParts(const MorphologicalTree& tree, NodeId nodeId, TValue* output, TValue value) {
        for (int pixel : tree.getProperParts(nodeId)) {
            output[pixel] = value;
        }
    }

    template <typename TValue>
    static void writeSubtreeProperParts(const MorphologicalTree& tree, NodeId nodeId, TValue* output, TValue value) {
        for (NodeId subtreeNodeId : tree.getNodeSubtree(nodeId)) {
            writeProperParts(tree, subtreeNodeId, output, value);
        }
    }

    static void filteringBySubtractiveScoreRuleImpl(AltitudeView view, std::vector<float>& prob, ImageFloatPtr imgOutputPtr) {
        view.requireTopologyUnchanged("AttributeFilters::filteringBySubtractiveScoreRule");
        const MorphologicalTree& tree = view.topology();
        requireScoreSize(tree, prob, "AttributeFilters::filteringBySubtractiveScoreRule");
        requireOutputImage(tree, imgOutputPtr, "AttributeFilters::filteringBySubtractiveScoreRule");
        std::unique_ptr<float[]> mapLevel(new float[tree.getNumInternalNodeSlots()]);

        const NodeId rootNodeId = tree.getRoot();
        mapLevel[rootNodeId] = static_cast<float>(altitudeOf(view, rootNodeId));

        for (NodeId nodeId : tree.getAliveNodeIds()) {
            if (!tree.isRoot(nodeId)) {
                const NodeId parentNodeId = tree.getNodeParent(nodeId);
                mapLevel[nodeId] = mapLevel[parentNodeId] + (static_cast<float>(residueOf(view, nodeId)) * prob[nodeId]);
            }
        }

        auto imgOutput = imgOutputPtr->rawData();
        for (NodeId nodeId : tree.getAliveNodeIds()) {
            writeProperParts(tree, nodeId, imgOutput, mapLevel[nodeId]);
        }
    }

    static void filteringBySubtractiveRuleImpl(AltitudeView view, std::vector<bool>& criterion, ImagePtr<T> imgOutputPtr) {
        const char* context = "AttributeFilters::filteringBySubtractiveRule";
        view.requireTopologyUnchanged(context);
        const MorphologicalTree& tree = view.topology();
        requireCriterionSize(tree, criterion, context);
        requireOutputImage(tree, imgOutputPtr, context);
        std::vector<AltitudeDiff<T>> mapLevel(static_cast<std::size_t>(tree.getNumInternalNodeSlots()));
        const NodeId rootNodeId = tree.getRoot();
        mapLevel[rootNodeId] = static_cast<AltitudeDiff<T>>(altitudeOf(view, rootNodeId));

        for (NodeId nodeId : tree.getAliveNodeIds()) {
            if (!tree.isRoot(nodeId)) {
                const NodeId parentNodeId = tree.getNodeParent(nodeId);
                mapLevel[nodeId] = criterion[nodeId]
                    ? static_cast<AltitudeDiff<T>>(mapLevel[parentNodeId] + residueOf(view, nodeId))
                    : mapLevel[parentNodeId];
            }
        }

        auto imgOutput = imgOutputPtr->rawData();
        for (NodeId nodeId : tree.getAliveNodeIds()) {
            writeProperParts(tree, nodeId, imgOutput, static_cast<T>(mapLevel[nodeId]));
        }
    }

    static void filteringByDirectRuleImpl(AltitudeView view, std::vector<bool>& criterion, ImagePtr<T> imgOutputPtr) {
        const char* context = "AttributeFilters::filteringByDirectRule";
        view.requireTopologyUnchanged(context);
        const MorphologicalTree& tree = view.topology();
        requireCriterionSize(tree, criterion, context);
        requireOutputImage(tree, imgOutputPtr, context);
        std::vector<T> mapLevel(static_cast<std::size_t>(tree.getNumInternalNodeSlots()));
        const NodeId rootNodeId = tree.getRoot();
        mapLevel[rootNodeId] = altitudeOf(view, rootNodeId);

        for (NodeId nodeId : tree.getAliveNodeIds()) {
            if (!tree.isRoot(nodeId)) {
                const NodeId parentNodeId = tree.getNodeParent(nodeId);
                mapLevel[nodeId] = criterion[nodeId] ? altitudeOf(view, nodeId) : mapLevel[parentNodeId];
            }
        }

        auto imgOutput = imgOutputPtr->rawData();
        for (NodeId nodeId : tree.getAliveNodeIds()) {
            writeProperParts(tree, nodeId, imgOutput, mapLevel[nodeId]);
        }
    }

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
                    writeSubtreeProperParts(tree, childNodeId, imgOutput, altitudeOf(view, childNodeId));
                }
            }
        }
    }

    static void filteringByPruningMaxCriterionImpl(AltitudeView view, std::vector<bool>& keepCriterion, ImagePtr<T> imgOutputPtr) {
        const char* context = "AttributeFilters::filteringByPruningMax";
        view.requireTopologyUnchanged(context);
        const MorphologicalTree& tree = view.topology();
        requireCriterionSize(tree, keepCriterion, context);
        requireOutputImage(tree, imgOutputPtr, context);
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

    static void filteringByPruningMinAttributeImpl(AltitudeView view, const float* attribute, float threshold, ImagePtr<T> imgOutputPtr) {
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

    static void filteringByPruningMaxAttributeImpl(AltitudeView view, const float* attribute, float threshold, ImagePtr<T> imgOutputPtr) {
        const char* context = "AttributeFilters::filteringByPruningMax";
        view.requireTopologyUnchanged(context);
        const MorphologicalTree& tree = view.topology();
        requireAttributePointer(attribute, context);
        requireOutputImage(tree, imgOutputPtr, context);
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

    static std::vector<bool> getAdaptiveCriterionImpl(const WeightedMorphologicalTree<T>& weighted, const float* attribute, float threshold, AltitudeDiff<T> delta) {
        requireAttributePointer(attribute, "AttributeFilters::getAdaptiveCriterion");

        const MorphologicalTree& tree = weighted.topology();
        ComputerMSER<T> mser(weighted);
        std::vector<uint8_t> isMSER = mser.computeMSER(delta);
        (void)isMSER;

        std::vector<float> stability = mser.getStabilities();
        std::vector<bool> isPruned(tree.getNumInternalNodeSlots(), false);
        for (NodeId nodeId : tree.getAliveNodeIds()) {
            if (attribute[nodeId] < threshold) {
                if (std::isnan(stability[nodeId])) {
                    isPruned[nodeId] = true;
                } else {
                    const float max = stability[nodeId];
                    const NodeId indexDescMaxStability = mser.descendantWithMaxStability(nodeId);
                    const NodeId indexAscMaxStability = mser.ascendantWithMaxStability(nodeId);
                    const float maxDesc = stability[indexDescMaxStability];
                    const float maxAnc = stability[indexAscMaxStability];

                    if (max >= maxDesc && max >= maxAnc) {
                        isPruned[nodeId] = true;
                    } else if (maxDesc >= max && maxDesc >= maxAnc) {
                        isPruned[indexDescMaxStability] = true;
                    } else {
                        isPruned[indexAscMaxStability] = true;
                    }
                }
            }
        }
        return isPruned;
    }

    static std::vector<bool> getAdaptiveCriterionImpl(const WeightedMorphologicalTree<T>& weighted, std::vector<bool>& criterion, AltitudeDiff<T> delta) {
        const MorphologicalTree& tree = weighted.topology();
        requireCriterionSize(tree, criterion, "AttributeFilters::getAdaptiveCriterion");
        ComputerMSER<T> mser(weighted);
        std::vector<uint8_t> isMSER = mser.computeMSER(delta);
        (void)isMSER;

        std::vector<float> stability = mser.getStabilities();
        std::vector<bool> isPruned(tree.getNumInternalNodeSlots(), false);
        for (NodeId nodeId : tree.getAliveNodeIds()) {
            if (!criterion[nodeId]) {
                if (std::isnan(stability[nodeId])) {
                    isPruned[nodeId] = true;
                } else {
                    const float max = stability[nodeId];
                    const NodeId indexDescMaxStability = mser.descendantWithMaxStability(nodeId);
                    const NodeId indexAscMaxStability = mser.ascendantWithMaxStability(nodeId);
                    const float maxDesc = stability[indexDescMaxStability];
                    const float maxAnc = stability[indexAscMaxStability];

                    if (max >= maxDesc && max >= maxAnc) {
                        isPruned[nodeId] = true;
                    } else if (maxDesc >= max && maxDesc >= maxAnc) {
                        isPruned[indexDescMaxStability] = true;
                    } else {
                        isPruned[indexAscMaxStability] = true;
                    }
                }
            }
        }
        return isPruned;
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
     * @brief Applies pruning-min filtering from an attribute buffer.
     *
     * Nodes with attribute values above `threshold` remain traversable; rejected
     * subtrees are reconstructed at the ancestor level selected by the pruning-min
     * rule.
     */
    [[nodiscard]] ImagePtr<T> filteringByPruningMin(const std::shared_ptr<float[]>& attr, float threshold) {
        return filteringByPruningMin(attr.get(), threshold);
    }

    /**
     * @brief Applies pruning-min filtering from a raw internal-node attribute buffer.
     */
    [[nodiscard]] ImagePtr<T> filteringByPruningMin(const float* attr, float threshold) {
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
    [[nodiscard]] ImagePtr<T> filteringByPruningMax(const std::shared_ptr<float[]>& attr, float threshold) {
        return filteringByPruningMax(attr.get(), threshold);
    }

    /**
     * @brief Applies pruning-max filtering from a raw internal-node attribute buffer.
     */
    [[nodiscard]] ImagePtr<T> filteringByPruningMax(const float* attr, float threshold) {
        requireStableTree("AttributeFilters::filteringByPruningMax");
        assert(attr != nullptr);
        ImagePtr<T> imgOutput = Image<T>::create(this->tree.getNumRowsOfImage(), this->tree.getNumColsOfImage());
        filteringByPruningMaxAttributeImpl(view(), attr, threshold, imgOutput);
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
    static void filteringByPruningMin(const WeightedTreeView<T>& weighted, const std::shared_ptr<float[]>& attribute, float threshold, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMin(weighted, attribute.get(), threshold, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-min filtering from a weighted owner and owned attribute buffer.
     */
    static void filteringByPruningMin(const WeightedMorphologicalTree<T>& weighted, const std::shared_ptr<float[]>& attribute, float threshold, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMin(weighted.asView(), attribute.get(), threshold, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-min filtering from a raw attribute buffer into an output image.
     */
    static void filteringByPruningMin(const WeightedTreeView<T>& weighted, const float* attribute, float threshold, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMinAttributeImpl(weighted, attribute, threshold, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-min filtering from a weighted owner and raw attribute buffer.
     */
    static void filteringByPruningMin(const WeightedMorphologicalTree<T>& weighted, const float* attribute, float threshold, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMin(weighted.asView(), attribute, threshold, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-max filtering from an owned attribute buffer into an output image.
     */
    static void filteringByPruningMax(const WeightedTreeView<T>& weighted, const std::shared_ptr<float[]>& attribute, float threshold, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMax(weighted, attribute.get(), threshold, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-max filtering from a weighted owner and owned attribute buffer.
     */
    static void filteringByPruningMax(const WeightedMorphologicalTree<T>& weighted, const std::shared_ptr<float[]>& attribute, float threshold, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMax(weighted.asView(), attribute.get(), threshold, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-max filtering from a raw attribute buffer into an output image.
     */
    static void filteringByPruningMax(const WeightedTreeView<T>& weighted, const float* attribute, float threshold, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMaxAttributeImpl(weighted, attribute, threshold, imgOutputPtr);
    }

    /**
     * @brief Writes pruning-max filtering from a weighted owner and raw attribute buffer.
     */
    static void filteringByPruningMax(const WeightedMorphologicalTree<T>& weighted, const float* attribute, float threshold, ImagePtr<T> imgOutputPtr) {
        filteringByPruningMax(weighted.asView(), attribute, threshold, imgOutputPtr);
    }

    /**
     * @brief Builds an MSER-assisted pruning criterion from an attribute threshold.
     */
    [[nodiscard]] static std::vector<bool> getAdaptiveCriterion(const WeightedMorphologicalTree<T>& weighted, const std::shared_ptr<float[]>& attribute, float threshold, AltitudeDiff<T> delta) {
        return getAdaptiveCriterionImpl(weighted, attribute.get(), threshold, delta);
    }

    /**
     * @brief Builds an MSER-assisted pruning criterion from a raw attribute buffer.
     */
    [[nodiscard]] static std::vector<bool> getAdaptiveCriterion(const WeightedMorphologicalTree<T>& weighted, const float* attribute, float threshold, AltitudeDiff<T> delta) {
        return getAdaptiveCriterionImpl(weighted, attribute, threshold, delta);
    }

    /**
     * @brief Builds an MSER-assisted pruning criterion from an existing criterion mask.
     */
    [[nodiscard]] static std::vector<bool> getAdaptiveCriterion(const WeightedMorphologicalTree<T>& weighted, std::vector<bool>& criterion, AltitudeDiff<T> delta) {
        return getAdaptiveCriterionImpl(weighted, criterion, delta);
    }
};


} // namespace mmcfilters
