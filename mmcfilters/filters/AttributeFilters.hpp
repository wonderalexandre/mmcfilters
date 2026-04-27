#pragma once

#include "../utils/Common.hpp"
#include "../trees/WeightedMorphologicalTree.hpp"
#include "../attributes/AttributeComputedIncrementally.hpp"
#include "../attributes/ComputerMSER.hpp"

#include <cassert>

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
class AttributeFilters {
protected:
    const MorphologicalTree& tree;
    const AltitudeBuffer* altitude_ = nullptr;

    static AltitudeType altitudeOf(const AltitudeBuffer* altitude, NodeId nodeId) {
        return WeightedMorphologicalTree::getAltitude(altitude, nodeId);
    }

    static AltitudeDiffType residueOf(const MorphologicalTree& tree, const AltitudeBuffer* altitude, NodeId nodeId) {
        return WeightedMorphologicalTree::getNodeResidue(tree, altitude, nodeId);
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

    static void filteringBySubtractiveScoreRuleImpl(const MorphologicalTree& tree, const AltitudeBuffer* altitude, std::vector<float>& prob, ImageFloatPtr imgOutputPtr) {
        std::unique_ptr<float[]> mapLevel(new float[tree.getNumInternalNodeSlots()]);

        const NodeId rootNodeId = tree.getRoot();
        mapLevel[rootNodeId] = static_cast<float>(altitudeOf(altitude, rootNodeId));

        for (NodeId nodeId : tree.getAliveNodeIds()) {
            if (!tree.isRoot(nodeId)) {
                const NodeId parentNodeId = tree.getNodeParent(nodeId);
                mapLevel[nodeId] = mapLevel[parentNodeId] + (static_cast<float>(residueOf(tree, altitude, nodeId)) * prob[nodeId]);
            }
        }

        auto imgOutput = imgOutputPtr->rawData();
        for (NodeId nodeId : tree.getAliveNodeIds()) {
            writeProperParts(tree, nodeId, imgOutput, mapLevel[nodeId]);
        }
    }

    static void filteringBySubtractiveRuleImpl(const MorphologicalTree& tree, const AltitudeBuffer* altitude, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr) {
        std::unique_ptr<AltitudeType[]> mapLevel(new AltitudeType[tree.getNumInternalNodeSlots()]);
        const NodeId rootNodeId = tree.getRoot();
        mapLevel[rootNodeId] = altitudeOf(altitude, rootNodeId);

        for (NodeId nodeId : tree.getAliveNodeIds()) {
            if (!tree.isRoot(nodeId)) {
                const NodeId parentNodeId = tree.getNodeParent(nodeId);
                mapLevel[nodeId] = criterion[nodeId]
                    ? static_cast<AltitudeType>(mapLevel[parentNodeId] + residueOf(tree, altitude, nodeId))
                    : mapLevel[parentNodeId];
            }
        }

        auto imgOutput = imgOutputPtr->rawData();
        for (NodeId nodeId : tree.getAliveNodeIds()) {
            writeProperParts(tree, nodeId, imgOutput, static_cast<uint8_t>(mapLevel[nodeId]));
        }
    }

    static void filteringByDirectRuleImpl(const MorphologicalTree& tree, const AltitudeBuffer* altitude, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr) {
        std::unique_ptr<AltitudeType[]> mapLevel(new AltitudeType[tree.getNumInternalNodeSlots()]);
        const NodeId rootNodeId = tree.getRoot();
        mapLevel[rootNodeId] = altitudeOf(altitude, rootNodeId);

        for (NodeId nodeId : tree.getAliveNodeIds()) {
            if (!tree.isRoot(nodeId)) {
                const NodeId parentNodeId = tree.getNodeParent(nodeId);
                mapLevel[nodeId] = criterion[nodeId] ? altitudeOf(altitude, nodeId) : mapLevel[parentNodeId];
            }
        }

        auto imgOutput = imgOutputPtr->rawData();
        for (NodeId nodeId : tree.getAliveNodeIds()) {
            writeProperParts(tree, nodeId, imgOutput, static_cast<uint8_t>(mapLevel[nodeId]));
        }
    }

    static void filteringByPruningMinCriterionImpl(const MorphologicalTree& tree, const AltitudeBuffer* altitude, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr) {
        std::stack<NodeId> stack;
        stack.push(tree.getRoot());
        auto imgOutput = imgOutputPtr->rawData();

        while (!stack.empty()) {
            const NodeId nodeId = stack.top();
            stack.pop();
            writeProperParts(tree, nodeId, imgOutput, static_cast<uint8_t>(altitudeOf(altitude, nodeId)));
            for (NodeId childNodeId : tree.getChildren(nodeId)) {
                if (criterion[childNodeId]) {
                    stack.push(childNodeId);
                } else {
                    writeSubtreeProperParts(tree, childNodeId, imgOutput, static_cast<uint8_t>(altitudeOf(altitude, childNodeId)));
                }
            }
        }
    }

    static void filteringByPruningMaxCriterionImpl(const MorphologicalTree& tree, const AltitudeBuffer* altitude, std::vector<bool>& keepCriterion, ImageUInt8Ptr imgOutputPtr) {
        std::vector<uint8_t> criterion(tree.getNumInternalNodeSlots(), false);
        AttributeComputedIncrementally::traversePostOrder(
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
            writeProperParts(tree, nodeId, imgOutput, static_cast<uint8_t>(altitudeOf(altitude, nodeId)));
            for (NodeId childNodeId : tree.getChildren(nodeId)) {
                if (!criterion[childNodeId]) {
                    stack.push(childNodeId);
                } else {
                    writeSubtreeProperParts(tree, childNodeId, imgOutput, static_cast<uint8_t>(altitudeOf(altitude, childNodeId)));
                }
            }
        }
    }

    static void filteringByPruningMinAttributeImpl(const MorphologicalTree& tree, const AltitudeBuffer* altitude, const float* attribute, float threshold, ImageUInt8Ptr imgOutputPtr) {
        assert(attribute != nullptr);
        auto imgOutput = imgOutputPtr->rawData();
        std::stack<NodeId> stack;
        stack.push(tree.getRoot());
        while (!stack.empty()) {
            const NodeId nodeId = stack.top();
            stack.pop();
            writeProperParts(tree, nodeId, imgOutput, static_cast<uint8_t>(altitudeOf(altitude, nodeId)));
            for (NodeId childNodeId : tree.getChildren(nodeId)) {
                if (attribute[childNodeId] > threshold) {
                    stack.push(childNodeId);
                } else {
                    writeSubtreeProperParts(tree, childNodeId, imgOutput, static_cast<uint8_t>(altitudeOf(altitude, nodeId)));
                }
            }
        }
    }

    static void filteringByPruningMaxAttributeImpl(const MorphologicalTree& tree, const AltitudeBuffer* altitude, const float* attribute, float threshold, ImageUInt8Ptr imgOutputPtr) {
        assert(attribute != nullptr);
        std::vector<uint8_t> criterion(tree.getNumInternalNodeSlots(), false);
        AttributeComputedIncrementally::traversePostOrder(
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
            writeProperParts(tree, nodeId, imgOutput, static_cast<uint8_t>(altitudeOf(altitude, nodeId)));
            for (NodeId childNodeId : tree.getChildren(nodeId)) {
                if (!criterion[childNodeId]) {
                    stack.push(childNodeId);
                } else {
                    writeSubtreeProperParts(tree, childNodeId, imgOutput, static_cast<uint8_t>(altitudeOf(altitude, nodeId)));
                }
            }
        }
    }

    static std::vector<bool> getAdaptiveCriterionImpl(const MorphologicalTree& tree, const AltitudeBuffer* altitude, const float* attribute, float threshold, int delta) {
        assert(attribute != nullptr);

        ComputerMSER mser(tree, altitude);
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

    static std::vector<bool> getAdaptiveCriterionImpl(const MorphologicalTree& tree, const AltitudeBuffer* altitude, std::vector<bool>& criterion, int delta) {
        ComputerMSER mser(tree, altitude);
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

public:
    explicit AttributeFilters(const MorphologicalTree& tree)
        : tree{tree}, altitude_{nullptr} {}

    explicit AttributeFilters(const WeightedMorphologicalTree& weighted)
        : tree{weighted.tree_}, altitude_{&weighted.altitude_} {}

    ~AttributeFilters() = default;

    std::vector<bool> getAdaptiveCriterion(std::vector<bool>& criterion, int delta) {
        return AttributeFilters::getAdaptiveCriterionImpl(this->tree, this->altitude_, criterion, delta);
    }

    ImageUInt8Ptr filteringByPruningMin(const std::shared_ptr<float[]>& attr, float threshold) {
        return filteringByPruningMin(attr.get(), threshold);
    }

    ImageUInt8Ptr filteringByPruningMin(const float* attr, float threshold) {
        assert(attr != nullptr);
        ImageUInt8Ptr imgOutput = ImageUInt8::create(this->tree.getNumRowsOfImage(), this->tree.getNumColsOfImage());
        filteringByPruningMinAttributeImpl(tree, altitude_, attr, threshold, imgOutput);
        return imgOutput;
    }

    ImageUInt8Ptr filteringByPruningMax(const std::shared_ptr<float[]>& attr, float threshold) {
        return filteringByPruningMax(attr.get(), threshold);
    }

    ImageUInt8Ptr filteringByPruningMax(const float* attr, float threshold) {
        assert(attr != nullptr);
        ImageUInt8Ptr imgOutput = ImageUInt8::create(this->tree.getNumRowsOfImage(), this->tree.getNumColsOfImage());
        filteringByPruningMaxAttributeImpl(tree, altitude_, attr, threshold, imgOutput);
        return imgOutput;
    }

    ImageUInt8Ptr filteringByPruningMin(std::vector<bool>& criterion) {
        ImageUInt8Ptr imgOutput = ImageUInt8::create(this->tree.getNumRowsOfImage(), this->tree.getNumColsOfImage());
        filteringByPruningMinCriterionImpl(tree, altitude_, criterion, imgOutput);
        return imgOutput;
    }

    ImageUInt8Ptr filteringByPruningMax(std::vector<bool>& criterion) {
        ImageUInt8Ptr imgOutput = ImageUInt8::create(this->tree.getNumRowsOfImage(), this->tree.getNumColsOfImage());
        filteringByPruningMaxCriterionImpl(tree, altitude_, criterion, imgOutput);
        return imgOutput;
    }

    ImageUInt8Ptr filteringByDirectRule(std::vector<bool>& criterion) {
        ImageUInt8Ptr imgOutput = ImageUInt8::create(this->tree.getNumRowsOfImage(), this->tree.getNumColsOfImage());
        filteringByDirectRuleImpl(tree, altitude_, criterion, imgOutput);
        return imgOutput;
    }

    ImageUInt8Ptr filteringBySubtractiveRule(std::vector<bool>& criterion) {
        ImageUInt8Ptr imgOutput = ImageUInt8::create(this->tree.getNumRowsOfImage(), this->tree.getNumColsOfImage());
        filteringBySubtractiveRuleImpl(tree, altitude_, criterion, imgOutput);
        return imgOutput;
    }

    ImageFloatPtr filteringBySubtractiveScoreRule(std::vector<float>& prob) {
        ImageFloatPtr imgOutput = ImageFloat::create(this->tree.getNumRowsOfImage(), this->tree.getNumColsOfImage());
        filteringBySubtractiveScoreRuleImpl(tree, altitude_, prob, imgOutput);
        return imgOutput;
    }

    static void filteringBySubtractiveScoreRule(const MorphologicalTree& tree, std::vector<float>& prob, ImageFloatPtr imgOutputPtr) {
        filteringBySubtractiveScoreRuleImpl(tree, nullptr, prob, imgOutputPtr);
    }

    static void filteringBySubtractiveScoreRule(const WeightedMorphologicalTree& weighted, std::vector<float>& prob, ImageFloatPtr imgOutputPtr) {
        filteringBySubtractiveScoreRuleImpl(weighted.tree_, &weighted.altitude_, prob, imgOutputPtr);
    }

    static void filteringBySubtractiveRule(const MorphologicalTree& tree, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr) {
        filteringBySubtractiveRuleImpl(tree, nullptr, criterion, imgOutputPtr);
    }

    static void filteringBySubtractiveRule(const WeightedMorphologicalTree& weighted, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr) {
        filteringBySubtractiveRuleImpl(weighted.tree_, &weighted.altitude_, criterion, imgOutputPtr);
    }

    static void filteringByDirectRule(const MorphologicalTree& tree, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr) {
        filteringByDirectRuleImpl(tree, nullptr, criterion, imgOutputPtr);
    }

    static void filteringByDirectRule(const WeightedMorphologicalTree& weighted, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr) {
        filteringByDirectRuleImpl(weighted.tree_, &weighted.altitude_, criterion, imgOutputPtr);
    }

    static void filteringByPruningMin(const MorphologicalTree& tree, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr) {
        filteringByPruningMinCriterionImpl(tree, nullptr, criterion, imgOutputPtr);
    }

    static void filteringByPruningMin(const WeightedMorphologicalTree& weighted, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr) {
        filteringByPruningMinCriterionImpl(weighted.tree_, &weighted.altitude_, criterion, imgOutputPtr);
    }

    static void filteringByPruningMax(const MorphologicalTree& tree, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr) {
        filteringByPruningMaxCriterionImpl(tree, nullptr, criterion, imgOutputPtr);
    }

    static void filteringByPruningMax(const WeightedMorphologicalTree& weighted, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr) {
        filteringByPruningMaxCriterionImpl(weighted.tree_, &weighted.altitude_, criterion, imgOutputPtr);
    }

    static void filteringByPruningMin(const MorphologicalTree& tree, const std::shared_ptr<float[]>& attribute, float threshold, ImageUInt8Ptr imgOutputPtr) {
        filteringByPruningMin(tree, attribute.get(), threshold, imgOutputPtr);
    }

    static void filteringByPruningMin(const WeightedMorphologicalTree& weighted, const std::shared_ptr<float[]>& attribute, float threshold, ImageUInt8Ptr imgOutputPtr) {
        filteringByPruningMin(weighted, attribute.get(), threshold, imgOutputPtr);
    }

    static void filteringByPruningMin(const MorphologicalTree& tree, const float* attribute, float threshold, ImageUInt8Ptr imgOutputPtr) {
        filteringByPruningMinAttributeImpl(tree, nullptr, attribute, threshold, imgOutputPtr);
    }

    static void filteringByPruningMin(const WeightedMorphologicalTree& weighted, const float* attribute, float threshold, ImageUInt8Ptr imgOutputPtr) {
        filteringByPruningMinAttributeImpl(weighted.tree_, &weighted.altitude_, attribute, threshold, imgOutputPtr);
    }

    static void filteringByPruningMax(const MorphologicalTree& tree, const std::shared_ptr<float[]>& attribute, float threshold, ImageUInt8Ptr imgOutputPtr) {
        filteringByPruningMax(tree, attribute.get(), threshold, imgOutputPtr);
    }

    static void filteringByPruningMax(const WeightedMorphologicalTree& weighted, const std::shared_ptr<float[]>& attribute, float threshold, ImageUInt8Ptr imgOutputPtr) {
        filteringByPruningMax(weighted, attribute.get(), threshold, imgOutputPtr);
    }

    static void filteringByPruningMax(const MorphologicalTree& tree, const float* attribute, float threshold, ImageUInt8Ptr imgOutputPtr) {
        filteringByPruningMaxAttributeImpl(tree, nullptr, attribute, threshold, imgOutputPtr);
    }

    static void filteringByPruningMax(const WeightedMorphologicalTree& weighted, const float* attribute, float threshold, ImageUInt8Ptr imgOutputPtr) {
        filteringByPruningMaxAttributeImpl(weighted.tree_, &weighted.altitude_, attribute, threshold, imgOutputPtr);
    }

    static std::vector<bool> getAdaptiveCriterion(const MorphologicalTree& tree, const std::shared_ptr<float[]>& attribute, float threshold, int delta) {
        return getAdaptiveCriterion(tree, attribute.get(), threshold, delta);
    }

    static std::vector<bool> getAdaptiveCriterion(const WeightedMorphologicalTree& weighted, const std::shared_ptr<float[]>& attribute, float threshold, int delta) {
        return getAdaptiveCriterionImpl(weighted.tree_, &weighted.altitude_, attribute.get(), threshold, delta);
    }

    static std::vector<bool> getAdaptiveCriterion(const MorphologicalTree& tree, const float* attribute, float threshold, int delta) {
        return getAdaptiveCriterionImpl(tree, nullptr, attribute, threshold, delta);
    }

    static std::vector<bool> getAdaptiveCriterion(const WeightedMorphologicalTree& weighted, const float* attribute, float threshold, int delta) {
        return getAdaptiveCriterionImpl(weighted.tree_, &weighted.altitude_, attribute, threshold, delta);
    }

    static std::vector<bool> getAdaptiveCriterion(const MorphologicalTree& tree, std::vector<bool>& criterion, int delta) {
        return getAdaptiveCriterionImpl(tree, nullptr, criterion, delta);
    }

    static std::vector<bool> getAdaptiveCriterion(const WeightedMorphologicalTree& weighted, std::vector<bool>& criterion, int delta) {
        return getAdaptiveCriterionImpl(weighted.tree_, &weighted.altitude_, criterion, delta);
    }
};

} // namespace mmcfilters
