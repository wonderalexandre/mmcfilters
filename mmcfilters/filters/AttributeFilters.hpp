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
    MorphologicalTree& tree;
    const AltitudeBuffer* altitude_ = nullptr;

    static AltitudeType altitudeOf(const AltitudeBuffer* altitude, NodeId nodeId) {
        return tree_altitude_ops::getAltitude(altitude, nodeId);
    }

    static AltitudeDiffType residueOf(MorphologicalTree& tree, const AltitudeBuffer* altitude, NodeId nodeId) {
        return tree_altitude_ops::getNodeResidue(tree, altitude, nodeId);
    }

    template <typename TValue>
    static void writeProperParts(MorphologicalTree& tree, NodeId nodeId, TValue* output, TValue value) {
        for (int pixel : tree.getProperParts(nodeId)) {
            output[pixel] = value;
        }
    }

    template <typename TValue>
    static void writeSubtreeProperParts(MorphologicalTree& tree, NodeId nodeId, TValue* output, TValue value) {
        for (NodeId subtreeNodeId : tree.getNodeSubtree(nodeId)) {
            writeProperParts(tree, subtreeNodeId, output, value);
        }
    }

    static void filteringBySubtractiveScoreRuleImpl(MorphologicalTree& tree, const AltitudeBuffer* altitude, std::vector<float>& prob, ImageFloatPtr imgOutputPtr) {
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

    static void filteringBySubtractiveRuleImpl(MorphologicalTree& tree, const AltitudeBuffer* altitude, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr) {
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

    static void filteringByDirectRuleImpl(MorphologicalTree& tree, const AltitudeBuffer* altitude, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr) {
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

    static void filteringByPruningMinCriterionImpl(MorphologicalTree& tree, const AltitudeBuffer* altitude, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr) {
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

    static void filteringByPruningMaxCriterionImpl(MorphologicalTree& tree, const AltitudeBuffer* altitude, std::vector<bool>& keepCriterion, ImageUInt8Ptr imgOutputPtr) {
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

    static void filteringByPruningMinAttributeImpl(MorphologicalTree& tree, const AltitudeBuffer* altitude, const float* attribute, float threshold, ImageUInt8Ptr imgOutputPtr) {
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

    static void filteringByPruningMaxAttributeImpl(MorphologicalTree& tree, const AltitudeBuffer* altitude, const float* attribute, float threshold, ImageUInt8Ptr imgOutputPtr) {
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

public:
    explicit AttributeFilters(MorphologicalTree& tree)
        : tree{tree}, altitude_{nullptr} {}

    explicit AttributeFilters(WeightedMorphologicalTree& weighted)
        : tree{weighted.tree}, altitude_{&weighted.altitude} {}

    ~AttributeFilters() = default;

    std::vector<bool> getAdaptiveCriterion(std::vector<bool>& criterion, int delta) {
        return AttributeFilters::getAdaptiveCriterion(this->tree, criterion, delta);
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

    static void filteringBySubtractiveScoreRule(MorphologicalTree& tree, std::vector<float>& prob, ImageFloatPtr imgOutputPtr) {
        filteringBySubtractiveScoreRuleImpl(tree, nullptr, prob, imgOutputPtr);
    }

    static void filteringBySubtractiveScoreRule(WeightedMorphologicalTree& weighted, std::vector<float>& prob, ImageFloatPtr imgOutputPtr) {
        filteringBySubtractiveScoreRuleImpl(weighted.tree, &weighted.altitude, prob, imgOutputPtr);
    }

    static void filteringBySubtractiveRule(MorphologicalTree& tree, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr) {
        filteringBySubtractiveRuleImpl(tree, nullptr, criterion, imgOutputPtr);
    }

    static void filteringBySubtractiveRule(WeightedMorphologicalTree& weighted, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr) {
        filteringBySubtractiveRuleImpl(weighted.tree, &weighted.altitude, criterion, imgOutputPtr);
    }

    static void filteringByDirectRule(MorphologicalTree& tree, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr) {
        filteringByDirectRuleImpl(tree, nullptr, criterion, imgOutputPtr);
    }

    static void filteringByDirectRule(WeightedMorphologicalTree& weighted, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr) {
        filteringByDirectRuleImpl(weighted.tree, &weighted.altitude, criterion, imgOutputPtr);
    }

    static void filteringByPruningMin(MorphologicalTree& tree, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr) {
        filteringByPruningMinCriterionImpl(tree, nullptr, criterion, imgOutputPtr);
    }

    static void filteringByPruningMin(WeightedMorphologicalTree& weighted, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr) {
        filteringByPruningMinCriterionImpl(weighted.tree, &weighted.altitude, criterion, imgOutputPtr);
    }

    static void filteringByPruningMax(MorphologicalTree& tree, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr) {
        filteringByPruningMaxCriterionImpl(tree, nullptr, criterion, imgOutputPtr);
    }

    static void filteringByPruningMax(WeightedMorphologicalTree& weighted, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr) {
        filteringByPruningMaxCriterionImpl(weighted.tree, &weighted.altitude, criterion, imgOutputPtr);
    }

    static void filteringByPruningMin(MorphologicalTree& tree, const std::shared_ptr<float[]>& attribute, float threshold, ImageUInt8Ptr imgOutputPtr) {
        filteringByPruningMin(tree, attribute.get(), threshold, imgOutputPtr);
    }

    static void filteringByPruningMin(WeightedMorphologicalTree& weighted, const std::shared_ptr<float[]>& attribute, float threshold, ImageUInt8Ptr imgOutputPtr) {
        filteringByPruningMin(weighted, attribute.get(), threshold, imgOutputPtr);
    }

    static void filteringByPruningMin(MorphologicalTree& tree, const float* attribute, float threshold, ImageUInt8Ptr imgOutputPtr) {
        filteringByPruningMinAttributeImpl(tree, nullptr, attribute, threshold, imgOutputPtr);
    }

    static void filteringByPruningMin(WeightedMorphologicalTree& weighted, const float* attribute, float threshold, ImageUInt8Ptr imgOutputPtr) {
        filteringByPruningMinAttributeImpl(weighted.tree, &weighted.altitude, attribute, threshold, imgOutputPtr);
    }

    static void filteringByPruningMax(MorphologicalTree& tree, const std::shared_ptr<float[]>& attribute, float threshold, ImageUInt8Ptr imgOutputPtr) {
        filteringByPruningMax(tree, attribute.get(), threshold, imgOutputPtr);
    }

    static void filteringByPruningMax(WeightedMorphologicalTree& weighted, const std::shared_ptr<float[]>& attribute, float threshold, ImageUInt8Ptr imgOutputPtr) {
        filteringByPruningMax(weighted, attribute.get(), threshold, imgOutputPtr);
    }

    static void filteringByPruningMax(MorphologicalTree& tree, const float* attribute, float threshold, ImageUInt8Ptr imgOutputPtr) {
        filteringByPruningMaxAttributeImpl(tree, nullptr, attribute, threshold, imgOutputPtr);
    }

    static void filteringByPruningMax(WeightedMorphologicalTree& weighted, const float* attribute, float threshold, ImageUInt8Ptr imgOutputPtr) {
        filteringByPruningMaxAttributeImpl(weighted.tree, &weighted.altitude, attribute, threshold, imgOutputPtr);
    }

    static std::vector<bool> getAdaptiveCriterion(MorphologicalTree& tree, const std::shared_ptr<float[]>& attribute, float threshold, int delta) {
        return getAdaptiveCriterion(tree, attribute.get(), threshold, delta);
    }

    static std::vector<bool> getAdaptiveCriterion(WeightedMorphologicalTree& weighted, const std::shared_ptr<float[]>& attribute, float threshold, int delta) {
        return getAdaptiveCriterion(weighted.tree, attribute.get(), threshold, delta);
    }

    static std::vector<bool> getAdaptiveCriterion(MorphologicalTree& tree, const float* attribute, float threshold, int delta) {
        assert(attribute != nullptr);

        ComputerMSER mser(tree);
        std::vector<uint8_t> isMSER = mser.computeMSER(delta);

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

    static std::vector<bool> getAdaptiveCriterion(WeightedMorphologicalTree& weighted, const float* attribute, float threshold, int delta) {
        return getAdaptiveCriterion(weighted.tree, attribute, threshold, delta);
    }

    static std::vector<bool> getAdaptiveCriterion(MorphologicalTree& tree, std::vector<bool>& criterion, int delta) {
        ComputerMSER mser(tree);
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

    static std::vector<bool> getAdaptiveCriterion(WeightedMorphologicalTree& weighted, std::vector<bool>& criterion, int delta) {
        return getAdaptiveCriterion(weighted.tree, criterion, delta);
    }
};

} // namespace mmcfilters
