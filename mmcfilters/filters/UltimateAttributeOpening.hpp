#pragma once

#include "../utils/Common.hpp"
#include "../trees/TreeAltitudeOps.hpp"
#include "../trees/WeightedMorphologicalTree.hpp"
#include "../attributes/ComputerMSER.hpp"
#include "../attributes/AttributeComputedIncrementally.hpp"

#include <cassert>

namespace mmcfilters {

/**
 * @brief Computes an Ultimate Attribute Opening by accumulating maximal contrasts.
 */
class UltimateAttributeOpening {
protected:
    int maxCriterion;
    const float* attrs_increasing = nullptr;
    std::shared_ptr<float[]> ownedAttrsIncreasing_;
    MorphologicalTree& tree;
    const AltitudeBuffer* altitude_ = nullptr;
    std::vector<uint8_t> maxContrastLUT;
    std::vector<int> associatedIndexLUT;
    std::vector<uint8_t> selectedForFiltering;

    AltitudeType altitudeOf(NodeId nodeId) const {
        return tree_altitude_ops::getAltitude(altitude_, nodeId);
    }

    void computeUAO(NodeId currentNodeId, AltitudeType altitudeNodeNotInNR, bool qPropag, bool isCalculateResidue) {
        const NodeId parentNodeId = tree.getNodeParent(currentNodeId);
        const AltitudeType altitudeNodeInNR = altitudeOf(currentNodeId);
        bool flagPropag = false;
        int contrast = 0;

        if (this->isSelectedForPruning(currentNodeId)) {
            altitudeNodeNotInNR = altitudeOf(parentNodeId);
            if (this->attrs_increasing[currentNodeId] <= this->maxCriterion) {
                isCalculateResidue = hasNodeSelectedInPrimitive(currentNodeId);
            }
        }

        if (this->attrs_increasing[currentNodeId] <= this->maxCriterion) {
            if (isCalculateResidue) {
                contrast = static_cast<int>(std::abs(altitudeNodeInNR - altitudeNodeNotInNR));
            }

            if (this->maxContrastLUT[parentNodeId] >= contrast) {
                this->maxContrastLUT[currentNodeId] = this->maxContrastLUT[parentNodeId];
                this->associatedIndexLUT[currentNodeId] = this->associatedIndexLUT[parentNodeId];
            } else {
                this->maxContrastLUT[currentNodeId] = static_cast<uint8_t>(contrast);
                this->associatedIndexLUT[currentNodeId] = !qPropag
                    ? static_cast<int>(this->attrs_increasing[currentNodeId] + 1)
                    : this->associatedIndexLUT[parentNodeId];
                flagPropag = true;
            }
        }

        for (NodeId childNodeId : tree.getChildren(currentNodeId)) {
            this->computeUAO(childNodeId, altitudeNodeNotInNR, flagPropag, isCalculateResidue);
        }
    }

    void executeImpl(int maxCriterion, const std::vector<uint8_t>& selectedForFiltering) {
        this->maxCriterion = maxCriterion;
        this->selectedForFiltering = selectedForFiltering;

        for (NodeId id : tree.getAliveNodeIds()) {
            maxContrastLUT[id] = 0;
            associatedIndexLUT[id] = 0;
        }

        const NodeId rootNodeId = tree.getRoot();
        const AltitudeType level = altitudeOf(rootNodeId);
        for (NodeId childNodeId : tree.getChildren(rootNodeId)) {
            computeUAO(childNodeId, level, false, false);
        }
    }

    bool isSelectedForPruning(NodeId currentNodeId) {
        const NodeId parentNodeId = tree.getNodeParent(currentNodeId);
        if (parentNodeId == InvalidNode) {
            return false;
        }
        return this->attrs_increasing[currentNodeId] != this->attrs_increasing[parentNodeId];
    }

    bool hasNodeSelectedInPrimitive(NodeId currentNodeId) {
        std::stack<NodeId> stack;
        stack.push(currentNodeId);
        while (!stack.empty()) {
            const NodeId nodeId = stack.top();
            stack.pop();
            if (selectedForFiltering[nodeId]) {
                return true;
            }

            for (NodeId childNodeId : tree.getChildren(nodeId)) {
                if (this->attrs_increasing[childNodeId] == this->attrs_increasing[nodeId]) {
                    stack.push(childNodeId);
                }
            }
        }
        return false;
    }

public:
    UltimateAttributeOpening(MorphologicalTree& tree, const AltitudeBuffer* altitude, const float* attrs_increasing)
        : tree(tree),
          altitude_(altitude),
          maxContrastLUT(this->tree.getNumInternalNodeSlots()),
          associatedIndexLUT(this->tree.getNumInternalNodeSlots()) {
        assert(attrs_increasing != nullptr);
        this->selectedForFiltering.assign(this->tree.getNumInternalNodeSlots(), true);
        this->attrs_increasing = attrs_increasing;
    }

    UltimateAttributeOpening(MorphologicalTree& tree, const std::shared_ptr<float[]>& attrs_increasing)
        : UltimateAttributeOpening(tree, attrs_increasing.get()) {
        this->ownedAttrsIncreasing_ = attrs_increasing;
    }

    UltimateAttributeOpening(MorphologicalTree& tree, const std::vector<float>& attrs_increasing)
        : UltimateAttributeOpening(tree, attrs_increasing.data()) {}

    UltimateAttributeOpening(MorphologicalTree& tree, const float* attrs_increasing)
        : UltimateAttributeOpening(tree, nullptr, attrs_increasing) {}

    UltimateAttributeOpening(WeightedMorphologicalTree& weighted, const std::shared_ptr<float[]>& attrs_increasing)
        : UltimateAttributeOpening(weighted, attrs_increasing.get()) {
        this->ownedAttrsIncreasing_ = attrs_increasing;
    }

    UltimateAttributeOpening(WeightedMorphologicalTree& weighted, const std::vector<float>& attrs_increasing)
        : UltimateAttributeOpening(weighted, attrs_increasing.data()) {}

    UltimateAttributeOpening(WeightedMorphologicalTree& weighted, const float* attrs_increasing)
        : UltimateAttributeOpening(weighted.tree, &weighted.altitude, attrs_increasing) {}

    ~UltimateAttributeOpening() = default;

    void execute(int maxCriterion) {
        std::vector<uint8_t> tmp(this->tree.getNumInternalNodeSlots(), true);
        executeImpl(maxCriterion, tmp);
    }

    void execute(int maxCriterion, const std::vector<uint8_t>& selectedForFiltering) {
        executeImpl(maxCriterion, selectedForFiltering);
    }

    void executeWithMSER(int maxCriterion, int deltaMSER) {
        ComputerMSER mser(this->tree);
        executeImpl(maxCriterion, mser.computeMSER(deltaMSER));
    }

    ImageUInt8Ptr getMaxContrastImage() {
        const int size = this->tree.getNumColsOfImage() * this->tree.getNumRowsOfImage();
        ImageUInt8Ptr imgOut = ImageUInt8::create(this->tree.getNumColsOfImage(), this->tree.getNumRowsOfImage());
        auto out = imgOut->rawData();

        for (int pidx = 0; pidx < size; pidx++) {
            out[pidx] = this->maxContrastLUT[tree.getSmallestComponent(pidx)];
        }
        return imgOut;
    }

    ImageInt32Ptr getAssociatedImage() {
        const int size = this->tree.getNumColsOfImage() * this->tree.getNumRowsOfImage();
        ImageInt32Ptr imgOut = ImageInt32::create(this->tree.getNumColsOfImage(), this->tree.getNumRowsOfImage());
        auto out = imgOut->rawData();

        for (int pidx = 0; pidx < size; pidx++) {
            out[pidx] = this->associatedIndexLUT[tree.getSmallestComponent(pidx)];
        }
        return imgOut;
    }

    ImageUInt8Ptr getAssociatedColorImage() {
        return ImageUtils::createRandomColor(this->getAssociatedImage()->rawData(), this->tree.getNumRowsOfImage(), this->tree.getNumColsOfImage());
    }
};

} // namespace mmcfilters
