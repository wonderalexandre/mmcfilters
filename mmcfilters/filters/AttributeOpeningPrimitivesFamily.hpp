#pragma once

#include "../trees/MorphologicalTree.hpp"
#include "../trees/WeightedMorphologicalTree.hpp"
#include "../filters/AttributeFilters.hpp"
#include "../attributes/ComputerMSER.hpp"
#include "../attributes/AttributeComputedIncrementally.hpp"
#include "../utils/Common.hpp"

#include <cassert>

namespace mmcfilters {

/**
 * @brief Manages families of attribute-opening primitives.
 *
 * The class stores the thresholds, selected nodes, and residual image needed
 * to derive primitive families from an increasing attribute.
 */
class AttributeOpeningPrimitivesFamily {
protected:
    const float* attrs_increasing = nullptr;
    std::shared_ptr<float[]> ownedAttrsIncreasing_;
    float maxCriterion;
    std::vector<float> thresholds;
    std::vector<NodeId> nodesWithMaximumCriterium;

    MorphologicalTree& tree;
    std::unique_ptr<AttributeFilters> filters_;
    std::vector<uint8_t> selectedForFiltering;
    ImageUInt8Ptr restOfImage;
    int numPrimitives;

    void initializeRestOfImage(float threshold);
    void initializeNodesWithMaximumCriterium();

    void make_unique_vector(std::vector<float>& v) {
        std::sort(v.begin(), v.end());
        auto new_end = std::unique(v.begin(), v.end());
        v.erase(new_end, v.end());
    }

    AttributeOpeningPrimitivesFamily(
        MorphologicalTree& tree,
        std::unique_ptr<AttributeFilters> filters,
        const AltitudeBuffer* altitude,
        const float* attrs_increasing,
        float maxCriterion,
        int deltaMSER)
        : attrs_increasing{attrs_increasing},
          maxCriterion{maxCriterion},
          tree{tree},
          filters_(std::move(filters)),
          numPrimitives{0} {
        assert(this->attrs_increasing != nullptr);
        assert(this->filters_ != nullptr);
        if (deltaMSER > 0) {
            ComputerMSER mser(this->tree, altitude);
            this->selectedForFiltering = mser.computeMSER(deltaMSER);
        } else {
            this->selectedForFiltering.assign(this->tree.getNumInternalNodeSlots(), true);
        }

        float maxThreshold = 0.0f;
        for (NodeId nodeId : this->tree.getAliveNodeIds()) {
            if (this->attrs_increasing[nodeId] <= this->maxCriterion && this->isSelectedForPruning(nodeId)) {
                this->numPrimitives++;
                if (this->attrs_increasing[nodeId] > maxThreshold) {
                    maxThreshold = this->attrs_increasing[nodeId];
                }
            }
        }
        this->initializeRestOfImage(maxThreshold);
        this->initializeNodesWithMaximumCriterium();
    }

public:
    AttributeOpeningPrimitivesFamily(MorphologicalTree& tree, const std::shared_ptr<float[]>& attr, float maxCriterion)
        : AttributeOpeningPrimitivesFamily(tree, attr, maxCriterion, 0) {}

    AttributeOpeningPrimitivesFamily(MorphologicalTree& tree, const std::vector<float>& attr, float maxCriterion)
        : AttributeOpeningPrimitivesFamily(tree, attr, maxCriterion, 0) {}

    AttributeOpeningPrimitivesFamily(MorphologicalTree& tree, const float* attrs_increasing, float maxCriterion)
        : AttributeOpeningPrimitivesFamily(tree, attrs_increasing, maxCriterion, 0) {}

    AttributeOpeningPrimitivesFamily(MorphologicalTree& tree, const std::shared_ptr<float[]>& attrs_increasing, float maxCriterion, int deltaMSER)
        : AttributeOpeningPrimitivesFamily(tree, attrs_increasing.get(), maxCriterion, deltaMSER) {
        this->ownedAttrsIncreasing_ = attrs_increasing;
    }

    AttributeOpeningPrimitivesFamily(MorphologicalTree& tree, const std::vector<float>& attrs_increasing, float maxCriterion, int deltaMSER)
        : AttributeOpeningPrimitivesFamily(tree, attrs_increasing.data(), maxCriterion, deltaMSER) {}

    AttributeOpeningPrimitivesFamily(MorphologicalTree& tree, const float* attrs_increasing, float maxCriterion, int deltaMSER)
        : AttributeOpeningPrimitivesFamily(
            tree,
            std::make_unique<AttributeFilters>(tree),
            nullptr,
            attrs_increasing,
            maxCriterion,
            deltaMSER) {}

    AttributeOpeningPrimitivesFamily(WeightedMorphologicalTree& weighted, const std::shared_ptr<float[]>& attr, float maxCriterion)
        : AttributeOpeningPrimitivesFamily(weighted, attr, maxCriterion, 0) {}

    AttributeOpeningPrimitivesFamily(WeightedMorphologicalTree& weighted, const std::vector<float>& attr, float maxCriterion)
        : AttributeOpeningPrimitivesFamily(weighted, attr, maxCriterion, 0) {}

    AttributeOpeningPrimitivesFamily(WeightedMorphologicalTree& weighted, const float* attrs_increasing, float maxCriterion)
        : AttributeOpeningPrimitivesFamily(weighted, attrs_increasing, maxCriterion, 0) {}

    AttributeOpeningPrimitivesFamily(WeightedMorphologicalTree& weighted, const std::shared_ptr<float[]>& attrs_increasing, float maxCriterion, int deltaMSER)
        : AttributeOpeningPrimitivesFamily(weighted, attrs_increasing.get(), maxCriterion, deltaMSER) {
        this->ownedAttrsIncreasing_ = attrs_increasing;
    }

    AttributeOpeningPrimitivesFamily(WeightedMorphologicalTree& weighted, const std::vector<float>& attrs_increasing, float maxCriterion, int deltaMSER)
        : AttributeOpeningPrimitivesFamily(weighted, attrs_increasing.data(), maxCriterion, deltaMSER) {}

    AttributeOpeningPrimitivesFamily(WeightedMorphologicalTree& weighted, const float* attrs_increasing, float maxCriterion, int deltaMSER)
        : AttributeOpeningPrimitivesFamily(
            weighted.tree,
            std::make_unique<AttributeFilters>(weighted),
            &weighted.altitude,
            attrs_increasing,
            maxCriterion,
            deltaMSER) {}

    ~AttributeOpeningPrimitivesFamily() = default;

    std::vector<float> getThresholdsPrimitive() {
        if (this->thresholds.empty()) {
            for (NodeId nodeId : this->tree.getAliveNodeIds()) {
                if (this->attrs_increasing[nodeId] <= this->maxCriterion && this->isSelectedForPruning(nodeId)) {
                    this->thresholds.push_back(this->attrs_increasing[nodeId]);
                }
            }
            this->make_unique_vector(this->thresholds);
        }
        return thresholds;
    }

    bool isSelectedForPruning(NodeId node) {
        const NodeId parentNodeId = tree.getNodeParent(node);
        if (parentNodeId == InvalidNode) {
            return false;
        }
        return this->attrs_increasing[node] != this->attrs_increasing[parentNodeId];
    }

    bool hasNodeSelectedInPrimitive(NodeId currentNode) {
        if (!this->selectedForFiltering[currentNode]) {
            std::stack<NodeId> s;
            s.push(currentNode);
            while (!s.empty()) {
                const NodeId node = s.top();
                s.pop();
                if (selectedForFiltering[node]) {
                    return true;
                }

                for (NodeId son : tree.getChildren(node)) {
                    if (this->attrs_increasing[son] == this->attrs_increasing[node]) {
                        s.push(son);
                    }
                }
            }
            return false;
        }
        return true;
    }

    std::vector<NodeId> getNodesWithMaximumCriterium() {
        return this->nodesWithMaximumCriterium;
    }

    ImageUInt8Ptr getRestOfImage() {
        return this->restOfImage;
    }

    int getNumPrimitives() {
        return this->numPrimitives;
    }

    MorphologicalTree& getTree() {
        return this->tree;
    }
};

} // namespace mmcfilters

inline void mmcfilters::AttributeOpeningPrimitivesFamily::initializeRestOfImage(float thrRestImage) {
    this->restOfImage = this->filters_->filteringByPruningMin(this->attrs_increasing, thrRestImage);
}

inline void mmcfilters::AttributeOpeningPrimitivesFamily::initializeNodesWithMaximumCriterium() {
    std::stack<NodeId> s;
    const NodeId rootId = tree.getRoot();
    for (NodeId childNodeId : tree.getChildren(rootId)) {
        s.push(childNodeId);
    }

    while (!s.empty()) {
        const NodeId nodeId = s.top();
        s.pop();
        if (this->attrs_increasing[rootId] != this->attrs_increasing[nodeId] && this->attrs_increasing[nodeId] <= this->maxCriterion) {
            this->nodesWithMaximumCriterium.push_back(nodeId);
        } else {
            for (NodeId childNodeId : tree.getChildren(nodeId)) {
                s.push(childNodeId);
            }
        }
    }
}
