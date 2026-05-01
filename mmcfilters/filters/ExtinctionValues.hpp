#pragma once

#include "../trees/WeightedMorphologicalTree.hpp"
#include "../utils/AdjacencyRelation.hpp"
#include "../utils/Common.hpp"
#include "../attributes/AttributeComputedIncrementally.hpp"
#include "../contours/ContoursComputedIncrementally.hpp"

#include <cassert>
#include <limits>

namespace mmcfilters {

/**
 * @brief Record describing one regional extremum and its extinction value.
 */
struct RegionalExtremaNode {
    NodeId leaf;
    NodeId cutoffNode;
    float extinction;

    RegionalExtremaNode(NodeId leaf, NodeId cutoffNode, float extinction)
        : leaf(leaf), cutoffNode(cutoffNode), extinction(extinction) {}
};

/**
 * @brief Computes and stores extinction values for regional extrema.
 */
class ExtinctionValues {
protected:
    std::vector<RegionalExtremaNode> regionalExtremaNodes;
    const MorphologicalTree& tree;
    const AltitudeBuffer* altitude_ = nullptr;

    AltitudeType altitudeOf(NodeId nodeId) const {
        return WeightedMorphologicalTree::getAltitude(altitude_, nodeId);
    }

    void initialize(const float* attr) {
        assert(attr != nullptr);
        std::vector<NodeId> leaves = this->tree.getLeaves();
        regionalExtremaNodes.reserve(leaves.size());
        std::vector<uint8_t> visited(this->tree.getNumInternalNodeSlots(), false);
        for (NodeId leafNodeId : leaves) {
            float extinction = std::numeric_limits<float>::max();
            NodeId cutoffNodeId = leafNodeId;
            NodeId parentNodeId = this->tree.getNodeParent(cutoffNodeId);
            bool flag = true;
            while (flag && !this->tree.isRoot(cutoffNodeId)) {
                if (this->tree.getNumChildren(parentNodeId) > 1) {
                    for (NodeId sonNodeId : this->tree.getChildren(parentNodeId)) {
                        if (flag) {
                            if (visited[sonNodeId] && sonNodeId != cutoffNodeId && attr[sonNodeId] == attr[cutoffNodeId]) {
                                flag = false;
                            } else if (sonNodeId != cutoffNodeId && attr[sonNodeId] > attr[cutoffNodeId]) {
                                flag = false;
                            }
                            visited[sonNodeId] = true;
                        }
                    }
                }
                if (flag) {
                    cutoffNodeId = parentNodeId;
                    parentNodeId = this->tree.getNodeParent(cutoffNodeId);
                }
            }
            if (!this->tree.isRoot(cutoffNodeId)) {
                extinction = attr[cutoffNodeId];
            }
            regionalExtremaNodes.emplace_back(leafNodeId, cutoffNodeId, extinction);
        }

        std::sort(regionalExtremaNodes.begin(), regionalExtremaNodes.end(), [](const auto& a, const auto& b) {
            return a.extinction > b.extinction;
        });
    }

public:
    ExtinctionValues(const MorphologicalTree& tree, const std::shared_ptr<float[]>& attr)
        : ExtinctionValues(tree, attr.get()) {}

    ExtinctionValues(const MorphologicalTree& tree, const std::vector<float>& attr)
        : ExtinctionValues(tree, attr.data()) {}

    ExtinctionValues(const MorphologicalTree& tree, const float* attr)
        : tree(tree), altitude_(nullptr) {
        initialize(attr);
    }

    ExtinctionValues(const WeightedMorphologicalTree& weighted, const std::shared_ptr<float[]>& attr)
        : ExtinctionValues(weighted, attr.get()) {}

    ExtinctionValues(const WeightedMorphologicalTree& weighted, const std::vector<float>& attr)
        : ExtinctionValues(weighted, attr.data()) {}

    ExtinctionValues(const WeightedMorphologicalTree& weighted, const float* attr)
        : tree(weighted.tree_), altitude_(&weighted.altitude_) {
        initialize(attr);
    }

    ImageFloatPtr saliencyMap(int extremaToKeep, bool unweighted = true) {
        std::vector<uint8_t> keep(tree.getNumInternalNodeSlots(), false);
        std::vector<float> extinctionByNode(tree.getNumInternalNodeSlots(), 0.0f);
        std::vector<NodeId> keptNodes;
        const int leafToKeep = std::min(extremaToKeep, static_cast<int>(regionalExtremaNodes.size()));
        for (int i = 0; i < leafToKeep; ++i) {
            const NodeId cutoffNode = this->regionalExtremaNodes[i].cutoffNode;
            if (!keep[cutoffNode]) {
                keptNodes.push_back(cutoffNode);
            }
            keep[cutoffNode] = true;
            extinctionByNode[cutoffNode] = unweighted ? static_cast<float>(leafToKeep - i) : this->regionalExtremaNodes[i].extinction;
        }

        ImageFloatPtr imgOutputPtr = ImageFloat::create(tree.getNumRowsOfImage(), tree.getNumColsOfImage(), 0);
        auto saliencyOutput = imgOutputPtr->rawData();

        auto contours = ContoursComputedIncrementally::extractCompactContours(tree);
        for (NodeId node : keptNodes) {
            for (int p : contours.getContour(node)) {
                saliencyOutput[p] = extinctionByNode[node];
            }
        }

        return imgOutputPtr;
    }

    ImageUInt8Ptr filtering(int extremaToKeep) {
        std::vector<uint8_t> criterion(tree.getNumInternalNodeSlots(), false);
        const int leafToKeep = std::min(extremaToKeep, static_cast<int>(regionalExtremaNodes.size()));
        for (int i = 0; i < leafToKeep; i++) {
            criterion[regionalExtremaNodes[i].leaf] = true;
        }
        for (NodeId nodeId : tree.getPostOrderNodes()) {
            if (!tree.isRoot(nodeId) && criterion[nodeId]) {
                criterion[tree.getNodeParent(nodeId)] = true;
            }
        }

        ImageUInt8Ptr imgOutputPtr = ImageUInt8::create(tree.getNumRowsOfImage(), tree.getNumColsOfImage(), 0);
        auto imgOutput = imgOutputPtr->rawData();
        std::stack<NodeId> stack;
        stack.push(tree.getRoot());
        while (!stack.empty()) {
            const NodeId nodeId = stack.top();
            stack.pop();
            const AltitudeType level = altitudeOf(nodeId);
            for (int pixel : tree.getProperParts(nodeId)) {
                imgOutput[pixel] = static_cast<uint8_t>(level);
            }
            for (NodeId childNodeId : tree.getChildren(nodeId)) {
                if (criterion[childNodeId]) {
                    stack.push(childNodeId);
                } else {
                    for (NodeId subtreeNodeId : tree.getNodeSubtree(childNodeId)) {
                        for (int pixel : tree.getProperParts(subtreeNodeId)) {
                            imgOutput[pixel] = static_cast<uint8_t>(level);
                        }
                    }
                }
            }
        }
        return imgOutputPtr;
    }

    std::vector<RegionalExtremaNode>& getExtinctionValues() { return regionalExtremaNodes; }
};

} // namespace mmcfilters
