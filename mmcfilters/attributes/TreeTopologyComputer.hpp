#pragma once

#include "AttributeComputer.hpp"
#include "AttributeComputedIncrementally.hpp"
#include "../trees/MorphologicalTree.hpp"


namespace mmcfilters {

/**
 * @brief Calcula atributos estruturais relacionados à topologia da árvore.
 */
class TreeTopologyComputer : public AttributeComputer {
public:
    std::vector<Attribute> attributes() const override {
        return {HEIGHT_NODE, DEPTH_NODE, IS_LEAF_NODE, IS_ROOT_NODE, NUM_CHILDREN_NODE, NUM_SIBLINGS_NODE, NUM_DESCENDANTS_NODE, NUM_LEAF_DESCENDANTS_NODE, LEAF_RATIO_NODE, BALANCE_NODE, AVG_CHILD_HEIGHT_NODE};
    }

    void compute(MorphologicalTree* tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<Attribute>& requestedAttributes, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>&) const override {
        if (PRINT_LOG) std::cout << "\n==== AttributeComputer: Computing STRUCTURE_TREE group" << std::endl;

        bool computeHeight = std::find(requestedAttributes.begin(), requestedAttributes.end(), HEIGHT_NODE) != requestedAttributes.end();
        bool computeDepth = std::find(requestedAttributes.begin(), requestedAttributes.end(), DEPTH_NODE) != requestedAttributes.end();
        bool computeIsLeaf = std::find(requestedAttributes.begin(), requestedAttributes.end(), IS_LEAF_NODE) != requestedAttributes.end();
        bool computeIsRoot = std::find(requestedAttributes.begin(), requestedAttributes.end(), IS_ROOT_NODE) != requestedAttributes.end();
        bool computeNumChildren = std::find(requestedAttributes.begin(), requestedAttributes.end(), NUM_CHILDREN_NODE) != requestedAttributes.end();
        bool computeNumSiblings = std::find(requestedAttributes.begin(), requestedAttributes.end(), NUM_SIBLINGS_NODE) != requestedAttributes.end();
        bool computeNumDescendants = std::find(requestedAttributes.begin(), requestedAttributes.end(), NUM_DESCENDANTS_NODE) != requestedAttributes.end();
        bool computeNumLeafDescendants = std::find(requestedAttributes.begin(), requestedAttributes.end(), NUM_LEAF_DESCENDANTS_NODE) != requestedAttributes.end();
        bool computeLeafRatio = std::find(requestedAttributes.begin(), requestedAttributes.end(), LEAF_RATIO_NODE) != requestedAttributes.end();
        bool computeBalance = std::find(requestedAttributes.begin(), requestedAttributes.end(), BALANCE_NODE) != requestedAttributes.end();
        bool computeAvgChildHeight = std::find(requestedAttributes.begin(), requestedAttributes.end(), AVG_CHILD_HEIGHT_NODE) != requestedAttributes.end();

        std::shared_ptr<float[]> bufferHeight = computeHeight ? buffer : std::shared_ptr<float[]>(new float[tree->getNumNodes()]);
        auto indexOfHeight = [&](NodeId idx) {
            return computeHeight ? attrNames->linearIndex(idx, HEIGHT_NODE) : idx;
        };

        std::shared_ptr<float[]> bufferNumDesc = computeNumDescendants ? buffer : std::shared_ptr<float[]>(new float[tree->getNumNodes()]);
        auto indexOfNumDescendants = [&](NodeId idx) {
            return computeNumDescendants ? attrNames->linearIndex(idx, NUM_DESCENDANTS_NODE) : idx;
        };

        std::shared_ptr<float[]> bufferNumLeafDesc = computeNumLeafDescendants ? buffer : std::shared_ptr<float[]>(new float[tree->getNumNodes()]);
        auto indexOfNumLeafDescendants = [&](NodeId idx) {
            return computeNumLeafDescendants ? attrNames->linearIndex(idx, NUM_LEAF_DESCENDANTS_NODE) : idx;
        };

        std::shared_ptr<float[]> bufferDepth = computeDepth ? buffer : std::shared_ptr<float[]>(new float[tree->getNumNodes()]);
        auto indexOfDepth = [&](NodeId idx) {
            return computeDepth ? attrNames->linearIndex(idx, DEPTH_NODE) : idx;
        };

        AttributeComputedIncrementally::computerAttribute(tree,
            tree->getRootById(),
            [&](NodeId node) {
                NodeId parent = tree->getParentById(node);

                float parentDepth = (parent != InvalidNode) ? bufferDepth[indexOfDepth(parent)] : InvalidNode;
                bufferDepth[indexOfDepth(node)] = parent != InvalidNode ? parentDepth + 1.0f : 0.0f;
                
                bufferHeight[indexOfHeight(node)] = 0.0f;
                bufferNumDesc[indexOfNumDescendants(node)] = 0.0f;
                bufferNumLeafDesc[indexOfNumLeafDescendants(node)] = tree->getNumChildrenById(node) == 0 ? 1.0f : 0.0f;

                if (computeHeight)
                    buffer[attrNames->linearIndex(node, HEIGHT_NODE)] = 0.0f;
                if (computeIsLeaf)
                    buffer[attrNames->linearIndex(node, IS_LEAF_NODE)] = tree->getNumChildrenById(node) == 0 ? 1.0f : 0.0f;
                if (computeIsRoot)
                    buffer[attrNames->linearIndex(node, IS_ROOT_NODE)] = parent!=InvalidNode ? 0.0f : 1.0f;
                if (computeNumChildren)
                    buffer[attrNames->linearIndex(node, NUM_CHILDREN_NODE)] = static_cast<float>(tree->getNumChildrenById(node));
                if (computeNumSiblings)
                    buffer[attrNames->linearIndex(node, NUM_SIBLINGS_NODE)] = parent!=InvalidNode ? static_cast<float>( tree->getNumChildrenById(parent) - 1) : 0.0f;
                if (computeLeafRatio)
                    buffer[attrNames->linearIndex(node, LEAF_RATIO_NODE)] = 0.0f;
                if (computeBalance)
                    buffer[attrNames->linearIndex(node, BALANCE_NODE)] = 0.0f;
                if (computeAvgChildHeight)
                    buffer[attrNames->linearIndex(node, AVG_CHILD_HEIGHT_NODE)] = 0.0f;
            },
            [&](NodeId parent, NodeId child) {
                
                bufferNumDesc[indexOfNumDescendants(parent)] += bufferNumDesc[indexOfNumDescendants(child)] + 1.0f;
                bufferNumLeafDesc[indexOfNumLeafDescendants(parent)] += bufferNumLeafDesc[indexOfNumLeafDescendants(child)];

                float childHeight = bufferHeight[indexOfHeight(child)];
                float& parentHeight = bufferHeight[indexOfHeight(parent)];
                parentHeight = std::max(parentHeight, childHeight + 1.0f);
                int numChildren = tree->getNumChildrenById(parent);

                if (computeBalance) {
                    float& minH = buffer[attrNames->linearIndex(parent, BALANCE_NODE)];
                    if (numChildren == 1)
                        minH = childHeight;
                    else
                        minH = std::min(minH, childHeight);
                }

                if (computeAvgChildHeight) {
                    float& sumH = buffer[attrNames->linearIndex(parent, AVG_CHILD_HEIGHT_NODE)];
                    if (numChildren == 1)
                        sumH = childHeight;
                    else
                        sumH += childHeight;
                }
            },
            [&](NodeId idx) {
                
                if (computeLeafRatio) {
                    float desc = bufferNumDesc[indexOfNumDescendants(idx)];
                    float folhas = bufferNumLeafDesc[indexOfNumLeafDescendants(idx)];
                    buffer[attrNames->linearIndex(idx, LEAF_RATIO_NODE)] = desc > 0.0f ? folhas / (desc + 1.0f) : 1.0f;
                }

                if (tree->getNumChildrenById(idx) > 0) {
                    if (computeBalance) {
                        float alturaMax = bufferHeight[indexOfHeight(idx)];
                        float alturaMin = buffer[attrNames->linearIndex(idx, BALANCE_NODE)];
                        buffer[attrNames->linearIndex(idx, BALANCE_NODE)] = alturaMax - alturaMin;
                    }

                    if (computeAvgChildHeight) {
                        buffer[attrNames->linearIndex(idx, AVG_CHILD_HEIGHT_NODE)] = buffer[attrNames->linearIndex(idx, AVG_CHILD_HEIGHT_NODE)] /
                                                                                     static_cast<float>(tree->getNumChildrenById(idx));
                    }
                }
            }
        );
    }
};

} // namespace mmcfilters

