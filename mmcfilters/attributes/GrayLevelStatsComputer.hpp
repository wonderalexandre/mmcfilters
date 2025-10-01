#pragma once

#include "AttributeComputer.hpp"
#include "AttributeComputedIncrementally.hpp"
#include "../trees/MorphologicalTree.hpp"


namespace mmcfilters {

/**
 * @brief Calcula estatísticas básicas de níveis de cinza (média, variância, altura).
 */
class GrayLevelStatsComputer : public AttributeComputer {
public:
    std::vector<Attribute> attributes() const override {
        return {LEVEL, MEAN_LEVEL, VARIANCE_LEVEL, GRAY_HEIGHT};
    }

    std::vector<AttributeOrGroup> requiredAttributes() const override {
        return {VOLUME};
    }

    void compute(
        MorphologicalTree* tree,
        std::shared_ptr<float[]> buffer,
        std::shared_ptr<AttributeNames> attrNames,
        const std::vector<Attribute>& requestedAttributes,
        const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources) const override
    {
        if (PRINT_LOG) std::cout << "\n==== AttributeComputer: Computing GrayLevelStatsComputer " << std::endl;

        auto indexOfMean = [&](NodeId idx) { return attrNames->linearIndex(idx, MEAN_LEVEL); };
        auto indexOfLevel = [&](NodeId idx) { return attrNames->linearIndex(idx, LEVEL); };
        auto indexOfVariance = [&](NodeId idx) { return attrNames->linearIndex(idx, VARIANCE_LEVEL); };
        auto indexOfGrayHeight = [&](NodeId idx) { return attrNames->linearIndex(idx, GRAY_HEIGHT); };

        bool computeMeanLevel = std::find(requestedAttributes.begin(), requestedAttributes.end(), MEAN_LEVEL) != requestedAttributes.end();
        bool computeVarianceLevel = std::find(requestedAttributes.begin(), requestedAttributes.end(), VARIANCE_LEVEL) != requestedAttributes.end();
        bool computeLevel = std::find(requestedAttributes.begin(), requestedAttributes.end(), LEVEL) != requestedAttributes.end();
        bool computeGrayHeight = std::find(requestedAttributes.begin(), requestedAttributes.end(), GRAY_HEIGHT) != requestedAttributes.end();

        auto [dependencyAttrNamesVol, bufferVol] = dependencySources[0];
        auto indexOfVol = [&](NodeId idx) { return dependencyAttrNamesVol->linearIndex(idx, VOLUME); };

        std::shared_ptr<long[]> sumGrayLevelSquare = nullptr;
        if (computeVarianceLevel) {
            sumGrayLevelSquare = std::shared_ptr<long[]>(new long[tree->getNumNodes()]);
        }

        AttributeComputedIncrementally::computerAttribute(
            tree,
            tree->getRootById(),
            [&](NodeId node) {
                if (computeVarianceLevel)
                    sumGrayLevelSquare[node] = static_cast<long>(tree->getNumCNPsById(node) * std::pow(tree->getLevelById(node), 2));
                if (computeLevel)
                    buffer[indexOfLevel(node)] = static_cast<float>(tree->getLevelById(node));
                if (computeGrayHeight)
                    buffer[indexOfGrayHeight(node)] = static_cast<float>(tree->getLevelById(node));
            },
            [&](NodeId parent, NodeId child) {
                if (computeVarianceLevel)
                    sumGrayLevelSquare[parent] += sumGrayLevelSquare[child];
                if (computeGrayHeight) {
                    float childValue = buffer[indexOfGrayHeight(child)];
                    float& parentValue = buffer[indexOfGrayHeight(parent)];
                    if (tree->isMaxtree() || tree->isMaxtreeNodeById(parent))
                        parentValue = std::max(parentValue, childValue);
                    else
                        parentValue = std::min(parentValue, childValue);
                }
            },
            [&](NodeId node) {
                float area = static_cast<float>(tree->getAreaById(node));
                if (computeMeanLevel)
                    buffer[indexOfMean(node)] = bufferVol[indexOfVol(node)] / area;
                if (computeVarianceLevel) {
                    float meanGrayLevel = bufferVol[indexOfVol(node)] / area;
                    double meanGrayLevelSquare = sumGrayLevelSquare[node] / area;
                    float var = static_cast<float>(meanGrayLevelSquare - (meanGrayLevel * meanGrayLevel));
                    buffer[indexOfVariance(node)] = var > 0.0f ? var : 0.0f;
                }
            });

        if (computeGrayHeight) {
            for (NodeId node : tree->getIteratorPostOrderTraversalById()) {
                if (tree->isLeafById(node))
                    buffer[indexOfGrayHeight(node)] = 0.0f;
                else
                    buffer[indexOfGrayHeight(node)] = std::abs(tree->getLevelById(node) - buffer[indexOfGrayHeight(node)]) + 1.0f;
            }
        }
    }
};

} // namespace mmcfilters

