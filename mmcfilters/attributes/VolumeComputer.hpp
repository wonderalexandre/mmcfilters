#pragma once

#include "AttributeComputer.hpp"
#include "AttributeComputedIncrementally.hpp"
#include "../trees/MorphologicalTree.hpp"


namespace mmcfilters {

/**
 * @brief Calcula volume e volume relativo acumulando níveis cinza sobre a árvore.
 */
class VolumeComputer : public AttributeComputer {
public:
    std::vector<Attribute> attributes() const override { return {VOLUME, RELATIVE_VOLUME}; }

    void compute(
        MorphologicalTree* tree,
        std::shared_ptr<float[]> buffer,
        std::shared_ptr<AttributeNames> attrNames,
        const std::vector<Attribute>& requestedAttributes,
        const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>&) const override
    {
        if (PRINT_LOG) std::cout << "\n==== AttributeComputer: Computing VOLUME" << std::endl;
        auto indexOfVol = [&](NodeId idx) { return attrNames->linearIndex(idx, VOLUME); };
        auto indexOfRel = [&](NodeId idx) { return attrNames->linearIndex(idx, RELATIVE_VOLUME); };

        bool computeVolume = std::find(requestedAttributes.begin(), requestedAttributes.end(), VOLUME) != requestedAttributes.end();
        bool computeRelative = std::find(requestedAttributes.begin(), requestedAttributes.end(), RELATIVE_VOLUME) != requestedAttributes.end();

        AttributeComputedIncrementally::computerAttribute(
            tree,
            tree->getRootById(),
            [&](NodeId node) {
                if (computeVolume)
                    buffer[indexOfVol(node)] = static_cast<float>(tree->getNumCNPsById(node) * tree->getLevelById(node));
                if (computeRelative)
                    buffer[indexOfRel(node)] = 0.0f;
            },
            [&](NodeId parent, NodeId child) {
                if (computeVolume)
                    buffer[indexOfVol(parent)] += buffer[indexOfVol(child)];
                if (computeRelative)
                    buffer[indexOfRel(parent)] += buffer[indexOfRel(child)] + static_cast<float>(tree->getAreaById(child) * std::abs(tree->getLevelById(child) - tree->getLevelById(parent)));
            },
            [&](NodeId node) {
                if (computeRelative)
                    buffer[indexOfRel(node)] += static_cast<float>(tree->getAreaById(node));
            });
    }
};

} // namespace mmcfilters

