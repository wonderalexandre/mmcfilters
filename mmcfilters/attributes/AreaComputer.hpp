#pragma once

#include "AttributeComputer.hpp"
#include "../trees/MorphologicalTree.hpp"

namespace mmcfilters {

/**
 * @brief Computa a área (número de pixels) de cada nó da árvore.
 */
class AreaComputer : public AttributeComputer {
public:
    std::vector<Attribute> attributes() const override { return {AREA}; }

    void compute(MorphologicalTree* tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<Attribute>&, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>&) const override{
        if (PRINT_LOG) std::cout << "\n==== AttributeComputer: Computing AREA" << std::endl;
        auto indexOf = [&](NodeId idx) { return attrNames->linearIndex(idx, AREA); };
        for (NodeId id : tree->getNodeIds()) {
            buffer[indexOf(id)] = static_cast<float>(tree->getAreaById(id));
        }
    }
};

} // namespace mmcfilters

