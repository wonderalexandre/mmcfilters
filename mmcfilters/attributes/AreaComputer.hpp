#pragma once

#include "AttributeComputer.hpp"
#include "AttributeComputedIncrementally.hpp"
#include "../trees/MorphologicalTree.hpp"

namespace mmcfilters {

/**
 * @brief Computes the canonical subtree area attribute.
 *
 * @details
 * `AREA` is the most basic size descriptor in the hierarchy and is defined as
 * the number of proper parts covered by the full support of a node, including
 * the supports contributed by all of its descendants. In image-domain trees,
 * this corresponds to the number of pixels represented by the connected
 * component associated with the node.
 *
 * The implementation follows a post-order accumulation:
 * - the local contribution of a node is the number of proper parts directly
 *   owned by that node;
 * - each child then contributes its already accumulated subtree area to the
 *   parent.
 *
 * Several other attribute computers depend on `AREA`, so this computer serves
 * as a foundational building block for the incremental attribute pipeline.
 */
class AreaComputer : public AttributeComputer {
public:
    /**
     * @brief Returns the attribute naturally produced by this computer.
     */
    std::vector<Attribute> attributes() const override { return {AREA}; }

    /**
     * @brief Materialises `AREA` for every live node of the tree.
     */
    void compute(const MorphologicalTree& tree, const AltitudeBuffer*, std::span<float> buffer, const AttributeNames& attrNames, std::span<const Attribute>, std::span<const DependencySource>) const override{
        if (PRINT_LOG) std::cout << "\n==== AttributeComputer: Computing AREA" << std::endl;
        auto indexOf = [&](NodeId idx) { return attrNames.linearIndex(idx, AREA); };
        AttributeComputedIncrementally::traversePostOrder(
            tree,
            tree.getRoot(),
            [&](NodeId nodeId) {
                buffer[indexOf(nodeId)] = static_cast<float>(tree.getNumProperParts(nodeId));
            },
            [&](NodeId parentNodeId, NodeId childNodeId) {
                buffer[indexOf(parentNodeId)] += buffer[indexOf(childNodeId)];
            },
            [](NodeId) {}
        );
    }
};

} // namespace mmcfilters
