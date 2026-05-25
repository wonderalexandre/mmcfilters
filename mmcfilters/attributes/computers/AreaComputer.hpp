#pragma once

#include "../AttributeComputer.hpp"
#include "../../trees/detail/TreeTraversalDetail.hpp"
#include "../../trees/MorphologicalTree.hpp"

namespace mmcfilters::attributes::computers {

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
    using AttributeComputer::compute;
    using AttributeComputer::computeUnitAttributes;

    /**
     * @brief Returns the attribute naturally produced by this computer.
     */
    [[nodiscard]] std::vector<Attribute> attributes() const override { return {AREA}; }

    /**
     * @brief Materialises `AREA` for every live node of the tree.
     *
     * `altitude`, `requestedAttributes`, and dependencies are ignored because
     * area is topology-only and this computer owns a single scalar attribute.
     */
    void compute(const MorphologicalTree& tree, AttributeAltitudeView, std::span<float> buffer, const AttributeNames& attrNames, std::span<const Attribute>, std::span<const DependencySource>) const override{
        requireAttributeBufferShape(tree, buffer, attrNames);
        computeAreaAttribute(tree, buffer, attrNames);
    }

    /**
     * @brief Computes area by summing direct proper-part counts bottom-up.
     *
     * @param tree Tree whose dense internal node ids index `buffer`.
     * @param buffer Flat internal-node output buffer.
     * @param attrNames Layout containing `AREA`.
     */
    static void computeAreaAttribute(const MorphologicalTree& tree, std::span<float> buffer, const AttributeNames& attrNames) {
        auto indexOfArea = [&](NodeId nodeId) { return attrNames.linearIndex(nodeId, AREA); };
        ::mmcfilters::detail::traversePostOrder(
            tree,
            tree.getRoot(),
            [&](NodeId nodeId) {
                buffer[indexOfArea(nodeId)] = static_cast<float>(tree.getNumProperParts(nodeId));
            },
            [&](NodeId parentNodeId, NodeId childNodeId) {
                buffer[indexOfArea(parentNodeId)] += buffer[indexOfArea(childNodeId)];
            },
            [](NodeId) {});
    }

    /**
     * @brief Materializes `AREA` for one-pixel unit supports.
     *
     * Every exported unit proper part has area `1`.
     */
    void computeUnitAttributes(const MorphologicalTree& tree, AttributeAltitudeView, std::span<const NodeId> unitProperParts, std::span<float> buffer, const AttributeNames& attrNames, std::span<const Attribute> requestedAttributes) const override {
        requireUnitAttributeBufferShape(tree, unitProperParts, buffer, attrNames);
        if (!requestsAttribute(requestedAttributes, AREA)) {
            return;
        }
        for (NodeId leafIndex = 0; leafIndex < static_cast<NodeId>(unitProperParts.size()); ++leafIndex) {
            buffer[attrNames.linearIndex(leafIndex, AREA)] = 1.0f;
        }
    }
};

} // namespace mmcfilters::attributes::computers
