#pragma once

#include "AttributeComputerDomain.hpp"
#include "../detail/AttributeKernelSupport.hpp"
#include "../../trees/detail/TreeTraversalDetail.hpp"
#include "../../trees/MorphologicalTree.hpp"

#include <array>
#include <concepts>
#include <string_view>

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
class AreaComputer {
public:
    /// Family name used in dependency-plan diagnostics.
    static constexpr std::string_view familyName = "area";

    /// Execution domain required by the computer.
    static constexpr AttributeComputerDomain domain = AttributeComputerDomain::Topology;

    /**
     * @brief Canonical list of attributes produced by this computer.
     */
    inline static constexpr std::array<Attribute, 1> producedAttributes{AREA};

    /**
     * @brief Computes area by summing direct proper-part counts bottom-up.
     *
     * @details
     * The context buffer is written in dense internal-node order. The request
     * subset is ignored because `AreaComputer` produces a single descriptor;
     * callers are responsible for passing an `AttributeNames` layout containing
     * `AREA`.
     *
     * @param context Non-owning compute context whose layout contains `AREA`.
     */
    template <std::floating_point Real>
    static void compute(const AttributeComputeContext<Real>& context) {
        computeImpl(context.tree, context.buffer, context.attrNames);
    }

private:
    template <std::floating_point Real>
    static void computeImpl(const MorphologicalTree& tree, std::span<Real> buffer, const AttributeNames& attrNames) {
        auto indexOfArea = [&](NodeId nodeId) { return attrNames.linearIndex(nodeId, AREA); };
        ::mmcfilters::detail::traversePostOrder(
            tree,
            tree.getRoot(),
            [&](NodeId nodeId) {
                buffer[indexOfArea(nodeId)] = static_cast<Real>(tree.getNumProperParts(nodeId));
            },
            [&](NodeId parentNodeId, NodeId childNodeId) {
                buffer[indexOfArea(parentNodeId)] += buffer[indexOfArea(childNodeId)];
            },
            [](NodeId) {});
    }

public:
    /**
     * @brief Materializes `AREA` for one-pixel unit supports.
     *
     * Every exported unit proper part has area `1`.
     */
    template <std::floating_point Real>
    static void computeUnitRows(const UnitAttributeComputeContext<Real>& context) {
        requireUnitAttributeBufferShape(context.tree, context.unitProperParts, context.buffer, context.attrNames);
        if (!requestsAttribute(context.requestedAttributes, AREA)) {
            return;
        }
        for (NodeId leafIndex = 0; leafIndex < static_cast<NodeId>(context.unitProperParts.size()); ++leafIndex) {
            context.buffer[context.attrNames.linearIndex(leafIndex, AREA)] = Real{1};
        }
    }

};

} // namespace mmcfilters::attributes::computers
