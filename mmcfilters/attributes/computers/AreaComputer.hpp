#pragma once

#include "AttributeComputerDomain.hpp"
#include "AttributeComputerFamily.hpp"
#include "../detail/AttributeKernelSupport.hpp"
#include "../../trees/detail/TreeTraversalDetail.hpp"
#include "../../trees/detail/CommittedTreeAccess.hpp"
#include "../../utils/Contract.hpp"
#include "../../trees/MorphologicalTree.hpp"

#include <array>
#include <concepts>
#include <stdexcept>
#include <string_view>

namespace mmcfilters::attributes::computers {

namespace detail {

template <std::floating_point Real> inline void validateAreaContext(const AttributeComputeContext<Real>& context) {
    requireAttributeBufferShape(context.tree, context.buffer, context.attrNames);
    if (!context.attrNames.contains(AREA)) {
        throw std::invalid_argument("AREA computation requires an AREA column in the output layout.");
    }
}

namespace kernel {

/** @brief Computes subtree area over an established tree. @param context Established tree, AREA column, and output buffer. */
template <std::floating_point Real> inline void computeArea(const AttributeComputeContext<Real>& context) {
    const int stride = context.attrNames.NUM_ATTRIBUTES;
    const int offset = context.attrNames.indexMap.find(AREA)->second;
    auto indexOfArea = [&](NodeId node) { return static_cast<std::size_t>(node * stride + offset); };
    ::mmcfilters::detail::kernel::traversePostOrder(
        context.tree, context.tree.getRoot(),
        [&](NodeId node) {
            context.buffer[indexOfArea(node)] = static_cast<Real>(::mmcfilters::detail::CommittedTreeAccess::numProperParts(context.tree, node));
        },
        [&](NodeId parent, NodeId child) { context.buffer[indexOfArea(parent)] += context.buffer[indexOfArea(child)]; }, [](NodeId) {});
}

} // namespace kernel
} // namespace detail

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

    /// Stable family id used by the scheduler.
    static constexpr AttributeComputerFamily family = AttributeComputerFamily::Area;

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
    template <std::floating_point Real> static void compute(const AttributeComputeContext<Real>& context) {
        MMCFILTERS_CONTRACT_CHECKED_ONLY(detail::validateAreaContext(context));
        detail::kernel::computeArea(context);
    }

  public:
    /**
     * @brief Materializes `AREA` for one-pixel unit supports.
     *
     * Every exported unit proper part has area `1`.
     *
     * @param context Operation context or diagnostic label.
     */
    template <std::floating_point Real> static void computeUnitRows(const UnitAttributeComputeContext<Real>& context) {
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
