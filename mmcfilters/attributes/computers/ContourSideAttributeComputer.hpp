#pragma once

#include "AttributeComputerDomain.hpp"
#include "AttributeComputerFamily.hpp"
#include "../detail/AttributeKernelSupport.hpp"
#include "detail/ContourSideAttributeMaterialization.hpp"
#include "detail/ContourSideLocalEventComputation.hpp"
#include "../../utils/Contract.hpp"

#include <array>
#include <concepts>
#include <span>
#include <string_view>

namespace mmcfilters::attributes::computers {

namespace detail {

namespace kernel {

/**
 * @brief Computes requested contour-side attributes over an established tree.
 * @param context Established tree, output layout, and output buffer.
 * @param request Contour-side columns to materialize.
 */
template <std::floating_point Real>
inline void computeContourSideAttributes(const AttributeComputeContext<Real>& context, const ContourSideRequest& request) {
    if (!request.any()) {
        return;
    }

    const std::vector<ContourSideCounts> sideCounts = computeContourSideCounts(context.tree);
    materializeContourSideAttributes(context, request, sideCounts);
}

} // namespace kernel

template <std::floating_point Real> inline void validateContourSideContext(const AttributeComputeContext<Real>& context) {
    requireAttributeBufferShape(context.tree, context.buffer, context.attrNames);
    requireRequestedAttributeColumns(context);
    constexpr std::array<local_events::WindowOffset, 5> contourWindow{{{0, 0}, {-1, 0}, {0, -1}, {0, 1}, {1, 0}}};
    local_events::detail::validateEventEngineInput(context.tree, contourWindow);
}

} // namespace detail

/**
 * @brief Scalar computer backed by local contour-side counts.
 *
 * @details
 * This computer projects scalar attributes from its local-event computation.
 * The public attribute pipeline/topology backend invokes it for `CONTOUR_*`
 * requests.
 *
 * The local event computation is performed once and the following attributes
 * are materialized as projections of the same bucket:
 *
 * - `CONTOUR_PIXELS`: pixels with at least one exposed 4-neighbour side;
 * - `CONTOUR_PERIMETER`: total exposed 4-neighbour sides;
 * - `CONTOUR_SIDE_NORTH`, `CONTOUR_SIDE_WEST`, `CONTOUR_SIDE_EAST`,
 *   `CONTOUR_SIDE_SOUTH`: directional exposed-side counters.
 */
class ContourSideAttributeComputer {
  public:
    /// Family name used in dependency-plan diagnostics.
    static constexpr std::string_view familyName = "contour-side";

    /// Stable family id used by the scheduler.
    static constexpr AttributeComputerFamily family = AttributeComputerFamily::ContourSide;

    /// Execution domain required by the computer.
    static constexpr AttributeComputerDomain domain = AttributeComputerDomain::Topology;

    /**
     * @brief Canonical list of scalar contour attributes materialized by this computer.
     */
    inline static constexpr std::array<Attribute, 6> producedAttributes{CONTOUR_PIXELS,    CONTOUR_PERIMETER, CONTOUR_SIDE_NORTH,
                                                                        CONTOUR_SIDE_WEST, CONTOUR_SIDE_EAST, CONTOUR_SIDE_SOUTH};

    /**
     * @brief Computes requested contour-side scalar attributes for live nodes.
     *
     * @details
     * The output buffer is indexed by dense internal node id and interpreted by
     * `context.attrNames`. The method computes the local-event contour-side
     * counts once, then projects only the requested scalar columns. Altitude and
     * dependencies are ignored because contour-side counts depend only on
     * topology and regular-grid support.
     *
     * @param context Operation context or diagnostic label.
     */
    template <std::floating_point Real> static void compute(const AttributeComputeContext<Real>& context) {
        const detail::ContourSideRequest request = detail::ContourSideRequest::from(context.requestedAttributes);
        MMCFILTERS_CONTRACT_CHECKED_ONLY(detail::validateContourSideContext(context));
        detail::kernel::computeContourSideAttributes(context, request);
    }

  public:
    /**
     * @brief Materializes contour-side attributes for one-pixel unit supports.
     *
     * A one-pixel unit support contributes one contour pixel and four exposed
     * sides, one in each cardinal direction.
     *
     * @param context Operation context or diagnostic label.
     */
    template <std::floating_point Real> static void computeUnitRows(const UnitAttributeComputeContext<Real>& context) {
        detail::ContourSideAttributeMaterialization::materializeUnitContourSideAttributes(context.tree, context.unitProperParts, context.buffer,
                                                                                          context.attrNames, context.requestedAttributes);
    }
};

} // namespace mmcfilters::attributes::computers
