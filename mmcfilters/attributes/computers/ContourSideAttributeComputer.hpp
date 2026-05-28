#pragma once

#include "../detail/AttributeKernelSupport.hpp"
#include "detail/ContourSideAttributeMaterialization.hpp"
#include "detail/ContourSideLocalEventComputation.hpp"

#include <concepts>
#include <span>
#include <vector>

namespace mmcfilters::attributes::computers {

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
    /**
     * @brief Returns the scalar contour attributes materialized by this computer.
     */
    [[nodiscard]] std::vector<Attribute> attributes() const {
        return {
            CONTOUR_PIXELS,
            CONTOUR_PERIMETER,
            CONTOUR_SIDE_NORTH,
            CONTOUR_SIDE_WEST,
            CONTOUR_SIDE_EAST,
            CONTOUR_SIDE_SOUTH,
        };
    }

    /**
     * @brief Computes requested contour-side scalar attributes for live nodes.
     *
     * @details
     * The output buffer is indexed by dense internal node id and interpreted by
     * `context.attrNames`. The method computes the local-event contour-side
     * counts once, then projects only the requested scalar columns. Altitude and
     * dependencies are ignored because contour-side counts are topology and
     * image-domain support descriptors.
     */
    template <std::floating_point Real>
    static void compute(const AttributeComputeContext<Real>& context) {
        computeImpl(
            context.tree,
            context.buffer,
            context.attrNames,
            context.requestedAttributes);
    }

private:
    template <std::floating_point Real>
    static void computeImpl(
        const MorphologicalTree& tree,
        std::span<Real> buffer,
        const AttributeNames& attrNames,
        std::span<const Attribute> requestedAttributes) {
        const auto sideCounts = detail::ContourSideLocalEventComputation::computeContourSideCounts(tree);
        detail::ContourSideAttributeMaterialization::materializeAttributesFromContourSideCounts(
            tree,
            sideCounts,
            buffer,
            attrNames,
            requestedAttributes);
    }

public:
    /**
     * @brief Materializes contour-side attributes for one-pixel unit supports.
     *
     * A one-pixel unit support contributes one contour pixel and four exposed
     * sides, one in each cardinal direction.
     */
    template <std::floating_point Real>
    static void computeUnitRows(const UnitAttributeComputeContext<Real>& context) {
        detail::ContourSideAttributeMaterialization::materializeUnitContourSideAttributes(
            context.tree,
            context.unitProperParts,
            context.buffer,
            context.attrNames,
            context.requestedAttributes);
    }
};

} // namespace mmcfilters::attributes::computers
