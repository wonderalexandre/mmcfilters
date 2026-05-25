#pragma once

#include "../AttributeComputer.hpp"
#include "ContourSideAttributeData.hpp"
#include "detail/ContourSideLocalEventComputation.hpp"

#include <span>
#include <stdexcept>
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
class ContourSideAttributeComputer : public AttributeComputer {
public:
    using AttributeComputer::compute;
    using AttributeComputer::computeUnitAttributes;

    /**
     * @brief Projected local contour-side bucket consumed by this computer.
     */
    using ContourSideCounts = ::mmcfilters::attributes::computers::ContourSideCounts;

    /**
     * @brief Returns the scalar contour attributes materialized by this computer.
     */
    [[nodiscard]] std::vector<Attribute> attributes() const override {
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
     * The output buffer is indexed by dense internal node id. Altitude and
     * dependencies are ignored because contour-side counts are topology and
     * image-domain support descriptors.
     */
    void compute(
        const MorphologicalTree& tree,
        AttributeAltitudeView,
        std::span<float> buffer,
        const AttributeNames& attrNames,
        std::span<const Attribute> requestedAttributes,
        std::span<const DependencySource>) const override {
        const auto sideCounts = detail::ContourSideLocalEventComputation::computeContourSideCounts(tree);
        materializeAttributesFromContourSideCounts(tree, sideCounts, buffer, attrNames, requestedAttributes);
    }

    /**
     * @brief Projects precomputed side counters into scalar contour attributes.
     *
     * @details
     * This helper makes the no-recomputation path explicit: callers can compute
     * `ContourSideCounts` once, then write any subset of scalar projections from
     * the same dense node-slot buffer.
     *
     * @throws std::invalid_argument If `sideCounts` does not cover every
     * internal node slot.
     */
    static void materializeAttributesFromContourSideCounts(
        const MorphologicalTree& tree,
        std::span<const ContourSideCounts> sideCounts,
        std::span<float> buffer,
        const AttributeNames& attrNames,
        std::span<const Attribute> requestedAttributes) {
        const std::size_t numNodeSlots = static_cast<std::size_t>(tree.getNumInternalNodeSlots());
        if (sideCounts.size() < numNodeSlots) {
            throw std::invalid_argument("Local-event contour side counts do not cover all tree node slots.");
        }

        const bool computeContourPixels = requestsAttribute(requestedAttributes, CONTOUR_PIXELS);
        const bool computeContourPerimeter = requestsAttribute(requestedAttributes, CONTOUR_PERIMETER);
        const bool computeNorth = requestsAttribute(requestedAttributes, CONTOUR_SIDE_NORTH);
        const bool computeWest = requestsAttribute(requestedAttributes, CONTOUR_SIDE_WEST);
        const bool computeEast = requestsAttribute(requestedAttributes, CONTOUR_SIDE_EAST);
        const bool computeSouth = requestsAttribute(requestedAttributes, CONTOUR_SIDE_SOUTH);

        for (NodeId nodeId : tree.getAliveNodeIds()) {
            materializeNode(
                nodeId,
                sideCounts[static_cast<std::size_t>(nodeId)],
                buffer,
                attrNames,
                computeContourPixels,
                computeContourPerimeter,
                computeNorth,
                computeWest,
                computeEast,
                computeSouth);
        }
    }

    /**
     * @brief Materializes contour-side attributes for one-pixel unit supports.
     *
     * A one-pixel unit support contributes one contour pixel and four exposed
     * sides, one in each cardinal direction.
     */
    void computeUnitAttributes(
        const MorphologicalTree& tree,
        AttributeAltitudeView,
        std::span<const NodeId> unitProperParts,
        std::span<float> buffer,
        const AttributeNames& attrNames,
        std::span<const Attribute> requestedAttributes) const override {
        requireUnitAttributeBufferShape(tree, unitProperParts, buffer, attrNames);

        const bool computeContourPixels = requestsAttribute(requestedAttributes, CONTOUR_PIXELS);
        const bool computeContourPerimeter = requestsAttribute(requestedAttributes, CONTOUR_PERIMETER);
        const bool computeNorth = requestsAttribute(requestedAttributes, CONTOUR_SIDE_NORTH);
        const bool computeWest = requestsAttribute(requestedAttributes, CONTOUR_SIDE_WEST);
        const bool computeEast = requestsAttribute(requestedAttributes, CONTOUR_SIDE_EAST);
        const bool computeSouth = requestsAttribute(requestedAttributes, CONTOUR_SIDE_SOUTH);

        const ContourSideCounts singleton{
            .contourPixels = 1,
            .exposedSides = 4,
            .north = 1,
            .west = 1,
            .east = 1,
            .south = 1,
        };

        for (NodeId leafIndex = 0; leafIndex < static_cast<NodeId>(unitProperParts.size()); ++leafIndex) {
            materializeNode(
                leafIndex,
                singleton,
                buffer,
                attrNames,
                computeContourPixels,
                computeContourPerimeter,
                computeNorth,
                computeWest,
                computeEast,
                computeSouth);
        }
    }

private:
    /**
     * @brief Writes the selected scalar projections for one output row.
     */
    static void materializeNode(
        NodeId outputIndex,
        const ContourSideCounts& counts,
        std::span<float> buffer,
        const AttributeNames& attrNames,
        bool computeContourPixels,
        bool computeContourPerimeter,
        bool computeNorth,
        bool computeWest,
        bool computeEast,
        bool computeSouth) {
        if (computeContourPixels) {
            buffer[attrNames.linearIndex(outputIndex, CONTOUR_PIXELS)] =
                static_cast<float>(counts.contourPixels);
        }
        if (computeContourPerimeter) {
            buffer[attrNames.linearIndex(outputIndex, CONTOUR_PERIMETER)] =
                static_cast<float>(counts.exposedSides);
        }
        if (computeNorth) {
            buffer[attrNames.linearIndex(outputIndex, CONTOUR_SIDE_NORTH)] =
                static_cast<float>(counts.north);
        }
        if (computeWest) {
            buffer[attrNames.linearIndex(outputIndex, CONTOUR_SIDE_WEST)] =
                static_cast<float>(counts.west);
        }
        if (computeEast) {
            buffer[attrNames.linearIndex(outputIndex, CONTOUR_SIDE_EAST)] =
                static_cast<float>(counts.east);
        }
        if (computeSouth) {
            buffer[attrNames.linearIndex(outputIndex, CONTOUR_SIDE_SOUTH)] =
                static_cast<float>(counts.south);
        }
    }
};

} // namespace mmcfilters::attributes::computers
