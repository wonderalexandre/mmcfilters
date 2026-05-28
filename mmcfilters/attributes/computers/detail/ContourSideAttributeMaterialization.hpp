#pragma once

#include "ContourSideAttributeData.hpp"
#include "../../detail/AttributeKernelSupport.hpp"
#include "../../../trees/MorphologicalTree.hpp"

#include <concepts>
#include <cstddef>
#include <span>
#include <stdexcept>

namespace mmcfilters::attributes::computers::detail {

/**
 * @brief Internal scalar projection for precomputed contour-side buckets.
 *
 * @details
 * `ContourSideCounts` is a local-event storage bucket. This helper is the only
 * layer that knows how to translate that bucket into public `CONTOUR_*`
 * scalar attributes, keeping the public computer API expressed in terms of
 * requested attributes and output buffers.
 */
class ContourSideAttributeMaterialization {
public:
    /**
     * @brief Projects precomputed side counters into scalar contour attributes.
     *
     * @details
     * This helper makes the no-recomputation path explicit for detail code and
     * tests: callers can compute `ContourSideCounts` once, then write any
     * subset of scalar projections from the same dense node-slot buffer.
     *
     * @throws std::invalid_argument If `sideCounts` does not cover every
     * internal node slot.
     */
    template <std::floating_point Real>
    static void materializeAttributesFromContourSideCounts(
        const MorphologicalTree& tree,
        std::span<const ContourSideCounts> sideCounts,
        std::span<Real> buffer,
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
    template <std::floating_point Real>
    static void materializeUnitContourSideAttributes(
        const MorphologicalTree& tree,
        std::span<const NodeId> unitProperParts,
        std::span<Real> buffer,
        const AttributeNames& attrNames,
        std::span<const Attribute> requestedAttributes) {
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
    template <std::floating_point Real>
    static void materializeNode(
        NodeId outputIndex,
        const ContourSideCounts& counts,
        std::span<Real> buffer,
        const AttributeNames& attrNames,
        bool computeContourPixels,
        bool computeContourPerimeter,
        bool computeNorth,
        bool computeWest,
        bool computeEast,
        bool computeSouth) {
        if (computeContourPixels) {
            buffer[attrNames.linearIndex(outputIndex, CONTOUR_PIXELS)] =
                static_cast<Real>(counts.contourPixels);
        }
        if (computeContourPerimeter) {
            buffer[attrNames.linearIndex(outputIndex, CONTOUR_PERIMETER)] =
                static_cast<Real>(counts.exposedSides);
        }
        if (computeNorth) {
            buffer[attrNames.linearIndex(outputIndex, CONTOUR_SIDE_NORTH)] =
                static_cast<Real>(counts.north);
        }
        if (computeWest) {
            buffer[attrNames.linearIndex(outputIndex, CONTOUR_SIDE_WEST)] =
                static_cast<Real>(counts.west);
        }
        if (computeEast) {
            buffer[attrNames.linearIndex(outputIndex, CONTOUR_SIDE_EAST)] =
                static_cast<Real>(counts.east);
        }
        if (computeSouth) {
            buffer[attrNames.linearIndex(outputIndex, CONTOUR_SIDE_SOUTH)] =
                static_cast<Real>(counts.south);
        }
    }
};

} // namespace mmcfilters::attributes::computers::detail
