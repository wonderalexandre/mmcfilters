#pragma once

#include "ContourSideAttributeData.hpp"
#include "../../detail/AttributeKernelSupport.hpp"
#include "../../../trees/MorphologicalTree.hpp"
#include "../../../trees/detail/TreeTraversalDetail.hpp"

#include <concepts>
#include <cstddef>
#include <span>
#include <stdexcept>

namespace mmcfilters::attributes::computers::detail {

/** @brief Selection mask for scalar attributes derived from contour-side counts. */
struct ContourSideRequest {
    bool contourPixels = false; ///< Whether contour-pixel count is requested.
    bool perimeter = false;     ///< Whether exposed-side perimeter is requested.
    bool north = false;         ///< Whether north-side count is requested.
    bool west = false;          ///< Whether west-side count is requested.
    bool east = false;          ///< Whether east-side count is requested.
    bool south = false;         ///< Whether south-side count is requested.

    /** @brief Reports whether at least one contour scalar is requested. @return True when any request flag is set. */
    [[nodiscard]] bool any() const noexcept { return contourPixels || perimeter || north || west || east || south; }

    /**
     * @brief Builds the selection mask from requested scalar attributes.
     * @param requestedAttributes Requested scalar attributes.
     * @return Contour-side selection mask.
     */
    [[nodiscard]] static ContourSideRequest from(std::span<const Attribute> requestedAttributes) {
        return {.contourPixels = requestsAttribute(requestedAttributes, ContourPixels),
                .perimeter = requestsAttribute(requestedAttributes, ContourPerimeter),
                .north = requestsAttribute(requestedAttributes, ContourSideNorth),
                .west = requestsAttribute(requestedAttributes, ContourSideWest),
                .east = requestsAttribute(requestedAttributes, ContourSideEast),
                .south = requestsAttribute(requestedAttributes, ContourSideSouth)};
    }
};

namespace kernel {

/**
 * @brief Materializes requested contour-side scalars from established counts.
 * @param context Established tree, output layout, and output buffer.
 * @param request Scalar columns to materialize.
 * @param sideCounts Per-node contour-side counts.
 */
template <std::floating_point Real>
inline void materializeContourSideAttributes(const AttributeComputeContext<Real>& context, const ContourSideRequest& request,
                                             std::span<const ContourSideCounts> sideCounts) {
    if (!request.any()) {
        return;
    }

    const int stride = context.attrNames.NUM_ATTRIBUTES;
    const auto offsetOf = [&](Attribute attribute) { return context.attrNames.indexMap.find(attribute)->second; };
    const int contourPixelsOffset = request.contourPixels ? offsetOf(ContourPixels) : 0;
    const int perimeterOffset = request.perimeter ? offsetOf(ContourPerimeter) : 0;
    const int northOffset = request.north ? offsetOf(ContourSideNorth) : 0;
    const int westOffset = request.west ? offsetOf(ContourSideWest) : 0;
    const int eastOffset = request.east ? offsetOf(ContourSideEast) : 0;
    const int southOffset = request.south ? offsetOf(ContourSideSouth) : 0;
    const auto outputIndex = [&](NodeId node, int offset) { return static_cast<std::size_t>(node * stride + offset); };

    ::mmcfilters::detail::kernel::traversePostOrder(
        context.tree, context.tree.root(), [](NodeId) {}, [](NodeId, NodeId) {},
        [&](NodeId node) {
            const ContourSideCounts& counts = sideCounts[static_cast<std::size_t>(node)];
            if (request.contourPixels)
                context.buffer[outputIndex(node, contourPixelsOffset)] = static_cast<Real>(counts.contourPixels);
            if (request.perimeter)
                context.buffer[outputIndex(node, perimeterOffset)] = static_cast<Real>(counts.exposedSides);
            if (request.north)
                context.buffer[outputIndex(node, northOffset)] = static_cast<Real>(counts.north);
            if (request.west)
                context.buffer[outputIndex(node, westOffset)] = static_cast<Real>(counts.west);
            if (request.east)
                context.buffer[outputIndex(node, eastOffset)] = static_cast<Real>(counts.east);
            if (request.south)
                context.buffer[outputIndex(node, southOffset)] = static_cast<Real>(counts.south);
        });
}

} // namespace kernel

/**
 * @brief Internal scalar projection for precomputed contour-side buckets.
 *
 * @details
 * `ContourSideCounts` is a finite-window storage bucket. This helper is the only
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
     * @param tree Tree topology.
     * @param sideCounts Contour-side counters.
     * @param buffer Buffer read or written by the operation.
     * @param attrNames Layout that maps attributes to buffer columns.
     * @param requestedAttributes Attributes requested for materialization.
     *
     * @throws std::invalid_argument If `sideCounts` does not cover every
     * internal node slot.
     *
     */
    template <std::floating_point Real>
    static void materializeAttributesFromContourSideCounts(const MorphologicalTree& tree, std::span<const ContourSideCounts> sideCounts, std::span<Real> buffer,
                                                           const AttributeNames& attrNames, std::span<const Attribute> requestedAttributes) {
        const std::size_t numNodeSlots = static_cast<std::size_t>(tree.numInternalNodeSlots());
        if (sideCounts.size() < numNodeSlots) {
            throw std::invalid_argument("Local-event contour side counts do not cover all tree node slots.");
        }

        const AttributeComputeContext<Real> context{tree, buffer, attrNames, requestedAttributes};
        const ContourSideRequest request = ContourSideRequest::from(requestedAttributes);
        kernel::materializeContourSideAttributes(context, request, sideCounts);
    }

    /**
     * @brief Materializes contour-side attributes for one-pixel unit supports.
     *
     * A one-pixel unit support contributes one contour pixel and four exposed
     * sides, one in each cardinal direction.
     *
     * @param tree Tree topology.
     * @param unitPixels Proper-part data.
     * @param buffer Buffer read or written by the operation.
     * @param attrNames Layout that maps attributes to buffer columns.
     * @param requestedAttributes Attributes requested for materialization.
     */
    template <std::floating_point Real>
    static void materializeUnitContourSideAttributes(const MorphologicalTree& tree, std::span<const PixelId> unitPixels, std::span<Real> buffer,
                                                     const AttributeNames& attrNames, std::span<const Attribute> requestedAttributes) {
        requireUnitAttributeBufferShape(tree, unitPixels, buffer, attrNames);

        const bool computeContourPixels = requestsAttribute(requestedAttributes, ContourPixels);
        const bool computeContourPerimeter = requestsAttribute(requestedAttributes, ContourPerimeter);
        const bool computeNorth = requestsAttribute(requestedAttributes, ContourSideNorth);
        const bool computeWest = requestsAttribute(requestedAttributes, ContourSideWest);
        const bool computeEast = requestsAttribute(requestedAttributes, ContourSideEast);
        const bool computeSouth = requestsAttribute(requestedAttributes, ContourSideSouth);

        const ContourSideCounts singleton{
            .contourPixels = 1,
            .exposedSides = 4,
            .north = 1,
            .west = 1,
            .east = 1,
            .south = 1,
        };

        for (NodeId leafIndex = 0; leafIndex < static_cast<NodeId>(unitPixels.size()); ++leafIndex) {
            materializeNode(leafIndex, singleton, buffer, attrNames, computeContourPixels, computeContourPerimeter, computeNorth, computeWest, computeEast,
                            computeSouth);
        }
    }

  private:
    /**
     * @brief Writes the selected scalar projections for one output row.
     *
     * @param outputIndex Index.
     * @param counts Counters.
     * @param buffer Buffer read or written by the operation.
     * @param attrNames Layout that maps attributes to buffer columns.
     * @param computeContourPixels Flag controlling compute contour pixels.
     * @param computeContourPerimeter Flag controlling compute contour perimeter.
     * @param computeNorth Flag controlling compute north.
     * @param computeWest Flag controlling compute west.
     * @param computeEast Flag controlling compute east.
     * @param computeSouth Flag controlling compute south.
     */
    template <std::floating_point Real>
    static void materializeNode(NodeId outputIndex, const ContourSideCounts& counts, std::span<Real> buffer, const AttributeNames& attrNames,
                                bool computeContourPixels, bool computeContourPerimeter, bool computeNorth, bool computeWest, bool computeEast,
                                bool computeSouth) {
        if (computeContourPixels) {
            buffer[attrNames.linearIndex(outputIndex, ContourPixels)] = static_cast<Real>(counts.contourPixels);
        }
        if (computeContourPerimeter) {
            buffer[attrNames.linearIndex(outputIndex, ContourPerimeter)] = static_cast<Real>(counts.exposedSides);
        }
        if (computeNorth) {
            buffer[attrNames.linearIndex(outputIndex, ContourSideNorth)] = static_cast<Real>(counts.north);
        }
        if (computeWest) {
            buffer[attrNames.linearIndex(outputIndex, ContourSideWest)] = static_cast<Real>(counts.west);
        }
        if (computeEast) {
            buffer[attrNames.linearIndex(outputIndex, ContourSideEast)] = static_cast<Real>(counts.east);
        }
        if (computeSouth) {
            buffer[attrNames.linearIndex(outputIndex, ContourSideSouth)] = static_cast<Real>(counts.south);
        }
    }
};

} // namespace mmcfilters::attributes::computers::detail
