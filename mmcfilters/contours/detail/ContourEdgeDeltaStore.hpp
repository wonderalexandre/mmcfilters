#pragma once

#include "../../utils/Common.hpp"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace mmcfilters::contours::detail {

/**
 * @brief Compact contour edge changes indexed by tree node.
 *
 * Each edge is stored as `4 * pixel + side`. Addition and removal buffers use
 * one contiguous slice per internal node.
 */
class ContourEdgeDeltaStore {
  public:
    /** @brief One distinct packed-edge change assigned to a tree node. */
    struct Event {
        NodeId node = InvalidNode; ///< Node receiving the edge change.
        int packedEdge = -1;       ///< Packed support-pixel edge.
    };

    /**
     * @brief Creates empty edge slices for every internal node.
     * @param numNodes Number of internal node slots.
     */
    explicit ContourEdgeDeltaStore(int numNodes = 0)
        : additionSlices_(static_cast<std::size_t>(numNodes)), removalSlices_(static_cast<std::size_t>(numNodes)) {}

    /**
     * @brief Returns the contour edges added at one node.
     * @param node Node identifier.
     * @return Borrowed span of packed edges added at `node`.
     */
    [[nodiscard]] std::span<const int> additions(NodeId node) const noexcept { return slice(additionEdges_, additionSlices_, node); }

    /**
     * @brief Returns the contour edges removed at one node.
     * @param node Node identifier.
     * @return Borrowed span of packed edges removed at `node`.
     */
    [[nodiscard]] std::span<const int> removals(NodeId node) const noexcept { return slice(removalEdges_, removalSlices_, node); }

    /**
     * @brief Groups distinct edge-change events into contiguous node slices.
     * @param numNodes Number of internal node slots.
     * @param additions Distinct additions in generation order.
     * @param removals Distinct removals in generation order.
     * @return Compact addition and removal slices.
     */
    [[nodiscard]] static ContourEdgeDeltaStore groupDistinct(int numNodes, std::span<const Event> additions,
                                                             std::span<const Event> removals) {
        ContourEdgeDeltaStore store(numNodes);
        std::vector<uint32_t> writePositions(static_cast<std::size_t>(numNodes));
        groupEvents(additions, store.additionEdges_, store.additionSlices_, writePositions);
        groupEvents(removals, store.removalEdges_, store.removalSlices_, writePositions);
        return store;
    }

  private:
    /** @brief Half-open location represented by an offset and count. */
    struct Slice {
        /// Offset of the first edge in the shared buffer.
        uint32_t offset = 0;
        /// Number of edges in the slice.
        uint32_t count = 0;
    };

    /**
     * @brief Returns one node slice from a shared edge buffer.
     * @param edges Shared packed-edge buffer.
     * @param slices Slice metadata indexed by node.
     * @param node Node identifier.
     * @return Borrowed packed-edge span for `node`.
     */
    [[nodiscard]] static std::span<const int> slice(const std::vector<int>& edges, const std::vector<Slice>& slices,
                                                    NodeId node) noexcept {
        const Slice& nodeSlice = slices[static_cast<std::size_t>(node)];
        if (nodeSlice.count == 0) {
            return {};
        }
        return std::span<const int>(edges.data() + nodeSlice.offset, static_cast<std::size_t>(nodeSlice.count));
    }

    /**
     * @brief Groups one event stream with a counting pass followed by a direct fill.
     * @param events Distinct changes in generation order.
     * @param edges Destination buffer grouped by node.
     * @param slices Destination slices indexed by node.
     * @param writePositions Reusable end positions for direct placement.
     */
    static void groupEvents(std::span<const Event> events, std::vector<int>& edges, std::vector<Slice>& slices,
                            std::vector<uint32_t>& writePositions) {
        for (Slice& nodeSlice : slices) {
            nodeSlice = {};
        }
        for (const Event& event : events) {
            assert(event.node >= 0 && event.node < static_cast<NodeId>(slices.size()));
            ++slices[static_cast<std::size_t>(event.node)].count;
        }

        std::size_t nextOffset = 0;
        for (std::size_t node = 0; node < slices.size(); ++node) {
            Slice& nodeSlice = slices[node];
            nodeSlice.offset = checkedUint32(nextOffset, "contour edge delta offset");
            nextOffset += nodeSlice.count;
            writePositions[node] = checkedUint32(nextOffset, "contour edge delta end");
        }

        edges.resize(events.size());
        for (const Event& event : events) {
            uint32_t& writePosition = writePositions[static_cast<std::size_t>(event.node)];
            edges[static_cast<std::size_t>(--writePosition)] = event.packedEdge;
        }
    }

    /**
     * @brief Converts a checked buffer size to `uint32_t`.
     * @param value Buffer offset or count.
     * @param context Operation name used in diagnostics.
     * @return `value` represented as `uint32_t`.
     */
    [[nodiscard]] static uint32_t checkedUint32(std::size_t value, const char* context) {
        if (value > static_cast<std::size_t>(std::numeric_limits<uint32_t>::max())) {
            throw std::overflow_error(std::string(context) + " exceeds uint32_t.");
        }
        return static_cast<uint32_t>(value);
    }

    std::vector<int> additionEdges_;       ///< Packed additions shared by all nodes.
    std::vector<int> removalEdges_;        ///< Packed removals shared by all nodes.
    std::vector<Slice> additionSlices_;    ///< Addition slice for each node slot.
    std::vector<Slice> removalSlices_;     ///< Removal slice for each node slot.
};

} // namespace mmcfilters::contours::detail
