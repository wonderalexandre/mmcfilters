#pragma once

#include "PendingEdgeLists.hpp"

#include <algorithm>
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
     * @brief Compacts temporary edge lists into contiguous node slices.
     * @param additions Pending additions indexed by node.
     * @param removals Pending removals indexed by node.
     * @param numPossibleEdges Size of the packed-edge domain.
     * @return Compact addition and removal slices.
     */
    [[nodiscard]] static ContourEdgeDeltaStore compact(const PendingEdgeLists& additions, const PendingEdgeLists& removals,
                                                       int numPossibleEdges) {
        const int numNodes = additions.numNodes();
        ContourEdgeDeltaStore store(numNodes);
        store.additionEdges_.reserve(additions.numEntries());
        store.removalEdges_.reserve(removals.numEntries());
        std::vector<uint16_t> edgeMarks(static_cast<std::size_t>(numPossibleEdges), 0);
        uint16_t markGeneration = 1;

        for (NodeId node = 0; node < numNodes; ++node) {
            appendNodeEdges(additions, node, store.additionEdges_, store.additionSlices_[static_cast<std::size_t>(node)], edgeMarks,
                            markGeneration);
            appendNodeEdges(removals, node, store.removalEdges_, store.removalSlices_[static_cast<std::size_t>(node)], edgeMarks,
                            markGeneration);
        }
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
     * @brief Appends one node's distinct pending edges to a compact buffer.
     * @param lists Pending edge lists.
     * @param node Source node identifier.
     * @param edges Destination packed-edge buffer.
     * @param nodeSlice Slice metadata written for `node`.
     * @param edgeMarks Deduplication generations indexed by packed edge.
     * @param markGeneration Active generation, advanced before appending.
     */
    static void appendNodeEdges(const PendingEdgeLists& lists, NodeId node, std::vector<int>& edges, Slice& nodeSlice,
                                std::vector<uint16_t>& edgeMarks, uint16_t& markGeneration) {
        advanceMarkGeneration(edgeMarks, markGeneration);
        nodeSlice.offset = checkedUint32(edges.size(), "contour edge delta offset");
        lists.appendDistinct(node, edges, edgeMarks, markGeneration);
        nodeSlice.count = checkedUint32(edges.size() - nodeSlice.offset, "contour edge delta count");
    }

    /**
     * @brief Advances the deduplication generation, clearing on wraparound.
     * @param edgeMarks Deduplication generations indexed by packed edge.
     * @param markGeneration Active generation.
     */
    static void advanceMarkGeneration(std::vector<uint16_t>& edgeMarks, uint16_t& markGeneration) {
        ++markGeneration;
        if (markGeneration == 0) {
            std::fill(edgeMarks.begin(), edgeMarks.end(), 0);
            markGeneration = 1;
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
