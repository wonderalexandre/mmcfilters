#pragma once

#include "PendingPixelLists.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace mmcfilters::detail {

/**
 * @brief Compact per-node delta representation for boundary trace edges.
 *
 * Values are packed side ids, using `4 * pixel + side_index`. The structure is
 * intentionally parallel to `ContourDeltaStore`, but it stores side-level
 * boundary events instead of contour pixels.
 */
struct ContourTraceDeltaStore {
    /**
     * @brief Offset and length of one node slice in a flat value vector.
     */
    struct Span {
        /** @brief Stores the offset. */
        uint32_t offset = 0;
        /** @brief Stores the size. */
        uint32_t size = 0;
    };

    /// Concatenated compacted local boundary-edge additions.
    std::vector<int> addValues;
    /// Concatenated compacted local boundary-edge removals.
    std::vector<int> removeValues;
    /// Per-node spans into `addValues`.
    std::vector<Span> addSpans;
    /// Per-node spans into `removeValues`.
    std::vector<Span> removeSpans;

    /**
     * @brief Constructs `ContourTraceDeltaStore` from the supplied inputs.
     *
     * @param numNodes Number of internal nodes.
     */
    explicit ContourTraceDeltaStore(int numNodes = 0) : addSpans(static_cast<std::size_t>(numNodes)), removeSpans(static_cast<std::size_t>(numNodes)) {}

    /**
     * @brief Returns read-only local boundary-edge additions for node.
     *
     * @param node Node identifier used by the operation.
     * @return Read-only local boundary-edge additions for `node`.
     *
     */
    std::span<const int> additions(NodeId node) const {
        const Span& span = addSpans[static_cast<std::size_t>(node)];
        if (span.size == 0) {
            return {};
        }
        return std::span<const int>(addValues.data() + span.offset, static_cast<std::size_t>(span.size));
    }

    /**
     * @brief Returns read-only local boundary-edge removals for node.
     *
     * @param node Node identifier used by the operation.
     * @return Read-only local boundary-edge removals for `node`.
     *
     */
    std::span<const int> removals(NodeId node) const {
        const Span& span = removeSpans[static_cast<std::size_t>(node)];
        if (span.size == 0) {
            return {};
        }
        return std::span<const int>(removeValues.data() + span.offset, static_cast<std::size_t>(span.size));
    }

    /**
     * @brief Builds the persistent compact store from transient edge lists.
     *
     * @param additions Values added to the local representation.
     * @param removals Values removed from the local representation.
     * @param numPackedEdges Number represented by `numPackedEdges`.
     * @return The resulting persistent compact store from transient edge lists.
     */
    static ContourTraceDeltaStore fromPendingEdgeLists(const PendingPixelLists& additions, const PendingPixelLists& removals, int numPackedEdges) {
        const int numNodes = additions.numLists();
        ContourTraceDeltaStore store(numNodes);
        store.addValues.reserve(additions.entryCount());
        store.removeValues.reserve(removals.entryCount());
        std::vector<uint16_t> edgeMark(static_cast<std::size_t>(numPackedEdges), 0);
        uint16_t markGeneration = 1;

        for (NodeId node = 0; node < numNodes; ++node) {
            store.appendCompacted(additions, node, store.addValues, store.addSpans[node], edgeMark, markGeneration);
            store.appendCompacted(removals, node, store.removeValues, store.removeSpans[node], edgeMark, markGeneration);
        }
        return store;
    }

  private:
    /**
     * @brief Appends one node contribution to the compact delta storage.
     *
     * @param lists Per-node lists that store the accumulated delta entries.
     * @param node Node identifier used by the operation.
     * @param values Values read or written by the operation.
     * @param span Per-node span metadata written for the compact storage.
     * @param edgeMark Generation-mark buffer used to avoid revisiting edges.
     * @param markGeneration Active mark generation.
     */
    static void appendCompacted(const PendingPixelLists& lists, NodeId node, std::vector<int>& values, Span& span, std::vector<uint16_t>& edgeMark,
                                uint16_t& markGeneration) {
        nextMarkGeneration(edgeMark, markGeneration);
        span.offset = checkedU32(values.size(), "contour trace delta offset");
        lists.appendUniqueValues(node, values, edgeMark, markGeneration);
        span.size = checkedU32(values.size() - span.offset, "contour trace delta size");
    }

    /**
     * @brief Advances mark generation.
     *
     * @param edgeMark Generation-mark buffer used to avoid revisiting edges.
     * @param markGeneration Active mark generation.
     */
    static void nextMarkGeneration(std::vector<uint16_t>& edgeMark, uint16_t& markGeneration) {
        ++markGeneration;
        if (markGeneration == 0) {
            std::fill(edgeMark.begin(), edgeMark.end(), 0);
            markGeneration = 1;
        }
    }

    /**
     * @brief Checks and converts u32.
     *
     * @param value Value used by the operation.
     * @param context Operation name used in diagnostics.
     * @return Value converted to `uint32_t` after range validation.
     */
    static uint32_t checkedU32(std::size_t value, const char* context) {
        if (value > static_cast<std::size_t>(std::numeric_limits<uint32_t>::max())) {
            throw std::overflow_error(std::string(context) + " exceeds uint32_t.");
        }
        return static_cast<uint32_t>(value);
    }
};

} // namespace mmcfilters::detail
