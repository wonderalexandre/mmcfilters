#pragma once

#include "PendingPixelLists.hpp"

namespace mmcfilters::detail {

/**
 * @brief Compact per-node delta representation for contour additions/removals.
 *
 * The extraction pass is mutation-heavy, so it uses `PendingPixelLists`. Once the
 * pass finishes, the persistent result only needs read-only per-node spans.
 *
 * Values are stored in two flat arrays. `addSpans[node]` and
 * `removeSpans[node]` identify each node slice, following the same idea as a CSR
 * representation. Per-node duplicates are removed during compaction.
 */
struct ContourDeltaStore {
    /**
     * @brief Offset and length of one node slice in a flat value vector.
     */
    struct Span {
        uint32_t offset = 0;
        uint32_t size = 0;
    };

    /// Concatenated compacted local contour additions.
    std::vector<int> addValues;
    /// Concatenated compacted local contour removals.
    std::vector<int> removeValues;
    /// Per-node spans into `addValues`.
    std::vector<Span> addSpans;
    /// Per-node spans into `removeValues`.
    std::vector<Span> removeSpans;

    explicit ContourDeltaStore(int numNodes = 0)
        : addSpans(static_cast<std::size_t>(numNodes)),
          removeSpans(static_cast<std::size_t>(numNodes)) {}

    /// @return Read-only local contour additions for `node`.
    std::span<const int> additions(NodeId node) const {
        const Span& span = addSpans[static_cast<std::size_t>(node)];
        if (span.size == 0) {
            return {};
        }
        return std::span<const int>(
            addValues.data() + span.offset,
            static_cast<std::size_t>(span.size));
    }

    /// @return Read-only local contour removals for `node`.
    std::span<const int> removals(NodeId node) const {
        const Span& span = removeSpans[static_cast<std::size_t>(node)];
        if (span.size == 0) {
            return {};
        }
        return std::span<const int>(
            removeValues.data() + span.offset,
            static_cast<std::size_t>(span.size));
    }

    /**
     * @brief Builds the persistent compact store from transient extraction lists.
     *
     * Each node list is appended once and deduplicated with generation marks.
     */
    static ContourDeltaStore fromPendingPixelLists(
        const PendingPixelLists& contours,
        const PendingPixelLists& removals,
        int numPixels) {
        const int numNodes = contours.numLists();
        ContourDeltaStore store(numNodes);
        store.addValues.reserve(contours.entryCount());
        store.removeValues.reserve(removals.entryCount());
        std::vector<uint16_t> pixelMark(static_cast<std::size_t>(numPixels), 0);
        uint16_t markGeneration = 1;

        for (NodeId node = 0; node < numNodes; ++node) {
            store.appendCompacted(contours, node, store.addValues, store.addSpans[node], pixelMark, markGeneration);
            store.appendCompacted(removals, node, store.removeValues, store.removeSpans[node], pixelMark, markGeneration);
        }
        return store;
    }

private:
    static void appendCompacted(
        const PendingPixelLists& lists,
        NodeId node,
        std::vector<int>& values,
        Span& span,
        std::vector<uint16_t>& pixelMark,
        uint16_t& markGeneration) {
        nextMarkGeneration(pixelMark, markGeneration);
        span.offset = checkedU32(values.size(), "contour delta offset");
        lists.appendUniqueValues(node, values, pixelMark, markGeneration);
        span.size = checkedU32(values.size() - span.offset, "contour delta size");
    }

    static void nextMarkGeneration(std::vector<uint16_t>& pixelMark, uint16_t& markGeneration) {
        ++markGeneration;
        if (markGeneration == 0) {
            std::fill(pixelMark.begin(), pixelMark.end(), 0);
            markGeneration = 1;
        }
    }

    static uint32_t checkedU32(std::size_t value, const char* context) {
        if (value > static_cast<std::size_t>(std::numeric_limits<uint32_t>::max())) {
            throw std::overflow_error(std::string(context) + " exceeds uint32_t.");
        }
        return static_cast<uint32_t>(value);
    }
};

} // namespace mmcfilters::detail
