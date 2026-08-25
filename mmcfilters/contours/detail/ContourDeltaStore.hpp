#pragma once

#include "PendingPixelLists.hpp"
#include "../../trees/MorphologicalTree.hpp"

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
        /** @brief Offset. */
        uint32_t offset = 0;
        /** @brief Size. */
        uint32_t size = 0;
    };

    /// Concatenated compacted local contour additions.
    std::vector<PixelId> addValues;
    /// Concatenated compacted local contour removals.
    std::vector<PixelId> removeValues;
    /// Per-node spans into `addValues`.
    std::vector<Span> addSpans;
    /// Per-node spans into `removeValues`.
    std::vector<Span> removeSpans;

    /**
     * @brief Constructs `ContourDeltaStore` from the supplied inputs.
     *
     * @param numNodes Number of internal nodes.
     */
    explicit ContourDeltaStore(int numNodes = 0) : addSpans(static_cast<std::size_t>(numNodes)), removeSpans(static_cast<std::size_t>(numNodes)) {}

    /**
     * @brief Returns read-only local contour additions for node.
     *
     * @param node Node identifier.
     * @return Read-only local contour additions for `node`.
     *
     */
    std::span<const PixelId> additions(NodeId node) const {
        const Span& span = addSpans[static_cast<std::size_t>(node)];
        if (span.size == 0) {
            return {};
        }
        return std::span<const PixelId>(addValues.data() + span.offset, static_cast<std::size_t>(span.size));
    }

    /**
     * @brief Returns read-only local contour removals for node.
     *
     * @param node Node identifier.
     * @return Read-only local contour removals for `node`.
     *
     */
    std::span<const PixelId> removals(NodeId node) const {
        const Span& span = removeSpans[static_cast<std::size_t>(node)];
        if (span.size == 0) {
            return {};
        }
        return std::span<const PixelId>(removeValues.data() + span.offset, static_cast<std::size_t>(span.size));
    }

    /**
     * @brief Builds the persistent compact store from transient extraction lists.
     *
     * Each node list is appended once and deduplicated with generation marks.
     *
     * @param contours Contour data.
     * @param removals Values removed from the local representation.
     * @param numPixels Number.
     * @return The resulting persistent compact store from transient extraction lists.
     */
    static ContourDeltaStore fromPendingPixelLists(const PendingPixelLists& contours, const PendingPixelLists& removals, int numPixels) {
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

    /**
     * @brief Copies one semantic contour-event source into compact persistent storage.
     *
     * The event source must expose `numAdditions()`, `numRemovals()`,
     * `additions(node)`, and `removals(node)`. This keeps contour
     * materialization independent from the algorithm used to derive boundary
     * lifetimes.
     *
     * @param tree Tree whose live nodes index the event source.
     * @param events Semantic contour birth/stop events.
     * @return Compact per-node additions and removals.
     */
    template <typename EventSource>
    static ContourDeltaStore fromEventSource(const MorphologicalTree& tree, const EventSource& events) {
        ContourDeltaStore store(tree.numInternalNodeSlots());
        store.addValues.reserve(events.numAdditions());
        store.removeValues.reserve(events.numRemovals());

        for (NodeId node = 0; node < tree.numInternalNodeSlots(); ++node) {
            if (!tree.isAlive(node)) {
                continue;
            }
            store.appendSpan(events.additions(node), store.addValues, store.addSpans[static_cast<std::size_t>(node)]);
            store.appendSpan(events.removals(node), store.removeValues, store.removeSpans[static_cast<std::size_t>(node)]);
        }
        return store;
    }

  private:
    /** @brief Appends one already-compacted event span. */
    static void appendSpan(std::span<const PixelId> source, std::vector<PixelId>& values, Span& span) {
        span.offset = checkedU32(values.size(), "contour event offset");
        values.insert(values.end(), source.begin(), source.end());
        span.size = checkedU32(source.size(), "contour event size");
    }

    /**
     * @brief Appends one node contribution to the compact delta storage.
     *
     * @param lists Per-node lists that store the accumulated delta entries.
     * @param node Node identifier.
     * @param values Values read or written by the operation.
     * @param span Per-node span metadata written for the compact storage.
     * @param pixelMark Generation-mark buffer used to avoid revisiting pixels.
     * @param markGeneration Active mark generation.
     */
    static void appendCompacted(const PendingPixelLists& lists, NodeId node, std::vector<PixelId>& values, Span& span, std::vector<uint16_t>& pixelMark,
                                uint16_t& markGeneration) {
        nextMarkGeneration(pixelMark, markGeneration);
        span.offset = checkedU32(values.size(), "contour delta offset");
        lists.appendUniqueValues(node, values, pixelMark, markGeneration);
        span.size = checkedU32(values.size() - span.offset, "contour delta size");
    }

    /**
     * @brief Advances mark generation.
     *
     * @param pixelMark Generation-mark buffer used to avoid revisiting pixels.
     * @param markGeneration Active mark generation.
     */
    static void nextMarkGeneration(std::vector<uint16_t>& pixelMark, uint16_t& markGeneration) {
        ++markGeneration;
        if (markGeneration == 0) {
            std::fill(pixelMark.begin(), pixelMark.end(), 0);
            markGeneration = 1;
        }
    }

    /**
     * @brief Checks and converts u32.
     *
     * @param value Value.
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
