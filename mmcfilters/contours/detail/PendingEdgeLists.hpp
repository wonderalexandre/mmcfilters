#pragma once

#include "../../utils/Common.hpp"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace mmcfilters::contours::detail {

/**
 * @brief Temporary packed-edge lists indexed by node.
 *
 * Each tree node owns one singly linked list. The contiguous entry buffer keeps
 * allocation overhead independent of the number of nodes that receive edges.
 */
class PendingEdgeLists {
  public:
    /**
     * @brief Creates one empty edge list for each node.
     * @param numNodes Number of internal node slots.
     * @param expectedNumEntries Optional capacity estimate for all entries.
     */
    explicit PendingEdgeLists(int numNodes, int expectedNumEntries = 0) : headByNode_(static_cast<std::size_t>(numNodes), -1) {
        if (expectedNumEntries > 0) {
            entries_.reserve(static_cast<std::size_t>(expectedNumEntries));
        }
    }

    /**
     * @brief Adds a packed contour edge to a node list.
     * @param node Destination node identifier.
     * @param packedEdge Packed contour edge.
     */
    void add(NodeId node, int packedEdge) {
        const int entryIndex = static_cast<int>(entries_.size());
        entries_.push_back(Entry{packedEdge, headByNode_[static_cast<std::size_t>(node)]});
        headByNode_[static_cast<std::size_t>(node)] = entryIndex;
    }

    /**
     * @brief Appends the distinct edges of one node to an output buffer.
     *
     * `edgeMarks` and `markGeneration` let callers reuse one deduplication
     * buffer across nodes.
     * @param node Source node identifier.
     * @param output Buffer that receives distinct packed edges.
     * @param edgeMarks Deduplication generation indexed by packed edge.
     * @param markGeneration Active generation in `edgeMarks`.
     */
    void appendDistinct(NodeId node, std::vector<int>& output, std::vector<uint16_t>& edgeMarks, uint16_t markGeneration) const {
        for (int entryIndex = headByNode_[static_cast<std::size_t>(node)]; entryIndex != -1;
             entryIndex = entries_[static_cast<std::size_t>(entryIndex)].nextIndex) {
            const int packedEdge = entries_[static_cast<std::size_t>(entryIndex)].packedEdge;
            assert(packedEdge >= 0 && packedEdge < static_cast<int>(edgeMarks.size()));
            if (edgeMarks[static_cast<std::size_t>(packedEdge)] != markGeneration) {
                edgeMarks[static_cast<std::size_t>(packedEdge)] = markGeneration;
                output.push_back(packedEdge);
            }
        }
    }

    /**
     * @brief Returns the number of node lists.
     * @return Number of indexed node slots.
     */
    [[nodiscard]] int numNodes() const noexcept { return static_cast<int>(headByNode_.size()); }

    /**
     * @brief Returns the number of stored list entries.
     * @return Total number of packed-edge entries.
     */
    [[nodiscard]] std::size_t numEntries() const noexcept { return entries_.size(); }

  private:
    /** @brief One packed edge and the next entry in its node list. */
    struct Entry {
        /// Packed contour edge.
        int packedEdge = -1;
        /// Previous list head, or -1 at the end.
        int nextIndex = -1;
    };

    std::vector<Entry> entries_; ///< Contiguous storage shared by all node lists.
    std::vector<int> headByNode_; ///< Entry index at the head of each node list.
};

} // namespace mmcfilters::contours::detail
