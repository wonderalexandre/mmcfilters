#pragma once

#include "../../utils/Common.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace mmcfilters::detail {

/**
 * @brief Temporary per-node pixel lists backed by a reusable contiguous buffer.
 *
 * This store is intentionally small and transient. Extraction uses it for many
 * append/consume operations, then `ContourDeltaStore` compacts the surviving
 * per-node values into persistent read-only spans.
 *
 * The structure owns one singly linked list per node. Freed entries are recycled
 * through `freeHead_`, so memory follows the peak number of simultaneously live
 * entries rather than the total number of list operations.
 */
struct PendingPixelLists {
    struct Entry {
        int value;
        int next;
    };

    /// Builds a store for `numNodes` lists with an optional capacity hint.
    explicit PendingPixelLists(int numNodes, int capacityHint = 0): head_(numNodes, -1) {
        if (capacityHint > 0) {
            entries_.reserve(capacityHint);
        }
    }

    /// Inserts `value` into the list owned by `node`.
    void add(NodeId node, int value) {
        const int slot = allocate();
        entries_[slot] = Entry{value, head_[node]};
        head_[node] = slot;
    }

    /// Moves every element of `node` into `out` while recycling the used slots.
    void consumeInto(NodeId node, std::vector<int>& out) {
        int idx = head_[node];
        while (idx != -1) {
            const int next = entries_[idx].next;
            out.push_back(entries_[idx].value);
            recycle(idx);
            idx = next;
        }
        head_[node] = -1;
    }

    /**
     * @brief Appends unique, in-domain values from `node` using caller-owned marks.
     *
     * The mark buffer is reused by the caller across nodes. `markGeneration`
     * identifies the current logical set and avoids sorting or allocating a
     * per-node temporary set.
     */
    void appendUniqueValues(
        NodeId node,
        std::vector<int>& out,
        std::vector<uint16_t>& pixelMark,
        uint16_t markGeneration) const {
        for (int idx = head_[node]; idx != -1; idx = entries_[idx].next) {
            const int value = entries_[idx].value;
            if (value < 0 || value >= static_cast<int>(pixelMark.size())) {
                continue;
            }
            if (pixelMark[static_cast<std::size_t>(value)] != markGeneration) {
                pixelMark[static_cast<std::size_t>(value)] = markGeneration;
                out.push_back(value);
            }
        }
    }

    /// @return Number of node-owned lists.
    int numLists() const {
        return static_cast<int>(head_.size());
    }

    /// @return Number of arena slots currently allocated, including free slots.
    std::size_t entryCount() const {
        return entries_.size();
    }

private:
    /// Returns a free slot, either from the free list or via `push_back`.
    int allocate() {
        if (freeHead_ == -1) {
            entries_.push_back(Entry{0, -1});
            return static_cast<int>(entries_.size() - 1);
        }
        const int idx = freeHead_;
        freeHead_ = entries_[idx].next;
        return idx;
    }

    /// Returns an index to the free list.
    void recycle(int idx) {
        entries_[idx].next = freeHead_;
        freeHead_ = idx;
    }

    std::vector<Entry> entries_;
    std::vector<int> head_;
    int freeHead_ = -1;
};

} // namespace mmcfilters::detail
