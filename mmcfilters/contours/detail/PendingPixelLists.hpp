#pragma once

#include "../../utils/Common.hpp"

#include <cassert>
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
    /** @brief Stores one pending-pixel list entry. */
    struct Entry {
        /** @brief Stores the value. */
        int value;
        /** @brief Stores the next. */
        int next;
    };

    /**
     * @brief Builds a store for `numNodes` lists with an optional capacity hint.
     *
     * @param numNodes Number represented by `numNodes`.
     * @param capacityHint Initial storage reservation hint.
     */
    explicit PendingPixelLists(int numNodes, int capacityHint = 0) : head_(numNodes, -1) {
        if (capacityHint > 0) {
            entries_.reserve(capacityHint);
        }
    }

    /**
     * @brief Inserts `value` into the list owned by `node`.
     *
     * @param node Node identifier used by the operation.
     * @param value Value used by the operation.
     */
    void add(NodeId node, int value) {
        const int slot = allocate();
        entries_[slot] = Entry{value, head_[node]};
        head_[node] = slot;
    }

    /**
     * @brief Moves every element of `node` into `out` while recycling the used slots.
     *
     * @param node Node identifier used by the operation.
     * @param out Destination receiving the result.
     */
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
     *
     * @param node Node identifier used by the operation.
     * @param out Destination receiving the result.
     * @param pixelMark Generation marks used to deduplicate pixels.
     * @param markGeneration Active pixel-mark generation.
     */
    void appendUniqueValues(NodeId node, std::vector<int>& out, std::vector<uint16_t>& pixelMark, uint16_t markGeneration) const {
        for (int idx = head_[node]; idx != -1; idx = entries_[idx].next) {
            const int value = entries_[idx].value;
            assert(value >= 0 && value < static_cast<int>(pixelMark.size()));
            if (pixelMark[static_cast<std::size_t>(value)] != markGeneration) {
                pixelMark[static_cast<std::size_t>(value)] = markGeneration;
                out.push_back(value);
            }
        }
    }

    /**
     * @brief Returns number of node-owned lists.
     *
     * @return Number of node-owned lists.
     */
    int numLists() const { return static_cast<int>(head_.size()); }

    /**
     * @brief Returns number of arena slots currently allocated, including free slots.
     *
     * @return Number of arena slots currently allocated, including free slots.
     */
    std::size_t entryCount() const { return entries_.size(); }

  private:
    /**
     * @brief Returns a free slot, either from the free list or via `push_back`.
     *
     * @return A free slot, either from the free list or via push_back.
     */
    int allocate() {
        if (freeHead_ == -1) {
            entries_.push_back(Entry{0, -1});
            return static_cast<int>(entries_.size() - 1);
        }
        const int idx = freeHead_;
        freeHead_ = entries_[idx].next;
        return idx;
    }

    /**
     * @brief Returns an index to the free list.
     *
     * @param idx Zero-based index used by the operation.
     */
    void recycle(int idx) {
        entries_[idx].next = freeHead_;
        freeHead_ = idx;
    }

    /** @brief Stores the entries. */
    std::vector<Entry> entries_;
    /** @brief Stores the head. */
    std::vector<int> head_;
    /** @brief Stores the free head. */
    int freeHead_ = -1;
};

} // namespace mmcfilters::detail
