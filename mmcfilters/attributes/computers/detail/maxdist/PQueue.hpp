#pragma once

#include <cassert>
#include <vector>
#include <algorithm>
#include <limits>

namespace mmcfilters::attributes::computers::detail::maxdist {
/**
 * @brief Bucket priority queue for non-negative integer costs.
 *
 * PQueue stores one doubly linked list per cost bucket and supports FIFO or
 * LIFO removal from the current minimum/maximum non-empty bucket. The
 * implementation is specialized for the squared-distance labels used by
 * EdtDIFT: costs must be finite integers in `[0, maxCost]`, where
 * `maxCost` is the `nbuckets` value passed to the constructor.
 *
 * The queue does not own element payloads. Elements are dense integer ids in
 * `[0, nelems)`, and each element has an externally visible cost and state.
 * The caller must avoid inserting the same element twice without removing
 * or updating it first.
 */
class PQueue // sPQueue
{
  public:
    /**
     * @brief Sentinel cost larger than any bucketed finite distance.
     */
    inline static constexpr int PINF = std::numeric_limits<int>::max();

    /**
     * @brief Sentinel cost smaller than any bucketed finite distance.
     */
    inline static constexpr int NINF = std::numeric_limits<int>::min();

    /**
     * @brief Sentinel index for absent pixels and empty bucket links.
     */
    inline static constexpr int NIL = -1;

    /**
     * @brief Per-element queue state used by the dynamic distance transform.
     */
    enum class State {
        /**
         * Element is not currently stored in any bucket.
         */
        NOT_PROCESSED,

        /**
         * Element is present in exactly one bucket.
         */
        QUEUED,

        /**
         * Element was removed from a bucket and its current label is final for
         * the active propagation pass.
         */
        POPPED
    };

    // ------------------------------------------------------------------------------------------
    // Public Methods
    // ------------------------------------------------------------------------------------------
    /**
     * @brief Creates a queue with buckets for every cost in `[0, nbuckets]`.
     *
     * `nbuckets` is an inclusive maximum finite cost, not a bucket count.
     * `nelems` defines the dense element id domain. Both values must be
     * non-negative.
     *
     * @param nbuckets Largest queue bucket index.
     * @param nelems Number of addressable queue elements.
     *
     * Invalid dimensions violate this internal kernel precondition and are
     * diagnosed only by debug assertions.
     */
    PQueue(int nbuckets, int nelems) : nadded_{0} {
        assert(nbuckets >= 0 && nelems >= 0);
        pixels_.nelems = nelems;
        pixels_.cost.resize(nelems);
        pixels_.elem.resize(nelems);

        std::fill(pixels_.cost.begin(), pixels_.cost.end(), 0);

        buckets_.nbuckets = nbuckets + 1;
        buckets_.first.resize(static_cast<std::size_t>(buckets_.nbuckets));
        buckets_.last.resize(static_cast<std::size_t>(buckets_.nbuckets));
        buckets_.maxvalue = NINF;
        buckets_.minvalue = PINF;

        for (int i = 0; i < buckets_.nbuckets; i++) {
            buckets_.first[i] = NIL;
            buckets_.last[i] = NIL;
        }

        for (int pidx = 0; pidx < pixels_.nelems; pidx++) {
            pixels_.elem[pidx].next = NIL;
            pixels_.elem[pidx].prev = NIL;
            pixels_.elem[pidx].state = State::NOT_PROCESSED;
        }
    }

    /**
     * @brief True when no element is currently queued.
     *
     * @return True when no element is currently queued; otherwise false.
     */
    inline bool isEmpty() const noexcept { return nadded_ == 0; }

    /**
     * @brief True when every element id in the configured domain is queued.
     *
     * @return True when every element id in the configured domain is queued; otherwise false.
     */
    inline bool isFull() const noexcept { return nadded_ == pixels_.nelems; }

    /**
     * @brief Updates the stored cost label of an element.
     *
     * This method only writes the cost array. If the element is already
     * queued, use update() so the bucket links remain consistent.
     *
     * @param elem Element identifier addressed by the queue.
     * @param newcost Replacement queue cost.
     */
    inline void setCost(int elem, int newcost) { pixels_.cost[elem] = newcost; }

    /**
     * @brief Updates the externally visible state of an element.
     *
     * @param elem Element identifier addressed by the queue.
     * @param state State read or updated by the operation.
     */
    inline void setState(int elem, State state) { pixels_.elem[elem].state = state; }

    /**
     * @brief Mutable access to an element cost label.
     *
     * @param elem Element identifier addressed by the queue.
     * @return Reference to the resulting object.
     */
    int& cost(int elem) { return pixels_.cost[elem]; }

    /**
     * @brief Read-only access to an element cost label.
     *
     * @param elem Element identifier addressed by the queue.
     * @return Read-only access to an element cost label.
     */
    int cost(int elem) const { return pixels_.cost[elem]; }

    /**
     * @brief Mutable access to all cost labels.
     *
     * Changing costs of queued elements can invalidate bucket membership. The
     * current MAX_DIST code uses this primarily for visualization/readback.
     *
     * @return Reference to the resulting object.
     */
    std::vector<int>& cost() { return pixels_.cost; }

    /**
     * @brief Read-only access to all cost labels.
     *
     * @return Reference to the resulting object.
     */
    const std::vector<int>& cost() const { return pixels_.cost; }

    /**
     * @brief Mutable access to an element state.
     *
     * @param elem Element identifier addressed by the queue.
     * @return Reference to the resulting object.
     */
    State& state(int elem) { return pixels_.elem[elem].state; }

    /**
     * @brief Read-only access to an element state.
     *
     * @param elem Element identifier addressed by the queue.
     * @return Read-only access to an element state.
     */
    State state(int elem) const { return pixels_.elem[elem].state; }

    /**
     * @brief Inserts an element in the bucket matching its current cost.
     *
     * @param elem Element identifier addressed by the queue.
     *
     * An out-of-domain cost violates the established EdtDIFT invariant and is
     * diagnosed only by a debug assertion.
     */
    void insert(int elem) {
        int bucket = pixels_.cost[elem];
        assert(bucket >= 0 && bucket < buckets_.nbuckets);

        ++nadded_;

        if (bucket < buckets_.minvalue)
            buckets_.minvalue = bucket;
        if (bucket > buckets_.maxvalue)
            buckets_.maxvalue = bucket;

        if (buckets_.first[bucket] == NIL) {
            buckets_.first[bucket] = elem;
            pixels_.elem[elem].prev = NIL;
        } else {
            pixels_.elem[buckets_.last[bucket]].next = elem;
            pixels_.elem[elem].prev = buckets_.last[bucket];
        }

        buckets_.last[bucket] = elem;
        pixels_.elem[elem].next = NIL;
        pixels_.elem[elem].state = State::QUEUED;
    }

    /**
     * @brief Removes an element from its current bucket.
     *
     * The method assumes `elem` is currently queued. It updates neighbouring
     * links and bucket head/tail pointers, then marks the element as POPPED.
     *
     * @param elem Element identifier addressed by the queue.
     */
    void remove(int elem) {
        --nadded_;
        int bucket = pixels_.cost[elem];
        int prev = pixels_.elem[elem].prev;
        int next = pixels_.elem[elem].next;

        if (buckets_.first[bucket] == elem) {
            buckets_.first[bucket] = next;
            if (next == NIL)
                buckets_.last[bucket] = NIL;
            else
                pixels_.elem[next].prev = NIL;
        } else {
            pixels_.elem[prev].next = next;
            if (next == NIL)
                buckets_.last[bucket] = prev;
            else
                pixels_.elem[next].prev = prev;
        }
        pixels_.elem[elem].state = State::POPPED;
    }

    /**
     * @brief Changes a queued element to a new cost bucket.
     *
     * This is the safe operation for distance relaxation: it removes the
     * element from the old bucket, changes the cost, and reinserts it into the
     * matching bucket.
     *
     * @param elem Element identifier addressed by the queue.
     * @param newcost Replacement queue cost.
     */
    void update(int elem, int newcost) {
        remove(elem);
        pixels_.cost[elem] = newcost;
        insert(elem);
    }

    /**
     * @brief Returns the largest non-empty bucket cost.
     *
     * @return The largest non-empty bucket cost.
     *
     */
    int maxValue() {
        requireNotEmpty();
        return findMaxBucket();
    }

    /**
     * @brief Returns the smallest non-empty bucket cost.
     *
     * @return The smallest non-empty bucket cost.
     *
     */
    int minValue() {
        requireNotEmpty();
        return findMinBucket();
    }

    /**
     * @brief Returns the oldest element in the smallest non-empty bucket.
     *
     * @return The oldest element in the smallest non-empty bucket.
     *
     */
    int minElemFIFO() {
        requireNotEmpty();
        int bucket = findMinBucket();
        return buckets_.first[bucket];
    }

    /**
     * @brief Pops the oldest element from the smallest non-empty bucket.
     *
     * @return The removed oldest element from the smallest non-empty bucket.
     *
     */
    int popMinFIFO() {
        requireNotEmpty();
        --nadded_;
        int bucket = findMinBucket();
        return bucketFIFO(bucket);
    }

    /**
     * @brief Pops the oldest element from the largest non-empty bucket.
     *
     * @return The removed oldest element from the largest non-empty bucket.
     *
     */
    int popMaxFIFO() {
        requireNotEmpty();
        --nadded_;
        int bucket = findMaxBucket();
        return bucketFIFO(bucket);
    }

    /**
     * @brief Pops the newest element from the smallest non-empty bucket.
     *
     * @return The removed newest element from the smallest non-empty bucket.
     *
     */
    int popMinLIFO() {
        requireNotEmpty();
        --nadded_;
        int bucket = findMinBucket();
        return bucketLIFO(bucket);
    }

    /**
     * @brief Pops the newest element from the largest non-empty bucket.
     *
     * @return The removed newest element from the largest non-empty bucket.
     *
     */
    int popMaxLIFO() {
        requireNotEmpty();
        --nadded_;
        int bucket = findMaxBucket();
        return bucketLIFO(bucket);
    }

  private:
    // ---------------------------------------------------------------------------------------------
    // Internal structures
    // ---------------------------------------------------------------------------------------------
    /**
     * @brief Intrusive doubly linked-list node for one queued element.
     */
    struct PixelListNode // sPNode
    {
        /**
         * @brief Next element id in the same bucket, or NIL.
         */
        int next;

        /**
         * @brief Previous element id in the same bucket, or NIL.
         */
        int prev;

        /**
         * @brief Current lifecycle state of the element.
         */
        State state;
    };

    /**
     * @brief Dense per-element storage: links, state, and current cost label.
     */
    struct PixelList // sPQDoublyLinkedLists
    {
        /**
         * @brief Intrusive link/state storage indexed by element id.
         */
        std::vector<PixelListNode> elem;

        /**
         * @brief Number of valid dense element ids.
         */
        int nelems;

        /**
         * @brief Cost label indexed by element id.
         */
        std::vector<int> cost;
    };

    /**
     * @brief Per-cost bucket heads/tails plus cached search bounds.
     */
    struct BucketList {
        /**
         * @brief Oldest element id per bucket, or NIL when the bucket is empty.
         */
        std::vector<int> first;

        /**
         * @brief Newest element id per bucket, or NIL when the bucket is empty.
         */
        std::vector<int> last;

        /**
         * @brief Number of allocated buckets, equal to max finite cost plus one.
         */
        int nbuckets;

        /**
         * @brief Cached lower bound for the next non-empty minimum-bucket scan.
         */
        int minvalue;

        /**
         * @brief Cached upper bound for the next non-empty maximum-bucket scan.
         */
        int maxvalue;
    };

    // ---------------------------------------------------------------------------
    //  Attributes (object members)
    // ---------------------------------------------------------------------------
    /**
     * @brief Number of elements currently present in buckets.
     */
    int nadded_;

    /**
     * @brief Bucket heads/tails and min/max scan caches.
     */
    BucketList buckets_;

    /**
     * @brief Dense per-element links, states, and cost labels.
     */
    PixelList pixels_;

    // ---------------------------------------------------------------------------
    // Private Methods
    // ---------------------------------------------------------------------------
    /**
     * @brief Rejects public min/max/pop operations on an empty queue.
     */
    void requireNotEmpty() const { assert(!isEmpty()); }

    /**
     * @brief Removes and returns the oldest element in `bucket`.
     *
     * The caller must guarantee that `bucket` is non-empty and has already
     * decremented `nadded_`.
     *
     * @param bucket Bucket read or updated by the operation.
     * @return The removed oldest element in bucket.
     */
    int bucketFIFO(int bucket) {
        int elem = buckets_.first[bucket];
        int next = pixels_.elem[elem].next;

        if (next == NIL) {
            buckets_.first[bucket] = NIL;
            buckets_.last[bucket] = NIL;
        } else {
            buckets_.first[bucket] = next;
            pixels_.elem[next].prev = NIL;
        }
        pixels_.elem[elem].state = State::POPPED;
        return elem;
    }

    /**
     * @brief Removes and returns the newest element in `bucket`.
     *
     * The caller must guarantee that `bucket` is non-empty and has already
     * decremented `nadded_`.
     *
     * @param bucket Bucket read or updated by the operation.
     * @return The removed newest element in bucket.
     */
    int bucketLIFO(int bucket) {
        int elem = buckets_.last[bucket];
        int prev = pixels_.elem[elem].prev;

        if (prev == NIL) {
            buckets_.first[bucket] = NIL;
            buckets_.last[bucket] = NIL;
        } else {
            buckets_.last[bucket] = prev;
            pixels_.elem[prev].next = NIL;
        }
        pixels_.elem[elem].state = State::POPPED;
        return elem;
    }

    /**
     * @brief Finds the smallest non-empty bucket and refreshes the cache.
     *
     * The queue must be non-empty. Public callers enforce this through
     * requireNotEmpty(); keeping the scan unchecked avoids duplicating the
     * branch inside the hot path.
     *
     * @return The located smallest non-empty bucket and refreshes the cache.
     */
    int findMinBucket() {
        int current = buckets_.minvalue;
        if (buckets_.first[current] == NIL) {
            do {
                ++current;
            } while ((current < buckets_.nbuckets) && (buckets_.first[current] == NIL));

            if (current < buckets_.nbuckets)
                buckets_.minvalue = current;
            else
                return NIL;
        }
        return current;
    }

    /**
     * @brief Finds the largest non-empty bucket and refreshes the cache.
     *
     * The queue must be non-empty. Public callers enforce this through
     * requireNotEmpty().
     *
     * @return The located largest non-empty bucket and refreshes the cache.
     */
    int findMaxBucket() {
        int current = buckets_.maxvalue;
        if (buckets_.first[current] == NIL) {
            do {
                --current;
            } while ((current >= 0) && (buckets_.first[current] == NIL));

            if (current >= 0)
                buckets_.maxvalue = current;
            else
                return NIL;
        }

        return current;
    }
};
} // namespace mmcfilters::attributes::computers::detail::maxdist
