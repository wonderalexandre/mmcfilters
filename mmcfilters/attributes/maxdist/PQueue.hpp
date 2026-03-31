#pragma once 
#include <vector>
#include <optional>
#include <cassert>
#include <algorithm>
#include <cstdint>

// ------------------------------------------------------------------------------------------
// IFT / DIFT Bucket Queue - Falcão's circular + growing Dial implementation
//
// Two structural improvements over the plain linear Dial queue:
// 1. CIRCULAR - the bucket array of capacity C is addressed as physical slot = cost % C
//    The scan pointer advances modulo C and wraps around, so the same C slots are reused 
//    across the entire algorithm. No slot is ever "past" the pointer: reinsert() with an 
//    arbitrarily low cost drops the pixel into the correct slot with no rewind needed.
//
// 2. GROWING - if a new cost exceeds the current capacity, the array doubles until it fits, 
//    and all ACTIVE pixels are rehashed into the new layout. The queue can therefore start 
//    small and expand lazily - useful when the true cost range is unknown at construction time,
//    or when DIFT seeds introduce costs that exceed the original IFT range.
//
//  3. INFINITY - PQueue::INF (= INT_MAX) is a first-class cost. INF pixels live in a dedicated 
//     overflow list that sits entirely outside the circular array - they never touch slot 
//     arithmentic or trigger a grow. All finite-cost pixels are always popped before any INF
//     pixel. decrease_key() / reinsert() from INF to a finite value moves the pixel into 
//     the correct circular slot ,growing if needed.
//  
// Correctness invariant (circular aliasing):
//  capacity > maxActiveCost - scan_ at all times.
//  ensureCapacity() grows before any insertion that would break it.
//
// Complexity:
//   insert / decrease_key / reinsert / reopen       - O(1) amortised
//   pop_min                                         - O(1) amortised
//   grow (rare, triggered by ensure_capacity)       - O(n + C_new)
//   storage                                         - O(n + C)

namespace mmcfilters {
  namespace maxdist
  {
    class PQueue
    {
      public:
        // Sentinel for "infinite" cost. Use freely in every public method - the queue 
        // routes INF pixels to a separate overflow list.
        static const int INF;

        enum class State : uint8_t { ABSENT, ACTIVE, DONE };

        // numPixels      : total number of pixels in the image.
        // initialCapacity: starting number of buckets. Does not need to be equal to the true 
        //                   max cost - the queue grows automatically. Rounded up to the next power
        //                   of two internally so modulo reduces to a cheap bitmask operation.
        PQueue(int numPixels, int initialCapacity = 256);

        // -------------- Core Operations -----------------------------------------------------------

        // Insert pixel pidx with the given cost (maybe INFO). pidx must be ABSENT 
        void insert(int pidx, int cost);

        // Lower the cost of an ACTIVE pixel. newCost must be <= old cost.
        // Lowering from INF to a finite value is explicitly supported.
        void decreaseCost(int pidx, int newCost);

        // Re-enqueue a DONE pixel as new seed (DIFT). Cost may be INF
        // No scan rewind needed - circular addressing handles any cost value.
        void reinsert(int pdix, int cost);

        // Reset a DONE pixel to ABSENT without enqueuing it (DIFT). Cost may be INF
        // The propagation loop will insert() when a neighbor relaxes it
        void reopen(int pdix, int cost);

        // Remove and return the pixel with the minimum cost.
        // Finite-cost pixels are always returned before INF-cost pixels.
        // Returns -1 if the queue is empty.
        int popMin();

        // Peek at the minimum cost without removing anything
        // Returns INF if only infinite-cost pixels remain, nullopt if empty.
        std::optional<int> minCost() const;

        // ------- Accessors ----------------------------------------------------------------
        inline bool empty() const noexcept { return size_ == 0; }
        inline int size() const noexcept { return size_; }
        inline int costOf(int pidx) const noexcept { return cost_[pidx]; }
        inline State stateOf(int pidx) const noexcept { return state_[pidx]; }
        int capacity() const noexcept { return static_cast<int>(buckets_.size()); }

        void reset();

      private:
        std::vector<int> buckets_;   // circular bucket array; size always a power of 2
        std::vector<int> cost_;      
        std::vector<int> next_;
        std::vector<int> prev_;
        std::vector<State> state_;
        int numPixels_;
        int infHead_;          // head of overflow list for INF-cost pixels 
        int mask_;             // buckets_.size() - 1 - bitmask for fast modulo
        int scan_;             // logical scan position; physical slot = scan _ & mask_
        int size_;   
        int finiteSize_;       // ACTIVE pixels with finite cost - O(1) INF check


        // Map a cost value to its physical slot in the circular array.
        inline int slot(int cost) const noexcept { return cost & mask_; }

        // Ensure capacity > (cost - scan_) so no two live finite costs alias
        void ensureCapacity(int cost);

        void grow();

        void pushFront(int& head, int pidx);
        void removeFromList(int& head, int pidx);

        // Walk scan_ forward until its slot is non-empty.
        // The circular invariant guarantees all ACTIVE pixels liew within
        // [scan_, scan_ + capacity + 1], so this terminates in O(C) steps 
        // amortised over all popMin() calls.
        void advanceScan();
    };
  }
}





