#pragma once 
#include <vector>
#include <optional>
#include <stdexcept>
#include <cassert>

// ----------------------------------------------------------------------------
// IFT Bucket Queue (Dial's Algorithm)
//
// A priority queue for integer costs in [0, max_cost].
// Supports O(1) insert, decrease-key, and amortised O(1) pop-min.
// Storage; O(n + C), where C = max_cost
//
// Pixels are identified by integer IDs in [0, num_pixels).
// Each pixel can be in one of three states:
//    - Not yet inserted    (state = ABSENT)
//    - In the queue        (state = ACTIVE)
//    - Already removed     (state = DONE)
// ----------------------------------------------------------------------------
namespace mmcfilters
{
  namespace maxdist
  {
    class PQueue
    {
    public:
      enum class State : u_int8_t { ABSENT, ACTIVE, DONE };

      // Construct the queue fpr 'numPixels' pixels and costs in [0, maxCost]
      PQueue(int numPixels, int maxCost);

      // Insert pixel pidx with the given cost. pidx must be ABSENT.
      void insert(int pidx, int cost);

      // Lower the cost of pixel p (already ACTIVE). New cost must be <= old cost.
      void decreaseCost(int pidx, int newCost);

      // Remove and return the pixel with the minimum cost.
      // Returns -1 if the queue is empty
      int popMin();
      
      // Peak at the minimum cost without removing anything
      std::optional<int> minCost() const;

  
    inline bool empty() const noexcept { return size_ == 0; }
    inline int costOf(int pidx) const noexcept { return cost_[pidx]; }
    inline stateOf(int pidx) const noexcept { return state_[pidx]; }
    
    void reset();


    private:
      // Prepend pixel pidx to the bucket's intrusive list
      void pushFront(int &head, int pidx);

      // Unlink pixel pidx from the given bucket's list.
      void removeFromBucket(int cost, int pidx);

      // Advance the scan pointer to the next non-empty bucket.
      // Because path costs are non-decreasing along optimal paths (the
      // IFT monotone condition), we never need to look back.
      void advanceScan();

    private:
      // Each bucket is the head of an intrusive doubly-linked list.
      // Sentinel value -1 means "no pixel".
      std::vector<int> buckets_;  // buckets_[cost] = head pixel list
      std::vector<int> cost_;     // current cost of each pixel 
      std::vector<int> next_;     // next pixel in the same bucket list
      std::vector<int> prev_;     // prev pixel in the same bucket list
      
      std::vector<State> state_;
      int numPixels_;
      int maxCost_;
      int scan_;                  // monotone scan pointer - never decremented
      int size_;
    };
  }
}