#include <iostream>
#include "PQueue.hpp"
#include <cassert>
#include <algorithm>
#include <limits>

namespace mmcfilters
{
  namespace maxdist
  {
    const int PQueue::INF = std::numeric_limits<int>::max();

    PQueue::PQueue(int numPixels, int initialCapacity)
      : cost_(numPixels, -1),
        next_(numPixels, -1),
        prev_(numPixels, -1),
        state_(numPixels, State::ABSENT),
        numPixels_{numPixels},
        infHead_{-1},
        scan_{0},
        size_{0},
        finiteSize_{0}
    {
      int cap = 1;
      while (cap < initialCapacity) cap <<= 1;
      buckets_.assign(cap, -1);
      mask_ = cap - 1;
    }

    void PQueue::insert(int pidx, int cost)
    {
      assert(state_[pidx] == State::ABSENT);
      assert(cost >= 0);
      cost_[pidx] = cost;
      state_[pidx] = State::ACTIVE;
      if (cost == INF) 
        pushFront(infHead_, pidx);
      else {
        ensureCapacity(cost);
        pushFront(buckets_[slot(cost)], pidx);
        ++finiteSize_;
      }
      ++size_;
    }

    void PQueue::decreaseCost(int pidx, int newCost)
    {
      assert(state_[pidx] == State::ACTIVE);
      assert(newCost >= 0 && newCost <= cost_[pidx]);
      
      // Remove from current location.
      if (cost_[pidx] == INF)
        removeFromList(infHead_, pidx);
      else {
        removeFromList(buckets_[slot(cost_[pidx])], pidx);
        --finiteSize_;
      }
      
      // Insert into new location.
      cost_[pidx] = newCost;
      if (newCost == INF) 
        pushFront(infHead_, pidx);
      else {
        ensureCapacity(newCost);
        pushFront(buckets_[slot(newCost)], pidx);
        ++finiteSize_;
      }      
     // size_ unchanged - pixel stays in the queue.
    }

    void PQueue::reinsert(int pidx, int cost)
    {
      assert(state_[pidx] == State::DONE);
      assert(cost >= 0);

      cost_[pidx] = cost;
      state_[pidx] = State::ACTIVE;

      if (cost == INF) 
        pushFront(infHead_, pidx);
      else {
        ensureCapacity(cost);
        pushFront(buckets_[slot(cost)], pidx);
        ++finiteSize_;
      }

      ++size_;
    }

    void PQueue::reopen(int pidx, int cost)
    {
      // assert(state_[pidx] == State::DONE);
      assert(cost >= 0);
      cost_[pidx] = cost;
      state_[pidx] = State::ABSENT;
    }

    int PQueue::popMin()
    {
      if (empty()) 
        return -1;

      int pidx;
      if (finiteSize_ > 0) {
        advanceScan();
        int &head = buckets_[scan_ & mask_];
        pidx = head;
        removeFromList(head, pidx);
        --finiteSize_;
      }
      else {
        pidx = infHead_;
        removeFromList(infHead_, pidx);
      }
      
      state_[pidx] = State::DONE;
      --size_;

      return pidx;
    }

    std::optional<int> PQueue::minCost() const
    {
      if (empty()) 
        return std::nullopt;
      if (finiteSize_ == 0) 
        return INF;

      int cap = static_cast<int>(buckets_.size());
      for (int i = 0; i < cap; ++i) {
        if (buckets_[(scan_ + i) & mask_] != -1)
          return scan_ + i;        
      }
      return INF; // unreacheable if finiteSize > 0, but satisfies the compiler
    }

    void PQueue::reset()
    {
      std::fill(cost_.begin(), cost_.end(), -1);
      std::fill(next_.begin(), next_.end(), -1);
      std::fill(prev_.begin(), prev_.end(), -1);
      std::fill(state_.begin(), state_.end(), State::ABSENT);
      std::fill(buckets_.begin(), buckets_.end(), -1);
      infHead_ = -1;
      scan_ = 0;
      size_ = 0;
      finiteSize_ = 0;
    }

    void PQueue::ensureCapacity(int cost) 
    {
      while (cost - scan_ >= static_cast<int>(buckets_.size()))
        grow();
    }

    void PQueue::grow()
    {
      int newCap = static_cast<int>(buckets_.size());
      int newMask = newCap - 1;

      std::vector<int> newBuckets(newCap, -1);

      for (int pidx = 0; pidx < numPixels_; pidx++) {
        if (state_[pidx] != State::ACTIVE || cost_[pidx] == INF)
          continue;
        
        int s = cost_[pidx] & newMask;
        next_[pidx] = newBuckets[s];
        prev_[pidx] = -1;
        if (newBuckets[s] != -1)
          prev_[newBuckets[s]] = pidx;
        newBuckets[s] = pidx;
      }

      buckets_ = std::move(newBuckets);
      mask_ = newMask;
    }

    void PQueue::pushFront(int &head, int pidx)
    {
      next_[pidx] = head;
      prev_[pidx] = -1;
      if (head != -1)
        prev_[head] = pidx;
      head = pidx;
    }

    void PQueue::removeFromList(int &head, int pidx)
    {
      if (prev_[pidx] != -1) 
        next_[prev_[pidx]] = next_[pidx];
      else 
        head = next_[pidx];

      if (next_[pidx] != -1) 
        prev_[next_[pidx]] = prev_[pidx];
      next_[pidx] = prev_[pidx] = -1;
    }

    void PQueue::advanceScan() 
    {
      while (buckets_[scan_ & mask_] == -1) 
        ++scan_;
    }
  }
}