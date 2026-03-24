#include "PQueue.hpp"

namespace mmcfilters
{
  namespace maxdist
  {
    PQueue::PQueue(int numPixels, int maxCost):
      buckets_(maxCost + 1),
      cost_(numPixels, -1),
      next_(numPixels, -1),
      prev_(numPixels, -1),
      state_(numPixels, State::ABSENT),
      numPixels_(numPixels),
      maxCost_(maxCost),
      scan_(0),
      size_(0)
    {}

    void PQueue::insert(int pidx, int cost)
    {
      assert(state_[pidx] == State::ABSENT);
      assert(cost >= 0 && cost <= maxCost_);
      cost_[pidx] = cost;
      state_[pidx] = State::ACTIVE;
      pushFront(buckets_[cost], pidx);
      ++size_;
    }

    void PQueue::decreaseCost(int pidx, int newCost)
    {
      assert(state_[pidx] == State::ACTIVE);
      assert(newCost >= 0 && newCost <= cost_[pidx]);
      removeFromBucket(cost_[pidx], pidx);
      cost_[pidx] = newCost;
      pushFront(buckets_[newCost], pidx);
      // pidx is still in the queue.
    }

    // Remove and return the pixel with the minimum cost.
    // Returns -1 if the queue is empty.
    int PQueue::popMin() 
    {
      if (empty())
        return -1;
      
      advanceScan();
      int pidx = buckets_[scan_];
      removeFromBucket(scan_, pidx);
      state_[pidx] = State::DONE;
      --size_;
      return pidx;
    }
    
    std::optional<int> PQueue::minCost() const 
    {
      if (empty()) 
        return std::nullopt;

      int s = scan_;
      while (s <= maxCost_ && buckets_[s] == -1)
        ++s;
      
      return (s <= maxCost_) ? std::optional<int>(s) : std::nullopt;
    }    
    
    void PQueue::reset()
    {
      std::fill(cost_.begin(), cost_.end(), -1);
      std::fill(next_.begin(), next_.end(), -1);
      std::fill(prev_.begin(), prev_.end(), -1);
      std::fill(state_.begin(), state_.end(), State::ABSENT);
      std::fill(buckets_.begin(), buckets_.end(), -1);
      scan_ = 0;
      size_ = 0;
    }

    void PQueue::pushFront(int &head, int pidx) 
    {
      next_[pidx] = head;  
      prev_[pidx] = -1;
      if (head != -1)
        prev_[head] = pidx;
      
      head = pidx;
    }

    void PQueue::removeFromBucket(int cost, int pidx)
    {
      int &head = buckets_[cost];
      if (prev_[pidx] != -1) 
        next_[prev_[pidx]] = next_[pidx];    // update next pointer of the previous node
      else  
        head = next_[pidx];  //pidx is the head, so update the new head.
      
      if (next_[pidx] != -1)  
        prev_[next_[pidx]] = prev_[pidx];  // update prev pointer of the next node
      
      // remove next and previous pointers of pidx
      next_[pidx] = prev_[pidx] = -1;
    }

    void PQueue::advanceScan()
    {
      while (scan_ <= maxCost_ && buckets_[scan_] == -1)
        ++scan_;
    }
  }
}