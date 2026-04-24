#pragma once 

#include <vector>
#include <algorithm>
#include <limits>


namespace mmcfilters
{
  namespace maxdist
  {
    class PQueue   // sPQueue 
    {
    public:
      inline static constexpr int PINF = std::numeric_limits<int>::max();   // positive infinity cost
      inline static constexpr int NINF = std::numeric_limits<int>::min();   // negative infinity cost
      inline static constexpr int NIL = -1;                                 // NULL value for indices (pixels and buckets)

      enum class State {
        NOT_PROCESSED,     // ABSENT = WHITE
        QUEUED,            // ACTIVE = GRAY
        POPPED             // DONE   = BLACK   
      };

      // ------------------------------------------------------------------------------------------
      // Public Methods
      // ------------------------------------------------------------------------------------------
      PQueue(int nbuckets, int nelems)
        : nadded_{0}
      {
        pixels_.nelems = nelems;
        pixels_.cost.resize(nelems);
        pixels_.elem.resize(nelems);

        std::fill(pixels_.cost.begin(), pixels_.cost.end(), 0);

        buckets_.nbuckets = nbuckets;
        buckets_.first.resize(nbuckets+1);
        buckets_.last.resize(nbuckets+1);
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

      inline bool isEmpty() const noexcept { return nadded_ == 0; }
      inline bool isFull() const noexcept { return nadded_ == pixels_.nelems; }

      inline void setCost(int elem, int newcost) { pixels_.cost[elem] = newcost; }
      inline void setState(int elem, State state) { pixels_.elem[elem].state = state; }

      int& cost(int elem) { return pixels_.cost[elem]; }
      int cost(int elem) const { return pixels_.cost[elem]; }
      std::vector<int>& cost() { return pixels_.cost; }
      const std::vector<int>& cost() const { return pixels_.cost; }

      State& state(int elem) { return pixels_.elem[elem].state; }
      State state(int elem) const { return pixels_.elem[elem].state; }

      void insert(int elem) {
        ++nadded_;
        int bucket = pixels_.cost[elem];

        if (bucket < buckets_.minvalue)
          buckets_.minvalue = bucket;
        if (bucket > buckets_.maxvalue)
          buckets_.maxvalue = bucket;

        if (buckets_.first[bucket] == NIL) {
          buckets_.first[bucket] = elem;
          pixels_.elem[elem].prev = NIL;
        }
        else {
          pixels_.elem[buckets_.last[bucket]].next = elem;
          pixels_.elem[elem].prev = buckets_.last[bucket];
        }

        buckets_.last[bucket] = elem;
        pixels_.elem[elem].next = NIL;
        pixels_.elem[elem].state = State::QUEUED;
      }

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
        }
        else {
          pixels_.elem[prev].next = next;
          if (next == NIL)
            buckets_.last[bucket] = prev;
          else
            pixels_.elem[next].prev = prev;
        }
        pixels_.elem[elem].state = State::POPPED;
      }

      void update(int elem, int newcost) {
        remove(elem);
        pixels_.cost[elem] = newcost;
        insert(elem);
      }

      int maxValue() {
        return findMaxBucket();
      }

      int minValue() {
        return findMinBucket();
      }

      int minElemFIFO() {
        int bucket = findMinBucket();
        return buckets_.first[bucket];
      }

      int popMinFIFO() {
        --nadded_;
        int bucket = findMinBucket();
        return bucketFIFO(bucket);
      }

      int popMaxFIFO() {
        --nadded_;
        int bucket = findMaxBucket();
        return bucketFIFO(bucket);
      }

      int popMinLIFO() {
        --nadded_;
        int bucket = findMinBucket();
        return bucketLIFO(bucket);
      }

      int popMaxLIFO() {
        --nadded_;
        int bucket = findMaxBucket();
        return bucketLIFO(bucket);
      }

    private:
      // ---------------------------------------------------------------------------------------------
      // Internal structures
      // ---------------------------------------------------------------------------------------------
      struct PixelListNode  // sPNode
      {
        int next;
        int prev;
        State state; 
      };

      struct PixelList   // sPQDoublyLinkedLists 
      {
        std::vector<PixelListNode> elem;
        int nelems;
        std::vector<int> cost;
      }; 

      struct BucketList 
      {
        std::vector<int> first;
        std::vector<int> last;
        int nbuckets;
        int minvalue;
        int maxvalue;
      };

      // ---------------------------------------------------------------------------
      //  Attributes (object members)
      // ---------------------------------------------------------------------------
      int nadded_;
      BucketList buckets_;
      PixelList pixels_;

      // ---------------------------------------------------------------------------
      // Private Methods
      // ---------------------------------------------------------------------------
      int bucketFIFO(int bucket) {
        int elem = buckets_.first[bucket];
        int next = pixels_.elem[elem].next;

        if (next == NIL) {
          buckets_.first[bucket] = NIL;
          buckets_.last[bucket] = NIL;
        }
        else {
          buckets_.first[bucket] = next;
          pixels_.elem[elem].prev = NIL;
        }
        pixels_.elem[elem].state = State::POPPED;
        return elem;
      }

      int bucketLIFO(int bucket) {
        int elem = buckets_.last[bucket];
        int prev = pixels_.elem[elem].prev;

        if (prev == NIL) {
          buckets_.last[bucket] = NIL;
          buckets_.last[bucket] = NIL;
        }
        else {
          buckets_.last[bucket] = prev;
          pixels_.elem[elem].next = NIL;
        }
        pixels_.elem[elem].state = State::POPPED;
        return elem;
      }
      
      int findMinBucket() {
        int current = buckets_.minvalue;
        if (buckets_.first[current] == NIL) {
          do
          {
            ++current;
          } while ((current < buckets_.nbuckets) && (buckets_.first[current] == NIL));

          if (current < buckets_.nbuckets)
            buckets_.minvalue = current;
          else
            return NIL;
        }
        return current;
      }

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
  }
}
