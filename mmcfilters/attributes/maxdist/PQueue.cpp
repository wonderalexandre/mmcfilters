#include "PQueue.hpp"
#include <limits>

namespace mmcfilters
{
  namespace maxdist
  {
    const int PQueue::PINF = std::numeric_limits<int>::max();   
    const int PQueue::NINF = std::numeric_limits<int>::min();   
    const int PQueue::NIL = -1;

    // nbuckets = maximum cost and nelems = number of pixels
    PQueue::PQueue(int nbuckets, int nelems)
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

    void PQueue::insert(int elem)
    {
      ++nadded_;
      int bucket = pixels_.cost[elem];
      
      // Update min and max value if necessary
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

    void PQueue::remove(int elem)
    {
      --nadded_;
      int bucket = pixels_.cost[elem];
      int prev = pixels_.elem[elem].prev;
      int next = pixels_.elem[elem].next;

      // elem is the head of the "bucket"
      if (buckets_.first[bucket] == elem) {
        buckets_.first[bucket] = next;
        if (next == NIL) // elem is also the last element
          buckets_.last[bucket] = NIL;
        else
          pixels_.elem[next].prev = NIL;
      }
      else {   // elem is either in the middle or in the last item in the bucket
        pixels_.elem[prev].next = next;
        if (next == NIL)  // if elem is the last element
          buckets_.last[bucket] = prev;
        else
          pixels_.elem[next].prev = prev;
      }
      pixels_.elem[elem].state = State::POPPED;
    }

    void PQueue::update(int elem, int newcost)
    {
      remove(elem);
      pixels_.cost[elem] = newcost;
      insert(elem);
    }

    int PQueue::maxValue()
    {
      return findMaxBucket();
    }

    int PQueue::minValue()
    {
      return findMinBucket();
    }

    int PQueue::minElemFIFO()
    {
      int bucket = findMinBucket();
      return buckets_.first[bucket];
    }

    int PQueue::popMinFIFO()
    {
      --nadded_;
      int bucket = findMinBucket();
      return bucketFIFO(bucket);
    }

    int PQueue::popMaxFIFO()
    {
      --nadded_;
      int bucket = findMaxBucket();
      return bucketFIFO(bucket);
    }

    int PQueue::popMinLIFO()
    {
      --nadded_;
      int bucket = findMinBucket();
      return bucketLIFO(bucket);
    }

    int PQueue::popMaxLIFO()
    {
      --nadded_;
      int bucket = findMaxBucket();
      return bucketLIFO(bucket);
    }

    int PQueue::bucketFIFO(int bucket)
    {
      int elem = buckets_.first[bucket];
      int next = pixels_.elem[elem].next;
      
      if (next == NIL) { // there was a single element in the list
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

    int PQueue::bucketLIFO(int bucket)
    {
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
    
    int PQueue::findMinBucket()
    {
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

    int PQueue::findMaxBucket()
    {
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
  }
}
