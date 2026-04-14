#pragma once 

#include <vector>


namespace mmcfilters
{
  namespace maxdist
  {
    class PQueue   // sPQueue 
    {
    public:
      static const int PINF;   // positive infinity cost
      static const int NINF;   // negative infinity cost
      static const int NIL;    // NULL value for indices (pixels and buckets)

      enum class State {
        NOT_PROCESSED,     // ABSENT = WHITE
        QUEUED,            // ACTIVE = GRAY
        POPPED             // DONE   = BLACK   
      };

      // ------------------------------------------------------------------------------------------
      // Public Methods
      // ------------------------------------------------------------------------------------------
      PQueue(int nbuckets, int nelems);

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

      void insert(int elem);
      void remove(int elem);
      void update(int elem, int newcost);

      int maxValue();
      int minValue();

      int minElemFIFO();

      int popMinFIFO();
      int popMaxFIFO();
      int popMinLIFO();
      int popMaxLIFO();

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
      int bucketFIFO(int bucket);
      int bucketLIFO(int bucket);
      
      int findMinBucket();
      int findMaxBucket();
    };
  }
}
