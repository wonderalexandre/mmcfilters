#pragma once 
#include "utils/Image.hpp"
#include "utils/AdjacencyRelation.hpp"
#include "PQueue.hpp"

#include "Geometry.hpp"


namespace mmcfilters
{
  namespace maxdist
  {    
    class AdaptiveAdj
    {
    public:
      class Neighbors
      {
      public:
        class Iterator
        {
        public:
          Iterator(const Neighbors &neighbors, int idx = 0);

          std::pair<Point2D, int>  operator*() const;
          inline Iterator& operator++() noexcept { idx_++; return *this; }

          bool operator==(const Iterator &other) const noexcept;
          inline bool operator!=(Iterator &other) const noexcept {  return !(*this == other); }

        private:
          int idx_;
          const Neighbors &neighbors_;
        };        
      
        Neighbors(const Point2D &p, const AdaptiveAdj &adj);

        Point2D point(int idx) const;
        inline int nextAdj(int idx) const { return adj_.nextAdj_[idx]; }
        inline Point2D operator()(int idx) const { return point(idx); }

        inline int size() const noexcept { return adj_.offset_.size(); }

        Iterator begin() const;
        Iterator end() const;

      private:
        const Point2D &p_;
        const AdaptiveAdj &adj_;
      };

      AdaptiveAdj(std::vector<Point2D> offset, std::vector<int> nextAdj);      
      AdaptiveAdj(std::vector<Point2D> &&offset, std::vector<int> nextAdj);
      AdaptiveAdj(std::vector<Point2D> offset, std::vector<int> &&nextAdj);
      AdaptiveAdj(std::vector<Point2D> &&offset, std::vector<int> &&nextAdj);

      AdaptiveAdj(std::vector<Point2D> offset, std::initializer_list<int> nextAdj);
      AdaptiveAdj(std::initializer_list<Point2D> offset, std::vector<int> nextAdj);
      AdaptiveAdj(std::initializer_list<Point2D> offset, std::initializer_list<int> nextAdj);
      
      AdaptiveAdj(std::vector<Point2D>&& offset, std::initializer_list<int> nextAdj);
      AdaptiveAdj(std::initializer_list<Point2D> offset, std::vector<int>&& nextAdj);

      Neighbors neighbors(const Point2D &p) const;
    
    private:
      std::vector<Point2D> offset_;
      std::vector<int> nextAdj_;
    };

    class AdaptiveAdjBank
    {
    public:
      AdaptiveAdjBank();

      const AdaptiveAdj &adj(int idx) const;
      const AdaptiveAdj &operator[](int idx) const;

    private:
      std::vector<AdaptiveAdj> bank_;
    };

    // class EdtDIFT
    // {
    // public:
    //   static const int NIL = -1;

    //   EdtDIFT(int nrows, int ncols);
    //   void run();

    // private:
    //   void setUpAdjMap();

    // private:
    //   ImageUInt8 bin_;
    //   ImageInt32 root_;
    //   ImageInt32 pred_;
    //   ImageInt32 cost_;
    //   ImageInt32 Bedt_;
    //   ImageUInt8 adjMap_;
    //   ImageUInt8 O_;

    //   PQueue Q_;
    //   AdjacencyRelation adj8_;
    // };
  }
}