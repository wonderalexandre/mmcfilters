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
      
        Neighbors(const Point2D &p, const AdaptiveAdj &adj, size_t end);

        Point2D point(int idx) const;
        inline int nextAdj(int idx) const { return adj_.nextAdj_[idx]; }
        inline Point2D operator()(int idx) const { return point(idx); }

        inline int size() const { return end_; }

        Iterator begin() const;
        Iterator end() const;

      private:
        Point2D p_;
        const AdaptiveAdj &adj_;
        size_t end_;
      };

      AdaptiveAdj(std::vector<Point2D> offset, std::vector<int> nextAdj, size_t npropagation);      
      AdaptiveAdj(std::vector<Point2D> &&offset, std::vector<int> nextAdj, size_t npropagation);
      AdaptiveAdj(std::vector<Point2D> offset, std::vector<int> &&nextAdj, size_t npropagation);
      AdaptiveAdj(std::vector<Point2D> &&offset, std::vector<int> &&nextAdj, size_t npropagation);

      AdaptiveAdj(std::vector<Point2D> offset, std::initializer_list<int> nextAdj, size_t npropagation);
      AdaptiveAdj(std::initializer_list<Point2D> offset, std::vector<int> nextAdj, size_t npropagation);
      AdaptiveAdj(std::initializer_list<Point2D> offset, std::initializer_list<int> nextAdj, size_t npropagation);
      
      AdaptiveAdj(std::vector<Point2D>&& offset, std::initializer_list<int> nextAdj, size_t npropagation);
      AdaptiveAdj(std::initializer_list<Point2D> offset, std::vector<int>&& nextAdj, size_t npropagation);

      Neighbors neighbors(const Point2D &p) const;
      Neighbors neighborsPropogation(const Point2D &p) const;

    private:
      std::vector<Point2D> offset_;
      std::vector<int> nextAdj_;
      size_t npropagation_;         // number of elements for propagation
    };

    class AdaptiveAdjBank
    {
    public:
      AdaptiveAdjBank();

      inline size_t size() const noexcept { return bank_.size(); }

      const AdaptiveAdj &adj(int idx) const;
      const AdaptiveAdj &operator[](int idx) const;

    private:
      std::vector<AdaptiveAdj> bank_;
    };

    class EdtDIFT
    {
    public:
      static const int NIL = -1;

      EdtDIFT(int nrows, int ncols);
      void run();

      inline void addPixelToBinaryImage(int pidx) { bin_[pidx] = 1; }

      void insertNeighborsPQueue(int pidx);
      void seed(int pidx);

      void open(int pidx);

      void treeRemoval(const std::vector<int> &toRemove);

      int maxBedt(const std::vector<int> &Ncontour) const;

      void saveDistanceTransform(const std::string &filename = "") const;
      void saveUnderlyingBinaryImage(const std::string &filename = "") const;
      void displayRootMapForSmallImages() const;

    private:
      void setUpAdjMap();

    private:
      ImageUInt8 bin_;
      ImageInt32 root_;
      ImageInt32 Bedt_;
      ImageUInt8 adjMap_;
      ImageUInt8 O_;

      PQueue Q_;
      AdjacencyRelation adj4_;
      AdaptiveAdjBank AAB_;
      Box2D domain_;
      std::vector<int> stack_;
    };
  }
}