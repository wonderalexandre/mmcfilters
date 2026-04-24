#pragma once

#include <iostream>
#include <iomanip>
#include <utility>
#include <vector>
#include <string>

#include "utils/Image.hpp"
#include "utils/AdjacencyRelation.hpp"
#include "PQueue.hpp"
#include "Geometry.hpp"


namespace mmcfilters
{
  namespace maxdist
  {
    namespace detail
    {
      inline int square(int value) noexcept
      {
        return value * value;
      }

      inline void printUInt8Image(const mmcfilters::ImageUInt8Ptr& image)
      {
        for (int row = 0; row < image->getNumRows(); ++row) {
          for (int col = 0; col < image->getNumCols(); ++col) {
            std::cout << std::setw(4) << static_cast<int>((*image)[row * image->getNumCols() + col]);
          }
          std::cout << "\n";
        }
      }
    }

    class AdaptiveAdj
    {
    public:
      class Neighbors
      {
      public:
        class Iterator
        {
        public:
          Iterator(const Neighbors &neighbors, int idx = 0)
            : idx_{idx}, neighbors_{neighbors}
          {}

          std::pair<Point2D, int>  operator*() const
          {
            Point2D p = neighbors_.point(idx_);
            int n = neighbors_.nextAdj(idx_);
            return std::make_pair(p, n);
          }
          inline Iterator& operator++() noexcept { idx_++; return *this; }

          bool operator==(const Iterator &other) const noexcept
          {
            return idx_ == other.idx_;
          }
          inline bool operator!=(Iterator &other) const noexcept {  return !(*this == other); }

        private:
          int idx_;
          const Neighbors &neighbors_;
        };        
      
        Neighbors(const Point2D &p, const AdaptiveAdj &adj, size_t end)
          : p_{p}, adj_{adj}, end_{end}
        {}

        Point2D point(int idx) const
        {
          return p_ + adj_.offset_[idx];
        }
        inline int nextAdj(int idx) const { return adj_.nextAdj_[idx]; }
        inline Point2D operator()(int idx) const { return point(idx); }

        inline int size() const { return end_; }

        Iterator begin() const
        {
          return Iterator(*this, 0);
        }

        Iterator end() const
        {
          return Iterator(*this, end_);
        }

      private:
        Point2D p_;
        const AdaptiveAdj &adj_;
        size_t end_;
      };

      AdaptiveAdj(std::vector<Point2D> offset, std::vector<int> nextAdj, size_t npropagation)
        : offset_(std::move(offset)), nextAdj_(std::move(nextAdj)), npropagation_{npropagation}
      {}
      AdaptiveAdj(std::vector<Point2D> &&offset, std::vector<int> nextAdj, size_t npropagation)
        : offset_(std::move(offset)), nextAdj_(std::move(nextAdj)), npropagation_{npropagation}
      {}
      AdaptiveAdj(std::vector<Point2D> offset, std::vector<int> &&nextAdj, size_t npropagation)
        : offset_(std::move(offset)), nextAdj_(std::move(nextAdj)), npropagation_{npropagation}
      {}
      AdaptiveAdj(std::vector<Point2D> &&offset, std::vector<int> &&nextAdj, size_t npropagation)
        : offset_(std::move(offset)), nextAdj_(std::move(nextAdj)), npropagation_{npropagation}
      {}

      AdaptiveAdj(std::vector<Point2D> offset, std::initializer_list<int> nextAdj, size_t npropagation)
        : offset_(std::move(offset)), nextAdj_(nextAdj), npropagation_{npropagation}
      {}
      AdaptiveAdj(std::initializer_list<Point2D> offset, std::vector<int> nextAdj, size_t npropagation)
        : offset_(offset), nextAdj_(std::move(nextAdj)), npropagation_{npropagation}
      {}
      AdaptiveAdj(std::initializer_list<Point2D> offset, std::initializer_list<int> nextAdj, size_t npropagation)
        : offset_{offset}, nextAdj_{nextAdj}, npropagation_{npropagation}
      {}
      
      AdaptiveAdj(std::vector<Point2D>&& offset, std::initializer_list<int> nextAdj, size_t npropagation)
        : offset_(std::move(offset)), nextAdj_{nextAdj}, npropagation_{npropagation}
      {}
      AdaptiveAdj(std::initializer_list<Point2D> offset, std::vector<int>&& nextAdj, size_t npropagation)
        : offset_(offset), nextAdj_(std::move(nextAdj)), npropagation_{npropagation}
      {}

      Neighbors neighbors(const Point2D &p) const
      {
        return Neighbors(p, *this, offset_.size());
      }

      Neighbors neighborsPropogation(const Point2D &p) const
      {
        return Neighbors(p, *this, npropagation_);
      }

    private:
      std::vector<Point2D> offset_;
      std::vector<int> nextAdj_;
      size_t npropagation_;         // number of elements for propagation
    };

    class AdaptiveAdjBank
    {
    public:
      AdaptiveAdjBank()
      {
        AdaptiveAdj adj1(
          { Point2D( 1, -1), Point2D( 1,  0), Point2D( 1,  1), Point2D(-1, -1),
            Point2D(-1,  0), Point2D(-1,  1), Point2D( 0, -1), Point2D( 0,  1)},
          { 5, 1, 6, 7, 2, 8, 3, 4}, 8);
        bank_.push_back(adj1);

        AdaptiveAdj adj2(
          { Point2D( 1, -1), Point2D( 1, 0), Point2D(1, 1) },
          {5, 1, 6}, 3);
        bank_.push_back(adj2);

        AdaptiveAdj adj3(
          { Point2D( -1, -1), Point2D( -1, 0), Point2D( -1, 1) },
          { 7, 2, 8}, 3);
        bank_.push_back(adj3);

        AdaptiveAdj adj4(
          { Point2D(-1, -1), Point2D(0, -1), Point2D(1, -1) },
          {7, 3, 5}, 3);
        bank_.push_back(adj4);

        AdaptiveAdj adj5(
          { Point2D(-1, 1), Point2D(0, 1), Point2D(1, 1) },
          {8, 4, 6}, 3);
        bank_.push_back(adj5);

        AdaptiveAdj adj6(
          { Point2D( 0, -1), Point2D(1, -1), Point2D(1, 0), Point2D(-1, -1),
            Point2D(1,  1)}, {3, 5, 1, 7, 6}, 3);
        bank_.push_back(adj6);

        AdaptiveAdj adj7(
          { Point2D( 0, 1), Point2D( 1, 1), Point2D( 1, 0), Point2D(-1, 1),
            Point2D( 1,-1) }, {4, 6, 1, 8, 5}, 3);
        bank_.push_back(adj7);

        AdaptiveAdj adj8(
          { Point2D( 0, -1), Point2D(-1, -1), Point2D(-1, 0), Point2D( 1,-1),
            Point2D(-1, 1) }, { 3, 7, 2, 5, 8}, 3);
        bank_.push_back(adj8);

        AdaptiveAdj adj9(
          { Point2D( 0, 1), Point2D(-1, 1), Point2D(-1, 0), Point2D(1, 1),
            Point2D(-1, -1) }, {4, 8, 2, 6, 7}, 3);
        bank_.push_back(adj9);
      }

      inline size_t size() const noexcept { return bank_.size(); }

      const AdaptiveAdj &adj(int idx) const
      {
        return bank_[idx];
      }

      const AdaptiveAdj &operator[](int idx) const
      {
        return adj(idx);
      }

    private:
      std::vector<AdaptiveAdj> bank_;
    };

    class EdtDIFT
    {
    public:
      inline static constexpr int NIL = -1;

      EdtDIFT(int nrows, int ncols) :
        bin_{nrows, ncols},
        root_{nrows, ncols},
        Bedt_{nrows, ncols},
        adjMap_{nrows, ncols},
        O_{nrows, ncols},
        Q_{detail::square(static_cast<int>(std::min(ncols, nrows) / 2.0 + 1)), nrows * ncols},
        adj4_{nrows, ncols, 1.0},
        domain_{ncols, nrows},
        stack_(nrows * ncols)
      {
        bin_.fill(0);
        root_.fill(0);
        Bedt_.fill(0);
        O_.fill(0);
        adjMap_.fill(0);

        for (int pidx = 0; pidx < bin_.getSize(); ++pidx) {
          root_[pidx] = pidx;
        }

        setUpAdjMap();
      }

      void run()
      {
        while (!Q_.isEmpty()) {
          int pidx = Q_.popMinFIFO();
          Point2D p = domain_.point(pidx);
          O_[pidx] = 0;

          int ridx = root_[pidx];
          Point2D r = domain_.point(ridx);
          Bedt_[ridx] = std::max(Bedt_[ridx], Q_.cost(pidx));

          const AdaptiveAdj &AA = AAB_[adjMap_[pidx]];
          for (const auto &[q, ai] : AA.neighborsPropogation(p)) {
            int qidx = domain_.index(q);
            int dx = q.x() - r.x();
            int dy = q.y() - r.y();
            int tmp = detail::square(dx) + detail::square(dy);

            if (tmp < Q_.cost(qidx) && O_[qidx] == 1) {
              if (Q_.state(qidx) != PQueue::State::QUEUED) {
                Q_.setCost(qidx, tmp);
                Q_.insert(qidx);
              }
              else {
                Q_.update(qidx, tmp);
              }

              root_[qidx] = ridx;
              adjMap_[qidx] = ai;
            }
          }
        }
      }

      inline void addPixelToBinaryImage(int pidx) { bin_[pidx] = 1; }

      void insertNeighborsPQueue(int pidx)
      {
        for (int qidx : adj4_.getNeighborPixels(pidx)) {
          if (bin_[qidx] > 0 && Q_.cost(qidx) != PQueue::PINF && Q_.state(qidx) != PQueue::State::QUEUED) {
            Q_.setState(qidx, PQueue::State::NOT_PROCESSED);
            Q_.insert(qidx);
          }
        }
      }

      void seed(int pidx)
      {
        root_[pidx] = pidx;
        Q_.setCost(pidx, 0);
        Q_.insert(pidx);
      }

      void open(int pidx)
      {
        O_[pidx] = 1;
        Q_.setCost(pidx, PQueue::PINF);
      }

      void treeRemoval(const std::vector<int> &toRemove)
      {
        int top = -1;
        for (int pidx : toRemove) {
          O_[pidx] = 1;
          Q_.setCost(pidx, PQueue::PINF);
          Q_.setState(pidx, PQueue::State::NOT_PROCESSED);
          ++top;
          stack_[top] = pidx;
        }

        while (top > -1) {
          int pidx = stack_[top];
          Point2D p = domain_.point(pidx);
          --top;

          const AdaptiveAdj &AA = AAB_[adjMap_[pidx]];
          for (const auto& [q, ai] : AA.neighbors(p)) {
            int qidx = domain_.index(q);

            if (Q_.cost(root_[qidx]) == PQueue::PINF) {
              if (O_[qidx] == 0) {
                O_[qidx] = 1;
                Q_.setCost(qidx, PQueue::PINF);
                Q_.setState(qidx, PQueue::State::NOT_PROCESSED);
                ++top;
                stack_[top] = qidx;
              }
            }
            else if (bin_[qidx] > 0 && Q_.state(qidx) != PQueue::State::QUEUED) {
              Q_.setState(qidx, PQueue::State::NOT_PROCESSED);
              Q_.insert(qidx);
            }
          }
        }
      }

      int maxBedt(const std::vector<int> &Ncontour) const
      {
        int maxValue = 0;
        for (int pidx : Ncontour) {
          int d = Bedt_[pidx];
          if (d > maxValue)
            maxValue = d;
        }

        return maxValue;
      }

      ImageUInt8Ptr distanceTransformImage() const
      {
        const std::vector<int> &cost = Q_.cost();

        ImageUInt8Ptr d = ImageUInt8::create(bin_.getNumRows(), bin_.getNumCols());
        int maxCost = 0;
        for (int pidx = 0; pidx < bin_.getSize(); ++pidx) {
          if (maxCost < cost[pidx]) {
            maxCost = cost[pidx];
          }
        }

        if (maxCost == 0) {
          d->fill(0);
          return d;
        }

        for (int pidx = 0; pidx < bin_.getSize(); ++pidx) {
          (*d)[pidx] = static_cast<uint8_t>((static_cast<float>(cost[pidx]) / static_cast<float>(maxCost) * 255));
        }

        return d;
      }

      ImageUInt8Ptr binaryImageForVisualisation() const
      {
        ImageUInt8Ptr b = ImageUInt8::create(bin_.getNumRows(), bin_.getNumCols());
        for (int pidx = 0; pidx < b->getSize(); pidx++) {
          (*b)[pidx] = bin_[pidx] == 0 ? 255 : 0;
        }
        return b;
      }

      void printDistanceTransform() const
      {
        detail::printUInt8Image(distanceTransformImage());
        std::cout << "\n ---- Queue Cost --- \n";
      }

      void printUnderlyingBinaryImage() const
      {
        detail::printUInt8Image(binaryImageForVisualisation());
      }

      void displayRootMapForSmallImages() const
      {
        for (int row = 0; row < root_.getNumRows(); row++) {
          std::cout << std::setw(4) << row;
          for (int col = 0; col < root_.getNumCols(); col++) {
            std::cout << std::setw(4) << root_[row * root_.getNumCols() + col];
          }
          std::cout << "\n";
        }
      }

    private:
      void setUpAdjMap()
      {
        int lastCol = adjMap_.getNumCols() - 1;
        for (int i = 1; i < adjMap_.getNumRows(); i++) {
          adjMap_[i * adjMap_.getNumCols()] = 1;
          adjMap_[(i * adjMap_.getNumCols()) + lastCol] = 2;
        }

        int lastRow = adjMap_.getNumCols() * (adjMap_.getNumRows()-1);
        for (int i = 1; i < adjMap_.getNumCols()-1; i++) {
          adjMap_[i] = 4;
          adjMap_[lastRow + i] = 3;
        }

        adjMap_[0] = 6;
        adjMap_[adjMap_.getNumCols()-1] = 8;
        adjMap_[adjMap_.getNumCols() * (adjMap_.getNumRows()-1)] = 5;
        adjMap_[adjMap_.getNumCols() * (adjMap_.getNumRows()-1) + adjMap_.getNumCols() - 1] = 7;
      }

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
