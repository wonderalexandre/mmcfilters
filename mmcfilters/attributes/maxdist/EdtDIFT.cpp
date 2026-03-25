#include "EdtDIFT.hpp"

namespace mmcfilters
{
  namespace maxdist
  {
    #define SQUARE(x) ((x)*(x))
    #define MIN(a, b) (((a) < (b)) ? (a) : (b)) 
    #define MAX(a, b) (((a) > (b)) ? (a) : (b))

    // ---------------------------------------------------------------------------------
    // Adaptive Adjacency Neighbors Iterator
    // ---------------------------------------------------------------------------------
    AdaptiveAdj::Neighbors::Iterator::Iterator(const AdaptiveAdj::Neighbors &neighbors, 
      int idx)
      :neighbors_{neighbors}, idx_{idx}
    {}

    bool AdaptiveAdj::Neighbors::Iterator::operator==(const Iterator &other) const noexcept
    {
      return idx_ == other.idx_;
    }

    std::pair<Point2D, int> AdaptiveAdj::Neighbors::Iterator::operator*() const
    {
      Point2D p = neighbors_.point(idx_);
      int n = neighbors_.nextAdj(idx_);
      return std::make_pair(p, n);
    }

    // ---------------------------------------------------------------------------------------
    // Adaptive Adjacency Neighbors 
    // ---------------------------------------------------------------------------------------
    AdaptiveAdj::Neighbors::Neighbors(const Point2D &p, const AdaptiveAdj &adj)
      : p_{p}, adj_{adj}
    {}

    Point2D AdaptiveAdj::Neighbors::point(int idx) const
    {
      return p_ + adj_.offset_[idx];
    }

    AdaptiveAdj::Neighbors::Iterator AdaptiveAdj::Neighbors::begin() const
    {
      return Iterator(*this, 0);
    }

    AdaptiveAdj::Neighbors::Iterator AdaptiveAdj::Neighbors::end() const
    {
      return Iterator(*this, adj_.offset_.size());
    }

    // -----------------------------------------------------------------------------------------
    // Adaptive Adjacency 
    // ------------------------------------------------------------------------------------------
    AdaptiveAdj::AdaptiveAdj(std::vector<Point2D> offset, std::vector<int> nextAdj)
      : offset_(std::move(offset)), nextAdj_(std::move(nextAdj))
    {}
 
    AdaptiveAdj::AdaptiveAdj(std::vector<Point2D> &&offset, std::vector<int> nextAdj)
      : offset_(std::move(offset)), nextAdj_(std::move(nextAdj))
    {}

    AdaptiveAdj::AdaptiveAdj(std::vector<Point2D> offset, std::vector<int> &&nextAdj)
      : offset_(std::move(offset)), nextAdj_(std::move(nextAdj))
    {}

    AdaptiveAdj::AdaptiveAdj(std::vector<Point2D> &&offset, std::vector<int> &&nextAdj)
      : offset_(std::move(offset)), nextAdj_(std::move(nextAdj))
    {}

    AdaptiveAdj::AdaptiveAdj(std::vector<Point2D> offset, std::initializer_list<int> nextAdj)
       : offset_(std::move(offset)), nextAdj_(nextAdj)
    {}

    AdaptiveAdj::AdaptiveAdj(std::initializer_list<Point2D> offset, std::vector<int> nextAdj)
       : offset_(offset), nextAdj_(std::move(nextAdj))
    {}

    AdaptiveAdj::AdaptiveAdj(std::initializer_list<Point2D> offset, std::initializer_list<int> nextAdj)
      : offset_{offset}, nextAdj_{nextAdj}
    {}
    
    AdaptiveAdj::AdaptiveAdj(std::vector<Point2D>&& offset, std::initializer_list<int> nextAdj)
      : offset_(std::move(offset)), nextAdj_{nextAdj}
    {}

    AdaptiveAdj::AdaptiveAdj(std::initializer_list<Point2D> offset, std::vector<int>&& nextAdj)
      : offset_(offset), nextAdj_(nextAdj)
    {}

    AdaptiveAdj::Neighbors AdaptiveAdj::neighbors(const Point2D &p) const
    {
      return Neighbors(p, *this);
    }

    // --------------------------------------------------------------------------------------------------
    // AdaptiveAdjBank
    // --------------------------------------------------------------------------------------------------
    AdaptiveAdjBank::AdaptiveAdjBank()
    {
      bank_.reserve(9);

      // All:
      AdaptiveAdj adj1(
        { Point2D( 1, -1), Point2D( 1,  0), Point2D( 1,  1), Point2D(-1, -1),
          Point2D(-1,  0), Point2D(-1,  1), Point2D( 0, -1), Point2D( 0,  1)}, 
        { 5, 1, 6, 7, 2, 8, 3, 4});
      bank_.push_back(adj1);
      
      // Right:
      AdaptiveAdj adj2(
        { Point2D( 1, -1), Point2D( 1, 0), Point2D(1, 1) },
        {5, 1, 6});
      bank_.push_back(adj2);

      // Left:
      AdaptiveAdj adj3(
        { Point2D( -1, -1), Point2D( -1, 0), Point2D( -1, 1) },
        { 7, 2, 8});
      bank_.push_back(adj3);

      // Top:
      AdaptiveAdj adj4(
        { Point2D(-1, -1), Point2D(0, -1), Point2D(1, -1) },
        {7, 3, 5});
      bank_.push_back(adj4);

      // Bottom:
      AdaptiveAdj adj5(
        { Point2D(-1, 1), Point2D(0, 1), Point2D(1, 1) },
        {8, 4, 6});
      bank_.push_back(adj5);

      // Top right:
      AdaptiveAdj adj6(
        { Point2D( 0, -1), Point2D(1, -1), Point2D(1, 0), Point2D(-1, -1),
          Point2D(-1, -1), Point2D(1,  1)}, {4, 6, 1, 8, 5});
      bank_.push_back(adj6);

      // bottom right:
      AdaptiveAdj adj7(
        { Point2D( 0, 1), Point2D( 1, 1), Point2D( 1, 0), Point2D(-1, 1),
          Point2D( 1,-1) }, {4, 6, 1, 8, 5});
      bank_.push_back(adj7);

      // Top left:
      AdaptiveAdj adj8(
        { Point2D( 0, -1), Point2D(-1, -1), Point2D(-1, 0), Point2D( 1,-1),
          Point2D(-1, 1) }, { 3, 7, 2, 5, 8});
      bank_.push_back(adj8);

      // Bottom left:
      AdaptiveAdj adj9(
        { Point2D( 0, 1), Point2D(-1, 1), Point2D(-1, 0), Point2D(1, 1),
          Point2D(-1, -1) }, {4, 8, 2, 6, 7});
      bank_.push_back(adj9);
    }

    const AdaptiveAdj &AdaptiveAdjBank::adj(int idx) const
    {
      return bank_[idx];
    }

    const AdaptiveAdj &AdaptiveAdjBank::operator[](int idx) const
    {
      return adj(idx);
    }

    // ---------------------------------------------------------------------------------------------------
    // EdtDIFT
    // ---------------------------------------------------------------------------------------------------
    EdtDIFT::EdtDIFT(int nrows, int ncols) :
      bin_{nrows, ncols},
      root_{nrows, ncols},
      cost_{nrows, ncols},
      Bedt_{nrows, ncols},
      O_{nrows, ncols},
      adjMap_{nrows, ncols},
      adj4_{nrows, ncols, 1.0},
      Q_{nrows * ncols, SQUARE(static_cast<int>(MIN(ncols, nrows) / 2.0 + 1))},
      domain_{ncols, ncols}
    { 
      for (int pidx = 0; pidx < bin_.getSize(); ++pidx) {
        root_[pidx] = pidx;
      }

      setUpAdjMap();
    }

    void EdtDIFT::setUpAdjMap()
    {
      int lastCol = Bedt_.getNumCols() - 1;
      for (int i = 1; i < Bedt_.getNumRows(); i++) {
        adjMap_[i * Bedt_.getNumCols()] = 1;            // fill first column
        adjMap_[i * Bedt_.getNumCols() + lastCol] = 2;  // fill last column
      }

      int lastRow = Bedt_.getNumCols() * (Bedt_.getNumRows()-1);
      for (int i = 1; i < Bedt_.getNumCols()-1; i++) {
        adjMap_[i] = 4;             // fill first row
        adjMap_[lastRow + i] = 3;   // fill last row
      }

      adjMap_[0] = 6;
      adjMap_[Bedt_.getNumCols()-1] = 8;
      adjMap_[Bedt_.getNumCols() * (Bedt_.getNumRows()-1)] = 5;
      adjMap_[Bedt_.getNumCols() * (Bedt_.getNumRows()-1) + Bedt_.getNumCols() - 1] = 7;
    }

    void EdtDIFT::treeRemoval(const std::vector<int> &toRemove,
        std::vector<int> &stack)
    {
      int top = -1;
      for (int pidx : toRemove) {
        O_[pidx] = 1;
        cost_[pidx] = Q_.maxCost();
        Q_.reopen(pidx, cost_[pidx]);
        stack[++top] = pidx;
      }

      while (top > -1) {
        int pidx = stack[top];
        Point2D p = domain_.point(pidx);
        --top;

        const AdaptiveAdj &AA = AAB_[adjMap_[pidx]];
        for (const auto& [q, ai] : AA.neighbors(p)) {
          int qidx = domain_.index(q);
          if (cost_[root_[qidx]] == Q_.maxCost()) {
            if (O_[qidx] == 0) {
              O_[qidx] = 1;
              cost_[qidx] = Q_.maxCost();
              Q_.reopen(qidx, cost_[qidx]);
              stack[++top] = qidx;
            }
          }
          else if (bin_[qidx] > 0 && Q_.stateOf(qidx) != PQueue::State::ACTIVE) {
            Q_.reinsert(qidx, cost_[qidx]);
          }
        }
      }
    }

    void EdtDIFT::run()
    {
      while (!Q_.empty()) {
        int pidx = Q_.popMin();
        Point2D p = domain_.point(pidx);

        int ridx = root_[pidx];

        Point2D r = domain_.point(ridx);
        Bedt_[ridx] = MAX(Bedt_[ridx], cost_[pidx]);

        const AdaptiveAdj &AA = AAB_[adjMap_[pidx]];
        for (const auto &[q, ai] : AA.neighbors(p)) {
          int qidx = domain_.index(q);
          int dx = q.x() - r.x();
          int dy = q.y() - r.y();
          int tmp = SQUARE(dx) + SQUARE(dy);

          if (tmp < cost_[qidx] && O_[qidx] == 1) {
            if (Q_.stateOf(qidx) != PQueue::State::ACTIVE) {
              cost_[qidx] = tmp;
              Q_.insert(qidx, cost_[qidx]);
            }
            else {
              Q_.decreaseCost(qidx, tmp);
            }

            root_[qidx] = ridx;
            adjMap_[qidx] = ai;
          }
        }
      }
    }
    
    void EdtDIFT::seed(int pidx)
    {
      root_[pidx] = pidx;
      cost_[pidx] = 0;
      Q_.insert(pidx, cost_[pidx]);
    }

    void EdtDIFT::insertNeighborsPQueue(int pidx)
    {
      Point2D p = domain_.point(pidx);

      for (int qidx : adj4_.getNeighborPixels(pidx)) {
        if (bin_[qidx] > 0 && cost_[qidx] != Q_.maxCost() 
              && Q_.stateOf(qidx) != PQueue::State::ACTIVE) {
          Q_.reopen(qidx, cost_[qidx]);
        }
      }
    }
  }
}