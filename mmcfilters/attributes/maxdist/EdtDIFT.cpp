#include "EdtDIFT.hpp"
#include "../../../tests/Tests.hpp"

#include "../../../tests/Tests.hpp"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../../../external/stb/stb_image_write.h"

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
    AdaptiveAdj::Neighbors::Neighbors(const Point2D &p, const AdaptiveAdj &adj, size_t end)
      : p_{p}, adj_{adj}, end_{end}
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
      return Iterator(*this, end_);
    }

    // -----------------------------------------------------------------------------------------
    // Adaptive Adjacency 
    // ------------------------------------------------------------------------------------------
    AdaptiveAdj::AdaptiveAdj(std::vector<Point2D> offset, std::vector<int> nextAdj, size_t npropagation)
      : offset_(std::move(offset)), nextAdj_(std::move(nextAdj)), npropagation_{npropagation}
    {}
 
    AdaptiveAdj::AdaptiveAdj(std::vector<Point2D> &&offset, std::vector<int> nextAdj, size_t npropagation)
      : offset_(std::move(offset)), nextAdj_(std::move(nextAdj)), npropagation_{npropagation}
    {}

    AdaptiveAdj::AdaptiveAdj(std::vector<Point2D> offset, std::vector<int> &&nextAdj, size_t npropagation)
      : offset_(std::move(offset)), nextAdj_(std::move(nextAdj)), npropagation_{npropagation}
    {}

    AdaptiveAdj::AdaptiveAdj(std::vector<Point2D> &&offset, std::vector<int> &&nextAdj, size_t npropagation)
      : offset_(std::move(offset)), nextAdj_(std::move(nextAdj)), npropagation_{npropagation}
    {}

    AdaptiveAdj::AdaptiveAdj(std::vector<Point2D> offset, std::initializer_list<int> nextAdj, size_t npropagation)
       : offset_(std::move(offset)), nextAdj_(nextAdj), npropagation_{npropagation}
    {}

    AdaptiveAdj::AdaptiveAdj(std::initializer_list<Point2D> offset, std::vector<int> nextAdj, size_t npropagation)
       : offset_(offset), nextAdj_(std::move(nextAdj)), npropagation_{npropagation}
    {}

    AdaptiveAdj::AdaptiveAdj(std::initializer_list<Point2D> offset, std::initializer_list<int> nextAdj, size_t npropagation)
      : offset_{offset}, nextAdj_{nextAdj}, npropagation_{npropagation}
    {}
    
    AdaptiveAdj::AdaptiveAdj(std::vector<Point2D>&& offset, std::initializer_list<int> nextAdj, size_t npropagation)
      : offset_(std::move(offset)), nextAdj_{nextAdj}, npropagation_{npropagation}
    {}

    AdaptiveAdj::AdaptiveAdj(std::initializer_list<Point2D> offset, std::vector<int>&& nextAdj, size_t npropagation)
      : offset_(offset), nextAdj_(std::move(nextAdj)), npropagation_{npropagation}
    {}

    AdaptiveAdj::Neighbors AdaptiveAdj::neighbors(const Point2D &p) const
    {
      return Neighbors(p, *this, offset_.size());
    }

    AdaptiveAdj::Neighbors AdaptiveAdj::neighborsPropogation(const Point2D &p) const
    {
      return Neighbors(p, *this, npropagation_);
    }

    // --------------------------------------------------------------------------------------------------
    // AdaptiveAdjBank
    // --------------------------------------------------------------------------------------------------
    AdaptiveAdjBank::AdaptiveAdjBank()
    {
      // bank_.reserve(9);

      // All:
      AdaptiveAdj adj1(
        { Point2D( 1, -1), Point2D( 1,  0), Point2D( 1,  1), Point2D(-1, -1),
          Point2D(-1,  0), Point2D(-1,  1), Point2D( 0, -1), Point2D( 0,  1)}, 
        { 5, 1, 6, 7, 2, 8, 3, 4}, 8);
      bank_.push_back(adj1);
      
      // Right:
      AdaptiveAdj adj2(
        { Point2D( 1, -1), Point2D( 1, 0), Point2D(1, 1) },
        {5, 1, 6}, 3);
      bank_.push_back(adj2);

      // Left:
      AdaptiveAdj adj3(
        { Point2D( -1, -1), Point2D( -1, 0), Point2D( -1, 1) },
        { 7, 2, 8}, 3);
      bank_.push_back(adj3);

      // Top:
      AdaptiveAdj adj4(
        { Point2D(-1, -1), Point2D(0, -1), Point2D(1, -1) },
        {7, 3, 5}, 3);
      bank_.push_back(adj4);

      // Bottom:
      AdaptiveAdj adj5(
        { Point2D(-1, 1), Point2D(0, 1), Point2D(1, 1) },
        {8, 4, 6}, 3);
      bank_.push_back(adj5);

      // Top right:
      AdaptiveAdj adj6(
        { Point2D( 0, -1), Point2D(1, -1), Point2D(1, 0), Point2D(-1, -1),
          Point2D(1,  1)}, {3, 5, 1, 7, 6}, 3);
      bank_.push_back(adj6);

      // bottom right:
      AdaptiveAdj adj7(
        { Point2D( 0, 1), Point2D( 1, 1), Point2D( 1, 0), Point2D(-1, 1),
          Point2D( 1,-1) }, {4, 6, 1, 8, 5}, 3);
      bank_.push_back(adj7);

      // Top left:
      AdaptiveAdj adj8(
        { Point2D( 0, -1), Point2D(-1, -1), Point2D(-1, 0), Point2D( 1,-1),
          Point2D(-1, 1) }, { 3, 7, 2, 5, 8}, 3);
      bank_.push_back(adj8);

      // Bottom left:
      AdaptiveAdj adj9(
        { Point2D( 0, 1), Point2D(-1, 1), Point2D(-1, 0), Point2D(1, 1),
          Point2D(-1, -1) }, {4, 8, 2, 6, 7}, 3);
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
      Bedt_{nrows, ncols},
      O_{nrows, ncols},
      adjMap_{nrows, ncols},
      adj4_{nrows, ncols, 1.0},
      Q_{SQUARE(static_cast<int>(MIN(ncols, nrows) / 2.0 + 1)), nrows * ncols},
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

    void EdtDIFT::setUpAdjMap()
    {
      int lastCol = adjMap_.getNumCols() - 1;
      for (int i = 1; i < adjMap_.getNumRows(); i++) {
        adjMap_[i * adjMap_.getNumCols()] = 1;            // fill first column
        adjMap_[(i * adjMap_.getNumCols()) + lastCol] = 2;  // fill last column
      }

      int lastRow = adjMap_.getNumCols() * (adjMap_.getNumRows()-1);
      for (int i = 1; i < adjMap_.getNumCols()-1; i++) {
        adjMap_[i] = 4;             // fill first row
        adjMap_[lastRow + i] = 3;   // fill last row
      }

      adjMap_[0] = 6;
      adjMap_[adjMap_.getNumCols()-1] = 8;
      adjMap_[adjMap_.getNumCols() * (adjMap_.getNumRows()-1)] = 5;
      adjMap_[adjMap_.getNumCols() * (adjMap_.getNumRows()-1) + adjMap_.getNumCols() - 1] = 7;
    }

    void EdtDIFT::treeRemoval(const std::vector<int> &toRemove)
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

    void EdtDIFT::run()
    {
      while (!Q_.isEmpty()) {
        int pidx = Q_.popMinFIFO();
        Point2D p = domain_.point(pidx);
        O_[pidx] = 0;

        int ridx = root_[pidx];

        Point2D r = domain_.point(ridx);

        Bedt_[ridx] = MAX(Bedt_[ridx], Q_.cost(pidx));

        const AdaptiveAdj &AA = AAB_[adjMap_[pidx]];
        for (const auto &[q, ai] : AA.neighborsPropogation(p)) {
          int qidx = domain_.index(q);
          int dx = q.x() - r.x();
          int dy = q.y() - r.y();
          int tmp = SQUARE(dx) + SQUARE(dy);

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
    
    void EdtDIFT::seed(int pidx)
    {
      root_[pidx] = pidx;
      Q_.setCost(pidx, 0);
      Q_.insert(pidx);
    }

    void EdtDIFT::open(int pidx)
    {
      O_[pidx] = 1;
      Q_.setCost(pidx, PQueue::PINF);
    }

    void EdtDIFT::insertNeighborsPQueue(int pidx)
    {
      for (int qidx : adj4_.getNeighborPixels(pidx)) {
        if (bin_[qidx] > 0 && Q_.cost(qidx) != PQueue::PINF && Q_.state(qidx) != PQueue::State::QUEUED) {
          Q_.setState(qidx, PQueue::State::NOT_PROCESSED);
          Q_.insert(qidx);
        }
      }
    }

    int EdtDIFT::maxBedt(const std::vector<int> &Ncontour) const
    {
      int max = 0;
      for (int pidx : Ncontour) {
        int d = Bedt_[pidx];
        if (d > max)
          max = d;
      }

      return max;
    }

    void EdtDIFT::saveDistanceTransform(const std::string &filename) const
    {
      const std::vector<int> &cost = Q_.cost();

      ImageUInt8Ptr d = ImageUInt8::create(bin_.getNumRows(), bin_.getNumCols());

      if (!filename.empty()) {
        int maxCost = 0;
        for (int pidx = 0; pidx < bin_.getSize(); ++pidx) {
          if (maxCost < cost[pidx])
            maxCost = cost[pidx];
        }

        for (int pidx = 0; pidx < bin_.getSize(); ++pidx) {
          (*d)[pidx] = static_cast<uint8_t>((static_cast<float>(cost[pidx]) / static_cast<float>(maxCost) * 255));
        }

        stbi_write_png(filename.c_str(), d->getNumCols(), d->getNumRows(), 1, d->rawData(), 0);     
      }
      else {
        ImageUInt8Ptr d = ImageUInt8::create(bin_.getNumRows(), bin_.getNumCols());
        for (int pidx = 0; pidx < bin_.getSize(); pidx++) {
          (*d)[pidx] = cost[pidx];
        }
        
        printImage(d);
        std::cout << "\n ---- Queue Cost --- \n";
      }
    }

    void EdtDIFT::saveUnderlyingBinaryImage(const std::string &filename) const
    {
      if (!filename.empty()) {
        ImageUInt8Ptr b = ImageUInt8::create(bin_.getNumRows(), bin_.getNumCols());
        for (int pidx = 0; pidx < b->getSize(); pidx++)
          (*b)[pidx] = bin_[pidx] == 0 ? 255 : 0; 

        stbi_write_png(filename.c_str(), b->getNumCols(), b->getNumRows(), 1, b->rawData(), 0);
      }
      else {
        ImageUInt8Ptr b = bin_.clone();
        printImage(b);
      }
    }

    void EdtDIFT::displayRootMapForSmallImages() const
    {
      for (int row = 0; row < root_.getNumRows(); row++) {
        std::cout << std::setw(4) << row;  // row index
        for (int col = 0; col < root_.getNumCols(); col++) {
          std::cout << std::setw(4) << root_[ImageUtils::to1D(row, col, root_.getNumCols())];
        }
        std::cout << "\n";
      }
    }
  }
}