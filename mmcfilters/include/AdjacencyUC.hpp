#ifndef ADJACENCY_UC_HPP
#define ADJACENCY_UC_HPP

#include <vector>
#include <cstdint>
#include <iterator>
#include "../include/ImageUtils.hpp"

enum class DiagonalConnection : uint8_t {
  None = 0,
  SW = 1 << 0,
  NE = 1 << 1,
  SE = 1 << 2,
  NW = 1 << 3
};

  inline DiagonalConnection operator|(DiagonalConnection a, DiagonalConnection b) {
    return static_cast<DiagonalConnection>(static_cast<uint8_t>(a) | static_cast<uint8_t>(b));
  }
  
  inline DiagonalConnection& operator|=(DiagonalConnection &a, DiagonalConnection b) {
    a = a | b;
    return a;
  }
  
  inline bool operator&(DiagonalConnection a, DiagonalConnection b) {
    return static_cast<uint8_t>(a) & static_cast<uint8_t>(b);
  }

class AdjacencyUC {

  private:
    int numRows, numCols;
    std::vector<DiagonalConnection> dconn_;

    // offsets: N, W, S, E, SW, NE, SE, NW
    const std::vector<int> dr_ = {-1, 0, 1, 0,  1, -1, 1, -1};
    const std::vector<int> dc_ = {0, -1, 0, 1, -1,  1, 1, -1};

    const std::vector<DiagonalConnection> required_diag_ = {
      DiagonalConnection::SW,
      DiagonalConnection::NE,
      DiagonalConnection::SE,
      DiagonalConnection::NW
    };

  public:
    AdjacencyUC(int numRows, int numCols) : numRows(numRows), numCols(numCols), dconn_(numRows * numCols, DiagonalConnection::None) {}
    void setDiagonalConnection(int row, int col, DiagonalConnection conn) {
      dconn_[ImageUtils::to1D(row, col, this->numCols)] |= conn;
    }

    void setDiagonalConnection(int idx, DiagonalConnection conn) {
      dconn_[idx] |= conn;
    }

    class NeighborIterator {
    private:
      AdjacencyUC &adj_;
      int row_, col_, id_;

      void advanceToValid() {
        while (id_ < adj_.dr_.size()) {
          int r = row_ + adj_.dr_[id_];
          int c = col_ + adj_.dc_[id_];
          if (r >= 0 && c >= 0 && r < adj_.numRows && c < adj_.numCols) {
            if (id_ < 4 || (adj_.dconn_[ImageUtils::to1D(row_, col_, adj_.numCols)] & adj_.required_diag_[id_ - 4])) {
              return;
            }
          }
          ++id_;
        }
      }
    public:
      NeighborIterator(AdjacencyUC &adj, int row, int col, int id)
        : adj_(adj), row_(row), col_(col), id_(id)
      {
        advanceToValid();
      }

      int operator*() const {
        int dr = adj_.dr_[id_];
        int dc = adj_.dc_[id_];
        return ImageUtils::to1D(row_ + dr, col_ + dc, adj_.numCols);
      }

      NeighborIterator& operator++() {
        ++id_;
        advanceToValid();
        return *this;
      }

      bool operator==(const NeighborIterator &other) const {
        return id_ == other.id_;
      }

      bool operator!=(const NeighborIterator &other) const {
        return !(*this == other);
      }
    };

    class NeighborRange {
      private:
        AdjacencyUC &adj_;
        int row_, col_;
      public:
        NeighborRange(AdjacencyUC &adj, int row, int col) : adj_(adj), row_(row), col_(col) {}
        NeighborIterator begin() { return NeighborIterator(adj_, row_, col_, 0); }
        NeighborIterator end() { return NeighborIterator(adj_, row_, col_, 8); }
      };
    
      NeighborRange getNeighboringPixels(int p) {
        auto [pRow, pCol] = ImageUtils::to2D(p, numCols);
        return NeighborRange(*this, pRow, pCol);
      }
      

      NeighborRange getNeighboringPixels(int row, int col) {
        return NeighborRange(*this, row, col);
      }

};

#endif // ADJACENCY_UC_HPP