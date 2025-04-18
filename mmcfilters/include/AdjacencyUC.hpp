#ifndef ADJACENCY_UC_HPP
#define ADJACENCY_UC_HPP

#include <cstdint>
#include <vector>
#include "../include/ImageUtils.hpp"
#include "../include/Common.hpp"

enum class DiagonalConnection : uint8_t {
  None = 0,
  SW = 1 << 0,
  NE = 1 << 1,
  SE = 1 << 2,
  NW = 1 << 3
};

// Operadores auxiliares
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
  uint8_t* dconnFlags = nullptr;     // 4-connect.  +  diag. connect.
                                    //  N, W, S, E,   SW, NE, SE, NW
  const std::vector<int> offsetRows = {-1, 0, 1, 0,    1, -1,  1, -1}; 
  const std::vector<int> offsetCols = {0, -1, 0, 1,    -1,  1,  1, -1};
  bool enableDiagonalConnection;
  const std::vector<DiagonalConnection> requiredDiagonal = {
    DiagonalConnection::SW, DiagonalConnection::NE,
    DiagonalConnection::SE, DiagonalConnection::NW
  };

public:
  AdjacencyUC(int rows, int cols, bool enableDiagonalConnection) : numRows(rows), numCols(cols), enableDiagonalConnection(enableDiagonalConnection){
    if(enableDiagonalConnection)
      dconnFlags = new uint8_t[rows * cols]();
  }

  ~AdjacencyUC() {
    delete[] dconnFlags;
  }

  void setDiagonalConnection(int row, int col, DiagonalConnection conn) {
    dconnFlags[ImageUtils::to1D(row, col, numCols)] |= static_cast<uint8_t>(conn);
  }

  void setDiagonalConnection(int idx, DiagonalConnection conn) {
    dconnFlags[idx] |= static_cast<uint8_t>(conn);
  }

  bool hasConnection(int row, int col, DiagonalConnection conn) const {
    return dconnFlags[ImageUtils::to1D(row, col, numCols)] & static_cast<uint8_t>(conn);
  }

  uint8_t getConnections(int row, int col) const {
    return dconnFlags[ImageUtils::to1D(row, col, numCols)];
  }

  class NeighborIterator {
  private:
    AdjacencyUC &instance;
    int row, col, id;

    void advanceToValid() {
      while (id< instance.offsetRows.size()) {
        int r = row + instance.offsetRows[id];
        int c = col + instance.offsetCols[id];
        if (r >= 0 && c >= 0 && r < instance.numRows && c < instance.numCols) {
          if (id < 4 || (instance.enableDiagonalConnection && instance.dconnFlags[ImageUtils::to1D(row, col, instance.numCols)] & static_cast<uint8_t>(instance.requiredDiagonal[id - 4]))) {
            return;
          }
        }
        ++id;
      }
    }

  public:
    NeighborIterator(AdjacencyUC &adj, int row, int col, int id): instance(adj), row(row), col(col), id(id){
      advanceToValid();
    }

    int operator*() const {
      int dr = instance.offsetRows[id];
      int dc = instance.offsetCols[id];
      return ImageUtils::to1D(row + dr, col + dc, instance.numCols);
    }

    NeighborIterator& operator++() {
      ++id;
      advanceToValid();
      return *this;
    }

    bool operator==(const NeighborIterator &other) const {
      return id == other.id;
    }

    bool operator!=(const NeighborIterator &other) const {
      return !(*this == other);
    }
  };

  class NeighborRange {
  private:
    AdjacencyUC &instance;
    int row, col;
    
  public:
    NeighborRange(AdjacencyUC &instance, int row, int col)
      : instance(instance), row(row), col(col) {}

    NeighborIterator begin() { return NeighborIterator(instance, row, col, 0); }
    NeighborIterator end() { return NeighborIterator(instance, row, col, 8); }
  };

  NeighborRange getNeighboringPixels(int p) {
    auto [row, col] = ImageUtils::to2D(p, numCols);
    return NeighborRange(*this, row, col);
  }

  NeighborRange getNeighboringPixels(int row, int col) {
    return NeighborRange(*this, row, col);
  }


};

#endif // ADJACENCY_UC_HPP
