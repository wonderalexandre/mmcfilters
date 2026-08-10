#pragma once

#include "RegularGridAdjacency2D.hpp"

namespace mmcfilters::detail {

/**
 * @brief Allocation-free grid traversal after the central pixel was established.
 *
 * The public adjacency API validates its central coordinate.  Internal kernels
 * whose loop bounds already prove `0 <= pixel < rows * cols` use this facade so
 * checked builds do not repeat that validation in the hot loop.
 */
class CommittedGridAccess {
  public:
    /**
     * @brief Constructs a radius adjacency from established parameters.
     * @param rows Established row count.
     * @param cols Established column count.
     * @param radius Established finite radius.
     * @return Constructed regular-grid adjacency.
     */
    [[nodiscard]] static RegularGridAdjacency2D radiusAdjacency(int rows, int cols, double radius) {
        return RegularGridAdjacency2D(rows, cols, radius, RegularGridAdjacency2D::EstablishedRadiusTag{});
    }

    /** @brief Iterates the full adjacent stencil of an established pixel. @param adjacency Established adjacency. @param pixel Established pixel id. @return Adjacent-index range. */
    [[nodiscard]] static RegularGridAdjacency2D::AdjacentIndexRange adjacent(const RegularGridAdjacency2D& adjacency, int pixel) noexcept {
        const int cols = adjacency.getNumCols();
        return RegularGridAdjacency2D::AdjacentIndexRange(adjacency, pixel / cols, pixel % cols);
    }

    /** @brief Iterates in-domain neighbors of an established pixel. @param adjacency Established adjacency. @param pixel Established pixel id. @return Neighbor-index range. */
    [[nodiscard]] static RegularGridAdjacency2D::NeighborIndexRange neighbors(const RegularGridAdjacency2D& adjacency, int pixel) noexcept {
        const int cols = adjacency.getNumCols();
        return RegularGridAdjacency2D::NeighborIndexRange(adjacency, pixel / cols, pixel % cols);
    }

    /** @brief Iterates forward neighbors of an established pixel. @param adjacency Established adjacency. @param pixel Established pixel id. @return Forward-neighbor range. */
    [[nodiscard]] static RegularGridAdjacency2D::ForwardNeighborIndexRange forwardNeighbors(const RegularGridAdjacency2D& adjacency, int pixel) noexcept {
        const int cols = adjacency.getNumCols();
        return RegularGridAdjacency2D::ForwardNeighborIndexRange(adjacency, pixel / cols, pixel % cols);
    }

    /** @brief Tests whether an established pixel lies on the grid boundary. @param adjacency Established adjacency. @param pixel Established pixel id. @return True for a boundary pixel. */
    [[nodiscard]] static bool isBoundary(const RegularGridAdjacency2D& adjacency, int pixel) noexcept {
        const int cols = adjacency.getNumCols();
        return adjacency.isGridBoundary(pixel / cols, pixel % cols);
    }
};

} // namespace mmcfilters::detail
