#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <iterator>
#include <limits>
#include <numbers>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcfilters {

/**
 * @brief Integer displacement in a row-major regular 2D grid.
 */
struct GridOffset2D {
    /// Signed row displacement.
    int row = 0;
    /// Signed column displacement.
    int col = 0;

    /**
     * @brief Compares both displacement coordinates.
     *
     * @return True when the documented condition holds; otherwise false.
     */
    bool operator==(const GridOffset2D&) const noexcept = default;
};

/**
 * @brief Construction family of an immutable regular-grid stencil.
 */
enum class RegularGridAdjacencyShape { EuclideanDisk, StructuringElement };

/**
 * @brief Immutable regular-grid 2D adjacency with allocation-free traversal.
 *
 * `RegularGridAdjacency2D` stores a reusable stencil of offsets for a regular
 * 2D proper-part grid. The immutable stencil can be traversed either fully or in
 * forward-only mode, which exposes only one directed half of the neighbourhood
 * and is therefore convenient when unique undirected edges are needed.
 *
 * Coordinates use `(row, col)` order. Linear grid indices use
 * `row * numCols + col`. Each traversal range owns its cursor state,
 * so ranges over the same relation are reentrant and can be nested safely.
 */
class RegularGridAdjacency2D {
  private:
    struct StructuringElementTag {};

    /** @brief Stores the num cols. */
    int numCols;
    /** @brief Stores the num rows. */
    int numRows;
    /** @brief Stores the radius. */
    double radius;
    /** @brief Stores the radius2. */
    double radius2;
    /** @brief Stores the n. */
    int n;
    /** @brief Stores the shape. */
    RegularGridAdjacencyShape shape = RegularGridAdjacencyShape::EuclideanDisk;

    /** @brief Stores the offset row. */
    std::vector<int> offsetRow;
    /** @brief Stores the offset col. */
    std::vector<int> offsetCol;
    /** @brief Stores the forward offset indices. */
    std::vector<int> forwardOffsetIndices;

    /**
     * @brief Validates domain dimensions.
     *
     * @param rows Number of rows in the domain.
     * @param cols Number of columns in the domain.
     */
    static void requireDomainDimensions(int rows, int cols) {
        if (rows < 0 || cols < 0) {
            throw std::invalid_argument("RegularGridAdjacency2D grid dimensions must be non-negative.");
        }
    }

    /**
     * @brief Tests whether offset holds.
     *
     * @param rowOffset Row displacement of the neighbor.
     * @param colOffset Column displacement of the neighbor.
     * @return True when offset; otherwise false.
     */
    [[nodiscard]] bool containsOffset(int rowOffset, int colOffset) const noexcept {
        for (int index = 0; index < n; ++index) {
            if (offsetRow[static_cast<std::size_t>(index)] == rowOffset && offsetCol[static_cast<std::size_t>(index)] == colOffset) {
                return true;
            }
        }
        return false;
    }

    /**
     * @brief Checks whether the relation uses the expected ordered offsets.
     *
     * @param expected Expected ordered adjacency offsets.
     * @return True when the documented condition holds; otherwise false.
     */
    [[nodiscard]] bool matchesOffsets(std::initializer_list<GridOffset2D> expected) const noexcept {
        if (static_cast<std::size_t>(n) != expected.size()) {
            return false;
        }
        for (const GridOffset2D offset : expected) {
            if (!containsOffset(offset.row, offset.col)) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Builds forward offset indices.
     */
    void buildForwardOffsetIndices() {
        forwardOffsetIndices.clear();
        forwardOffsetIndices.reserve(static_cast<std::size_t>(n / 2));
        for (int index = 1; index < n; ++index) {
            const int dx = offsetCol[static_cast<std::size_t>(index)];
            const int dy = offsetRow[static_cast<std::size_t>(index)];
            if (dy > 0 || (dy == 0 && dx > 0)) {
                forwardOffsetIndices.push_back(index);
            }
        }
    }

    /**
     * @brief Constructs `RegularGridAdjacency2D` from the supplied inputs.
     *
     * @param rows Number of rows in the domain.
     * @param cols Number of columns in the domain.
     * @param offsets Adjacency offsets used to construct the regular grid relation.
     */
    RegularGridAdjacency2D(int rows, int cols, std::vector<GridOffset2D> offsets, StructuringElementTag)
        : numCols(cols), numRows(rows), radius(0.0), radius2(0.0), n(0), shape(RegularGridAdjacencyShape::StructuringElement) {
        requireDomainDimensions(rows, cols);
        if (offsets.empty()) {
            throw std::invalid_argument("A structuring-element adjacency requires at least the origin offset.");
        }

        for (const GridOffset2D offset : offsets) {
            if (offset.row == std::numeric_limits<int>::min() || offset.col == std::numeric_limits<int>::min()) {
                throw std::invalid_argument("Structuring-element offsets must be safely negatable.");
            }
        }

        std::sort(offsets.begin(), offsets.end(), [](const GridOffset2D& lhs, const GridOffset2D& rhs) {
            if (lhs.row != rhs.row) {
                return lhs.row < rhs.row;
            }
            return lhs.col < rhs.col;
        });
        const auto duplicate = std::adjacent_find(offsets.begin(), offsets.end());
        if (duplicate != offsets.end()) {
            throw std::invalid_argument("A structuring-element adjacency cannot contain duplicate offsets.");
        }

        const auto origin = std::find(offsets.begin(), offsets.end(), GridOffset2D{0, 0});
        if (origin == offsets.end()) {
            throw std::invalid_argument("A structuring-element adjacency must contain the origin offset.");
        }
        for (const GridOffset2D offset : offsets) {
            if (!std::binary_search(offsets.begin(), offsets.end(), GridOffset2D{-offset.row, -offset.col},
                                    [](const GridOffset2D& lhs, const GridOffset2D& rhs) {
                                        if (lhs.row != rhs.row) {
                                            return lhs.row < rhs.row;
                                        }
                                        return lhs.col < rhs.col;
                                    })) {
                throw std::invalid_argument("An adjacency-inducing structuring element must be centrally symmetric.");
            }
        }

        std::vector<GridOffset2D> ordered;
        ordered.reserve(offsets.size());
        ordered.push_back({0, 0});
        offsets.erase(origin);
        std::sort(offsets.begin(), offsets.end(), [](const GridOffset2D& lhs, const GridOffset2D& rhs) {
            auto angle = [](GridOffset2D offset) {
                double value = std::atan2(-static_cast<double>(offset.row), -static_cast<double>(offset.col));
                if (value < 0.0) {
                    value += 2.0 * std::numbers::pi;
                }
                return value;
            };
            const double lhsAngle = angle(lhs);
            const double rhsAngle = angle(rhs);
            if (lhsAngle != rhsAngle) {
                return lhsAngle < rhsAngle;
            }
            const std::int64_t lhsRadius = static_cast<std::int64_t>(lhs.row) * lhs.row + static_cast<std::int64_t>(lhs.col) * lhs.col;
            const std::int64_t rhsRadius = static_cast<std::int64_t>(rhs.row) * rhs.row + static_cast<std::int64_t>(rhs.col) * rhs.col;
            return lhsRadius < rhsRadius;
        });
        ordered.insert(ordered.end(), offsets.begin(), offsets.end());

        if (ordered.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            throw std::length_error("Structuring-element stencil exceeds the supported offset count.");
        }
        n = static_cast<int>(ordered.size());
        offsetRow.reserve(ordered.size());
        offsetCol.reserve(ordered.size());
        for (const GridOffset2D offset : ordered) {
            offsetRow.push_back(offset.row);
            offsetCol.push_back(offset.col);
            const double squaredDistance = static_cast<double>(offset.row) * offset.row + static_cast<double>(offset.col) * offset.col;
            radius2 = std::max(radius2, squaredDistance);
        }
        radius = std::sqrt(radius2);
        buildForwardOffsetIndices();
    }

    /**
     * @brief Validates coordinates.
     *
     * @param row Zero-based row coordinate.
     * @param col Zero-based column coordinate.
     */
    void requireCoordinates(int row, int col) const {
        if (row < 0 || row >= numRows || col < 0 || col >= numCols) {
            throw std::out_of_range("Index out of bounds.");
        }
    }

    /**
     * @brief Validates linear index.
     *
     * @param index Zero-based index used by the operation.
     */
    void requireLinearIndex(int index) const {
        const std::int64_t domainSize = static_cast<std::int64_t>(numRows) * static_cast<std::int64_t>(numCols);
        if (index < 0 || static_cast<std::int64_t>(index) >= domainSize) {
            throw std::out_of_range("Index out of bounds.");
        }
    }

  public:
    /**
     * @brief Builds an adjacency relation for a `numRows` by `numCols` grid.
     *
     * @param numRows Number of grid rows.
     * @param numCols Number of grid columns.
     * @param radius Radius of the neighbourhood stencil. `1.0` gives 4-connectivity
     * and `1.5` gives 8-connectivity on the integer grid.
     */
    RegularGridAdjacency2D(int numRows, int numCols, double radius) {
        requireDomainDimensions(numRows, numCols);
        const double maxSafeRadius = (std::sqrt(static_cast<double>(std::numeric_limits<int>::max())) - 1.0) / 2.0;
        if (!std::isfinite(radius) || radius < 0.0 || radius > maxSafeRadius) {
            throw std::invalid_argument("RegularGridAdjacency2D radius must be finite, non-negative, and representable by the integer stencil.");
        }
        this->numRows = numRows;
        this->numCols = numCols;
        this->radius = radius;
        this->radius2 = radius * radius;

        int i, j, k, dx, dy, r0, r2, i0 = 0;
        this->n = 0;
        r0 = (int)radius;
        r2 = (int)radius2;
        for (dy = -r0; dy <= r0; dy++)
            for (dx = -r0; dx <= r0; dx++)
                if (((dx * dx) + (dy * dy)) <= r2)
                    this->n++;

        i = 0;
        this->offsetCol.resize(this->n);
        this->offsetRow.resize(this->n);

        for (dy = -r0; dy <= r0; dy++) {
            for (dx = -r0; dx <= r0; dx++) {
                if (((dx * dx) + (dy * dy)) <= r2) {
                    this->offsetCol[i] = dx;
                    this->offsetRow[i] = dy;
                    if ((dx == 0) && (dy == 0))
                        i0 = i;
                    i++;
                }
            }
        }

        float aux;
        std::vector<float> da(n);
        std::vector<float> dr(n);

        /* Set clockwise */
        for (i = 0; i < n; i++) {
            dx = this->offsetCol[i];
            dy = this->offsetRow[i];
            dr[i] = std::sqrt((dx * dx) + (dy * dy));
            if (i != i0) {
                da[i] = (std::atan2(-dy, -dx) * 180.0 / std::numbers::pi);
                if (da[i] < 0.0)
                    da[i] += 360.0;
            }
        }
        da[i0] = 0.0;
        dr[i0] = 0.0;

        /* Place the central grid index first. */
        aux = da[i0];
        da[i0] = da[0];
        da[0] = aux;

        aux = dr[i0];
        dr[i0] = dr[0];
        dr[0] = aux;

        int auxX, auxY;
        auxX = this->offsetCol[i0];
        auxY = this->offsetRow[i0];
        this->offsetCol[i0] = this->offsetCol[0];
        this->offsetRow[i0] = this->offsetRow[0];

        this->offsetCol[0] = auxX;
        this->offsetRow[0] = auxY;

        /* sort by angle */
        for (i = 1; i < n - 1; i++) {
            k = i;
            for (j = i + 1; j < n; j++)
                if (da[j] < da[k]) {
                    k = j;
                }
            aux = da[i];
            da[i] = da[k];
            da[k] = aux;
            aux = dr[i];
            dr[i] = dr[k];
            dr[k] = aux;

            auxX = this->offsetCol[i];
            auxY = this->offsetRow[i];
            this->offsetCol[i] = this->offsetCol[k];
            this->offsetRow[i] = this->offsetRow[k];

            this->offsetCol[k] = auxX;
            this->offsetRow[k] = auxY;
        }

        /* sort by radius for each angle */
        for (i = 1; i < n - 1; i++) {
            k = i;
            for (j = i + 1; j < n; j++)
                if ((dr[j] < dr[k]) && (da[j] == da[k])) {
                    k = j;
                }
            aux = dr[i];
            dr[i] = dr[k];
            dr[k] = aux;

            auxX = this->offsetCol[i];
            auxY = this->offsetRow[i];
            this->offsetCol[i] = this->offsetCol[k];
            this->offsetRow[i] = this->offsetRow[k];

            this->offsetCol[k] = auxX;
            this->offsetRow[k] = auxY;
        }

        // Forward-only stencil keeps one orientation of every undirected edge
        // while preserving the complete clockwise stencil order.
        buildForwardOffsetIndices();
    }

    /**
     * @brief Builds adjacency induced by a symmetric structuring element.
     *
     * The origin must appear exactly once and every `(dr, dc)` offset must have
     * its opposite `(-dr, -dc)`. This undirectedness contract is required by
     * algorithms that traverse each grid edge once through
     * `getForwardNeighborIndices()`.
     *
     * @param numRows Number of rows in the domain.
     * @param numCols Number of columns in the domain.
     * @param offsets Symmetric structuring-element offsets.
     * @return The resulting adjacency induced by a symmetric structuring element.
     */
    [[nodiscard]] static RegularGridAdjacency2D fromStructuringElement(int numRows, int numCols, std::span<const GridOffset2D> offsets) {
        return RegularGridAdjacency2D(numRows, numCols, std::vector<GridOffset2D>(offsets.begin(), offsets.end()), StructuringElementTag{});
    }

    /**
     * @brief Builds a centered rectangular structuring-element adjacency.
     *
     * @param numRows Number of rows in the domain.
     * @param numCols Number of columns in the domain.
     * @param rowRadius Vertical radius of the rectangular stencil.
     * @param colRadius Horizontal radius of the rectangular stencil.
     * @return The resulting centered rectangular structuring-element adjacency.
     */
    [[nodiscard]] static RegularGridAdjacency2D rectangular(int numRows, int numCols, int rowRadius, int colRadius) {
        if (rowRadius < 0 || colRadius < 0) {
            throw std::invalid_argument("Rectangular adjacency radii must be non-negative.");
        }
        const std::int64_t height = 2 * static_cast<std::int64_t>(rowRadius) + 1;
        const std::int64_t width = 2 * static_cast<std::int64_t>(colRadius) + 1;
        const std::int64_t count = height * width;
        if (count > static_cast<std::int64_t>(std::numeric_limits<int>::max())) {
            throw std::length_error("Rectangular adjacency exceeds the supported offset count.");
        }

        std::vector<GridOffset2D> offsets;
        offsets.reserve(static_cast<std::size_t>(count));
        for (int rowOffset = -rowRadius; rowOffset <= rowRadius; ++rowOffset) {
            for (int colOffset = -colRadius; colOffset <= colRadius; ++colOffset) {
                offsets.push_back({rowOffset, colOffset});
            }
        }
        return fromStructuringElement(numRows, numCols, offsets);
    }

    /**
     * @brief Builds a centered digital line from `(-dr,-dc)` to `(dr,dc)`.
     *
     * Sampling follows the longest axis and rounds symmetrically, producing a
     * centrally symmetric digital segment suitable as an undirected adjacency.
     *
     * @param numRows Number of rows in the domain.
     * @param numCols Number of columns in the domain.
     * @param rowExtent Vertical extent of the digital line.
     * @param colExtent Horizontal extent of the digital line.
     * @return The resulting centered digital line from (-dr,-dc) to (dr,dc).
     */
    [[nodiscard]] static RegularGridAdjacency2D line(int numRows, int numCols, int rowExtent, int colExtent) {
        if (rowExtent == std::numeric_limits<int>::min() || colExtent == std::numeric_limits<int>::min()) {
            throw std::invalid_argument("Line adjacency extents must be safely negatable.");
        }
        const std::int64_t absoluteRow = std::abs(static_cast<std::int64_t>(rowExtent));
        const std::int64_t absoluteCol = std::abs(static_cast<std::int64_t>(colExtent));
        const std::int64_t steps = std::max(absoluteRow, absoluteCol);
        if (2 * steps + 1 > static_cast<std::int64_t>(std::numeric_limits<int>::max())) {
            throw std::length_error("Line adjacency exceeds the supported offset count.");
        }
        if (steps == 0) {
            const GridOffset2D origin{0, 0};
            return fromStructuringElement(numRows, numCols, std::span<const GridOffset2D>(&origin, 1));
        }

        auto roundedRatio = [](std::int64_t numerator, std::int64_t denominator) {
            if (numerator >= 0) {
                return (numerator + denominator / 2) / denominator;
            }
            return -((-numerator + denominator / 2) / denominator);
        };

        std::vector<GridOffset2D> offsets;
        offsets.reserve(static_cast<std::size_t>(2 * steps + 1));
        for (std::int64_t sample = -steps; sample <= steps; ++sample) {
            offsets.push_back({static_cast<int>(roundedRatio(sample * rowExtent, steps)), static_cast<int>(roundedRatio(sample * colExtent, steps))});
        }
        return fromStructuringElement(numRows, numCols, offsets);
    }

    /**
     * @brief Builds a centered horizontal-line adjacency.
     *
     * @param numRows Number of rows in the domain.
     * @param numCols Number of columns in the domain.
     * @param halfLength Half-length of the centred line.
     * @return The resulting centered horizontal-line adjacency.
     */
    [[nodiscard]] static RegularGridAdjacency2D horizontalLine(int numRows, int numCols, int halfLength) {
        if (halfLength < 0) {
            throw std::invalid_argument("Horizontal-line half-length must be non-negative.");
        }
        return line(numRows, numCols, 0, halfLength);
    }

    /**
     * @brief Builds a centered vertical-line adjacency.
     *
     * @param numRows Number of rows in the domain.
     * @param numCols Number of columns in the domain.
     * @param halfLength Half-length of the centred line.
     * @return The resulting centered vertical-line adjacency.
     */
    [[nodiscard]] static RegularGridAdjacency2D verticalLine(int numRows, int numCols, int halfLength) {
        if (halfLength < 0) {
            throw std::invalid_argument("Vertical-line half-length must be non-negative.");
        }
        return line(numRows, numCols, halfLength, 0);
    }

    /**
     * @brief Returns the number of offsets in the current stencil.
     *
     * The count includes the central origin offset at stencil position `0`.
     *
     * @return The number of offsets in the current stencil.
     */
    int getSize() const noexcept { return this->n; }

    /**
     * @brief Returns the number of rows in the attached grid domain.
     *
     * @return The number of rows in the attached grid domain.
     */
    int getNumRows() const noexcept { return numRows; }

    /**
     * @brief Returns the number of columns in the attached grid domain.
     *
     * @return The number of columns in the attached grid domain.
     */
    int getNumCols() const noexcept { return numCols; }

    /**
     * @brief Returns how the immutable stencil was constructed.
     *
     * @return How the immutable stencil was constructed.
     */
    RegularGridAdjacencyShape getShape() const noexcept { return shape; }

    /**
     * @brief Tests adjacency between two linear grid indices.
     *
     * A grid index is adjacent to itself because every supported stencil
     * includes the origin.
     *
     * @param p Point used by the operation.
     * @param q Second point used by the operation.
     * @return True when the documented condition holds; otherwise false.
     */
    inline bool isAdjacent(int p, int q) const noexcept {
        if (numCols <= 0) {
            return false;
        }
        int py = p / numCols, px = p % numCols;
        int qy = q / numCols, qx = q % numCols;

        return isAdjacent(px, py, qx, qy);
    }

    /**
     * @brief Tests adjacency between two grid coordinates.
     *
     * Coordinates are given as `(x, y)` pairs, where `x` is the column and `y`
     * is the row. The method applies the configured stencil but does not check
     * whether either coordinate lies inside the grid domain.
     *
     * @param px Row coordinate of the first grid point.
     * @param py Column coordinate of the first grid point.
     * @param qx Row coordinate of the second grid point.
     * @param qy Column coordinate of the second grid point.
     * @return True when the documented condition holds; otherwise false.
     */
    inline bool isAdjacent(int px, int py, int qx, int qy) const noexcept {
        const std::int64_t dx = static_cast<std::int64_t>(px) - qx;
        const std::int64_t dy = static_cast<std::int64_t>(py) - qy;
        if (shape == RegularGridAdjacencyShape::EuclideanDisk) {
            return static_cast<double>(dx) * dx + static_cast<double>(dy) * dy <= radius2;
        }
        if (dx < std::numeric_limits<int>::min() || dx > std::numeric_limits<int>::max() || dy < std::numeric_limits<int>::min() ||
            dy > std::numeric_limits<int>::max()) {
            return false;
        }
        return containsOffset(static_cast<int>(dy), static_cast<int>(dx));
    }

    /**
     * @brief Returns the configured or bounding Euclidean radius.
     *
     * Radius-built relations preserve the original input radius. For a
     * structuring element this is the maximum distance from its origin.
     *
     * @return The configured or bounding Euclidean radius.
     */
    double getRadius() const noexcept { return this->radius; }

    /**
     * @brief Returns true when the stencil represents canonical 4-connectivity.
     *
     * @return True when the stencil represents canonical 4-connectivity.
     */
    bool is4connectivity() const noexcept {
        if (shape == RegularGridAdjacencyShape::EuclideanDisk) {
            return radius == 1.0;
        }
        return matchesOffsets({{0, 0}, {-1, 0}, {0, -1}, {1, 0}, {0, 1}});
    }

    /**
     * @brief Returns true when the stencil represents canonical 8-connectivity.
     *
     * @return True when the stencil represents canonical 8-connectivity.
     */
    bool is8connectivity() const noexcept {
        if (shape == RegularGridAdjacencyShape::EuclideanDisk) {
            return radius == 1.5;
        }
        return matchesOffsets({{0, 0}, {-1, -1}, {-1, 0}, {-1, 1}, {0, -1}, {0, 1}, {1, -1}, {1, 0}, {1, 1}});
    }

    /**
     * @brief Tests whether grid-topology formulas may interpret this as 4/8 connectivity.
     *
     * @return True if grid-topology formulas may interpret this as 4/8 connectivity; otherwise false.
     */
    bool isCanonical4Or8Connectivity() const noexcept { return is4connectivity() || is8connectivity(); }

    /**
     * @brief Returns whether a linear grid index lies on the grid boundary.
     *
     * The index is interpreted in row-major order. No explicit bounds check is
     * performed before converting the index to `(row, col)`.
     *
     * @param index Zero-based index used by the operation.
     * @return Whether a linear grid index lies on the grid boundary.
     */
    bool isGridBoundary(int index) const {
        requireLinearIndex(index);
        return isGridBoundary(index / numCols, index % numCols);
    }

    /**
     * @brief Returns whether `(row, col)` lies on the grid boundary.
     *
     * The method assumes coordinates are inside the grid domain.
     *
     * @param row Zero-based row coordinate.
     * @param col Zero-based column coordinate.
     * @return Whether (row, col) lies on the grid boundary.
     */
    bool isGridBoundary(int row, int col) const noexcept { return row == 0 || col == 0 || row == this->numRows - 1 || col == this->numCols - 1; }

    /**
     * @brief Returns the row offset stored at stencil position `index`.
     *
     * The method does not perform bounds checking.
     *
     * @param index Zero-based index used by the operation.
     * @return The row offset stored at stencil position index.
     */
    int getOffsetRow(int index) const noexcept { return offsetRow[index]; }

    /**
     * @brief Returns the column offset stored at stencil position `index`.
     *
     * The method does not perform bounds checking.
     *
     * @param index Zero-based index used by the operation.
     * @return The column offset stored at stencil position index.
     */
    int getOffsetCol(int index) const noexcept { return offsetCol[index]; }

    /**
     * @brief Allocation-free iterator over one independent grid traversal.
     *
     * The iterator yields linear grid indices and respects both grid bounds
     * and the optional forward-only mask. All cursor state belongs to the
     * iterator, not to the adjacency relation.
     */
    template <bool ForwardOnly> class IteratorAdjacencyT {
      private:
        /** @brief Stores the relation. */
        const RegularGridAdjacency2D* relation_ = nullptr;
        /** @brief Stores the row. */
        int row_ = 0;
        /** @brief Stores the col. */
        int col_ = 0;
        /** @brief Stores the index. */
        int index_ = 0;

        /**
         * @brief Returns size.
         *
         * @return Size.
         */
        [[nodiscard]] int stencilSize() const noexcept {
            if constexpr (ForwardOnly) {
                return static_cast<int>(relation_->forwardOffsetIndices.size());
            }
            return relation_->n;
        }

        /**
         * @brief Returns the current adjacency-stencil index.
         *
         * @return Index of the active adjacency offset.
         */
        [[nodiscard]] int stencilIndex() const noexcept {
            if constexpr (ForwardOnly) {
                return relation_->forwardOffsetIndices[static_cast<std::size_t>(index_)];
            }
            return index_;
        }

        /**
         * @brief Advances to valid.
         */
        void seekValid() noexcept {
            const int size = stencilSize();
            while (index_ < size) {
                const int offsetIndex = stencilIndex();
                const std::int64_t neighborRow = static_cast<std::int64_t>(row_) + relation_->offsetRow[static_cast<std::size_t>(offsetIndex)];
                const std::int64_t neighborCol = static_cast<std::int64_t>(col_) + relation_->offsetCol[static_cast<std::size_t>(offsetIndex)];
                if (neighborRow >= 0 && neighborRow < relation_->numRows && neighborCol >= 0 && neighborCol < relation_->numCols) {
                    return;
                }
                ++index_;
            }
        }

      public:
        /// Forward-iterator category for independent multi-pass traversal.
        using iterator_category = std::forward_iterator_tag;
        /// C++20 iterator concept for the traversal cursor.
        using iterator_concept = std::forward_iterator_tag;

        /// Linear grid index yielded by the iterator.
        using value_type = int;
        /// Signed distance type used by iterator algorithms.
        using difference_type = std::ptrdiff_t;
        /// Value-returning reference type.
        using reference = int;
        /// No pointer type is exposed because dereference returns a value.
        using pointer = void;

        /**
         * @brief Constructs an empty adjacency iterator.
         */
        IteratorAdjacencyT() noexcept = default;

        /**
         * @brief Creates an iterator over one traversal context.
         *
         * @param relation Adjacency relation traversed by the iterator.
         * @param row Zero-based row coordinate.
         * @param col Zero-based column coordinate.
         * @param index Zero-based index used by the operation.
         */
        IteratorAdjacencyT(const RegularGridAdjacency2D* relation, int row, int col, int index) noexcept
            : relation_(relation), row_(row), col_(col), index_(index) {
            seekValid();
        }

        /**
         * @brief Advances this iterator without modifying the relation or other ranges.
         *
         * @return Reference to the resulting object.
         */
        IteratorAdjacencyT& operator++() noexcept {
            ++index_;
            seekValid();
            return *this;
        }

        /**
         * @brief Advances the iterator and returns its previous position.
         *
         * @return Iterator position before the advancement.
         */
        IteratorAdjacencyT operator++(int) noexcept {
            IteratorAdjacencyT previous = *this;
            ++(*this);
            return previous;
        }

        /**
         * @brief Returns true when both iterators identify the same traversal position.
         *
         * @param other Object to compare with or transfer from.
         * @return True when both iterators identify the same traversal position.
         */
        bool operator==(const IteratorAdjacencyT& other) const noexcept {
            return relation_ == other.relation_ && row_ == other.row_ && col_ == other.col_ && index_ == other.index_;
        }

        /**
         * @brief Returns true when the traversal positions differ.
         *
         * @param other Object to compare with or transfer from.
         * @return True when the traversal positions differ.
         */
        bool operator!=(const IteratorAdjacencyT& other) const noexcept { return !(*this == other); }

        /**
         * @brief Returns the current neighbour as a linear grid index.
         *
         * @return The current neighbour as a linear grid index.
         */
        int operator*() const noexcept {
            const int offsetIndex = stencilIndex();
            const std::int64_t neighborRow = static_cast<std::int64_t>(row_) + relation_->offsetRow[static_cast<std::size_t>(offsetIndex)];
            const std::int64_t neighborCol = static_cast<std::int64_t>(col_) + relation_->offsetCol[static_cast<std::size_t>(offsetIndex)];
            return static_cast<int>(neighborRow * relation_->numCols + neighborCol);
        }
    };

    /// Iterator over the complete stencil, including its configured origin.
    using IteratorAdjacency = IteratorAdjacencyT<false>;
    /// Iterator over one directed half of the stencil.
    using ForwardIteratorAdjacency = IteratorAdjacencyT<true>;

    /**
     * @brief Small value range carrying one immutable traversal context.
     */
    template <bool ForwardOnly, int FirstOffset> class GridIndexRangeT {
      private:
        /** @brief Stores the relation. */
        const RegularGridAdjacency2D* relation_;
        /** @brief Stores the row. */
        int row_;
        /** @brief Stores the col. */
        int col_;

      public:
        /**
         * @brief Creates a range for one validated grid coordinate.
         *
         * @param relation Adjacency relation traversed by the iterator.
         * @param row Zero-based row coordinate.
         * @param col Zero-based column coordinate.
         */
        GridIndexRangeT(const RegularGridAdjacency2D& relation, int row, int col) noexcept : relation_(&relation), row_(row), col_(col) {}

        /// Iterator type returned by this range.
        using iterator = IteratorAdjacencyT<ForwardOnly>;

        /**
         * @brief Returns an iterator positioned at the first valid neighbour.
         *
         * @return An iterator positioned at the first valid neighbour.
         */
        [[nodiscard]] iterator begin() const noexcept { return iterator(relation_, row_, col_, FirstOffset); }

        /**
         * @brief Returns the traversal sentinel iterator.
         *
         * @return The traversal sentinel iterator.
         */
        [[nodiscard]] iterator end() const noexcept {
            const int endIndex = [&] {
                if constexpr (ForwardOnly) {
                    return static_cast<int>(relation_->forwardOffsetIndices.size());
                }
                return relation_->n;
            }();
            return iterator(relation_, row_, col_, endIndex);
        }
    };

    /// Range over adjacent indices including the origin.
    using AdjacentIndexRange = GridIndexRangeT<false, 0>;
    /// Range over neighbouring indices excluding the origin.
    using NeighborIndexRange = GridIndexRangeT<false, 1>;
    /// Range over one directed half of the neighbouring indices.
    using ForwardNeighborIndexRange = GridIndexRangeT<true, 0>;

    /**
     * @brief Returns adjacent grid indices including the origin.
     *
     * @param row Zero-based row coordinate.
     * @param col Zero-based column coordinate.
     * @return Adjacent grid indices including the origin.
     */
    [[nodiscard]] AdjacentIndexRange getAdjacentIndices(int row, int col) const {
        requireCoordinates(row, col);
        return AdjacentIndexRange(*this, row, col);
    }

    /**
     * @brief Returns adjacent grid indices including the origin.
     *
     * This overload accepts a validated row-major linear index.
     *
     * @param gridIndex Index represented by `gridIndex`.
     * @return Adjacent grid indices including the origin.
     */
    [[nodiscard]] AdjacentIndexRange getAdjacentIndices(int gridIndex) const {
        requireLinearIndex(gridIndex);
        return getAdjacentIndices(gridIndex / numCols, gridIndex % numCols);
    }

    /**
     * @brief Returns valid neighbouring grid indices excluding the origin.
     *
     * @param row Zero-based row coordinate.
     * @param col Zero-based column coordinate.
     * @return Valid neighbouring grid indices excluding the origin.
     */
    [[nodiscard]] NeighborIndexRange getNeighborIndices(int row, int col) const {
        requireCoordinates(row, col);
        return NeighborIndexRange(*this, row, col);
    }

    /**
     * @brief Returns neighbouring grid indices excluding the origin.
     *
     * This overload accepts a validated row-major linear index.
     *
     * @param gridIndex Index represented by `gridIndex`.
     * @return Neighbouring grid indices excluding the origin.
     */
    [[nodiscard]] NeighborIndexRange getNeighborIndices(int gridIndex) const {
        requireLinearIndex(gridIndex);
        return getNeighborIndices(gridIndex / numCols, gridIndex % numCols);
    }

    /**
     * @brief Returns the directed positive half of the neighbourhood.
     *
     * @param row Zero-based row coordinate.
     * @param col Zero-based column coordinate.
     * @return The directed positive half of the neighbourhood.
     */
    [[nodiscard]] ForwardNeighborIndexRange getForwardNeighborIndices(int row, int col) const {
        requireCoordinates(row, col);
        return ForwardNeighborIndexRange(*this, row, col);
    }

    /**
     * @brief Returns the directed positive half of the neighbourhood.
     *
     * This overload accepts a validated row-major linear index.
     *
     * @param gridIndex Index represented by `gridIndex`.
     * @return The directed positive half of the neighbourhood.
     */
    [[nodiscard]] ForwardNeighborIndexRange getForwardNeighborIndices(int gridIndex) const {
        requireLinearIndex(gridIndex);
        return getForwardNeighborIndices(gridIndex / numCols, gridIndex % numCols);
    }
};

} // namespace mmcfilters
