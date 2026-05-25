#pragma once

#include "Image.hpp"

#include <cmath>
#include <cstdint>
#include <iterator>
#include <numbers>
#include <stdexcept>
#include <vector>


namespace mmcfilters {

/**
 * @brief Two-dimensional adjacency relation with configurable radius and efficient iteration.
 *
 * `AdjacencyRelation` stores a reusable stencil of offsets for a regular image
 * grid. The stencil can be traversed either fully or in forward-only mode,
 * which exposes only one directed half of the neighbourhood and is therefore
 * convenient when unique undirected edges are needed.
 *
 * Coordinates use `(row, col)` order for row-major image domains. Linear pixel
 * ids use `row * numCols + col`. Iteration methods mutate an internal cursor,
 * so one `AdjacencyRelation` object should not be traversed by multiple active
 * ranges at the same time.
 */
class AdjacencyRelation {
private:
    int id;
    
    int row;
    int col;    
    int numCols;
    int numRows;
    double radius;
    double radius2;
    int n;
    bool forwardOnly = false;

    std::vector<int> offsetRow;
    std::vector<int> offsetCol;
    std::vector<uint8_t> forwardMask;


public:
    /**
     * @brief Builds an adjacency relation for a `numRows` by `numCols` image.
     *
     * @param numRows Number of image rows.
     * @param numCols Number of image columns.
     * @param radius Radius of the neighbourhood stencil. `1.0` gives 4-connectivity
     * and `1.5` gives 8-connectivity on the integer grid.
     */    
    AdjacencyRelation(int numRows, int numCols, double radius){
        this->numRows = numRows;
        this->numCols = numCols;
        this->radius = radius;
        this->radius2 = radius * radius;

        int i, j, k, dx, dy, r0, r2, i0 = 0;
        this->n = 0;
        r0 = (int) radius;
        r2 = (int) radius2;
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
                    this->offsetCol[i] =dx;
                    this->offsetRow[i] =dy;
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

        /* place central pixel at first */
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

        // Forward-only mask: keep only offsets in the positive half-plane.
        forwardMask.resize(n, 0);
        for (int k = 1; k < n; ++k) {
            int dx = offsetCol[k], dy = offsetRow[k];
            forwardMask[k] = (dy > 0 || (dy == 0 && dx > 0)) ? 1 : 0;
        }
        forwardMask[0] = 0;
    }


    /**
     * @brief Advances to the next valid offset under the current bounds and mode.
     * @return Index of the next valid offset, or `getSize()` if none remain.
     */
    int nextValid() {
        id += 1;
        while (id < n) {

            // Apply the forward-only mask when required.
            if (forwardOnly && !forwardMask[id]) { id += 1; continue; }

            // Candidate neighbour coordinates.
            const int newRow = row + offsetRow[id];
            const int newCol = col + offsetCol[id];

            if (newRow >= 0 && newRow < numRows && newCol >= 0 && newCol < numCols) {
                return id;
            }
            id += 1;
        }
        return n;
    }

    /**
     * @brief Returns the number of offsets in the current stencil.
     *
     * The count includes the central origin offset at stencil position `0`.
     */
    int getSize() const {
        return this->n;
    }

    /**
     * @brief Returns the number of rows in the attached image domain.
     */
    int getNumRows() const noexcept {
        return numRows;
    }

    /**
     * @brief Returns the number of columns in the attached image domain.
     */
    int getNumCols() const noexcept {
        return numCols;
    }

    /**
     * @brief Prepares iteration over adjacent pixels, including the origin.
     *
     * @param row Zero-based row coordinate inside the image domain.
     * @param col Zero-based column coordinate inside the image domain.
     * @return This adjacency relation configured as a range.
     * @throws std::out_of_range If `(row, col)` is outside the image domain.
     */
    AdjacencyRelation& getAdjPixels(int row, int col){
        if (row < 0 || row >= this->numRows || col < 0 || col >= this->numCols) {
            throw std::out_of_range("Index out of bounds.");
        }
        this->row = row;
        this->col = col;
        this->id = -1;
        this->forwardOnly = false;

        return *this;
    }

    /**
     * @brief Prepares iteration over adjacent pixels from a linear index, including the origin.
     *
     * @param indexVector Row-major linear pixel index.
     * @return This adjacency relation configured as a range.
     * @throws std::out_of_range If `indexVector` maps outside the image domain.
     */
    AdjacencyRelation& getAdjPixels(int indexVector){
        return getAdjPixels(indexVector / this->numCols, indexVector % this->numCols);
    }
    

    /**
     * @brief Prepares iteration over valid neighbouring pixels, excluding the origin.
     *
     * @param row Zero-based row coordinate inside the image domain.
     * @param col Zero-based column coordinate inside the image domain.
     * @return This adjacency relation configured as a range.
     * @throws std::out_of_range If `(row, col)` is outside the image domain.
     */
    AdjacencyRelation& getNeighborPixels(int row, int col){
        if (row < 0 || row >= this->numRows || col < 0 || col >= this->numCols) {
            throw std::out_of_range("Index out of bounds.");
        }
        this->row = row;
        this->col = col;
        this->id = 0;
        this->forwardOnly = false;
        return *this;
    }

    /**
     * @brief Prepares iteration over valid neighbouring pixels from a linear index, excluding the origin.
     *
     * @param indexVector Row-major linear pixel index.
     * @return This adjacency relation configured as a range.
     * @throws std::out_of_range If `indexVector` maps outside the image domain.
     */
    AdjacencyRelation& getNeighborPixels(int indexVector){
        return getNeighborPixels(indexVector / this->numCols, indexVector % this->numCols);
    }

    /**
     * @brief Prepares forward-only neighbour iteration at `(row, col)`.
     *
     * @details Only one directed half of the neighbourhood is emitted:
     * neighbours whose offset has `dy > 0`, or `dy == 0 && dx > 0`.
     *
     * @param row Zero-based row coordinate inside the image domain.
     * @param col Zero-based column coordinate inside the image domain.
     * @return This adjacency relation configured as a range.
     * @throws std::out_of_range If `(row, col)` is outside the image domain.
     */
    AdjacencyRelation& getNeighborPixelsForward(int row, int col){
        if (row < 0 || row >= this->numRows || col < 0 || col >= this->numCols) {
            throw std::out_of_range("Index out of bounds.");
        }
        this->row = row;
        this->col = col;
        this->id = 0;
        this->forwardOnly = true;
        return *this;
    }
    
    /**
     * @brief Prepares forward-only neighbour iteration from a linear index.
     *
     * @param indexVector Row-major linear pixel index.
     * @return This adjacency relation configured as a range.
     * @throws std::out_of_range If `indexVector` maps outside the image domain.
     */
    AdjacencyRelation& getNeighborPixelsForward(int indexVector){
        return getNeighborPixelsForward(indexVector / this->numCols, indexVector % this->numCols);
    }
    

    /**
     * @brief Tests adjacency between two linear pixel indices.
     *
     * The central pixel is adjacent to itself because the radius test includes a
     * zero-length offset.
     */
    inline bool isAdjacent(int p, int q) const noexcept {
        int py = p / numCols, px = p % numCols;
        int qy = q / numCols, qx = q % numCols;

        return isAdjacent(px, py, qx, qy);
    }

    /**
     * @brief Tests adjacency between two pixel coordinates.
     *
     * Coordinates are given as `(x, y)` pairs, where `x` is the column and `y`
     * is the row. The method only applies the radius test; it does not check
     * whether either coordinate lies inside the configured image domain.
     */
    inline bool isAdjacent(int px, int py, int qx, int qy) const noexcept {
        int dx = px - qx;
        int dy = py - qy;
        return double(dx)*dx + double(dy)*dy <= radius2;
    }

    /**
     * @brief Returns the configured neighbourhood radius.
     */
    double getRadius() const {
        return this->radius;
    }

    /**
     * @brief Returns true when the stencil represents 4-connectivity.
     */
    bool is4connectivity() const { return this->radius == 1;}

    /**
     * @brief Returns true when the stencil represents 8-connectivity.
     */
    bool is8connectivity() const {return this->radius == 1.5;}

    /**
     * @brief Returns whether a linear pixel index lies on the image border.
     *
     * The index is interpreted in row-major order. No explicit bounds check is
     * performed before converting the index to `(row, col)`.
     */
    bool isBorderDomainImage(int index){
        auto[row, col] = ImageUtils::to2D(index, this->numCols);
        return isBorderDomainImage(row, col);
    }

    /**
     * @brief Returns whether `(row, col)` lies on the image border.
     *
     * The method assumes coordinates are inside the image domain.
     */
    bool isBorderDomainImage(int row, int col){
        return row == 0 || col == 0 || row == this->numRows - 1 || col == this->numCols - 1;
    }

    /**
     * @brief Returns the row offset stored at stencil position `index`.
     *
     * The method does not perform bounds checking.
     */
    int getOffsetRow(int index){
        return offsetRow[index];
    }

    /**
     * @brief Returns the column offset stored at stencil position `index`.
     *
     * The method does not perform bounds checking.
     */
    int getOffsetCol(int index){
        return offsetCol[index];
    }


    /**
     * @brief Lightweight iterator over the currently configured traversal.
     *
     * The iterator yields linear pixel indices and respects both image bounds
     * and the optional forward-only mask.
     *
     * Iterators hold a pointer to the parent `AdjacencyRelation`; incrementing
     * any iterator advances the parent's mutable traversal cursor.
     */
    class IteratorAdjacency { 
    private:
        int index;
        AdjacencyRelation* instance; 

    public:
        /// Input-iterator category for range-based traversal.
        using iterator_category = std::input_iterator_tag;

        /// Linear pixel index yielded by the iterator.
        using value_type = int;

        /**
         * @brief Creates an iterator over `obj` at stencil position `id`.
         */
        IteratorAdjacency(AdjacencyRelation* obj, int id) :  index(id), instance(obj) { }

        /**
         * @brief Returns the adjacency relation currently being traversed.
         */
        AdjacencyRelation* getInstance() { return instance; } 

        /**
         * @brief Advances to the next valid neighbour in the current traversal.
         */
        IteratorAdjacency& operator++() { 
            this->index = instance->nextValid();  
            return *this; 
        }

        /**
         * @brief Returns true when both iterators point to the same stencil position.
         */
        bool operator==(const IteratorAdjacency& other) const { 
            return index == other.index; 
        }

        /**
         * @brief Returns true when the iterators point to different stencil positions.
         */
        bool operator!=(const IteratorAdjacency& other) const { 
            return !(*this == other);
        }

        /**
         * @brief Returns the current neighbour as a linear pixel index.
         */
        int operator*() const { 
            return (instance->row + instance->offsetRow[index]) * instance->numCols + (instance->col + instance->offsetCol[index]); 
        }
    };
    /**
     * @brief Returns the beginning of the current traversal.
     *
     * Call one of the `get*Pixels` methods before requesting `begin()`.
     */
    IteratorAdjacency begin() { 
        return IteratorAdjacency(this, nextValid()); 
    }

    /**
     * @brief Returns the end sentinel of the current traversal.
     */
    IteratorAdjacency end() { 
        return IteratorAdjacency(this, this->n); 
    }
};

} // namespace mmcfilters
