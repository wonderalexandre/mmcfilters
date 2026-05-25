#pragma once

#include "../utils/AdjacencyRelation.hpp"
#include "../utils/Image.hpp"
#include "../utils/Common.hpp"

#include <algorithm>
#include <cstdint>
#include <deque>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <vector>

namespace mmcfilters {

/**
 * @brief Selects the interpolation/connectivity convention used to build a tree of shapes.
 */
enum class ToSInterpolation {
    SelfDual,
    Min4cMax8c,
    Min8cMax4c
};

using ToSLevel = uint16_t;

inline constexpr int ToSInterpolationScale = 2;
inline constexpr int ToSInterpolationPadding = 1;
inline constexpr int ToSUInt8Depth = 8;
inline constexpr int ToSSelfDualDepth = 9;
inline constexpr int ToSDefaultInfinityRow = 0;
inline constexpr int ToSDefaultInfinityCol = 0;

/**
 * @brief Common interface for union-find builders of morphological trees.
 */
class IMorphologicalTreeBuilder {
public:
    virtual ~IMorphologicalTreeBuilder() = default;

    /**
     * @brief Builds the union-find parent representation for an 8-bit image.
     *
     * @return Tuple containing the parent image, the processed pixel order, and
     * the number of created tree nodes.
     */
    virtual std::tuple<std::vector<int>, std::vector<int>, int> createTreeByUnionFind(const ImageUInt8Ptr& imgPtr) const = 0;
};


/**
 * @brief Union-find builder for max-trees and min-trees.
 */
class BuilderComponentTree : public IMorphologicalTreeBuilder{
private:
    AdjacencyRelation* adj;
    bool isMaxtree;

public:
    /**
     * @brief Creates a component-tree builder using `adj` as image adjacency.
     *
     * `isMaxtree` selects max-tree ordering when true and min-tree ordering when
     * false. The adjacency object is borrowed and must outlive the builder.
     */
    explicit BuilderComponentTree(AdjacencyRelation* adj, bool isMaxtree) : adj(adj), isMaxtree(isMaxtree) { }
    ~BuilderComponentTree() { }

    /// @cond INTERNAL
    template <typename PixelType>
    std::vector<int> sort(ImagePtr<PixelType> imgPtr) const {
        const int n = imgPtr->getSize();
        if (n <= 0) {
            throw std::invalid_argument("BuilderComponentTree requires a non-empty image.");
        }
        std::vector<int> orderedPixels(n);
        PixelType* img = imgPtr->rawData();

        if constexpr (!std::is_same_v<PixelType, std::uint8_t>) {
            if (PRINT_LOG) std::cout << "Sorting image with comparison sort, size: " << n << std::endl;
            std::iota(orderedPixels.begin(), orderedPixels.end(), 0);
            if (isMaxtree) {
                std::stable_sort(orderedPixels.begin(), orderedPixels.end(),
                                 [&](int a, int b) { return img[a] < img[b]; });
            } else {
                std::stable_sort(orderedPixels.begin(), orderedPixels.end(),
                                 [&](int a, int b) { return img[a] > img[b]; });
            }
        } else {
            if (PRINT_LOG) std::cout << "Sorting uint8 image with counting sort, size: " << n << std::endl;
        
            // counting sort com faixa [0..maxvalue]; 
            int maxvalue =  static_cast<int>(img[0]);
            for (int i = 1; i < n; i++) if(maxvalue < img[i]) maxvalue = img[i];
            std::vector<uint32_t> counter(static_cast<size_t>(maxvalue) + 1, 0);
            if(isMaxtree){
                for (int i = 0; i < n; i++)
                    counter[img[i]]++;
                for (int i = 1; i <= maxvalue; i++)
                    counter[i] += counter[i - 1];
                for (int i = n - 1; i >= 0; --i)
                    orderedPixels[--counter[img[i]]] = i;	

            }else{
                for (int i = 0; i < n; i++)
                    counter[maxvalue - img[i]]++;
                for (int i = 1; i <= maxvalue; i++)
                    counter[i] += counter[i - 1];
                for (int i = n - 1; i >= 0; --i)
                    orderedPixels[--counter[maxvalue - img[i]]] = i;
            }
        }
        return orderedPixels;
    }
    /// @endcond

    /**
     * @brief Builds the parent image of a uint8 max-tree or min-tree.
     */
    std::tuple<std::vector<int>, std::vector<int>, int> createTreeByUnionFind(const ImageUInt8Ptr& imgPtr) const override {
        return createTreeByUnionFind<uint8_t>(imgPtr);
    }


    /**
     * @brief Builds the parent image of a typed max-tree or min-tree.
     */
    template <typename PixelType>
    std::tuple<std::vector<int>, std::vector<int>, int> createTreeByUnionFind(const ImagePtr<PixelType>& imgPtr) const {
        std::vector<int> orderedPixels = sort(imgPtr);
        auto img = imgPtr->rawData();

        int numPixels = imgPtr->getSize();
        std::vector<int> zPar(numPixels, InvalidNode);
        std::vector<int> parent(numPixels, InvalidNode);
        auto findRoot = [&](int p) {
            while (zPar[p] != p) { zPar[p] = zPar[zPar[p]]; p = zPar[p]; }
            return p;
        };

        for (int i = numPixels - 1; i >= 0; i--) {
            int p = orderedPixels[i];
            parent[p] = p;
            zPar[p] = p;
            for (int q : adj->getNeighborPixels(p)) {
                if (zPar[q] != InvalidNode) {
                    int r = findRoot(q);
                    if (p != r) { parent[r] = p; zPar[r] = p; }
                }
            }
        }

        int numNodes = 0;
        for (int i = 0; i < numPixels; i++) {
            int p = orderedPixels[i];
            int q = parent[p];
            if (img[parent[q]] == img[q]) parent[p] = parent[q];
            if (parent[p] == p || img[parent[p]] != img[p]) ++numNodes;
        }
        return std::make_tuple(std::move(parent), std::move(orderedPixels), std::move(numNodes));
    }
};



	/************************ Tree of Shapes support ************************/

    /// @cond INTERNAL
	/*
	 * Adaptive adjacency backend used by the tree-of-shapes construction.
 *
 * Diagonal links are activated on demand so that the interpolated grid can
 * emulate the required 4/8-connectivity behaviour during the union-find pass.
 * Each pixel may carry four diagonal flags: SW, NE, SE, and NW.
 */
enum class DiagonalConnection : uint8_t {
    None = 0,
    SW = 1 << 0,
    NE = 1 << 1,
    SE = 1 << 2,
    NW = 1 << 3
};

	// Helper operators for diagonal-connection flags.
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

/**
 * @brief Adaptive adjacency used during 4/8-connected tree-of-shapes construction.
 */
class AdjacencyUC {
private:
    int numRows, numCols;
    std::vector<uint8_t> dconnFlags;     // 4-connect.  +  diag. connect.
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
        dconnFlags.resize(rows * cols, 0);
    }

    ~AdjacencyUC() {
        
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

	/**
	 * @brief Iterator over valid neighbours including enabled diagonal links.
	 */
	class NeighborIterator {
    private:
        AdjacencyUC &instance;
        int row, col;
        std::size_t id;

        void advanceToValid() {
        while (id < instance.offsetRows.size()) {
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

	/**
	 * @brief Range helper that produces valid-neighbour iterators.
	 */
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

    NeighborRange getNeighborPixels(int p) {
        auto [row, col] = ImageUtils::to2D(p, numCols);
        return NeighborRange(*this, row, col);
    }

    NeighborRange getNeighborPixels(int row, int col) {
        return NeighborRange(*this, row, col);
    }

};


	/**
	 * @brief Discrete priority queue used during tree-of-shapes construction.
	 */
	class PriorityQueueToS {
private:
    std::vector<std::deque<int>> buckets;
    int currentPriority;
    int numElements;
    int maxPriorityLevels;
    

public:
    PriorityQueueToS(int depthOfImage=ToSUInt8Depth) : currentPriority(0), numElements(0), maxPriorityLevels(1 << depthOfImage){
        buckets.resize(maxPriorityLevels);
    }

    void initial(int element, int priority) {
        currentPriority = priority;
        buckets[priority].push_back(element);
        numElements++;
    }
    int getCurrentPriority()  {return currentPriority;}
    bool isEmpty()  {return numElements == 0;}

    void priorityPush(int element, int lower, int upper) {
        int priority;
        if (lower > currentPriority) {
            priority = lower;
        } else if (upper < currentPriority) {
            priority = upper;
        } else {
            priority = currentPriority;
        }
        numElements++;
        buckets[priority].push_back(element);
    }

    int priorityPop() {
        if (buckets[currentPriority].empty()) {
            int nextPriority = -1;
            for (int priority = currentPriority + 1; priority < maxPriorityLevels; ++priority) {
                if (!buckets[priority].empty()) {
                    nextPriority = priority;
                    break;
                }
            }

            if (nextPriority == -1) {
                for (int priority = currentPriority - 1; priority >= 0; --priority) {
                    if (!buckets[priority].empty()) {
                        nextPriority = priority;
                        break;
                    }
                }
            }

            if (nextPriority == -1) {
                throw std::runtime_error("PriorityQueueToS is empty.");
            }
            currentPriority = nextPriority;
        }

        int element = buckets[currentPriority].front(); 
        buckets[currentPriority].pop_front();           

        numElements--;  
        return element;
	    }
	};
    /// @endcond




	/**
	 * @brief Builds trees of shapes (ToS) using a union-find construction.
	 *
	 * The builder interpolates the input image according to `ToSInterpolation`,
	 * constructs a max-tree on the interpolated domain, and projects the resulting
	 * hierarchy back to the original image domain.
	 */
	class BuilderTreeOfShape: public IMorphologicalTreeBuilder {
private:
    ToSInterpolation interpolation;
    int infinitySeedRow;
    int infinitySeedCol;

    inline bool usesConnectivityMap() const noexcept {
        return interpolation == ToSInterpolation::Min4cMax8c
            || interpolation == ToSInterpolation::Min8cMax4c;
    }

    inline bool usesHighDiagonalAtSaddle() const noexcept {
        return interpolation == ToSInterpolation::Min4cMax8c;
    }

    inline int interpolatedNumRows(int numRows) const noexcept {
        return ToSInterpolationScale * numRows + ToSInterpolationPadding;
    }

    inline int interpolatedNumCols(int numCols) const noexcept {
        return ToSInterpolationScale * numCols + ToSInterpolationPadding;
    }

    inline int originalPointRow(int row) const noexcept {
        return ToSInterpolationScale * row + ToSInterpolationPadding;
    }

    inline int originalPointCol(int col) const noexcept {
        return ToSInterpolationScale * col + ToSInterpolationPadding;
    }

    inline ToSLevel scaledOriginalLevel(uint8_t value) const noexcept {
        return static_cast<ToSLevel>(ToSInterpolationScale * static_cast<int>(value));
    }

    inline int infinitySeedIndex(int interpNumRows, int interpNumCols) const {
        if (infinitySeedRow < 0 || infinitySeedCol < 0 || infinitySeedRow >= interpNumRows || infinitySeedCol >= interpNumCols) {
            throw std::invalid_argument("Tree-of-shapes infinity seed must be inside the interpolated domain.");
        }
        return ImageUtils::to1D(infinitySeedRow, infinitySeedCol, interpNumCols);
    }

    inline void setDiagonal0Connection(AdjacencyUC& adj, int row, int col) const {
        adj.setDiagonalConnection(row, col-1, DiagonalConnection::SE);
        adj.setDiagonalConnection(row+1, col, DiagonalConnection::NW);

        adj.setDiagonalConnection(row - 1, col - 1, DiagonalConnection::SE);
        adj.setDiagonalConnection(row, col, DiagonalConnection::SE | DiagonalConnection::NW);
        adj.setDiagonalConnection(row + 1, col + 1, DiagonalConnection::NW);

        adj.setDiagonalConnection(row-1, col, DiagonalConnection::SE);
        adj.setDiagonalConnection(row, col+1, DiagonalConnection::NW);
    }

    inline void setDiagonal1Connection(AdjacencyUC& adj, int row, int col) const {
        adj.setDiagonalConnection(row, col-1, DiagonalConnection::NE);
        adj.setDiagonalConnection(row-1, col, DiagonalConnection::SW);

        adj.setDiagonalConnection(row-1, col+1, DiagonalConnection::SW);
        adj.setDiagonalConnection(row, col, DiagonalConnection::SW | DiagonalConnection::NE);
        adj.setDiagonalConnection(row + 1, col - 1, DiagonalConnection::NE);

        adj.setDiagonalConnection(row+1, col, DiagonalConnection::NE);
        adj.setDiagonalConnection(row, col+1, DiagonalConnection::SW);
    }




	public:

	    /**
	     * @brief Creates a tree-of-shapes builder.
	     *
	     * @param interpolation Interpolation/connectivity convention.
	     * @param infinitySeedRow Row of the propagation seed in the interpolated domain.
	     * @param infinitySeedCol Column of the propagation seed in the interpolated domain.
	     */
	    explicit BuilderTreeOfShape(
	        ToSInterpolation interpolation,
        int infinitySeedRow = ToSDefaultInfinityRow,
        int infinitySeedCol = ToSDefaultInfinityCol)
	        : interpolation(interpolation), infinitySeedRow(infinitySeedRow), infinitySeedCol(infinitySeedCol) {}
	    ~BuilderTreeOfShape() { }

        /// @cond INTERNAL
	     std::tuple<std::vector<ToSLevel>, std::vector<ToSLevel>, AdjacencyUC> interpolateImage(const ImageUInt8Ptr& imgPtr) const{
        // Implements the self-dual span-based immersion used by Boutry's PhD
        // thesis: ISpan(u) followed by front propagation from a median-valued
        // outer boundary. Internal levels are stored in Z/2 by scaling values
        // by 2, so odd integers represent half gray levels.
        auto img = imgPtr->rawData();
        int numRows = imgPtr->getNumRows();
        int numCols = imgPtr->getNumCols();
        if (numRows <= 0 || numCols <= 0) {
            throw std::invalid_argument("BuilderTreeOfShape requires a non-empty image.");
        }
        if (numRows == 1 && numCols == 1) {
            const int interpNumRows = interpolatedNumRows(numRows);
            const int interpNumCols = interpolatedNumCols(numCols);
            std::vector<ToSLevel> interpolationMin(interpNumRows * interpNumCols, scaledOriginalLevel(img[0]));
            std::vector<ToSLevel> interpolationMax(interpNumRows * interpNumCols, scaledOriginalLevel(img[0]));
            AdjacencyUC adj(interpNumRows, interpNumCols, false);
            return std::make_tuple(std::move(interpolationMin), std::move(interpolationMax), std::move(adj));
        }

        constexpr int adjCircleCol[] = {-1, +1, -1, +1};
        constexpr int adjCircleRow[] = {-1, -1, +1, +1};

        constexpr int adjRetHorCol[] = {0, 0};
        constexpr int adjRetHorRow[] = {-1, +1};

        constexpr int adjRetVerCol[] = {+1, -1};
        constexpr int adjRetVerRow[] = {0, 0};

        int interpNumCols = interpolatedNumCols(numCols);
        int interpNumRows = interpolatedNumRows(numRows);
        int size = interpNumCols * interpNumRows;

        // Allocate interpolation result buffers for minimum and maximum levels.
        std::vector<ToSLevel> interpolationMin(size);
        std::vector<ToSLevel> interpolationMax(size);

        int numBoundary = 2 * (numRows + numCols) - 4;

        std::vector<int> pixels(numBoundary);

        int pT, i = 0; // Boundary-pixel counter for the median buffer.
        
        for (int p = 0; p < numCols * numRows; p++) {
            auto [row, col] = ImageUtils::to2D(p, numCols);

            // Check whether the pixel lies on the image border.
            if (row == 0 || row == numRows - 1 || col == 0 || col == numCols - 1) {
                pixels[i++] = static_cast<int>(img[p]); // Add the pixel to the median buffer.
            }

            // Compute the interpolated-image index.
            pT = ImageUtils::to1D(originalPointRow(row), originalPointCol(col), interpNumCols);

            // Assign interpolation values.
            interpolationMin[pT] = interpolationMax[pT] = scaledOriginalLevel(img[p]);
        }

        std::sort(pixels.begin(), pixels.end());
        int median;
        if (numBoundary % 2 == 0) {
            median = pixels[numBoundary / 2 - 1] + pixels[numBoundary / 2];
        } else {
            median = ToSInterpolationScale * pixels[numBoundary / 2];
        }
        //std::cout << "Interpolation (Median): " << median << std::endl;

        
        int qT, qCol, qRow, min, max;
        const int* adjCol = nullptr;
        const int* adjRow = nullptr;
        int adjSize;
        AdjacencyUC adj(interpNumRows, interpNumCols, false);

        for (int row=0; row < interpNumRows; row++){
            for (int col=0; col < interpNumCols; col++){
                if (col % 2 == 1 && row % 2 == 1) continue;
                pT = ImageUtils::to1D(row, col, interpNumCols);
                if(col == 0 || col == interpNumCols - 1 || row == 0 || row == interpNumRows - 1){
                    max = median;
                    min = median;
                }else{
                    if (col % 2 == 0 && row % 2 == 0) { 
                        adjCol = adjCircleCol;
                        adjRow = adjCircleRow;
                        adjSize = 4;
                    } else if (col % 2 == 0 && row % 2 == 1) {
                        adjCol = adjRetVerCol;
                        adjRow = adjRetVerRow;
                        adjSize = 2;
                    } else if (col % 2 == 1 && row % 2 == 0) {
                        adjCol = adjRetHorCol;
                        adjRow = adjRetHorRow;
                        adjSize = 2;
                    } 

                    min = std::numeric_limits<int>::max();
                    max = std::numeric_limits<int>::min();
                    for (int i = 0; i < adjSize; i++) {
                        qRow = row + adjRow[i];
                        qCol = col + adjCol[i];

                        if (qRow >= 0 && qCol >= 0 && qRow < interpNumRows && qCol < interpNumCols) {
                            qT = ImageUtils::to1D(qRow, qCol, interpNumCols);

                            if (interpolationMax[qT] > max) {
                                max = interpolationMax[qT];
                            }
                            if (interpolationMin[qT] < min) {
                                min = interpolationMin[qT];
                            }
                        } else {
                            if (median > max) {
                                max = median;
                            }
                            if (median < min) {
                                min = median;
                            }
                        }
                    }
                }
                interpolationMin[pT] = static_cast<ToSLevel>(min);
                interpolationMax[pT] = static_cast<ToSLevel>(max);
            }
        }
        return std::make_tuple(std::move(interpolationMin), std::move(interpolationMax), std::move(adj));
       
    }

    std::tuple<std::vector<ToSLevel>, std::vector<ToSLevel>, AdjacencyUC> interpolateImage4c8c(const ImageUInt8Ptr& imgPtr) const{
        // Implements the optimized 2D immersion/connectivity-map rules from
        // Carlinet, Crozet, and Geraud, "The Tree of Shapes Turned into a
        // Max-Tree: A Simple and Efficient Linear Algorithm", ICIP 2018.
        auto img = imgPtr->rawData();
        int numRows = imgPtr->getNumRows();
        int numCols = imgPtr->getNumCols();
        if (numRows <= 0 || numCols <= 0) {
            throw std::invalid_argument("BuilderTreeOfShape requires a non-empty image.");
        }

        int interpNumCols = interpolatedNumCols(numCols);
        int interpNumRows = interpolatedNumRows(numRows);
        int size = interpNumCols * interpNumRows;
        AdjacencyUC adj(interpNumRows, interpNumCols, true);


        // Allocate interpolation result buffers for minimum and maximum levels.
        std::vector<ToSLevel> interpolationMin(size);
        std::vector<ToSLevel> interpolationMax(size);

        int pT;
         // Compute interval from 2-faces.
        for (int p = 0; p < numCols * numRows; p++) {
            auto [row, col] = ImageUtils::to2D(p, numCols);

            // Compute the interpolated-image index.
            pT = ImageUtils::to1D(originalPointRow(row), originalPointCol(col), interpNumCols);

            // Assign interpolation values.
            interpolationMin[pT] = interpolationMax[pT] = img[p];
        }

        auto getValue = [&](int row, int col) -> int {
            int origRow = (row - 1) / 2;
            int origCol = (col - 1) / 2;
            return img[ImageUtils::to1D(origRow, origCol, numCols)];
        };

        // Borders.
        for (int row=0; row < interpNumRows; row++){
            int col;
            if(row % 2 == 1){ //horizontal e vertical
                col = 0;
                int v1 = getValue(row, col+1);
                interpolationMin[ImageUtils::to1D(row, col, interpNumCols)] = v1;
                interpolationMax[ImageUtils::to1D(row, col, interpNumCols)] = v1;

                col = interpNumCols - 1;
                v1 = getValue(row, col -1);
                interpolationMin[ImageUtils::to1D(row, col, interpNumCols)] = v1;
                interpolationMax[ImageUtils::to1D(row, col, interpNumCols)] = v1;
            }else{ //circulos
                if(row == 0){
                    col = 0;
                    int v1 = getValue(row+1, col+1);
                    interpolationMin[ImageUtils::to1D(row, col, interpNumCols)] = v1;
                    interpolationMax[ImageUtils::to1D(row, col, interpNumCols)] = v1;

                    col = interpNumCols - 1;
                    v1 = getValue(row+1, col -1);
                    interpolationMin[ImageUtils::to1D(row, col, interpNumCols)] = v1;
                    interpolationMax[ImageUtils::to1D(row, col, interpNumCols)] = v1;

                }else if(row == interpNumRows-1){
                    col = 0;
                    int v1 = getValue(row-1, 1);
                    interpolationMin[ImageUtils::to1D(row, col, interpNumCols)] = v1;
                    interpolationMax[ImageUtils::to1D(row, col, interpNumCols)] = v1;

                    col = interpNumCols - 1;
                    v1 = getValue(row-1, col - 1);
                    interpolationMin[ImageUtils::to1D(row, col, interpNumCols)] = v1;
                    interpolationMax[ImageUtils::to1D(row, col, interpNumCols)] = v1;
                }else{
                    col = 0;
                    int v1 = getValue(row-1, col+1);
                    int v2 = getValue(row+1, col+1);
                    interpolationMin[ImageUtils::to1D(row, 0, interpNumCols)] = std::min(v1, v2);
                    interpolationMax[ImageUtils::to1D(row, 0, interpNumCols)] = std::max(v1, v2);

                    col = interpNumCols - 1;
                    v1 = getValue(row-1, col-1);
                    v2 = getValue(row+1, col-1);
                    interpolationMin[ImageUtils::to1D(row, col, interpNumCols)] = std::min(v1, v2);
                    interpolationMax[ImageUtils::to1D(row, col, interpNumCols)] = std::max(v1, v2);
                }
            }
        }
        
        for (int col=1; col < interpNumCols-1; col++){
            int row;
            if(col % 2 == 1){ //horizontal e vertical
                row = 0;
                int v1 = getValue(row+1, col);
                interpolationMin[ImageUtils::to1D(row, col, interpNumCols)] = v1;
                interpolationMax[ImageUtils::to1D(row, col, interpNumCols)] = v1;

                row = interpNumRows - 1;
                v1 = getValue(row-1, col);
                interpolationMin[ImageUtils::to1D(row, col, interpNumCols)] = v1;
                interpolationMax[ImageUtils::to1D(row, col, interpNumCols)] = v1;
            }else{ //circulos
                row = 0;
                int v1 = getValue(row+1, col-1);
                int v2 = getValue(row+1, col+1);
                interpolationMin[ImageUtils::to1D(row, col, interpNumCols)] = std::min(v1, v2);
                interpolationMax[ImageUtils::to1D(row, col, interpNumCols)] = std::max(v1, v2);

                row = interpNumRows - 1;
                v1 = getValue(row-1, col-1);
                v2 = getValue(row-1, col+1);
                interpolationMin[ImageUtils::to1D(row, col, interpNumCols)] = std::min(v1, v2);
                interpolationMax[ImageUtils::to1D(row, col, interpNumCols)] = std::max(v1, v2);
            }
        }

        // Compute interval from 1-faces 
        for (int row=1; row < interpNumRows-1; row++){
            for (int col=1; col < interpNumCols-1; col++){
                if (row % 2 == 1 && col % 2 == 1) continue;  // Already defined.

                pT = ImageUtils::to1D(row, col, interpNumCols);
                if (col % 2 == 0 && row % 2 == 1) {
                    int v1 = getValue(row, col+1);
                    int v2 = getValue(row, col-1);
                    interpolationMin[pT] = std::min(v1, v2);
                    interpolationMax[pT] = std::max(v1, v2);
                } else if (col % 2 == 1 && row % 2 == 0) {
                    int v1 = getValue(row+1, col);
                    int v2 = getValue(row-1, col);
                    interpolationMin[pT] = std::min(v1, v2);
                    interpolationMax[pT] = std::max(v1, v2);
                } 
            }
        }
         // Compute interval from 0-faces 
         for (int row=1; row < interpNumRows-1; row++){
            for (int col=1; col < interpNumCols-1; col++){
                if (row % 2 == 1 && col % 2 == 1) continue;  // Already defined.
                pT = ImageUtils::to1D(row, col, interpNumCols);
                if (row % 2 == 0 && col % 2 == 0) {
                    // | v0 | v1 |
                    // | v2 | v3 |
                    int v0 = getValue(row - 1, col - 1);
                    int v1 = getValue(row + 1, col - 1);
                    int v2 = getValue(row - 1, col + 1);
                    int v3 = getValue(row + 1, col + 1);


                    const int min_v0v3 = std::min(v0, v3);
                    const int max_v0v3 = std::max(v0, v3);
                    const int min_v1v2 = std::min(v1, v2);
                    const int max_v1v2 = std::max(v1, v2);
                    const bool diagonal0IsHigh = min_v0v3 > max_v1v2;
                    const bool diagonal1IsHigh = min_v1v2 > max_v0v3;

                    if (diagonal0IsHigh || diagonal1IsHigh) {
                        const bool chooseDiagonal0 = usesHighDiagonalAtSaddle() ? diagonal0IsHigh : !diagonal0IsHigh;
                        if (chooseDiagonal0) {
                            setDiagonal0Connection(adj, row, col);
                            interpolationMin[pT] = min_v0v3;
                            interpolationMax[pT] = max_v0v3;
                        } else {
                            setDiagonal1Connection(adj, row, col);
                            interpolationMin[pT] = min_v1v2;
                            interpolationMax[pT] = max_v1v2;
                        }
                    }else{
                        // Non-critical configuration.
                        interpolationMin[pT] = std::min({v0, v1, v2, v3});
                        interpolationMax[pT] = std::max({v0, v1, v2, v3});
                    }
                }

            }
        }
        return std::make_tuple(std::move(interpolationMin), std::move(interpolationMax), std::move(adj));

       
    }

    std::tuple<std::vector<ToSLevel>, std::vector<int>, AdjacencyUC> sort(const ImageUInt8Ptr& imgPtr) const{

        int numRows = imgPtr->getNumRows();
        int numCols = imgPtr->getNumCols();

        int interpNumCols = interpolatedNumCols(numCols);
        int interpNumRows = interpolatedNumRows(numRows);
        int size = interpNumCols * interpNumRows;
        auto [interpolationMin, interpolationMax, adj] = usesConnectivityMap()? interpolateImage4c8c(imgPtr): interpolateImage(imgPtr);
        
        std::vector<uint8_t> dejavu(size, 0);  
        std::vector<int> imgR(size);  // Ordered pixels.
        std::vector<ToSLevel> imgU(size);        // Interpolated image levels.
        
        PriorityQueueToS queue(usesConnectivityMap() ? ToSUInt8Depth : ToSSelfDualDepth);  // Priority queue.
        int infinityPixel = infinitySeedIndex(interpNumRows, interpNumCols);
        int priorityQueueOld = interpolationMin[infinityPixel];
        queue.initial(infinityPixel, priorityQueueOld);
        dejavu[infinityPixel] = true;

        int order = 0; 
        int depth = 0;
        while (!queue.isEmpty()) {
            int h = queue.priorityPop();  // Pop the element with highest priority.
            int priorityQueue = queue.getCurrentPriority(); // Current priority.
            if(usesConnectivityMap()){
                if(priorityQueue != priorityQueueOld) depth++;
                imgU[h] = depth;
            }else{
                imgU[h] = priorityQueue;
            }
            
            // Store h in the correct output order.
            imgR[order++] = h;
            
            // Adjacencies.
            for(int n: adj.getNeighborPixels(h)){
                if (!dejavu[n]) {
                    queue.priorityPush(n, interpolationMin[n], interpolationMax[n]);
                    dejavu[n] = true;  // Mark as processed.
                }
            }
            priorityQueueOld = priorityQueue;
        }
        return std::make_tuple(std::move(imgU), std::move(imgR), std::move(adj));
    }


    // Tests whether an interpolated pixel corresponds to an original pixel.
    inline bool isOriginal1D(int p, int interpNumCols) const{
        /*
        p = row x  numCols + col
        numCols is odd. If row is odd, row * numCols is odd.
        Therefore row and col are odd iff p is even and row is odd.
        */
        int row = p / interpNumCols;
        
        // original <=> p is even and row is odd.
        return ((p & 1) == 0) && ((row & 1) == 1);
    }

    // Maps an interpolated pixel to the original image domain.
    inline int toOriginal1D(int pStar, int interNumCols, int numCols) const{
        int r = pStar / interNumCols;
        int c = pStar - r * interNumCols;         // Avoid the modulo operator.
        return ((r - 1) >> 1) * numCols + ((c - 1) >> 1);
    }
        /// @endcond

	    /**
	     * @brief Builds the projected tree-of-shapes parent representation.
	     *
	     * @return Tuple containing the parent image on the original domain, the
	     * projected pixel order, and the number of created tree nodes.
	     */
	    std::tuple<std::vector<int>, std::vector<int>, int> createTreeByUnionFind(const ImageUInt8Ptr& imgPtr) const override {

        int numRows = imgPtr->getNumRows();
        int numCols = imgPtr->getNumCols();
        int numPixels = numRows * numCols;

        int interpNumCols = interpolatedNumCols(numCols);
        int interpNumRows = interpolatedNumRows(numRows);
        int numPixelsInterp = interpNumCols * interpNumRows;
        auto [imgInterpolate, orderedPixelsInterpolete, adj] = sort(imgPtr);

        // Union-find over the interpolated image.
        std::vector<int> zPar(numPixelsInterp, InvalidNode);
        std::vector<int> parentInterpolate(numPixelsInterp, InvalidNode);
        auto findRoot = [&](int pStar) {
            while (zPar[pStar] != pStar) { zPar[pStar] = zPar[zPar[pStar]]; pStar = zPar[pStar]; }
            return pStar;
        };

        for (int i = numPixelsInterp - 1; i >= 0; --i) { // Leaves to root.
            int pStar = orderedPixelsInterpolete[i];
            parentInterpolate[pStar] = pStar;
            zPar[pStar] = pStar;
            for (int qStar : adj.getNeighborPixels(pStar)) {
                if (zPar[qStar] != -1) {
                    int rStar = findRoot(qStar);
                    if (pStar != rStar) { parentInterpolate[rStar] = pStar; zPar[rStar] = pStar; }
                }
            }
        }

        auto sameLevel = [&](int aStar, int bStar){ 
            return imgInterpolate[aStar] == imgInterpolate[bStar]; 
        };
        auto repOf = [&](int pStar) { // Plateau representative: parent at same level, otherwise the point itself.
            int parStar = parentInterpolate[pStar];
            return (parStar == pStar || imgInterpolate[parStar] == imgInterpolate[pStar]) ? parStar : pStar;
        };
        
        // Step 1 - canonicalization and representative marking.
        int numNodes=0;
        std::vector<int> representativeOriginalPixels(numPixelsInterp, InvalidNode); // Valid only for plateau representatives.
        for (int i = 0; i < numPixelsInterp; ++i) {  // Root to leaves.
            int pStar = orderedPixelsInterpolete[i];
            int qStar = parentInterpolate[pStar];
            
            // Canonicalize one step inside the plateau.
            if (sameLevel(parentInterpolate[qStar], qStar))
                parentInterpolate[pStar] = parentInterpolate[qStar];
            
            if (parentInterpolate[pStar] == pStar || imgInterpolate[parentInterpolate[pStar]] != imgInterpolate[pStar])
                ++numNodes;

            if (isOriginal1D(pStar, interpNumCols)) {
                int rep = repOf(pStar);              // Computed on demand.
                if (representativeOriginalPixels[rep] == InvalidNode) representativeOriginalPixels[rep] = pStar;  // Select the first original point in the plateau.
            }
        }

        auto representativeOriginalAbove = [&](int repStar) {
            int currentRep = repStar;
            while (currentRep != InvalidNode) {
                if (representativeOriginalPixels[currentRep] != InvalidNode) {
                    return representativeOriginalPixels[currentRep];
                }
                int parentStar = parentInterpolate[currentRep];
                if (parentStar == currentRep) {
                    break;
                }
                currentRep = repOf(parentStar);
            }
            return InvalidNode;
        };


        std::vector<int> parent(numPixels, InvalidNode);
        std::vector<int> orderedPixels; orderedPixels.reserve(numPixels);
        int projectedRootOriginal = InvalidNode;
        for (int i = 0; i < numPixelsInterp; ++i) { // Root to leaves.
            int pStar = orderedPixelsInterpolete[i];
            if (!isOriginal1D(pStar, interpNumCols)) continue;

            int parStar   = parentInterpolate[pStar];
            int qStar;
            if (parStar == pStar) { // Tree root.
                qStar = pStar;
            }
            else if (sameLevel(parStar, pStar)) {
                int repP = repOf(pStar);
                if (representativeOriginalPixels[repP] == pStar) {
                    // pStar is the original representative of its plateau; its parent comes from the plateau above.
                    int repAbove = repOf(parentInterpolate[repP]);
                    qStar = representativeOriginalAbove(repAbove);
                } else {
                    // Same plateau: parent is the original representative of that plateau.
                    qStar = representativeOriginalPixels[repP];
                }
            }
            else {
                // Different level: parent is the original representative of the parent's plateau.
                int repPar = repOf(parStar);
                qStar = representativeOriginalAbove(repPar);
            }

            if (qStar == InvalidNode) {
                if (projectedRootOriginal == InvalidNode) {
                    projectedRootOriginal = pStar;
                }
                qStar = projectedRootOriginal;
            }

            // Projection to the original parent and ordered-pixel buffers.
            int p = toOriginal1D(pStar, interpNumCols, numCols);
            int q = toOriginal1D(qStar,  interpNumCols, numCols);
            parent[p] = q;
            orderedPixels.push_back(p);
        }
        

        return std::make_tuple(std::move(parent), std::move(orderedPixels), numNodes);
    }

};

} // namespace mmcfilters
