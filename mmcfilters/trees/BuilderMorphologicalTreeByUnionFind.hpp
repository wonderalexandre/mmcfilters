#pragma once

#include "../utils/AdjacencyRelation.hpp"
#include "../utils/Common.hpp"


namespace mmcfilters {

/**
 * @brief Interface para construtores de árvores morfológicas baseados em union-find.
 */
class IMorphologicalTreeBuilder {
public:
    virtual ~IMorphologicalTreeBuilder() = default;
    
    virtual std::tuple<std::vector<int>, std::vector<int>, int> createTreeByUnionFind(const ImageUInt8Ptr& imgPtr) const = 0;
};


/**
 * @brief Implementa construção de component trees via algoritmo union-find.
 */
class BuilderComponentTree : public IMorphologicalTreeBuilder{
private:
    AdjacencyRelation* adj;
    bool isMaxtree;

public:
    explicit BuilderComponentTree(AdjacencyRelation* adj, bool isMaxtree) : adj(adj), isMaxtree(isMaxtree) { }
    ~BuilderComponentTree() { }

    template <typename PixelType>
    std::vector<int> sort(ImagePtr<PixelType> imgPtr) const {
        const int n = imgPtr->getSize();
        std::vector<int> orderedPixels(n);
        PixelType* img = imgPtr->rawData();

        if constexpr (std::is_floating_point_v<PixelType>) {
            if (PRINT_LOG) std::cout << "Sorting floating point image with size: " << n << std::endl;
            std::iota(orderedPixels.begin(), orderedPixels.end(), 0);
            if (isMaxtree) {
                std::stable_sort(orderedPixels.begin(), orderedPixels.end(),
                                 [&](int a, int b) { return img[a] < img[b]; });
            } else {
                std::stable_sort(orderedPixels.begin(), orderedPixels.end(),
                                 [&](int a, int b) { return img[a] > img[b]; });
            }
        } else {
            if (PRINT_LOG) std::cout << "Sorting integer image with size: " << n << std::endl;
        
            // counting sort com faixa [0..maxvalue]; 
            int maxvalue =  static_cast<int>(img[0]);
            for (int i = 1; i < n; i++) if(maxvalue < img[i]) maxvalue = img[i];
            std::vector<uint32_t> counter(static_cast<size_t>(maxvalue) + 1, 0);
            if(isMaxtree){
                for (int i = 0; i < n; i++)
                    counter[img[i]]++;
                for (int i = 1; i < maxvalue; i++) 
                    counter[i] += counter[i - 1];
                counter[maxvalue] += counter[maxvalue-1];
                for (int i = n - 1; i >= 0; --i)
                    orderedPixels[--counter[img[i]]] = i;	

            }else{
                for (int i = 0; i < n; i++)
                    counter[maxvalue - img[i]]++;
                for (int i = 1; i < maxvalue; i++) 
                    counter[i] += counter[i - 1];
                counter[maxvalue] += counter[maxvalue-1];
                for (int i = n - 1; i >= 0; --i)
                    orderedPixels[--counter[maxvalue - img[i]]] = i;
            }
        }
        return orderedPixels;
    }

    std::tuple<std::vector<int>, std::vector<int>, int> createTreeByUnionFind(const ImageUInt8Ptr& imgPtr) const override {
        return createTreeByUnionFind<uint8_t>(imgPtr);
    }

    
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



/************************ToS******************* */

/*
Adjacência adaptativa para construção da ToS com 4 e 8-conectividade.
As conexões diagonais são ativadas conforme necessário durante a construção da ToS para garantir a conectividade adequada dos componentes.
As conexões diagonais são representadas usando flags de bits em um vetor. Cada pixel pode ter até quatro conexões diagonais: SW, NE, SE, NW
Essas conexões são armazenadas em um vetor de uint8_t, onde cada bit representa uma conexão diagonal específica.
Por exemplo, se o bit 0 estiver definido, significa que há uma conexão diagonal SW do pixel atual para o pixel ao sul-oeste.
*/
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

/**
 * @brief Representa adjacência adaptativa para construção da ToS de 4/8-conectividade.
 *
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
     * @brief Iterador sobre os vizinhos válidos considerando conexões diagonais.
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
     * @brief Range helper que produz iteradores de vizinhança válidos.
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
 * @brief Fila de prioridades discreta utilizada durante a construção da ToS.
 */
class PriorityQueueToS {
private:
    std::vector<std::deque<int>> buckets;
    int currentPriority;
    int numElements;
    int maxPriorityLevels;
    

public:
    PriorityQueueToS(int depthOfImage=8) : currentPriority(0), numElements(0), maxPriorityLevels(1 << depthOfImage){
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
        // Se o bucket atual estiver vazio, precisamos ajustar a prioridade
        if (buckets[currentPriority].empty()) {
            int i = currentPriority;
            int j = currentPriority;
            while (true) {

                // Tentar diminuir a prioridade
                if (j > 0 && buckets[j].empty()) {
                    j--;
                }
                if (!buckets[j].empty()) { // Encontrou o próximo bucket não vazio diminuindo a prioridade
                    currentPriority = j;
                    break;
                }

                // Tentar aumentar a prioridade
                if (i < maxPriorityLevels && buckets[i].empty()) {
                    i++;
                }
                if (i < maxPriorityLevels && !buckets[i].empty()) { // Encontrou o próximo bucket não vazio aumentando a prioridade
                    currentPriority = i;
                    break;
                }
            }
        }

        int element = buckets[currentPriority].front(); 
        buckets[currentPriority].pop_front();           

        numElements--;  
        return element;
    }
};




/**
 * @brief Constrói árvores de formas (ToS) usando o algoritmo de Thierry Géraud et al.
 */
class BuilderTreeOfShape: public IMorphologicalTreeBuilder {
private:
    bool is4c8cConnectivity;




public:

    explicit BuilderTreeOfShape(): BuilderTreeOfShape(true) {}
    explicit BuilderTreeOfShape(bool is4c8cConnectivity): is4c8cConnectivity(is4c8cConnectivity) {}
    ~BuilderTreeOfShape() { }

     /**
      * Implementation based on the paper: 
      *  - Thesi of the N.Boutry
      * - T. Géraud, E. Carlinet, and S. Crozet, Self-Duality and Digital Topology: Links Between the Morphological Tree of Shapes and Well-Composed Gray-Level Images, ISMM 2015
      * - N.Boutry, T.Géraud, L.Najman, "How to Make nD Functions Digitally Well-Composed in a Self-dual Way", ISMM 2015.
      * - N.Boutry, T.Géraud, L.Najman, "On Making {$n$D} Images Well-Composed by a Self-Dual Local Interpolation", DGCI 2014
      */
     std::tuple<std::vector<uint8_t>, std::vector<uint8_t>, AdjacencyUC> interpolateImage(const ImageUInt8Ptr& imgPtr) const{
        auto img = imgPtr->rawData();
        int numRows = imgPtr->getNumRows();
        int numCols = imgPtr->getNumCols();

        constexpr int adjCircleCol[] = {-1, +1, -1, +1};
        constexpr int adjCircleRow[] = {-1, -1, +1, +1};

        constexpr int adjRetHorCol[] = {0, 0};
        constexpr int adjRetHorRow[] = {-1, +1};

        constexpr int adjRetVerCol[] = {+1, -1};
        constexpr int adjRetVerRow[] = {0, 0};

        int interpNumCols = numCols * 2 + 1;
        int interpNumRows = numRows * 2 + 1;
        int size = interpNumCols * interpNumRows;

        // Aloca memória para os resultados de interpolação (mínimo e máximo)
        std::vector<uint8_t> interpolationMin(size);
        std::vector<uint8_t> interpolationMax(size);
        //this->imgU = std::make_unique<uint8_t[]>(size);
        //this->interpolationMax = std::make_unique<uint8_t[]>(size);

        int numBoundary = 2 * (numRows + numCols) - 4;
        
        //std::unique_ptr<uint8_t[]> pixelsPtr(new uint8_t[numBoundary]);// Para calcular a mediana
        std::vector<uint8_t> pixels(numBoundary);
        //uint8_t* pixels = pixelsPtr.get();

        int pT, i = 0; // i é um contador para o array pixels
        
        for (int p = 0; p < numCols * numRows; p++) {
            auto [row, col] = ImageUtils::to2D(p, numCols);

            // Verifica se o pixel está na borda
            if (row == 0 || row == numRows - 1 || col == 0 || col == numCols - 1) {
                pixels[i++] = img[p]; // Adiciona o pixel ao array pixels
            }

            // Calcula o índice para imagem interpolada
            pT = ImageUtils::to1D(2 * row + 1, 2 * col + 1, interpNumCols);

            // Define os valores de interpolação
            interpolationMin[pT] = interpolationMax[pT] = img[p];
        }

        //std::sort(pixels, pixels + numBoundary);
        std::sort(pixels.begin(), pixels.end());
        int median;
        if (numBoundary % 2 == 0) {
            median = (pixels[numBoundary / 2 - 1] + pixels[numBoundary / 2]) / 2;
        } else {
            median = pixels[numBoundary / 2];
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
                interpolationMin[pT] = min;
                interpolationMax[pT] = max;
            }
        }
        return std::make_tuple(std::move(interpolationMin), std::move(interpolationMax), std::move(adj));
       
    }

    std::tuple<std::vector<uint8_t>, std::vector<uint8_t>, AdjacencyUC> interpolateImage4c8c(const ImageUInt8Ptr& imgPtr) const{
        auto img = imgPtr->rawData();
        int numRows = imgPtr->getNumRows();
        int numCols = imgPtr->getNumCols();

        int interpNumCols = numCols * 2 + 1;
        int interpNumRows = numRows * 2 + 1;
        int size = interpNumCols * interpNumRows;
        AdjacencyUC adj(interpNumRows, interpNumCols, true);


        // Aloca memória para os resultados de interpolação (mínimo e máximo)
        std::vector<uint8_t> interpolationMin(size);
        std::vector<uint8_t> interpolationMax(size);

        int pT;
         // Compute interval from 2-faces.
        for (int p = 0; p < numCols * numRows; p++) {
            auto [row, col] = ImageUtils::to2D(p, numCols);

            // Calcula o índice para imagem interpolada
            pT = ImageUtils::to1D(2 * row + 1, 2 * col + 1, interpNumCols);

            // Define os valores de interpolação
            interpolationMin[pT] = interpolationMax[pT] = img[p];
        }

        auto getValue = [&](int row, int col) -> int {
            int origRow = (row - 1) / 2;
            int origCol = (col - 1) / 2;
            return img[ImageUtils::to1D(origRow, origCol, numCols)];
        };

        // Bordas
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
                if (row % 2 == 1 && col % 2 == 1) continue;  // já definido

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
                if (row % 2 == 1 && col % 2 == 1) continue;  // já definido
                pT = ImageUtils::to1D(row, col, interpNumCols);
                if (row % 2 == 0 && col % 2 == 0) {
                    // | v0 | v1 |
                    // | v2 | v3 |
                    int v0 = getValue(row - 1, col - 1);
                    int v1 = getValue(row + 1, col - 1);
                    int v2 = getValue(row - 1, col + 1);
                    int v3 = getValue(row + 1, col + 1);


                    int min_v0v3 = std::min(v0, v3);
                    int max_v0v3 = std::max(v0, v3);
                    int min_v1v2 = std::min(v1, v2);
                    int max_v1v2 = std::max(v1, v2);
                    if (max_v1v2 > min_v0v3) {
                        
                        // Saddle point configuration 1
                        adj.setDiagonalConnection(row, col-1, DiagonalConnection::SE);
                        adj.setDiagonalConnection(row+1, col, DiagonalConnection::NW);
                        
                        adj.setDiagonalConnection(row - 1, col - 1, DiagonalConnection::SE);
                        adj.setDiagonalConnection(row, col, DiagonalConnection::SE | DiagonalConnection::NW);
                        adj.setDiagonalConnection(row + 1, col + 1, DiagonalConnection::NW);

                        adj.setDiagonalConnection(row-1, col, DiagonalConnection::SE);
                        adj.setDiagonalConnection(row, col+1, DiagonalConnection::NW);

                        interpolationMin[pT] = min_v0v3;
                        interpolationMax[pT] = max_v0v3;
                    }
                    else if (max_v0v3 > min_v1v2) {
                        // Saddle point configuration 2
                        adj.setDiagonalConnection(row, col-1, DiagonalConnection::NE);
                        adj.setDiagonalConnection(row-1, col, DiagonalConnection::SW);

                        adj.setDiagonalConnection(row-1, col+1, DiagonalConnection::SW);
                        adj.setDiagonalConnection(row, col, DiagonalConnection::SW | DiagonalConnection::NE);
                        adj.setDiagonalConnection(row + 1, col - 1, DiagonalConnection::NE);

                        adj.setDiagonalConnection(row+1, col, DiagonalConnection::NE);
                        adj.setDiagonalConnection(row, col+1, DiagonalConnection::SW);

                        interpolationMin[pT] = min_v1v2;
                        interpolationMax[pT] = max_v1v2;
                    }else{
                        // Non-critical configuration.
                        interpolationMin[pT] = std::min(min_v0v3, min_v1v2);
                        interpolationMax[pT] = std::min(max_v0v3, max_v1v2);
                    }
                }

            }
        }
        return std::make_tuple(std::move(interpolationMin), std::move(interpolationMax), std::move(adj));

       
    }

    std::tuple<std::vector<uint8_t>, std::vector<int>, AdjacencyUC> sort(const ImageUInt8Ptr& imgPtr) const{

        int numRows = imgPtr->getNumRows();
        int numCols = imgPtr->getNumCols();

        int interpNumCols = numCols * 2 + 1;
        int interpNumRows = numRows * 2 + 1;
        int size = interpNumCols * interpNumRows;
        auto [interpolationMin, interpolationMax, adj] = is4c8cConnectivity? interpolateImage4c8c(imgPtr): interpolateImage(imgPtr);
        
        std::vector<uint8_t> dejavu(size, 0);  
        std::vector<int> imgR(size);  // Pixels ordenados
        std::vector<uint8_t> imgU(size);        // Níveis de cinza da imagem
        
        PriorityQueueToS queue;  // Fila de prioridade
        int pInfinito = ImageUtils::to1D(1, 1, interpNumCols);
        int priorityQueueOld = interpolationMin[pInfinito];
        queue.initial(pInfinito, priorityQueueOld);  
        dejavu[pInfinito] = true;

        int order = 0; 
        int depth = 0;
        while (!queue.isEmpty()) {
            int h = queue.priorityPop();  // Retirar o elemento com maior prioridade
            int priorityQueue = queue.getCurrentPriority(); // Prioridade corrente
            if(is4c8cConnectivity){
                if(priorityQueue != priorityQueueOld) depth++;
                imgU[h] = depth;
            }else{
                imgU[h] = priorityQueue;
            }
            
            // Armazenar o índice h em imgR na ordem correta
            imgR[order++] = h;
            
            // Adjacências
            for(int n: adj.getNeighborPixels(h)){
                if (!dejavu[n]) {
                    queue.priorityPush(n, interpolationMin[n], interpolationMax[n]);
                    dejavu[n] = true;  // Marcar como processado
                }
            }
            priorityQueueOld = priorityQueue;
        }
        return std::make_tuple(std::move(imgU), std::move(imgR), std::move(adj));
    }


    //Testa se é um pixel é original 
    inline bool isOriginal1D(int p, int interpNumCols) const{
        /*
        p = row x  numCols + col
        Sabemos que numCols é ímpar. 
        - Testaremos se row e col são ímpares assim:
        Se row é ímpar, então row x numCols é ímpar. 
        Logo, row e col são ímpares sse p é par E col é ímpar
        */
        int row = p / interpNumCols;
        
        // original <=> (p é par) ∧ (row é ímpar)
        return ((p & 1) == 0) && ((row & 1) == 1);
    }

    // mapeia pixel interpolado para pixel original
    inline int toOriginal1D(int pStar, int interNumCols, int numCols) const{
        int r = pStar / interNumCols;
        int c = pStar - r * interNumCols;         // evita operador %
        return ((r - 1) >> 1) * numCols + ((c - 1) >> 1);
    }

    std::tuple<std::vector<int>, std::vector<int>, int> createTreeByUnionFind(const ImageUInt8Ptr& imgPtr) const override {

        int numRows = imgPtr->getNumRows();
        int numCols = imgPtr->getNumCols();
        int numPixels = numRows * numCols;

        int interpNumCols = numCols * 2 + 1;
        int interpNumRows = numRows * 2 + 1;
        int numPixelsInterp = interpNumCols * interpNumRows;
        auto [imgInterpolate, orderedPixelsInterpolete, adj] = sort(imgPtr);

        // ---------- UF em imagem interpolada ----------
        std::vector<int> zPar(numPixelsInterp, InvalidNode);
        std::vector<int> parentInterpolate(numPixelsInterp, InvalidNode);
        auto findRoot = [&](int pStar) {
            while (zPar[pStar] != pStar) { zPar[pStar] = zPar[zPar[pStar]]; pStar = zPar[pStar]; }
            return pStar;
        };

        for (int i = numPixelsInterp - 1; i >= 0; --i) { // folhas -> raiz
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
        auto repOf = [&](int pStar) { // Representante do platô de s: se o pai está no mesmo nível, o rep é o pai; senão, é ele próprio.
            int parStar = parentInterpolate[pStar];
            return (parStar == pStar || imgInterpolate[parStar] == imgInterpolate[pStar]) ? parStar : pStar;
        };
        
        // Passo 1 — canonização + marcar apenas representantes
        int numNodes=0;
        std::vector<int> repPixelsOriginais(numPixelsInterp, InvalidNode); // válido só quando o índice é representante de platô
        for (int i = 0; i < numPixelsInterp; ++i) {  // raiz -> folhas
            int pStar = orderedPixelsInterpolete[i];
            int qStar = parentInterpolate[pStar];
            
            // canoniza 1 passo dentro do platô
            if (sameLevel(parentInterpolate[qStar], qStar))
                parentInterpolate[pStar] = parentInterpolate[qStar];
            
            if (parentInterpolate[pStar] == pStar || imgInterpolate[parentInterpolate[pStar]] != imgInterpolate[pStar])
                ++numNodes;

            if (isOriginal1D(pStar, interpNumCols)) {
                int rep = repOf(pStar);              // calculado na hora
                if (repPixelsOriginais[rep] == InvalidNode) repPixelsOriginais[rep] = pStar;  // elege o o primeiro original do platô
            }
        }


        std::vector<int> parent(numPixels, InvalidNode);
        std::vector<int> orderedPixels; orderedPixels.reserve(numPixels);
        for (int i = 0; i < numPixelsInterp; ++i) { // raiz -> folhas
            int pStar = orderedPixelsInterpolete[i];
            if (!isOriginal1D(pStar, interpNumCols)) continue;

            int parStar   = parentInterpolate[pStar];
            int qStar;
            if (parStar == pStar) { // raiz da arvore
                qStar = pStar;
            }
            else if (sameLevel(parStar, pStar)) {
                int repP = repOf(pStar);
                if (repPixelsOriginais[repP] == pStar) {
                    // pStar é o representante ORIGINAL do seu platô. Seu pai vem do platô ACIMA
                    int repAbove = repOf(parentInterpolate[repP]);
                    qStar = repPixelsOriginais[repAbove];
                } else {
                    // mesmo platô. Seu pai é o representante ORIGINAL do mesmo platô
                    qStar = repPixelsOriginais[repP];
                }
            }
            else {
                // nível diferente. Seu pai é o representante ORIGINAL do platô do pai
                int repPar = repOf(parStar);
                qStar = repPixelsOriginais[repPar];
            }

            //projeção para o parent e orderedPixels originais
            int p = toOriginal1D(pStar, interpNumCols, numCols);
            int q = toOriginal1D(qStar,  interpNumCols, numCols);
            parent[p] = q;
            orderedPixels.push_back(p);
        }
        

        return std::make_tuple(std::move(parent), std::move(orderedPixels), numNodes);
    }

};

} // namespace mmcfilters

