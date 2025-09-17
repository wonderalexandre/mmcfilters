#ifndef BUILDER_MORPHOLOGICAL_TREE_BY_UNION_FIND_HPP
#define BUILDER_MORPHOLOGICAL_TREE_BY_UNION_FIND_HPP

#include "../include/ComponentTree.hpp"
#include "../include/AdjacencyRelation.hpp"
#include "../include/Common.hpp"
#include "../include/AdjacencyUC.hpp"

// Builder externo para construir ComponentTree via Union-Find (Pixels-only).
class BuilderComponentTreeByUnionFind {
private:
    

public:
    explicit BuilderComponentTreeByUnionFind() {}

    // Ordenação estável dos pixels por nível de cinza
    std::vector<int> countingSort(const ImageUInt8Ptr& imgPtr, bool isMaxtree) {
        int n = imgPtr->getSize();
        auto img = imgPtr->rawData();
        int maxvalue = img[0];
        for (int i = 1; i < n; i++) if (maxvalue < img[i]) maxvalue = img[i];

        std::vector<uint32_t> counter(maxvalue + 1, 0);
        std::vector<int> orderedPixels(n);

        if (isMaxtree) {
            for (int i = 0; i < n; i++) counter[img[i]]++;
            for (int i = 1; i < maxvalue; i++) counter[i] += counter[i - 1];
            counter[maxvalue] += counter[maxvalue - 1];
            for (int i = n - 1; i >= 0; --i) orderedPixels[--counter[img[i]]] = i;
        } else {
            for (int i = 0; i < n; i++) counter[maxvalue - img[i]]++;
            for (int i = 1; i < maxvalue; i++) counter[i] += counter[i - 1];
            counter[maxvalue] += counter[maxvalue - 1];
            for (int i = n - 1; i >= 0; --i) orderedPixels[--counter[maxvalue - img[i]]] = i;
        }
        return orderedPixels;
    }

    // Pixels: fluxo padrão por UF
    std::tuple<std::vector<int>, std::vector<int>, int> createTreeByUnionFind(const ImageUInt8Ptr& imgPtr, bool isMaxtree, AdjacencyRelation* adj) {
        std::vector<int> orderedPixels = countingSort(imgPtr, isMaxtree);
        auto img = imgPtr->rawData();

        int numPixels = imgPtr->getSize();
        std::vector<int> zPar(numPixels, -1);
        std::vector<int> parent(numPixels, -1);
        auto findRoot = [&](int p) {
            while (zPar[p] != p) { zPar[p] = zPar[zPar[p]]; p = zPar[p]; }
            return p;
        };
        

        for (int i = numPixels - 1; i >= 0; i--) {
            int p = orderedPixels[i];
            parent[p] = p;
            zPar[p] = p;
            for (int q : adj->getNeighborPixels(p)) {
                if (zPar[q] != -1) {
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
        return std::make_tuple(parent, orderedPixels, numNodes);
    }
};






class BuilderTreeOfShapeByUnionFind {
private:

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
    

public:


    BuilderTreeOfShapeByUnionFind(){}
    ~BuilderTreeOfShapeByUnionFind() { }

     /**
      * Implementation based on the paper: 
      *  - Thesi of the N.Boutry
      * - T. Géraud, E. Carlinet, and S. Crozet, Self-Duality and Digital Topology: Links Between the Morphological Tree of Shapes and Well-Composed Gray-Level Images, ISMM 2015
      * - N.Boutry, T.Géraud, L.Najman, "How to Make nD Functions Digitally Well-Composed in a Self-dual Way", ISMM 2015.
      * - N.Boutry, T.Géraud, L.Najman, "On Making {$n$D} Images Well-Composed by a Self-Dual Local Interpolation", DGCI 2014
      */
     std::tuple<std::vector<uint8_t>, std::vector<uint8_t>, AdjacencyUC> interpolateImage(const ImageUInt8Ptr& imgPtr) {
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

                    min = INT_MAX;
                    max = INT_MIN;
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
        return std::make_tuple(interpolationMin, interpolationMax, adj);
       
    }

    std::tuple<std::vector<uint8_t>, std::vector<uint8_t>, AdjacencyUC> interpolateImage4c8c(const ImageUInt8Ptr&  imgPtr) {
        auto img = imgPtr->rawData();
        int numRows = imgPtr->getNumRows();
        int numCols = imgPtr->getNumCols();

        bool is4c8cConnectivity = true;
        int interpNumCols = numCols * 2 + 1;
        int interpNumRows = numRows * 2 + 1;
        int size = interpNumCols * interpNumRows;
        AdjacencyUC adj(interpNumRows, interpNumCols, true);


        // Aloca memória para os resultados de interpolação (mínimo e máximo)
        //this->interpolationMin = std::make_unique<uint8_t[]>(size);
        //this->interpolationMax = std::make_unique<uint8_t[]>(size);
        std::vector<uint8_t> interpolationMin(size);
        std::vector<uint8_t> interpolationMax(size);


        int pT, i = 0; // i é um contador para o array pixels
        
         // Compute interval from 2-faces.
        for (int p = 0; p < numCols * numRows; p++) {
            auto [row, col] = ImageUtils::to2D(p, numCols);

            // Calcula o índice para imagem interpolada
            pT = ImageUtils::to1D(2 * row + 1, 2 * col + 1, interpNumCols);

            // Define os valores de interpolação
            interpolationMin[pT] = interpolationMax[pT] = img[p];
        }

        int qT, qCol, qRow, min, max;
        const int* adjCol = nullptr;
        const int* adjRow = nullptr;
        int adjSize;

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
        return std::make_tuple(interpolationMin, interpolationMax, adj);

       
    }

    std::tuple<std::vector<uint8_t>, std::vector<int>, AdjacencyUC> sort(const ImageUInt8Ptr& imgPtr, bool is4c8cConnectivity) {

        auto img = imgPtr->rawData();
        int numRows = imgPtr->getNumRows();
        int numCols = imgPtr->getNumCols();

        int interpNumCols = numCols * 2 + 1;
        int interpNumRows = numRows * 2 + 1;
        int size = interpNumCols * interpNumRows;
        auto [interpolationMin, interpolationMax, adj] = is4c8cConnectivity? interpolateImage4c8c(imgPtr): interpolateImage(imgPtr);

        
        //std::unique_ptr<bool[]> dejavu(new bool[size]());  // Vetor de booleanos, inicializado com false
        //this->imgR = std::make_unique<int[]>(size);  // Pixels ordenados
        //this->imgU = std::make_unique<uint8_t[]>(size);        // Níveis de cinza da imagem
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
    inline bool isOriginal1D(int p, int interpNumCols) {
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
    inline int toOriginal1D(int pStar, int interNumCols, int numCols) {
        int r = pStar / interNumCols;
        int c = pStar - r * interNumCols;         // evita operador %
        return ((r - 1) >> 1) * numCols + ((c - 1) >> 1);
    }


    std::tuple<std::vector<int>, std::vector<int>, int> createTreeByUnionFind(const ImageUInt8Ptr& imgPtr, bool is4c8cConnectivity) {

        auto img = imgPtr->rawData();
        int numRows = imgPtr->getNumRows();
        int numCols = imgPtr->getNumCols();
        int numPixels = numRows * numCols;

        int interpNumCols = numCols * 2 + 1;
        int interpNumRows = numRows * 2 + 1;
        int numPixelsInterp = interpNumCols * interpNumRows;
        auto [imgInterpolate, imgOrderedInterpolete, adj] = sort(imgPtr, is4c8cConnectivity);

        std::vector<int> zPar(numPixelsInterp, -1);
        std::vector<int> parentInterpolate(numPixelsInterp, -1);
        auto findRoot = [&](int p) {
            while (zPar[p] != p) { zPar[p] = zPar[zPar[p]]; p = zPar[p]; }
            return p;
        };

        //Construção do parentInterpolate por union-find
        for (int i = numPixelsInterp - 1; i >= 0; --i) { //processamento da construção: folhas para raiz
            int pStar = imgOrderedInterpolete[i];
            parentInterpolate[pStar] = pStar;
            zPar[pStar]   = pStar;
            for (int qStar : adj.getNeighborPixels(pStar)) {
                if (zPar[qStar] != -1) {
                    int rStar = findRoot(qStar);
                    if (pStar != rStar) { 
                        parentInterpolate[rStar] = pStar; 
                        zPar[rStar] = pStar; 
                    }
                }
            }
        }
        
        // Passo 1 — canonização + marcar apenas representantes
        std::vector<int> plateauRep(numPixelsInterp, -1);
        int numNodesInterpolateTree=0;
        for (int i = 0; i < numPixelsInterp; ++i) { // raiz -> folhas
            int p = imgOrderedInterpolete[i];
            int q = parentInterpolate[p];
            if (imgInterpolate[parentInterpolate[q]] == imgInterpolate[q])
                parentInterpolate[p] = parentInterpolate[q];   // canoniza 1 passo    
            // representante do platô de p (com 1 passo já vale em raiz->folhas):
            plateauRep[p] = (parentInterpolate[p] == p || imgInterpolate[parentInterpolate[p]] != imgInterpolate[p]) ? p : parentInterpolate[p];
            if (parentInterpolate[p] == p || imgInterpolate[parentInterpolate[p]] != imgInterpolate[p]) ++numNodesInterpolateTree;
        }

        std::vector<int> repOrig(numPixelsInterp, -1); // só posições == reps serão usadas
        for (int i = 0; i < numPixelsInterp; ++i) { // raiz -> folhas
            int p = imgOrderedInterpolete[i];
            if (!isOriginal1D(p, interpNumCols)) continue;

            int rep = plateauRep[p];
            if (repOrig[rep] == -1) repOrig[rep] = p;         // “primeiro original vence”
            // (se preferir “último vence”, basta sempre sobrescrever)
        }


        std::vector<int> parent(numPixels, -1);
        std::vector<int> imgOrdered;
        imgOrdered.reserve(numPixels);

        for (int i = 0; i < numPixelsInterp; ++i) { // raiz -> folhas
            int pStar = imgOrderedInterpolete[i];
            if (!isOriginal1D(pStar, interpNumCols)) continue;

            int repP   = plateauRep[pStar];        // rep do platô de p*
            int par    = parentInterpolate[pStar]; // pai em ToS interpolada
            int repPar = plateauRep[par];          // rep do platô do pai

            int qStar;                             // pai ORIGINAL de p*
            bool sameLevel = (imgInterpolate[par] == imgInterpolate[pStar]);
            bool isRoot    = (par == pStar);

            if (isRoot) {
                // raiz global
                qStar = pStar;
            } else if (sameLevel) {
                if (repOrig[repP] == pStar) {
                    // p* é o representante-ORIGINAL do seu platô -> o pai deve ser o rep-ORIGINAL do platô ACIMA,
                    // não o mesmo. Subimos do representante do platô até o representante do platô do pai.
                    int tRep = repP;
                    // sobe até o platô acima (nível diferente)
                    if (tRep != parentInterpolate[tRep]) {
                        tRep = parentInterpolate[tRep];
                    }
                    int repAbove = plateauRep[tRep];

                    // Se necessário, resolvemos para cima até achar um repOrig definido
                    while (repAbove != -1 && repOrig[repAbove] == -1 && tRep != parentInterpolate[tRep]) {
                        tRep = parentInterpolate[tRep];
                        repAbove = plateauRep[tRep];
                    }

                    qStar = (repAbove != -1 && repOrig[repAbove] != -1) ? repOrig[repAbove] : pStar; // fallback raro
                } else {
                    // nó não-representante no mesmo platô -> pai é o rep-ORIGINAL do mesmo platô
                    qStar = repOrig[repP];
                }
            } else {
                // nível diferente -> pai é o rep-ORIGINAL do platô de cima
                qStar = repOrig[repPar];
            }

            // (asserts opcionais para garantir existência)
            // assert(qStar != -1 && isOriginal1D(qStar, interpNumCols));

            int p = toOriginal1D(pStar, interpNumCols, numCols);
            int q = toOriginal1D(qStar,  interpNumCols, numCols);

            parent[p] = q;
            imgOrdered.push_back(p);
        }
        int roots = 0;
        for (int p : imgOrdered) if (parent[p] == p) ++roots;
        std::cout << "roots:"<<roots<<std::endl;
        assert(roots == 1);

        return std::make_tuple(parent, imgOrdered, numNodesInterpolateTree);

    }

};

#endif // BUILDER_MORPHOLOGICAL_TREE_BY_UNION_FIND_HPP
