#ifndef BUILDER_MORPHOLOGICAL_TREE_BY_UNION_FIND_HPP
#define BUILDER_MORPHOLOGICAL_TREE_BY_UNION_FIND_HPP

#include "../include/ComponentTree.hpp"
#include "../include/AdjacencyRelation.hpp"
#include "../include/Common.hpp"
#include "../include/AdjacencyUC.hpp"

// Builder externo para construir ComponentTree via Union-Find (Pixels-only).
class BuilderComponentTreeByUnionFind {
private:
    ComponentTree* tree;

public:
    explicit BuilderComponentTreeByUnionFind(ComponentTree* t) : tree(t) {}

    // Ordenação estável dos pixels por nível de cinza
    std::vector<int> countingSort(ImageUInt8Ptr imgPtr) {
        int n = tree->getNumRowsOfImage() * tree->getNumColsOfImage();
        auto img = imgPtr->rawData();
        int maxvalue = img[0];
        for (int i = 1; i < n; i++) if (maxvalue < img[i]) maxvalue = img[i];

        std::vector<uint32_t> counter(maxvalue + 1, 0);
        std::vector<int> orderedPixels(n);

        if (tree->isMaxtree()) {
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
    void createTreeByUnionFind(ImageUInt8Ptr imgPtr) {
        std::vector<int> orderedPixels = countingSort(imgPtr);
        auto& pixelToNodeId = tree->pixelToNodeId;
        auto* adj = tree->adj.get();

        int numPixels = tree->getNumRowsOfImage() * tree->getNumColsOfImage();
        std::vector<int> zPar(numPixels, -1);
        std::vector<int> parent(numPixels, -1);
        auto findRoot = [&](int p) {
            while (zPar[p] != p) { zPar[p] = zPar[zPar[p]]; p = zPar[p]; }
            return p;
        };
        auto img = imgPtr->rawData();

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

        tree->reserveNodes(numNodes);
        tree->pixelBuffer = std::make_shared<PixelSetManager>(numPixels, numNodes);
        tree->pixelView = tree->pixelBuffer->view();
        auto& pixelView = tree->pixelView;
        int indice = 0;
        for (int i = 0; i < numPixels; i++) {
            int p = orderedPixels[i];

            //Construção da árvore e arena
            if (p == parent[p]) {
                int threshold1 = tree->maxtreeTreeType ? 0 : 255;
                int threshold2 = img[p];
                pixelToNodeId[p] = tree->root = tree->makeNode(p, -1, threshold1, threshold2);
            } else if (img[p] != img[parent[p]]) {
                int threshold1 = tree->maxtreeTreeType ? img[parent[p]] + 1 : img[parent[p]] - 1;
                int threshold2 = img[p];
                pixelToNodeId[p] = tree->makeNode(p, pixelToNodeId[parent[p]], threshold1, threshold2);
            } else {
                pixelToNodeId[p] = pixelToNodeId[parent[p]];
            }

            //Construção de PixelSetManager
            if (p == parent[p] || img[p] != img[parent[p]]) {
                pixelView.indexToPixel[indice] = p;
                pixelView.pixelToIndex[p] = indice;
                pixelView.sizeSets[indice] = 1;
                pixelView.pixelsNext[p] = p;
                indice++;
            } else {
                pixelView.pixelsNext[p] = pixelView.pixelsNext[parent[p]];
                pixelView.pixelsNext[parent[p]] = p;
                int idx = pixelView.pixelToIndex[parent[p]];
                pixelView.sizeSets[idx]++;
            }
        }

        assert((indice == numNodes) && "Erro na contagem de sets");
    }
};






class BuilderTreeOfShapeByUnionFind {
private:
    int interpNumRows;
    int interpNumCols;
    std::unique_ptr<uint8_t[]> interpolationMin;
    std::unique_ptr<uint8_t[]> interpolationMax;
    std::unique_ptr<uint8_t[]> imgU;
    std::unique_ptr<int[]> parent;
    std::unique_ptr<int[]> imgR; 
    
    std::unique_ptr<AdjacencyUC> adj;
    bool is4c8cConnectivity;

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


    
    int getInterpNumRows() {return this->interpNumRows;}
    int getInterpNumCols() {return this->interpNumCols;}
    uint8_t* getImgU() {return this->imgU.get();}
    int* getParent() {return this->parent.get();}
    int* getImgR() {return this->imgR.get();}
    AdjacencyUC* getAdjacency() { return adj.get(); }
    uint8_t* getInterpolationMin() { return interpolationMin.get(); }
    uint8_t* getInterpolationMax() { return interpolationMax.get(); }

    BuilderTreeOfShapeByUnionFind(){ }
    ~BuilderTreeOfShapeByUnionFind() { }

     /**
      * Implementation based on the paper: 
      *  - Thesi of the N.Boutry
      * - T. Géraud, E. Carlinet, and S. Crozet, Self-Duality and Digital Topology: Links Between the Morphological Tree of Shapes and Well-Composed Gray-Level Images, ISMM 2015
      * - N.Boutry, T.Géraud, L.Najman, "How to Make nD Functions Digitally Well-Composed in a Self-dual Way", ISMM 2015.
      * - N.Boutry, T.Géraud, L.Najman, "On Making {$n$D} Images Well-Composed by a Self-Dual Local Interpolation", DGCI 2014
      */
     void interpolateImage(ImageUInt8Ptr imgPtr) {
        auto img = imgPtr->rawData();
        int numRows = imgPtr->getNumRows();
        int numCols = imgPtr->getNumCols();

        this->is4c8cConnectivity = false;
        constexpr int adjCircleCol[] = {-1, +1, -1, +1};
        constexpr int adjCircleRow[] = {-1, -1, +1, +1};

        constexpr int adjRetHorCol[] = {0, 0};
        constexpr int adjRetHorRow[] = {-1, +1};

        constexpr int adjRetVerCol[] = {+1, -1};
        constexpr int adjRetVerRow[] = {0, 0};

        this->interpNumCols = numCols * 2 + 1;
        this->interpNumRows = numRows * 2 + 1;
        int size = interpNumCols * interpNumRows;

        // Aloca memória para os resultados de interpolação (mínimo e máximo)
        this->interpolationMin = std::make_unique<uint8_t[]>(size);
        this->interpolationMax = std::make_unique<uint8_t[]>(size);

        int numBoundary = 2 * (numRows + numCols) - 4;
        
        std::unique_ptr<uint8_t[]> pixelsPtr(new uint8_t[numBoundary]);// Para calcular a mediana
        uint8_t* pixels = pixelsPtr.get();

        int pT, i = 0; // i é um contador para o array pixels
        
        for (int p = 0; p < numCols * numRows; p++) {
            auto [row, col] = ImageUtils::to2D(p, numCols);

            // Verifica se o pixel está na borda
            if (row == 0 || row == numRows - 1 || col == 0 || col == numCols - 1) {
                pixels[i++] = img[p]; // Adiciona o pixel ao array pixels
            }

            // Calcula o índice para imagem interpolada
            pT = ImageUtils::to1D(2 * row + 1, 2 * col + 1, this->interpNumCols);

            // Define os valores de interpolação
            this->interpolationMin[pT] = this->interpolationMax[pT] = img[p];
        }

        std::sort(pixels, pixels + numBoundary);
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
        this->adj = std::make_unique<AdjacencyUC>(interpNumRows, interpNumCols, false);

        for (int row=0; row < this->interpNumRows; row++){
            for (int col=0; col < this->interpNumCols; col++){
                if (col % 2 == 1 && row % 2 == 1) continue;
                pT = ImageUtils::to1D(row, col, this->interpNumCols);
                if(col == 0 || col == this->interpNumCols - 1 || row == 0 || row == this->interpNumRows - 1){
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

                        if (qRow >= 0 && qCol >= 0 && qRow < this->interpNumRows && qCol < this->interpNumCols) {
                            qT = ImageUtils::to1D(qRow, qCol, this->interpNumCols);

                            if (interpolationMax[qT] > max) {
                                max = this->interpolationMax[qT];
                            }
                            if (interpolationMin[qT] < min) {
                                min = this->interpolationMin[qT];
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
                this->interpolationMin[pT] = min;
                this->interpolationMax[pT] = max;
            }
        }
       
    }

    void interpolateImage4c8c(ImageUInt8Ptr imgPtr) {
        auto img = imgPtr->rawData();
        int numRows = imgPtr->getNumRows();
        int numCols = imgPtr->getNumCols();

        this->is4c8cConnectivity = true;
        this->interpNumCols = numCols * 2 + 1;
        this->interpNumRows = numRows * 2 + 1;
        int size = interpNumCols * interpNumRows;
        this->adj = std::make_unique<AdjacencyUC>(interpNumRows, interpNumCols, true);


        // Aloca memória para os resultados de interpolação (mínimo e máximo)
        this->interpolationMin = std::make_unique<uint8_t[]>(size);
        this->interpolationMax = std::make_unique<uint8_t[]>(size);
        int pT, i = 0; // i é um contador para o array pixels
        
         // Compute interval from 2-faces.
        for (int p = 0; p < numCols * numRows; p++) {
            auto [row, col] = ImageUtils::to2D(p, numCols);

            // Calcula o índice para imagem interpolada
            pT = ImageUtils::to1D(2 * row + 1, 2 * col + 1, this->interpNumCols);

            // Define os valores de interpolação
            this->interpolationMin[pT] = this->interpolationMax[pT] = img[p];
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
        for (int row=0; row < this->interpNumRows; row++){
            int col;
            if(row % 2 == 1){ //horizontal e vertical
                col = 0;
                int v1 = getValue(row, col+1);
                this->interpolationMin[ImageUtils::to1D(row, col, this->interpNumCols)] = v1;
                this->interpolationMax[ImageUtils::to1D(row, col, this->interpNumCols)] = v1;

                col = this->interpNumCols - 1;
                v1 = getValue(row, col -1);
                this->interpolationMin[ImageUtils::to1D(row, col, this->interpNumCols)] = v1;
                this->interpolationMax[ImageUtils::to1D(row, col, this->interpNumCols)] = v1;
            }else{ //circulos
                if(row == 0){
                    col = 0;
                    int v1 = getValue(row+1, col+1);
                    this->interpolationMin[ImageUtils::to1D(row, col, this->interpNumCols)] = v1;
                    this->interpolationMax[ImageUtils::to1D(row, col, this->interpNumCols)] = v1;

                    col = this->interpNumCols - 1;
                    v1 = getValue(row+1, col -1);
                    this->interpolationMin[ImageUtils::to1D(row, col, this->interpNumCols)] = v1;
                    this->interpolationMax[ImageUtils::to1D(row, col, this->interpNumCols)] = v1;

                }else if(row == this->interpNumRows-1){
                    col = 0;
                    int v1 = getValue(row-1, 1);
                    this->interpolationMin[ImageUtils::to1D(row, col, this->interpNumCols)] = v1;
                    this->interpolationMax[ImageUtils::to1D(row, col, this->interpNumCols)] = v1;

                    col = this->interpNumCols - 1;
                    v1 = getValue(row-1, col - 1);
                    this->interpolationMin[ImageUtils::to1D(row, col, this->interpNumCols)] = v1;
                    this->interpolationMax[ImageUtils::to1D(row, col, this->interpNumCols)] = v1;
                }else{
                    col = 0;
                    int v1 = getValue(row-1, col+1);
                    int v2 = getValue(row+1, col+1);
                    this->interpolationMin[ImageUtils::to1D(row, 0, this->interpNumCols)] = std::min(v1, v2);
                    this->interpolationMax[ImageUtils::to1D(row, 0, this->interpNumCols)] = std::max(v1, v2);

                    col = this->interpNumCols - 1;
                    v1 = getValue(row-1, col-1);
                    v2 = getValue(row+1, col-1);
                    this->interpolationMin[ImageUtils::to1D(row, col, this->interpNumCols)] = std::min(v1, v2);
                    this->interpolationMax[ImageUtils::to1D(row, col, this->interpNumCols)] = std::max(v1, v2);
                }
            }
        }
        
        for (int col=1; col < this->interpNumCols-1; col++){
            int row;
            if(col % 2 == 1){ //horizontal e vertical
                row = 0;
                int v1 = getValue(row+1, col);
                this->interpolationMin[ImageUtils::to1D(row, col, this->interpNumCols)] = v1;
                this->interpolationMax[ImageUtils::to1D(row, col, this->interpNumCols)] = v1;

                row = this->interpNumRows - 1;
                v1 = getValue(row-1, col);
                this->interpolationMin[ImageUtils::to1D(row, col, this->interpNumCols)] = v1;
                this->interpolationMax[ImageUtils::to1D(row, col, this->interpNumCols)] = v1;
            }else{ //circulos
                row = 0;
                int v1 = getValue(row+1, col-1);
                int v2 = getValue(row+1, col+1);
                this->interpolationMin[ImageUtils::to1D(row, col, this->interpNumCols)] = std::min(v1, v2);
                this->interpolationMax[ImageUtils::to1D(row, col, this->interpNumCols)] = std::max(v1, v2);

                row = this->interpNumRows - 1;
                v1 = getValue(row-1, col-1);
                v2 = getValue(row-1, col+1);
                this->interpolationMin[ImageUtils::to1D(row, col, this->interpNumCols)] = std::min(v1, v2);
                this->interpolationMax[ImageUtils::to1D(row, col, this->interpNumCols)] = std::max(v1, v2);
            }
        }

        // Compute interval from 1-faces 
        for (int row=1; row < this->interpNumRows-1; row++){
            for (int col=1; col < this->interpNumCols-1; col++){
                if (row % 2 == 1 && col % 2 == 1) continue;  // já definido

                pT = ImageUtils::to1D(row, col, this->interpNumCols);
                if (col % 2 == 0 && row % 2 == 1) {
                    int v1 = getValue(row, col+1);
                    int v2 = getValue(row, col-1);
                    this->interpolationMin[pT] = std::min(v1, v2);
                    this->interpolationMax[pT] = std::max(v1, v2);
                } else if (col % 2 == 1 && row % 2 == 0) {
                    int v1 = getValue(row+1, col);
                    int v2 = getValue(row-1, col);
                    this->interpolationMin[pT] = std::min(v1, v2);
                    this->interpolationMax[pT] = std::max(v1, v2);
                } 
            }
        }
         // Compute interval from 0-faces 
         for (int row=1; row < this->interpNumRows-1; row++){
            for (int col=1; col < this->interpNumCols-1; col++){
                if (row % 2 == 1 && col % 2 == 1) continue;  // já definido
                pT = ImageUtils::to1D(row, col, this->interpNumCols);
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
                        this->adj->setDiagonalConnection(row, col-1, DiagonalConnection::SE);
                        this->adj->setDiagonalConnection(row+1, col, DiagonalConnection::NW);
                        
                        this->adj->setDiagonalConnection(row - 1, col - 1, DiagonalConnection::SE);
                        this->adj->setDiagonalConnection(row, col, DiagonalConnection::SE | DiagonalConnection::NW);
                        this->adj->setDiagonalConnection(row + 1, col + 1, DiagonalConnection::NW);

                        this->adj->setDiagonalConnection(row-1, col, DiagonalConnection::SE);
                        this->adj->setDiagonalConnection(row, col+1, DiagonalConnection::NW);

                        this->interpolationMin[pT] = min_v0v3;
                        this->interpolationMax[pT] = max_v0v3;
                    }
                    else if (max_v0v3 > min_v1v2) {
                        // Saddle point configuration 2
                        this->adj->setDiagonalConnection(row, col-1, DiagonalConnection::NE);
                        this->adj->setDiagonalConnection(row-1, col, DiagonalConnection::SW);

                        this->adj->setDiagonalConnection(row-1, col+1, DiagonalConnection::SW);
                        this->adj->setDiagonalConnection(row, col, DiagonalConnection::SW | DiagonalConnection::NE);
                        this->adj->setDiagonalConnection(row + 1, col - 1, DiagonalConnection::NE);

                        this->adj->setDiagonalConnection(row+1, col, DiagonalConnection::NE);
                        this->adj->setDiagonalConnection(row, col+1, DiagonalConnection::SW);

                        this->interpolationMin[pT] = min_v1v2;
                        this->interpolationMax[pT] = max_v1v2;
                    }else{
                        // Non-critical configuration.
                        this->interpolationMin[pT] = std::min(min_v0v3, min_v1v2);
                        this->interpolationMax[pT] = std::min(max_v0v3, max_v1v2);
                    }
                }

            }
        }
       
    }

    void sort() {
        int size = this->interpNumCols * this->interpNumRows;
        std::unique_ptr<bool[]> dejavu(new bool[size]());  // Vetor de booleanos, inicializado com false
        this->imgR = std::make_unique<int[]>(size);  // Pixels ordenados
        this->imgU = std::make_unique<uint8_t[]>(size);        // Níveis de cinza da imagem
        
        PriorityQueueToS queue;  // Fila de prioridade
        int pInfinito = ImageUtils::to1D(0, 0, interpNumCols);
        int priorityQueueOld = this->interpolationMin[pInfinito];
        queue.initial(pInfinito, priorityQueueOld);  
        dejavu[pInfinito] = true;

        int order = 0; 
        int depth = 0;
        while (!queue.isEmpty()) {
            int h = queue.priorityPop();  // Retirar o elemento com maior prioridade
            int priorityQueue = queue.getCurrentPriority(); // Prioridade corrente
            if(this->is4c8cConnectivity){
                if(priorityQueue != priorityQueueOld) depth++;
                imgU[h] = depth;
            }else{
                imgU[h] = priorityQueue;
            }
            
            // Armazenar o índice h em imgR na ordem correta
            this->imgR[order++] = h;
            
            // Adjacências
            for(int n: adj->getNeighborPixels(h)){
                if (!dejavu[n]) {
                    queue.priorityPush(n, this->interpolationMin[n], this->interpolationMax[n]);
                    dejavu[n] = true;  // Marcar como processado
                }
            }
            priorityQueueOld = priorityQueue;
        }
        //delete[] dejavu;
    }

    int findRoot(int* zPar, int p) {
        if (zPar[p] == p) {
            return p;
        } else {
            zPar[p] = findRoot(zPar, zPar[p]);
            return zPar[p];
        }
    }

    void createTreeByUnionFind() {
        int size = this->interpNumCols * this->interpNumRows;
        this->parent =  std::make_unique<int[]>(size);

        //this->parent = new int[interpNumCols * interpNumRows];
        std::unique_ptr<int[]> zParPtr(new int[size]);
        int* zPar = zParPtr.get(); // Pega o ponteiro do std::unique_ptr
        const int NIL = -1;
        for (int p = 0; p < size; p++) {
            zPar[p] = NIL; // Assumindo que NIL é uma constante definida em outro lugar
        }
        for (int i = size - 1; i >= 0; i--) {
            int p = this->imgR[i];
            this->parent[p] = p;
            zPar[p] = p;

            for(int n: adj->getNeighborPixels(p)){
                if (zPar[n] != NIL) {
                    int r = findRoot(zPar, n);
                    if (p != r) {
                        this->parent[r] = p;
                        zPar[r] = p;
                    }
                }
            }
        }

        // Canonização da árvore
        for (int i = 0; i < size; i++) {
            int p = this->imgR[i];
            int q = this->parent[p];
            if (this->imgU[parent[q]] == this->imgU[q]) { 
                this->parent[p] = this->parent[q];
            }
        }

        //delete[] zPar; // Liberar memória de zPar

        
    }


};

#endif // BUILDER_MORPHOLOGICAL_TREE_BY_UNION_FIND_HPP
