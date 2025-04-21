#include "../include/BuilderTreeOfShapeByUnionFind.hpp"
#include "../include/ImageUtils.hpp"
#include "../include/AdjacencyRelation.hpp"
#include <iostream>
    
    int BuilderTreeOfShapeByUnionFind::getInterpNumRows() {return this->interpNumRows;}
    int BuilderTreeOfShapeByUnionFind::getInterpNumCols() {return this->interpNumCols;}
    PixelType* BuilderTreeOfShapeByUnionFind::getInterpolationMin() {return this->interpolationMin;}
    PixelType* BuilderTreeOfShapeByUnionFind::getInterpolationMax() {return this->interpolationMax;}
    int* BuilderTreeOfShapeByUnionFind::getImgR() {return this->imgR;}
    PixelType* BuilderTreeOfShapeByUnionFind::getImgU() {return this->imgU;}
    int* BuilderTreeOfShapeByUnionFind::getParent() {return this->parent;}
    AdjacencyUC* BuilderTreeOfShapeByUnionFind::getAdjacency(){ return this->adj;}

    BuilderTreeOfShapeByUnionFind::BuilderTreeOfShapeByUnionFind(){
        
    }

    BuilderTreeOfShapeByUnionFind::~BuilderTreeOfShapeByUnionFind() {
        delete[] interpolationMin;
        delete[] interpolationMax;
        delete[] parent;
        delete[] imgR;
        delete[] imgU;
        delete adj;
    }

     /**
      * Implementation based on the paper: 
      *  - Thesi of the N.Boutry
      * - T. Géraud, E. Carlinet, and S. Crozet, Self-Duality and Digital Topology: Links Between the Morphological Tree of Shapes and Well-Composed Gray-Level Images, ISMM 2015
      * - N.Boutry, T.Géraud, L.Najman, "How to Make nD Functions Digitally Well-Composed in a Self-dual Way", ISMM 2015.
      * - N.Boutry, T.Géraud, L.Najman, "On Making {$n$D} Images Well-Composed by a Self-Dual Local Interpolation", DGCI 2014
      */
     void BuilderTreeOfShapeByUnionFind::interpolateImage(ImagePtr imgPtr) {
        PixelType* img = imgPtr->rawData();
        int numRows = imgPtr->numRows;
        int numCols = imgPtr->numCols;

        this->is4c8cConnectivity = false;
        constexpr int adjCircleCol[] = {-1, +1, -1, +1};
        constexpr int adjCircleRow[] = {-1, -1, +1, +1};

        constexpr int adjRetHorCol[] = {0, 0};
        constexpr int adjRetHorRow[] = {-1, +1};

        constexpr int adjRetVerCol[] = {+1, -1};
        constexpr int adjRetVerRow[] = {0, 0};

        this->interpNumCols = numCols * 2 + 1;
        this->interpNumRows = numRows * 2 + 1;

        // Aloca memória para os resultados de interpolação (mínimo e máximo)
        this->interpolationMin = new PixelType[interpNumCols * interpNumRows];
        this->interpolationMax = new PixelType[interpNumCols * interpNumRows];

        int numBoundary = 2 * (numRows + numCols) - 4;
        int* pixels = new int[numBoundary];  // Para calcular a mediana

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
        delete[] pixels;

        
        int qT, qCol, qRow, min, max;
        const int* adjCol = nullptr;
        const int* adjRow = nullptr;
        int adjSize;
        this->adj = new AdjacencyUC(interpNumRows, interpNumCols, false);

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

    void BuilderTreeOfShapeByUnionFind::interpolateImage4c8c(ImagePtr imgPtr) {
        PixelType* img = imgPtr->rawData();
        int numRows = imgPtr->numRows;
        int numCols = imgPtr->numCols;

        this->is4c8cConnectivity = true;
        this->interpNumCols = numCols * 2 + 1;
        this->interpNumRows = numRows * 2 + 1;
        this->adj = new AdjacencyUC(interpNumRows, interpNumCols, true);


        // Aloca memória para os resultados de interpolação (mínimo e máximo)
        this->interpolationMin = new PixelType[interpNumCols * interpNumRows];
        this->interpolationMax = new PixelType[interpNumCols * interpNumRows];
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

    void BuilderTreeOfShapeByUnionFind::sort() {
        int size = this->interpNumCols * this->interpNumRows;
        bool* dejavu = new bool[size]();  // Vetor de booleanos, inicializado com false
        this->imgR = new int[size];        // Pixels ordenados
        this->imgU = new PixelType[size];        // Níveis de cinza da imagem
        
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
            for(int n: adj->getNeighboringPixels(h)){
                if (!dejavu[n]) {
                    queue.priorityPush(n, this->interpolationMin[n], this->interpolationMax[n]);
                    dejavu[n] = true;  // Marcar como processado
                }
            }
            priorityQueueOld = priorityQueue;
        }
        delete[] dejavu;
    }

    int BuilderTreeOfShapeByUnionFind::findRoot(int zPar[], int p) {
        if (zPar[p] == p) {
            return p;
        } else {
            zPar[p] = findRoot(zPar, zPar[p]);
            return zPar[p];
        }
    }

    void BuilderTreeOfShapeByUnionFind::createTreeByUnionFind() {
        this->parent = new int[interpNumCols * interpNumRows];
        int* zPar = new int[interpNumCols * interpNumRows];
        const int NIL = -1;
        for (int p = 0; p < interpNumCols * interpNumRows; p++) {
            zPar[p] = NIL; // Assumindo que NIL é uma constante definida em outro lugar
        }
        for (int i = this->interpNumCols * this->interpNumRows - 1; i >= 0; i--) {
            int p = this->imgR[i];
            this->parent[p] = p;
            zPar[p] = p;

            for(int n: adj->getNeighboringPixels(p)){
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
        for (int i = 0; i < this->interpNumCols * this->interpNumRows; i++) {
            int p = this->imgR[i];
            int q = this->parent[p];
            if (this->imgU[parent[q]] == this->imgU[q]) { 
                this->parent[p] = this->parent[q];
            }
        }

        delete[] zPar; // Liberar memória de zPar

        
    }
