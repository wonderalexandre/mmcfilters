#include "../include/BuilderTreeOfShapeByUnionFind.hpp"
#include "../include/ImageUtils.hpp"
#include <iostream>
    
    int BuilderTreeOfShapeByUnionFind::getInterpNumRows() {return this->interpNumRows;}
    int BuilderTreeOfShapeByUnionFind::getInterpNumCols() {return this->interpNumCols;}
    int* BuilderTreeOfShapeByUnionFind::getInterpolationMin() {return this->interpolationMin;}
    int* BuilderTreeOfShapeByUnionFind::getInterpolationMax() {return this->interpolationMax;}
    int* BuilderTreeOfShapeByUnionFind::getImgR() {return this->imgR;}
    int* BuilderTreeOfShapeByUnionFind::getImgU() {return this->imgU;}
    int* BuilderTreeOfShapeByUnionFind::getParent() {return this->parent;}


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
      * - N.Boutry, T.Géraud, L.Najman, "How to Make nD Functions Digitally Well-Composed in a Self-dual Way", ISMM 2015.
      * - N.Boutry, T.Géraud, L.Najman, "On Making {$n$D} Images Well-Composed by a Self-Dual Local Interpolation", DGCI 2014
      */
    void BuilderTreeOfShapeByUnionFind::interpolateImage(int* img, int numRows, int numCols) {
        constexpr int adjCircleCol[] = {-1, +1, -1, +1};
        constexpr int adjCircleRow[] = {-1, -1, +1, +1};

        constexpr int adjRetHorCol[] = {0, 0};
        constexpr int adjRetHorRow[] = {-1, +1};

        constexpr int adjRetVerCol[] = {+1, -1};
        constexpr int adjRetVerRow[] = {0, 0};

        this->interpNumCols = numCols * 2 + 1;
        this->interpNumRows = numRows * 2 + 1;

        // Aloca memória para os resultados de interpolação (mínimo e máximo)
        this->interpolationMin = new int[interpNumCols * interpNumRows];
        this->interpolationMax = new int[interpNumCols * interpNumRows];

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
                    } else {
                        continue;
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

    

    void BuilderTreeOfShapeByUnionFind::sort() {
        int size = this->interpNumCols * this->interpNumRows;
        bool* dejavu = new bool[size]();  // Vetor de booleanos, inicializado com false
        this->imgR = new int[size];        // Pixels ordenados
        this->imgU = new int[size];        // Níveis de cinza da imagem
        
        PriorityQueueToS queue;  // Fila de prioridade
        int pInfinito = 0;
        queue.initial(pInfinito, this->interpolationMin[pInfinito]);  
        dejavu[pInfinito] = true;

        this->adj = new AdjacencyRelation(interpNumRows, interpNumCols, 1);

        int i = 0;  // Contador para preencher imgR na ordem correta
        while (!queue.isEmpty()) {
            //queue.printCurrentPriority();
            int priorityQueue = queue.getCurrentPriority();
            int h = queue.priorityPop();  // Retirar o elemento com maior prioridade

            // Preencher imgU com o valor da prioridade corrente da fila
            imgU[h] = queue.getCurrentPriority();  // Prioridade corrente

            // Armazenar o índice h em imgR na ordem correta
            this->imgR[i] = h;
            
            // Adjacências
            for(int n: adj->getAdjPixels(h)){
                if (!dejavu[n]) {
                    queue.priorityPush(n, this->interpolationMin[n], this->interpolationMax[n]);
                    dejavu[n] = true;  // Marcar como processado
                }
            }
            i++;
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

            for(int n: adj->getAdjPixels(p)){
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
