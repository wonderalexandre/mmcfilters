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
      * Implementation based on the paper: N.Boutry, T.Géraud, L.Najman, "How to Make nD Functions Digitally Well-Composed in a Self-dual Way", ISMM 2015.
      */
    void BuilderTreeOfShapeByUnionFind::interpolateImage(int* img, int num_rows, int num_cols) {
        constexpr int adjCircleX[] = {-1, +1, -1, +1};
        constexpr int adjCircleY[] = {-1, -1, +1, +1};

        constexpr int adjRetHorX[] = {0, 0};
        constexpr int adjRetHorY[] = {-1, +1};

        constexpr int adjRetVerX[] = {+1, -1};
        constexpr int adjRetVerY[] = {0, 0};

        this->interpNumCols = num_cols * 2 + 1;
        this->interpNumRows = num_rows * 2 + 1;

        // Aloca memória para os resultados de interpolação (mínimo e máximo)
        this->interpolationMin = new int[interpNumCols * interpNumRows];
        this->interpolationMax = new int[interpNumCols * interpNumRows];

        int numBoundary = 2 * (num_rows + num_cols) - 4;
        int* pixels = new int[numBoundary];  // Para calcular a mediana

        int x, y, pT, i = 0; // i é um contador para o array pixels
        
        for (int p = 0; p < num_cols * num_rows; p++) {
            auto [x, y] = ImageUtils::to2D(p, num_cols);

            // Verifica se o pixel está na borda
            if (x == 0 || x == num_cols - 1 || y == 0 || y == num_rows - 1) {
                pixels[i++] = img[p]; // Adiciona o pixel ao array pixels
            }

            // Calcula o índice para imagem interpolada
            pT = ImageUtils::to1D(2 * x + 1, 2 * y + 1, this->interpNumCols);

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

        
        int qT, qX, qY, min, max;
        const int* adjX = nullptr;
        const int* adjY = nullptr;
        int adjSize;

        for (y=0; y < this->interpNumRows; y++){
            for (x=0; x < this->interpNumCols; x++){
                if (x % 2 == 1 && y % 2 == 1) continue;
                pT = ImageUtils::to1D(x, y, this->interpNumCols);
                if(x == 0 || x == this->interpNumCols - 1 || y == 0 || y == this->interpNumRows - 1){
                    max = median;
                    min = median;
                }else{
                    if (x % 2 == 0 && y % 2 == 0) { 
                        adjX = adjCircleX;
                        adjY = adjCircleY;
                        adjSize = 4;
                    } else if (x % 2 == 0 && y % 2 == 1) {
                        adjX = adjRetVerX;
                        adjY = adjRetVerY;
                        adjSize = 2;
                    } else if (x % 2 == 1 && y % 2 == 0) {
                        adjX = adjRetHorX;
                        adjY = adjRetHorY;
                        adjSize = 2;
                    } else {
                        continue;
                    }

                    min = INT_MAX;
                    max = INT_MIN;
                    for (int i = 0; i < adjSize; i++) {
                        qY = y + adjY[i];
                        qX = x + adjX[i];

                        if (qY >= 0 && qX >= 0 && qY < this->interpNumRows && qX < this->interpNumCols) {
                            qT = ImageUtils::to1D(qX, qY, this->interpNumCols);

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
