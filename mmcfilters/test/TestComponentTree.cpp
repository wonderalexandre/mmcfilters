#include "../tests/Tests.hpp"
#include "../mmcfilters/include/AdjacencyRelation.hpp"

int main() {
    // Definição da imagem e parâmetros
    int numRows, numCols;
    int* img = getPassatImage(numRows, numCols);
    int n = numRows * numCols;
    double radioAdj = 1.5;

    // Criação das Component Trees
    MorphologicalTreePtr maxtree = std::make_shared<MorphologicalTree>(img, numRows, numCols, true, radioAdj);
    MorphologicalTreePtr mintree = std::make_shared<MorphologicalTree>(img, numRows, numCols, false, radioAdj);

    // Executar testes
    testComponentTree(mintree, "Min-Tree", mintree->reconstructionImage(), numRows, numCols);
    testComponentTree(maxtree, "Max-Tree", maxtree->reconstructionImage(), numRows, numCols);


    // Liberação de memória
    delete[] imgMaxtree;
    delete[] imgMintree;
    //delete maxtree;
    //delete mintree;
    delete[] img;

    return 0;
}