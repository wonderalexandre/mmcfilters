#include "Tests.hpp"
#include "../include/AdjacencyRelation.hpp"
#include "../include/Common.hpp"

int main() {
    // Definição da imagem e parâmetros
    ImageUInt8Ptr img = getSimpleImage();
    
    double radioAdj = 1.5;

    // Criação das Component Trees
    MorphologicalTreePtr maxtree = std::make_shared<MorphologicalTree>(img, true, radioAdj);
    MorphologicalTreePtr mintree = std::make_shared<MorphologicalTree>(img, false, radioAdj);

    
    printTree(mintree->getRoot());
    testComponentTree(mintree, "Min-Tree", mintree->reconstructionImage());

    printTree(maxtree->getRoot());
    testComponentTree(maxtree, "Max-Tree", maxtree->reconstructionImage());


    
    return 0;
}