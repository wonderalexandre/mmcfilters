#include "Tests.hpp"
#include "../include/AdjacencyRelation.hpp"
#include "../include/Common.hpp"

int main() {
    // Definição da imagem e parâmetros
    ImageUInt8Ptr img = getPassatImage();
    
    double radioAdj = 1.5;

    // Criação das Component Trees
    MorphologicalTreePtr maxtree = std::make_shared<MorphologicalTree>(img, true, radioAdj);
    MorphologicalTreePtr mintree = std::make_shared<MorphologicalTree>(img, false, radioAdj);

    // Executar testes
    testComponentTree(mintree, "Min-Tree", mintree->reconstructionImage());
    testComponentTree(maxtree, "Max-Tree", maxtree->reconstructionImage());


    
    return 0;
}