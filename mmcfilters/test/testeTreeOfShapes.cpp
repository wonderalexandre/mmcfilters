#include <iostream>
#include "../include/BuilderTreeOfShapeByUnionFind.hpp"
#include "../include/NodeCT.hpp"
#include "../include/ImageUtils.hpp"
#include "../include/ComponentTree.hpp"

void printTree(NodeCT* root, int indent = 0) {
    
    // Imprime o nó atual com indentação
    for (int i = 0; i < indent; ++i) {
        std::cout << "|-";
    }
    std::cout << "Node: " << root->getIndex() <<  ", Level: " << root->getLevel()<< ", cnps: " << root->getCNPs().size() << std::endl;

    // Chama recursivamente a função para cada filho
    for (NodeCT* child : root->getChildren()) {
        printTree(child, indent + 1);
    }
}


int main() {
    // Exemplo de teste com imagem representada como ponteiro 1D
    
    int image[] = {
        2, 2, 0, 0,
        2, 1, 1, 0,
        2, 2, 0, 0
    };
    int num_rows = 3;
    int num_cols = 4;
    /*
   int image[] = {
        1, 3,
        4, 1
    };
    int num_rows = 2;
    int num_cols = 2;
    
    int image[] = {
        9, 11, 15,
        7, 1,  13,
        3, 5,  3
    };
    int num_rows = 3;
    int num_cols = 3;
   
    
   int image[] = {
    128, 124, 150, 137, 106,
    116, 128, 156, 165, 117,
    117, 90,  131, 108, 151,
    107, 87,  118, 109, 167,
    107, 73,  125, 157, 117
    };
    int num_rows = 5;
    int num_cols = 5;
 */

    BuilderTreeOfShapeByUnionFind builder;

    // Receber os ponteiros de interpolação (mínimo e máximo)
    builder.interpolateImage(image, num_rows, num_cols);
    int* interpolationMin = builder.getInterpolationMin();
    int* interpolationMax = builder.getInterpolationMax();

    // Imprimir os resultados da interpolação
    std::cout << "Interpolação: " << builder.getInterpNumRows() << " x " << builder.getInterpNumCols() << std::endl;
    for (int y = 0; y < builder.getInterpNumRows(); ++y) {
        for (int x = 0; x < builder.getInterpNumCols(); ++x) {
            int index = ImageUtils::to1D(x, y, builder.getInterpNumCols());
            std::cout << "[" << interpolationMin[index] << ", " << interpolationMax[index] << "] ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    
    
    // Ordenar a interpolação mínima
    builder.sort();
    int* imgR = builder.getImgR();
    int* imgU = builder.getImgU();

    std::cout << "imgU: " << builder.getInterpNumRows() << " x " << builder.getInterpNumCols() << std::endl;
    // Imprimir os resultados da interpolação ordenada
    for (int y = 0; y < builder.getInterpNumRows(); ++y) {    
        for (int x = 0; x < builder.getInterpNumCols(); ++x) {
            int index = ImageUtils::to1D(x, y, builder.getInterpNumCols());
            std::cout  << imgU[index] << ", ";
        }
        std::cout << std::endl;
    }
    std::cout << "imgR: " << builder.getInterpNumRows() << " x " << builder.getInterpNumCols() << std::endl;
    for (int y = 0; y < builder.getInterpNumRows(); ++y) {
        for (int x = 0; x < builder.getInterpNumCols(); ++x) {
            int index = ImageUtils::to1D(x, y, builder.getInterpNumCols());
            std::cout << imgR[index] << ", ";
        }
        std::cout << std::endl;
    }

    builder.createTreeByUnionFind();
    int* parent = builder.getParent();
    std::cout << "parent: " << builder.getInterpNumRows() << " x " << builder.getInterpNumCols() << std::endl;
    for (int y = 0; y < builder.getInterpNumRows(); ++y) {
        for (int x = 0; x < builder.getInterpNumCols(); ++x) {
            int index = ImageUtils::to1D(x, y, builder.getInterpNumCols());
            std::cout << parent[index] << ", ";
        }
        std::cout << std::endl;
    }


    printTree( ComponentTree(image, num_rows, num_cols).getRoot() );

     
    return 0;
}
