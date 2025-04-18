#include <iomanip>
#include <fstream>
#include <iostream>
#include <sstream>

#include "../include/BuilderTreeOfShapeByUnionFind.hpp"
#include "../include/NodeMT.hpp"
#include "../include/ImageUtils.hpp"
#include "../include/MorphologicalTree.hpp"
#include "Tests.hpp"

/*
void printTree(NodeMTPtr root, int indent = 0) {
    
    // Imprime o nó atual com indentação
    for (int i = 0; i < indent; ++i) {
        std::cout << "|-";
    }
    std::cout << "Node: " << root->getIndex() <<  ", Level: " << root->getLevel()<< ", cnps: " << root->getCNPs().size() << ", |children|:" << root->getChildren().size() << ", Area:" << root->getAreaCC() << std::endl;

    // Chama recursivamente a função para cada filho
    for (NodeMTPtr child : root->getChildren()) {
        printTree(child, indent + 1);
    }
}
*/

int main() {
    // Exemplo de teste com imagem representada como ponteiro 1D
    /*
    int image[] = {
        2, 2, 0, 0,
        2, 1, 1, 0,
        2, 2, 0, 0
    };
    int num_rows = 3;
    int num_cols = 4;
    

    int image[] = {
        5, 5, 5, 5, 5, 5, 5, 9, 9,
        5, 5, 0, 0, 5, 5, 9, 5, 5,
        5, 0, 5, 0, 5, 5, 9, 9, 5,
        5, 0, 0, 5, 5, 5, 9, 9, 5,
        5, 0, 0, 5, 5, 5, 5, 5, 5
    };
    int num_rows = 5;
    int num_cols = 9;

    
   int image[] = {
        3, 3, 3, 3,
        3, 0, 8, 3,
        3, 9, 4, 3,
        3, 3, 3, 3,
    };
    int num_rows = 4;
    int num_cols = 4;
    
    int image[] = {
        9, 11, 15,
        7, 1,  13,
        3, 5,  3
    };
    int num_rows = 3;
    int num_cols = 3;
    */
    /*
    int image[] = {
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1,
    1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1,
    1, 2, 2, 1, 1, 2, 2, 1, 1, 0, 0, 1, 1, 0, 0, 1,
    1, 2, 2, 1, 1, 2, 2, 1, 1, 0, 0, 1, 1, 0, 0, 1,
    1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1,
    1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
    };
    int num_rows = 8;
    int num_cols = 16;
    */
   int num_rows, num_cols;
    int* image = getSimpleImage(num_rows, num_cols);
/*
    int image[] = {
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 1, 1, 1, 7, 7, 7, 4, 4,
        4, 4, 1, 1, 1, 1, 7, 7, 7, 4, 4,
        4, 4, 1, 1, 1, 1, 7, 7, 7, 4, 4,
        4, 4, 1, 1, 4, 4, 4, 7, 7, 4, 4,
        4, 4, 1, 1, 1, 1, 7, 7, 7, 4, 4,
        4, 4, 1, 1, 1, 1, 7, 7, 7, 4, 4,
        4, 4, 4, 1, 1, 1, 7, 7, 7, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
    };
    int num_cols = 11;
    int num_rows = 10;

    
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
    builder.interpolateImage4c8c(image, num_rows, num_cols);
    int* interpolationMin = builder.getInterpolationMin();
    int* interpolationMax = builder.getInterpolationMax();

    // Imprimir os resultados da interpolação
    std::cout << "\nInterpolação: " << builder.getInterpNumRows() << " x " << builder.getInterpNumCols() << std::endl;
    for (int r = 0; r < builder.getInterpNumRows(); ++r) {
        for (int c = 0; c < builder.getInterpNumCols(); ++c) {
            int index = ImageUtils::to1D(r, c, builder.getInterpNumCols());
    
            std::ostringstream cell;
            if (r % 2 == 1 && c % 2 == 1) {
                cell << "  " << interpolationMax[index] << "  ";
            } else {
                cell << "[" << interpolationMin[index] << "," << interpolationMax[index] << "]";
            }
    
            std::string cellStr = cell.str();
    
            // Garante que a célula tenha exatamente 8 caracteres (alinha à direita)
            std::cout << std::setw(8) << cellStr;
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    
    AdjacencyUC* adj = builder.getAdjacency();
    for(int index : adj->getNeighboringPixels(3, 3)){
        auto [r, c] = ImageUtils::to2D(index, builder.getInterpNumCols());
        std::cout << "(" << r << ", " << c << ") = " <<  "[" << interpolationMin[index] << "," << interpolationMax[index] << "]" << std::endl;
    }


    
    // Ordenar a interpolação mínima
    builder.sort();
    int* imgR = builder.getImgR();
    int* imgU = builder.getImgU();

    std::cout << "\nimgU: " << builder.getInterpNumRows() << " x " << builder.getInterpNumCols() << std::endl;
    // Imprimir os resultados da interpolação ordenada
    for (int row = 0; row < builder.getInterpNumRows(); ++row) {    
        for (int col = 0; col < builder.getInterpNumCols(); ++col) {
            int index = ImageUtils::to1D(row, col, builder.getInterpNumCols());
            std::cout << std::setw(2) << imgU[index] << ", ";
        }
        std::cout << std::endl;
    }
    std::cout << "\nimgR: " << builder.getInterpNumRows() << " x " << builder.getInterpNumCols() << std::endl;
    int index = 0;
    int* imgRInv = new int[builder.getInterpNumRows() * builder.getInterpNumCols()];
    for(int i=0 ; i < builder.getInterpNumRows() * builder.getInterpNumCols(); i++){
        imgRInv[imgR[i]] = i;
    }

    for (int row = 0; row < builder.getInterpNumRows(); ++row) {
        for (int col = 0; col < builder.getInterpNumCols(); ++col) {
            int index = ImageUtils::to1D(row, col, builder.getInterpNumCols());

            std::cout << std::setw(3) << imgRInv[index] << ", ";
        }
        std::cout << std::endl;
    }
    delete[] imgRInv;

    builder.createTreeByUnionFind();
    int* parent = builder.getParent();
    std::cout << "\nparent: " << builder.getInterpNumRows() << " x " << builder.getInterpNumCols() << std::endl;
    for (int row = 0; row < builder.getInterpNumRows(); ++row) {
        for (int col = 0; col < builder.getInterpNumCols(); ++col) {
            int index = ImageUtils::to1D(row, col, builder.getInterpNumCols());
            std::cout << std::setw(3) << parent[index] << ", ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
    MorphologicalTree tree(image, num_rows, num_cols, "4c8c");
    std::cout << "Depth:" << tree.getDepth() << ", |nodes|:" << tree.getNumNodes() << std::endl;
    printTree( tree.getRoot() );

     
    return 0;
}
