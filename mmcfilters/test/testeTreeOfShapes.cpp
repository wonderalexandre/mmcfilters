#include <iomanip>
#include <fstream>
#include <iostream>
#include <sstream>

#include "../include/BuilderTreeOfShapeByUnionFind.hpp"
#include "../include/NodeMT.hpp"
#include "../include/ImageUtils.hpp"
#include "../include/MorphologicalTree.hpp"
#include "Tests.hpp"


int main() {
   
    ImagePtr image = getSimpleImage();

    BuilderTreeOfShapeByUnionFind builder;

    // Receber os ponteiros de interpolação (mínimo e máximo)
    builder.interpolateImage4c8c(image);
    PixelType* interpolationMin = builder.getInterpolationMin();
    PixelType* interpolationMax = builder.getInterpolationMax();

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
    PixelType* imgU = builder.getImgU();

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
    MorphologicalTree tree(image, "4c8c");
    std::cout << "Depth:" << tree.getDepth() << ", |nodes|:" << tree.getNumNodes() << std::endl;
    printTree( tree.getRoot() );

     
    return 0;
}
