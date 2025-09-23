#include <iomanip>
#include <fstream>
#include <iostream>
#include <sstream>

#include "../mmcfilters/trees/BuilderMorphologicalTreeByUnionFind.hpp"
#include "../mmcfilters/utils/Common.hpp"
#include "Tests.hpp"

using namespace mmcfilters;

int main() {
   
    ImageUInt8Ptr image = getSimpleImage();

    BuilderTreeOfShape builder;

    // Receber os ponteiros de interpolação (mínimo e máximo)
    auto [interpolationMin, interpolationMax, _] = builder.interpolateImage4c8c(image);
    int interpNumRows = image->getNumRows() * 2 + 1;
    int interpNumCols = image->getNumCols() * 2 + 1;
    
        
    // Imprimir os resultados da interpolação
    std::cout << "\nInterpolação: " << interpNumRows << " x " << interpNumCols << std::endl;
    for (int r = 0; r < interpNumRows; ++r) {
        for (int c = 0; c < interpNumCols; ++c) {
            int index = ImageUtils::to1D(r, c, interpNumCols);
    
            std::ostringstream cell;
            if (r % 2 == 1 && c % 2 == 1) {
                cell << "  " << static_cast<int>(interpolationMax[index]) << "  ";
            } else {
                cell << "[" << static_cast<int>(interpolationMin[index]) << "," << static_cast<int>(interpolationMax[index]) << "]";
            }
    
            std::string cellStr = cell.str();
    
            // Garante que a célula tenha exatamente 8 caracteres (alinha à direita)
            std::cout << std::setw(10) << cellStr;
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    

    // Ordenar a interpolação mínima
    auto [imgU, imgR, adjUC] = builder.sort(image);

    std::cout << "Adjacency Relation" << std::endl;
    for(int index : adjUC.getNeighborPixels(3, 3)){
        auto [r, c] = ImageUtils::to2D(index, interpNumCols);
        std::cout << "(" << r << ", " << c << ") = " <<  "[" << static_cast<int>(interpolationMin[index]) << "," << static_cast<int>(interpolationMax[index]) << "]" << std::endl;
    }
    
    std::cout << "\nimgU: " << interpNumRows << " x " << interpNumCols << std::endl;
    // Imprimir os resultados da interpolação ordenada
    for (int row = 0; row < interpNumRows; ++row) {    
        for (int col = 0; col < interpNumCols; ++col) {
            int index = ImageUtils::to1D(row, col, interpNumCols);
            std::cout << std::setw(3) << static_cast<int>(imgU[index]) << ", ";
        }
        std::cout << std::endl;
    }
    std::cout << "\nimgR: " << interpNumRows << " x " << interpNumCols << std::endl;
    // Imprimir os resultados da interpolação ordenada
    for (int row = 0; row < interpNumRows; ++row) {    
        for (int col = 0; col < interpNumCols; ++col) {
            int index = ImageUtils::to1D(row, col, interpNumCols);
            std::cout << std::setw(4) << imgR[index] << ", ";
        }
        std::cout << std::endl;
    }

    auto [parent, orderedPixels, numNodes] = builder.createTreeByUnionFind(image);
    std::cout << "\nparent: " << interpNumRows << " x " << interpNumCols << std::endl;
    for (int row = 0; row < interpNumRows; ++row) {
        for (int col = 0; col < interpNumCols; ++col) {
            int index = ImageUtils::to1D(row, col, interpNumCols);
            std::cout << std::setw(3) << parent[index] << ", ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
    MorphologicalTree tree(image, "4c8c");
    std::cout << "|nodes|:" << tree.getNumNodes() << std::endl;
    printTree( tree.getRoot() );

     
    return 0;
}
