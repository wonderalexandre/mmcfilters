#include <iostream>
#include "../include/BuilderTreeOfShapeByUnionFind.hpp"
#include "../include/NodeMT.hpp"
#include "../include/ImageUtils.hpp"
#include "../include/MorphologicalTree.hpp"
#include "../include/AdjacencyUC.hpp"
#include "../include/ImageUtils.hpp"

int main() {

    int numRows = 100;
    int numCols = 100;
    AdjacencyUC adj(numRows, numCols, true);

    int pRow = 2;
    int pCol = 2;
    adj.setDiagonalConnection(pRow, pCol, DiagonalConnection::NE | DiagonalConnection::SW);
    for (int q : adj.getNeighboringPixels(0,0)) {
        auto [qRow, qCol] = ImageUtils::to2D(q, numCols);
        std::cout << "Vizinho: (" << qRow << ", " << qCol << ")"<< std::endl;
    }

}