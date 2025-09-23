#include <iostream>
#include "../mmcfilters/utils/AdjacencyRelation.hpp"
#include "../mmcfilters/trees/BuilderMorphologicalTreeByUnionFind.hpp"
#include "../mmcfilters/utils/Common.hpp"

using namespace mmcfilters;

int main() {

    int numRows = 100;
    int numCols = 100;
    AdjacencyRelation adj4(numRows, numCols, 1);
    AdjacencyRelation adj8(numRows, numCols, 1.5);
    std::cout << "Adjacency 4-connectivity:" << adj4.getSize() << std::endl;
    std::cout << "Adjacency 8-connectivity:" << adj8.getSize() << std::endl;

    AdjacencyUC adj(numRows, numCols, true);

    int pRow = 2;
    int pCol = 2;
    adj.setDiagonalConnection(pRow, pCol, DiagonalConnection::NE | DiagonalConnection::SW);
    for (int q : adj.getNeighborPixels(0,0)) {
        auto [qRow, qCol] = ImageUtils::to2D(q, numCols);
        std::cout << "Vizinho: (" << qRow << ", " << qCol << ")"<< std::endl;
    }

}
