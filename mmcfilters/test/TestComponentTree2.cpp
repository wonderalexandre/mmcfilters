
#include "../include/AdjacencyRelation.hpp"
#include "../include/ComponentTree.hpp"
#include "../include/Common.hpp"

#include "Tests2.hpp"

#include <cassert>
#include <vector>
#include <iostream>


int main(){
    auto img = getLenaCropImage();
    double radioAdj = 1.5;
    printImage(img);

    // Criação das Component Trees
    AdjacencyRelationPtr adj =std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    ComponentTreePtr maxtree = std::make_shared<ComponentTree>(img, false, adj);

    auto imgMaxtree = maxtree->reconstructionImage();
    printTree(maxtree->getRoot());
    testComponentTree(maxtree, "maxtreeFZ sem grafo", imgMaxtree);

    NodeId nodeId = maxtree->getLeaves().front();// maxtree->getSC(28);
    NodeCT node = maxtree->proxy(7);
    std::cout << "\nNode - ID: " << node.getIndex() << ", Level: " << node.getLevel() << ", Area: " << node.getArea() << "\n" << std::endl;
    maxtree->prunning(node);


    auto imgPrunned = maxtree->reconstructionImage();
    printTree(maxtree->getRoot());
    printImage(imgPrunned);
    testComponentTree(maxtree, "maxtreeFZ sem grafo", imgPrunned);
    return 0;
}
