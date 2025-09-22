
#include "../include/AdjacencyRelation.hpp"
#include "../include/MorphologicalTree.hpp"
#include "../include/Common.hpp"
#include "../include/NodeMT.hpp"

#include "Tests.hpp"

#include <cassert>
#include <vector>
#include <iostream>


int main(){
    auto img = getSimpleImage();
    double radioAdj = 1.5;
    printImage(img);

    // Criação das Component Trees
    MorphologicalTreePtr maxtree = std::make_shared<MorphologicalTree>(img, true);
    
    
    auto imgMaxtree = maxtree->reconstructionImage();
    printTree(maxtree->getRoot());
    testComponentTree(maxtree, "maxtreeFZ sem grafo", imgMaxtree);

    NodeId nodeId = maxtree->getLeaves().front();// maxtree->getSC(28);
    NodeMT node = maxtree->proxy(5);
    
    std::cout << "\nNode - ID: " << node.getIndex() << ", Level: " << node.getLevel() << ", Area: " << node.getArea() << "\n" << std::endl;
    maxtree->prunning(node);
    

    auto imgPrunned = maxtree->reconstructionImage();
    printTree(maxtree->getRoot());
    printImage(imgPrunned);
    testComponentTree(maxtree, "maxtreeFZ sem grafo", imgPrunned);
    return 0;
}
