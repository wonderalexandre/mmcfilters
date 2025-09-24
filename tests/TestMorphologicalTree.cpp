
#include "../mmcfilters/utils/AdjacencyRelation.hpp"
#include "../mmcfilters/trees/MorphologicalTree.hpp"
#include "../mmcfilters/utils/Common.hpp"
#include "../mmcfilters/trees/NodeMT.hpp"

#include "Tests.hpp"

#include <cassert>
#include <vector>
#include <iostream>

using namespace mmcfilters;


int main(){
    auto img = getSimpleImage();
    printImage(img);

    // Criação das Component Trees
    MorphologicalTreePtr maxtree = std::make_shared<MorphologicalTree>(img, true);
    
    
    auto imgMaxtree = maxtree->reconstructionImage();
    printTree(maxtree->getRoot());
    testComponentTree(maxtree, "maxtreeFZ sem grafo", imgMaxtree);

    NodeMT node = maxtree->proxy(5);
    
    std::cout << "\nNode - ID: " << node.getIndex() << ", Level: " << node.getLevel() << ", Area: " << node.getArea() << "\n" << std::endl;
    maxtree->prunning(node);
    

    auto imgPrunned = maxtree->reconstructionImage();
    printTree(maxtree->getRoot());
    printImage(imgPrunned);
    testComponentTree(maxtree, "maxtreeFZ sem grafo", imgPrunned);
    return 0;
}
