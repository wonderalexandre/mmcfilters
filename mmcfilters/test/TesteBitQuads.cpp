
#include "Tests.hpp"
#include "../include/AdjacencyRelation.hpp"

#include "../include/NodeMT.hpp"
#include "../include/ImageUtils.hpp"
#include "../include/Common.hpp"

#include "../include/MorphologicalTree.hpp"
#include "../include/AttributeComputedIncrementally.hpp"
#include "../include/ComputerAttributeBasedBitQuads.hpp"


int main(int argc, char* argv[]) {
    // Definição da imagem e parâmetros
    ImageUInt8Ptr image = getICIP14Image();
    
    
    printImage(image);
    


    // Criação das Component Trees
    MorphologicalTreePtr tree = nullptr;
    std::string treeType = "maxtree"; 
    if(treeType=="mintree"){
        tree = std::make_shared<MorphologicalTree>(image, false);
        //std::cout << "mintree" << std::endl;
    }else if(treeType=="maxtree"){
        tree = std::make_shared<MorphologicalTree>(image, true);
        //std::cout << "maxtree" << std::endl;
    }else{
        treeType = treeType=="ToS-4c8c"? "4c8c":"self-dual";
        tree = std::make_shared<MorphologicalTree>(image, treeType);
        //std::cout << "tree of shapes - "<< treeType << std::endl;
    }
    
    std::cout << "--- Tree --- "<< std::endl;
    printTree(tree->getRoot());
    std::cout << std::endl;
    

    ComputerAttributeBasedBitQuads computer(tree);
    std::vector<AttributeBasedBitQuads> attr = computer.getAttributes();
    std::cout << "Patterns: " << std::endl;
    for (NodeMTPtr node : tree->getRoot()->getNodesDescendants()) {
        std::cout << "Node ID: " << node->getIndex() << ",\tLevel: " << node->getLevel() << ",\tPatterns: " << attr[node->getIndex()].printPattern() << std::endl;
    }

    return 0;
}