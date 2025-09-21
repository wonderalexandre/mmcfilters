
#include "Tests2.hpp"
#include "../include/AdjacencyRelation.hpp"

#include "../include/Common.hpp"

#include "../include/ComponentTree.hpp"
#include "../include/AttributeComputedIncrementallyCT.hpp"
#include "../include/ComputerAttributeBasedBitQuadsCT.hpp"


int main(int argc, char* argv[]) {
    // Definição da imagem e parâmetros
    ImageUInt8Ptr image = getICIP14Image();
    
    
    printImage(image);
    

    ComponentTreePtr tree = nullptr;
    std::string treeType = "mintree";
    if(treeType=="mintree"){
        tree = std::make_shared<ComponentTree>(image, false);
        //std::cout << "mintree" << std::endl;
    }else if(treeType=="maxtree"){
        tree = std::make_shared<ComponentTree>(image, true);
        //std::cout << "maxtree" << std::endl;
    }else{
        treeType = treeType=="ToS-4c8c"? "4c8c":"self-dual";
        tree = std::make_shared<ComponentTree>(image, treeType);
        //std::cout << "tree of shapes - "<< treeType << std::endl;
    }
    
    std::cout << "--- Tree: " << treeType << " --- ["<<tree->getTreeType() << "]" << std::endl;
    printTree(tree->getRoot());
    std::cout << std::endl;
    

    ComputerAttributeBasedBitQuadsCT computer(tree);
    std::vector<AttributeBasedBitQuadsCT> attr = computer.getAttributes();
    std::cout << "Patterns: " << std::endl;
    for (NodeId node : tree->getNodeIds()) {
        std::cout << "Node ID: " << node << ",\tLevel: " << tree->getLevelById(node) << ",\tPatterns: " << attr[node].printPattern() << std::endl;
    }

    return 0;
}