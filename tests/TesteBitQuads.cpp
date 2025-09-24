
#include "Tests.hpp"
#include "../mmcfilters/utils/AdjacencyRelation.hpp"

#include "../mmcfilters/utils/Common.hpp"

#include "../mmcfilters/trees/MorphologicalTree.hpp"
#include "../mmcfilters/attributes/AttributeComputedIncrementally.hpp"
#include "../mmcfilters/attributes/ComputerAttributeBasedBitQuads.hpp"

using namespace mmcfilters;

int main() {
    // Definição da imagem e parâmetros
    ImageUInt8Ptr image = getICIP14Image();
    
    
    printImage(image);
    

    MorphologicalTreePtr tree = nullptr;
    std::string treeType = "mintree";
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
    
    std::cout << "--- Tree: " << treeType << " --- ["<<tree->getTreeType() << "]" << std::endl;
    printTree(tree->getRoot());
    std::cout << std::endl;
    

    ComputerAttributeBasedBitQuads computer(tree);
    std::vector<AttributeBasedBitQuads> attr = computer.getAttributes();
    std::cout << "Patterns: " << std::endl;
    for (NodeId node : tree->getNodeIds()) {
        std::cout << "Node ID: " << node << ",\tLevel: " << tree->getLevelById(node) << ",\tPatterns: " << attr[node].printPattern() << std::endl;
    }

    return 0;
}
