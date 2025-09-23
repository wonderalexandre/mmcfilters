#include "../mmcfilters/trees/MorphologicalTree.hpp"
#include "../mmcfilters/attributes/AttributeComputedIncrementally.hpp"
#include "../mmcfilters/attributes/ComputerMSER.hpp"
#include "../mmcfilters/filters/AttributeFilters.hpp"
#include "./Tests.hpp"

#include <iomanip> 

#include <iostream>
#include <fstream>
#include <stdexcept>

#include <vector>

using namespace mmcfilters;

int main(int argc, char const *argv[])
{
    ImageUInt8Ptr image = getSimpleImage();

   // printImage(img_pointer, numRows, numCols);
    std::cout << "img_pointer ok" << std::endl;
    
    // Criar um MorphologicalTree
    MorphologicalTreePtr tree = std::make_shared<MorphologicalTree>(image, true);
    std::cout << "tree ok" << std::endl;

    printTree(tree->getRoot());

    // Criar um AttributeComputedIncrementally::computerArea
    int n = tree->getNumNodes();	
    
    int delta = 3;
    auto [attrNames, attrsPtr] = AttributeComputedIncrementally::computerBasicAttributes(tree);
    std::cout << "attributes ok" << std::endl;

    // Depuração do mapeamento em `attributeNamesDelta`
    std::cout << "Mapeamento de índices para atributos em attributeNames:" << std::endl;
    for (const auto& pair : attrNames.indexMap) {
        
        Attribute attribute = pair.first;
        int offset = pair.second;
        std::cout << "\nAtributo: " << attrNames.toString(attribute) << ", Offset: " << offset << std::endl;

        // Exibir os valores dos atributos para cada nó
        for (NodeMTPtr node : tree->getIndexNode()) {
            int nodeIndex = node->getIndex();
            if (attrNames.indexMap.count(attribute)) {
                std::cout << "Node " << nodeIndex << " - " << attrNames.toString(attribute) << ": " 
                        << attrsPtr[nodeIndex + offset] 
                        << std::endl;
            }
        }

    }

    AttributeFilters filter(tree);
    float* attrDinamics = AttributeComputedIncrementally::computerAttributeByIndex(tree, Attribute::DYNAMICS);    
    ImageUInt8Ptr imgOut = filter.filteringByExtinctionValue(tree, attrDinamics, 1);
    printImage(imgOut);

    return 0;
}
 