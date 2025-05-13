#include "../include/MorphologicalTree.hpp"
#include "../include/AttributeComputedIncrementally.hpp"
#include "../include/ComputerMSER.hpp"
#include "../include/AttributeFilters.hpp"
#include "./Tests.hpp"

#include <iomanip> 

#include <iostream>
#include <fstream>
#include <stdexcept>

#include <vector>

int main(int argc, char const *argv[])
{
    ImagePtr image = getWonderImage();
    printImage(image);
    
    // Criar um ComponentTree
    MorphologicalTreePtr tree = std::make_shared<MorphologicalTree>(image, false);
    printTree(tree->getRoot());

    // Criar um AttributeComputedIncrementally::computerArea
    int n = tree->getNumNodes();	
    
    
    auto [attrsNamesArea, attrsArea] = AttributeComputedIncrementally::computeSingleAttribute(tree, GeometricAttribute::AREA);
    auto [attrNamesVolume, attrsVolume] = AttributeComputedIncrementally::computeSingleAttribute(tree, GeometricAttribute::VOLUME);
    auto [attrNamesLevel, attrsLevel] = AttributeComputedIncrementally::computeSingleAttribute(tree, GeometricAttribute::LEVEL);
    auto [attrNamesMeanLevel, attrsMeanLevel] = AttributeComputedIncrementally::computeSingleAttribute(tree, GeometricAttribute::MEAN_LEVEL);
    auto [attrNamesVarianceLevel, attrsVarianceLevel] = AttributeComputedIncrementally::computeSingleAttribute(tree, GeometricAttribute::VARIANCE_LEVEL);
    auto [attrsNamesDynamics, attrsDynamics] = AttributeComputedIncrementally::computeSingleAttribute(tree, GeometricAttribute::DYNAMICS);
    auto [attrNamesRectangularity, attrsRectangularity] = AttributeComputedIncrementally::computeSingleAttribute(tree, GeometricAttribute::RECTANGULARITY);
    auto [attrNamesRatio, attrsRatio] = AttributeComputedIncrementally::computeSingleAttribute(tree, GeometricAttribute::RATIO_WH);
    auto [attrNamesBoxWidth, attrsBoxWidth] = AttributeComputedIncrementally::computeSingleAttribute(tree, GeometricAttribute::BOX_WIDTH);
    auto [attrNamesBoxHeight, attrsBoxHeight] = AttributeComputedIncrementally::computeSingleAttribute(tree, GeometricAttribute::BOX_HEIGHT);
    auto [attrNamesOrientation, attrsOrientation] = AttributeComputedIncrementally::computeSingleAttribute(tree, GeometricAttribute::AXIS_ORIENTATION);
    auto [attrNamesInertia, attrsInertia] = AttributeComputedIncrementally::computeSingleAttribute(tree, GeometricAttribute::INERTIA);
    auto [attrNamesLength, attrsLength] = AttributeComputedIncrementally::computeSingleAttribute(tree, GeometricAttribute::LENGTH_MINOR_AXIS);
    auto [attrNamesEccentricity, attrsEccentricity] = AttributeComputedIncrementally::computeSingleAttribute(tree, GeometricAttribute::ECCENTRICITY);
    auto [attrNamesCompactness, attrsCompactness] = AttributeComputedIncrementally::computeSingleAttribute(tree, GeometricAttribute::COMPACTNESS);
    
    for(NodeMTPtr node : tree->getIndexNode()){
        int nodeIndex = node->getIndex();
        std::cout << "Atributo AREA do nó " << nodeIndex << ": " << attrsArea[nodeIndex] << std::endl;
        std::cout << "Atributo DYNAMICS do nó " << nodeIndex << ": " << attrsDynamics[nodeIndex] << std::endl;
        std::cout << "Atributo ORIENTATION do nó " << nodeIndex << ": " << attrsOrientation[nodeIndex] << std::endl;
        std::cout << "Atributo VOLUME do nó " << nodeIndex << ": " << attrsVolume[nodeIndex] << std::endl;
        std::cout << "Atributo LEVEL do nó " << nodeIndex << ": " << attrsLevel[nodeIndex] << std::endl;
        std::cout << "Atributo MEAN_LEVEL do nó " << nodeIndex << ": " << attrsMeanLevel[nodeIndex] << std::endl;
        std::cout << "Atributo VARIANCE_LEVEL do nó " << nodeIndex << ": " << attrsVarianceLevel[nodeIndex] << std::endl;
        std::cout << "Atributo RECTANGULARITY do nó " << nodeIndex << ": " << attrsRectangularity[nodeIndex] << std::endl;
        std::cout << "Atributo RATIO_WH do nó " << nodeIndex << ": " << attrsRatio[nodeIndex] << std::endl;
        std::cout << "Atributo BOX_WIDTH do nó " << nodeIndex << ": " << attrsBoxWidth[nodeIndex] << std::endl;
        std::cout << "Atributo BOX_HEIGHT do nó " << nodeIndex << ": " << attrsBoxHeight[nodeIndex] << std::endl;
        std::cout << "Atributo AXIS_ORIENTATION do nó " << nodeIndex << ": " << attrsOrientation[nodeIndex] << std::endl;
        std::cout << "Atributo INERTIA do nó " << nodeIndex << ": " << attrsInertia[nodeIndex] << std::endl;
        std::cout << "Atributo LENGTH_MINOR_AXIS do nó " << nodeIndex << ": " << attrsLength[nodeIndex] << std::endl;
        std::cout << "Atributo ECCENTRICITY do nó " << nodeIndex << ": " << attrsEccentricity[nodeIndex] << std::endl;
        std::cout << "Atributo COMPACTNESS do nó " << nodeIndex << ": " << attrsCompactness[nodeIndex] << std::endl;
       // printConnectedComponent(node, tree);
    }


    auto [attrNames, attrsPtr] = AttributeComputedIncrementally::computeAttributes(tree, {  GeometricAttribute::AREA, 
                                                                                            GeometricAttribute::DYNAMICS, 
                                                                                            GeometricAttribute::AXIS_ORIENTATION,
                                                                                            GeometricAttribute::VOLUME,
                                                                                            GeometricAttribute::LEVEL,
                                                                                            GeometricAttribute::MEAN_LEVEL,
                                                                                            GeometricAttribute::VARIANCE_LEVEL,
                                                                                            GeometricAttribute::RECTANGULARITY,
                                                                                            GeometricAttribute::RATIO_WH,
                                                                                            GeometricAttribute::LENGTH_MAJOR_AXIS,
                                                                                            GeometricAttribute::LENGTH_MINOR_AXIS,
                                                                                            GeometricAttribute::ECCENTRICITY,
                                                                                            GeometricAttribute::COMPACTNESS,
                                                                                            GeometricAttribute::BOX_WIDTH,
                                                                                            GeometricAttribute::BOX_HEIGHT,
                                                                                            GeometricAttribute::HU_MOMENT_1,
                                                                                            GeometricAttribute::HU_MOMENT_2,
                                                                                            GeometricAttribute::HU_MOMENT_3,
                                                                                            GeometricAttribute::HU_MOMENT_4,
                                                                                            GeometricAttribute::HU_MOMENT_5,
                                                                                            GeometricAttribute::HU_MOMENT_6,
                                                                                            GeometricAttribute::HU_MOMENT_7,
                                                                                            GeometricAttribute::CENTRAL_MOMENT_20,
                                                                                            GeometricAttribute::CENTRAL_MOMENT_02,
                                                                                            GeometricAttribute::CENTRAL_MOMENT_11,
                                                                                            GeometricAttribute::CENTRAL_MOMENT_30,
                                                                                            GeometricAttribute::CENTRAL_MOMENT_03,
                                                                                            GeometricAttribute::CENTRAL_MOMENT_21,
                                                                                            GeometricAttribute::CENTRAL_MOMENT_12} );

    
    // Depuração do mapeamento em `attributeNamesDelta`
    std::cout << "\n\nMapeamento de índices para atributos em attributeNames:" << std::endl;
    std::cout << "Número de atributos: " << attrNames->NUM_ATTRIBUTES << std::endl;
    for (const auto& pair : attrNames->indexMap) {
        
        Attribute attribute = pair.first;
        int offset = pair.second;
        std::cout << "\nAtributo: " << attrNames->toString(attribute) << ", Offset: " << offset << std::endl;

        // Exibir os valores dos atributos para cada nó
        for (NodeMTPtr node : tree->getIndexNode()) {
            int nodeIndex = node->getIndex();
            if (attrNames->indexMap.count(attribute)) {
                std::cout << "Node " << nodeIndex << " - " << attrNames->toString(attribute) << ": " 
                        //<< attrsPtr[nodeIndex * attrNames->NUM_ATTRIBUTES + offset]
                        << attrsPtr[attrNames->linearIndex(nodeIndex, attribute)]
                        << std::endl;
            }
        }

    }
    


    /*
    AttributeFilters filter(tree);
    float* attrDinamics = AttributeComputedIncrementally::computerAttributeByIndex(tree, GeometricAttribute::DYNAMICS);    
    ImagePtr imgOut = filter.filteringByExtinctionValue(tree, attrDinamics, 1);
    printImage(imgOut);
    */


    return 0;
}
 