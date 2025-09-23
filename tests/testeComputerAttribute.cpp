
#include "../mmcfilters/utils/AdjacencyRelation.hpp"
#include "../mmcfilters/trees/MorphologicalTree.hpp"
#include "../mmcfilters/utils/Common.hpp"
#include "../mmcfilters/trees/NodeMT.hpp"
#include "../mmcfilters/attributes/AttributeComputedIncrementally.hpp"

#include "Tests.hpp"


#include <iomanip>
#include <iostream>
#include <vector>

using namespace mmcfilters;

int main(int argc, char const* argv[]) {
    (void)argc;
    (void)argv;

    ImageUInt8Ptr image = getWonderImage();
    printImage(image);

    MorphologicalTreePtr treePtr = std::make_shared<MorphologicalTree>(image, false);
    MorphologicalTree* tree = treePtr.get();
    printTree(treePtr->getRoot());


    auto [attrNamesRectangularity, attrsRectangularity] = AttributeComputedIncrementally::computeSingleAttribute(tree, Attribute::RECTANGULARITY);
    auto [attrsNamesArea, attrsArea] = AttributeComputedIncrementally::computeSingleAttribute(tree, Attribute::AREA);
    auto [attrNamesVolume, attrsVolume] = AttributeComputedIncrementally::computeSingleAttribute(tree, Attribute::VOLUME);
    auto [attrNamesLevel, attrsLevel] = AttributeComputedIncrementally::computeSingleAttribute(tree, Attribute::LEVEL);
    auto [attrNamesMeanLevel, attrsMeanLevel] = AttributeComputedIncrementally::computeSingleAttribute(tree, Attribute::MEAN_LEVEL);
    auto [attrNamesVarianceLevel, attrsVarianceLevel] = AttributeComputedIncrementally::computeSingleAttribute(tree, Attribute::VARIANCE_LEVEL);
    auto [attrsNamesGrayHeight, attrsGrayHeight] = AttributeComputedIncrementally::computeSingleAttribute(tree, Attribute::GRAY_HEIGHT);
    auto [attrNamesRatio, attrsRatio] = AttributeComputedIncrementally::computeSingleAttribute(tree, Attribute::RATIO_WH);
    auto [attrNamesBoxWidth, attrsBoxWidth] = AttributeComputedIncrementally::computeSingleAttribute(tree, Attribute::BOX_WIDTH);
    auto [attrNamesBoxHeight, attrsBoxHeight] = AttributeComputedIncrementally::computeSingleAttribute(tree, Attribute::BOX_HEIGHT);
    auto [attrNamesOrientation, attrsOrientation] = AttributeComputedIncrementally::computeSingleAttribute(tree, Attribute::AXIS_ORIENTATION);
    auto [attrNamesInertia, attrsInertia] = AttributeComputedIncrementally::computeSingleAttribute(tree, Attribute::INERTIA);
    auto [attrNamesLength, attrsLength] = AttributeComputedIncrementally::computeSingleAttribute(tree, Attribute::LENGTH_MINOR_AXIS);
    auto [attrNamesEccentricity, attrsEccentricity] = AttributeComputedIncrementally::computeSingleAttribute(tree, Attribute::ECCENTRICITY);
    auto [attrNamesCompactness, attrsCompactness] = AttributeComputedIncrementally::computeSingleAttribute(tree, Attribute::COMPACTNESS);

    for (NodeId nodeIndex: tree->getNodeIds()) {
        std::cout << "Atributo AREA do nó " << nodeIndex << ": " << attrsArea[nodeIndex] << std::endl;
        std::cout << "Atributo GRAY_HEIGHT do nó " << nodeIndex << ": " << attrsGrayHeight[nodeIndex] << std::endl;
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
    }

    auto [attrNames, attrsPtr] = AttributeComputedIncrementally::computeAttributes(tree, {
                                                                                            AttributeGroup::ALL
    });

    // Depuração do mapeamento em `attributeNamesDelta`
    std::cout << "\n\nMapeamento de índices para atributos em attributeNames:" << std::endl;
    std::cout << "Número de atributos: " << attrNames->NUM_ATTRIBUTES << std::endl;
    for (const auto& pair : attrNames->indexMap) {
        
        Attribute attribute = pair.first;
        int offset = pair.second;
        std::cout << "\nAtributo: " << attrNames->toString(attribute) << ", Offset: " << offset << std::endl;

        // Exibir os valores dos atributos para cada nó
        for (NodeId nodeIndex: tree->getNodeIds()) {
            if (attrNames->indexMap.count(attribute)) {
                std::cout << "Node " << nodeIndex << " - " << attrNames->toString(attribute) << ": " 
                        //<< attrsPtr[nodeIndex * attrNames->NUM_ATTRIBUTES + offset]
                        << attrsPtr[attrNames->linearIndex(nodeIndex, attribute)]
                        << std::endl;
            }
        }
    }


    auto [attrsNamesDelta, attrsAreaDelta] = AttributeComputedIncrementally::computeSingleAttributeWithDelta(tree, Attribute::GRAY_HEIGHT, 0, "zero-padding");

    for (const auto& pair : attrsNamesDelta->indexMap) {
        
        AttributeKey attributeKey = pair.first;
        int offset = pair.second;
        std::cout << "\nAtributo: " << attrsNamesDelta->toString(attributeKey) << ", Offset: " << offset << std::endl;
        // Exibir os valores dos atributos para cada nó
        for (NodeId nodeIndex: tree->getNodeIds()) {
            std::cout << "Node " << nodeIndex << " - " << attrsNamesDelta->toString(attributeKey) << ": " << attrsAreaDelta[attrsNamesDelta->linearIndex(nodeIndex, attributeKey)] << std::endl;    
        }
    }


    ImageFloatPtr attributeMapping = ImageFloat::fromExternal(static_cast<float*>(attrsRectangularity.get()), image->getNumRows(), image->getNumCols());
    MorphologicalTreePtr treeAttrMap = MorphologicalTree::createFromAttributeMapping(attributeMapping, image, false, 1.5);
    printTree(treeAttrMap->getRoot());
    return 0;
}
