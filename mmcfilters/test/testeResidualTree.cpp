#include "../include/MorphologicalTree.hpp"
#include "../include/AttributeComputedIncrementally.hpp"
#include "../include/AttributeOpeningPrimitivesFamily.hpp"
#include "../include/ResidualTree.hpp"
#include "../include/UltimateAttributeOpening.hpp"
#include "Tests.hpp"


#include <iostream>
#include <iomanip> 

#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>

int main(int argc, char const *argv[])
{
    
    
    ImagePtr img = getPassatImage();
    
    printImage(img);
    std::cout << "img_pointer ok" << std::endl;
    
    // Criar um ComponentTree
    MorphologicalTreePtr tree = std::make_shared<MorphologicalTree>(img, true);
    std::cout << "tree ok" << std::endl;
    

    // Criar um AttributeComputedIncrementally::computerArea
    const int n = tree->getNumNodes();	
    float* attr = AttributeComputedIncrementally::computerAttributeByIndex(tree, GeometricAttribute::BOX_HEIGHT); //size: n * numAttribute
    std::cout << "attributes ok" << std::endl;

    // Criar um AttributeOpeningPrimitivesFamily
    int maxCriterion = img->numRows; 
    AttributeOpeningPrimitivesFamily* primitives = new AttributeOpeningPrimitivesFamily(tree, attr, maxCriterion);
    std::cout << "primitives ok" << std::endl;
    
    ResidualTree* residualTree = new ResidualTree(primitives);
    std::cout << "residualTree ok" << std::endl;


    ImagePtr imgPos = residualTree->getPositiveResidues();
    printImage(imgPos);
    std::cout << "imgPos ok" << std::endl;
    
    ImagePtr imgNeg = residualTree->getNegativeResidues();
    printImage(imgNeg);
    std::cout << "imgNeg ok" << std::endl;
    

    ImagePtr imgRec = residualTree->reconstruction();
    printImage(imgRec);

    ImagePtr contrast = residualTree->getMaxConstrastImage();
    //printImage(contrast);
    std::cout << "contrast ok" << std::endl;
    
    UltimateAttributeOpening *uao = new UltimateAttributeOpening(tree, attr);
    uao->execute(maxCriterion);
    ImagePtr contrastUAO = residualTree->getMaxConstrastImage();
    //printImage(contrastUAO);

    
    printf("isEquals: %d\n", contrast->isEqual(contrastUAO));



    

    return 0;
}
