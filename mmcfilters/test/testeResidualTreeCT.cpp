#include "../include/ComponentTree.hpp"
#include "../include/AttributeComputedIncrementallyCT.hpp"
#include "../include/AttributeOpeningPrimitivesFamilyCT.hpp"
#include "../include/ResidualTreeCT.hpp"
#include "../include/UltimateAttributeOpeningCT.hpp"
#include "Tests2.hpp"


#include <iostream>
#include <iomanip> 

#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>

int main(int argc, char const *argv[])
{
    
    
    ImageUInt8Ptr img = getPassatImage();
    
    printImage(img);
    std::cout << "img_pointer ok" << std::endl;
    
    // Criar um ComponentTree
    ComponentTreePtr tree = std::make_shared<ComponentTree>(img, true);
    std::cout << "tree ok" << std::endl;
    

    // Criar um AttributeComputedIncrementally::computerArea
    const int n = tree->getNumNodes();	
    auto [names, attr] = AttributeComputedIncrementallyCT::computeSingleAttribute(tree, Attribute::BOX_HEIGHT); //size: n * numAttribute
    std::cout << "attributes ok" << std::endl;

    // Criar um AttributeOpeningPrimitivesFamily
    int maxCriterion = img->getNumRows(); 
    std::shared_ptr<AttributeOpeningPrimitivesFamilyCT> primitives = std::make_shared<AttributeOpeningPrimitivesFamilyCT>(tree, attr, maxCriterion);
    std::cout << "primitives ok" << std::endl;
    
    ResidualTreeCT* residualTree = new ResidualTreeCT(primitives);
    std::cout << "residualTree ok" << std::endl;


    ImageUInt8Ptr imgPos = residualTree->getPositiveResidues();
    printImage(imgPos);
    std::cout << "imgPos ok" << std::endl;
    
    ImageUInt8Ptr imgNeg = residualTree->getNegativeResidues();
    printImage(imgNeg);
    std::cout << "imgNeg ok" << std::endl;
    

    ImageUInt8Ptr imgRec = residualTree->reconstruction();
    printImage(imgRec);

    ImageUInt8Ptr contrast = residualTree->getMaxConstrastImage();
    //printImage(contrast);
    std::cout << "contrast ok" << std::endl;
    
    UltimateAttributeOpeningCT *uao = new UltimateAttributeOpeningCT(tree, attr);
    uao->execute(maxCriterion);
    ImageUInt8Ptr contrastUAO = residualTree->getMaxConstrastImage();
    //printImage(contrastUAO);

    
    printf("isEquals: %d\n", contrast->isEqual(contrastUAO));



    

    return 0;
}
