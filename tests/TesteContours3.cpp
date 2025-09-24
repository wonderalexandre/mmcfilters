
#include "Tests.hpp"
#include "../mmcfilters/utils/AdjacencyRelation.hpp"

#include "../mmcfilters/trees/NodeMT.hpp"
#include "../mmcfilters/utils/Common.hpp"

#include "../mmcfilters/trees/MorphologicalTree.hpp"
#include "../mmcfilters/attributes/AttributeComputedIncrementally.hpp"
#include "../mmcfilters/contours/ContoursComputedIncrementally.hpp"

using namespace mmcfilters;

int main(int argc, char* argv[]) {
    // Definição da imagem e parâmetros
    
    ImageUInt8Ptr image = getSimpleImage();
    
    if(argc != 2){
        std::cout << "Execute assim: " << argv[0] << " <ToS_type>" << std::endl;
        return 1;
    }
   // ImagePtr image = openImage(argv[2]);
    
    int numRows = image->getNumRows();
    int numCols = image->getNumCols();

    AdjacencyRelationPtr adj = std::make_shared<AdjacencyRelation>(numRows, numCols, 1);

    //std::cout << "\nImage:"<< argv[2] << " \tResolution (cols x rows): " << numCols << " x " << numRows << std::endl;
    
    printImage(image);

    // Criação das Component Trees
    
    MorphologicalTreePtr tree = nullptr;
    std::string treeType = std::string(argv[1]);
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
    
    printImage(tree->reconstructionImage());
    
    //std::cout << "Depth:" << tree->getDepth() << ", |nodes|:" << tree->getNumNodes() << std::endl;
    
    
    //testComponentTree(tree, "ToS", img, numRows, numCols);
    //std::cout << std::endl;
    
    printConnectedComponents(tree);
    std::cout << std::endl;

    std::cout << "--- Tree --- "<< std::endl;
    printTree(tree->getRoot());
    std::cout << std::endl;
    printMappingSC(tree, 3);
    
    auto contoursCT = ContoursComputedIncrementally::extractCompactContours(tree.get());
    //ImagePtr imgContours = Image::create(numRows, numCols, 0);
    //bool isEquals = true;
    
    for (auto [node, contourPixels] : contoursCT.contoursLazy()) {
        ImageUInt8Ptr imgContours = ImageUInt8::create(numRows, numCols, 0);
        
        //contorno incremental
        for(int p: contourPixels){
            (*imgContours)[p] = 1;
        }
        
        std::cout << "\nNode:"<< node << std::endl;
        printConnectedComponent(tree->proxy(node), tree);
        printImage(imgContours, 3);
    }

    /*
    contoursMT.visitContours(tree, [&](NodeMTPtr node, const std::unordered_set<int>& contourNode) {
       
        ImageUInt8Ptr imgContours = ImageUInt8::create(numRows, numCols, 0);
        //contorno incremental
        for(int p: contourNode){
            (*imgContours)[p] = 1;
        }
        
        std::cout << "\nNode:"<< node->getIndex() << std::endl;
        printImage(imgContours, 3);
                
    });
    

    ImageUInt8Ptr imgContours2 = ImageUInt8::create(numRows, numCols, 0);
    for(int p: contoursMT.getContour(tree->getNodeByIndex(3))){
        (*imgContours2)[p] = 1;
    }
    
   std::cout << "\nNode 3" << std::endl;
    printImage(imgContours2, 3);
*/
    return 0;
}
