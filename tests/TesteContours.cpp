
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
    
    if(argc != 3){
        std::cout << "Execute assim: " << argv[0] << " <ToS_type> <filename>" << std::endl;
        return 1;
    }
   // ImagePtr image = openImage(argv[2]);
    
    int numRows = image->getNumRows();
    int numCols = image->getNumCols();

    AdjacencyRelationPtr adj = std::make_shared<AdjacencyRelation>(numRows, numCols, 1);

    //std::cout << "\nImage:"<< argv[2] << " \tResolution (cols x rows): " << numCols << " x " << numRows << std::endl;
    
    //printImage(img, numRows, numCols);

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
    
    
    
    //std::cout << "Depth:" << tree->getDepth() << ", |nodes|:" << tree->getNumNodes() << std::endl;
    
    
    //testComponentTree(tree, "ToS", img, numRows, numCols);
    //std::cout << std::endl;
    
    //printConnectedComponents(tree);
    //std::cout << std::endl;

    //std::cout << "--- Tree --- "<< std::endl;
    //printTree(tree->getRoot());
    //std::cout << std::endl;
    //printMappingSC(tree, 3);
    
    auto contoursCT = ContoursComputedIncrementally::extractCompactContours(tree.get());
    //std::vector<std::unordered_set<int>> countors =  ContoursComputedIncrementally::extractCountors(tree);
    
    ImageUInt8Ptr imgBin = ImageUInt8::create(numRows, numCols, 0);
    ImageUInt8Ptr contoursInc = ImageUInt8::create(numRows, numCols, 0);
    ImageUInt8Ptr contoursNonInc = ImageUInt8::create(numRows, numCols, 0);

    bool isEquals = true;
    for (auto&& [node, contourNode] : contoursCT.contoursLazy()) {
        std::list<int> cc;
        for (int p : tree->getCNPsById(node)) {
            cc.push_back(p);
        }

        // contorno incremental
        for (int p : contourNode) {
            (*contoursInc)[p] = 1;
        }

        // contorno não incremental
        for (int p : cc) {
            if (tree->isDescendant(tree->getSCById(p), node)) {
                (*imgBin)[p] = 1;
            }
        }

        for (int p : cc) {
            auto [row, col] = ImageUtils::to2D(p, numCols);
            if ((*imgBin)[p] == 1 && (row == 0 || col == 0 || col == numCols - 1 || row == numRows - 1)) {
                (*contoursNonInc)[p] = 1;
            } else {
                for (int q : adj->getAdjPixels(p)) {
                    if ((*imgBin)[p] == 1 && (*imgBin)[q] == 0) {
                        (*contoursNonInc)[p] = 1;
                    }
                }
            }
        }

        bool isEqualsCC = true;
        for (int p : cc) {
            if ((*contoursNonInc)[p] != (*contoursInc)[p]) {
                isEquals = false;
                isEqualsCC = false;
                auto point = ImageUtils::to2D(p, numCols);
                std::cout << "(row, col) = (" << point.first << ", " << point.second << ")\n";
            }
        }
        std::cout << "CC(nodeID):" << node << "\t is equals? " << isEqualsCC << std::endl;

        if (!isEqualsCC) {
            std::cout << "\nContorno não incremental" << std::endl;
            printImage(contoursNonInc, 3);

            std::cout << "\nContorno incremental" << std::endl;
            printImage(contoursInc, 3);
        }

        // limpa as estruturas para o próximo nó
        for (int p : cc) {
            (*imgBin)[p] = 0;
            (*contoursInc)[p] = 0;
            (*contoursNonInc)[p] = 0;
        }
    }
  
/*
    
        for(NodeMTPtr node: tree->getRoot()->getIteratorPostOrderTraversal()){

            for(int p=0; p < n; p++){
                imgBin[p] = 0;
                contoursInc[p] = 0;
                contoursNonInc[p] = 0;
            }
            
            
            std::unordered_set<int> contourNode = contoursMT->getContour(node);//countors[node->getIndex()];
            for(int p: contourNode){
                contoursInc[p] = 1;
            }
            
            
            for (int row = 0; row < numRows; ++row) {
                for (int col = 0; col < numCols; ++col) {
                    if(tree->isDescendant(tree->getSC(ImageUtils::to1D(row, col, numCols)), node)){
                        imgBin[ImageUtils::to1D(row, col, numCols)] = 1;
                    }
                }
            }

            for (int p=0; p < n; p++ ) {
                auto [row, col] = ImageUtils::to2D(p, numCols);
                if(imgBin[p]==1 && (row ==0 || col ==0 || col == numCols-1 || row == numRows-1)){
                    contoursNonInc[p] = 1;
                }else{
                    for (int q : adj->getAdjPixels(p)) {
                        if (imgBin[p]==1 && imgBin[q]==0) {
                            contoursNonInc[p] = 1;
                        }
                    }
                }
            }

            
        
            bool isEqualsDepth = true;
            for (int p=0; p < n; p++ ) {
                if(contoursNonInc[p] != contoursInc[p]){
                    isEquals = false;
                    isEqualsDepth = false;
                    std::pair<int, int> point = ImageUtils::to2D(p, numCols);
                    std::cout << "(row, col) = (" << point.first << ", " << point.second << ")\n";
                }
            }
            std::cout << "nodeID:"<< node->getIndex() << "\tSão iguais:" << isEqualsDepth << std::endl;
            if(!isEqualsDepth){
                std::cout << "\nImagem binaria" << std::endl;
                printImage(imgBin, numRows, numCols, 3);

                std::cout << "\nContorno não incremental" << std::endl;
                printImage(contoursNonInc, numRows, numCols, 3);
        
                std::cout << "\nContorno incremental" << std::endl;
                printImage(contoursInc, numRows, numCols, 3);
                
                
                break;
            }
           
        }*/
       
        //std::cout << "Depth:" << tree->getDepth() << ", |nodes|:" << tree->getNumNodes() << std::endl;
    if(isEquals){
        std::cout <<"Filename:" << argv[2] <<", TreeType: "<< treeType <<  ", |nodes|:" << tree->getNumNodes() << "\tIs equals? Yes" << std::endl;
    }else{
        std::cout <<"Filename:" << argv[2] << "\tIs equals? No" << std::endl;
    }
    return 0;
}
