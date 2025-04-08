
#include "Tests.hpp"
#include "../include/AdjacencyRelation.hpp"

#include "../include/NodeMT.hpp"
#include "../include/ImageUtils.hpp"

#include "../include/MorphologicalTree.hpp"
#include "../include/AttributeComputedIncrementally.hpp"


int main() {
    // Definição da imagem e parâmetros
    int numRows, numCols;
    int* img = getSimpleImage(numRows, numCols);
    int n = numRows * numCols;
    double radioAdj = 1.5;
    AdjacencyRelationPtr adj = std::make_shared<AdjacencyRelation>(numRows, numCols, 1);


    std::cout << "Resolution (cols x rows): " << numCols << " x " << numRows << std::endl;

    //imagem binaria
    
    printImage(img, numRows, numCols);

    // Criação das Component Trees
    MorphologicalTreePtr tree = std::make_shared<MorphologicalTree>(img, numRows, numCols);

    testComponentTree(tree, "ToS", img, numRows, numCols);
    std::cout << "Depth:" << tree->getDepth() << std::endl;
    std::cout << "--- Tree --- "<< std::endl;
    printTree(tree->getRoot());
    std::cout << std::endl;
    
    printConnectedComponents(tree);
    std::cout << std::endl;
    printMappingSC(tree, 3);

    std::vector<std::unordered_set<int>> countors = AttributeComputedIncrementally::extractCountors(tree);
    std::vector<std::vector<NodeMTPtr>> nodesByDepth = tree->getNodesByDepth();

    bool isEquals = true;
    for(int depth=tree->getDepth(); depth >= 0; depth--){
        std::vector<NodeMTPtr> nodesDepth = nodesByDepth[depth];

        int* contoursInc = new int[n]();
        int* contoursNonInc = new int[n]();
        for(NodeMTPtr node: nodesDepth){
            std::unordered_set<int> contourNode = countors[node->getIndex()];
            for(int p: contourNode){
                contoursInc[p] = 1;
            }
            int* imgBin = new int[n]();
            for(int p: node->getPixelsOfCC()){
                imgBin[p] = 1;
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
           // std::cout << "\nCBinaria" << std::endl;
           // printImage(imgBin, numRows, numCols);
            delete[] imgBin;
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
        std::cout << "Depth:"<< depth << "\tSão iguais:" << isEqualsDepth << std::endl;
        
        if(!isEqualsDepth){
            

            std::cout << "\nContorno não incremental" << std::endl;
            printImage(contoursNonInc, numRows, numCols, 3);
    
            std::cout << "\nContorno incremental" << std::endl;
            printImage(contoursInc, numRows, numCols, 3);
            break;
        }
        
        delete[] contoursInc;
        delete[] contoursNonInc;
    }
    
    //delete[] imgBin;

    if(isEquals){
        std::cout << "\nSão iguais" << std::endl;
    }else{
        std::cout << "\nSão diferentes" << std::endl;
    }
}