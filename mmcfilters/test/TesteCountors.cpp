
#include "Tests.hpp"
#include "../include/AdjacencyRelation.hpp"

#include "../include/NodeMT.hpp"
#include "../include/ImageUtils.hpp"

#include "../include/MorphologicalTree.hpp"
#include "../include/AttributeComputedIncrementally.hpp"


int main() {
    // Definição da imagem e parâmetros
    int numRows, numCols;
    int* img = getPassatImage(numRows, numCols);
    int n = numRows * numCols;
    double radioAdj = 1.5;
    AdjacencyRelationPtr adj = std::make_shared<AdjacencyRelation>(numRows, numCols, 1);


    std::cout << "Resolution: " << numCols << " x " << numRows << std::endl;

    //imagem binaria
    bool* imgBin = new bool[n]();


    // Criação das Component Trees
    MorphologicalTreePtr tree = std::make_shared<MorphologicalTree>(img, numRows, numCols);
    std::vector<std::unordered_set<int>> countors = AttributeComputedIncrementally::extractCountors(tree);
    std::vector<std::vector<NodeMTPtr>> nodesByDepth = tree->getNodesByDepth();
    bool isEquals = true;
    for(int depth=tree->getDepth(); depth >= 0; depth--){
        std::vector<NodeMTPtr> nodesDepth = nodesByDepth[depth];

        std::vector<bool> contoursInc(numRows * numCols, false);
        for(NodeMTPtr node: nodesDepth){
            std::unordered_set<int> contourNode = countors[node->getIndex()];
            for(int p: contourNode){
                contoursInc[p] = true;
            }
            for(int p: node->getCNPs()){
                imgBin[p] = true;
            }
        }
        
        std::vector<bool> contoursNonInc(numRows * numCols, false);
        for (int p=0; p < numRows*numCols; p++ ) {
            for (int q : adj->getAdjPixels(p)) {
                if (imgBin[p] && !imgBin[q]) {
                    contoursNonInc[p] = true;
                }
            }
        }
        bool isEqualsDepth = true;
        for (int p=0; p < (numRows*numCols); p++ ) {
            if(contoursNonInc[p] != contoursInc[p]){
                isEquals = false;
                isEqualsDepth = false;
                std::pair<int, int> point = ImageUtils::to2D(p, numCols);
                std::cout << "(" << point.first << ", " << point.second << ")\n";
            }
        }
        std::cout << "Depth:"<< depth << "\tSão iguais:" << isEqualsDepth << std::endl;
        
    }

    if(isEquals){
        std::cout << "\nSão iguais" << std::endl;
    }else{
        std::cout << "\nSão diferentes" << std::endl;
    }
}