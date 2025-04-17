
#include "Tests.hpp"
#include "../include/AdjacencyRelation.hpp"

#include "../include/NodeMT.hpp"
#include "../include/ImageUtils.hpp"

#include "../include/MorphologicalTree.hpp"
#include "../include/AttributeComputedIncrementally.hpp"
#include "../../external/stb/stb_image.h"
#include "../../external/stb/stb_image_write.h"

#include <chrono>
#include <iostream>

int* openImage(std::string filename, int& numRows, int& numCols){
    //std::cout << "filename:" << filename << std::endl;
    int nchannels;
    unsigned char* data = stbi_load(filename.c_str(), &numCols, &numRows, &nchannels, 1);
    
    int* img = new int[numCols * numRows];
    for (int i = 0; i < numCols * numRows; i++) {
        img[i] = static_cast<int>(data[i]);  // Converte de `unsigned char` para `int`
    }
    stbi_image_free(data); // Liberar a memória da imagem carregada

    return img;
}

int main(int argc, char* argv[]) {
    // Definição da imagem e parâmetros
    int numRows, numCols;
    //int* img = getSimpleImage(numRows, numCols);
    
    if(argc != 3){
        std::cout << "Execute assim: " << argv[0] << " <ToS_type> <filename>" << std::endl;
        return 1;
    }
    int* img = openImage(argv[2], numRows, numCols);
    std::string tosType = (std::string(argv[1])=="self-dual"? "self-dual":"4c8c");

    int n = numRows * numCols;
    AdjacencyRelationPtr adj = std::make_shared<AdjacencyRelation>(numRows, numCols, 1);

    std::cout << "\nImage:"<< argv[2] << " \tResolution (cols x rows): " << numCols << " x " << numRows << std::endl;
    
    //printImage(img, numRows, numCols);

    // Criação das Component Trees
    auto start = std::chrono::high_resolution_clock::now();
    MorphologicalTreePtr tree = std::make_shared<MorphologicalTree>(img, numRows, numCols, tosType);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Tempo para construir a árvore: " << elapsed.count() << " segundos" << std::endl;
    
    std::cout << "Depth:" << tree->getDepth() << ", |nodes|:" << tree->getNumNodes() << std::endl;
    
    
    //testComponentTree(tree, "ToS", img, numRows, numCols);
    //std::cout << std::endl;
    
    //printConnectedComponents(tree);
    //std::cout << std::endl;

    //std::cout << "--- Tree --- "<< std::endl;
    //printTree(tree->getRoot());
    //std::cout << std::endl;
    //printMappingSC(tree, 3);
    

    std::vector<std::unordered_set<int>> countors = AttributeComputedIncrementally::extractCountors(tree);
    std::vector<std::vector<NodeMTPtr>> nodesByDepth = tree->getNodesByDepth();
    int* imgBin = new int[n];
    int* contoursInc = new int[n];
    int* contoursNonInc = new int[n];
    bool isEquals = true;
    for(int depth=tree->getDepth(); depth >= 0; depth--){
        std::vector<NodeMTPtr> nodesDepth = nodesByDepth[depth];

        for(NodeMTPtr node: nodesDepth){

            for(int p=0; p < n; p++){
                imgBin[p] = 0;
                contoursInc[p] = 0;
                contoursNonInc[p] = 0;
            }
            
            std::unordered_set<int> contourNode = countors[node->getIndex()];
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
        
        
    }
    
    if(isEquals){
        std::cout <<"\tFilename:" << argv[1] << "\tIs equals? Yes" << std::endl;
    }else{
        std::cout <<"\tFilename:" << argv[1] << "\tIs equals? No" << std::endl;
    }
    delete[] imgBin;
    delete[] contoursInc;
    delete[] contoursNonInc;
    return 0;
}