#include "../include/ComponentTree.hpp"
#include "../include/AttributeComputedIncrementally.hpp"
#include "../include/AttributeOpeningPrimitivesFamily.hpp"
#include "../include/ResidualTree.hpp"
#include "../include/UltimateAttributeOpening.hpp"

/*



#include "../include/PrimitivesFamily.hpp"
#include "../include/NodeCT.hpp"
*/


#include <iostream>
#include <iomanip> 

#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>

int* readImageFromFile(const std::string& filename, int numRows, int numCols) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Não foi possível abrir o arquivo");
    }

    std::vector<int> data;
    int value;
    while (file >> value) {
        data.push_back(value);
    }

    if (data.size() != numRows * numCols) {
        throw std::runtime_error("O tamanho dos dados não corresponde ao número de linhas e colunas fornecidos");
    }

    int* img_pointer = new int[data.size()];
    std::copy(data.begin(), data.end(), img_pointer);

    return img_pointer;
}

void printImage(int* img, int numRows, int numCols) {
    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            int pos = i * numCols + j;
            std::cout << std::setw(3) << img[pos] << " "; // Use std::setw(3) para formatar a saída
        }
        std::cout << std::endl;
    }
}

int main(int argc, char const *argv[])
{
    
    
    const int numRows = 25;
    const int numCols = 25;
    const int size = numRows * numCols;

    // Definir os dados da matriz
    int data[size] = {
        203, 203, 203, 171, 153, 157, 165, 165, 161, 160, 165, 165, 165, 161, 163, 181, 203, 203, 203, 203, 203, 203, 203, 203, 203,
        203, 203, 128,  78,  78,  78,  78,  78,  78,  78,  78,  78,  78,  78,  78,  78, 105, 190, 203, 203, 203, 203, 203, 203, 203,
        203, 177,  78,  79, 101, 113,  91,  76, 107, 113,  96,  78,  78,  78,  78,  78,  78, 134, 203, 203, 185, 158, 185, 203, 203,
        203, 138,  78, 100, 117,  67,  38,  38,  69, 126, 126,  78,  78,  78,  78,  78,  78, 116, 203, 153,  56,  54,  68, 200, 203,
        203, 105,  78, 109,  58,  38,  38,  38,  71, 126, 123,  78,  78,  78,  78,  78,  78, 116, 184,  63,  54,  54,  54, 181, 203,
        203, 102,  78, 106,  38,  38,  38,  41, 105, 104,  87,  78,  78,  78,  78,  78,  78, 143,  97,  54,  54,  54,  54, 189, 203,
        203,  89,  78,  98,  38,  38,  49,  92,  84,  78,  78, 119, 153, 148, 131, 129, 151, 180,  55,  54,  54,  54,  67, 203, 203,
        189,  78,  78, 119, 102,  97,  90,  78,  78,  90, 180, 201, 203, 203, 203, 203, 203, 118,  54,  54,  54,  54, 105, 203, 203,
        180,  78,  78, 123, 126, 107,  78,  78,  78,  97, 203, 203, 203, 203, 203, 203, 203,  75,  54,  54,  54,  54, 137, 203, 203,
        178,  78,  78, 126, 126, 120,  79,  78,  78,  96, 203, 203, 203, 203, 203, 203, 203,  54,  54,  54,  54,  79, 194, 203, 203,
        179,  78,  78, 126, 126, 126,  96,  78,  78,  86, 203, 203, 203, 203, 203, 203, 203,  74,  54,  54,  87, 190, 203, 203, 203,
        183,  78,  78, 126, 126, 126, 114,  78,  78,  78, 199, 203, 203, 203, 203, 203, 203, 173,  86, 103, 196, 203, 203, 203, 203,
        178,  78,  78, 123, 126, 126, 126,  85,  78,  78, 170, 203, 203, 185, 165, 165, 170, 190, 203, 203, 203, 203, 203, 203, 203,
        174,  78,  78, 116, 126, 126, 126,  92,  78,  78, 142, 203, 198, 130, 126, 126, 126, 126, 150, 193, 203, 203, 203, 203, 203,
        170,  78,  78, 101, 126, 126, 126,  84,  78,  78, 115, 203, 170, 126, 126, 126, 126, 126, 122, 106, 153, 193, 203, 203, 203,
        165,  78,  78,  79, 105, 125, 115, 112, 157, 139,  84, 195, 161, 126, 126, 126, 126, 126,  93,  72, 109, 129, 169, 203, 203,
        165,  78,  78,  78,  78,  80, 112, 161, 161, 161, 123, 160, 165, 126, 126, 126, 126, 119,  73,  72,  74, 116, 126, 168, 203,
        165,  78,  78,  78,  78,  84, 154, 121,  89, 150, 157, 130, 190, 126, 126, 126, 126,  86,  72,  72,  72,  96, 126, 132, 199,
        165,  78,  78,  78,  78, 106, 161,  72,  48,  90, 161, 104, 190, 142, 126, 126, 103,  72,  72,  72,  72,  83, 126, 126, 192,
        165,  78,  78,  78,  78, 114, 161,  58,  46,  89, 161, 107, 133, 195, 136, 119,  75,  72,  72,  72,  72,  72, 114, 134, 202,
        165,  78,  78,  78,  78, 115, 161, 113,  87, 144, 161, 106,  80, 192, 192, 135, 118, 118, 118, 120, 121, 121, 125, 180, 203,
        179,  78,  78,  78,  78,  94, 161, 161, 161, 161, 161,  94,  78, 143, 203, 200, 156, 126, 126, 126, 126, 126, 173, 203, 203,
        203,  98,  78,  78,  78,  78, 137, 161, 161, 161, 146,  79,  78,  83, 197, 203, 203, 191, 182, 177, 171, 175, 203, 203, 203,
        203, 170, 113,  95,  79,  78,  85, 133, 140, 130,  86,  78,  78,  83, 197, 203, 203, 203, 203, 203, 203, 203, 203, 203, 203,
        203, 203, 203, 203, 198, 184, 177, 164, 159, 154, 153, 143, 141, 181, 203, 203, 203, 203, 203, 203, 203, 203, 203, 203, 203
    };
    
    int* img_pointer = data;
    
    
    /*
    const int numRows = 500;
    const int numCols = 500;
    int* img_pointer = readImageFromFile("/Users/wonderalexandre/GitHub/ComponentTreeLearn/dat/imgTeste.txt", numRows, numCols);
    */

    printImage(img_pointer, numRows, numCols);
    std::cout << "img_pointer ok" << std::endl;
    
    // Criar um ComponentTree
    ComponentTree* tree = new ComponentTree(img_pointer, numRows, numCols, true);
    std::cout << "tree ok" << std::endl;
    

    // Criar um AttributeComputedIncrementally::computerArea
    const int n = tree->getNumNodes();	
    int indexAttr = 6; //height
    float* attr = AttributeComputedIncrementally::computerAttribute(tree, indexAttr); //size: n * numAttribute
    std::cout << "attributes ok" << std::endl;

    // Criar um AttributeOpeningPrimitivesFamily
    int maxCriterion = numRows; 
    AttributeOpeningPrimitivesFamily* primitives = new AttributeOpeningPrimitivesFamily(tree, attr, maxCriterion);
    std::cout << "primitives ok" << std::endl;
    
    ResidualTree* residualTree = new ResidualTree(primitives);
    std::cout << "residualTree ok" << std::endl;


    int* imgPos = residualTree->getPositiveResidues();
    printImage(imgPos, numRows, numCols);
    std::cout << "imgPos ok" << std::endl;
    
    int* imgNeg = residualTree->getNegativeResidues();
    printImage(imgNeg, numRows, numCols);
    std::cout << "imgNeg ok" << std::endl;
    

    int* imgRec = residualTree->reconstruction();
    printImage(imgRec, numRows, numCols);

    int* contrast = residualTree->getMaxConstrastImage();
    //printImage(contrast, numRows, numCols);
    std::cout << "contrast ok" << std::endl;
    
     std::vector<float> attrVector(attr, attr + n);
    UltimateAttributeOpening *uao = new UltimateAttributeOpening(tree, attrVector);
    uao->execute(maxCriterion);
    int* contrastUAO = residualTree->getMaxConstrastImage();
    //printImage(contrastUAO, numRows, numCols);

    bool isEquals = true;
    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            int pos = i * numCols + j;
            if(contrast[pos] != contrastUAO[pos]){
                isEquals = false;
                break;
            }
        }
    }
    printf("isEquals: %d\n", isEquals);



    

    return 0;
}
