#pragma once

#include "../trees/MorphologicalTree.hpp"
#include "../trees/NodeMT.hpp"
#include "../attributes/AttributeComputedIncrementally.hpp"

#include <vector>
#include <tuple>
namespace mmcfilters {


/**
 * @brief Calcula derivadas auxiliares para pipelines de aprendizado sobre a árvore.
 */
class ComputerDerivatives {
    
    private:
        

    public:

        static std::tuple<float*, float*, float*> gradients(MorphologicalTreePtr tree, float* attributes, std::vector<float>& sigmoid, float* gradLoss) {
            int rows = tree->getNumRowsOfImage();
            int cols = tree->getNumColsOfImage();
            
            float *dWeight = new float[rows * cols];
            float *dBias = new float[rows];
            for(NodeMTPtr node: tree->getListNodes()){
                int id = node->getIndex();
                dBias[id] = (sigmoid[id] * (1 - sigmoid[id])) * node->getResidue();
                for (int j = 0; j < cols; j++)
                    dWeight[id + (rows * j)] = (sigmoid[id] * (1 - sigmoid[id])) * attributes[id + (rows * j)] * node->getResidue();

                if (id > 0) {
                    int idParent = node->getParent()->getIndex();
                    dBias[id] += dBias[idParent];
                    for (int j = 0; j < cols; j++)
                        dWeight[id + (rows * j)] += dWeight[idParent + rows * j];
                }
            }

            float *gradWeight = new float[cols];
            float *gradBias = new float[1];
            gradBias[0] = 0;
            for (int j = 0; j < cols; j++)
                gradWeight[j] = 0;
            

            std::unique_ptr<float[]> summationGrad(new float[tree->getNumNodes()]);
            float *gradInput = new float[tree->getNumRowsOfImage() * tree->getNumColsOfImage()];
            AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
                    [&summationGrad, &gradLoss, &sigmoid](NodeMTPtr node) -> void { //pre-processing
                        summationGrad[node->getIndex()] = 0;
                        for(int p: node->getCNPs()){
                            summationGrad[node->getIndex()] += gradLoss[p];
                        }
                        summationGrad[node->getIndex()] = summationGrad[node->getIndex()] * sigmoid[node->getIndex()]; 
                    },
                    [&summationGrad](NodeMTPtr parent, NodeMTPtr child) -> void { //merge-processing
                        summationGrad[parent->getIndex()] += summationGrad[child->getIndex()];
                    },
                    [&summationGrad, &gradInput, &gradWeight,&gradBias, &dBias, &dWeight, &rows, &cols, &gradLoss ](NodeMTPtr node) -> void { //post-processing	
                        for(int p: node->getCNPs()){
                            gradInput[p] = summationGrad[node->getIndex()]; 
                            gradBias[0] += gradLoss[p] * dBias[node->getIndex()];
                            for (int j = 0; j < cols; j++)
                                gradWeight[j] += gradLoss[p] * dWeight[node->getIndex() + (rows * j)];
                            }		
                    }
            );
            
            return {gradWeight, gradBias, gradInput};
        }

};

} // namespace mmcfilters

