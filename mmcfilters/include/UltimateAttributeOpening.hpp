#include <array>
#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"
#include "../include/ComputerMSER.hpp"
#include "../include/AttributeComputedIncrementally.hpp"

#ifndef ULTIMATE_ATTR_OPENING_H
#define ULTIMATE_ATTR_OPENING_H


class UltimateAttributeOpening{

  protected:
    int maxCriterion;
    std::vector<float> attrs_increasing;
    ComponentTree* tree;
    int* maxContrastLUT;
    int* associatedIndexLUT;
    
    void computeUAO(NodeCT* currentNode, int levelNodeNotInNR, bool qPropag, bool isCalculateResidue);
    void execute(int maxCriterion, std::vector<bool> selectedForFiltering);
    
    bool isSelectedForPruning(NodeCT* currentNode); //first Node in Nr(i)
    bool hasNodeSelectedInPrimitive(NodeCT* currentNode); //has node selected inside Nr(i)
    std::vector<bool> selectedForFiltering; //mappping between nodes and selected nodes
   
    

  public:

    UltimateAttributeOpening(ComponentTree* tree,  std::vector<float> attrs_increasing);

    ~UltimateAttributeOpening();

    int* getMaxConstrastImage();

    int* getAssociatedImage();

    int* getAssociatedColorImage();    

    void execute(int maxCriterion);
    
    void executeWithMSER(int maxCriterion, int deltaMSER);
    
};

#endif





	

