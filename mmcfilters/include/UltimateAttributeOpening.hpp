#include <array>
#include "../include/NodeMT.hpp"
#include "../include/MorphologicalTree.hpp"
#include "../include/ComputerMSER.hpp"
#include "../include/AttributeComputedIncrementally.hpp"
#include "../include/Common.hpp"

#ifndef ULTIMATE_ATTR_OPENING_H
#define ULTIMATE_ATTR_OPENING_H


class UltimateAttributeOpening{

  protected:
    int maxCriterion;
    std::vector<float> attrs_increasing;
    MorphologicalTreePtr tree;
    int* maxContrastLUT;
    int* associatedIndexLUT;
    
    void computeUAO(NodeMTPtr currentNode, int levelNodeNotInNR, bool qPropag, bool isCalculateResidue);
    void execute(int maxCriterion, std::vector<bool> selectedForFiltering);
    
    bool isSelectedForPruning(NodeMTPtr currentNode); //first Node in Nr(i)
    bool hasNodeSelectedInPrimitive(NodeMTPtr currentNode); //has node selected inside Nr(i)
    std::vector<bool> selectedForFiltering; //mappping between nodes and selected nodes
   
    

  public:

    UltimateAttributeOpening(MorphologicalTreePtr tree,  std::vector<float> attrs_increasing);

    ~UltimateAttributeOpening();

    int* getMaxConstrastImage();

    int* getAssociatedImage();

    int* getAssociatedColorImage();    

    void execute(int maxCriterion);
    
    void executeWithMSER(int maxCriterion, int deltaMSER);
    
};

#endif





	

