#include <array>
#include "../include/NodeMT.hpp"
#include "../include/MorphologicalTree.hpp"
#include "../include/ComputerMSER.hpp"
#include "../include/AttributeComputedIncrementally.hpp"
#include "../include/Common.hpp"

#ifndef ULTIMATE_ATTR_OPENING_H
#define ULTIMATE_ATTR_OPENING_H


class UltimateAttributeOpening;
using UltimateAttributeOpeningPtr = std::shared_ptr<UltimateAttributeOpening>;
class UltimateAttributeOpening{

  protected:
    int maxCriterion;
    std::shared_ptr<float[]> attrs_increasing;
    MorphologicalTreePtr tree;
    std::shared_ptr<uint8_t[]> maxContrastLUT;
    std::shared_ptr<int[]> associatedIndexLUT;
    
    void computeUAO(NodeMTPtr currentNode, int levelNodeNotInNR, bool qPropag, bool isCalculateResidue);
    void execute(int maxCriterion, std::vector<bool> selectedForFiltering);
    
    bool isSelectedForPruning(NodeMTPtr currentNode); //first Node in Nr(i)
    bool hasNodeSelectedInPrimitive(NodeMTPtr currentNode); //has node selected inside Nr(i)
    std::vector<bool> selectedForFiltering; //mappping between nodes and selected nodes
   
    

  public:

    UltimateAttributeOpening(MorphologicalTreePtr tree,  std::shared_ptr<float[]> attrs_increasing);

    ~UltimateAttributeOpening();

    ImageUInt8Ptr getMaxConstrastImage();

    ImageInt32Ptr getAssociatedImage();

    ImageUInt8Ptr getAssociatedColorImage();    

    void execute(int maxCriterion);
    
    void executeWithMSER(int maxCriterion, int deltaMSER);
    
};

#endif





	

