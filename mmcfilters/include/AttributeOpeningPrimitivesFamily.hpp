#pragma once


#include "../include/MorphologicalTree.hpp"
#include "../include/ComputerMSER.hpp"
#include "../include/AttributeComputedIncrementally.hpp"
#include "../include/Common.hpp"

#include <vector>
#include <list>

class AttributeOpeningPrimitivesFamily;
using AttributeOpeningPrimitivesFamilyPtr = std::shared_ptr<AttributeOpeningPrimitivesFamily>;

class AttributeOpeningPrimitivesFamily{
  
  protected:
    std::shared_ptr<float[]> attrs_increasing;
    float maxCriterion;
    std::list<float> thresholds;
    std::list<NodeId> nodesWithMaximumCriterium;

    //PrimitivesFamily
    MorphologicalTree* tree;
    std::vector<uint8_t> selectedForFiltering; //mappping between index nodes and selected nodes
    ImageUInt8Ptr restOfImage;
    int numPrimitives;
    

    void initializeRestOfImage(float threshold);
    void initializeNodesWithMaximumCriterium();
    
  public:
    AttributeOpeningPrimitivesFamily(ComponentTreePtr tree,   std::shared_ptr<float[]> attr, float maxCriterion): AttributeOpeningPrimitivesFamily(tree.get(), attr, maxCriterion) {}
    AttributeOpeningPrimitivesFamily(MorphologicalTree* tree,   std::shared_ptr<float[]> attr, float maxCriterion);

    AttributeOpeningPrimitivesFamily(ComponentTreePtr tree,   std::shared_ptr<float[]> attrs_increasing, float maxCriterion, int deltaMSER): AttributeOpeningPrimitivesFamily(tree.get(), attrs_increasing, maxCriterion, deltaMSER) {}
    AttributeOpeningPrimitivesFamily(MorphologicalTree* tree,   std::shared_ptr<float[]> attrs_increasing, float maxCriterion, int deltaMSER);

    ~AttributeOpeningPrimitivesFamily();

    std::list<float> getThresholdsPrimitive();

    //PrimitivesFamily
    bool isSelectedForPruning(NodeId node) ; //first Node in Nr(i)

    bool hasNodeSelectedInPrimitive(NodeId node) ; //has node selected inside Nr(i)

    std::list<NodeId> getNodesWithMaximumCriterium() ; 

    ImageUInt8Ptr getRestOfImage() ;

    int getNumPrimitives() ;
    
    MorphologicalTree* getTree() ;
};





