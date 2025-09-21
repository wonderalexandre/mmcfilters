#pragma once


#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"
#include "../include/ComputerMSERCT.hpp"
#include "../include/AttributeComputedIncrementallyCT.hpp"
#include "../include/Common.hpp"

#include <vector>
#include <list>

class AttributeOpeningPrimitivesFamilyCT;
using AttributeOpeningPrimitivesFamilyCTPtr = std::shared_ptr<AttributeOpeningPrimitivesFamilyCT>;

class AttributeOpeningPrimitivesFamilyCT{
  
  protected:
    std::shared_ptr<float[]> attrs_increasing;
    float maxCriterion;
    std::list<float> thresholds;
    std::list<NodeId> nodesWithMaximumCriterium;

    //PrimitivesFamily
    ComponentTree* tree;
    std::vector<uint8_t> selectedForFiltering; //mappping between index nodes and selected nodes
    ImageUInt8Ptr restOfImage;
    int numPrimitives;
    

    void initializeRestOfImage(float threshold);
    void initializeNodesWithMaximumCriterium();
    
  public:
    AttributeOpeningPrimitivesFamilyCT(ComponentTreePtr tree,   std::shared_ptr<float[]> attr, float maxCriterion): AttributeOpeningPrimitivesFamilyCT(tree.get(), attr, maxCriterion) {}
    AttributeOpeningPrimitivesFamilyCT(ComponentTree* tree,   std::shared_ptr<float[]> attr, float maxCriterion);

    AttributeOpeningPrimitivesFamilyCT(ComponentTreePtr tree,   std::shared_ptr<float[]> attrs_increasing, float maxCriterion, int deltaMSER): AttributeOpeningPrimitivesFamilyCT(tree.get(), attrs_increasing, maxCriterion, deltaMSER) {}
    AttributeOpeningPrimitivesFamilyCT(ComponentTree* tree,   std::shared_ptr<float[]> attrs_increasing, float maxCriterion, int deltaMSER);

    ~AttributeOpeningPrimitivesFamilyCT();

    std::list<float> getThresholdsPrimitive();

    //PrimitivesFamily
    bool isSelectedForPruning(NodeId node) ; //first Node in Nr(i)

    bool hasNodeSelectedInPrimitive(NodeId node) ; //has node selected inside Nr(i)

    std::list<NodeId> getNodesWithMaximumCriterium() ; 

    ImageUInt8Ptr getRestOfImage() ;

    int getNumPrimitives() ;
    
    ComponentTree* getTree() ;
};





