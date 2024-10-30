#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"
#include "../include/ComputerMSER.hpp"
#include "../include/AttributeComputedIncrementally.hpp"

#include <vector>
#include <list>

#ifndef ATTRIBUTE_OPENING_PRIMITIVES_FAMILY_H
#define ATTRIBUTE_OPENING_PRIMITIVES_FAMILY_H

class AttributeOpeningPrimitivesFamily{
  
  protected:
    float* attrs_increasing;
    float maxCriterion;
    std::list<float> thresholds;
    std::list<NodeCT*> nodesWithMaximumCriterium;

    //PrimitivesFamily
    ComponentTree* tree;
    std::vector<bool> selectedForFiltering; //mappping between index nodes and selected nodes
    int* restOfImage;
    int numPrimitives;
    

    void initializeRestOfImage(float threshold);
    void initializeNodesWithMaximumCriterium();
    
  public:
    AttributeOpeningPrimitivesFamily(ComponentTree* tree,  float* attr, float maxCriterion);

    AttributeOpeningPrimitivesFamily(ComponentTree* tree,  float* attrs_increasing, float maxCriterion, int deltaMSER);
    
    ~AttributeOpeningPrimitivesFamily();

    std::list<float> getThresholdsPrimitive();

    //PrimitivesFamily
    bool isSelectedForPruning(NodeCT* node) ; //first Node in Nr(i)

    bool hasNodeSelectedInPrimitive(NodeCT* node) ; //has node selected inside Nr(i)

    std::list<NodeCT*> getNodesWithMaximumCriterium() ; 

    ComponentTree* getTree() ;

    int* getRestOfImage() ;

    int getNumPrimitives() ;
    

};

#endif





