#include "../include/NodeMT.hpp"
#include "../include/MorphologicalTree.hpp"
#include "../include/ComputerMSER.hpp"
#include "../include/AttributeComputedIncrementally.hpp"
#include "../include/Common.hpp"

#include <vector>
#include <list>

#ifndef ATTRIBUTE_OPENING_PRIMITIVES_FAMILY_H
#define ATTRIBUTE_OPENING_PRIMITIVES_FAMILY_H

class AttributeOpeningPrimitivesFamily;
using AttributeOpeningPrimitivesFamilyPtr = std::shared_ptr<AttributeOpeningPrimitivesFamily>;

class AttributeOpeningPrimitivesFamily{
  
  protected:
    std::shared_ptr<float[]> attrs_increasing;
    float maxCriterion;
    std::list<float> thresholds;
    std::list<NodeMTPtr> nodesWithMaximumCriterium;

    //PrimitivesFamily
    MorphologicalTreePtr tree;
    std::vector<bool> selectedForFiltering; //mappping between index nodes and selected nodes
    ImageUInt8Ptr restOfImage;
    int numPrimitives;
    

    void initializeRestOfImage(float threshold);
    void initializeNodesWithMaximumCriterium();
    
  public:
    AttributeOpeningPrimitivesFamily(MorphologicalTreePtr tree,   std::shared_ptr<float[]> attr, float maxCriterion);

    AttributeOpeningPrimitivesFamily(MorphologicalTreePtr tree,   std::shared_ptr<float[]> attrs_increasing, float maxCriterion, int deltaMSER);
    
    ~AttributeOpeningPrimitivesFamily();

    std::list<float> getThresholdsPrimitive();

    //PrimitivesFamily
    bool isSelectedForPruning(NodeMTPtr node) ; //first Node in Nr(i)

    bool hasNodeSelectedInPrimitive(NodeMTPtr node) ; //has node selected inside Nr(i)

    std::list<NodeMTPtr> getNodesWithMaximumCriterium() ; 

    MorphologicalTreePtr getTree() ;

    ImageUInt8Ptr getRestOfImage() ;

    int getNumPrimitives() ;
    

};

#endif





