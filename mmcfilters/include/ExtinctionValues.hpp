

#ifndef EXTINCTION_VALUES_H
#define EXTINCTION_VALUES_H

#include "../include/ImageUtils.hpp"
#include "../include/MorphologicalTree.hpp"
#include "../include/NodeMT.hpp"
#include "../include/AdjacencyRelation.hpp"
#include "../include/Common.hpp"

#define PI 3.14159265358979323846


struct RegionalExtremaNode{
    NodeMTPtr leaf;
    NodeMTPtr cutoffNode;
    float extinction;
    
    RegionalExtremaNode(NodeMTPtr leaf, NodeMTPtr cutoffNode, float extinction) : leaf(leaf), cutoffNode(cutoffNode), extinction(extinction) { }
};
using RegionalExtremaNodePtr = std::shared_ptr<RegionalExtremaNode>;

class ExtinctionValues; // Forward declaration 
using ExtinctionValuesPtr = std::shared_ptr<ExtinctionValues>;

class ExtinctionValues{

    private:
        std::vector<RegionalExtremaNodePtr> regionalExtremaNodes;
        MorphologicalTreePtr tree;
        std::shared_ptr<float[]> attribute;

    public:
	    ExtinctionValues(MorphologicalTreePtr tree, std::shared_ptr<float[]> attr);

        ImageFloatPtr saliencyMap(int leafToKeep, bool unweighted = true);    
        ImageUInt8Ptr filtering(int leafToKeep);

        std::vector<RegionalExtremaNodePtr>& getExtinctionValues() { return regionalExtremaNodes; }
};

#endif // EXTINCTION_VALUES_H