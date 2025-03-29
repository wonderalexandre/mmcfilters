#include <list>

#include "../include/NodeMT.hpp"
#include "../include/NodeRes.hpp"
#include "../include/AttributeOpeningPrimitivesFamily.hpp"

#ifndef RESIDUAL_TREE_H
#define RESIDUAL_TREE_H


class ResidualTree{

    protected:
      NodeRes* root;
      AttributeOpeningPrimitivesFamily* primitivesFamily;
      MorphologicalTreePtr tree;
      int* maxContrastLUT;
      int* associatedIndexesLUT;
      int numNodes;
      int* restOfImage;
      //std::list<NodeRes*> listNodes;
      NodeRes** nodes;

    public:
        ResidualTree(AttributeOpeningPrimitivesFamily* primitivesFamily);

        //void computerNodeRes(NodeCT *currentNode);

        void computerMaximumResidues();

        void createTree();

        int* reconstruction();

        ~ResidualTree();

        //std::list<NodeRes*> getListNodes();

        NodeRes* getRoot();

        NodeRes* getNodeRes(NodeMTPtr node);

        int* getMaxConstrastImage();

        int* filtering(std::vector<bool> criterion, int* imgOutput);

        int* getAssociatedImage();

        int* getAssociatedColorImage();   

        int* getRestOfImage();

        int* getPositiveResidues();

        int* getNegativeResidues();

        MorphologicalTreePtr getCTree();

};


#endif