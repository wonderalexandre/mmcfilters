#include <list>

#include "../include/NodeMT.hpp"
#include "../include/NodeRes.hpp"
#include "../include/AttributeOpeningPrimitivesFamily.hpp"
#include "../include/Common.hpp"

#ifndef RESIDUAL_TREE_H
#define RESIDUAL_TREE_H


class ResidualTree{

    protected:
      NodeRes* root;
      AttributeOpeningPrimitivesFamily* primitivesFamily;
      MorphologicalTreePtr tree;
      PixelValueType* maxContrastLUT;
      int* associatedIndexesLUT;
      int numNodes;
      ImagePtr restOfImage;
      //std::list<NodeRes*> listNodes;
      NodeRes** nodes;

    public:
        ResidualTree(AttributeOpeningPrimitivesFamily* primitivesFamily);

        //void computerNodeRes(NodeCT *currentNode);

        void computerMaximumResidues();

        void createTree();

        ImagePtr reconstruction();

        ~ResidualTree();

        //std::list<NodeRes*> getListNodes();

        NodeRes* getRoot();

        NodeRes* getNodeRes(NodeMTPtr node);

        ImagePtr getMaxConstrastImage();

        ImagePtr filtering(std::vector<bool> criterion);

        int* getAssociatedImage();

        int* getAssociatedColorImage();   

        ImagePtr getRestOfImage();

        ImagePtr getPositiveResidues();

        ImagePtr getNegativeResidues();

        MorphologicalTreePtr getCTree();

};


#endif