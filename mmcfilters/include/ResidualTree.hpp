#include <list>

#include "../include/NodeMT.hpp"
#include "../include/NodeRes.hpp"
#include "../include/AttributeOpeningPrimitivesFamily.hpp"
#include "../include/Common.hpp"

#ifndef RESIDUAL_TREE_H
#define RESIDUAL_TREE_H


class ResidualTree{

    protected:
      NodeResPtr root;
      AttributeOpeningPrimitivesFamily* primitivesFamily;
      MorphologicalTreePtr tree;
      ImagePtr maxContrastLUT;
      int* associatedIndexesLUT;
      int numNodes;
      ImagePtr restOfImage;
      std::vector<NodeResPtr> nodes;

    public:
        ResidualTree(AttributeOpeningPrimitivesFamily* primitivesFamily);

        //void computerNodeRes(NodeCT *currentNode);

        void computerMaximumResidues();

        void createTree();

        ImagePtr reconstruction();

        ~ResidualTree();

        //std::list<NodeRes*> getListNodes();

        NodeResPtr getRoot();

        NodeResPtr getNodeRes(NodeMTPtr node);

        ImagePtr getMaxConstrastImage();

        ImagePtr filtering(std::vector<bool> criterion);

        int* getAssociatedImage();

        PixelType* getAssociatedColorImage();   

        ImagePtr getRestOfImage();

        ImagePtr getPositiveResidues();

        ImagePtr getNegativeResidues();

        MorphologicalTreePtr getCTree();

};


#endif