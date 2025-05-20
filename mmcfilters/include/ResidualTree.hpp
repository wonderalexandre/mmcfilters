#include <list>

#include "../include/NodeMT.hpp"
#include "../include/NodeRes.hpp"
#include "../include/AttributeOpeningPrimitivesFamily.hpp"
#include "../include/Common.hpp"

#ifndef RESIDUAL_TREE_H
#define RESIDUAL_TREE_H


class ResidualTree;
using ResidualTreePtr = std::shared_ptr<ResidualTree>;
class ResidualTree{

    protected:
      NodeResPtr root;
      std::shared_ptr<AttributeOpeningPrimitivesFamily> primitivesFamily;
      MorphologicalTreePtr tree;
      ImageUInt8Ptr maxContrastLUT;
      std::shared_ptr<int[]> associatedIndexesLUT;
      int numNodes;
      ImageUInt8Ptr restOfImage;
      std::vector<NodeResPtr> nodes;

    public:
        ResidualTree(std::shared_ptr<AttributeOpeningPrimitivesFamily> primitivesFamily);

        //void computerNodeRes(NodeCT *currentNode);

        void computerMaximumResidues();

        void createTree();

        ImageUInt8Ptr reconstruction();

        ~ResidualTree();

        //std::list<NodeRes*> getListNodes();

        NodeResPtr getRoot();

        NodeResPtr getNodeRes(NodeMTPtr node);

        ImageUInt8Ptr getMaxConstrastImage();

        ImageUInt8Ptr filtering(std::vector<bool> criterion);

        ImageInt32Ptr getAssociatedImage();

        ImageUInt8Ptr getAssociatedColorImage();   

        ImageUInt8Ptr getRestOfImage();

        ImageUInt8Ptr getPositiveResidues();

        ImageUInt8Ptr getNegativeResidues();

        MorphologicalTreePtr getCTree();

};


#endif