#include <list>
#include <vector>

#include "../include/NodeCT.hpp"
#include "../include/AdjacencyRelation.hpp"

#ifndef COMPONENT_TREE_H
#define COMPONENT_TREE_H


class ComponentTree {

protected:
	int numCols;
	int numRows;
	int treeType; //0-mintree, 1-maxtree, 2-tree of shapes
	NodeCT* root;
	int numNodes;
	std::list<NodeCT*> listNodes;
	NodeCT** nodes;

	void reconstruction(NodeCT* node, int* imgOut);

public:
   	static const int MAX_TREE = 0;
	static const int MIN_TREE = 1;
	static const int TREE_OF_SHAPES = 2;

	ComponentTree(int* img, int numRows, int numCols, bool isMaxtree, double radiusOfAdjacencyRelation);

	ComponentTree(int* img, int numRows, int numCols, bool isMaxtree);

	ComponentTree(int* img, int numRows, int numCols);

    ~ComponentTree();

	int* getInputImage();
	
	NodeCT* getRoot();

	bool isMaxtree();

	int getTreeType();

	NodeCT* getSC(int pixel);

	std::list<NodeCT*> getListNodes();

	int getNumNodes();

	int getNumRowsOfImage();

	int getNumColsOfImage();

	int* reconstructionImage();

	int* getImageAferPruning(NodeCT* node);

	void pruning(NodeCT* node);
	
};

#endif