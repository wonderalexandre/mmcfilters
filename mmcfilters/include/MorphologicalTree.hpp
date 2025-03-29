#include <list>
#include <vector>

#include "../include/NodeMT.hpp"
#include "../include/AdjacencyRelation.hpp"
#include "../include/Common.hpp"

#ifndef COMPONENT_TREE_H
#define COMPONENT_TREE_H


class MorphologicalTree {

protected:
	int numCols;
	int numRows;
	int treeType; //0-mintree, 1-maxtree, 2-tree of shapes
	NodeMTPtr root;
	int numNodes;
	std::list<NodeMTPtr> listNodes;
	std::vector<NodeMTPtr> nodes;
	AdjacencyRelationPtr adj;
	
	void reconstruction(NodeMTPtr node, int* imgOut);

public:
   	static const int MAX_TREE = 0;
	static const int MIN_TREE = 1;
	static const int TREE_OF_SHAPES = 2;

	MorphologicalTree(int* img, int numRows, int numCols, bool isMaxtree, double radiusOfAdjacencyRelation);

	MorphologicalTree(int* img, int numRows, int numCols, bool isMaxtree);

	MorphologicalTree(int* img, int numRows, int numCols);

    ~MorphologicalTree();

	int* getInputImage();
	
	NodeMTPtr getRoot();

	bool isMaxtree();

	int getTreeType();

	NodeMTPtr getSC(int pixel);

	std::list<NodeMTPtr> getListNodes();

	int getNumNodes();

	int getNumRowsOfImage();

	int getNumColsOfImage();

	int* reconstructionImage();

	int* getImageAferPruning(NodeMTPtr node);

	void pruning(NodeMTPtr node);

	bool isAncestor(NodeMTPtr u, NodeMTPtr v);
	
	bool isDescendant(NodeMTPtr u, NodeMTPtr v);

	bool isComparable(NodeMTPtr u, NodeMTPtr v);

	bool isStrictAncestor(NodeMTPtr u, NodeMTPtr v);
	
	bool isStrictDescendant(NodeMTPtr u, NodeMTPtr v);
	
	bool isStrictComparable(NodeMTPtr u, NodeMTPtr v);
	
	
};

#endif