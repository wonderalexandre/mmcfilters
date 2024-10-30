#include <algorithm>
#include <climits>
#include <vector>
#include <utility>
#include <array>
#include <list>

#include "../include/AdjacencyRelation.hpp"


#ifndef BUILDER_COMPONENT_TREE_BY_UNION_FIND_H
#define BUILDER_COMPONENT_TREE_BY_UNION_FIND_H

class BuilderComponentTreeByUnionFind {
private:
	int* parent;
	int *orderedPixels;
    
public:
    
    void sort(int* img, int numRows, int numCols, bool isMaxtree);
	void createTreeByUnionFind(int* img, int numRows, int numCols, bool isMaxtree, AdjacencyRelation* adj);
	int findRoot(int *zPar, int x);
    int* getParent();
    int* getOrderedPixels();
    BuilderComponentTreeByUnionFind(int* img, int numRows, int numCols, bool isMaxtree, AdjacencyRelation* adj);
    ~BuilderComponentTreeByUnionFind();
};

#endif
