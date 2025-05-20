#include <algorithm>
#include <climits>
#include <vector>
#include <utility>
#include <array>
#include <list>

#include "../include/AdjacencyRelation.hpp"
#include "../include/Common.hpp"


#ifndef BUILDER_COMPONENT_TREE_BY_UNION_FIND_H
#define BUILDER_COMPONENT_TREE_BY_UNION_FIND_H

class BuilderComponentTreeByUnionFind {
private:
	int* parent;
	int *orderedPixels;
    
public:
    
    void sort(ImageUInt8Ptr img, bool isMaxtree);
	void createTreeByUnionFind(ImageUInt8Ptr img, bool isMaxtree, AdjacencyRelationPtr adj);
	int findRoot(int *zPar, int x);
    int* getParent();
    int* getOrderedPixels();
    BuilderComponentTreeByUnionFind(ImageUInt8Ptr img, bool isMaxtree, AdjacencyRelationPtr adj);
    ~BuilderComponentTreeByUnionFind();
};

#endif
