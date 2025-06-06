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

template <typename PixelType>
class BuilderComponentTreeByUnionFind {
private:
	
    std::unique_ptr<int[]> parent; 
	std::unique_ptr<int[]> orderedPixels; 
    
public:

    void sort(ImagePtr<PixelType> img, bool isMaxtree);

	void createTreeByUnionFind(ImagePtr<PixelType> img, bool isMaxtree, AdjacencyRelationPtr adj);

	int findRoot(int* zPar, int x);
    int* getParent();
    int* getOrderedPixels();

    BuilderComponentTreeByUnionFind(ImagePtr<PixelType> img, bool isMaxtree, AdjacencyRelationPtr adj);

    ~BuilderComponentTreeByUnionFind();
};

#include "BuilderComponentTreeByUnionFind.tpp"

#endif
