#include <algorithm>
#include <climits>
#include <vector>
#include <utility>
#include <array>
#include <list>

#include "../include/AdjacencyRelation.hpp"
#include "../include/BuilderComponentTreeByUnionFind.hpp"



BuilderComponentTreeByUnionFind::BuilderComponentTreeByUnionFind(ImagePtr img, bool isMaxtree, AdjacencyRelationPtr adj){
    this->sort(img, isMaxtree);
    this->createTreeByUnionFind(img, isMaxtree, adj);
}

BuilderComponentTreeByUnionFind::~BuilderComponentTreeByUnionFind() {
    delete[] orderedPixels;
    delete[] parent;
}

int* BuilderComponentTreeByUnionFind::getParent(){
	return this->parent;
}

int* BuilderComponentTreeByUnionFind::getOrderedPixels(){
	return this->orderedPixels;
}


void BuilderComponentTreeByUnionFind::sort(ImagePtr imgPtr, bool isMaxtree){
	const int n = imgPtr->numRows * imgPtr->numCols;
	PixelType* img = imgPtr->rawData();
	
	int maxvalue = img[0];
	for (int i = 1; i < n; i++)
		if(maxvalue < img[i]) maxvalue = img[i];
			
	std::vector<int> counter(maxvalue + 1, 0); 
	this->orderedPixels = new int[n];
		
	if(isMaxtree){
		for (int i = 0; i < n; i++)
			counter[img[i]]++;

		for (int i = 1; i < maxvalue; i++) 
			counter[i] += counter[i - 1];
		counter[maxvalue] += counter[maxvalue-1];
		
		for (int i = n - 1; i >= 0; --i)
			this->orderedPixels[--counter[img[i]]] = i;	

	}else{
		for (int i = 0; i < n; i++)
			counter[maxvalue - img[i]]++;

		for (int i = 1; i < maxvalue; i++) 
			counter[i] += counter[i - 1];
		counter[maxvalue] += counter[maxvalue-1];

		for (int i = n - 1; i >= 0; --i)
			this->orderedPixels[--counter[maxvalue - img[i]]] = i;
	}
	
}

int BuilderComponentTreeByUnionFind::findRoot(int *zPar, int x) {
	if (zPar[x] == x)
		return x;
	else {
		zPar[x] = findRoot(zPar, zPar[x]);
		return zPar[x];
	}
}

void BuilderComponentTreeByUnionFind::createTreeByUnionFind(ImagePtr imgPtr, bool isMaxtree, AdjacencyRelationPtr adj) {
	const int n = imgPtr->numRows * imgPtr->numCols;
	PixelType* img = imgPtr->rawData();
	int *zPar = new int[n];
	this->parent = new int[n];
	
	for (int p = 0; p < n; p++) {
		zPar[p] =  -1;
	}
		
	for(int i=n-1; i >= 0; i--){
		int p = orderedPixels[i];
		parent[p] = p;
		zPar[p] = p;
		for (int n : adj->getAdjPixels(p)) {
			if(zPar[n] != -1){
				int r = this->findRoot(zPar, n);
				if(p != r){
					parent[r] = p;
					zPar[r] = p;
				}
			}
		}
	}
			
	// canonizacao da arvore
	for (int i = 0; i < n; i++) {
		int p = orderedPixels[i];
		int q = parent[p];
				
		if(img[parent[q]] == img[q]){
			parent[p] = parent[q];
		}
	}
		
	delete[] zPar;
    zPar = nullptr;
}

