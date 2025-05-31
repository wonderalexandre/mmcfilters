#include <algorithm>
#include <climits>
#include <vector>
#include <utility>
#include <array>
#include <list>
#include <type_traits>
#include <numeric>
#include <iostream>

#include "../include/AdjacencyRelation.hpp"
#include "../include/BuilderComponentTreeByUnionFind.hpp"


template <typename PixelType>
BuilderComponentTreeByUnionFind<PixelType>::BuilderComponentTreeByUnionFind(ImagePtr<PixelType> img, bool isMaxtree, AdjacencyRelationPtr adj){
    this->sort(img, isMaxtree);
    this->createTreeByUnionFind(img, isMaxtree, adj);
}

template <typename PixelType>
BuilderComponentTreeByUnionFind<PixelType>::~BuilderComponentTreeByUnionFind() {
    //delete[] orderedPixels;
    //delete[] parent;
}

template <typename PixelType>
int* BuilderComponentTreeByUnionFind<PixelType>::getParent(){
	return this->parent.get();
}

template <typename PixelType>
int* BuilderComponentTreeByUnionFind<PixelType>::getOrderedPixels(){
	return this->orderedPixels.get();
}

template <typename PixelType>
void BuilderComponentTreeByUnionFind<PixelType>::sort(ImagePtr<PixelType> imgPtr, bool isMaxtree){
	const int n = imgPtr->getSize();
	this->orderedPixels = std::make_unique<int[]>(n);
	PixelType* img = imgPtr->rawData();
	if constexpr (std::is_floating_point<PixelType>::value){
		if(PRINT_LOG) std::cout << "Sorting floating point image with size: " << n << std::endl;
		int* orderedPixelsPtr = this->orderedPixels.get();
        std::iota(orderedPixelsPtr, orderedPixelsPtr + n, 0);
        if (isMaxtree) {
			std::sort(orderedPixelsPtr, orderedPixelsPtr + n,
				[&](int a, int b) { return img[a] > img[b]; });
		} else {
			std::sort(orderedPixelsPtr, orderedPixelsPtr + n,
				[&](int a, int b) { return img[a] < img[b]; });
		}	
	}else{
		if(PRINT_LOG) std::cout << "Sorting integer image with size: " << n << std::endl;
		int maxvalue =  img[0];
		for (int i = 1; i < n; i++)
			if(maxvalue < img[i]) maxvalue = img[i];
				
		std::unique_ptr<int[]> counter(new int[maxvalue + 1]());;	
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
	
}

template <typename PixelType>
int BuilderComponentTreeByUnionFind<PixelType>::findRoot(int* zPar, int x) {
	if (zPar[x] == x)
		return x;
	else {
		zPar[x] = findRoot(zPar, zPar[x]);
		return zPar[x];
	}
}

template <typename PixelType>
void BuilderComponentTreeByUnionFind<PixelType>::createTreeByUnionFind(ImagePtr<PixelType> imgPtr, bool isMaxtree, AdjacencyRelationPtr adj) {
	if constexpr (std::is_floating_point<PixelType>::value){
		std::cout << "Creating tree by union-find for floating point image." << std::endl;
	} else {
		std::cout << "Creating tree by union-find for integer image." << std::endl;
	}


	const int n = imgPtr->getSize();
	auto img = imgPtr->rawData();
	
	std::unique_ptr<int[]> zParPtr = std::make_unique<int[]>(n);
	int* zPar = zParPtr.get();
	this->parent = std::make_unique<int[]>(n);
	float epsilon = 1e-5f; // Tolerância para comparação de pixels flutuantes
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
				
		if constexpr (std::is_floating_point<PixelType>::value) {
			if (std::fabs(img[parent[q]] - img[q]) <= epsilon) {
				parent[p] = parent[q];
			}
		} else {
			if (img[parent[q]] == img[q]) {
				parent[p] = parent[q];
			}
		}


	}
		
}

