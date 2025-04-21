#include <list>
#include <vector>
#include <stack>


#include "../include/NodeMT.hpp"
#include "../include/MorphologicalTree.hpp"
#include "../include/AdjacencyRelation.hpp"
#include "../include/AttributeComputedIncrementally.hpp"
#include "../include/BuilderTreeOfShapeByUnionFind.hpp"
#include "../include/ImageUtils.hpp"
#include "../include/BuilderComponentTreeByUnionFind.hpp"




 MorphologicalTree::~MorphologicalTree(){
	
 }

MorphologicalTree::MorphologicalTree(ImagePtr imgPtr, std::string ToSInperpolation){
	this->treeType = TREE_OF_SHAPES;
	this->numRows = imgPtr->numRows;
	this->numCols = imgPtr->numCols;
	this->nodes.resize(numRows * numCols, nullptr);

	BuilderTreeOfShapeByUnionFind* builder = new BuilderTreeOfShapeByUnionFind();
	if(ToSInperpolation == "4c8c")
		builder->interpolateImage4c8c(imgPtr);
	else
		builder->interpolateImage(imgPtr);
	
	builder->sort();
	int* imgR = builder->getImgR();
	PixelType* imgU = builder->getImgU();
	
	builder->createTreeByUnionFind();
	int* parent = builder->getParent();
	
	int size = builder->getInterpNumCols() * builder->getInterpNumRows();
    std::vector<NodeMTPtr> nodesTmp(size);
    
	
    this->numNodes = 0;
    for (int i = 0; i < size; i++) {
		int p = imgR[i];
        auto [row, col] = ImageUtils::to2D(p, builder->getInterpNumCols());
		    
		if (p == parent[p]) { //representante do node raiz
            this->root = nodesTmp[p] = std::make_shared<NodeMT>(this->numNodes++, nullptr, imgU[p]);
		}
		else if (imgU[p] != imgU[parent[p]]) { //representante de um node
			nodesTmp[p] = std::make_shared<NodeMT>(this->numNodes++, nodesTmp[parent[p]], imgU[p]);
			nodesTmp[parent[p]]->addChild(nodesTmp[p]);
		}
		else if (imgU[p] == imgU[parent[p]]) {
			nodesTmp[p] = nodesTmp[parent[p]];
		}

		if(row % 2 == 1 && col % 2 == 1){
			int pixelUnterpolate =  ImageUtils::to1D(row/2, col/2, numCols);
			this->nodes[pixelUnterpolate] = nodesTmp[p];					
			this->nodes[pixelUnterpolate]->addCNPs(pixelUnterpolate);
		}
	}
	
	if(ToSInperpolation == "4c8c"){
		AttributeComputedIncrementally::computerAttribute(this->root,
			[&](NodeMTPtr node) -> void { //pre-processing
				node->setLevel((imgPtr->data)[node->getCNPs().front()]);
			},
			[](NodeMTPtr parent, NodeMTPtr child) -> void { },
			[](NodeMTPtr node) -> void {}
		);
	}

	computerTreeAttribute();
	delete builder;
	imgR = nullptr;
	imgU = nullptr;
	parent = nullptr;
	
} 

void MorphologicalTree::computerTreeAttribute(){
	this->indexToNode.resize(this->numNodes);
	this->numNodes = 0;
	AttributeComputedIncrementally::computerAttribute(this->root,
		[this](NodeMTPtr node) -> void { //pre-processing
			this->indexToNode[this->numNodes] = node;
			node->setIndex(this->numNodes++);
		},
		[](NodeMTPtr parent, NodeMTPtr child) -> void { //merge-processing
			
		},
		[](NodeMTPtr node) -> void { // post-processing
		}
	);

	int timer = 0;
	int maxDepth = 0;
	int* depth = new int[this->numNodes];
	AttributeComputedIncrementally::computerAttribute(this->root,
		[this, &timer, depth](NodeMTPtr node) -> void { //pre-processing
			node->setAreaCC( node->getCNPs().size() );
			node->setTimePreOrder(timer++);
			depth[node->getIndex()] =  node->getParent() == nullptr ? 0 : depth[node->getParent()->getIndex()] + 1;
		},
		[](NodeMTPtr parent, NodeMTPtr child) -> void { //merge-processing
			parent->setAreaCC( parent->getAreaCC() + child->getAreaCC() );
		},
		[&timer, depth, &maxDepth](NodeMTPtr node) -> void { // post-processing
			node->setTimePostOrder(timer++);
			maxDepth = std::max(maxDepth, depth[node->getIndex()]);
		}
	);
	this->depth = maxDepth;
	delete[] depth;
}

 
MorphologicalTree::MorphologicalTree(ImagePtr imgPtr, bool isMaxtree, double radiusOfAdjacencyRelation){
	this->numRows = imgPtr->numRows;
	this->numCols = imgPtr->numCols;
	PixelType* img = imgPtr->rawData();
	
	this->treeType = isMaxtree? MAX_TREE : MIN_TREE;

	this->adj = std::make_shared<AdjacencyRelation>(numRows, numCols, radiusOfAdjacencyRelation);	
	BuilderComponentTreeByUnionFind* builder = new BuilderComponentTreeByUnionFind(imgPtr, isMaxtree, adj);
	
	int n = numRows * numCols;
	int* orderedPixels = builder->getOrderedPixels();
	int* parent = builder->getParent();
		
	this->nodes.resize(n, nullptr);

	this->numNodes = 0;
	for (int i = 0; i < n; i++) {
		int p = orderedPixels[i];
		if (p == parent[p]) { //representante do node raiz
			this->root = this->nodes[p] = std::make_shared<NodeMT>(this->numNodes++, nullptr, img[p]);
			this->nodes[p]->addCNPs(p);
		}
		else if (img[p] != img[parent[p]]) { //representante de um node
			this->nodes[p] = std::make_shared<NodeMT>(this->numNodes++, this->nodes[parent[p]], img[p]);
			this->nodes[p]->addCNPs(p);
			this->nodes[parent[p]]->addChild(this->nodes[p]);
		}
		else if (img[p] == img[parent[p]]) {
			this->nodes[parent[p]]->addCNPs(p);
			this->nodes[p] = this->nodes[parent[p]];
		}
	}
	
	
	computerTreeAttribute();
	delete builder;
	builder = nullptr;
	orderedPixels = nullptr;
	parent = nullptr;
}

int MorphologicalTree::getDepth(){
	return this->depth;
}


NodeMTPtr MorphologicalTree::getSC(int pixel){
	return this->nodes[pixel];
}
	
NodeMTPtr MorphologicalTree::getNodeByIndex(int index){
	return this->indexToNode[index];
}

NodeMTPtr MorphologicalTree::getRoot() {
	return this->root;
}

bool MorphologicalTree::isMaxtree(){
	return this->treeType == MAX_TREE;
}

int MorphologicalTree::getTreeType(){
	return this->treeType;
}

std::vector<NodeMTPtr>& MorphologicalTree::getIndexNode(){
	return this->indexToNode;
}

int MorphologicalTree::getNumNodes(){
	return this->numNodes;
}

int MorphologicalTree::getNumRowsOfImage(){
	return this->numRows;
}

int MorphologicalTree::getNumColsOfImage(){
	return this->numCols;
}

bool MorphologicalTree::isAncestor(NodeMTPtr u, NodeMTPtr v) {
    return u->getTimePreOrder() <= v->getTimePreOrder() && u->getTimePostOrder() >= v->getTimePostOrder();
}

bool MorphologicalTree::isDescendant(NodeMTPtr u, NodeMTPtr v) {
    return v->getTimePreOrder() <= u->getTimePreOrder() && v->getTimePostOrder() >= u->getTimePostOrder();
}

bool MorphologicalTree::isComparable(NodeMTPtr u, NodeMTPtr v) {
    return isAncestor(u, v) || isAncestor(v, u);
}

bool MorphologicalTree::isStrictAncestor(NodeMTPtr u, NodeMTPtr v) {
    return u != v &&
           u->getTimePreOrder() <= v->getTimePreOrder() &&
           u->getTimePostOrder() >= v->getTimePostOrder();
}

bool MorphologicalTree::isStrictDescendant(NodeMTPtr u, NodeMTPtr v) {
    return u != v &&
           v->getTimePreOrder() <= u->getTimePreOrder() &&
           v->getTimePostOrder() >= u->getTimePostOrder();
}

bool MorphologicalTree::isStrictComparable(NodeMTPtr u, NodeMTPtr v) {
    return isStrictAncestor(u, v) || isStrictAncestor(v, u);
}

NodeMTPtr MorphologicalTree::findLowestCommonAncestor(NodeMTPtr u, NodeMTPtr v){
    // Troca para garantir que u não é mais profundo que v
    if (u->getTimePreOrder() > v->getTimePreOrder())
        std::swap(u, v);

    while (!isAncestor(u, v)) {
        u = u->getParent();
    }
    return u;
}




ImagePtr MorphologicalTree::getImageAferPruning(NodeMTPtr nodePruning){
	ImagePtr imgOut = Image::create(getNumRowsOfImage(), getNumColsOfImage());
	std::stack<NodeMTPtr> s;
	s.push(this->root);
	while(!s.empty()){
		NodeMTPtr node = s.top();s.pop();
		if(node->getIndex() == nodePruning->getIndex()){
			for(int p: node->getPixelsOfCC()){
				if(node->getParent() != nullptr)
					imgOut->data[p] = node->getParent()->getLevel();
				else
					imgOut->data[p] = node->getLevel();
			}
		}
		else{
			for(int p: node->getCNPs()){
				imgOut->data[p] = node->getLevel();
			}
			for(NodeMTPtr child: node->getChildren()){
				s.push(child);
			}
		}
	}
	return imgOut;
}

void MorphologicalTree::pruning(NodeMTPtr nodePruning){
	if(nodePruning->getParent() != nullptr){
		for(int p: nodePruning->getPixelsOfCC()){
			nodePruning->getParent()->addCNPs(p);
			this->nodes[p] = nodePruning->getParent()->getParent();
		}
		nodePruning->getParent()->getChildren().remove(nodePruning);
		nodePruning->setParent(nullptr);
		nodePruning = nullptr;
		this->computerTreeAttribute();
	}
}

ImagePtr MorphologicalTree::reconstructionImage(){
	ImagePtr imgOut = Image::create(getNumRowsOfImage(), getNumColsOfImage());
	this->reconstruction(this->root, imgOut->rawData());
	return imgOut;
}


void MorphologicalTree::reconstruction(NodeMTPtr node, PixelType* dataOut){
	for (int p : node->getCNPs()){
		dataOut[p] = node->getLevel();
	}
	for(NodeMTPtr child: node->getChildren()){
		reconstruction(child, dataOut);
	}
}

std::vector<std::vector<NodeMTPtr>> MorphologicalTree::getNodesByDepth(){
	std::vector<std::vector<NodeMTPtr>> nodesByDepth(this->depth + 1);
	MorphologicalTree::extractDepthMap(this->root, 0, nodesByDepth);
	return nodesByDepth;
}