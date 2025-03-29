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
	//nodes = nullptr;
 }

MorphologicalTree::MorphologicalTree(int* img, int numRows, int numCols, bool isMaxtree) 
	: MorphologicalTree(img, numRows, numCols, isMaxtree, 1.5){ }


MorphologicalTree::MorphologicalTree(int* img, int numRows, int numCols){
	this->numRows = numRows;
	this->numCols = numCols;
	this->treeType = TREE_OF_SHAPES;
	this->nodes.resize(this->numRows * this->numCols, nullptr);

	BuilderTreeOfShapeByUnionFind* builder = new BuilderTreeOfShapeByUnionFind();
	builder->interpolateImage(img, numRows, numCols);
	int* interpolationMin = builder->getInterpolationMin();
	int* interpolationMax = builder->getInterpolationMax();
	
	builder->sort();
	int* imgR = builder->getImgR();
	int* imgU = builder->getImgU();
	
	builder->createTreeByUnionFind();
	int* parent = builder->getParent();
	
	int size = builder->getInterpNumCols() * builder->getInterpNumRows();
    std::vector<NodeMTPtr> nodes(size);
    
	
    this->numNodes = 0;
    for (int i = 0; i < size; i++) {
		int p = imgR[i];
        auto [px, py] = ImageUtils::to2D(p, builder->getInterpNumCols());
		int pixelUnterpolate = (px/2) + (py/2) * numCols;
            
		if (p == parent[p]) { //representante do node raiz
            this->root = nodes[p] = std::make_shared<NodeMT>(this->numNodes, pixelUnterpolate, nullptr, imgU[p]);
		}
		else if (imgU[p] != imgU[parent[p]]) { //representante de um node
			nodes[p] = std::make_shared<NodeMT>(this->numNodes, pixelUnterpolate, nodes[parent[p]], imgU[p]);
			nodes[parent[p]]->addChild(nodes[p]);
		}
		else if (imgU[p] == imgU[parent[p]]) {
			nodes[p] = nodes[parent[p]];
		}

		if(px % 2 == 1 && py % 2 == 1){
			nodes[p]->addCNPs(pixelUnterpolate);
			this->nodes[pixelUnterpolate] = nodes[p];					
		}
	}
	if(this->root->getCNPs().size() == 0){
		this->root->setResidue(0);
	}
	AttributeComputedIncrementally::computerAttribute(this->root,
		[this](NodeMTPtr _node) -> void { //pre-processing
			_node->setAreaCC( _node->getCNPs().size() );
			_node->setNumDescendants( _node->getChildren().size() );
			_node->setIndex(this->numNodes++);
			this->listNodes.push_back(_node);
		},
		[](NodeMTPtr _root, NodeMTPtr _child) -> void { //merge-processing
			_root->setAreaCC( _root->getAreaCC() + _child->getAreaCC() );
			_root->setNumDescendants( _root->getNumDescendants() + _child->getNumDescendants() );
		},
		[](NodeMTPtr node) -> void { //post-processing
									
		}
	);
	
	delete builder;
	imgR = nullptr;
	imgU = nullptr;
	parent = nullptr;
	interpolationMin = nullptr;
	interpolationMax = nullptr;
} 

 
MorphologicalTree::MorphologicalTree(int* img, int numRows, int numCols, bool isMaxtree, double radiusOfAdjacencyRelation){
	this->numRows = numRows;
	this->numCols = numCols;
	this->treeType = isMaxtree? MAX_TREE : MIN_TREE;

	this->adj = std::make_shared<AdjacencyRelation>(numRows, numCols, radiusOfAdjacencyRelation);	
	BuilderComponentTreeByUnionFind* builder = new BuilderComponentTreeByUnionFind(img, numRows, numCols, isMaxtree, adj);
	
	int n = this->numRows * this->numCols;
	int* orderedPixels = builder->getOrderedPixels();
	int* parent = builder->getParent();
		
	this->nodes.resize(n, nullptr);

	this->numNodes = 0;
	for (int i = 0; i < n; i++) {
		int p = orderedPixels[i];
		if (p == parent[p]) { //representante do node raiz
			this->root = this->nodes[p] = std::make_shared<NodeMT>(this->numNodes++, p, nullptr, img[p]);
			this->listNodes.push_back(this->nodes[p]);
			this->nodes[p]->addCNPs(p);
		}
		else if (img[p] != img[parent[p]]) { //representante de um node
			this->nodes[p] = std::make_shared<NodeMT>(this->numNodes++, p, this->nodes[parent[p]], img[p]);
			this->listNodes.push_back(this->nodes[p]);
			this->nodes[p]->addCNPs(p);
			this->nodes[parent[p]]->addChild(this->nodes[p]);
		}
		else if (img[p] == img[parent[p]]) {
			this->nodes[parent[p]]->addCNPs(p);
			this->nodes[p] = this->nodes[parent[p]];
		}
	}
	
	int timer = 0;
	AttributeComputedIncrementally::computerAttribute(this->root,
		[&timer](NodeMTPtr _node) -> void { // pre-processing
			_node->setAreaCC(_node->getCNPs().size());
			_node->setNumDescendants(_node->getChildren().size());
			_node->setTimePreOrder(timer++);
		},
		[](NodeMTPtr _root, NodeMTPtr _child) -> void { // merge-processing
			_root->setAreaCC(_root->getAreaCC() + _child->getAreaCC());
			_root->setNumDescendants(_root->getNumDescendants() + _child->getNumDescendants());
		},
		[&timer](NodeMTPtr node) -> void { // post-processing
			node->setTimePostOrder(timer++);
		}
	);
	
	delete builder;
	builder = nullptr;
	orderedPixels = nullptr;
	parent = nullptr;
}

NodeMTPtr MorphologicalTree::getSC(int pixel){
	return this->nodes[pixel];
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

std::list<NodeMTPtr> MorphologicalTree::getListNodes(){
	return this->listNodes;
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





int* MorphologicalTree::getImageAferPruning(NodeMTPtr nodePruning){
	int n = this->numRows * this->numCols;
	int* imgOut = new int[n];
	std::stack<NodeMTPtr> s;
	s.push(this->root);
	while(!s.empty()){
		NodeMTPtr node = s.top();s.pop();
		if(node->getIndex() == nodePruning->getIndex()){
			for(int p: node->getPixelsOfCC()){
				if(node->getParent() != nullptr)
					imgOut[p] = node->getParent()->getLevel();
				else
					imgOut[p] = node->getLevel();
			}
		}
		else{
			for(int p: node->getCNPs()){
				imgOut[p] = node->getLevel();
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
		int numDescendants = nodePruning->getParent()->getNumDescendants();
		int numDescendantsChild = nodePruning->getNumDescendants() + 1;
		nodePruning->getParent()->setNumDescendants(numDescendants - numDescendantsChild); 
		nodePruning->getParent()->getChildren().remove(nodePruning);
		nodePruning->setParent(nullptr);
		nodePruning = nullptr;
		
	}
}

int* MorphologicalTree::reconstructionImage(){
	int n = this->numRows * this->numCols;
	int *imgOut = new int[n];
	this->reconstruction(this->root, imgOut);
	return imgOut;
}

int* MorphologicalTree::getInputImage(){
	int n = this->numRows * this->numCols;
	int* img = new int[n];
	this->reconstruction(this->root, img);
	return img;
}
	

void MorphologicalTree::reconstruction(NodeMTPtr node, int* imgOut){
	for (int p : node->getCNPs()){
		imgOut[p] = node->getLevel();
	}
	for(NodeMTPtr child: node->getChildren()){
		reconstruction(child, imgOut);
	}
}