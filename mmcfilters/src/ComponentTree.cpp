#include <list>
#include <vector>
#include <stack>


#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"
#include "../include/AdjacencyRelation.hpp"
#include "../include/AttributeComputedIncrementally.hpp"
#include "../include/BuilderTreeOfShapeByUnionFind.hpp"
#include "../include/ImageUtils.hpp"
#include "../include/BuilderComponentTreeByUnionFind.hpp"


void ComponentTree::reconstruction(NodeCT* node, int* imgOut){
	for (int p : node->getCNPs()){
		imgOut[p] = node->getLevel();
	}
	for(NodeCT* child: node->getChildren()){
		reconstruction(child, imgOut);
	}
}


 ComponentTree::~ComponentTree(){
	for (NodeCT *node: this->listNodes){
		delete node;
		node = nullptr;
	}
	delete[] nodes;
	root = nullptr;
	nodes = nullptr;
 }

ComponentTree::ComponentTree(int* img, int numRows, int numCols, bool isMaxtree) 
	: ComponentTree(img, numRows, numCols, isMaxtree, 1.5){ }


ComponentTree::ComponentTree(int* img, int numRows, int numCols){
	this->numRows = numRows;
	this->numCols = numCols;
	this->treeType = TREE_OF_SHAPES;
	this->nodes = new NodeCT*[this->numRows * this->numCols];

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
    NodeCT** nodes = new NodeCT*[size];
    
	
    this->numNodes = 0;
    for (int i = 0; i < size; i++) {
		int p = imgR[i];
        auto [px, py] = ImageUtils::to2D(p, builder->getInterpNumCols());
		int pixelUnterpolate = (px/2) + (py/2) * numCols;
            
		if (p == parent[p]) { //representante do node raiz
            this->root = nodes[p] = new NodeCT(this->numNodes, pixelUnterpolate, nullptr, imgU[p]);
		}
		else if (imgU[p] != imgU[parent[p]]) { //representante de um node
			nodes[p] = new NodeCT(this->numNodes, pixelUnterpolate, nodes[parent[p]], imgU[p]);
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
		[this](NodeCT* _node) -> void { //pre-processing
			_node->setAreaCC( _node->getCNPs().size() );
			_node->setNumDescendants( _node->getChildren().size() );
			_node->setIndex(this->numNodes++);
			this->listNodes.push_back(_node);
		},
		[](NodeCT* _root, NodeCT* _child) -> void { //merge-processing
			_root->setAreaCC( _root->getAreaCC() + _child->getAreaCC() );
			_root->setNumDescendants( _root->getNumDescendants() + _child->getNumDescendants() );
		},
		[](NodeCT* node) -> void { //post-processing
									
		}
	);
	
	delete builder;
	imgR = nullptr;
	imgU = nullptr;
	parent = nullptr;
	interpolationMin = nullptr;
	interpolationMax = nullptr;
} 

 
ComponentTree::ComponentTree(int* img, int numRows, int numCols, bool isMaxtree, double radiusOfAdjacencyRelation){
	this->numRows = numRows;
	this->numCols = numCols;
	this->treeType = isMaxtree? MAX_TREE : MIN_TREE;

	AdjacencyRelation* adj = new AdjacencyRelation(numRows, numCols, radiusOfAdjacencyRelation);	
	BuilderComponentTreeByUnionFind* builder = new BuilderComponentTreeByUnionFind(img, numRows, numCols, isMaxtree, adj);
	
	int n = this->numRows * this->numCols;
	int* orderedPixels = builder->getOrderedPixels();
	int* parent = builder->getParent();
		
	this->nodes = new NodeCT*[n];
	this->numNodes = 0;
	for (int i = 0; i < n; i++) {
		int p = orderedPixels[i];
		if (p == parent[p]) { //representante do node raiz
			this->root = this->nodes[p] = new NodeCT(this->numNodes++, p, nullptr, img[p]);
			this->listNodes.push_back(this->nodes[p]);
			this->nodes[p]->addCNPs(p);
		}
		else if (img[p] != img[parent[p]]) { //representante de um node
			this->nodes[p] = new NodeCT(this->numNodes++, p, this->nodes[parent[p]], img[p]);
			this->listNodes.push_back(this->nodes[p]);
			this->nodes[p]->addCNPs(p);
			this->nodes[parent[p]]->addChild(this->nodes[p]);
		}
		else if (img[p] == img[parent[p]]) {
			this->nodes[parent[p]]->addCNPs(p);
			this->nodes[p] = this->nodes[parent[p]];
		}
	}
	
	AttributeComputedIncrementally::computerAttribute(this->root,
		[](NodeCT* _node) -> void { //pre-processing
			_node->setAreaCC( _node->getCNPs().size() );
			_node->setNumDescendants( _node->getChildren().size() );
		},
		[](NodeCT* _root, NodeCT* _child) -> void { //merge-processing
			_root->setAreaCC( _root->getAreaCC() + _child->getAreaCC() );
			_root->setNumDescendants( _root->getNumDescendants() + _child->getNumDescendants() );
		},
		[](NodeCT* node) -> void { //post-processing
									
		}
	);
	delete builder;
	delete adj;
	adj = nullptr;
	builder = nullptr;
	orderedPixels = nullptr;
	parent = nullptr;
}

NodeCT* ComponentTree::getSC(int pixel){
	return this->nodes[pixel];
}
	
NodeCT* ComponentTree::getRoot() {
	return this->root;
}

bool ComponentTree::isMaxtree(){
	return this->treeType == MAX_TREE;
}

int ComponentTree::getTreeType(){
	return this->treeType;
}

std::list<NodeCT*> ComponentTree::getListNodes(){
	return this->listNodes;
}

int ComponentTree::getNumNodes(){
	return this->numNodes;
}

int ComponentTree::getNumRowsOfImage(){
	return this->numRows;
}

int ComponentTree::getNumColsOfImage(){
	return this->numCols;
}

int* ComponentTree::getImageAferPruning(NodeCT* nodePruning){
	int n = this->numRows * this->numCols;
	int* imgOut = new int[n];
	std::stack<NodeCT*> s;
	s.push(this->root);
	while(!s.empty()){
		NodeCT* node = s.top();s.pop();
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
			for(NodeCT* child: node->getChildren()){
				s.push(child);
			}
		}
	}
	return imgOut;
}

void ComponentTree::pruning(NodeCT* nodePruning){
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
		free(nodePruning);

	}
}

int* ComponentTree::reconstructionImage(){
	int n = this->numRows * this->numCols;
	int *imgOut = new int[n];
	this->reconstruction(this->root, imgOut);
	return imgOut;
}

int* ComponentTree::getInputImage(){
	int n = this->numRows * this->numCols;
	int* img = new int[n];
	this->reconstruction(this->root, img);
	return img;
}
	