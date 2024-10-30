
#include <list>
#include <vector>
#include <stack>
#include <limits.h>

#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"
#include "../include/ComputerMSER.hpp"
#include "../include/AttributeComputedIncrementally.hpp"



    NodeCT* ComputerMSER::getNodeAscendant(NodeCT* node, int h){
		NodeCT* n = node;
		for(int i=0; i <= h; i++){
			if(this->tree->isMaxtree()){
				if(node->getLevel() >= n->getLevel() + h)
					return n;
			}else{
				if(node->getLevel() <= n->getLevel() - h)
					return n;
			}
			if(n->getParent() != nullptr)
				n = n->getParent();
			else 
				return n;
		}
		return n;
	}

	void ComputerMSER::maxAreaDescendants(NodeCT* nodeAsc, NodeCT* nodeDes){
		if(this->descendants[nodeAsc->getIndex()] == nullptr)
			this->descendants[nodeAsc->getIndex()] = nodeDes;
		
		if( this->attr_area[ this->descendants[nodeAsc->getIndex()]->getIndex() ] < this->attr_area[ nodeDes->getIndex() ])
			this->descendants[nodeAsc->getIndex()] = nodeDes;
		
	}
	
	double ComputerMSER::getStability(NodeCT* node){
		return (this->attr_mser[this->getAscendant(node)->getIndex()] - this->attr_mser[this->getDescendant(node)->getIndex()]) / this->attr_mser[node->getIndex()]  ;
	}

	ComputerMSER::~ComputerMSER(){
		delete [] this->attr_area;
	}

	ComputerMSER::ComputerMSER(ComponentTree* tree, double* attr_increasing) : ComputerMSER(tree) {
		this->attr_mser = attr_increasing;
	}
	
	ComputerMSER::ComputerMSER(ComponentTree* tree){
		this->tree = tree;
		this->maxVariation = 10.0;
		this->minArea = 0;
		this->maxArea = tree->getNumColsOfImage() * tree->getNumRowsOfImage();

		double *_attribute = new double[this->tree->getNumNodes()]; 
		AttributeComputedIncrementally::computerAttribute(this->tree->getRoot(),
						[&_attribute](NodeCT* node) -> void {
							_attribute[node->getIndex()] = node->getCNPs().size();
						},
						[&_attribute](NodeCT* root, NodeCT* child) -> void {
							_attribute[root->getIndex()] += _attribute[child->getIndex()];
						},
						[](NodeCT* node) -> void { 							
						});

		this->attr_area = _attribute;
		this->attr_mser = _attribute;
		
	}

	std::vector<bool> ComputerMSER::computerMSER(int delta){

		std::vector<NodeCT*> tmp_asc (this->tree->getNumNodes(), nullptr);
		this->ascendants = tmp_asc;

		std::vector<NodeCT*> tmp_des (this->tree->getNumNodes(), nullptr);
		this->descendants = tmp_des;

		std::vector<double> tmp_stab (this->tree->getNumNodes(), UNDEF);
		this->stability = tmp_stab;
		
		for(NodeCT *node: tree->getListNodes()){
			NodeCT *nodeAsc = this->getNodeAscendant(node, delta);
			this->maxAreaDescendants(nodeAsc, node);
			this->ascendants[node->getIndex()] = nodeAsc;
		}
		
		for(NodeCT *node: tree->getListNodes()){
			if(this->ascendants[node->getIndex()] != nullptr && this->descendants[node->getIndex()] != nullptr){
				this->stability[node->getIndex()] = this->getStability(node);
			}
		}
		
		this->num = 0;
		double maxStabilityDesc, maxStabilityAsc;
		std::vector<bool> mser(this->tree->getNumNodes(), false);
		for(NodeCT *node: tree->getListNodes()){
			if(this->stability[node->getIndex()] != UNDEF && this->stability[this->getAscendant(node)->getIndex()] != UNDEF && this->stability[this->getDescendant(node)->getIndex()] != UNDEF){
				maxStabilityDesc = this->stability[this->getDescendant(node)->getIndex()];
				maxStabilityAsc = this->stability[this->getAscendant(node)->getIndex()];
				if(this->stability[node->getIndex()] < maxStabilityDesc && this->stability[node->getIndex()] < maxStabilityAsc){
					if(stability[node->getIndex()] < this->maxVariation && this->attr_mser[node->getIndex()] >= this->minArea && this->attr_mser[node->getIndex()] <= this->maxArea){
						mser[node->getIndex()] = true;
						this->num++;
					}
				}
			}
		}
		return mser;
	}

	NodeCT* ComputerMSER::getNodeInPathWithMaxStability(NodeCT* node, std::vector<bool> isMSER){
		NodeCT* nodeAsc = this->ascendantWithMaxStability(node);
		NodeCT* nodeDes = this->descendantWithMaxStability(node);
		NodeCT* nodeMax = node;


		double max = stability[node->getIndex()];
        double maxDesc = stability[nodeDes->getIndex()];
        double maxAnc = stability[nodeAsc->getIndex()];
                    
                    if(max <= maxDesc && max <= maxAnc) {
                        return node;
                    }else if (maxDesc <= maxAnc) {
                        return nodeDes;
                    }else {
                       return nodeAsc;
                    }
		
	}


	NodeCT* ComputerMSER::descendantWithMaxStability(NodeCT* node) {
		return this->descendants[node->getIndex()];
	}
	
	std::vector<double> ComputerMSER::getStabilities(){
		return this->stability;
	}

    NodeCT* ComputerMSER::ascendantWithMaxStability(NodeCT* node) {
		return this->ascendants[node->getIndex()];
	}

	int ComputerMSER::getNumNodes() {
		return  num;
	}

	std::vector<NodeCT*> ComputerMSER::getAscendants(){
		return this->ascendants;
	}

	NodeCT* ComputerMSER::getAscendant(NodeCT* node){
		return this->ascendants[node->getIndex()];
	}
	
	NodeCT* ComputerMSER::getDescendant(NodeCT* node){
		return this->descendants[node->getIndex()];
	}

	std::vector<NodeCT*> ComputerMSER::getDescendants(){
		return this->descendants;
	}

	void ComputerMSER::setMaxVariation(double maxVariation) { this->maxVariation = maxVariation; }
	void ComputerMSER::setMinArea(int minArea) { this->minArea = minArea; }
	void ComputerMSER::setMaxArea(int maxArea) { this->maxArea = maxArea; }