
#include <list>
#include <vector>
#include <stack>
#include <limits.h>

#include "../include/NodeMT.hpp"
#include "../include/MorphologicalTree.hpp"
#include "../include/ComputerMSER.hpp"
#include "../include/AttributeComputedIncrementally.hpp"



    NodeMTPtr ComputerMSER::getNodeAscendant(NodeMTPtr node, int h){
		NodeMTPtr n = node;
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

	void ComputerMSER::maxAreaDescendants(NodeMTPtr nodeAsc, NodeMTPtr nodeDes){
		if(this->descendants[nodeAsc->getIndex()] == nullptr)
			this->descendants[nodeAsc->getIndex()] = nodeDes;
		
		if( this->descendants[nodeAsc->getIndex()]->getAreaCC() < nodeDes->getAreaCC() )
			this->descendants[nodeAsc->getIndex()] = nodeDes;
		
	}
	
	double ComputerMSER::getStability(NodeMTPtr node){
		return (this->getAttrMSER(this->getAscendant(node)) - this->getAttrMSER(this->getDescendant(node))) / this->getAttrMSER(node)  ;
	}

	ComputerMSER::~ComputerMSER(){ }

	ComputerMSER::ComputerMSER(MorphologicalTreePtr tree, std::shared_ptr<float[]> attr_increasing) : ComputerMSER(tree) {
		this->attr_mser = attr_increasing;
	}
	
	ComputerMSER::ComputerMSER(MorphologicalTreePtr tree): tree(tree), attr_mser(nullptr), maxVariation(10.0), minArea(0), maxArea(tree->getNumColsOfImage() * tree->getNumRowsOfImage()) { }

	float ComputerMSER::getAttrMSER(NodeMTPtr node){
		if(node == nullptr)
			return node->getAreaCC();
		else
			return this->attr_mser[node->getIndex()];
	}

	std::vector<bool> ComputerMSER::computerMSER(int delta){

		std::vector<NodeMTPtr> tmp_asc (this->tree->getNumNodes(), nullptr);
		this->ascendants = tmp_asc;

		std::vector<NodeMTPtr> tmp_des (this->tree->getNumNodes(), nullptr);
		this->descendants = tmp_des;

		std::vector<double> tmp_stab (this->tree->getNumNodes(), UNDEF);
		this->stability = tmp_stab;
		
		for(NodeMTPtr node: tree->getIndexNode()){
			NodeMTPtr nodeAsc = this->getNodeAscendant(node, delta);
			this->maxAreaDescendants(nodeAsc, node);
			this->ascendants[node->getIndex()] = nodeAsc;
		}
		
		for(NodeMTPtr node: tree->getIndexNode()){
			if(this->ascendants[node->getIndex()] != nullptr && this->descendants[node->getIndex()] != nullptr){
				this->stability[node->getIndex()] = this->getStability(node);
			}
		}
		
		this->num = 0;
		double maxStabilityDesc, maxStabilityAsc;
		std::vector<bool> mser(this->tree->getNumNodes(), false);
		for(NodeMTPtr node: tree->getIndexNode()){
			if(this->stability[node->getIndex()] != UNDEF && this->stability[this->getAscendant(node)->getIndex()] != UNDEF && this->stability[this->getDescendant(node)->getIndex()] != UNDEF){
				maxStabilityDesc = this->stability[this->getDescendant(node)->getIndex()];
				maxStabilityAsc = this->stability[this->getAscendant(node)->getIndex()];
				if(this->stability[node->getIndex()] < maxStabilityDesc && this->stability[node->getIndex()] < maxStabilityAsc){
					if(stability[node->getIndex()] < this->maxVariation && this->getAttrMSER(node) >= this->minArea && this->getAttrMSER(node) <= this->maxArea){
						mser[node->getIndex()] = true;
						this->num++;
					}
				}
			}
		}
		return mser;
	}

	NodeMTPtr ComputerMSER::getNodeInPathWithMaxStability(NodeMTPtr node, std::vector<bool> isMSER){
		NodeMTPtr nodeAsc = this->ascendantWithMaxStability(node);
		NodeMTPtr nodeDes = this->descendantWithMaxStability(node);
		NodeMTPtr nodeMax = node;


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


	NodeMTPtr ComputerMSER::descendantWithMaxStability(NodeMTPtr node) {
		return this->descendants[node->getIndex()];
	}
	
	std::vector<double> ComputerMSER::getStabilities(){
		return this->stability;
	}

    NodeMTPtr ComputerMSER::ascendantWithMaxStability(NodeMTPtr node) {
		return this->ascendants[node->getIndex()];
	}

	int ComputerMSER::getNumNodes() {
		return  num;
	}

	std::vector<NodeMTPtr> ComputerMSER::getAscendants(){
		return this->ascendants;
	}

	NodeMTPtr ComputerMSER::getAscendant(NodeMTPtr node){
		return this->ascendants[node->getIndex()];
	}
	
	NodeMTPtr ComputerMSER::getDescendant(NodeMTPtr node){
		return this->descendants[node->getIndex()];
	}

	std::vector<NodeMTPtr> ComputerMSER::getDescendants(){
		return this->descendants;
	}

	void ComputerMSER::setMaxVariation(double maxVariation) { this->maxVariation = maxVariation; }
	void ComputerMSER::setMinArea(int minArea) { this->minArea = minArea; }
	void ComputerMSER::setMaxArea(int maxArea) { this->maxArea = maxArea; }