#pragma once

#include "../include/ComponentTree.hpp"
#include "../include/Common.hpp"


#define UNDEF -999999999999


class ComputerMSERCT {
private:
	
	ComponentTree* tree;
	float maxVariation;
	float minAttr;
	float maxAttr;
	int num;
	std::vector<float> stability;
	std::vector<NodeId> ascendants;
	std::vector<NodeId> descendants;

	std::shared_ptr<float[]> attr_mser;
	
	
public:
	ComputerMSERCT(ComponentTreePtr tree): ComputerMSERCT(tree.get()) {}
	ComputerMSERCT(ComponentTreePtr tree, std::shared_ptr<float[]> attr_increasing){ ComputerMSERCT(tree.get(), attr_increasing); }
	
	ComputerMSERCT(ComponentTree* tree, std::shared_ptr<float[]> attr_increasing) : ComputerMSERCT(tree) {
		this->attr_mser = attr_increasing;
	}
	
	ComputerMSERCT(ComponentTree* tree): tree(tree), attr_mser(nullptr), maxVariation(10.0), minAttr(0), maxAttr(tree->getNumColsOfImage() * tree->getNumRowsOfImage()) { }

	~ComputerMSERCT(){}

	std::vector<uint8_t> computerMSER(int delta){
		std::pair<std::vector<NodeId>, std::vector<NodeId>> ascDesc = tree->computerAscendantsAndDescendants(delta);
		this->ascendants = std::move(ascDesc.first);
		this->descendants = std::move(ascDesc.second);
		this->stability.assign(tree->getNumNodes(), UNDEF);
		
		for(NodeId node: tree->getNodeIds()){
			if(this->ascendants[node] != InvalidNode && this->descendants[node] != InvalidNode){
				this->stability[node] = this->getStability(node);
			}
		}
		
		this->num = 0;
		double maxStabilityDesc, maxStabilityAsc;
		std::vector<uint8_t> mser(this->tree->getNumNodes(), false);
		for(NodeId node: tree->getNodeIds()){
			if(this->stability[node] != UNDEF && this->stability[this->ascendants[node]] != UNDEF && this->stability[this->descendants[node]] != UNDEF){
				maxStabilityDesc = this->stability[this->descendants[node]];
				maxStabilityAsc = this->stability[this->ascendants[node]];
				if(this->stability[node] < maxStabilityDesc && this->stability[node] < maxStabilityAsc){
					if(stability[node] < this->maxVariation && this->getAttrMSER(node) >= this->minAttr && this->getAttrMSER(node) <= this->maxAttr){
						mser[node] = true;
						this->num++;
					}
				}
			}
		}
		return mser;
	}

	
	double getStability(NodeId node){
		return (this->getAttrMSER(this->ascendants[node]) - this->getAttrMSER(this->descendants[node])) / this->getAttrMSER(node)  ;
	}

	float getAttrMSER(NodeId node){
		if(attr_mser == nullptr)
			return tree->getAreaById(node); 
		else
			return this->attr_mser[node];
	}

	NodeId getNodeInPathWithMaxStability(NodeId node, std::vector<uint8_t> isMSER){
		NodeId nodeAsc = this->ascendants[node];
		NodeId nodeDes = this->descendants[node];
		NodeId nodeMax = node;
                
        if(stability[node] <= stability[nodeDes] && stability[node] <= stability[nodeAsc]) {
            return node;
        }else if (stability[nodeDes] <= stability[nodeAsc]) {
            return nodeDes;
        }else {
            return nodeAsc;
        }
		
	}


	NodeId ascendantWithMaxStability(NodeId node){return this->ascendants[node];}
	NodeId descendantWithMaxStability(NodeId node) {return descendants[node];}
	std::vector<float> getStabilities(){return stability;}
	int getNumNodes() {return  num;}
	void setMaxVariation(float maxVariation) { this->maxVariation = maxVariation; }
	void setMinAttribute(int minAttr) { this->minAttr = minAttr; }
	void setMaxAttribute(int maxAttr) { this->maxAttr = maxAttr; }
};



#include <list>
#include <vector>
#include <stack>
#include <limits.h>

#include "../include/NodeMT.hpp"
#include "../include/MorphologicalTree.hpp"
#include "../include/ComputerMSER.hpp"
#include "../include/AttributeComputedIncrementally.hpp"



