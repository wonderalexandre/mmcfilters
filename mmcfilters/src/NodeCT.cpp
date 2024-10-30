#include "../include/NodeCT.hpp"
#include "../include/AdjacencyRelation.hpp"

#include <list>
#include <stdlib.h>

NodeCT::NodeCT(){}

NodeCT::NodeCT(int index, int rep, NodeCT* parent, int level) {
		this->index = index;
        this->rep = rep;
        this->parent = parent;
        this->level = level;
        if(parent == nullptr)
            this->residue = this->level;
        else{
            this->isMaxtree = level > parent->level;
            this->residue = abs(this->level - parent->level);
        }

}

void NodeCT::addCNPs(int p) {
    this->cnps.push_back(p);
}

void NodeCT::addChild(NodeCT* child) {
	this->children.push_back(child);
}

int NodeCT::getRep(){ return this->rep; }

int NodeCT::getIndex(){ return this->index; }

void NodeCT::setIndex(int index) {this->index = index;}

bool NodeCT::isMaxtreeNode(){ return this->isMaxtree; }

int NodeCT::getResidue(){ return this->residue; }

void NodeCT::setResidue(int residue){ this->residue = residue; }

int NodeCT::getLevel(){ return this->level; }

int NodeCT::getAreaCC() { return this->areaCC; }

void NodeCT::setAreaCC(int area) { this->areaCC = area; }

int NodeCT::getNumDescendants() { return this->numDescendants; }

void NodeCT::setNumDescendants(int num) { this->numDescendants = num; }

NodeCT* NodeCT::getParent(){  return this->parent; }

void NodeCT::setParent(NodeCT* parent){ this->parent = parent; }

std::list<int> NodeCT::getCNPs(){  return this->cnps; }

std::list<NodeCT*> NodeCT::getChildren(){  return this->children; }


int NodeCT::getNumSiblings() {
    if(this->parent != nullptr)
		return this->parent->getChildren().size();
	else
		return 0;
}
