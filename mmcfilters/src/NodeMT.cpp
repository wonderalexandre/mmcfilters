#include "../include/NodeMT.hpp"
#include "../include/AdjacencyRelation.hpp"

#include <list>
#include <stdlib.h>

NodeMT::NodeMT(){}

NodeMT::NodeMT(int index, int rep, NodeMTPtr parent, int level) {
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

void NodeMT::addCNPs(int p) {
    this->cnps.push_back(p);
}

void NodeMT::addChild(NodeMTPtr child) {
	this->children.push_back(child);
}

int NodeMT::getRep(){ return this->rep; }

int NodeMT::getIndex(){ return this->index; }

void NodeMT::setIndex(int index) {this->index = index;}

bool NodeMT::isMaxtreeNode(){ return this->isMaxtree; }

int NodeMT::getResidue(){ return this->residue; }

void NodeMT::setResidue(int residue){ this->residue = residue; }

int NodeMT::getLevel(){ return this->level; }

int NodeMT::getAreaCC() { return this->areaCC; }

void NodeMT::setAreaCC(int area) { this->areaCC = area; }

int NodeMT::getNumDescendants() { return this->numDescendants; }

void NodeMT::setNumDescendants(int num) { this->numDescendants = num; }

NodeMTPtr NodeMT::getParent(){  return this->parent; }

void NodeMT::setParent(NodeMTPtr parent){ this->parent = parent; }

std::list<int>& NodeMT::getCNPs()  { return this->cnps; }

std::list<NodeMTPtr>& NodeMT::getChildren(){  return this->children; }


int NodeMT::getNumSiblings() {
    if(this->parent != nullptr)
		return this->parent->getChildren().size();
	else
		return 0;
}

int NodeMT::getTimePostOrder() { return this->timePostOrder; }

void NodeMT::setTimePostOrder(int time) { this->timePostOrder = time; } 

int NodeMT::getTimePreOrder() { return this->timePreOrder; }

void NodeMT::setTimePreOrder(int time) { this->timePreOrder = time; }

bool NodeMT::isLeaf() {
    return this->children.empty();
}