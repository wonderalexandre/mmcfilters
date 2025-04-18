#include "../include/NodeMT.hpp"
#include "../include/AdjacencyRelation.hpp"

#include <list>
#include <stdlib.h>

NodeMT::NodeMT(){}

NodeMT::NodeMT(int index, NodeMTPtr parent, int level) {
		this->index = index;
        this->parent = parent;
        this->level = level;
}

void NodeMT::addCNPs(int p) {
    this->cnps.push_back(p);
}

void NodeMT::addChild(NodeMTPtr child) {
	this->children.push_back(child);
}

void NodeMT::setLevel(int level){
    this->level = level;
}

int NodeMT::getIndex(){ return this->index; }

void NodeMT::setIndex(int index) {this->index = index;}

bool NodeMT::isMaxtreeNode(){ 
    return parent != nullptr && level > parent->level;
}

int NodeMT::getResidue(){ 
    if(parent == nullptr)
        return this->level;
    else    
        return abs(this->level - parent->level);
 }

int NodeMT::getLevel(){ return this->level; }

int NodeMT::getAreaCC() { return this->areaCC; }

void NodeMT::setAreaCC(int area) { this->areaCC = area; }

int NodeMT::getNumDescendants() { 
    return (this->getTimePostOrder() - this->getTimePreOrder() - 1) / 2;
 }

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