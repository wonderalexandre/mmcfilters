#include "../include/NodeRes.hpp"

NodeRes::NodeRes(NodeMTPtr rootNr, int associeatedIndex, bool desirableResidue){
    this->rootNr = rootNr;
    this->associeatedIndex = associeatedIndex;
    this->desirableResidue = desirableResidue;
}

void NodeRes::addNodeInNr(NodeMTPtr node){
    this->nodes.push_back(node);
}

void NodeRes::addChild(NodeResPtr child){
    this->children.push_back(child);
}

void NodeRes::setParent(NodeResPtr parent){
    this->parent = parent;
}

int NodeRes::getAssocieatedIndex(){
    return this->associeatedIndex;
}

bool NodeRes::isDesirableResidue(){
    return this->desirableResidue;
}

std::list<NodeMTPtr> NodeRes::getNodeInNr(){
    return this->nodes;
}

std::list<NodeResPtr> NodeRes::getChildren(){
    return this->children;
}

NodeMTPtr NodeRes::getRootNr(){
    return this->rootNr;
}

NodeResPtr NodeRes::getParent(){
    return this->parent;
}

int NodeRes::getLevelNodeNotInNR(){
    return this->levelNodeNotInNR;
}

void NodeRes::setLevelNodeNotInNR(int level){
    this->levelNodeNotInNR = level;
}