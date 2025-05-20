#include "../include/AttributeOpeningPrimitivesFamily.hpp"
#include "../include/AttributeFilters.hpp"
#include "../include/NodeMT.hpp"
#include "../include/ComputerMSER.hpp"

#include <vector>
#include <stack>
#include <algorithm>
#include <iostream>

AttributeOpeningPrimitivesFamily::~AttributeOpeningPrimitivesFamily(){
    //delete[] this->restOfImage;
}

AttributeOpeningPrimitivesFamily::AttributeOpeningPrimitivesFamily(MorphologicalTreePtr tree, std::shared_ptr<float[]> attrs_increasing, float maxCriterion, int deltaMSER){
  this->tree = tree;
  this->attrs_increasing = attrs_increasing;
  this->maxCriterion = maxCriterion;

  if(deltaMSER > 0){
    ComputerMSER mser(this->tree);
    this->selectedForFiltering = mser.computerMSER(deltaMSER);
  }
  else{
    this->selectedForFiltering.assign(this->tree->getNumNodes(), true);
  }
  
  this->numPrimitives = 0;
  float maxThreshold = 0;
  for(NodeMTPtr node: this->tree->getIndexNode()){
    if(this->attrs_increasing[node->getIndex()] <= this->maxCriterion && this->isSelectedForPruning(node)){
      this->numPrimitives++;
      if(this->attrs_increasing[node->getIndex()] > maxThreshold)
        maxThreshold = this->attrs_increasing[node->getIndex()];
    }
  }
  this->initializeRestOfImage(maxThreshold);
  this->initializeNodesWithMaximumCriterium();
}

AttributeOpeningPrimitivesFamily::AttributeOpeningPrimitivesFamily(MorphologicalTreePtr tree, std::shared_ptr<float[]> attrs_increasing, float maxCriterion): AttributeOpeningPrimitivesFamily(tree, attrs_increasing, maxCriterion, 0){ }

int AttributeOpeningPrimitivesFamily::getNumPrimitives(){
  return this->numPrimitives;
}

std::list<float> AttributeOpeningPrimitivesFamily::getThresholdsPrimitive(){
  if(this->thresholds.size() == 0){
    for(NodeMTPtr node: this->tree->getIndexNode()){
      if(this->attrs_increasing[node->getIndex()] <= this->maxCriterion && this->isSelectedForPruning(node)){
        this->thresholds.push_back(this->attrs_increasing[node->getIndex()]);
      }
    }
    this->thresholds.sort();
    this->thresholds.unique();
  }
  return thresholds;
}

bool AttributeOpeningPrimitivesFamily::hasNodeSelectedInPrimitive(NodeMTPtr currentNode){
  if(!this->selectedForFiltering[currentNode->getIndex()]){
    std::stack<NodeMTPtr> s;
    s.push(currentNode);
    while (!s.empty()){
      NodeMTPtr node = s.top();
      s.pop();
      if (selectedForFiltering[node->getIndex()]){
        return true;
      }

      for (NodeMTPtr son : node->getChildren()){
        if (this->attrs_increasing[son->getIndex()] == this->attrs_increasing[son->getParent()->getIndex()]){ //same primitive?
          s.push(son);
        }
      }
    }
    return false;
  }
  return true;
}

bool AttributeOpeningPrimitivesFamily::isSelectedForPruning(NodeMTPtr node){
  return node->getParent() != nullptr && this->attrs_increasing[node->getIndex()] != this->attrs_increasing[node->getParent()->getIndex()];
}

ImageUInt8Ptr AttributeOpeningPrimitivesFamily::getRestOfImage(){
  return this->restOfImage;
}



void AttributeOpeningPrimitivesFamily::initializeRestOfImage(float thrRestImage){
  this->restOfImage = ImageUInt8::create(this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage());
  AttributeFilters::filteringByPruningMin(this->tree, this->attrs_increasing, thrRestImage, restOfImage);
}

void AttributeOpeningPrimitivesFamily::initializeNodesWithMaximumCriterium(){
  std::stack<NodeMTPtr> s;
  for(NodeMTPtr child: this->tree->getRoot()->getChildren()){
    s.push(child);
  }

  while(!s.empty()){
    NodeMTPtr node = s.top();s.pop();
    if(this->attrs_increasing[this->tree->getRoot()->getIndex()] != this->attrs_increasing[node->getIndex()] && this->attrs_increasing[node->getIndex()] <= this->maxCriterion){
      this->nodesWithMaximumCriterium.push_back(node);
    }
    else{
      for(NodeMTPtr child: node->getChildren()){
        s.push(child);
      }
    }
  }

}

std::list<NodeMTPtr> AttributeOpeningPrimitivesFamily::getNodesWithMaximumCriterium(){
  return this->nodesWithMaximumCriterium;
}



MorphologicalTreePtr AttributeOpeningPrimitivesFamily::getTree(){
  return this->tree;
}