#include "../include/AttributeOpeningPrimitivesFamily.hpp"
#include "../include/AttributeFilters.hpp"
#include "../include/NodeCT.hpp"
#include "../include/ComputerMSER.hpp"

#include <vector>
#include <stack>
#include <algorithm>
#include <iostream>

AttributeOpeningPrimitivesFamily::~AttributeOpeningPrimitivesFamily(){
    delete[] this->restOfImage;
}

AttributeOpeningPrimitivesFamily::AttributeOpeningPrimitivesFamily(ComponentTree* tree,  float* attrs_increasing, float maxCriterion, int deltaMSER){
  this->tree = tree;
  this->attrs_increasing = attrs_increasing;
  this->maxCriterion = maxCriterion;

  if(deltaMSER > 0){
    ComputerMSER *mser = new ComputerMSER(this->tree);
    this->selectedForFiltering = mser->computerMSER(deltaMSER);
    delete mser;
  }
  else{
    std::vector<bool> tmp(this->tree->getNumNodes(), true);
    this->selectedForFiltering = tmp;
  }
  
  this->numPrimitives = 0;
  float maxThreshold = 0;
  for(NodeCT* node: this->tree->getListNodes()){
    if(this->attrs_increasing[node->getIndex()] <= this->maxCriterion && this->isSelectedForPruning(node)){
      this->numPrimitives++;
      if(this->attrs_increasing[node->getIndex()] > maxThreshold)
        maxThreshold = this->attrs_increasing[node->getIndex()];
    }
  }
  this->initializeRestOfImage(maxThreshold);
  this->initializeNodesWithMaximumCriterium();
}

AttributeOpeningPrimitivesFamily::AttributeOpeningPrimitivesFamily(ComponentTree* tree,  float* attrs_increasing, float maxCriterion): AttributeOpeningPrimitivesFamily(tree, attrs_increasing, maxCriterion, 0){ }

int AttributeOpeningPrimitivesFamily::getNumPrimitives(){
  return this->numPrimitives;
}

std::list<float> AttributeOpeningPrimitivesFamily::getThresholdsPrimitive(){
  if(this->thresholds.size() == 0){
    for(NodeCT* node: this->tree->getListNodes()){
      if(this->attrs_increasing[node->getIndex()] <= this->maxCriterion && this->isSelectedForPruning(node)){
        this->thresholds.push_back(this->attrs_increasing[node->getIndex()]);
      }
    }
    this->thresholds.sort();
    this->thresholds.unique();
  }
  return thresholds;
}

bool AttributeOpeningPrimitivesFamily::hasNodeSelectedInPrimitive(NodeCT *currentNode){
  if(!this->selectedForFiltering[currentNode->getIndex()]){
    std::stack<NodeCT *> s;
    s.push(currentNode);
    while (!s.empty()){
      NodeCT *node = s.top();
      s.pop();
      if (selectedForFiltering[node->getIndex()]){
        return true;
      }

      for (NodeCT* son : node->getChildren()){
        if (this->attrs_increasing[son->getIndex()] == this->attrs_increasing[son->getParent()->getIndex()]){ //same primitive?
          s.push(son);
        }
      }
    }
    return false;
  }
  return true;
}

bool AttributeOpeningPrimitivesFamily::isSelectedForPruning(NodeCT* node){
  return node->getParent() != nullptr && this->attrs_increasing[node->getIndex()] != this->attrs_increasing[node->getParent()->getIndex()];
}

int* AttributeOpeningPrimitivesFamily::getRestOfImage(){
  return this->restOfImage;
}



void AttributeOpeningPrimitivesFamily::initializeRestOfImage(float thrRestImage){
  this->restOfImage = new int[this->tree->getNumRowsOfImage() * this->tree->getNumColsOfImage()];
  AttributeFilters::filteringByPruningMin(this->tree, this->attrs_increasing, thrRestImage, restOfImage);
}

void AttributeOpeningPrimitivesFamily::initializeNodesWithMaximumCriterium(){
  std::stack<NodeCT*> s;
  for(NodeCT* child: this->tree->getRoot()->getChildren()){
    s.push(child);
  }

  while(!s.empty()){
    NodeCT* node = s.top();s.pop();
    if(this->attrs_increasing[this->tree->getRoot()->getIndex()] != this->attrs_increasing[node->getIndex()] && this->attrs_increasing[node->getIndex()] <= this->maxCriterion){
      this->nodesWithMaximumCriterium.push_back(node);
    }
    else{
      for(NodeCT* child: node->getChildren()){
        s.push(child);
      }
    }
  }

}

std::list<NodeCT*> AttributeOpeningPrimitivesFamily::getNodesWithMaximumCriterium(){
  return this->nodesWithMaximumCriterium;
}



ComponentTree* AttributeOpeningPrimitivesFamily::getTree(){
  return this->tree;
}