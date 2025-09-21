#include "../include/AttributeOpeningPrimitivesFamilyCT.hpp"
#include "../include/AttributeFiltersCT.hpp"
#include "../include/NodeCT.hpp"
#include "../include/ComputerMSERCT.hpp"

#include <vector>
#include <stack>
#include <algorithm>
#include <iostream>

AttributeOpeningPrimitivesFamilyCT::~AttributeOpeningPrimitivesFamilyCT(){
    //delete[] this->restOfImage;
}

AttributeOpeningPrimitivesFamilyCT::AttributeOpeningPrimitivesFamilyCT(ComponentTree* tree, std::shared_ptr<float[]> attrs_increasing, float maxCriterion, int deltaMSER){
  this->tree = tree;
  this->attrs_increasing = attrs_increasing;
  this->maxCriterion = maxCriterion;

  if(deltaMSER > 0){
    ComputerMSERCT mser(this->tree);
    this->selectedForFiltering = mser.computerMSER(deltaMSER);
  }
  else{
    this->selectedForFiltering.assign(this->tree->getNumNodes(), true);
  }
  
  this->numPrimitives = 0;
  float maxThreshold = 0;
  for(NodeId node: this->tree->getNodeIds()){
    if(this->attrs_increasing[node] <= this->maxCriterion && this->isSelectedForPruning(node)){
      this->numPrimitives++;
      if(this->attrs_increasing[node] > maxThreshold)
        maxThreshold = this->attrs_increasing[node];
    }
  }
  this->initializeRestOfImage(maxThreshold);
  this->initializeNodesWithMaximumCriterium();
}

AttributeOpeningPrimitivesFamilyCT::AttributeOpeningPrimitivesFamilyCT(ComponentTree* tree, std::shared_ptr<float[]> attrs_increasing, float maxCriterion): AttributeOpeningPrimitivesFamilyCT(tree, attrs_increasing, maxCriterion, 0){ }

int AttributeOpeningPrimitivesFamilyCT::getNumPrimitives(){
  return this->numPrimitives;
}

std::list<float> AttributeOpeningPrimitivesFamilyCT::getThresholdsPrimitive(){
  if(this->thresholds.size() == 0){
    for(NodeId node: this->tree->getNodeIds()){
      if(this->attrs_increasing[node] <= this->maxCriterion && this->isSelectedForPruning(node)){
        this->thresholds.push_back(this->attrs_increasing[node]);
      }
    }
    this->thresholds.sort();
    this->thresholds.unique();
  }
  return thresholds;
}

bool AttributeOpeningPrimitivesFamilyCT::hasNodeSelectedInPrimitive(NodeId currentNode){
  if(!this->selectedForFiltering[currentNode]){
    std::stack<NodeId> s;
    s.push(currentNode);
    while (!s.empty()){
      NodeId node = s.top();
      s.pop();
      if (selectedForFiltering[node]){
        return true;
      }

      for (NodeId son : tree->getChildrenById(node)){
        if (this->attrs_increasing[son] == this->attrs_increasing[tree->getParentById(son)]){ //same primitive?
          s.push(son);
        }
      }
    }
    return false;
  }
  return true;
}

bool AttributeOpeningPrimitivesFamilyCT::isSelectedForPruning(NodeId node){
  return tree->getParentById(node) != InvalidNode && this->attrs_increasing[node] != this->attrs_increasing[tree->getParentById(node)];
}

ImageUInt8Ptr AttributeOpeningPrimitivesFamilyCT::getRestOfImage(){
  return this->restOfImage;
}



void AttributeOpeningPrimitivesFamilyCT::initializeRestOfImage(float thrRestImage){
  this->restOfImage = ImageUInt8::create(this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage());
  AttributeFiltersCT::filteringByPruningMin(this->tree, this->attrs_increasing, thrRestImage, restOfImage);
}

void AttributeOpeningPrimitivesFamilyCT::initializeNodesWithMaximumCriterium(){
  std::stack<NodeId> s;
  for(NodeId child: tree->getChildrenById(tree->getRootById())){
    s.push(child);
  }

  while(!s.empty()){
    NodeId node = s.top();s.pop();
    if(this->attrs_increasing[tree->getRootById()] != this->attrs_increasing[node] && this->attrs_increasing[node] <= this->maxCriterion){
      this->nodesWithMaximumCriterium.push_back(node);
    }
    else{
      for(NodeId child: tree->getChildrenById(node)){
        s.push(child);
      }
    }
  }

}

std::list<NodeId> AttributeOpeningPrimitivesFamilyCT::getNodesWithMaximumCriterium(){
  return this->nodesWithMaximumCriterium;
}



ComponentTree* AttributeOpeningPrimitivesFamilyCT::getTree(){
  return this->tree;
}