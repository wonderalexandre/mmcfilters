#include "../filters/AttributeOpeningPrimitivesFamily.hpp"
#include "../filters/AttributeFilters.hpp"

namespace mmcfilters {

AttributeOpeningPrimitivesFamily::~AttributeOpeningPrimitivesFamily(){ }

AttributeOpeningPrimitivesFamily::AttributeOpeningPrimitivesFamily(MorphologicalTree* tree, std::shared_ptr<float[]> attrs_increasing, float maxCriterion, int deltaMSER){
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

AttributeOpeningPrimitivesFamily::AttributeOpeningPrimitivesFamily(MorphologicalTree* tree, std::shared_ptr<float[]> attrs_increasing, float maxCriterion): AttributeOpeningPrimitivesFamily(tree, attrs_increasing, maxCriterion, 0){ }

int AttributeOpeningPrimitivesFamily::getNumPrimitives(){
  return this->numPrimitives;
}

std::vector<float> AttributeOpeningPrimitivesFamily::getThresholdsPrimitive(){
  if(this->thresholds.size() == 0){
    for(NodeId node: this->tree->getNodeIds()){
      if(this->attrs_increasing[node] <= this->maxCriterion && this->isSelectedForPruning(node)){
        this->thresholds.push_back(this->attrs_increasing[node]);
      }
    }
    this->make_unique_vector(this->thresholds);
  }
  return thresholds;
}

bool AttributeOpeningPrimitivesFamily::hasNodeSelectedInPrimitive(NodeId currentNode){
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

bool AttributeOpeningPrimitivesFamily::isSelectedForPruning(NodeId node){
  return tree->getParentById(node) != InvalidNode && this->attrs_increasing[node] != this->attrs_increasing[tree->getParentById(node)];
}

ImageUInt8Ptr AttributeOpeningPrimitivesFamily::getRestOfImage(){
  return this->restOfImage;
}



void AttributeOpeningPrimitivesFamily::initializeRestOfImage(float thrRestImage){
  this->restOfImage = ImageUInt8::create(this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage());
  AttributeFilters::filteringByPruningMin(this->tree, this->attrs_increasing, thrRestImage, restOfImage);
}

void AttributeOpeningPrimitivesFamily::initializeNodesWithMaximumCriterium(){
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

std::vector<NodeId> AttributeOpeningPrimitivesFamily::getNodesWithMaximumCriterium(){
  return this->nodesWithMaximumCriterium;
}



MorphologicalTree* AttributeOpeningPrimitivesFamily::getTree(){
  return this->tree;
}

} // namespace mmcfilters
