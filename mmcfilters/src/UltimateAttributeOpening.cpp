

#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"
#include "../include/AttributeComputedIncrementally.hpp"
#include "../include/UltimateAttributeOpening.hpp"
#include "../include/ComputerMSER.hpp"

#include "../include/ImageUtils.hpp"
#include <cstdlib>
#include <vector>

void UltimateAttributeOpening::execute(int maxCriterion){
  std::vector<bool> tmp(this->tree->getNumNodes(), true);
  execute(maxCriterion, tmp);
}

void UltimateAttributeOpening::executeWithMSER(int maxCriterion, int deltaMSER){
  ComputerMSER *mser = new ComputerMSER(this->tree);
  execute(maxCriterion, mser->computerMSER(deltaMSER));
}

void UltimateAttributeOpening::execute(int maxCriterion, std::vector<bool> selectedForFiltering){
  this->maxCriterion = maxCriterion;
  
  this->selectedForFiltering = selectedForFiltering;

  for (int id = 0; id < this->tree->getNumNodes(); id++){
    maxContrastLUT[id] = 0;
    associatedIndexLUT[id] = 0;
  }

  for (NodeCT *son : this->tree->getRoot()->getChildren()){
    computeUAO(son, this->tree->getRoot()->getLevel(), false, false);
  }
}

UltimateAttributeOpening::UltimateAttributeOpening(ComponentTree *tree, std::vector<float> attrs_increasing){
  this->tree = tree;
  this->maxContrastLUT = new int[this->tree->getNumNodes()];
  this->associatedIndexLUT = new int[this->tree->getNumNodes()];
  this->attrs_increasing = attrs_increasing;
}

UltimateAttributeOpening::~UltimateAttributeOpening(){
  free(maxContrastLUT);
  free(associatedIndexLUT);
}


int* UltimateAttributeOpening::getMaxConstrastImage(){
  int size = this->tree->getNumColsOfImage() * this->tree->getNumRowsOfImage();
  int *out = new int[size];

  for (int pidx = 0; pidx < size; pidx++){
    out[pidx] = this->maxContrastLUT[this->tree->getSC(pidx)->getIndex()];
  }
  return out;
}

int *UltimateAttributeOpening::getAssociatedImage(){
  int size = this->tree->getNumColsOfImage() * this->tree->getNumRowsOfImage();
  int *out = new int[size];

  for (int pidx = 0; pidx < size; pidx++){
    out[pidx] = this->associatedIndexLUT[this->tree->getSC(pidx)->getIndex()];
  }
  return out;
}
 
int* UltimateAttributeOpening::getAssociatedColorImage(){
  return ImageUtils::createRandomColor(this->getAssociatedImage(), this->tree->getNumColsOfImage(), this->tree->getNumRowsOfImage());
}

void UltimateAttributeOpening::computeUAO(NodeCT *currentNode, int levelNodeNotInNR, bool qPropag, bool isCalculateResidue){
  NodeCT *parentNode = currentNode->getParent();
  int levelNodeInNR = currentNode->getLevel();
  bool flagPropag = false;
  int contrast = 0;
  if (this->isSelectedForPruning(currentNode)){ // new primitive?
    levelNodeNotInNR = parentNode->getLevel();

    if (this->attrs_increasing[currentNode->getIndex()] <= this->maxCriterion){ // node selected for pruning = first node in Nr
      isCalculateResidue = hasNodeSelectedInPrimitive(currentNode);
    }
  }

  if (this->attrs_increasing[currentNode->getIndex()] <= this->maxCriterion){

    if (isCalculateResidue) // non Filter?
      contrast = (int)std::abs(levelNodeInNR - levelNodeNotInNR);

    if (this->maxContrastLUT[parentNode->getIndex()] >= contrast){
      this->maxContrastLUT[currentNode->getIndex()] = this->maxContrastLUT[parentNode->getIndex()];
      this->associatedIndexLUT[currentNode->getIndex()] = this->associatedIndexLUT[parentNode->getIndex()];
    }
    else{
      this->maxContrastLUT[currentNode->getIndex()] = contrast;
      if (!qPropag){                                                                                                      // new primitive with max contrast?
        this->associatedIndexLUT[currentNode->getIndex()] = this->attrs_increasing[currentNode->getIndex()] + 1;
      }
      else{
        this->associatedIndexLUT[currentNode->getIndex()] = this->associatedIndexLUT[parentNode->getIndex()];
      }
      flagPropag = true;
    }
  }

  for (NodeCT *son : currentNode->getChildren()){
    this->computeUAO(son, levelNodeNotInNR, flagPropag, isCalculateResidue);
  }
}

bool UltimateAttributeOpening::isSelectedForPruning(NodeCT *currentNode){
  // primitiva: attribute opening
  return this->attrs_increasing[currentNode->getIndex()] != this->attrs_increasing[currentNode->getParent()->getIndex()];
}

bool UltimateAttributeOpening::hasNodeSelectedInPrimitive(NodeCT *currentNode){
  std::stack<NodeCT *> s;
  s.push(currentNode);
  while (!s.empty()){
    NodeCT *node = s.top();
    s.pop();
    if (selectedForFiltering[node->getIndex()]){
      return true;
    }

    for (NodeCT *n : node->getChildren()){
      if (this->attrs_increasing[n->getIndex()] == this->attrs_increasing[n->getParent()->getIndex()]){ // if n in Nr?
        s.push(n);
      }
    }
  }
  return false;
}
