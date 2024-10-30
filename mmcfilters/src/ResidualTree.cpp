
#include <list>
#include <vector>
#include <stack>
#include <utility> // Para std::pair e std::make_pair

#include "../include/NodeCT.hpp"
#include "../include/NodeRes.hpp"
#include "../include/ImageUtils.hpp"
#include "../include/ResidualTree.hpp"
#include "../include/AttributeOpeningPrimitivesFamily.hpp"


ResidualTree::ResidualTree(AttributeOpeningPrimitivesFamily* primitivesFamily) {
  this->primitivesFamily = primitivesFamily;
  this->tree = primitivesFamily->getTree();
  this->nodes = new NodeRes*[this->tree->getNumNodes()];
  for(int i = 0; i < this->tree->getNumNodes(); i++){
    this->nodes[i] = nullptr;
  }
  this->maxContrastLUT = new int[this->tree->getNumNodes()];
  this->associatedIndexesLUT = new int[this->tree->getNumNodes()];
  this->createTree();
}

void ResidualTree::createTree(){
  this->numNodes = 0;
  
  this->restOfImage = this->primitivesFamily->getRestOfImage();
  std::list<NodeCT*> nodesWithMaximumCriterium = this->primitivesFamily->getNodesWithMaximumCriterium();
  bool isDesirableResidue = false;
  this->root = new NodeRes(nullptr, this->numNodes++, isDesirableResidue);
  //this->listNodes.push_back(this->root);  
  for (NodeCT *nodeMaxCriterion : nodesWithMaximumCriterium){
    this->nodes[nodeMaxCriterion->getParent()->getIndex()] = this->root; 
    this->nodes[nodeMaxCriterion->getParent()->getIndex()]->setLevelNodeNotInNR(nodeMaxCriterion->getParent()->getLevel());
    this->nodes[nodeMaxCriterion->getParent()->getIndex()]->setParent(nullptr);
    
    //computerNodeRes(nodeMaxCriterion);
    
    //connect the nodes in the residual tree
    for (NodeCT* currentNode: nodeMaxCriterion->getNodesDescendants()){
      NodeRes* parent = this->nodes[currentNode->getParent()->getIndex()];
      if (this->primitivesFamily->isSelectedForPruning(currentNode)){ // first node in Nr(i)
      
        bool isDesirableResidue = this->primitivesFamily->hasNodeSelectedInPrimitive(currentNode);
        this->nodes[currentNode->getIndex()] = new NodeRes(currentNode, this->numNodes++, isDesirableResidue);
        this->nodes[currentNode->getIndex()]->addNodeInNr(currentNode); 
        
        this->nodes[currentNode->getIndex()]->setLevelNodeNotInNR( currentNode->getParent()->getLevel() );

        this->nodes[currentNode->getIndex()]->setParent(parent);
        parent->addChild(nodes[currentNode->getIndex()]);
       // this->listNodes.push_back(this->nodes[currentNode->getIndex()]);
      }
      else{//Node in Nr(i)
        this->nodes[currentNode->getIndex()] = parent;
        this->nodes[currentNode->getIndex()]->addNodeInNr(currentNode); 
      }
    }
     
  }

  //this->computerMaximumResidues();
  
}

NodeRes* ResidualTree::getRoot(){
  return this->root;
}

NodeRes* ResidualTree::getNodeRes(NodeCT* node){
  return this->nodes[node->getIndex()];
}

ComponentTree* ResidualTree::getCTree(){
  return this->tree;
}

int* ResidualTree::getRestOfImage(){
  return this->restOfImage;
}

void ResidualTree::computerMaximumResidues(){

  for (int id = 0; id < this->tree->getNumNodes(); id++){
    this->maxContrastLUT[id] = 0;
    this->associatedIndexesLUT[id] = 0;
  }
  std::stack<NodeRes*> s;
  s.push( this->root );
  while(!s.empty()){
    NodeRes* nodeRes = s.top(); s.pop();

    for (NodeCT *nodeCT : nodeRes->getNodeInNr()){
      int levelNodeNotInNR = nodeRes->getLevelNodeNotInNR();
      int levelNodeInNR = nodeCT->getLevel();
      int contrast = 0;
      NodeCT *parentNodeCT = nodeCT->getParent();
      if (nodeRes->isDesirableResidue()) // is desirable residue?
        contrast = (int)std::abs(levelNodeInNR - levelNodeNotInNR);  
      
      if (this->maxContrastLUT[parentNodeCT->getIndex()] >= contrast){ //propagate max contrast e associeated index
        this->maxContrastLUT[nodeCT->getIndex()] = this->maxContrastLUT[parentNodeCT->getIndex()];
        this->associatedIndexesLUT[nodeCT->getIndex()] =  this->associatedIndexesLUT[parentNodeCT->getIndex()];
      }
      else{ //new max contrast

        this->maxContrastLUT[nodeCT->getIndex()] = contrast;

        bool regionWithMaxContrastIsPropagated = false;
        if(parentNodeCT->getParent() != nullptr){
          regionWithMaxContrastIsPropagated = this->maxContrastLUT[parentNodeCT->getParent()->getIndex()] < this->maxContrastLUT[parentNodeCT->getIndex()];
        }

        if (regionWithMaxContrastIsPropagated){   
          this->associatedIndexesLUT[nodeCT->getIndex()] = this->associatedIndexesLUT[parentNodeCT->getIndex()];
        }
        else{
          // new primitive with max contrast?
          this->associatedIndexesLUT[nodeCT->getIndex()] = nodeRes->getAssocieatedIndex();
        }
      }
    }
    
    for (NodeRes *child : nodeRes->getChildren()){
      s.push( child ); 
    }
  }


}

/*
void ResidualTree::computerNodeRes(NodeCT *currentNode){
  NodeCT *parentNode = currentNode->getParent();
  NodeRes* parent = this->nodes[currentNode->getParent()->getIndex()];

  if (this->primitivesFamily->isSelectedForPruning(currentNode)){ // first node in Nr(i)
    
    bool isDesirableResidue = this->primitivesFamily->hasNodeSelectedInPrimitive(currentNode);
    this->nodes[currentNode->getIndex()] = new NodeRes(currentNode, this->numNodes++, isDesirableResidue);
    this->nodes[currentNode->getIndex()]->addNodeInNr(currentNode); 
    
    this->nodes[currentNode->getIndex()]->setLevelNodeNotInNR( parentNode->getLevel() );

    this->nodes[currentNode->getIndex()]->setParent(parent);
    parent->addChild(nodes[currentNode->getIndex()]);
    this->listNodes.push_back(this->nodes[currentNode->getIndex()]);
  }
  else{//Node in Nr(i)
     this->nodes[currentNode->getIndex()] = parent;
     this->nodes[currentNode->getIndex()]->addNodeInNr(currentNode); 
  }

  int levelNodeNotInNR = this->nodes[currentNode->getIndex()]->getLevelNodeNotInNR();
  int levelNodeInNR = currentNode->getLevel();
  int contrast = 0; 
  
  if (this->nodes[currentNode->getIndex()]->isDesirableResidue()) // is desirable residue?
    contrast = (int)std::abs(levelNodeInNR - levelNodeNotInNR);  
    
  if (this->maxContrastLUT[parentNode->getIndex()] >= contrast){
    this->maxContrastLUT[currentNode->getIndex()] = this->maxContrastLUT[parentNode->getIndex()];
    this->associatedIndexesLUT[currentNode->getIndex()] =  this->associatedIndexesLUT[parentNode->getIndex()];
  }
  else{
    this->maxContrastLUT[currentNode->getIndex()] = contrast;
    
    bool regionWithMaxContrastIsPropagated = false;
    if(parentNode->getParent() != nullptr){
      regionWithMaxContrastIsPropagated = this->maxContrastLUT[parentNode->getParent()->getIndex()] < this->maxContrastLUT[parentNode->getIndex()];
    }

    if (regionWithMaxContrastIsPropagated){   
      this->associatedIndexesLUT[currentNode->getIndex()] = this->associatedIndexesLUT[parentNode->getIndex()];
    }
    else{
      // new primitive with max contrast?
      this->associatedIndexesLUT[currentNode->getIndex()] = this->nodes[currentNode->getIndex()]->getAssocieatedIndex();
    }
    
  }

  for (NodeCT *son : currentNode->getChildren()){
    this->computerNodeRes(son);
  }
}
*/    

int* ResidualTree::filtering(std::vector<bool> criterion, int* imgOutput){
  std::stack<NodeRes*> s;
  for (NodeRes *node : this->root->getChildren()){
    s.push(node);
  }
  
  int mapLevel[this->tree->getNumNodes()];
  for(NodeCT* nodeCT: this->tree->getListNodes()){
    mapLevel[nodeCT->getIndex()] = 0;
  } 

  while (!s.empty()){
    NodeRes *node = s.top(); s.pop();
    for (NodeCT *nodeCT : node->getNodeInNr()){
      if(nodeCT->getParent() != nullptr){
        if(criterion[node->getRootNr()->getIndex()]){
            mapLevel[nodeCT->getIndex()] =  mapLevel[nodeCT->getParent()->getIndex()] + nodeCT->getResidue();
        }else{
            mapLevel[nodeCT->getIndex()] =  mapLevel[nodeCT->getParent()->getIndex()];
        }
      }
    }
    
    for (NodeRes *child : node->getChildren()){
      s.push(child);
    }
  }

  for(NodeCT* node: tree->getListNodes()){
    for (int pixel : node->getCNPs()){
      if(this->tree->isMaxtree())
        imgOutput[pixel] = this->restOfImage[pixel] + mapLevel[node->getIndex()];
      else
        imgOutput[pixel] = this->restOfImage[pixel] - mapLevel[node->getIndex()];
    }
  }
  
  return imgOutput;

}

int* ResidualTree::getPositiveResidues(){

  std::stack<NodeRes*> s;
  for (NodeRes *node : this->root->getChildren()){
    s.push(node);
  }
  
  int mapLevelPos[this->tree->getNumNodes()];
  for(NodeCT* nodeCT: this->tree->getListNodes()){
    mapLevelPos[nodeCT->getIndex()] = 0;
  } 

  while (!s.empty()){
    NodeRes *node = s.top(); s.pop();
    for (NodeCT *nodeCT : node->getNodeInNr()){
      if(nodeCT->getParent() != nullptr){
        if(nodeCT->isMaxtreeNode()){
          mapLevelPos[nodeCT->getIndex()] =  mapLevelPos[nodeCT->getParent()->getIndex()] + nodeCT->getResidue();
        }else{
          mapLevelPos[nodeCT->getIndex()] =  mapLevelPos[nodeCT->getParent()->getIndex()];
        }
      }
    }
    for (NodeRes *child : node->getChildren()){
      s.push(child);
    }
  }

  int* imgOutput = new int[this->tree->getNumRowsOfImage() * this->tree->getNumColsOfImage()];
  //for(int p = 0; p < this->tree->getNumRowsOfImage() * this->tree->getNumColsOfImage(); p++)
  //  imgOutput[p] = 0;
  for(NodeCT* node: tree->getListNodes()){
    for (int pixel : node->getCNPs()){
      if(this->tree->getTreeType() != ComponentTree::MIN_TREE)
        imgOutput[pixel] = mapLevelPos[node->getIndex()];
      else
        imgOutput[pixel] = 0;
    }
  }
  
  return imgOutput;

}

int* ResidualTree::getNegativeResidues(){

  std::stack<NodeRes*> s;
  for (NodeRes *node : this->root->getChildren()){
    s.push(node);
  }
  
  int mapLevelNeg[this->tree->getNumNodes()];
  for(NodeCT* nodeCT: this->tree->getListNodes()){
    mapLevelNeg[nodeCT->getIndex()] = 0;
  } 

  while (!s.empty()){
    NodeRes *node = s.top(); s.pop();
    for (NodeCT *nodeCT : node->getNodeInNr()){
      if(nodeCT->getParent() != nullptr){
        if(!nodeCT->isMaxtreeNode()){
          mapLevelNeg[nodeCT->getIndex()] =  mapLevelNeg[nodeCT->getParent()->getIndex()] + nodeCT->getResidue();
        }else{
          mapLevelNeg[nodeCT->getIndex()] =  mapLevelNeg[nodeCT->getParent()->getIndex()];
        }
      }
    }
    for (NodeRes *child : node->getChildren()){
      s.push(child);
    }
  }

  int* imgOutput = new int[this->tree->getNumRowsOfImage() * this->tree->getNumColsOfImage()];
  //for(int p = 0; p < this->tree->getNumRowsOfImage() * this->tree->getNumColsOfImage(); p++)
  //  imgOutput[p] = 0;
  for(NodeCT* node: tree->getListNodes()){
    for (int pixel : node->getCNPs()){
      if(this->tree->getTreeType() != ComponentTree::MAX_TREE)
        imgOutput[pixel] = mapLevelNeg[node->getIndex()];
      else
        imgOutput[pixel] = 0;
      
    }
  }

  return imgOutput;

}

int* ResidualTree::reconstruction(){

  std::stack<NodeRes*> s;
  for (NodeRes *node : this->root->getChildren()){
    s.push(node);
  }
  
  int mapLevelPos[this->tree->getNumNodes()];
  int mapLevelNeg[this->tree->getNumNodes()];
  for(NodeCT* nodeCT: this->tree->getListNodes()){
    mapLevelPos[nodeCT->getIndex()] = 0;
    mapLevelNeg[nodeCT->getIndex()] = 0;
  } 

  while (!s.empty()){
    NodeRes *node = s.top(); s.pop();
    for (NodeCT *nodeCT : node->getNodeInNr()){
      if(nodeCT->getParent() != nullptr){
        if(nodeCT->isMaxtreeNode()){
          mapLevelPos[nodeCT->getIndex()] =  mapLevelPos[nodeCT->getParent()->getIndex()] + nodeCT->getResidue();
          mapLevelNeg[nodeCT->getIndex()] =  mapLevelNeg[nodeCT->getParent()->getIndex()];
        }else{
          mapLevelNeg[nodeCT->getIndex()] =  mapLevelNeg[nodeCT->getParent()->getIndex()] + nodeCT->getResidue();
          mapLevelPos[nodeCT->getIndex()] =  mapLevelPos[nodeCT->getParent()->getIndex()];
        }
      }
    }

    for (NodeRes *child : node->getChildren()){
      s.push(child);
    }
  }

  int* imgOutput = new int[this->tree->getNumRowsOfImage() * this->tree->getNumColsOfImage()];
  for(NodeCT* node: tree->getListNodes()){
    for (int pixel : node->getCNPs()){
       imgOutput[pixel] = this->restOfImage[pixel] - mapLevelNeg[node->getIndex()] + mapLevelPos[node->getIndex()];
    }
  }
  
  return imgOutput;
}


ResidualTree::~ResidualTree(){
  delete[] this->nodes;
  delete[] this->maxContrastLUT;
  delete[] this->associatedIndexesLUT;
  delete this->root;
}

/*std::list<NodeRes*> ResidualTree::getListNodes(){
  return this->listNodes;
}*/

int* ResidualTree::getMaxConstrastImage(){
  int size = this->tree->getNumColsOfImage() * this->tree->getNumRowsOfImage();
  int *out = new int[size];
  for (int pidx = 0; pidx < size; pidx++){
    out[pidx] = this->maxContrastLUT[this->tree->getSC(pidx)->getIndex()];
  }
  return out;
}

int* ResidualTree::getAssociatedImage(){
  int size = this->tree->getNumColsOfImage() * this->tree->getNumRowsOfImage();
  int *out = new int[size];
  for (int pidx = 0; pidx < size; pidx++){
    out[pidx] = this->associatedIndexesLUT[this->tree->getSC(pidx)->getIndex()];
  }
  return out;
}

int* ResidualTree::getAssociatedColorImage(){
  return ImageUtils::createRandomColor(this->getAssociatedImage(), this->tree->getNumColsOfImage(), this->tree->getNumRowsOfImage());
}
