#pragma once

#include "../trees/NodeRes.hpp"
#include "../trees/MorphologicalTree.hpp"
#include "../filters/AttributeOpeningPrimitivesFamily.hpp"
#include "../utils/Common.hpp"

namespace mmcfilters {

class ResidualTree;
using ResidualTreePtr = std::shared_ptr<ResidualTree>;
/**
 * @brief Estrutura derivada que organiza resíduos a partir de uma família de primitivas.
 */
class ResidualTree{

    protected:
      NodeResPtr root;
      std::shared_ptr<AttributeOpeningPrimitivesFamily> primitivesFamily;
      MorphologicalTree* tree;
      ImageUInt8Ptr maxContrastLUT;
      std::shared_ptr<int[]> associatedIndexesLUT;
      int numNodes;
      ImageUInt8Ptr restOfImage;
      std::vector<NodeResPtr> nodes;

    public:


    ResidualTree(std::shared_ptr<AttributeOpeningPrimitivesFamily> primitivesFamily) {
      this->primitivesFamily = primitivesFamily;
      this->tree = primitivesFamily->getTree();
      this->nodes.resize(this->tree->getNumNodes());
      this->maxContrastLUT = ImageUInt8::create(this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage()); 
      this->associatedIndexesLUT = std::shared_ptr<int[]>(new int[this->tree->getNumNodes()]);
      this->createTree();
    }

    void createTree(){
      this->numNodes = 0;
      
      this->restOfImage = this->primitivesFamily->getRestOfImage();
      std::vector<NodeId > nodesWithMaximumCriterium = this->primitivesFamily->getNodesWithMaximumCriterium();
      bool isDesirableResidue = false;
      this->root = std::make_shared<NodeRes>(tree, InvalidNode, this->numNodes++, isDesirableResidue);
      for (NodeId nodeMaxCriterion : nodesWithMaximumCriterium){
        this->nodes[tree->getParentById(nodeMaxCriterion)] = this->root; 
        this->nodes[tree->getParentById(nodeMaxCriterion)]->setLevelNodeNotInNR( tree->getLevelById( tree->getParentById(nodeMaxCriterion) ));
        this->nodes[tree->getParentById(nodeMaxCriterion)]->setParent(nullptr);
        
        //computerNodeRes(nodeMaxCriterion);
        
        //connect the nodes in the residual tree
        for (NodeId  currentNode: tree->getNodesOfSubtree( nodeMaxCriterion  )){
          NodeResPtr parent = this->nodes[tree->getParentById(currentNode)];
          if (this->primitivesFamily->isSelectedForPruning(currentNode)){ // first node in Nr(i)
          
            bool isDesirableResidue = this->primitivesFamily->hasNodeSelectedInPrimitive(currentNode);
            this->nodes[currentNode] = std::make_shared<NodeRes>(tree, currentNode, this->numNodes++, isDesirableResidue);
            this->nodes[currentNode]->addNodeInNr(currentNode); 
            
            this->nodes[currentNode]->setLevelNodeNotInNR( tree->getLevelById( tree->getParentById(currentNode) ) );

            this->nodes[currentNode]->setParent(parent);
            parent->addChild(nodes[currentNode]);
          // this->listNodes.push_back(this->nodes[currentNode]);
          }
          else{//Node in Nr(i)
            this->nodes[currentNode] = parent;
            this->nodes[currentNode]->addNodeInNr(currentNode); 
          }
        }
        
      }

      //this->computerMaximumResidues();
      
    }

    NodeResPtr getRoot(){
      return this->root;
    }

    NodeResPtr getNodeRes(NodeId  node){
      return this->nodes[node];
    }

    MorphologicalTree* getCTree(){
      return this->tree;
    }

    ImageUInt8Ptr getRestOfImage(){
      return this->restOfImage;
    }

    void computerMaximumResidues(){
      auto maxContrastLUTRaw = this->maxContrastLUT->rawData();
      for (int id = 0; id < this->tree->getNumNodes(); id++){
        maxContrastLUTRaw[id] = 0;
        this->associatedIndexesLUT[id] = 0;
      }
      std::stack<NodeResPtr> s;
      s.push( this->root );
      while(!s.empty()){
        NodeResPtr nodeRes = s.top(); s.pop();

        for (NodeId nodeCT : nodeRes->getNodeInNr()){
          int levelNodeNotInNR = nodeRes->getLevelNodeNotInNR();
          int levelNodeInNR = tree->getLevelById(nodeCT);
          int contrast = 0;
          NodeId parentNodeCT = tree->getParentById(nodeCT);
          if (nodeRes->isDesirableResidue()) // is desirable residue?
            contrast = (int)std::abs(levelNodeInNR - levelNodeNotInNR);  
          
          if (maxContrastLUTRaw[parentNodeCT] >= contrast){ //propagate max contrast e associeated index
            maxContrastLUTRaw[nodeCT] = maxContrastLUTRaw[parentNodeCT];
            this->associatedIndexesLUT[nodeCT] =  this->associatedIndexesLUT[parentNodeCT];
          }
          else{ //new max contrast

            maxContrastLUTRaw[nodeCT] = contrast;

            bool regionWithMaxContrastIsPropagated = false;
            if(tree->getParentById(nodeCT) != InvalidNode){
              regionWithMaxContrastIsPropagated = maxContrastLUTRaw[tree->getParentById(nodeCT)] < maxContrastLUTRaw[parentNodeCT];
            }

            if (regionWithMaxContrastIsPropagated){   
              this->associatedIndexesLUT[nodeCT] = this->associatedIndexesLUT[parentNodeCT];
            }
            else{
              // new primitive with max contrast?
              this->associatedIndexesLUT[nodeCT] = nodeRes->getAssocieatedIndex();
            }
          }
        }
        
        for (NodeResPtr child : nodeRes->getChildren()){
          s.push( child ); 
        }
      }


    }

    /*
    void computerNodeRes(NodeId currentNode){
      NodeId parentNode = tree->getParentById(currentNode);
      NodeResPtr parent = this->nodes[tree->getParentById(currentNode)];

      if (this->primitivesFamily->isSelectedForPruning(currentNode)){ // first node in Nr(i)
        
        bool isDesirableResidue = this->primitivesFamily->hasNodeSelectedInPrimitive(currentNode);
        this->nodes[currentNode] = new NodeRes(currentNode, this->numNodes++, isDesirableResidue);
        this->nodes[currentNode]->addNodeInNr(currentNode); 
        
        this->nodes[currentNode]->setLevelNodeNotInNR( parentNode->getLevel() );

        this->nodes[currentNode]->setParent(parent);
        parent->addChild(nodes[currentNode]);
        this->listNodes.push_back(this->nodes[currentNode]);
      }
      else{//Node in Nr(i)
        this->nodes[currentNode] = parent;
        this->nodes[currentNode]->addNodeInNr(currentNode); 
      }

      int levelNodeNotInNR = this->nodes[currentNode]->getLevelNodeNotInNR();
      int levelNodeInNR = currentNode->getLevel();
      int contrast = 0; 
      
      if (this->nodes[currentNode]->isDesirableResidue()) // is desirable residue?
        contrast = (int)std::abs(levelNodeInNR - levelNodeNotInNR);  
        
      if (this->maxContrastLUT[parentNode] >= contrast){
        this->maxContrastLUT[currentNode] = this->maxContrastLUT[parentNode];
        this->associatedIndexesLUT[currentNode] =  this->associatedIndexesLUT[parentNode];
      }
      else{
        this->maxContrastLUT[currentNode] = contrast;
        
        bool regionWithMaxContrastIsPropagated = false;
        if(parentNode->getParent() != nullptr){
          regionWithMaxContrastIsPropagated = this->maxContrastLUT[parentNode->getParent()] < this->maxContrastLUT[parentNode];
        }

        if (regionWithMaxContrastIsPropagated){   
          this->associatedIndexesLUT[currentNode] = this->associatedIndexesLUT[parentNode];
        }
        else{
          // new primitive with max contrast?
          this->associatedIndexesLUT[currentNode] = this->nodes[currentNode]->getAssocieatedIndex();
        }
        
      }

      for (NodeId son : currentNode->getChildren()){
        this->computerNodeRes(son);
      }
    }
    */    

    ImageUInt8Ptr filtering(std::vector<uint8_t>& criterion){
      std::stack<NodeResPtr> s;
      for (NodeResPtr node : this->root->getChildren()){
        s.push(node);
      }
      
      std::unique_ptr<int[]> mapLevel(new int[this->tree->getNumNodes()]);
      for(NodeId  nodeCT: this->tree->getNodeIds()){
        mapLevel[nodeCT] = 0;
      } 

      while (!s.empty()){
        NodeResPtr node = s.top(); s.pop();
        for (NodeId nodeCT : node->getNodeInNr()){
          if(tree->getParentById(nodeCT) != InvalidNode){
            if(criterion[node->getRootNr()]){
                mapLevel[nodeCT] =  mapLevel[tree->getParentById(nodeCT)] + tree->getResidueById(nodeCT);
            }else{
                mapLevel[nodeCT] =  mapLevel[tree->getParentById(nodeCT)];
            }
          }
        }
        
        for (NodeResPtr child : node->getChildren()){
          s.push(child);
        }
      }
      ImageUInt8Ptr imgOut = ImageUInt8::create(this->getCTree()->getNumRowsOfImage(), this->getCTree()->getNumColsOfImage());
      for(NodeId  node: tree->getNodeIds()){
        for (int pixel : tree->getCNPsById(node)){
          if(this->tree->isMaxtree())
            (*imgOut)[pixel] = (*this->restOfImage)[pixel] + mapLevel[node];
          else
            (*imgOut)[pixel] = (*this->restOfImage)[pixel] - mapLevel[node];
        }
      }
      
      return imgOut;

    }

    ImageUInt8Ptr getPositiveResidues(){

      std::stack<NodeResPtr> s;
      for (NodeResPtr node : this->root->getChildren()){
        s.push(node);
      }
      
      std::unique_ptr<int[]> mapLevelPos(new int[this->tree->getNumNodes()]);
      for(NodeId  nodeCT: this->tree->getNodeIds()){
        mapLevelPos[nodeCT] = 0;
      } 

      while (!s.empty()){
        NodeResPtr node = s.top(); s.pop();
        for (NodeId nodeCT : node->getNodeInNr()){
          if(tree->getParentById(nodeCT) != InvalidNode){
            if(tree->isMaxtreeNodeById(nodeCT)){
              mapLevelPos[nodeCT] =  mapLevelPos[tree->getParentById(nodeCT)] + tree->getResidueById(nodeCT);
            }else{
              mapLevelPos[nodeCT] =  mapLevelPos[tree->getParentById(nodeCT)];
            }
          }
        }
        for (NodeResPtr child : node->getChildren()){
          s.push(child);
        }
      }

      ImageUInt8Ptr imgOut = ImageUInt8::create(this->getCTree()->getNumRowsOfImage(), this->getCTree()->getNumColsOfImage());
      for(NodeId  node: tree->getNodeIds()){
        for (int pixel : tree->getCNPsById(node)){
          if(this->tree->getTreeType() != MorphologicalTree::MIN_TREE)
            (*imgOut)[pixel] = mapLevelPos[node];
          else
            (*imgOut)[pixel] = 0;
        }
      }
      
      return imgOut;

    }

    ImageUInt8Ptr getNegativeResidues(){

      std::stack<NodeResPtr> s;
      for (NodeResPtr node : this->root->getChildren()){
        s.push(node);
      }
      
      std::unique_ptr<int[]> mapLevelNeg(new int[this->tree->getNumNodes()]);
      for(NodeId  nodeCT: this->tree->getNodeIds()){
        mapLevelNeg[nodeCT] = 0;
      } 

      while (!s.empty()){
        NodeResPtr node = s.top(); s.pop();
        for (NodeId nodeCT : node->getNodeInNr()){
          if(tree->getParentById(nodeCT) != InvalidNode){
            if(!tree->isMaxtreeNodeById(nodeCT)){
              mapLevelNeg[nodeCT] =  mapLevelNeg[tree->getParentById(nodeCT)] + tree->getResidueById(nodeCT);
            }else{
              mapLevelNeg[nodeCT] =  mapLevelNeg[tree->getParentById(nodeCT)];
            }
          }
        }
        for (NodeResPtr child : node->getChildren()){
          s.push(child);
        }
      }

      ImageUInt8Ptr imgOut = ImageUInt8::create(this->getCTree()->getNumRowsOfImage(), this->getCTree()->getNumColsOfImage());
      for(NodeId  node: tree->getNodeIds()){
        for (int pixel : tree->getCNPsById(node)){
          if(this->tree->getTreeType() != MorphologicalTree::MAX_TREE)
            (*imgOut)[pixel] = mapLevelNeg[node];
          else
            (*imgOut)[pixel] = 0;
          
        }
      }

      return imgOut;

    }

    ImageUInt8Ptr reconstruction(){

      std::stack<NodeResPtr> s;
      for (NodeResPtr node : this->root->getChildren()){
        s.push(node);
      }
      
      std::unique_ptr<int[]> mapLevelNeg(new int[this->tree->getNumNodes()]);
      std::unique_ptr<int[]> mapLevelPos(new int[this->tree->getNumNodes()]);
      for(NodeId  nodeCT: this->tree->getNodeIds()){
        mapLevelPos[nodeCT] = 0;
        mapLevelNeg[nodeCT] = 0;
      } 

      while (!s.empty()){
        NodeResPtr node = s.top(); s.pop();
        for (NodeId nodeCT : node->getNodeInNr()){
          if(tree->getParentById(nodeCT) != InvalidNode){
            if(tree->isMaxtreeNodeById(nodeCT)){
              mapLevelPos[nodeCT] =  mapLevelPos[tree->getParentById(nodeCT)] + tree->getResidueById(nodeCT);
              mapLevelNeg[nodeCT] =  mapLevelNeg[tree->getParentById(nodeCT)];
            }else{
              mapLevelNeg[nodeCT] =  mapLevelNeg[tree->getParentById(nodeCT)] + tree->getResidueById(nodeCT);
              mapLevelPos[nodeCT] =  mapLevelPos[tree->getParentById(nodeCT)];
            }
          }
        }

        for (NodeResPtr child : node->getChildren()){
          s.push(child);
        }
      }

      ImageUInt8Ptr imgOut = ImageUInt8::create(this->getCTree()->getNumRowsOfImage(), this->getCTree()->getNumColsOfImage());
      for(NodeId  node: tree->getNodeIds()){
        for (int pixel : tree->getCNPsById(node)){
          (*imgOut)[pixel] = (*this->restOfImage)[pixel] - mapLevelNeg[node] + mapLevelPos[node];
        }
      }
      
      return imgOut;
    }


    ~ResidualTree(){

    }

    /*std::list<NodeResPtr> getListNodes(){
      return this->listNodes;
    }*/

    ImageUInt8Ptr getMaxConstrastImage(){
      return this->maxContrastLUT;
    }

    ImageInt32Ptr getAssociatedImage(){
      int size = this->tree->getNumColsOfImage() * this->tree->getNumRowsOfImage();
      ImageInt32Ptr out = ImageInt32::create(this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage());
      for (int pidx = 0; pidx < size; pidx++){
        (*out)[pidx] = this->associatedIndexesLUT[this->tree->getSCById(pidx)];
      }
      return out;
    }

    ImageUInt8Ptr getAssociatedColorImage(){
      return ImageUtils::createRandomColor(this->getAssociatedImage()->rawData(), this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage());
    }
  };

} // namespace mmcfilters

