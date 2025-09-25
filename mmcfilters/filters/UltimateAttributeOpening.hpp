#pragma once

#include "../utils/Common.hpp"
#include "../trees/MorphologicalTree.hpp"
#include "../attributes/ComputerMSER.hpp"
#include "../attributes/AttributeComputedIncrementally.hpp"


namespace mmcfilters {


class UltimateAttributeOpening;
using UltimateAttributeOpeningPtr = std::shared_ptr<UltimateAttributeOpening>;

/**
 * @brief Realiza Ultimate Attribute Opening acumulando contrastes máximos.
 */
class UltimateAttributeOpening{

  protected:
    int maxCriterion;
    std::shared_ptr<float[]> attrs_increasing;
    MorphologicalTree* tree;
    std::shared_ptr<uint8_t[]> maxContrastLUT;
    std::shared_ptr<int[]> associatedIndexLUT;
    std::vector<uint8_t> selectedForFiltering; //mappping between nodes and selected nodes

    void computeUAO(NodeId currentNode, int levelNodeNotInNR, bool qPropag, bool isCalculateResidue){
      NodeId parentNode = tree->getParentById(currentNode);
      int levelNodeInNR = tree->getLevelById(currentNode);
      bool flagPropag = false;
      int contrast = 0;
      if (this->isSelectedForPruning(currentNode)){ // new primitive?
        levelNodeNotInNR = tree->getLevelById(parentNode);

        if (this->attrs_increasing[currentNode] <= this->maxCriterion){ // node selected for pruning = first node in Nr
          isCalculateResidue = hasNodeSelectedInPrimitive(currentNode);
        }
      }

      if (this->attrs_increasing[currentNode] <= this->maxCriterion){

        if (isCalculateResidue) // non Filter?
          contrast = (int)std::abs(levelNodeInNR - levelNodeNotInNR);

        if (this->maxContrastLUT[parentNode] >= contrast){
          this->maxContrastLUT[currentNode] = this->maxContrastLUT[parentNode];
          this->associatedIndexLUT[currentNode] = this->associatedIndexLUT[parentNode];
        }
        else{
          this->maxContrastLUT[currentNode] = contrast;
          if (!qPropag){                                                                                                      // new primitive with max contrast?
            this->associatedIndexLUT[currentNode] = this->attrs_increasing[currentNode] + 1;
          }
          else{
            this->associatedIndexLUT[currentNode] = this->associatedIndexLUT[parentNode];
          }
          flagPropag = true;
        }
      }

      for (NodeId son : tree->getChildrenById(currentNode)){
        this->computeUAO(son, levelNodeNotInNR, flagPropag, isCalculateResidue);
      }
    }

    void execute(int maxCriterion, const std::vector<uint8_t>& selectedForFiltering){
      this->maxCriterion = maxCriterion;
      
      this->selectedForFiltering = selectedForFiltering;

      for (NodeId id: tree->getNodeIds()){
        maxContrastLUT[id] = 0;
        associatedIndexLUT[id] = 0;
      }

      NodeId rootId = tree->getRootById();
      int level = tree->getLevelById(rootId);
      for (NodeId son : tree->getChildrenById(rootId)){
        computeUAO(son, level, false, false);
      }
    }
    
    //first Node in Nr(i)
    bool isSelectedForPruning(NodeId currentNode){
      // primitiva: attribute opening
      return this->attrs_increasing[currentNode] != this->attrs_increasing[ tree->getParentById(currentNode)];
    }

    //has node selected inside Nr(i)
    bool hasNodeSelectedInPrimitive(NodeId currentNode){
      std::stack<NodeId> s;
      s.push(currentNode);
      while (!s.empty()){
        NodeId node = s.top();
        s.pop();
        if (selectedForFiltering[node]){
          return true;
        }

        for (NodeId n : tree->getChildrenById(node)){
          if (this->attrs_increasing[n] == this->attrs_increasing[tree->getParentById(n)]){ // if n in Nr?
            s.push(n);
          }
        }
      }
      return false;
    }

  public:

    UltimateAttributeOpening(MorphologicalTreePtr tree,  std::shared_ptr<float[]> attrs_increasing): UltimateAttributeOpening(tree.get(), std::move(attrs_increasing)) {}
    UltimateAttributeOpening(MorphologicalTree* tree, std::shared_ptr<float[]> attrs_increasing){
      this->tree = tree;
      this->maxContrastLUT = std::shared_ptr<uint8_t[]>(new uint8_t[this->tree->getNumNodes()]);
      this->associatedIndexLUT = std::shared_ptr<int[]>(new int[this->tree->getNumNodes()]);
      this->selectedForFiltering.assign(this->tree->getNumNodes(), true);
      this->attrs_increasing = attrs_increasing;
    }

    ~UltimateAttributeOpening(){ }


    void execute(int maxCriterion){
      std::vector<uint8_t> tmp(this->tree->getNumNodes(), true);
      execute(maxCriterion, tmp);
    }

    void executeWithMSER(int maxCriterion, int deltaMSER){
      ComputerMSER mser(this->tree);
      execute(maxCriterion, mser.computerMSER(deltaMSER));
    }


    ImageUInt8Ptr getMaxConstrastImage(){
      int size = this->tree->getNumColsOfImage() * this->tree->getNumRowsOfImage();
      ImageUInt8Ptr imgOut = ImageUInt8::create(this->tree->getNumColsOfImage(), this->tree->getNumRowsOfImage());
      auto out = imgOut->rawData();

      for (int pidx = 0; pidx < size; pidx++){
        out[pidx] = this->maxContrastLUT[tree->getSCById(pidx)];
      }
      return imgOut;
    }

    ImageInt32Ptr getAssociatedImage(){
      int size = this->tree->getNumColsOfImage() * this->tree->getNumRowsOfImage();
      ImageInt32Ptr imgOut = ImageInt32::create(this->tree->getNumColsOfImage(), this->tree->getNumRowsOfImage());
      auto out = imgOut->rawData();


      for (int pidx = 0; pidx < size; pidx++){
        out[pidx] = this->associatedIndexLUT[tree->getSCById(pidx)];
      }
      return imgOut;
    }
    
    ImageUInt8Ptr getAssociatedColorImage(){
      return ImageUtils::createRandomColor(this->getAssociatedImage()->rawData(), this->tree->getNumRowsOfImage(), this->tree->getNumColsOfImage());
    }
};







	
} // namespace mmcfilters

