#pragma once


#include "../trees/MorphologicalTree.hpp"
#include "../attributes/ComputerMSER.hpp"
#include "../attributes/AttributeComputedIncrementally.hpp"
#include "../utils/Common.hpp"

namespace mmcfilters {

class AttributeOpeningPrimitivesFamily;
using AttributeOpeningPrimitivesFamilyPtr = std::shared_ptr<AttributeOpeningPrimitivesFamily>;

/**
 * @brief Gerencia famílias de primitivas de abertura por atributos.
 *
 * Mantém thresholds, nós selecionados e imagem residual necessários para
 * construir uma árvore residual a partir de uma família de primitivas
 */
class AttributeOpeningPrimitivesFamily{
  
  protected:
    std::shared_ptr<float[]> attrs_increasing;
    float maxCriterion;
    std::vector<float> thresholds;
    std::vector<NodeId> nodesWithMaximumCriterium;

    //PrimitivesFamily
    MorphologicalTree* tree;
    std::vector<uint8_t> selectedForFiltering; //mappping between index nodes and selected nodes
    ImageUInt8Ptr restOfImage;
    int numPrimitives;
    

    void initializeRestOfImage(float threshold);
    void initializeNodesWithMaximumCriterium();
    void make_unique_vector(std::vector<float>& v) {
          std::sort(v.begin(), v.end());
          auto new_end = std::unique(v.begin(), v.end());
          v.erase(new_end, v.end());
    }

  public:
    AttributeOpeningPrimitivesFamily(MorphologicalTreePtr tree,   std::shared_ptr<float[]> attr, float maxCriterion): AttributeOpeningPrimitivesFamily(tree.get(), attr, maxCriterion) {}
    AttributeOpeningPrimitivesFamily(MorphologicalTree* tree,   std::shared_ptr<float[]> attr, float maxCriterion);

    AttributeOpeningPrimitivesFamily(MorphologicalTreePtr tree,   std::shared_ptr<float[]> attrs_increasing, float maxCriterion, int deltaMSER): AttributeOpeningPrimitivesFamily(tree.get(), attrs_increasing, maxCriterion, deltaMSER) {}
    AttributeOpeningPrimitivesFamily(MorphologicalTree* tree,   std::shared_ptr<float[]> attrs_increasing, float maxCriterion, int deltaMSER);

    ~AttributeOpeningPrimitivesFamily();

    std::vector<float> getThresholdsPrimitive();

    //PrimitivesFamily
    bool isSelectedForPruning(NodeId node) ; //first Node in Nr(i)

    bool hasNodeSelectedInPrimitive(NodeId node) ; //has node selected inside Nr(i)

    std::vector<NodeId> getNodesWithMaximumCriterium() ; 

    ImageUInt8Ptr getRestOfImage() ;

    int getNumPrimitives() ;
    
    MorphologicalTree* getTree() ;
};




} // namespace mmcfilters
