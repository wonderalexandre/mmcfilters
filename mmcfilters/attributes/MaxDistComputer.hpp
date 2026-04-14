#pragma once 

#include "AttributeComputer.hpp"
#include "maxdist/MaxDistAlg.hpp"
#include "../trees/MorphologicalTree.hpp"

namespace mmcfilters
{
  /**
   * @brief Compute maximum squared euclidean distance transform for each n
   *        morphological tree node.
   */
  class MaxDistComputer : public AttributeComputer 
  {
  public:
    std::vector<Attribute> attributes() const override 
    {
      return {MAX_DIST};
    }  

    void compute(MorphologicalTree *tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames,
      const std::vector<Attribute> &requestedAttributes, 
      const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>&) const override
    {
      if (PRINT_LOG)
        std::cout << "\n===== AttributeComputer: Computing MaxDist" << std::endl;

      maxdist::MaxDistAlg maxDistAlg(tree);
      std::vector<float> maxDist = maxDistAlg.getAttributes();

      auto indexOf = [&](NodeId idx) { return attrNames->linearIndex(idx, MAX_DIST); };
      for (NodeId nodeId : tree->getNodeIds()) {
        buffer[indexOf(nodeId)] = maxDist[nodeId];
      }
    }
  };
}