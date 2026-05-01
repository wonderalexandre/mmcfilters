#pragma once 

#include "AttributeComputer.hpp"
#include "maxdist/MaxDistAlg.hpp"
namespace mmcfilters{

  /**
   * @brief Computes the maximum Euclidean distance transform value per node.
   *
   * @details
   * `MAX_DIST` measures, for each component represented by the tree, the
   * largest Euclidean distance from a point in the component to its boundary.
   * The heavy lifting is delegated to `maxdist::MaxDistAlg`, which operates on
   * the image-domain support induced by the tree.
   *
   * Because the distance transform is defined with respect to image-domain
   * neighbourhoods, the tree must carry a valid adjacency relation. Importing
   * a hierarchy without connectivity metadata is therefore insufficient for
   * this attribute.
   */
  class MaxDistComputer : public AttributeComputer {

    public:

    /**
     * @brief Returns the single descriptor produced by this computer.
     */
    std::vector<Attribute> attributes() const override {
      return {MAX_DIST};
    }  

    /**
     * @brief Computes `MAX_DIST` for each live node of the tree.
     */
    void compute(const MorphologicalTree& tree, const AltitudeBuffer* altitude, std::span<float> buffer, const AttributeNames& attrNames, [[maybe_unused]] std::span<const Attribute> requestedAttributes,std::span<const DependencySource>) const override{
      if (!tree.hasAdjacencyRelation()) {
        throw std::invalid_argument("MAX_DIST requires an adjacency relation.");
      }
      if (altitude == nullptr) {
        throw std::invalid_argument("MAX_DIST requires an explicit altitude buffer. Use WeightedMorphologicalTree or provide an explicit altitude buffer.");
      }

      if (PRINT_LOG)
        std::cout << "\n===== AttributeComputer: Computing MaxDist" << std::endl;

      maxdist::MaxDistAlg maxDistAlg(tree, altitude);
      std::vector<float> maxDist = maxDistAlg.getAttributes();

      auto indexOf = [&](NodeId idx) { return attrNames.linearIndex(idx, MAX_DIST); };
      for (NodeId nodeId : tree.getAliveNodeIds()) {
        buffer[indexOf(nodeId)] = maxDist[nodeId];
      }
    }

    void computeUnitAttributes(
        const MorphologicalTree& tree,
        const AltitudeBuffer*,
        std::span<const NodeId> unitProperParts,
        std::span<float> buffer,
        const AttributeNames& attrNames,
        std::span<const Attribute> requestedAttributes) const override
    {
      requireUnitAttributeBufferShape(tree, unitProperParts, buffer, attrNames);
      if (!requestsAttribute(requestedAttributes, MAX_DIST)) {
        return;
      }
      for (NodeId leafIndex = 0; leafIndex < static_cast<NodeId>(unitProperParts.size()); ++leafIndex) {
        buffer[attrNames.linearIndex(leafIndex, MAX_DIST)] = 0.0f;
      }
    }

  };

}
