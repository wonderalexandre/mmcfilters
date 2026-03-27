#pragma once 

#include <vector>

#include "../../utils/Common.hpp"
#include "../../trees/MorphologicalTree.hpp"

#include "Geometry.hpp"

namespace mmcfilters 
{
  namespace maxdist
  {
    class MaxDistComputer
    {
    public:
      MaxDistComputer(MorphologicalTree *tree);

      std::vector<float> getAttributes() const;

    private:
      std::array<std::vector<NodeId>, 256> extractLevelMap() const;
      ImageInt32 erode(const Box2D& domain, ImageUInt8Ptr f) const;     // erosion by a cross in-place


    private:
      MorphologicalTree *tree_;
    };
  }
}