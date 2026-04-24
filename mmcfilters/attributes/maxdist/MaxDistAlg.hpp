#pragma once 

#include <vector>

#include "../../utils/Common.hpp"
#include "../../trees/MorphologicalTree.hpp"
#include "../../trees/TreeAltitudeOps.hpp"

#include "EdtDIFT.hpp"
#include "Geometry.hpp"

namespace mmcfilters 
{
  namespace maxdist
  {
    class MaxDistAlg
    {
    public:
      MaxDistAlg(const MorphologicalTree& tree, const AltitudeBuffer* altitude)
        :tree_{tree}, altitude_{altitude}
      {}

      std::vector<float> getAttributes() const
      {
        std::vector<float> maxDist(tree_.getNumInternalNodeSlots(), 0.0f);
        std::array<std::vector<NodeId>, 256> levelToNodes = extractLevelMap();

        ImageUInt8Ptr f = tree_altitude_ops::reconstructImage(tree_, altitude_);
        Box2D domain(f->getNumCols(), f->getNumRows());
        EdtDIFT edtDIFT(f->getNumRows(), f->getNumCols());

        std::unordered_map<NodeId, std::vector<int>> contours;
        ImageInt32 feroded = erode(domain, f);

        for (int level = 255; level >= 0; --level) {
          const std::vector<NodeId> &nodes = levelToNodes[level];

          if (nodes.empty())
            continue;

          for (NodeId nodeId : nodes) {
            contours.insert({nodeId, std::vector<int>()});

            std::vector<int> toRemove;
            std::vector<int> &Ncontour = contours[nodeId];

            for (NodeId childNodeId : tree_.getChildren(nodeId)) {
              for (int pidx : contours[childNodeId]) {
                if (feroded[pidx] < static_cast<int>(tree_altitude_ops::getAltitude(altitude_, nodeId)))
                  Ncontour.push_back(pidx);
                else
                  toRemove.push_back(pidx);
              }
              contours.erase(childNodeId);
            }

            if (!toRemove.empty())
              edtDIFT.treeRemoval(toRemove);

            for (int pidx : tree_.getProperParts(nodeId)) {
              edtDIFT.addPixelToBinaryImage(pidx);

              if (feroded[pidx] < static_cast<int>(tree_altitude_ops::getAltitude(altitude_, nodeId))) {
                Ncontour.push_back(pidx);
                edtDIFT.seed(pidx);
              }
              else {
                edtDIFT.open(pidx);
                edtDIFT.insertNeighborsPQueue(pidx);
              }
            }
          }

          edtDIFT.run();

          for (NodeId nodeId : nodes) {
            maxDist[nodeId] = edtDIFT.maxBedt(contours[nodeId]);
          }
        }

        return maxDist;
      }

    private:
      std::array<std::vector<NodeId>, 256> extractLevelMap() const
      {
        std::array<std::vector<NodeId>, 256> levelToNodes;
        FastStack<int> stack{static_cast<size_t>(tree_.getNumInternalNodeSlots())};

        stack.push(tree_.getRoot());

        while (!stack.empty()) {
          const NodeId nodeId = stack.pop();
          levelToNodes[static_cast<size_t>(tree_altitude_ops::getAltitude(altitude_, nodeId))].push_back(nodeId);

          for (NodeId childNodeId : tree_.getChildren(nodeId)) {
            stack.push(childNodeId);
          }
        }

        return levelToNodes;
      }

      ImageInt32 erode(const Box2D& domain, ImageUInt8Ptr f) const     // erosion by a cross in-place
      {
        const Point2D cross[] = { Point2D(-1, 0), Point2D(0, -1), Point2D(1 ,0), Point2D(0, 1) };
        const int OUT_OF_DOMAIN_VAL = -1;

        ImageInt32 feroded(domain.width(), domain.height());

        for (int pidx = 0; pidx < feroded.getSize(); ++pidx) {
          int minVal = static_cast<int>((*f)[pidx]);
          Point2D p = domain.point(pidx);

          for (const Point2D& q : cross) {
            Point2D n = p + q;
            if (domain.contains(n)) {
              int nval = (*f)[domain.index(n)];
              if (nval < minVal)
                minVal = nval;
            }
            else {
              minVal = OUT_OF_DOMAIN_VAL;
              break;
            }
          }
          feroded[pidx] = minVal;
        }

        return feroded;
      }


    private:
      const MorphologicalTree& tree_;
      const AltitudeBuffer* altitude_;
    };
  }
}
