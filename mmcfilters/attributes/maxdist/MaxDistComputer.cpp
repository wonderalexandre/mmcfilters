#include "MaxDistComputer.hpp"
#include "EdtDIFT.hpp"
#include "../../trees/NodeMT.hpp"
#include "../../../tests/Tests.hpp"

#include <map>

namespace mmcfilters
{
  namespace maxdist
  {
    MaxDistComputer::MaxDistComputer(MorphologicalTree *tree)
      :tree_{tree}
    {}

    std::array<std::vector<NodeId>, 256> MaxDistComputer::extractLevelMap() const
    {
      std::array<std::vector<NodeId>, 256> levelToNodes;
      FastStack<int> stack{static_cast<size_t>(tree_->getNumNodes())};
      
      stack.push(tree_->getRoot().getIndex());

      while (!stack.empty()) {
        NodeId nid = stack.pop();
        NodeMT node = tree_->proxy(nid);

        levelToNodes[node.getLevel()].push_back(node.getIndex());

        for (NodeMT c : node.getChildren()) {
          stack.push(c.getIndex());
        }
      }
      
      return levelToNodes;
    }

    std::vector<float> MaxDistComputer::getAttributes() const 
    {
      std::vector<float> maxDist(tree_->getNumNodes(), 0.0f);
      std::array<std::vector<NodeId>, 256> levelToNodes = extractLevelMap();
      
      ImageUInt8Ptr f = tree_->reconstructionImage();
      Box2D domain(f->getNumCols(), f->getNumRows());
      EdtDIFT edtDIFT(f->getNumRows(), f->getNumCols());

      std::unordered_map<NodeId, std::vector<int>> contours;
      ImageInt32 feroded = erode(domain, f);

      // process the level sets from 255 down to 0
      for (int level = 255; level >= 0; --level) {
        const std::vector<NodeId> &nodes = levelToNodes[level];

        // skip level that does not contain nodes
        if (nodes.empty())
          continue;

        // There exist at least one node in "level", so, we have to process them
        for (NodeId nid : nodes) {
          NodeMT node = tree_->proxy(nid);

          // define the contour for the node
          contours.insert({nid, std::vector<int>()});

          // contour pixels to be removed
          std::vector<int> toRemove;

          // define current node contour
          std::vector<int> &Ncontour = contours[nid];

          // resue children pixels and collect pixels which need to be removed 
          // in the DIFT
          for (const NodeMT &c : node.getChildren()) {
            for (int pidx : contours[c.getIndex()]) {
              if (feroded[pidx] < static_cast<int>(node.getLevel())) 
                Ncontour.push_back(pidx);
              else
                toRemove.push_back(pidx);
            }
            contours.erase(c.getIndex());
          }

         if (!toRemove.empty())
            edtDIFT.treeRemoval(toRemove);

          // compute new contour points
          for (int pidx : node.getCNPs()) {            
            // Incrementally create level-set binary image
            edtDIFT.addPixelToBinaryImage(pidx);

            // Check if CNP is contour pixel 
            if (feroded[pidx] < static_cast<int>(node.getLevel())) {
              // pidx us a contour pixel
              // add to contour               
              Ncontour.push_back(pidx);

              // set pidx as DIFT seed
              edtDIFT.seed(pidx);
            }
            else {
              edtDIFT.open(pidx);
              edtDIFT.insertNeighborsPQueue(pidx);
            }
          } // end of cnps       
        } // end of node in the level loop

        // The roots, cost and borders are set.
        // Compute maximum distance attribute for 
        // the nodes in the level
        edtDIFT.run();

        const ImageInt32 c = edtDIFT.cost();
        ImageUInt8Ptr fc = ImageUInt8::create(c.getNumRows(), c.getNumCols());

        for (int pidx = 0; pidx < c.getSize(); pidx++) {
          (*fc)[pidx] = static_cast<unsigned char>(c[pidx]);
        }

        printImage(fc);

        for (NodeId nid : nodes) 
          maxDist[nid] = edtDIFT.maxBedt(contours[nid]);
      }  // level loop

      return maxDist;
    }

    ImageInt32 MaxDistComputer::erode(const Box2D& domain, ImageUInt8Ptr f) const
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
  }
}
