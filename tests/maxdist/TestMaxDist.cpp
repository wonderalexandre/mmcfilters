#include "../../mmcfilters/attributes/maxdist/MaxDistComputer.hpp"
#include "../../mmcfilters/utils/Common.hpp"
#include "../../mmcfilters/utils/Image.hpp"
#include "../../mmcfilters/dataStructure/FastStack.hpp"

#include "Tests.hpp"

using namespace mmcfilters;
using namespace mmcfilters::maxdist;


ImageUInt8Ptr recContour(const std::vector<int> &contour, int nrows, int ncols)
{
  ImageUInt8Ptr img = ImageUInt8::create(nrows, ncols, 0);
  for (int pidx : contour) {
    (*img)[pidx] = 1;
  }

  return img;
}



int main()
{
  uint8_t *data = new u_int8_t[] {
    0, 0, 0, 0, 0, 0, 0,
    0, 4, 4, 4, 7, 7, 7,
    0, 4, 4 ,4, 7, 4, 7,
    0, 4, 4 ,4, 7, 4, 7,
    0, 4, 4 ,4, 7, 4, 7,
    0, 4, 4 ,4, 7, 7, 7,
    0, 0, 0, 0, 0, 0, 0
  };

  ImageUInt8Ptr f = ImageUInt8::fromRaw(data, 7, 7);
  MorphologicalTreePtr tree = std::make_shared<MorphologicalTree>(f, true, 1.0);      // 4-connected max-tree

  printImage(f);
  printTree(tree->getRoot());

  MaxDistComputer computer(tree.get());
  std::vector<float> maxDist = computer.getAttributes();

  FastStack<NodeId> stack(tree->getNumNodes());
  stack.push(tree->getRootById());

  while (!stack.empty()) {
    NodeId nid = stack.pop();
    NodeMT node = tree->proxy(nid);

    std::cout << "maxDist[" << nid << "] = " << maxDist[nid] << "\n";

    for (NodeId cid : node.getChildren()) {
      stack.push(cid);
    }
  }


  return 0;
}