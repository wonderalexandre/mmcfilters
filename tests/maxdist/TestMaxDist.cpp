#include "../../mmcfilters/attributes/maxdist/MaxDistAlg.hpp"
#include "../../mmcfilters/utils/Common.hpp"
#include "../../mmcfilters/utils/Image.hpp"
#include "../../mmcfilters/dataStructure/FastStack.hpp"
#include "../../mmcfilters/attributes/AttributeComputedIncrementally.hpp"
#include "../../mmcfilters/filters/ExtinctionValues.hpp"

#include "Tests.hpp"

#define STB_IMAGE_IMPLEMENTATION
#include "../../external/stb/stb_image.h"


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
  // uint8_t *data = new u_int8_t[] {
  //   0, 0, 0, 0, 0, 0, 0,
  //   0, 4, 4, 4, 7, 7, 7,
  //   0, 4, 4 ,4, 7, 4, 7,
  //   0, 4, 4 ,4, 7, 4, 7,
  //   0, 4, 4 ,4, 7, 4, 7,
  //   0, 4, 4 ,4, 7, 7, 7,
  //   0, 0, 0, 0, 0, 0, 0
  // };

  // uint8_t *data = new u_int8_t[] {
  //   7, 7, 7, 7, 7, 7, 7,
  //   7, 7, 7, 7, 7, 7, 7,
  //   7, 7, 7, 7, 7, 7, 7,
  //   7, 7, 7, 7, 7, 7, 7,
  //   7, 7, 7, 7, 7, 7, 7,
  //   7, 7, 7, 7, 7, 7, 7,
  //   7, 7, 7, 7, 7, 7, 7,
  // };


  // int width, height, nchannels;
  // uint8_t *data = stbi_load("../../dat/maxdist-test.png", &width, &height, &nchannels, 0);

  // ImageUInt8Ptr f = ImageUInt8::fromExternal(data, height, width);


  ImageUInt8Ptr f = getWonderImage();
  MorphologicalTreePtr tree = std::make_shared<MorphologicalTree>(f, true, 1.5); 

  auto [attrName, maxdist] = AttributeComputedIncrementally::computeSingleAttribute(tree, Attribute::MAX_DIST);
  ExtinctionValues extval{tree, maxdist};
  
  std::cout << "number of leaves: " << tree->getLeaves().size() << "\n\n";

  const std::vector<RegionalExtremaNode> leaves = extval.getExtinctionValues();
  for (const RegionalExtremaNode &node : leaves) {
    std::cout << "leaf.id: " << node.leaf << "\n";
    std::cout << "cutoff: " << node.cutoffNode << "\n";
    std::cout << "extinction value: " << node.extinction << "\n";
    std::cout << "\n";
  }

  // printImage(f);
  // printTree(tree->getRoot());

  // MaxDistAlg computer(tree.get());
  // std::vector<float> maxDist = computer.getAttributes();

  // FastStack<NodeId> stack(tree->getNumNodes());
  // stack.push(tree->getRootById());

  // while (!stack.empty()) {
  //   NodeId nid = stack.pop();
  //   NodeMT node = tree->proxy(nid);

  //   std::cout << "maxDist[" << nid << "] = " << maxDist[nid] << "\n";

  //   for (NodeId cid : node.getChildren()) {
  //     stack.push(cid);
  //   }
  // }

  // stbi_image_free(data);
  // data = nullptr;
  return 0;
}