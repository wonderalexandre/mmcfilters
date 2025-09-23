#pragma once

#include "../mmcfilters/trees/MorphologicalTree.hpp"
#include "../mmcfilters/trees/NodeMT.hpp"
#include "../mmcfilters/utils/Common.hpp"

#include "PybindUtils.hpp"

#include <stdexcept>
#include <pybind11/numpy.h>

namespace py = pybind11;

class MorphologicalTreePybind;
using MorphologicalTreePybindPtr = std::shared_ptr<MorphologicalTreePybind>;

/**
 * @brief Interface Pybind da árvore morfológica com utilidades de reconstrução.
 */
class MorphologicalTreePybind : public MorphologicalTree {


 public:
    using MorphologicalTree::MorphologicalTree;

    MorphologicalTreePybind(py::array_t<uint8_t, py::array::c_style | py::array::forcecast> input, std::string ToSInperpolation = "self-dual")
        : MorphologicalTree(
              [&]() {
                  auto buf = input.request();
                  if (buf.ndim != 2) {
                      throw std::invalid_argument("input must be a 2D uint8 array");
                  }
                  int rows = static_cast<int>(buf.shape[0]);
                  int cols = static_cast<int>(buf.shape[1]);
                  return ImageUInt8::fromExternal(static_cast<uint8_t*>(buf.ptr), rows, cols);
              }(),
              ToSInperpolation) { }

    MorphologicalTreePybind(py::array_t<uint8_t, py::array::c_style | py::array::forcecast> input, bool isMaxtree, double radiusOfAdjacencyRelation = 1.5)
        : MorphologicalTree(
              [&]() {
                  auto buf = input.request();
                  if (buf.ndim != 2) {
                      throw std::invalid_argument("input must be a 2D uint8 array");
                  }
                  int rows = static_cast<int>(buf.shape[0]);
                  int cols = static_cast<int>(buf.shape[1]);
                  return ImageUInt8::fromExternal(static_cast<uint8_t*>(buf.ptr), rows, cols);
              }(),
              isMaxtree,
              radiusOfAdjacencyRelation) { }

              
    MorphologicalTreePybind() = delete;

    py::array_t<uint8_t> reconstructionImage(){
        ImageUInt8Ptr imgOut = ImageUInt8::create(this->numRows, this->numCols);
        MorphologicalTree::reconstruction(this->root, imgOut->rawData());
        return PybindUtils::toNumpy(imgOut);
    }


    static MorphologicalTreePybindPtr createTreeFromAttributeMapping(py::array_t<float> attrMapping, py::array_t<uint8_t> input, bool isMaxtree, double radius=1.5) {
        auto buf_attr = attrMapping.request();
        if (buf_attr.ndim != 2) {
            throw std::invalid_argument("input must be a 2D float array");
        }
        int numRows = static_cast<int>(buf_attr.shape[0]);
        int numCols = static_cast<int>(buf_attr.shape[1]);

        ImageFloatPtr attributeMapping = ImageFloat::fromExternal(static_cast<float*>(buf_attr.ptr), numRows, numCols);

        auto buf_input = input.request();
        ImageUInt8Ptr img = ImageUInt8::fromExternal(static_cast<uint8_t*>(buf_input.ptr), numRows, numCols);

        MorphologicalTreePtr tree = MorphologicalTree::createFromAttributeMapping(attributeMapping, img, isMaxtree, radius);

        return std::static_pointer_cast<MorphologicalTreePybind>(tree);

    }

    static py::array_t<uint8_t> recNode(NodeMT node) {
        if (!node) {
            throw std::invalid_argument("NodeMT inválido para reconstrução");
        }

        int totalPixels = node.getArea();
        NodeMT parent = node.getParent();
        while (parent) {
            totalPixels = parent.getArea();
            parent = parent.getParent();
        }

        ImageUInt8Ptr imgOut = ImageUInt8::create(totalPixels, 1);
        imgOut->fill(0);
        for (int p : node.getPixelsOfCC()) {
            (*imgOut)[p] = 255;
        }
        return PybindUtils::toNumpy(imgOut);
    }

};

