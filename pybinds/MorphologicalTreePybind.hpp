#pragma once

#include "../mmcfilters/trees/MorphologicalTree.hpp"
#include "../mmcfilters/trees/NodeMT.hpp"
#include "../mmcfilters/utils/Common.hpp"

#include "PybindUtils.hpp"

#include <stdexcept>
#include <pybind11/numpy.h>
namespace mmcfilters {

namespace py = pybind11;

class MorphologicalTreePybind;
using MorphologicalTreePybindPtr = std::shared_ptr<MorphologicalTreePybind>;

/**
 * @brief Interface Pybind da árvore morfológica com utilidades de reconstrução.
 */
class MorphologicalTreePybind : public MorphologicalTree {


 public:
    using MorphologicalTree::MorphologicalTree;

    MorphologicalTreePybind(int r, int c, bool m, AdjacencyRelationPtr a) {
        MorphologicalTree::numRows = r;
        MorphologicalTree::numCols = c;
        MorphologicalTree::treeType = m ? MAX_TREE : MIN_TREE;
        MorphologicalTree::adj = a;
        MorphologicalTree::pixelToNodeId.resize(r * c, -1);
    }

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


    static MorphologicalTreePybindPtr createTreeFromAttributeMapping(
        py::array_t<float, py::array::c_style | py::array::forcecast> attrMapping,
        py::array_t<uint8_t, py::array::c_style | py::array::forcecast> input,
        bool isMaxtree, double radius = 1.5) {
        auto buf_attr = attrMapping.request();
        if (buf_attr.ndim != 2) {
            throw std::invalid_argument("attrMapping must be a 2D float32 array");
        }
        auto buf_input = input.request();
        if (buf_input.ndim != 2) {
            throw std::invalid_argument("input image must be a 2D uint8 array");
        }
        int numRows = static_cast<int>(buf_input.shape[0]);
        int numCols = static_cast<int>(buf_input.shape[1]);
        if (buf_attr.shape[0] != buf_input.shape[0] || buf_attr.shape[1] != buf_input.shape[1]) {
            throw std::invalid_argument("attrMapping and input must have identical shapes");
        }
        ImageUInt8Ptr img = ImageUInt8::fromExternal(static_cast<uint8_t*>(buf_input.ptr), numRows, numCols);
        ImageFloatPtr attributeMapping = ImageFloat::fromExternal(static_cast<float*>(buf_attr.ptr), numRows, numCols);
        
        auto tree = std::make_shared<MorphologicalTreePybind>(
            numRows, numCols, isMaxtree,
            std::make_shared<AdjacencyRelation>(numRows, numCols, radius));
        // call the overload that RECEIVES an existing instance (by reference) and builds it:
        MorphologicalTree::createFromAttributeMapping(tree, attributeMapping, img, isMaxtree, radius);
        return tree;
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


    static py::list repCNPsByFlood(NodeMT node)  {
        if (!node) {
            throw std::invalid_argument("NodeMT inválido.");
        }
        MorphologicalTree* tree = node.getTree();
        AdjacencyRelation* adj = tree->getAdjacencyRelation(); 
        const int N = tree->getNumRowsOfImage() * tree->getNumColsOfImage();
        const int target = node.getLevel(); // nível da flatzone

        auto levelOf = [&](int p) -> int {
            // Acesso ao nível do pixel p via SC-id (conforme seu item 3)
            return tree->getLevelById(tree->getSCById(p));
        };
        auto inNode = [&](int p) -> int {
            return tree->getSCById(p) == node.getIndex();
        };
        
        
        py::list reps;
        FastQueue<int> Q(1024);
        std::vector<uint8_t> visited(N, 0);
        for (int p : node.getPixelsOfCC()) {
            if (!inNode(p) || visited[p]) continue;

            // Encontramos uma nova flatzone; s será o representante
            reps.append(p);
            visited[p] = true;
            Q.push(p);

            while (!Q.empty()) {
                int cnp = Q.pop();
                
                // Ajuste esta chamada se sua adjacência tiver outra API
                for (int q : adj->getAdjPixels(cnp)) {
                    if (!inNode(q) || visited[q]) continue;
                    if (levelOf(q) != target) continue; // zona plana: mesmo nível
                    visited[q] = true;
                    Q.push(q);
                }
            }
        }

        return reps;
    }

};

} // namespace mmcfilters
