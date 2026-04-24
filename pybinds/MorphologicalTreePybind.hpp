#pragma once

#include "../mmcfilters/trees/MorphologicalTree.hpp"
#include "../mmcfilters/utils/Common.hpp"

#include "PybindUtils.hpp"

#include <stdexcept>
#include <pybind11/numpy.h>
namespace mmcfilters {

namespace py = pybind11;

class MorphologicalTreePybind;
using MorphologicalTreePybindPtr = std::shared_ptr<MorphologicalTreePybind>;

/**
 * @brief Pybind-facing wrapper for `MorphologicalTree`.
 *
 * This subclass exposes Python-friendly constructors and a few helper methods
 * that return NumPy arrays or Python lists while preserving the ownership
 * model of the underlying C++ tree.
 */
class MorphologicalTreePybind : public MorphologicalTree {
    static std::vector<int> collectPixelsOfConnectedComponent(const MorphologicalTree& tree, NodeId nodeId) {
        std::vector<int> pixels;
        for (NodeId subtreeNodeId : tree.getNodeSubtree(nodeId)) {
            for (int properPart : tree.getProperParts(subtreeNodeId)) {
                pixels.push_back(properPart);
            }
        }
        return pixels;
    }

 public:
    using MorphologicalTree::MorphologicalTree;

    MorphologicalTreePybind(int r, int c, bool m, std::optional<AdjacencyRelation> a = std::nullopt) {
        MorphologicalTree::treeType_ = m ? MAX_TREE : MIN_TREE;
        MorphologicalTree::adj_ = std::move(a);
        MorphologicalTree::numRows_ = r;
        MorphologicalTree::numCols_ = c;
        MorphologicalTree::resetEmptyStorage(static_cast<size_t>(r * c));
    }

    MorphologicalTreePybind(
        py::array_t<uint8_t, py::array::c_style | py::array::forcecast> input,
        ToSInterpolation interpolation = ToSInterpolation::SelfDual)
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
              interpolation) { }

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

    static py::array_t<uint8_t> reconstructNode(const MorphologicalTree& tree, NodeId nodeId) {
        if (!tree.isNode(nodeId) || !tree.isAlive(nodeId)) {
            throw std::invalid_argument("invalid NodeId for reconstruction");
        }

        ImageUInt8Ptr imgOut = ImageUInt8::create(tree.getNumRowsOfImage(), tree.getNumColsOfImage());
        imgOut->fill(0);
        for (int p : collectPixelsOfConnectedComponent(tree, nodeId)) {
            (*imgOut)[p] = 255;
        }
        return PybindUtils::toNumpy(imgOut);
    }


    static py::list representativeProperPartsByFlood(MorphologicalTree& tree, const AltitudeBuffer& altitude, NodeId nodeId) {
        if (!tree.isNode(nodeId) || !tree.isAlive(nodeId)) {
            throw std::invalid_argument("invalid NodeId");
        }

        AdjacencyRelation* adjacency = tree.getAdjacencyRelation();
        if (adjacency == nullptr) {
            throw std::invalid_argument("adjacency relation is unavailable for this tree type");
        }
        const int numPixels = tree.getNumRowsOfImage() * tree.getNumColsOfImage();
        const AltitudeType targetAltitude = tree_altitude_ops::getAltitude(altitude, nodeId);

        auto levelOf = [&](int p) -> int {
            const NodeId smallestComponent = tree.getSmallestComponent(p);
            return tree_altitude_ops::getAltitude(altitude, smallestComponent);
        };
        auto inNode = [&](int p) -> bool {
            return tree.getSmallestComponent(p) == nodeId;
        };

        py::list reps;
        FastQueue<int> queue(1024);
        std::vector<uint8_t> visited(numPixels, 0);

        for (int p : collectPixelsOfConnectedComponent(tree, nodeId)) {
            if (!inNode(p) || visited[p]) {
                continue;
            }

            reps.append(p);
            visited[p] = true;
            queue.push(p);

            while (!queue.empty()) {
                const int properPart = queue.pop();
                for (int q : adjacency->getAdjPixels(properPart)) {
                    if (!inNode(q) || visited[q]) {
                        continue;
                    }
                    if (levelOf(q) != targetAltitude) {
                        continue;
                    }
                    visited[q] = true;
                    queue.push(q);
                }
            }
        }

        return reps;
    }

};

} // namespace mmcfilters
