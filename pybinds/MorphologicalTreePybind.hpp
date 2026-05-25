#pragma once

#include "../mmcfilters/trees/MorphologicalTree.hpp"
#include "../mmcfilters/utils/Image.hpp"

#include "PybindUtils.hpp"

#include <cstdint>
#include <stdexcept>
#include <vector>

#include <pybind11/numpy.h>

namespace mmcfilters {

namespace py = pybind11;

/**
 * @brief Pybind-facing wrapper for `MorphologicalTree`.
 *
 * This subclass exposes helper methods that return NumPy arrays while
 * preserving the ownership model of the underlying C++ tree. Public Python
 * construction is centralized in `MorphologicalTreeFactory`; this wrapper
 * remains as the Python topology type returned by factory/import paths.
 */
class MorphologicalTreePybind : public MorphologicalTree {
    static std::vector<int> collectPixelsOfConnectedComponent(const MorphologicalTree& tree, NodeId nodeId) {
        std::vector<int> pixels;
        for (int properPart : tree.getConnectedComponent(nodeId)) {
            pixels.push_back(properPart);
        }
        return pixels;
    }

 public:
    explicit MorphologicalTreePybind(MorphologicalTree&& tree)
        : MorphologicalTree(std::move(tree)) {}

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
};

using MorphologicalTreePybindPtr = std::shared_ptr<MorphologicalTreePybind>;

} // namespace mmcfilters
