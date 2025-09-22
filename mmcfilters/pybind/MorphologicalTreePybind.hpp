#ifndef COMPONENT_TREE_PYBIND_H
#define COMPONENT_TREE_PYBIND_H


#include "../include/MorphologicalTree.hpp"
#include "../include/NodeMT.hpp"
#include "../include/Common.hpp"

#include "../pybind/PybindUtils.hpp"

#include <stdexcept>
#include <pybind11/numpy.h>



namespace py = pybind11;

class MorphologicalTreePybind;
using MorphologicalTreePybindPtr = std::shared_ptr<MorphologicalTreePybind>;

class MorphologicalTreePybind : public MorphologicalTree {


 public:
    using MorphologicalTree::MorphologicalTree;

    MorphologicalTreePybind(py::array_t<uint8_t> input, int numRows, int numCols, std::string ToSInperpolation="self-dual")
        : MorphologicalTree(ImageUInt8::fromExternal(static_cast<uint8_t*>(input.request().ptr), numRows, numCols), ToSInperpolation) { }

    MorphologicalTreePybind(py::array_t<uint8_t> input, int numRows, int numCols, bool isMaxtree, double radiusOfAdjacencyRelation=1.5)
        : MorphologicalTree(ImageUInt8::fromExternal(static_cast<uint8_t*>(input.request().ptr), numRows, numCols), isMaxtree, radiusOfAdjacencyRelation) { }

    MorphologicalTreePybind() = delete;

    py::array_t<uint8_t> reconstructionImage(){
        ImageUInt8Ptr imgOut = ImageUInt8::create(this->numRows, this->numCols);
        MorphologicalTree::reconstruction(this->root, imgOut->rawData());
        return PybindUtils::toNumpy(imgOut);
    }


    static MorphologicalTreePybindPtr createTreeFromAttributeMapping(py::array_t<float> attrMapping, py::array_t<uint8_t> input, int numRows, int numCols, bool isMaxtree, double radius=1.5) {
        auto buf_attr = attrMapping.request();
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



#endif