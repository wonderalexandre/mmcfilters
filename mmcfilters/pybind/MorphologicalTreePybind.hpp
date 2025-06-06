#ifndef COMPONENT_TREE_PYBIND_H
#define COMPONENT_TREE_PYBIND_H


#include "../include/MorphologicalTree.hpp"
#include "../include/Common.hpp"

#include "../pybind/PybindUtils.hpp"

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

    py::array_t<uint8_t> getImageAferPruning(NodeMTPtr node){
        int n = this->numRows * this->numCols;
        ImageUInt8Ptr imgOut = MorphologicalTree::getImageAferPruning(node); // Chamar método da superclasse
        return PybindUtils::toNumpy(imgOut);
    }

    py::array_t<uint8_t> reconstructionImage(){
        int n = this->numRows * this->numCols;
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

    static py::array_t<uint8_t> recNode(NodeMTPtr _node) {
        int n = _node->getAreaCC();
        NodeMTPtr parent = _node->getParent();
        while (parent != nullptr) {
            n = parent->getAreaCC();
            parent = parent->getParent();
        }

        ImageUInt8Ptr imgOut = ImageUInt8::create(n, 1);
        for (int p = 0; p < n; p++)
            (*imgOut)[p] = 0;
        for(int p: _node->getPixelsOfCC()){
            (*imgOut)[p] = 255;
        }
        return PybindUtils::toNumpy(imgOut);
    }

};



#endif