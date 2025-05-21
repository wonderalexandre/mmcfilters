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
   

    /*
    py::array_t<int> getOrderedPixels(){
        int n = this->numRows * this->numCols;
        return PybindUtils::toNumpy(this->orderedPixels, n);
    }

    py::array_t<int> getParent(){
        int n = this->numRows * this->numCols;
        return PybindUtils::toNumpy(this->parent, n);
    }*/

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

    /*static py::array_t<int> computerParent(py::array_t<int> input, int numRows, int numCols, bool isMaxtree){
		auto buf_input = input.request();
		int* img = (int *) buf_input.ptr;
		ComponentTree tree(img, numRows, numCols, isMaxtree);
		return PybindUtils::toNumpy(tree.getParent(), numRows * numCols);;
	}*/


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