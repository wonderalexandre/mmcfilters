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
    
    MorphologicalTreePybind(py::array_t<PixelValueType> input, int numRows, int numCols, std::string ToSInperpolation="self-dual")
        : MorphologicalTree(Image::fromExternal(static_cast<PixelValueType*>(input.request().ptr), numRows, numCols), ToSInperpolation) { }

	MorphologicalTreePybind(py::array_t<PixelValueType> input, int numRows, int numCols, bool isMaxtree, double radiusOfAdjacencyRelation=1.5)
        : MorphologicalTree(Image::fromExternal(static_cast<PixelValueType*>(input.request().ptr), numRows, numCols), isMaxtree, radiusOfAdjacencyRelation) { }
   

    /*
    py::array_t<int> getOrderedPixels(){
        int n = this->numRows * this->numCols;
        return PybindUtils::toNumpy(this->orderedPixels, n);
    }

    py::array_t<int> getParent(){
        int n = this->numRows * this->numCols;
        return PybindUtils::toNumpy(this->parent, n);
    }*/

    py::array_t<PixelValueType> getImageAferPruning(NodeMTPtr node){
        int n = this->numRows * this->numCols;
        ImagePtr imgOut = MorphologicalTree::getImageAferPruning(node); // Chamar método da superclasse
        return PybindUtils::toNumpy(imgOut->rawData(), n);
    }

    py::array_t<PixelValueType> reconstructionImage(){
        int n = this->numRows * this->numCols;
        PixelValueType* imgOut = new PixelValueType[n];
        MorphologicalTree::reconstruction(this->root, imgOut);
        return PybindUtils::toNumpy(imgOut, n);
    }

    /*static py::array_t<int> computerParent(py::array_t<int> input, int numRows, int numCols, bool isMaxtree){
		auto buf_input = input.request();
		int* img = (int *) buf_input.ptr;
		ComponentTree tree(img, numRows, numCols, isMaxtree);
		return PybindUtils::toNumpy(tree.getParent(), numRows * numCols);;
	}*/


    static py::array_t<PixelValueType> recNode(NodeMTPtr _node) {
        int n = _node->getAreaCC();
        NodeMTPtr parent = _node->getParent();
        while (parent != nullptr) {
            n = parent->getAreaCC();
            parent = parent->getParent();
        }

        auto img_numpy = py::array(py::buffer_info(
            nullptr, sizeof(PixelValueType), py::format_descriptor<PixelValueType>::value,
            1, {n}, {sizeof(PixelValueType)}
        ));
        auto buf_img = img_numpy.request();
        PixelValueType* imgOut = (PixelValueType*) buf_img.ptr;
        for (int p = 0; p < n; p++)
            imgOut[p] = 0;
        for(int p: _node->getPixelsOfCC()){
            imgOut[p] = 255;
        }
        return img_numpy;
    }

};



#endif