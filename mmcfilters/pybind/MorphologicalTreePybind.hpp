#ifndef COMPONENT_TREE_PYBIND_H
#define COMPONENT_TREE_PYBIND_H


#include "../include/MorphologicalTree.hpp"
#include "../pybind/PybindUtils.hpp"

#include <pybind11/numpy.h>



namespace py = pybind11;

class MorphologicalTreePybind;
using MorphologicalTreePybindPtr = std::shared_ptr<MorphologicalTreePybind>;

class MorphologicalTreePybind : public MorphologicalTree {


 public:
    using MorphologicalTree::MorphologicalTree;

    MorphologicalTreePybind(py::array_t<int> input, int numRows, int numCols)
        : MorphologicalTree(static_cast<int*>(input.request().ptr), numRows, numCols) { }

    MorphologicalTreePybind(py::array_t<int> input, int numRows, int numCols, bool isMaxtree)
        : MorphologicalTree(static_cast<int*>(input.request().ptr), numRows, numCols, isMaxtree) { }

	MorphologicalTreePybind(py::array_t<int> input, int numRows, int numCols, bool isMaxtree, double radiusOfAdjacencyRelation)
        : MorphologicalTree(static_cast<int*>(input.request().ptr), numRows, numCols, isMaxtree, radiusOfAdjacencyRelation) { }
   

    /*
    py::array_t<int> getOrderedPixels(){
        int n = this->numRows * this->numCols;
        return PybindUtils::toNumpy(this->orderedPixels, n);
    }

    py::array_t<int> getParent(){
        int n = this->numRows * this->numCols;
        return PybindUtils::toNumpy(this->parent, n);
    }*/

    py::array_t<int> getImageAferPruning(NodeMTPtr node){
        int n = this->numRows * this->numCols;
        int* imgOut = MorphologicalTree::getImageAferPruning(node); // Chamar método da superclasse
        return PybindUtils::toNumpy(imgOut, n);
    }

    py::array_t<int> reconstructionImage(){
        int n = this->numRows * this->numCols;
        int* imgOut = new int[n];
        MorphologicalTree::reconstruction(this->root, imgOut);
        return PybindUtils::toNumpy(imgOut, n);
    }

    /*static py::array_t<int> computerParent(py::array_t<int> input, int numRows, int numCols, bool isMaxtree){
		auto buf_input = input.request();
		int* img = (int *) buf_input.ptr;
		ComponentTree tree(img, numRows, numCols, isMaxtree);
		return PybindUtils::toNumpy(tree.getParent(), numRows * numCols);;
	}*/

};



#endif