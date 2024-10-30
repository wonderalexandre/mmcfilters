#include <list>

#include "../include/NodeCT.hpp"
#include "../include/NodeRes.hpp"
#include "../include/AttributeOpeningPrimitivesFamily.hpp"
#include "../include/ResidualTree.hpp"

#include "../pybind/PybindUtils.hpp"

#ifndef RESIDUAL_TREE_PYBIND_H
#define RESIDUAL_TREE_PYBIND_H


class ResidualTreePybind: public ResidualTree{

    public:
    using ResidualTree::ResidualTree;

        ResidualTreePybind(AttributeOpeningPrimitivesFamily* primitivesFamily): ResidualTree(primitivesFamily){}

        py::array_t<int> reconstruction(){
            int* imgOutput = ResidualTree::reconstruction();
            return PybindUtils::toNumpy(imgOutput, this->tree->getNumRowsOfImage() * this->tree->getNumColsOfImage());
        }

        py::array_t<int> filtering(std::vector<bool> criterion){
            int n = this->tree->getNumRowsOfImage() * this->tree->getNumColsOfImage();
            int* imgOutput = new int[n];

            ResidualTree::filtering(criterion, imgOutput);
            return PybindUtils::toNumpy(imgOutput, n);
        }

        py::array_t<int> getMaxConstrastImage(){
            return PybindUtils::toNumpy(ResidualTree::getMaxConstrastImage(), this->tree->getNumColsOfImage() * this->tree->getNumRowsOfImage());
        }       

        py::array_t<int> getAssociatedImage(){
            return PybindUtils::toNumpy(ResidualTree::getAssociatedImage(), this->tree->getNumColsOfImage() * this->tree->getNumRowsOfImage());
        }

        py::array_t<int> getAssociatedColoredImage(){
            return PybindUtils::toNumpy(ResidualTree::getAssociatedColorImage(), this->tree->getNumColsOfImage() * this->tree->getNumRowsOfImage() * 3);
        }

        py::array_t<int> getNegativeResidues(){
            return PybindUtils::toNumpy(ResidualTree::getNegativeResidues(),  this->tree->getNumColsOfImage() * this->tree->getNumRowsOfImage());
        }

        py::array_t<int> getPositiveResidues(){
            return PybindUtils::toNumpy(ResidualTree::getPositiveResidues(),  this->tree->getNumColsOfImage() * this->tree->getNumRowsOfImage());
        }
};

#endif