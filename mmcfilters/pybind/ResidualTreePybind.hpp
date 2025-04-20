#include <list>

#include "../include/NodeMT.hpp"
#include "../include/NodeRes.hpp"
#include "../include/AttributeOpeningPrimitivesFamily.hpp"
#include "../include/ResidualTree.hpp"
#include "../include/Common.hpp"

#include "../pybind/PybindUtils.hpp"

#ifndef RESIDUAL_TREE_PYBIND_H
#define RESIDUAL_TREE_PYBIND_H


class ResidualTreePybind: public ResidualTree{

    public:
    using ResidualTree::ResidualTree;

        ResidualTreePybind(AttributeOpeningPrimitivesFamily* primitivesFamily): ResidualTree(primitivesFamily){}

        py::array_t<PixelValueType> reconstruction(){
            PixelValueType* imgOutput = ResidualTree::reconstruction()->rawData();
            return PybindUtils::toNumpy(imgOutput, this->tree->getNumRowsOfImage() * this->tree->getNumColsOfImage());
        }

        py::array_t<PixelValueType> filtering(std::vector<bool> criterion){
            int n = this->tree->getNumRowsOfImage() * this->tree->getNumColsOfImage();
            PixelValueType* imgOutput = ResidualTree::filtering(criterion)->rawData();
            return PybindUtils::toNumpy(imgOutput, n);
        }

        py::array_t<int> getMaxConstrastImage(){
            return PybindUtils::toNumpy(ResidualTree::getMaxConstrastImage()->rawData(), this->tree->getNumColsOfImage() * this->tree->getNumRowsOfImage());
        }       

        py::array_t<int> getAssociatedImage(){
            return PybindUtils::toNumpyInt(ResidualTree::getAssociatedImage(), this->tree->getNumColsOfImage() * this->tree->getNumRowsOfImage());
        }

        py::array_t<int> getAssociatedColoredImage(){
            return PybindUtils::toNumpyInt(ResidualTree::getAssociatedColorImage(), this->tree->getNumColsOfImage() * this->tree->getNumRowsOfImage() * 3);
        }

        py::array_t<int> getNegativeResidues(){
            return PybindUtils::toNumpy(ResidualTree::getNegativeResidues()->rawData(),  this->tree->getNumColsOfImage() * this->tree->getNumRowsOfImage());
        }

        py::array_t<int> getPositiveResidues(){
            return PybindUtils::toNumpy(ResidualTree::getPositiveResidues()->rawData(),  this->tree->getNumColsOfImage() * this->tree->getNumRowsOfImage());
        }
};

#endif