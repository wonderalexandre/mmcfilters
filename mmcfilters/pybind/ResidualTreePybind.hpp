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

        py::array_t<PixelType> reconstruction(){
            return PybindUtils::toNumpy(ResidualTree::reconstruction()->rawData(), this->tree->getNumRowsOfImage() * this->tree->getNumColsOfImage());
        }

        py::array_t<PixelType> filtering(std::vector<bool> criterion){
            return PybindUtils::toNumpy(ResidualTree::filtering(criterion)->rawData(), this->tree->getNumRowsOfImage() * this->tree->getNumColsOfImage());
        }

        py::array_t<PixelType> getMaxConstrastImage(){
            return PybindUtils::toNumpy(ResidualTree::getMaxConstrastImage()->rawData(), this->tree->getNumColsOfImage() * this->tree->getNumRowsOfImage());
        }       

        py::array_t<int> getAssociatedImage(){
            return PybindUtils::toNumpyInt(ResidualTree::getAssociatedImage(), this->tree->getNumColsOfImage() * this->tree->getNumRowsOfImage());
        }

        py::array_t<PixelType> getAssociatedColoredImage(){
            return PybindUtils::toNumpy(ResidualTree::getAssociatedColorImage(), this->tree->getNumColsOfImage() * this->tree->getNumRowsOfImage() * 3);
        }

        py::array_t<PixelType> getNegativeResidues(){
            return PybindUtils::toNumpy(ResidualTree::getNegativeResidues()->rawData(),  this->tree->getNumColsOfImage() * this->tree->getNumRowsOfImage());
        }

        py::array_t<PixelType> getPositiveResidues(){
            return PybindUtils::toNumpy(ResidualTree::getPositiveResidues()->rawData(),  this->tree->getNumColsOfImage() * this->tree->getNumRowsOfImage());
        }
};

#endif