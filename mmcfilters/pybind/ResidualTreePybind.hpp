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

        ResidualTreePybind(AttributeOpeningPrimitivesFamilyPtr primitivesFamily): ResidualTree(primitivesFamily){}

        py::array_t<uint8_t> reconstruction(){
            return PybindUtils::toNumpy(ResidualTree::reconstruction());
        }

        py::array_t<uint8_t> filtering(std::vector<bool> criterion){
            return PybindUtils::toNumpy(ResidualTree::filtering(criterion));
        }

        py::array_t<uint8_t> getMaxConstrastImage(){
            return PybindUtils::toNumpy(ResidualTree::getMaxConstrastImage());
        }       

        py::array_t<int32_t> getAssociatedImage(){
            auto imgOut = ResidualTree::getAssociatedImage();
            return PybindUtils::toNumpyInt(imgOut->rawData(), imgOut->getSize());
        }

        py::array_t<uint8_t> getAssociatedColoredImage(){
            return PybindUtils::toNumpy(ResidualTree::getAssociatedColorImage());
        }

        py::array_t<uint8_t> getNegativeResidues(){
            return PybindUtils::toNumpy(ResidualTree::getNegativeResidues());
        }

        py::array_t<uint8_t> getPositiveResidues(){
            return PybindUtils::toNumpy(ResidualTree::getPositiveResidues());
        }
};

#endif