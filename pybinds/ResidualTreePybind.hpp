#pragma once

#include <list>

#include "../mmcfilters/trees/NodeMT.hpp"
#include "../mmcfilters/trees/NodeRes.hpp"
#include "../mmcfilters/filters/AttributeOpeningPrimitivesFamily.hpp"
#include "../mmcfilters/trees/ResidualTree.hpp"
#include "../mmcfilters/utils/Common.hpp"

#include "PybindUtils.hpp"
namespace mmcfilters {


/**
 * @brief Interface Pybind para reconstruções baseadas em árvores residuais.
 */
class ResidualTreePybind: public ResidualTree{

    public:
    using ResidualTree::ResidualTree;

        ResidualTreePybind(AttributeOpeningPrimitivesFamilyPtr primitivesFamily): ResidualTree(primitivesFamily){}

        py::array_t<uint8_t> reconstruction(){
            return PybindUtils::toNumpy(ResidualTree::reconstruction());
        }

        py::array_t<uint8_t> filtering(std::vector<uint8_t>& criterion){
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

} // namespace mmcfilters
