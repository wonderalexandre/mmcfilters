
#include "include/NodeCT.hpp"
#include "include/AdjacencyRelation.hpp"

#include "pybind/AttributeComputedIncrementallyPybind.hpp"
#include "pybind/ComponentTreePybind.hpp"
#include "pybind/AttributeFiltersPybind.hpp"
#include "pybind/UltimateAttributeOpeningPybind.hpp"
#include "pybind/AttributeOpeningPrimitivesFamilyPybind.hpp"
#include "pybind/ResidualTreePybind.hpp"


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <iterator>
#include <utility>


namespace py = pybind11;


void init_NodeCT(py::module &m){
    py::class_<NodeCT>(m, "NodeCT")
		.def(py::init<>())
		.def_property_readonly("id", &NodeCT::getIndex )
        .def_property_readonly("rep", &NodeCT::getRep )
		.def_property_readonly("cnps", &NodeCT::getCNPs )
		.def_property_readonly("level", &NodeCT::getLevel )
		.def_property_readonly("children", &NodeCT::getChildren )
		.def_property_readonly("parent", &NodeCT::getParent )
        .def_property_readonly("areaCC", &NodeCT::getAreaCC )
        .def_property_readonly("numDescendants", &NodeCT::getNumDescendants )
        .def_property_readonly("isMaxtree", &NodeCT::isMaxtreeNode )
        .def_property_readonly("numSiblings", &NodeCT::getNumSiblings )
        .def_property_readonly("residue", &NodeCT::getResidue ) 
        .def("pixelsOfCC",&NodeCT::getPixelsOfCC )
        .def("nodesOfPathToRoot",&NodeCT::getNodesOfPathToRoot )
        .def("nodesDescendants",&NodeCT::getNodesDescendants );
        
    py::class_<NodeCT::IteratorPixelsOfCC>(m, "IteratorPixelsOfCC")
		.def(py::init<NodeCT *, int>())
		.def_property_readonly("begin", &NodeCT::IteratorPixelsOfCC::begin )
        .def_property_readonly("end", &NodeCT::IteratorPixelsOfCC::end )
        .def("__iter__", [](NodeCT::IteratorPixelsOfCC &iter) {
            return py::make_iterator(iter.begin(), iter.end());
            }, py::keep_alive<0, 1>()); /* Keep vector alive while iterator is used */

    py::class_<NodeCT::IteratorNodesOfPathToRoot>(m, "IteratorNodesOfPathToRoot")
		.def(py::init<NodeCT *>())
		.def_property_readonly("begin", &NodeCT::IteratorNodesOfPathToRoot::begin )
        .def_property_readonly("end", &NodeCT::IteratorNodesOfPathToRoot::end )
        .def("__iter__", [](NodeCT::IteratorNodesOfPathToRoot &iter) {
            return py::make_iterator(iter.begin(), iter.end());
            }, py::keep_alive<0, 1>()); /* Keep vector alive while iterator is used */
            
    py::class_<NodeCT::IteratorNodesDescendants>(m, "IteratorNodesDescendants")
		.def(py::init<NodeCT *, int>())
		.def_property_readonly("begin", &NodeCT::IteratorNodesDescendants::begin )
        .def_property_readonly("end", &NodeCT::IteratorNodesDescendants::end )
        .def("__iter__", [](NodeCT::IteratorNodesDescendants &iter) {
            return py::make_iterator(iter.begin(), iter.end());
            }, py::keep_alive<0, 1>()); /* Keep vector alive while iterator is used */
}



void init_ComponentTree(py::module &m){
      py::class_<ComponentTreePybind>(m, "ComponentTree")
        .def(py::init<py::array_t<int>, int, int, bool, double>())
        .def(py::init<py::array_t<int>, int, int, bool>())
        .def(py::init<py::array_t<int>, int, int>())
        .def("reconstructionImage", &ComponentTreePybind::reconstructionImage )
		.def_property_readonly("numNodes", &ComponentTreePybind::getNumNodes )
        .def_property_readonly("listNodes", &ComponentTreePybind::getListNodes )
        .def_property_readonly("root", &ComponentTreePybind::getRoot )
        //.def_static("computerParent", &ComponentTreePybind::computerParent)
		//.def_property_readonly("parent", &ComponentTreePybind::getParent )
        //.def_property_readonly("orderedPixels", &ComponentTreePybind::getOrderedPixels )
        .def_property_readonly("treeType", &ComponentTreePybind::getTreeType)
        .def("getSC", &ComponentTreePybind::getSC );
        
        
        //.def("prunningMin", py::overload_cast<py::array_t<double> &, double>(&ComponentTree::prunningMin))
        //.def("prunningMin", &ComponentTree::prunningMin)
        //.def("computerArea", &ComponentTree::computerArea)
        //.def("prunningMin", py::overload_cast<py::array_t<double> &, double>(&ComponentTree::prunningMin))
}


void init_AttributeComputedIncrementally(py::module &m){
    	py::class_<AttributeComputedIncrementallyPybind>(m, "Attribute")
        //.def_static("computerAttribute", &AttributeComputedIncrementally::computerAttribute)
        .def_static("computerBasicAttributes", &AttributeComputedIncrementallyPybind::computerBasicAttributes)
        .def_static("computerArea", &AttributeComputedIncrementallyPybind::computerArea);
}

void init_AttributeFilters(py::module &m){
    py::class_<AttributeFiltersPybind>(m, "AttributeFilters")
    .def(py::init<ComponentTreePybind *>())
    .def("filteringMin", py::overload_cast<py::array_t<float> &, float>(&AttributeFiltersPybind::filteringByPruningMin))
    .def("filteringMin", py::overload_cast<std::vector<bool>>(&AttributeFiltersPybind::filteringByPruningMin))
    .def("filteringMax", py::overload_cast<std::vector<bool>>(&AttributeFiltersPybind::filteringByPruningMax))
    .def("filteringDirectRule", py::overload_cast<std::vector<bool>>(&AttributeFiltersPybind::filteringByDirectRule))
    .def("filteringSubtractiveRule", py::overload_cast<std::vector<bool>>(&AttributeFiltersPybind::filteringBySubtractiveRule))
    .def("filteringSubtractiveScoreRule", py::overload_cast<std::vector<float>>(&AttributeFiltersPybind::filteringBySubtractiveScoreRule))
    .def("filteringMax", py::overload_cast<py::array_t<float> &, float>(&AttributeFiltersPybind::filteringByPruningMax));   
}


void init_AdjacencyRelation(py::module &m){
    	py::class_<AdjacencyRelation>(m, "AdjacencyRelation")
        .def(py::init<int, int, double>())
        .def_property_readonly("size", &AdjacencyRelation::getSize )
        .def("getAdjPixels", py::overload_cast<int, int>( &AdjacencyRelation::getAdjPixels ));
}


void init_UltimateAttributeOpening(py::module &m){
    	py::class_<UltimateAttributeOpeningPybind>(m, "UltimateAttributeOpening")
        .def(py::init<ComponentTreePybind *, std::vector<float>>())
        .def("execute", py::overload_cast<int>(&UltimateAttributeOpeningPybind::execute))
        .def("executeWithMSER", &UltimateAttributeOpeningPybind::executeWithMSER)
        .def("getMaxConstrastImage", &UltimateAttributeOpeningPybind::getMaxConstrastImage)
        .def("getAssociatedImage", &UltimateAttributeOpeningPybind::getAssociatedImage)
        .def("getAssociatedColoredImage", &UltimateAttributeOpeningPybind::getAssociatedColorImage);
}

void init_ResidualTree(py::module &m){
    	py::class_<ResidualTreePybind>(m, "ResidualTree")
        .def(py::init<AttributeOpeningPrimitivesFamilyPybind *>())
        .def("reconstruction", &ResidualTreePybind::reconstruction)
        .def("filtering", &ResidualTreePybind::filtering)
        .def("computerMaximumResidues", &ResidualTreePybind::computerMaximumResidues)
        .def("getMaxConstrastImage", &ResidualTreePybind::getMaxConstrastImage)
        .def("getAssociatedImage", &ResidualTreePybind::getAssociatedImage)
        .def("getAssociatedColoredImage", &ResidualTreePybind::getAssociatedColoredImage)
        .def("getNegativeResidues", &ResidualTreePybind::getNegativeResidues)
        .def("getPositiveResidues", &ResidualTreePybind::getPositiveResidues);

}

void init_AttributeOpeningPrimitivesFamily(py::module &m){
    	py::class_<AttributeOpeningPrimitivesFamilyPybind>(m, "AttributeOpeningPrimitivesFamily")
        .def(py::init<ComponentTreePybind *, py::array_t<float>, float>())
        .def(py::init<ComponentTreePybind *, py::array_t<float>, float, int>())
        .def_property_readonly("numPrimitives", &AttributeOpeningPrimitivesFamilyPybind::getNumPrimitives)
        .def("getPrimitive", &AttributeOpeningPrimitivesFamilyPybind::getPrimitive)
        .def_property_readonly("restOfImage", &AttributeOpeningPrimitivesFamilyPybind::getRestOfNumpyImage)
        .def("getNodesWithMaximumCriterium", &AttributeOpeningPrimitivesFamilyPybind::getNodesWithMaximumCriterium)
        .def("getThresholdsPrimitive", &AttributeOpeningPrimitivesFamilyPybind::getThresholdsPrimitive);

}


PYBIND11_MODULE(mmcfilters, m) {
    // Optional docstring
    m.doc() = "A simple library for connected filters based on morphological trees";
    
    init_NodeCT(m);
    init_ComponentTree(m);
    init_AttributeComputedIncrementally(m);
    init_AttributeFilters(m);
    init_AdjacencyRelation(m);

    init_UltimateAttributeOpening(m);
    init_ResidualTree(m);
    init_AttributeOpeningPrimitivesFamily(m);

}
