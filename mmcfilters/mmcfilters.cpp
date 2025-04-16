
#include "include/NodeMT.hpp"
#include "include/AdjacencyRelation.hpp"
#include "include/Common.hpp"

#include "pybind/AttributeComputedIncrementallyPybind.hpp"
#include "pybind/MorphologicalTreePybind.hpp"
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
using namespace pybind11::literals;

void init_NodeCT(py::module &m){
    py::class_<NodeMT, std::shared_ptr<NodeMT>>(m, "NodeMT")
		.def(py::init<>())
		.def_property_readonly("id", &NodeMT::getIndex )
        .def("__str__", [](NodeMT &node) {
            std::ostringstream oss;
            oss << "NodeCT(id=" << node.getIndex() 
                << ", level=" << node.getLevel() 
                << ", numCNPs=" << node.getCNPs().size() 
                << ", area=" << node.getAreaCC(); 
            return oss.str();
        })
        .def("__repr__", [](NodeMT &node) { 
            std::ostringstream oss;
            oss << "NodeCT(id=" << node.getIndex() << ", level=" << node.getLevel() << ")";
            return oss.str();
        })
 		.def_property_readonly("cnps", &NodeMT::getCNPs )
		.def_property_readonly("level", &NodeMT::getLevel )
		.def_property_readonly("children", &NodeMT::getChildren )
		.def_property_readonly("parent", &NodeMT::getParent )
        .def_property_readonly("areaCC", &NodeMT::getAreaCC )
        .def_property_readonly("numDescendants", &NodeMT::getNumDescendants )
        .def_property_readonly("isMaxtree", &NodeMT::isMaxtreeNode )
        .def_property_readonly("numSiblings", &NodeMT::getNumSiblings )
        .def_property_readonly("residue", &NodeMT::getResidue ) 
        .def("pixelsOfCC",&NodeMT::getPixelsOfCC )
        .def("nodesOfPathToRoot",&NodeMT::getNodesOfPathToRoot )
        .def("nodesDescendants",&NodeMT::getNodesDescendants )
        .def("bfsTraversal", &NodeMT::getIteratorBreadthFirstTraversal)
        .def("postOrderTraversal", &NodeMT::getIteratorPostOrderTraversal)
        .def("recNode", [](NodeMTPtr node) {
            return MorphologicalTreePybind::recNode(node);
        });

        
}


void init_NodeCT_Iterators(py::module &m) {

    py::class_<typename NodeMT::IteratorPixelsOfCC>(m, "IteratorPixelsOfCC")
        .def(py::init<std::shared_ptr<NodeMT>, int>())
        .def("__iter__", [](typename NodeMT::IteratorPixelsOfCC &iter) {
            return py::make_iterator(iter.begin(), iter.end());
        }, py::keep_alive<0, 1>());


    py::class_<typename NodeMT::IteratorNodesOfPathToRoot>(m, "IteratorNodesOfPathToRoot")
        .def(py::init<std::shared_ptr<NodeMT>>())
        .def("__iter__", [](typename NodeMT::IteratorNodesOfPathToRoot &iter) {
            return py::make_iterator(iter.begin(), iter.end());
        }, py::keep_alive<0, 1>());

    py::class_<typename NodeMT::IteratorPostOrderTraversal>(m, "IteratorPostOrderTraversal")
        .def(py::init<std::shared_ptr<NodeMT>>())
        .def("__iter__", [](typename NodeMT::IteratorPostOrderTraversal &iter) {
            return py::make_iterator(iter.begin(), iter.end());
        }, py::keep_alive<0, 1>());

    py::class_<typename NodeMT::IteratorBreadthFirstTraversal>(m, "IteratorBreadthFirstTraversal")
        .def(py::init<std::shared_ptr<NodeMT>>())
        .def("__iter__", [](typename NodeMT::IteratorBreadthFirstTraversal &iter) {
            return py::make_iterator(iter.begin(), iter.end());
        }, py::keep_alive<0, 1>());

         
    py::class_<typename NodeMT::IteratorNodesDescendants>(m, "IteratorNodesDescendants")
    .def(py::init<std::shared_ptr<NodeMT>, int>())
    .def("__iter__", [](NodeMT::IteratorNodesDescendants &iter) {
        return py::make_iterator(iter.begin(), iter.end());
        }, py::keep_alive<0, 1>()); /* Keep vector alive while iterator is used */

}

void init_MorphologicalTree(py::module &m){
      py::class_<MorphologicalTreePybind, std::shared_ptr<MorphologicalTreePybind>>(m, "MorphologicalTree")
        .def(py::init<py::array_t<int>, int, int, bool, double>(),
            "input"_a, "rows"_a, "cols"_a, "isMaxtree"_a, "radius"_a = 1.5)
        .def(py::init<py::array_t<int>, int, int, std::string>(),
            "input"_a, "rows"_a, "cols"_a, "ToSInperpolation"_a = "self-dual")
        .def("reconstructionImage", &MorphologicalTreePybind::reconstructionImage )
		.def_property_readonly("numNodes", &MorphologicalTreePybind::getNumNodes )
        .def_property_readonly("listNodes", &MorphologicalTreePybind::getIndexNode )
        .def_property_readonly("root", &MorphologicalTreePybind::getRoot )
        .def_property_readonly("depth", &MorphologicalTreePybind::getDepth )
        //.def_static("computerParent", &ComponentTreePybind::computerParent)
		//.def_property_readonly("parent", &ComponentTreePybind::getParent )
        //.def_property_readonly("orderedPixels", &ComponentTreePybind::getOrderedPixels )
        .def_property_readonly("treeType", &MorphologicalTreePybind::getTreeType)
        .def("getSC", &MorphologicalTreePybind::getSC );
        
        
        //.def("prunningMin", py::overload_cast<py::array_t<double> &, double>(&ComponentTree::prunningMin))
        //.def("prunningMin", &ComponentTree::prunningMin)
        //.def("computerArea", &ComponentTree::computerArea)
        //.def("prunningMin", py::overload_cast<py::array_t<double> &, double>(&ComponentTree::prunningMin))
}


void init_AttributeComputedIncrementally(py::module &m){
    	py::class_<AttributeComputedIncrementallyPybind>(m, "Attribute")
        .def_static("computerAttribute", static_cast<void(*)(NodeMTPtr, 
                                                             std::function<void(NodeMTPtr)>, 
                                                             std::function<void(NodeMTPtr, NodeMTPtr)>, 
                                                             std::function<void(NodeMTPtr)>)>(&AttributeComputedIncrementally::computerAttribute))
        .def_static("computerBasicAttributes", &AttributeComputedIncrementallyPybind::computerBasicAttributes)
        .def_static("extractCountors", &AttributeComputedIncrementallyPybind::extractCountors)
        .def_static("computerArea", &AttributeComputedIncrementallyPybind::computerArea);
}

void init_AttributeFilters(py::module &m){
    py::class_<AttributeFiltersPybind>(m, "AttributeFilters")
    .def(py::init<MorphologicalTreePybindPtr>())
    .def("filteringMin", py::overload_cast<py::array_t<float> &, float>(&AttributeFiltersPybind::filteringByPruningMin))
    .def("filteringMin", py::overload_cast<std::vector<bool>&>(&AttributeFiltersPybind::filteringByPruningMin))
    .def("filteringMax", py::overload_cast<std::vector<bool>&>(&AttributeFiltersPybind::filteringByPruningMax))
    .def("filteringDirectRule", py::overload_cast<std::vector<bool>&>(&AttributeFiltersPybind::filteringByDirectRule))
    .def("filteringSubtractiveRule", py::overload_cast<std::vector<bool>&>(&AttributeFiltersPybind::filteringBySubtractiveRule))
    .def("filteringSubtractiveScoreRule", py::overload_cast<std::vector<float>&>(&AttributeFiltersPybind::filteringBySubtractiveScoreRule))
    .def("filteringMax", py::overload_cast<py::array_t<float> &, float>(&AttributeFiltersPybind::filteringByPruningMax))
    .def("getAdaptativeCriterion", &AttributeFiltersPybind::getAdaptativeCriterion);   
}


void init_AdjacencyRelation(py::module &m){
    	py::class_<AdjacencyRelation>(m, "AdjacencyRelation")
        .def(py::init<int, int, double>())
        .def_property_readonly("size", &AdjacencyRelation::getSize )
        .def("getAdjPixels", py::overload_cast<int, int>( &AdjacencyRelation::getAdjPixels ));
}


void init_UltimateAttributeOpening(py::module &m){
    	py::class_<UltimateAttributeOpeningPybind>(m, "UltimateAttributeOpening")
        .def(py::init<MorphologicalTreePybindPtr, std::vector<float>>())
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
        .def(py::init<MorphologicalTreePybindPtr, py::array_t<float>, float>())
        .def(py::init<MorphologicalTreePybindPtr, py::array_t<float>, float, int>())
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
    init_NodeCT_Iterators(m);
    init_MorphologicalTree(m);
    init_AttributeComputedIncrementally(m);
    init_AttributeFilters(m);
    init_AdjacencyRelation(m);

    init_UltimateAttributeOpening(m);
    init_ResidualTree(m);
    init_AttributeOpeningPrimitivesFamily(m);

}
