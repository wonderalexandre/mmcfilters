
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
        .def_property_readonly("area", &NodeMT::getAreaCC )
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
        .def_property_readonly("treeType", &MorphologicalTreePybind::getTreeType)
        .def_property_readonly("numRows", &MorphologicalTreePybind::getNumRowsOfImage )
        .def_property_readonly("numCols", &MorphologicalTreePybind::getNumColsOfImage )
        .def_property_readonly("depth", &MorphologicalTreePybind::getDepth )
        .def_property_readonly("leaves", &MorphologicalTreePybind::getLeaves )
        .def("getSC", &MorphologicalTreePybind::getSC );

        
        
}


void init_AttributeComputedIncrementally(py::module &m){
        auto cls = py::class_<AttributeComputedIncrementallyPybind>(m, "Attribute")	
        .def_static("computerAttribute", static_cast<void(*)(NodeMTPtr, 
                                                             std::function<void(NodeMTPtr)>, 
                                                             std::function<void(NodeMTPtr, NodeMTPtr)>, 
                                                             std::function<void(NodeMTPtr)>)>(&AttributeComputedIncrementally::computerAttribute))
        .def_static("computeAttributes", &AttributeComputedIncrementallyPybind::computeAttributesFromList)
        .def_static("computeSingleAttribute", &AttributeComputedIncrementallyPybind::computeSingleAttribute)
        .def_static("computeSingleAttributeWithDelta", &AttributeComputedIncrementallyPybind::computeSingleAttributeWithDelta)
        .def_static("describe", &AttributeComputedIncrementallyPybind::describeAttribute)
        .def_static("extractCountors", &AttributeComputedIncrementallyPybind::extractCompactCountors)
        .def_static("extractCountorsNonCompact", &AttributeComputedIncrementallyPybind::extractCountors)
        .def_static("extractionExtinctionValues", &AttributeComputedIncrementallyPybind::extractionExtinctionValues);
        

        py::class_<ContoursMT, std::shared_ptr<ContoursMT>>(m, "ContoursMT")
            .def("contours", &ContoursMT::contoursLazy);

        py::class_<ContoursMT::ContourPostOrderRange>(m, "ContourPostOrderRange")
            .def("__iter__", [](ContoursMT::ContourPostOrderRange &self) { return self.begin(); });

        py::class_<ContoursMT::ContourPostOrderIterator>(m, "ContourPostOrderIterator")
            .def("__iter__", [](ContoursMT::ContourPostOrderIterator &self) -> ContoursMT::ContourPostOrderIterator& { return self; })
            .def("__next__", [](ContoursMT::ContourPostOrderIterator &self) {
                if (self == ContoursMT::ContourPostOrderIterator(nullptr, nullptr))
                    throw py::stop_iteration();
                auto val = *self;
                ++self;
                return val;
            });

        py::enum_<AttributeGroup>(cls, "Group")
            .value("ALL", AttributeGroup::ALL)
            .value("GEOMETRIC", AttributeGroup::GEOMETRIC)
            .value("BOUNDING_BOX", AttributeGroup::BOUNDING_BOX)
            .value("CENTRAL_MOMENTS", AttributeGroup::CENTRAL_MOMENTS)
            .value("HU_MOMENTS", AttributeGroup::HU_MOMENTS)
            .value("MOMENT_BASED", AttributeGroup::MOMENT_BASED)
            .value("TREE_TOPOLOGY", AttributeGroup::TREE_TOPOLOGY)
            .export_values();

         py::enum_<Attribute>(cls, "Type")
            .value("AREA", Attribute::AREA)
            .value("VOLUME", Attribute::VOLUME)
            .value("RELATIVE_VOLUME", Attribute::RELATIVE_VOLUME)
            .value("LEVEL", Attribute::LEVEL)
            .value("GRAY_HEIGHT", Attribute::GRAY_HEIGHT)
            .value("MEAN_LEVEL", Attribute::MEAN_LEVEL)
            .value("VARIANCE_LEVEL", Attribute::VARIANCE_LEVEL)
            .value("BOX_WIDTH", Attribute::BOX_WIDTH)
            .value("BOX_HEIGHT", Attribute::BOX_HEIGHT)
            .value("RECTANGULARITY", Attribute::RECTANGULARITY)
            .value("DIAGONAL_LENGTH", Attribute::DIAGONAL_LENGTH)
            .value("BOX_COL_MIN", Attribute::BOX_COL_MIN)
            .value("BOX_COL_MAX", Attribute::BOX_COL_MAX)
            .value("BOX_ROW_MIN", Attribute::BOX_ROW_MIN)
            .value("BOX_ROW_MAX", Attribute::BOX_ROW_MAX)
            .value("RATIO_WH", Attribute::RATIO_WH)
            .value("CENTRAL_MOMENT_20", Attribute::CENTRAL_MOMENT_20)
            .value("CENTRAL_MOMENT_02", Attribute::CENTRAL_MOMENT_02)
            .value("CENTRAL_MOMENT_11", Attribute::CENTRAL_MOMENT_11)
            .value("CENTRAL_MOMENT_30", Attribute::CENTRAL_MOMENT_30)
            .value("CENTRAL_MOMENT_03", Attribute::CENTRAL_MOMENT_03)
            .value("CENTRAL_MOMENT_21", Attribute::CENTRAL_MOMENT_21)
            .value("CENTRAL_MOMENT_12", Attribute::CENTRAL_MOMENT_12)
            .value("AXIS_ORIENTATION", Attribute::AXIS_ORIENTATION)
            .value("LENGTH_MAJOR_AXIS", Attribute::LENGTH_MAJOR_AXIS)
            .value("LENGTH_MINOR_AXIS", Attribute::LENGTH_MINOR_AXIS)
            .value("ECCENTRICITY", Attribute::ECCENTRICITY)
            .value("COMPACTNESS", Attribute::COMPACTNESS)
            .value("INERTIA", Attribute::INERTIA)
            .value("HU_MOMENT_1", Attribute::HU_MOMENT_1)
            .value("HU_MOMENT_2", Attribute::HU_MOMENT_2)
            .value("HU_MOMENT_3", Attribute::HU_MOMENT_3)
            .value("HU_MOMENT_4", Attribute::HU_MOMENT_4)
            .value("HU_MOMENT_5", Attribute::HU_MOMENT_5)
            .value("HU_MOMENT_6", Attribute::HU_MOMENT_6)
            .value("HU_MOMENT_7", Attribute::HU_MOMENT_7)
            .value("HEIGHT_NODE", Attribute::HEIGHT_NODE)
            .value("DEPTH_NODE", Attribute::DEPTH_NODE)
            .value("IS_LEAF_NODE", Attribute::IS_LEAF_NODE)
            .value("IS_ROOT_NODE", Attribute::IS_ROOT_NODE)
            .value("NUM_CHILDREN_NODE", Attribute::NUM_CHILDREN_NODE)
            .value("NUM_SIBLINGS_NODE", Attribute::NUM_SIBLINGS_NODE)
            .value("NUM_DESCENDANTS_NODE", Attribute::NUM_DESCENDANTS_NODE)
            .value("NUM_LEAF_DESCENDANTS_NODE", Attribute::NUM_LEAF_DESCENDANTS_NODE)
            .value("LEAF_RATIO_NODE", Attribute::LEAF_RATIO_NODE)
            .value("BALANCE_NODE", Attribute::BALANCE_NODE)
            .value("AVG_CHILD_HEIGHT_NODE", Attribute::AVG_CHILD_HEIGHT_NODE)
            .export_values();
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
    .def("filteringByExtinctionValue", py::overload_cast<py::array_t<float> &, int>(&AttributeFiltersPybind::filteringByExtinctionValue))
    .def("saliencyMapByExtinction", py::overload_cast<py::array_t<float> &, int>(&AttributeFiltersPybind::saliencyMapByExtinction))
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
        .def(py::init<MorphologicalTreePybindPtr, py::array_t<float>&>())
        .def("execute", py::overload_cast<int>(&UltimateAttributeOpeningPybind::execute))
        .def("executeWithMSER", &UltimateAttributeOpeningPybind::executeWithMSER)
        .def("getMaxConstrastImage", &UltimateAttributeOpeningPybind::getMaxConstrastImage)
        .def("getAssociatedImage", &UltimateAttributeOpeningPybind::getAssociatedImage)
        .def("getAssociatedColoredImage", &UltimateAttributeOpeningPybind::getAssociatedColorImage);
}

void init_ResidualTree(py::module &m){
    	py::class_<ResidualTreePybind>(m, "ResidualTree")
        .def(py::init<std::shared_ptr<AttributeOpeningPrimitivesFamilyPybind>>())
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
    	py::class_<AttributeOpeningPrimitivesFamilyPybind, std::shared_ptr<AttributeOpeningPrimitivesFamilyPybind>>(m, "AttributeOpeningPrimitivesFamily")
        .def(py::init<MorphologicalTreePybindPtr, py::array_t<float>&, float>())
        .def(py::init<MorphologicalTreePybindPtr, py::array_t<float>&, float, int>())
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
