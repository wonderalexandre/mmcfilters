
#include "../mmcfilters/trees/NodeMT.hpp"
#include "../mmcfilters/utils/AdjacencyRelation.hpp"
#include "../mmcfilters/utils/Common.hpp"
#include "../mmcfilters/contours/ContoursComputedIncrementally.hpp"

#include "AttributeComputedIncrementallyPybind.hpp"
#include "ContoursComputedIncrementallyPybind.hpp"
#include "MorphologicalTreePybind.hpp"
#include "ExtinctionValuesPybind.hpp"
#include "AttributeFiltersPybind.hpp"
#include "UltimateAttributeOpeningPybind.hpp"
#include "AttributeOpeningPrimitivesFamilyPybind.hpp"
#include "ResidualTreePybind.hpp"


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>

using namespace mmcfilters;

#include <optional>
#include <sstream>
#include <utility>


namespace py = pybind11;
using namespace pybind11::literals;

void init_NodeCT(py::module &m){
    py::class_<NodeMT>(m, "NodeMT", py::module_local(false))
        .def(py::init<>())
        .def("__bool__", [](const NodeMT& node) { return static_cast<bool>(node); })
        .def_property_readonly("id", &NodeMT::getIndex)
        .def("__str__", [](const NodeMT &node) {
            std::ostringstream oss;
            oss << "NodeMT(id=" << node.getIndex()
                << ", level=" << node.getLevel()
                << ", numCNPs=" << node.getNumCNPs()
                << ", area=" << node.getArea() << ")";
            return oss.str();
        })
        .def("__repr__", [](const NodeMT &node) {
            std::ostringstream oss;
            oss << "NodeMT(id=" << node.getIndex() << ", level=" << node.getLevel() << ")";
            return oss.str();
        })
        .def_property_readonly("level", &NodeMT::getLevel)
        .def_property_readonly("area", &NodeMT::getArea)
        .def_property_readonly("repNode", &NodeMT::getRepNode)
        .def_property_readonly("numDescendants", &NodeMT::getNumDescendants)
        .def_property_readonly("isLeaf", &NodeMT::isLeaf)
        .def_property_readonly("residue", &NodeMT::getResidue)
        .def_property_readonly("isMaxtree", &NodeMT::isMaxtreeNode)
        .def_property_readonly("numSiblings", &NodeMT::getNumSiblings)
        .def_property_readonly("residue", &NodeMT::getResidue)
        .def_property_readonly("numCNPs", &NodeMT::getNumCNPs)
        .def_property_readonly("cnps", [](const NodeMT &node) {
            py::list result;
            for (int value : node.getCNPs()) {
                result.append(value);
            }
            return result;
        })
        .def_property_readonly("children", [](NodeMT &node) {
            py::list children;
            for (auto child : node.getChildren()) {
                if (child) {
                    children.append(child);
                }
            }
            return children;
        })
        .def_property_readonly("parent", [](NodeMT &node) -> py::object {
            NodeMT parent = node.getParent();
            if (!parent) {
                return py::none();
            }
            return py::cast(parent);
        })
        .def("pixelsOfCC", [](const NodeMT &node) {
            py::list pixels;
            for (int p : node.getPixelsOfCC()) {
                pixels.append(p);
            }
            return pixels;
        })
        .def("nodesOfPathToRoot", [](NodeMT &node) {
            py::list nodes;
            auto range = node.getNodesOfPathToRoot();
            for (auto it = range.begin(); it != range.end(); ++it) {
                NodeMT current = *it;
                if (current) {
                    nodes.append(current);
                }
            }
            return nodes;
        })
        .def("nodesDescendants", [](NodeMT &node) {
            auto range = node.getNodesDescendants();
            return py::make_iterator(range.begin(), range.end());
        })
        .def("bfsTraversal", [](NodeMT &node) {
            auto traversal = node.getIteratorBreadthFirstTraversal();
            return py::make_iterator(traversal.begin(), traversal.end());
        }, py::keep_alive<0, 1>())
        .def("postOrderTraversal", [](NodeMT &node) {
            auto traversal = node.getIteratorPostOrderTraversal();
            return py::make_iterator(traversal.begin(), traversal.end());
        }, py::keep_alive<0, 1>())
        .def("recNode", [](NodeMT node) {
            return MorphologicalTreePybind::recNode(node);
        })
        .def_property_readonly("repCNPs", [](NodeMT &node) {
            return MorphologicalTreePybind::repCNPsByFlood(node);
        });
}
void init_MorphologicalTree(py::module &m){
      py::class_<MorphologicalTree, std::shared_ptr<MorphologicalTree>>(m, "MorphologicalTreeBase", py::module_local(false));
      py::class_<MorphologicalTreePybind, std::shared_ptr<MorphologicalTreePybind>>(m, "MorphologicalTree", py::module_local(false))
        .def(py::init<py::array_t<int>, bool, double>(), "input"_a, "isMaxtree"_a, "radius"_a = 1.5)
        .def(py::init<py::array_t<int>,  std::string>(), "input"_a, "ToSInperpolation"_a = "self-dual")
        .def("reconstructionImage", &MorphologicalTreePybind::reconstructionImage )
        .def_property_readonly("listNodes", [](MorphologicalTreePybind &tree) {
            py::list nodes;
            for (NodeId id : tree.getNodeIds()) {
                nodes.append(tree.proxy(id));
            }
            return nodes;
        })
        .def_property_readonly("root", [](MorphologicalTreePybind &tree) {
            return tree.getRoot();
        })
        .def_property_readonly("treeType", [](MorphologicalTreePybind& self) { return self.getTreeType(); })
        .def_property_readonly("numRows", [](MorphologicalTreePybind& self) { return self.getNumRowsOfImage(); })
        .def_property_readonly("numCols", [](MorphologicalTreePybind& self) { return self.getNumColsOfImage(); })
        .def_property_readonly("numNodes", [](MorphologicalTreePybind& self) { return self.getNumNodes(); })
        .def_property_readonly("leaves", [](MorphologicalTreePybind &tree) {
            py::list leaves;
            for (NodeId id : tree.getLeaves()) {
                leaves.append(tree.proxy(id));
            }
            return leaves;
        })
        .def("getSC", [](MorphologicalTreePybind &tree, int pixel) {
            return tree.getSC(pixel);
        })
        .def("nodeIds", [](MorphologicalTreePybind &tree) {
            auto ids = tree.getNodeIds();
            return py::make_iterator(ids.begin(), ids.end());
        }, py::keep_alive<0, 1>())
        .def("getNode", [](MorphologicalTreePybind &tree, int nodeId) {
            return tree.proxy(nodeId);
        })
        .def_static("createFromAttributeMapping", &MorphologicalTreePybind::createTreeFromAttributeMapping );
        
}

void init_ContoursComputedIncrementally(py::module &m){
    // Alias locais
    using Contours     = ContoursComputedIncrementally::IncrementalContours;
    using ContourProxy = Contours::ContourProxy;
    using Range = decltype(std::declval<Contours&>().contoursLazy());
    using Iter  = decltype(std::declval<Range&>().begin());
    
    // Torna ContourProxy iterável em Python
    py::class_<ContourProxy>(m, "ContourProxy", py::module_local(false))
        .def("__iter__", [](const ContourProxy& p) {
            // usa os iteradores já existentes do proxy
            return py::make_iterator(p.begin(), p.end());
        }, py::keep_alive<0, 1>())
        .def("empty", &ContourProxy::empty);

    struct ContoursIterator {
        Contours* owner;
        Range range;
        Iter it, itEnd;
        ContoursIterator(Contours& self)
            : owner(&self), range(self.contoursLazy()), it(range.begin()), itEnd(range.end()) {}
    };

    py::class_<ContoursIterator>(m, "ContoursIterator", py::module_local(false))
        .def(py::init<Contours&>())
        .def("__iter__", [](ContoursIterator& self) -> ContoursIterator& { return self; }, py::return_value_policy::reference_internal)
        .def("__next__", [](ContoursIterator& self) -> py::object {
            if (self.it == self.itEnd) throw py::stop_iteration();
            auto entry = *self.it++;  // <— cópia do par temporário
            auto nodeId = std::get<0>(entry);
            auto proxy  = std::get<1>(entry);  // ContourProxy
            return py::make_tuple(nodeId, proxy);
        });


    py::class_<Contours, std::shared_ptr<Contours>>(m, "Contours", py::module_local(false))
        .def("contours", [](Contours &self) {
            return ContoursIterator(self);
        }, py::keep_alive<0, 1>())
        .def("getContour", &Contours::contour);

    py::class_<ContoursComputedIncrementallyPybind>(m, "ContourComputation", py::module_local(false))
        .def_static("extraction", &ContoursComputedIncrementallyPybind::extraction);
}

void init_AttributeComputedIncrementally(py::module &m){
        auto cls = py::class_<AttributeComputedIncrementallyPybind>(m, "Attribute", py::module_local(false))
        .def_static(
            "computerAttribute",
            [](MorphologicalTreePybind &tree,
               std::function<void(NodeMT)> preProcessing,
               std::function<void(NodeMT, NodeMT)> mergeProcessing,
               std::function<void(NodeMT)> postProcessing,
               std::optional<NodeMT> rootOpt) {
                NodeMT root = rootOpt.value_or(NodeMT());
                if (!root) {
                    root = tree.getRoot();
                }
                NodeId rootId = root ? root.getIndex() : tree.getRoot().getIndex();

                AttributeComputedIncrementally::computerAttribute(
                    &tree,
                    rootId,
                    [&](NodeId nodeId) { preProcessing(tree.proxy(nodeId)); },
                    [&](NodeId parentId, NodeId childId) { mergeProcessing(tree.proxy(parentId), tree.proxy(childId)); },
                    [&](NodeId nodeId) { postProcessing(tree.proxy(nodeId)); }
                );
            },
            py::arg("tree"),
            py::arg("preProcessing"),
            py::arg("mergeProcessing"),
            py::arg("postProcessing"),
            py::arg("root") = std::optional<NodeMT>{}
        )
        .def_static("computeAttributes", &AttributeComputedIncrementallyPybind::computeAttributesFromList)
        .def_static("computeSingleAttribute", &AttributeComputedIncrementallyPybind::computeSingleAttribute)
        .def_static("computeSingleAttributeWithDelta", &AttributeComputedIncrementallyPybind::computeSingleAttributeWithDelta)
        .def_static("describe", &AttributeComputedIncrementallyPybind::describeAttribute)
        .def_static("computerAttributeMapping", &AttributeComputedIncrementallyPybind::computerAttributeMapping);

        py::enum_<AttributeGroup>(cls, "Group", py::module_local(false))
            .value("ALL", AttributeGroup::ALL)
            .value("GEOMETRIC", AttributeGroup::GEOMETRIC)
            .value("BOUNDING_BOX", AttributeGroup::BOUNDING_BOX)
            .value("CENTRAL_MOMENTS", AttributeGroup::CENTRAL_MOMENTS)
            .value("HU_MOMENTS", AttributeGroup::HU_MOMENTS)
            .value("MOMENT_BASED", AttributeGroup::MOMENT_BASED)
            .value("TREE_TOPOLOGY", AttributeGroup::TREE_TOPOLOGY)
            .value("BITQUADS", AttributeGroup::BITQUADS)
            .export_values();

         py::enum_<Attribute>(cls, "Type", py::module_local(false))
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
            .value("CIRCULARITY", Attribute::CIRCULARITY)
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
            .value("BITQUADS_AREA", Attribute::BITQUADS_AREA)
            .value("BITQUADS_NUMBER_EULER", Attribute::BITQUADS_NUMBER_EULER)
            .value("BITQUADS_NUMBER_HOLES", Attribute::BITQUADS_NUMBER_HOLES)
            .value("BITQUADS_PERIMETER", Attribute::BITQUADS_PERIMETER)
            .value("BITQUADS_PERIMETER_CONTINUOUS", Attribute::BITQUADS_PERIMETER_CONTINUOUS)
            .value("BITQUADS_CIRCULARITY", Attribute::BITQUADS_CIRCULARITY)
            .value("BITQUADS_PERIMETER_AVERAGE", Attribute::BITQUADS_PERIMETER_AVERAGE)
            .value("BITQUADS_LENGTH_AVERAGE", Attribute::BITQUADS_LENGTH_AVERAGE)
            .value("BITQUADS_WIDTH_AVERAGE", Attribute::BITQUADS_WIDTH_AVERAGE)
            
            .export_values();
}

void init_AttributeFilters(py::module &m){
    py::class_<AttributeFiltersPybind>(m, "AttributeFilters", py::module_local(false))
    .def(py::init<MorphologicalTreePybindPtr>())
    .def("filteringMin", py::overload_cast<py::array_t<float> &, float>(&AttributeFiltersPybind::filteringByPruningMin))
    .def("filteringMin", py::overload_cast<std::vector<bool>&>(&AttributeFiltersPybind::filteringByPruningMin))
    .def("filteringMax", py::overload_cast<std::vector<bool>&>(&AttributeFiltersPybind::filteringByPruningMax))
    .def("filteringDirectRule", py::overload_cast<std::vector<bool>&>(&AttributeFiltersPybind::filteringByDirectRule))
    .def("filteringSubtractiveRule", py::overload_cast<std::vector<bool>&>(&AttributeFiltersPybind::filteringBySubtractiveRule))
    .def("filteringSubtractiveScoreRule", py::overload_cast<std::vector<float>&>(&AttributeFiltersPybind::filteringBySubtractiveScoreRule))
    .def("filteringMax", py::overload_cast<py::array_t<float> &, float>(&AttributeFiltersPybind::filteringByPruningMax))
    .def("filteringByExtinction", py::overload_cast<py::array_t<float> &, int>(&AttributeFiltersPybind::filteringByExtinctionValue))
    .def("saliencyMapByExtinction", py::overload_cast<py::array_t<float> &, int, bool>(&AttributeFiltersPybind::saliencyMapByExtinctionValue), "attr"_a, "leafToKeep"_a, "unweighted"_a = false)
    .def("getAdaptativeCriterion", &AttributeFiltersPybind::getAdaptativeCriterion);       
}

void init_ExtinctionValues(py::module &m){
    py::class_<ExtinctionValuesPybind>(m, "ExtinctionValues", py::module_local(false))
    .def(py::init<MorphologicalTreePybindPtr, py::array_t<float>&>())
    .def("filtering", &ExtinctionValuesPybind::filtering)
    .def("saliencyMap", &ExtinctionValuesPybind::saliencyMap, "leafToKeep"_a, "unweighted"_a = true)
    .def("getExtinctionValues", &ExtinctionValuesPybind::getExtinctionValuesPy);
    
}


void init_AdjacencyRelation(py::module &m){
    	py::class_<AdjacencyRelation>(m, "AdjacencyRelation", py::module_local(false))
        .def(py::init<int, int, double>())
        .def_property_readonly("size", &AdjacencyRelation::getSize )
        .def("getAdjPixels", py::overload_cast<int, int>( &AdjacencyRelation::getAdjPixels ));
}


void init_UltimateAttributeOpening(py::module &m){
    	py::class_<UltimateAttributeOpeningPybind>(m, "UltimateAttributeOpening", py::module_local(false))
        .def(py::init<MorphologicalTreePybindPtr, py::array_t<float>&>())
        .def("execute", py::overload_cast<int>(&UltimateAttributeOpeningPybind::execute))
        .def("executeWithMSER", &UltimateAttributeOpeningPybind::executeWithMSER)
        .def("getMaxConstrastImage", &UltimateAttributeOpeningPybind::getMaxConstrastImage)
        .def("getAssociatedImage", &UltimateAttributeOpeningPybind::getAssociatedImage)
        .def("getAssociatedColoredImage", &UltimateAttributeOpeningPybind::getAssociatedColorImage);
}

void init_ResidualTree(py::module &m){
    	py::class_<ResidualTreePybind>(m, "ResidualTree", py::module_local(false))
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
    	py::class_<AttributeOpeningPrimitivesFamilyPybind, std::shared_ptr<AttributeOpeningPrimitivesFamilyPybind>>(m, "AttributeOpeningPrimitivesFamily", py::module_local(false))
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
    init_MorphologicalTree(m);
    init_AttributeComputedIncrementally(m);
    init_ContoursComputedIncrementally(m);
    init_AttributeFilters(m);
    init_ExtinctionValues(m);
    init_AdjacencyRelation(m);

    init_UltimateAttributeOpening(m);
    init_ResidualTree(m);
    init_AttributeOpeningPrimitivesFamily(m);

}
