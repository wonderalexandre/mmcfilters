#include "../mmcfilters/utils/AdjacencyRelation.hpp"
#include "../mmcfilters/utils/Common.hpp"
#include "../mmcfilters/contours/ContoursComputedIncrementally.hpp"
#include "../mmcfilters/trees/WeightedMorphologicalTree.hpp"
#include "../mmcfilters/trees/adjust/CasfComponentTrees.hpp"
#include "../mmcfilters/trees/adjust/DualMinMaxTreeIncrementalFilter.hpp"

#include "AttributeComputedIncrementallyPybind.hpp"
#include "ContoursComputedIncrementallyPybind.hpp"
#include "MorphologicalTreePybind.hpp"
#include "ExtinctionValuesPybind.hpp"
#include "AttributeFiltersPybind.hpp"
#include "UltimateAttributeOpeningPybind.hpp"
#include "AttributeOpeningPrimitivesFamilyPybind.hpp"


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>

using namespace mmcfilters;

#include <optional>
#include <sstream>
#include <type_traits>
#include <utility>
#include <vector>


namespace py = pybind11;
using namespace pybind11::literals;

using UInt8InputArray = py::array_t<uint8_t, py::array::c_style | py::array::forcecast>;

template <class Range>
std::vector<NodeId> collectNodeIds(const Range &range) {
    std::vector<NodeId> ids;
    for (NodeId id : range) {
        ids.push_back(id);
    }
    return ids;
}

namespace {

ImageUInt8Ptr imageFromArray(const UInt8InputArray& input) {
    auto buf = input.request();
    if (buf.ndim != 2) {
        throw std::invalid_argument("input must be a 2D uint8 array");
    }
    const int rows = static_cast<int>(buf.shape[0]);
    const int cols = static_cast<int>(buf.shape[1]);
    return ImageUInt8::fromExternal(static_cast<uint8_t*>(buf.ptr), rows, cols);
}

MorphologicalTree& topology(MorphologicalTreePybind& tree) {
    return tree;
}

const MorphologicalTree& topology(const MorphologicalTreePybind& tree) {
    return tree;
}

const MorphologicalTree& topology(const WeightedMorphologicalTree& weighted) {
    return weighted.topology();
}

AltitudeType altitudeOf(const WeightedMorphologicalTree& weighted, NodeId nodeId) {
    return weighted.getAltitude(nodeId);
}

AltitudeDiffType residueOf(const WeightedMorphologicalTree& weighted, NodeId nodeId) {
    return weighted.getNodeResidue(nodeId);
}

py::array_t<uint8_t> reconstructionImageOf(const WeightedMorphologicalTree& weighted) {
    return PybindUtils::toNumpy(weighted.reconstructionImage());
}

std::pair<std::vector<NodeId>, std::vector<AltitudeType>> exportHigraHierarchyOf(const WeightedMorphologicalTree& weighted) {
    return weighted.exportHigraHierarchy();
}

class DualMinMaxTreeIncrementalFilterPybind {
private:
    std::shared_ptr<WeightedMorphologicalTree> minTree_;
    std::shared_ptr<WeightedMorphologicalTree> maxTree_;
    AdjacencyRelation adjacency_;
    mmcfilters::adjust::DualMinMaxTreeIncrementalFilter<AltitudeType> adjust_;

    static AdjacencyRelation makeAdjacency(const WeightedMorphologicalTree& minTree, const WeightedMorphologicalTree& maxTree) {
        const AdjacencyRelation* minAdjacency = minTree.topology().getAdjacencyRelation();
        const AdjacencyRelation* maxAdjacency = maxTree.topology().getAdjacencyRelation();
        if (minAdjacency == nullptr || maxAdjacency == nullptr) {
            throw std::invalid_argument("DualMinMaxTreeIncrementalFilter requires weighted component trees with adjacency information.");
        }
        if (minTree.topology().getNumRowsOfImage() != maxTree.topology().getNumRowsOfImage() ||
            minTree.topology().getNumColsOfImage() != maxTree.topology().getNumColsOfImage() ||
            minTree.topology().getNumTotalProperParts() != maxTree.topology().getNumTotalProperParts()) {
            throw std::invalid_argument("DualMinMaxTreeIncrementalFilter requires trees built on the same image domain.");
        }
        return *minAdjacency;
    }

public:
    DualMinMaxTreeIncrementalFilterPybind(std::shared_ptr<WeightedMorphologicalTree> minTree, std::shared_ptr<WeightedMorphologicalTree> maxTree)
        : minTree_(std::move(minTree)),
          maxTree_(std::move(maxTree)),
          adjacency_(makeAdjacency(*minTree_, *maxTree_)),
          adjust_(minTree_.get(), maxTree_.get(), adjacency_) {}

    void pruneMaxTreeAndUpdateMinTree(const std::vector<NodeId>& nodesToPrune) {
        adjust_.pruneMaxTreeAndUpdateMinTree(nodesToPrune);
    }

    void pruneMinTreeAndUpdateMaxTree(const std::vector<NodeId>& nodesToPrune) {
        adjust_.pruneMinTreeAndUpdateMaxTree(nodesToPrune);
    }

    std::shared_ptr<WeightedMorphologicalTree> minTree() const {
        return minTree_;
    }

    std::shared_ptr<WeightedMorphologicalTree> maxTree() const {
        return maxTree_;
    }
};

template <class TreeLike, class PyClass>
void bindTreeQueryApi(PyClass& cls) {
    cls.def_property_readonly("numInternalNodeSlots", [](TreeLike &self) {
            return topology(self).getNumInternalNodeSlots();
        })
        .def_property_readonly("numTotalProperParts", [](TreeLike &self) {
            return topology(self).getNumTotalProperParts();
        })
        .def_property_readonly("numHigraNodes", [](TreeLike &self) {
            return topology(self).getNumHigraNodes();
        }, "Size of the preserved imported Higra node-id domain.")
        .def("getRoot", [](TreeLike &self) {
            return topology(self).getRoot();
        }, "Return the current root node id.")
        .def_property_readonly("root", [](TreeLike &self) {
            return topology(self).getRoot();
        }, "Current root node id.")
        .def_property_readonly("numFreeNodeSlots", [](TreeLike &self) {
            return topology(self).getNumFreeNodeSlots();
        })
        .def_property_readonly("numLeafNodes", [](TreeLike &self) {
            return topology(self).getNumLeafNodes();
        })
        .def("getAliveNodeIds", [](TreeLike &self) {
            return collectNodeIds(topology(self).getAliveNodeIds());
        }, "Return all alive internal-node ids in the dense node-id domain.")
        .def_property_readonly("aliveNodeIds", [](TreeLike &self) {
            return collectNodeIds(topology(self).getAliveNodeIds());
        }, "Alive internal-node ids in the dense node-id domain.")
        .def_property_readonly("alive_node_ids", [](TreeLike &self) {
            return collectNodeIds(topology(self).getAliveNodeIds());
        }, "Alive internal-node ids in the dense node-id domain.")
        .def("getLeafNodeIds", [](TreeLike &self) {
            return topology(self).getLeaves();
        }, "Return alive leaf node ids in the dense node-id domain.")
        .def_property_readonly("leafNodeIds", [](TreeLike &self) {
            return topology(self).getLeaves();
        }, "Alive leaf node ids in the dense node-id domain.")
        .def_property_readonly("leaf_node_ids", [](TreeLike &self) {
            return topology(self).getLeaves();
        }, "Alive leaf node ids in the dense node-id domain.")
        .def("getChildren", [](TreeLike &self, NodeId nodeId) {
            return collectNodeIds(topology(self).getChildren(nodeId));
        }, "nodeId"_a, "Return the direct children of a node in the dense node-id domain.")
        .def("childrenOf", [](TreeLike &self, NodeId nodeId) {
            return collectNodeIds(topology(self).getChildren(nodeId));
        }, "nodeId"_a, "Return the direct children of a node in the dense node-id domain.")
        .def("children_of", [](TreeLike &self, NodeId nodeId) {
            return collectNodeIds(topology(self).getChildren(nodeId));
        }, "nodeId"_a, "Return the direct children of a node in the dense node-id domain.")
        .def("getNodeNumDescendants", [](TreeLike &self, NodeId nodeId) {
            return topology(self).getNodeNumDescendants(nodeId);
        }, "nodeId"_a, "Return the number of descendants of nodeId.")
        .def("getNodeNumSiblings", [](TreeLike &self, NodeId nodeId) {
            return topology(self).getNodeNumSiblings(nodeId);
        }, "nodeId"_a, "Return the number of siblings of nodeId.")
        .def("getNumProperParts", [](TreeLike &self, NodeId nodeId) {
            return topology(self).getNumProperParts(nodeId);
        }, "nodeId"_a, "Return the number of direct proper parts owned by nodeId.")
        .def("getNodeTimePreOrder", [](TreeLike &self, NodeId nodeId) {
            return topology(self).getNodeTimePreOrder(nodeId);
        }, "nodeId"_a, "Return the preorder timestamp of nodeId.")
        .def("getNodeTimePostOrder", [](TreeLike &self, NodeId nodeId) {
            return topology(self).getNodeTimePostOrder(nodeId);
        }, "nodeId"_a, "Return the postorder timestamp of nodeId.")
        .def("getProperParts", [](TreeLike &self, NodeId nodeId) {
            return collectNodeIds(topology(self).getProperParts(nodeId));
        }, "nodeId"_a, "Return the proper parts owned directly by a node.")
        .def("properPartsOf", [](TreeLike &self, NodeId nodeId) {
            return collectNodeIds(topology(self).getProperParts(nodeId));
        }, "nodeId"_a, "Return the proper parts owned directly by a node.")
        .def("proper_parts_of", [](TreeLike &self, NodeId nodeId) {
            return collectNodeIds(topology(self).getProperParts(nodeId));
        }, "nodeId"_a, "Return the proper parts owned directly by a node.")
        .def("reconstructNode", [](TreeLike &self, NodeId nodeId) {
            return MorphologicalTreePybind::reconstructNode(topology(self), nodeId);
        }, "nodeId"_a, "Reconstruct a binary mask for the connected component represented by nodeId.")
        .def("getPostOrderNodes", [](TreeLike &self, std::optional<NodeId> rootNodeId) {
            return rootNodeId.has_value()
                ? collectNodeIds(topology(self).getPostOrderNodes(*rootNodeId))
                : collectNodeIds(topology(self).getPostOrderNodes());
        }, "rootNodeId"_a = std::nullopt)
        .def("getIteratorBreadthFirstTraversal", [](TreeLike &self, std::optional<NodeId> rootNodeId) {
            return rootNodeId.has_value()
                ? collectNodeIds(topology(self).getIteratorBreadthFirstTraversal(*rootNodeId))
                : collectNodeIds(topology(self).getIteratorBreadthFirstTraversal());
        }, "rootNodeId"_a = std::nullopt)
        .def("getPathToRootNodes", [](TreeLike &self, NodeId nodeId) {
            return collectNodeIds(topology(self).getPathToRootNodes(nodeId));
        }, "nodeId"_a)
        .def("getPathBetweenNodes", [](TreeLike &self, NodeId sourceNodeId, NodeId targetNodeId) {
            return collectNodeIds(topology(self).getPathBetweenNodes(sourceNodeId, targetNodeId));
        }, "sourceNodeId"_a, "targetNodeId"_a)
        .def("getNodeSubtree", [](TreeLike &self, NodeId nodeId) {
            return collectNodeIds(topology(self).getNodeSubtree(nodeId));
        }, "nodeId"_a)
        .def("nodeSubtreeOf", [](TreeLike &self, NodeId nodeId) {
            return collectNodeIds(topology(self).getNodeSubtree(nodeId));
        }, "nodeId"_a)
        .def("node_subtree_of", [](TreeLike &self, NodeId nodeId) {
            return collectNodeIds(topology(self).getNodeSubtree(nodeId));
        }, "nodeId"_a)
        .def("getDescendants", [](TreeLike &self, NodeId nodeId) {
            return collectNodeIds(topology(self).getDescendants(nodeId));
        }, "nodeId"_a)
        .def("descendantsOf", [](TreeLike &self, NodeId nodeId) {
            return collectNodeIds(topology(self).getDescendants(nodeId));
        }, "nodeId"_a)
        .def("descendants_of", [](TreeLike &self, NodeId nodeId) {
            return collectNodeIds(topology(self).getDescendants(nodeId));
        }, "nodeId"_a)
        .def("getNodeParent", [](TreeLike &self, NodeId nodeId) {
            return topology(self).getNodeParent(nodeId);
        }, "nodeId"_a)
        .def("parentOf", [](TreeLike &self, NodeId nodeId) {
            return topology(self).getNodeParent(nodeId);
        }, "nodeId"_a)
        .def("parent_of", [](TreeLike &self, NodeId nodeId) {
            return topology(self).getNodeParent(nodeId);
        }, "nodeId"_a)
        .def("getSmallestComponent", [](TreeLike &self, int pixelId) {
            return topology(self).getSmallestComponent(pixelId);
        }, "pixelId"_a)
        .def("smallestComponentOf", [](TreeLike &self, int pixelId) {
            return topology(self).getSmallestComponent(pixelId);
        }, "pixelId"_a)
        .def("smallest_component_of", [](TreeLike &self, int pixelId) {
            return topology(self).getSmallestComponent(pixelId);
        }, "pixelId"_a)
        .def("getHigraNodeId", [](TreeLike &self, NodeId nodeId) {
            return topology(self).getHigraNodeId(nodeId);
        }, "nodeId"_a, "Return the preserved imported Higra node id for a live internal NodeId, or InvalidNode.")
        .def("getNumChildren", [](TreeLike &self, NodeId nodeId) {
            return topology(self).getNumChildren(nodeId);
        }, "nodeId"_a)
        .def("getFirstChild", [](TreeLike &self, NodeId nodeId) {
            return topology(self).getFirstChild(nodeId);
        }, "nodeId"_a)
        .def("getNextSibling", [](TreeLike &self, NodeId nodeId) {
            return topology(self).getNextSibling(nodeId);
        }, "nodeId"_a)
        .def("isNode", [](TreeLike &self, NodeId nodeId) {
            return topology(self).isNode(nodeId);
        }, "nodeId"_a)
        .def("isProperPart", [](TreeLike &self, NodeId nodeId) {
            return topology(self).isProperPart(nodeId);
        }, "nodeId"_a)
        .def("isAlive", [](TreeLike &self, NodeId nodeId) {
            return topology(self).isAlive(nodeId);
        }, "nodeId"_a)
        .def("isRoot", [](TreeLike &self, NodeId nodeId) {
            return topology(self).isRoot(nodeId);
        }, "nodeId"_a)
        .def("isLeaf", [](TreeLike &self, NodeId nodeId) {
            return topology(self).isLeaf(nodeId);
        }, "nodeId"_a)
        .def("hasChild", [](TreeLike &self, NodeId parentId, NodeId childId) {
            return topology(self).hasChild(parentId, childId);
        }, "parentId"_a, "childId"_a)
        .def("pruneNode", [](TreeLike &self, NodeId nodeId) {
            self.pruneNode(nodeId);
        }, "nodeId"_a)
        .def("mergeNodeIntoParent", [](TreeLike &self, NodeId nodeId) {
            self.mergeNodeIntoParent(nodeId);
        }, "nodeId"_a)
        .def_property_readonly("treeType", [](TreeLike& self) { return topology(self).getTreeType(); })
        .def_property_readonly("hasAdjacencyRelation", [](TreeLike& self) { return topology(self).hasAdjacencyRelation(); })
        .def_property_readonly("hasTreeOfShapesAdjacencyPolicy", [](TreeLike& self) { return topology(self).hasTreeOfShapesAdjacencyPolicy(); })
        .def("getTreeOfShapesMinTreeAdjacencyRadius", [](TreeLike& self) {
            return topology(self).getTreeOfShapesMinTreeAdjacencyRadius();
        })
        .def("getTreeOfShapesMaxTreeAdjacencyRadius", [](TreeLike& self) {
            return topology(self).getTreeOfShapesMaxTreeAdjacencyRadius();
        })
        .def("getTreeOfShapesMinTreeAdjacencyRelation", [](TreeLike& self) -> const AdjacencyRelation& {
            const AdjacencyRelation* adjacency = topology(self).getTreeOfShapesMinTreeAdjacencyRelation();
            if (!adjacency) {
                throw std::runtime_error("Tree-of-shapes adjacency policy is not available.");
            }
            return *adjacency;
        }, py::return_value_policy::reference_internal)
        .def("getTreeOfShapesMaxTreeAdjacencyRelation", [](TreeLike& self) -> const AdjacencyRelation& {
            const AdjacencyRelation* adjacency = topology(self).getTreeOfShapesMaxTreeAdjacencyRelation();
            if (!adjacency) {
                throw std::runtime_error("Tree-of-shapes adjacency policy is not available.");
            }
            return *adjacency;
        }, py::return_value_policy::reference_internal)
        .def_property_readonly("numRows", [](TreeLike& self) { return topology(self).getNumRowsOfImage(); })
        .def_property_readonly("numCols", [](TreeLike& self) { return topology(self).getNumColsOfImage(); })
        .def_property_readonly("numNodes", [](TreeLike& self) { return topology(self).getNumNodes(); });

    if constexpr (std::is_same_v<TreeLike, WeightedMorphologicalTree>) {
        cls.def("getAltitude", [](TreeLike &self, NodeId nodeId) {
                return altitudeOf(self, nodeId);
            }, "nodeId"_a, "Return the altitude associated with nodeId.")
            .def("altitudeOf", [](TreeLike &self, NodeId nodeId) {
                return altitudeOf(self, nodeId);
            }, "nodeId"_a, "Return the altitude associated with nodeId.")
            .def("altitude_of", [](TreeLike &self, NodeId nodeId) {
                return altitudeOf(self, nodeId);
            }, "nodeId"_a, "Return the altitude associated with nodeId.")
            .def("getNodeResidue", [](TreeLike &self, NodeId nodeId) {
                return residueOf(self, nodeId);
            }, "nodeId"_a, "Return the residue between nodeId and its parent.")
            .def("residueOf", [](TreeLike &self, NodeId nodeId) {
                return residueOf(self, nodeId);
            }, "nodeId"_a, "Return the residue between nodeId and its parent.")
            .def("residue_of", [](TreeLike &self, NodeId nodeId) {
                return residueOf(self, nodeId);
            }, "nodeId"_a, "Return the residue between nodeId and its parent.")
            .def("getRepresentativeProperPartsByFlood", [](TreeLike &self, NodeId nodeId) {
                return MorphologicalTreePybind::representativeProperPartsByFlood(
                    self.topology(),
                    self.getAltitudeBuffer(),
                    nodeId);
            }, "nodeId"_a, "Return one derived representative proper part per flat zone inside the connected component of nodeId.")
            .def("reconstructionImage", [](TreeLike &self) {
                return reconstructionImageOf(self);
            }, "Reconstruct the current tree into a 2D image using the attached image domain.")
            .def("exportHigraHierarchy", [](TreeLike &self) {
                return exportHigraHierarchyOf(self);
            }, "Export the current rooted tree to a new compact Higra (parent, altitude) representation.");
    }
}

} // namespace

void init_MorphologicalTree(py::module &m){
      py::enum_<ToSInterpolation>(m, "ToSInterpolation", py::module_local(false))
        .value("SelfDual", ToSInterpolation::SelfDual)
        .value("Min4cMax8c", ToSInterpolation::Min4cMax8c)
        .value("Min8cMax4c", ToSInterpolation::Min8cMax4c)
        .export_values();

      py::enum_<NodeIdSpace>(m, "NodeIdSpace", py::module_local(false))
        .value("MORPHOLOGICAL_TREE", NodeIdSpace::MORPHOLOGICAL_TREE)
        .value("HIGRA", NodeIdSpace::HIGRA)
        .export_values();

      py::class_<MorphologicalTree, std::shared_ptr<MorphologicalTree>>(m, "MorphologicalTreeBase", py::module_local(false));
      auto treeCls = py::class_<MorphologicalTreePybind, std::shared_ptr<MorphologicalTreePybind>>(m, "MorphologicalTree", py::module_local(false),
        "Morphological tree with a NodeId-first public API. Prefer getRoot/getAliveNodeIds/getChildren/getProperParts and related NodeId-based operations for new code. "
        "Weighted quantities such as altitude, image reconstruction, and compact Higra parent/altitude export live on WeightedMorphologicalTree.");
      treeCls
        .def_static("createComponentTree", [](UInt8InputArray input, bool isMaxtree, double radius) {
            return std::make_shared<MorphologicalTreePybind>(input, isMaxtree, radius);
        }, "input"_a, "isMaxtree"_a, "radius"_a = 1.5)
        .def_static("createMaxTree", [](UInt8InputArray input, double radius) {
            return std::make_shared<MorphologicalTreePybind>(input, true, radius);
        }, "input"_a, "radius"_a = 1.5)
        .def_static("createMinTree", [](UInt8InputArray input, double radius) {
            return std::make_shared<MorphologicalTreePybind>(input, false, radius);
        }, "input"_a, "radius"_a = 1.5)
        .def_static("createTreeOfShapes", [](UInt8InputArray input, ToSInterpolation interpolation, int infinitySeedRow, int infinitySeedCol) {
            return std::make_shared<MorphologicalTreePybind>(input, interpolation, infinitySeedRow, infinitySeedCol);
        },
            "input"_a,
            "interpolation"_a = ToSInterpolation::SelfDual,
            "infinitySeedRow"_a = ToSDefaultInfinityRow,
            "infinitySeedCol"_a = ToSDefaultInfinityCol)
        .def_static("createFromHigraParent", [](const std::vector<NodeId>& parent, int rows, int cols, int treeType, std::optional<double> radius) {
            return std::make_shared<MorphologicalTreePybind>(
                MorphologicalTree::createFromHigraParent(
                    parent,
                    rows,
                    cols,
                    treeType,
                    radius ? std::optional<AdjacencyRelation>(std::in_place, rows, cols, *radius) : std::nullopt));
        },
            "parent"_a,
            "rows"_a,
            "cols"_a,
            "treeType"_a,
            "radius"_a = py::none(),
            "Create a tree from an imported static Higra parent representation [leaves | internal nodes]. "
            "The imported Higra node-id domain is preserved until the tree is edited.");
      bindTreeQueryApi<MorphologicalTreePybind>(treeCls);

      auto weightedCls = py::class_<WeightedMorphologicalTree, std::shared_ptr<WeightedMorphologicalTree>>(m, "WeightedMorphologicalTree", py::module_local(false),
        "Wrapper pairing MorphologicalTree topology with an external dense altitude buffer. "
        "Imports can preserve an original Higra node-id domain; exports always create a new compact Higra domain.");
      weightedCls
        .def_static("createComponentTree", [](UInt8InputArray input, bool isMaxtree, double radius) {
            return std::make_shared<WeightedMorphologicalTree>(
                WeightedMorphologicalTree::createComponentTree(imageFromArray(input), isMaxtree, radius));
        }, "input"_a, "isMaxtree"_a, "radius"_a = 1.5)
        .def_static("createMaxTree", [](UInt8InputArray input, double radius) {
            return std::make_shared<WeightedMorphologicalTree>(
                WeightedMorphologicalTree::createComponentTree(imageFromArray(input), true, radius));
        }, "input"_a, "radius"_a = 1.5)
        .def_static("createMinTree", [](UInt8InputArray input, double radius) {
            return std::make_shared<WeightedMorphologicalTree>(
                WeightedMorphologicalTree::createComponentTree(imageFromArray(input), false, radius));
        }, "input"_a, "radius"_a = 1.5)
        .def_static("createTreeOfShapes", [](UInt8InputArray input, ToSInterpolation interpolation, int infinitySeedRow, int infinitySeedCol) {
            return std::make_shared<WeightedMorphologicalTree>(
                WeightedMorphologicalTree::createTreeOfShapes(imageFromArray(input), interpolation, infinitySeedRow, infinitySeedCol));
        },
            "input"_a,
            "interpolation"_a = ToSInterpolation::SelfDual,
            "infinitySeedRow"_a = ToSDefaultInfinityRow,
            "infinitySeedCol"_a = ToSDefaultInfinityCol)
        .def("setAltitude", [](WeightedMorphologicalTree &tree, NodeId nodeId, AltitudeType altitude) {
            if (!tree.topology().isNode(nodeId)) {
                throw std::invalid_argument("invalid NodeId for altitude update");
            }
            tree.setAltitude(nodeId, altitude);
        }, "nodeId"_a, "altitude"_a, "Set one node altitude inside the external dense altitude buffer.")
        .def("setAltitudeBuffer", [](WeightedMorphologicalTree &tree, const std::vector<AltitudeType>& altitude) {
            PybindUtils::requireVectorSize(altitude, static_cast<std::size_t>(tree.topology().getNumInternalNodeSlots()), "altitude");
            tree.setAltitudeBuffer(altitude);
        }, "altitude"_a, "Replace the dense altitude buffer indexed by internal NodeId.")
        .def_property(
            "altitude",
            [](const WeightedMorphologicalTree &tree) {
                return AltitudeBuffer(tree.getAltitudeBuffer());
            },
            [](WeightedMorphologicalTree &tree, const std::vector<AltitudeType>& altitude) {
                PybindUtils::requireVectorSize(altitude, static_cast<std::size_t>(tree.topology().getNumInternalNodeSlots()), "altitude");
                tree.setAltitudeBuffer(altitude);
            },
            "Dense altitude buffer indexed by internal NodeId.")
        .def("validateAltitudeBufferShape", static_cast<void (WeightedMorphologicalTree::*)() const>(&WeightedMorphologicalTree::validateAltitudeBufferShape))
        .def("validateMonotoneAltitude", static_cast<void (WeightedMorphologicalTree::*)() const>(&WeightedMorphologicalTree::validateMonotoneAltitude))
        .def_static("createFromHigraParent", [](const std::vector<NodeId>& parent, const std::vector<AltitudeType>& altitude, int rows, int cols, int treeType, std::optional<double> radius) {
            if (parent.size() != altitude.size()) {
                throw std::invalid_argument("parent and altitude must have the same size");
            }

            return std::make_shared<WeightedMorphologicalTree>(
                WeightedMorphologicalTree::createFromHigraParent(
                    parent,
                    altitude,
                    rows,
                    cols,
                    treeType,
                    radius ? std::optional<AdjacencyRelation>(std::in_place, rows, cols, *radius) : std::nullopt));
        },
            "parent"_a,
            "altitude"_a,
            "rows"_a,
            "cols"_a,
            "treeType"_a,
            "radius"_a = py::none(),
            "Create a weighted tree from an imported static Higra parent/altitude representation [leaves | internal nodes]. "
            "The imported Higra node-id domain is preserved until the tree is edited.");
      bindTreeQueryApi<WeightedMorphologicalTree>(weightedCls);
	}

void init_ContoursComputedIncrementally(py::module &m){
    // Alias locais
    using Contours     = ContoursComputedIncrementally::IncrementalContours;
    using ContourRange = Contours::ContourRange;
    using Range = decltype(std::declval<Contours&>().contoursByNode());
    using Iter  = decltype(std::declval<Range&>().begin());
    
    // Torna ContourRange iterável em Python
    py::class_<ContourRange>(m, "ContourRange", py::module_local(false))
        .def("__iter__", [](const ContourRange& p) {
            // usa os iteradores já existentes do range
            return py::make_iterator(p.begin(), p.end());
        }, py::keep_alive<0, 1>())
        .def("empty", &ContourRange::empty);

    struct ContoursIterator {
        Range range;
        Iter it, itEnd;
        ContoursIterator(Contours& self)
            : range(self.contoursByNode()), it(range.begin()), itEnd(range.end()) {}
    };

    py::class_<ContoursIterator>(m, "ContoursIterator", py::module_local(false))
        .def(py::init<Contours&>())
        .def("__iter__", [](ContoursIterator& self) -> ContoursIterator& { return self; }, py::return_value_policy::reference_internal)
        .def("__next__", [](ContoursIterator& self) -> py::object {
            if (self.it == self.itEnd) throw py::stop_iteration();
            auto entry = *self.it++;  // <— cópia do par temporário
            auto nodeId = std::get<0>(entry);
            auto proxy  = std::get<1>(entry);  // ContourRange
            return py::make_tuple(nodeId, proxy);
        });


    py::class_<Contours, std::shared_ptr<Contours>>(m, "Contours", py::module_local(false))
        .def("contoursByNode", [](Contours &self) {
            return ContoursIterator(self);
        }, py::keep_alive<0, 1>())
        .def("getContour", &Contours::getContour)
        .def("materializeAll", &Contours::materializeAll)
        .def_property_readonly("isMaterialized", &Contours::isMaterialized)
        .def("isContourMaterialized", &Contours::isContourMaterialized);

    py::class_<ContoursComputedIncrementallyPybind>(m, "ContourComputation", py::module_local(false))
        .def_static("extraction", py::overload_cast<MorphologicalTreePybindPtr>(&ContoursComputedIncrementallyPybind::extraction))
        .def_static("extraction", py::overload_cast<std::shared_ptr<WeightedMorphologicalTree>>(&ContoursComputedIncrementallyPybind::extraction));
}

void init_AttributeComputedIncrementally(py::module &m){
        auto cls = py::class_<AttributeComputedIncrementallyPybind>(m, "Attribute", py::module_local(false),
        "Incremental attribute utilities based on dense node-id traversal.")
        .def_static(
            "traversePostOrder",
            [](MorphologicalTreePybind &tree,
               std::function<void(NodeId)> preProcessing,
               std::function<void(NodeId, NodeId)> mergeProcessing,
               std::function<void(NodeId)> postProcessing,
               std::optional<NodeId> rootNodeIdOpt) {
                const NodeId rootNodeId = rootNodeIdOpt.value_or(tree.getRoot());
                if (!tree.isNode(rootNodeId) || !tree.isAlive(rootNodeId)) {
                    throw std::invalid_argument("rootNodeId inválido");
                }

                AttributeComputedIncrementally::traversePostOrder(
                    tree,
                    rootNodeId,
                    std::move(preProcessing),
                    std::move(mergeProcessing),
                    std::move(postProcessing)
                );
            },
            py::arg("tree"),
            py::arg("preProcessing"),
            py::arg("mergeProcessing"),
            py::arg("postProcessing"),
            py::arg("rootNodeId") = std::optional<NodeId>{},
            "Callback traversal using dense node-id values."
        )
        .def_static(
            "traversePostOrder",
            [](WeightedMorphologicalTree &weighted,
               std::function<void(NodeId)> preProcessing,
               std::function<void(NodeId, NodeId)> mergeProcessing,
               std::function<void(NodeId)> postProcessing,
               std::optional<NodeId> rootNodeIdOpt) {
                const auto& tree = weighted.topology();
                const NodeId rootNodeId = rootNodeIdOpt.value_or(tree.getRoot());
                if (!tree.isNode(rootNodeId) || !tree.isAlive(rootNodeId)) {
                    throw std::invalid_argument("rootNodeId inválido");
                }

                AttributeComputedIncrementally::traversePostOrder(
                    tree,
                    rootNodeId,
                    std::move(preProcessing),
                    std::move(mergeProcessing),
                    std::move(postProcessing)
                );
            },
            py::arg("tree"),
            py::arg("preProcessing"),
            py::arg("mergeProcessing"),
            py::arg("postProcessing"),
            py::arg("rootNodeId") = std::optional<NodeId>{},
            "Callback traversal using dense node-id values."
        )
	        .def_static("computeAttributes", py::overload_cast<MorphologicalTreePybindPtr, const std::vector<AttributeOrGroup>&, NodeIdSpace>(&AttributeComputedIncrementallyPybind::computeAttributesFromList),
                py::arg("tree"),
                py::arg("attributes"),
                py::arg("outputSpace") = NodeIdSpace::MORPHOLOGICAL_TREE)
            .def_static("computeAttributes", py::overload_cast<std::shared_ptr<WeightedMorphologicalTree>, const std::vector<AttributeOrGroup>&, NodeIdSpace>(&AttributeComputedIncrementallyPybind::computeAttributesFromList),
                py::arg("tree"),
                py::arg("attributes"),
                py::arg("outputSpace") = NodeIdSpace::MORPHOLOGICAL_TREE)
	        .def_static("computeSingleAttribute", py::overload_cast<MorphologicalTreePybindPtr, Attribute, NodeIdSpace>(&AttributeComputedIncrementallyPybind::computeSingleAttribute),
                py::arg("tree"),
                py::arg("attribute"),
                py::arg("outputSpace") = NodeIdSpace::MORPHOLOGICAL_TREE)
            .def_static("computeSingleAttribute", py::overload_cast<std::shared_ptr<WeightedMorphologicalTree>, Attribute, NodeIdSpace>(&AttributeComputedIncrementallyPybind::computeSingleAttribute),
                py::arg("tree"),
                py::arg("attribute"),
                py::arg("outputSpace") = NodeIdSpace::MORPHOLOGICAL_TREE)
	        .def_static("computeSingleAttributeWithDelta", py::overload_cast<MorphologicalTreePybindPtr, Attribute, int, std::string, NodeIdSpace>(&AttributeComputedIncrementallyPybind::computeSingleAttributeWithDelta),
                py::arg("tree"),
                py::arg("attribute"),
                py::arg("delta"),
                py::arg("padding") = "last-padding",
                py::arg("outputSpace") = NodeIdSpace::MORPHOLOGICAL_TREE)
            .def_static("computeSingleAttributeWithDelta", py::overload_cast<std::shared_ptr<WeightedMorphologicalTree>, Attribute, int, std::string, NodeIdSpace>(&AttributeComputedIncrementallyPybind::computeSingleAttributeWithDelta),
                py::arg("tree"),
                py::arg("attribute"),
                py::arg("delta"),
                py::arg("padding") = "last-padding",
                py::arg("outputSpace") = NodeIdSpace::MORPHOLOGICAL_TREE)
        .def_static("describe", &AttributeComputedIncrementallyPybind::describeAttribute)
        .def_static("computeAttributeMapping", py::overload_cast<MorphologicalTreePybindPtr, Attribute>(&AttributeComputedIncrementallyPybind::computeAttributeMapping))
        .def_static("computeAttributeMapping", py::overload_cast<std::shared_ptr<WeightedMorphologicalTree>, Attribute>(&AttributeComputedIncrementallyPybind::computeAttributeMapping));

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
            .value("MAX_DIST", Attribute::MAX_DIST)
            
            .export_values();
}

void init_AttributeFilters(py::module &m){
    using FloatArray = py::array_t<float, py::array::c_style | py::array::forcecast>;

    py::class_<AttributeFiltersPybind>(m, "AttributeFilters", py::module_local(false))
    .def(py::init<MorphologicalTreePybindPtr>())
    .def(py::init<std::shared_ptr<WeightedMorphologicalTree>>())
    .def("filteringMin", [](AttributeFiltersPybind &self, FloatArray attr, float threshold) {
        return self.filteringByPruningMin(std::move(attr), threshold);
    })
    .def("filteringMin", py::overload_cast<std::vector<bool>&>(&AttributeFiltersPybind::filteringByPruningMin))
    .def("filteringMax", py::overload_cast<std::vector<bool>&>(&AttributeFiltersPybind::filteringByPruningMax))
    .def("filteringDirectRule", py::overload_cast<std::vector<bool>&>(&AttributeFiltersPybind::filteringByDirectRule))
    .def("filteringSubtractiveRule", py::overload_cast<std::vector<bool>&>(&AttributeFiltersPybind::filteringBySubtractiveRule))
    .def("filteringSubtractiveScoreRule", py::overload_cast<std::vector<float>&>(&AttributeFiltersPybind::filteringBySubtractiveScoreRule))
    .def("filteringMax", [](AttributeFiltersPybind &self, FloatArray attr, float threshold) {
        return self.filteringByPruningMax(std::move(attr), threshold);
    })
    .def("filteringByExtinction", [](AttributeFiltersPybind &self, FloatArray attr, int leafToKeep) {
        return self.filteringByExtinctionValue(std::move(attr), leafToKeep);
    })
    .def("saliencyMapByExtinction", [](AttributeFiltersPybind &self, FloatArray attr, int leafToKeep, bool unweighted) {
        return self.saliencyMapByExtinctionValue(std::move(attr), leafToKeep, unweighted);
    }, "attr"_a, "leafToKeep"_a, "unweighted"_a = false)
    .def("getAdaptiveCriterion", &AttributeFiltersPybind::getAdaptiveCriterion);
}

void init_ExtinctionValues(py::module &m){
    using FloatArray = py::array_t<float, py::array::c_style | py::array::forcecast>;

    py::class_<ExtinctionValuesPybind>(m, "ExtinctionValues", py::module_local(false),
        "Extinction value utilities returning dense node-id tuples.")
    .def(py::init<MorphologicalTreePybindPtr, FloatArray>())
    .def(py::init<std::shared_ptr<WeightedMorphologicalTree>, FloatArray>())
    .def("filtering", &ExtinctionValuesPybind::filtering)
    .def("saliencyMap", &ExtinctionValuesPybind::saliencyMap, "leafToKeep"_a, "unweighted"_a = true)
    .def("getExtinctionValues", &ExtinctionValuesPybind::getExtinctionValuesPy,
        "Return extinction tuples as (leafNodeId, cutoffNodeId, value).");
    
}


void init_AdjacencyRelation(py::module &m){
    	py::class_<AdjacencyRelation>(m, "AdjacencyRelation", py::module_local(false))
        .def(py::init<int, int, double>())
        .def_property_readonly("size", &AdjacencyRelation::getSize )
        .def_property_readonly("radius", &AdjacencyRelation::getRadius )
        .def("getAdjPixels", py::overload_cast<int, int>( &AdjacencyRelation::getAdjPixels ));
}


void init_UltimateAttributeOpening(py::module &m){
    using FloatArray = py::array_t<float, py::array::c_style | py::array::forcecast>;

    	py::class_<UltimateAttributeOpeningPybind>(m, "UltimateAttributeOpening", py::module_local(false))
        .def(py::init<MorphologicalTreePybindPtr, FloatArray>())
        .def(py::init<std::shared_ptr<WeightedMorphologicalTree>, FloatArray>())
        .def("execute", py::overload_cast<int>(&UltimateAttributeOpeningPybind::execute))
        .def("executeWithMSER", &UltimateAttributeOpeningPybind::executeWithMSER)
        .def("getMaxContrastImage", &UltimateAttributeOpeningPybind::getMaxContrastImage)
        .def("getAssociatedImage", &UltimateAttributeOpeningPybind::getAssociatedImage)
        .def("getAssociatedColoredImage", &UltimateAttributeOpeningPybind::getAssociatedColorImage);
}

void init_AttributeOpeningPrimitivesFamily(py::module &m){
    using FloatArray = py::array_t<float, py::array::c_style | py::array::forcecast>;

    	py::class_<AttributeOpeningPrimitivesFamilyPybind, std::shared_ptr<AttributeOpeningPrimitivesFamilyPybind>>(m, "AttributeOpeningPrimitivesFamily", py::module_local(false))
        .def(py::init<MorphologicalTreePybindPtr, FloatArray, float>())
        .def(py::init<MorphologicalTreePybindPtr, FloatArray, float, int>())
        .def(py::init<std::shared_ptr<WeightedMorphologicalTree>, FloatArray, float>())
        .def(py::init<std::shared_ptr<WeightedMorphologicalTree>, FloatArray, float, int>())
        .def_property_readonly("numPrimitives", &AttributeOpeningPrimitivesFamilyPybind::getNumPrimitives)
        .def("getPrimitive", &AttributeOpeningPrimitivesFamilyPybind::getPrimitive)
        .def_property_readonly("restOfImage", &AttributeOpeningPrimitivesFamilyPybind::getRestOfNumpyImage)
        .def("getNodesWithMaximumCriterium", &AttributeOpeningPrimitivesFamilyPybind::getNodesWithMaximumCriterium)
        .def("getThresholdsPrimitive", &AttributeOpeningPrimitivesFamilyPybind::getThresholdsPrimitive);

}

void init_ComponentTreeAdjust(py::module &m) {
    namespace adjust = mmcfilters::adjust;

    py::enum_<adjust::CasfComponentTreesAttribute>(m, "CasfComponentTreesAttribute", py::module_local(false))
        .value("AREA", adjust::CasfComponentTreesAttribute::AREA)
        .value("BOUNDING_BOX_WIDTH", adjust::CasfComponentTreesAttribute::BOUNDING_BOX_WIDTH)
        .value("BOUNDING_BOX_HEIGHT", adjust::CasfComponentTreesAttribute::BOUNDING_BOX_HEIGHT)
        .value("BOUNDING_BOX_DIAGONAL", adjust::CasfComponentTreesAttribute::BOUNDING_BOX_DIAGONAL)
        .export_values();

    py::class_<DualMinMaxTreeIncrementalFilterPybind, std::shared_ptr<DualMinMaxTreeIncrementalFilterPybind>>(
        m,
        "DualMinMaxTreeIncrementalFilter",
        py::module_local(false),
        "Incremental updater for paired weighted min-tree/max-tree component trees.")
        .def(py::init<std::shared_ptr<WeightedMorphologicalTree>, std::shared_ptr<WeightedMorphologicalTree>>(),
            "minTree"_a,
            "maxTree"_a)
        .def("pruneMaxTreeAndUpdateMinTree",
            &DualMinMaxTreeIncrementalFilterPybind::pruneMaxTreeAndUpdateMinTree,
            "nodesToPrune"_a,
            "Update the min-tree incrementally and then prune the requested max-tree subtrees.")
        .def("pruneMinTreeAndUpdateMaxTree",
            &DualMinMaxTreeIncrementalFilterPybind::pruneMinTreeAndUpdateMaxTree,
            "nodesToPrune"_a,
            "Update the max-tree incrementally and then prune the requested min-tree subtrees.")
        .def_property_readonly("minTree", &DualMinMaxTreeIncrementalFilterPybind::minTree)
        .def_property_readonly("maxTree", &DualMinMaxTreeIncrementalFilterPybind::maxTree);

    py::class_<adjust::CasfComponentTrees, std::shared_ptr<adjust::CasfComponentTrees>>(
        m,
        "CasfComponentTrees",
        py::module_local(false),
        "Connected alternating sequential filter based on paired component-tree adjustment.")
        .def(py::init([](UInt8InputArray input, adjust::CasfComponentTreesAttribute attribute, double radius) {
                return std::make_shared<adjust::CasfComponentTrees>(imageFromArray(input), attribute, radius);
            }),
            "input"_a,
            "attribute"_a = adjust::CasfComponentTreesAttribute::AREA,
            "radius"_a = 1.5)
        .def("filter", [](adjust::CasfComponentTrees& self, const std::vector<double>& thresholds) {
            return PybindUtils::toNumpy(self.filter(thresholds));
        }, "thresholds"_a)
        .def_property_readonly("minTree", [](adjust::CasfComponentTrees& self) -> const WeightedMorphologicalTree& {
            return self.minTree();
        }, py::return_value_policy::reference_internal)
        .def_property_readonly("maxTree", [](adjust::CasfComponentTrees& self) -> const WeightedMorphologicalTree& {
            return self.maxTree();
        }, py::return_value_policy::reference_internal)
        .def_property_readonly("attribute", &adjust::CasfComponentTrees::attribute)
        .def("exportMinTree", &adjust::CasfComponentTrees::exportMinTree)
        .def("exportMaxTree", &adjust::CasfComponentTrees::exportMaxTree);

    m.attr("ComponentTreeCasf") = m.attr("CasfComponentTrees");
    m.attr("ComponentTreeCasfAttribute") = m.attr("CasfComponentTreesAttribute");
}

PYBIND11_MODULE(mmcfilters, m) {
    m.doc() = "Morphological tree filters with a NodeId-first Python API.";
    
    init_MorphologicalTree(m);
    init_AttributeComputedIncrementally(m);
    init_ContoursComputedIncrementally(m);
    init_AttributeFilters(m);
    init_ExtinctionValues(m);
    init_AdjacencyRelation(m);

    init_UltimateAttributeOpening(m);
    init_AttributeOpeningPrimitivesFamily(m);
    init_ComponentTreeAdjust(m);

}
