#include "ModuleBindings.hpp"

#include "MorphologicalTreePybind.hpp"
#include "PybindUtils.hpp"

#include "../mmcfilters/attributes/AttributeComputation.hpp"
#include "../mmcfilters/attributes/AttributeNames.hpp"
#include "../mmcfilters/trees/MorphologicalTreeFactory.hpp"
#include "../mmcfilters/trees/TreeAltitudeAlgorithms.hpp"
#include "../mmcfilters/trees/WeightedMorphologicalTree.hpp"
#include "../mmcfilters/utils/AdjacencyRelation.hpp"
#include "../mmcfilters/utils/Common.hpp"
#include "../mmcfilters/utils/Image.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <memory>
#include <optional>
#include <set>
#include <span>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace mmcfilters::pybindings {

namespace py = pybind11;
using namespace pybind11::literals;

using UInt8InputArray = py::array;
using FloatInputArray = py::array_t<float, py::array::c_style | py::array::forcecast>;

namespace {

template <class Range>
std::vector<NodeId> collectNodeIds(const Range& range) {
    std::vector<NodeId> ids;
    for (NodeId id : range) {
        ids.push_back(id);
    }
    return ids;
}

ImageUInt8Ptr imageFromArray(const UInt8InputArray& input) {
    if (!input.dtype().is(py::dtype::of<uint8_t>())) {
        throw std::invalid_argument("input must be a 2D uint8 array");
    }
    auto buf = input.request();
    if (buf.ndim != 2) {
        throw std::invalid_argument("input must be a 2D uint8 array");
    }
    const int rows = static_cast<int>(buf.shape[0]);
    const int cols = static_cast<int>(buf.shape[1]);
    if (buf.strides[1] != static_cast<py::ssize_t>(sizeof(uint8_t)) ||
        buf.strides[0] != static_cast<py::ssize_t>(cols * sizeof(uint8_t))) {
        throw std::invalid_argument("input must be a C-contiguous 2D uint8 array");
    }
    return ImageUInt8::fromExternal(static_cast<uint8_t*>(buf.ptr), rows, cols);
}

std::uint8_t pythonUInt8AltitudeValue(py::handle value, NodeId nodeId, const char* context) {
    if (py::isinstance<py::bool_>(value) || !PyIndex_Check(value.ptr())) {
        throw std::invalid_argument(std::string(context) + " must be an integer altitude in [0, 255].");
    }

    py::object indexValue = py::reinterpret_steal<py::object>(PyNumber_Index(value.ptr()));
    const long long rawAltitude = py::cast<long long>(indexValue);
    if (rawAltitude < 0 || rawAltitude > 255) {
        std::ostringstream oss;
        oss << context << " requires node altitudes in the uint8 domain [0, 255]; node "
            << nodeId << " has altitude " << rawAltitude << ".";
        throw std::invalid_argument(oss.str());
    }
    const auto altitude = static_cast<std::uint8_t>(rawAltitude);
    (void)TreeAltitudeAlgorithms::requireUInt8AltitudeValue(altitude, nodeId, context);
    return altitude;
}

std::vector<std::uint8_t> pythonUInt8AltitudeVector(py::handle values, const char* context) {
    if (py::isinstance<py::array>(values)) {
        py::array array = py::reinterpret_borrow<py::array>(values);
        if (!array.dtype().is(py::dtype::of<uint8_t>())) {
            throw std::invalid_argument(std::string(context) + " numpy arrays must have dtype uint8.");
        }
        auto buffer = array.request();
        if (buffer.ndim != 1) {
            throw std::invalid_argument(std::string(context) + " must be a 1D uint8 array.");
        }
        if (buffer.strides[0] != static_cast<py::ssize_t>(sizeof(uint8_t))) {
            throw std::invalid_argument(std::string(context) + " must be a C-contiguous 1D uint8 array.");
        }
        const auto* data = static_cast<const uint8_t*>(buffer.ptr);
        std::vector<std::uint8_t> altitude(static_cast<std::size_t>(buffer.shape[0]));
        for (std::size_t i = 0; i < altitude.size(); ++i) {
            altitude[i] = static_cast<std::uint8_t>(data[i]);
        }
        return altitude;
    }

    if (!py::isinstance<py::sequence>(values)) {
        throw std::invalid_argument(std::string(context) + " must be a sequence of integer altitudes in [0, 255].");
    }

    py::sequence sequence = py::reinterpret_borrow<py::sequence>(values);
    std::vector<std::uint8_t> altitude;
    altitude.reserve(static_cast<std::size_t>(py::len(sequence)));
    NodeId nodeId = 0;
    for (py::handle item : sequence) {
        altitude.push_back(pythonUInt8AltitudeValue(item, nodeId, context));
        ++nodeId;
    }
    return altitude;
}

MorphologicalTree& topology(MorphologicalTreePybind& tree) {
    return tree;
}

const MorphologicalTree& topology(const MorphologicalTreePybind& tree) {
    return tree;
}

const MorphologicalTree& topology(const WeightedMorphologicalTree<std::uint8_t>& weighted) {
    return weighted.topology();
}

std::uint8_t getAltitudeValue(const WeightedMorphologicalTree<std::uint8_t>& weighted, NodeId nodeId) {
    return weighted.getAltitude(nodeId);
}

AltitudeDiff<std::uint8_t> getResidueValue(const WeightedMorphologicalTree<std::uint8_t>& weighted, NodeId nodeId) {
    return weighted.getNodeResidue(nodeId);
}

py::array_t<uint8_t> reconstructionImageOf(const WeightedMorphologicalTree<std::uint8_t>& weighted) {
    return PybindUtils::toNumpy(weighted.reconstructionImage());
}

std::pair<std::vector<NodeId>, std::vector<std::uint8_t>> exportHigraHierarchyOf(const WeightedMorphologicalTree<std::uint8_t>& weighted) {
    return weighted.exportHigraHierarchy();
}

std::vector<Attribute> parseProjectionAttributes(py::handle attributes, int valuesPerNode) {
    std::vector<Attribute> parsed;

    try {
        parsed.push_back(py::cast<Attribute>(attributes));
    } catch (const py::cast_error&) {
        if (!py::isinstance<py::sequence>(attributes)) {
            throw std::invalid_argument("attributes must be an Attribute or a sequence of Attribute values");
        }
        py::sequence seq = py::reinterpret_borrow<py::sequence>(attributes);
        parsed.reserve(static_cast<size_t>(py::len(seq)));
        for (py::handle item : seq) {
            parsed.push_back(py::cast<Attribute>(item));
        }
    }

    if (static_cast<int>(parsed.size()) != valuesPerNode) {
        std::ostringstream oss;
        oss << "attributes must contain " << valuesPerNode << " item(s), got " << parsed.size();
        throw std::invalid_argument(oss.str());
    }
    if (std::set<Attribute>(parsed.begin(), parsed.end()).size() != parsed.size()) {
        throw std::invalid_argument("attributes must not contain duplicates");
    }
    return parsed;
}

AttributeNames makeProjectionAttributeNames(const std::vector<Attribute>& attributes) {
    std::unordered_map<Attribute, int> offsets;
    for (int i = 0; i < static_cast<int>(attributes.size()); ++i) {
        offsets[attributes[static_cast<size_t>(i)]] = i;
    }
    return AttributeNames(std::move(offsets));
}

py::array_t<float> projectNodeValuesToExportedHigraOf(
    const WeightedMorphologicalTree<std::uint8_t>& weighted,
    const FloatInputArray& nodeValues,
    py::object attributes) {
    const auto buffer = nodeValues.request();
    const int numNodeSlots = weighted.topology().getNumInternalNodeSlots();
    const int numHigraVertices = weighted.topology().getNumTotalProperParts() + weighted.topology().getNumNodes();

    if (buffer.ndim == 1) {
        PybindUtils::require1DArray(buffer, numNodeSlots, "nodeValues");
        const auto parsedAttributes = parseProjectionAttributes(attributes, 1);
        const AttributeNames attrNames = makeProjectionAttributeNames(parsedAttributes);
        auto projected = AttributeComputation::projectNodeValuesToExportedHigra(
            weighted,
            attrNames,
            std::span<const float>(static_cast<const float*>(buffer.ptr), static_cast<size_t>(numNodeSlots)));
        return PybindUtils::toNumpyOwned(std::move(projected), numHigraVertices);
    }

    if (buffer.ndim == 2) {
        if (buffer.shape[0] != numNodeSlots) {
            std::ostringstream oss;
            oss << "nodeValues must have " << numNodeSlots << " rows, got " << buffer.shape[0];
            throw std::invalid_argument(oss.str());
        }
        const int valuesPerNode = static_cast<int>(buffer.shape[1]);
        const auto parsedAttributes = parseProjectionAttributes(attributes, valuesPerNode);
        const AttributeNames attrNames = makeProjectionAttributeNames(parsedAttributes);
        auto projected = AttributeComputation::projectNodeValuesToExportedHigra(
            weighted,
            attrNames,
            std::span<const float>(
                static_cast<const float*>(buffer.ptr),
                static_cast<size_t>(numNodeSlots) * static_cast<size_t>(valuesPerNode)));
        return PybindUtils::toNumpyOwned2D(
            std::move(projected),
            numHigraVertices,
            valuesPerNode);
    }

    throw std::invalid_argument("nodeValues must be a 1D or 2D float array");
}

template <class TreeLike, class PyClass>
void bindTreeQueryApi(PyClass& cls) {
    cls.def_property_readonly("numInternalNodeSlots", [](TreeLike& self) {
            return topology(self).getNumInternalNodeSlots();
        }, "Size of the dense internal NodeId slot domain.")
        .def_property_readonly("numTotalProperParts", [](TreeLike& self) {
            return topology(self).getNumTotalProperParts();
        }, "Number of proper parts in the image-domain support.")
        .def_property_readonly("numHigraNodes", [](TreeLike& self) {
            return topology(self).getNumHigraNodes();
        }, "Size of the preserved imported Higra node-id domain.")
        .def("getRoot", [](TreeLike& self) {
            return topology(self).getRoot();
        }, "Return the current root node id.")
        .def_property_readonly("root", [](TreeLike& self) {
            return topology(self).getRoot();
        }, "Current root node id.")
        .def_property_readonly("numFreeNodeSlots", [](TreeLike& self) {
            return topology(self).getNumFreeNodeSlots();
        }, "Number of currently free internal-node slots.")
        .def_property_readonly("numLeafNodes", [](TreeLike& self) {
            return topology(self).getNumLeafNodes();
        }, "Number of alive leaf nodes.")
        .def("getAliveNodeIds", [](TreeLike& self) {
            return collectNodeIds(topology(self).getAliveNodeIds());
        }, "Return all alive internal-node ids in the dense node-id domain.")
        .def_property_readonly("aliveNodeIds", [](TreeLike& self) {
            return collectNodeIds(topology(self).getAliveNodeIds());
        }, "Alive internal-node ids in the dense node-id domain.")
        .def_property_readonly("alive_node_ids", [](TreeLike& self) {
            return collectNodeIds(topology(self).getAliveNodeIds());
        }, "Alive internal-node ids in the dense node-id domain.")
        .def("getLeafNodeIds", [](TreeLike& self) {
            return topology(self).getLeaves();
        }, "Return alive leaf node ids in the dense node-id domain.")
        .def_property_readonly("leafNodeIds", [](TreeLike& self) {
            return topology(self).getLeaves();
        }, "Alive leaf node ids in the dense node-id domain.")
        .def_property_readonly("leaf_node_ids", [](TreeLike& self) {
            return topology(self).getLeaves();
        }, "Alive leaf node ids in the dense node-id domain.")
        .def("getChildren", [](TreeLike& self, NodeId nodeId) {
            return collectNodeIds(topology(self).getChildren(nodeId));
        }, "nodeId"_a, "Return the direct children of a node in the dense node-id domain.")
        .def("getNodeNumDescendants", [](TreeLike& self, NodeId nodeId) {
            return topology(self).getNodeNumDescendants(nodeId);
        }, "nodeId"_a, "Return the number of descendants of nodeId.")
        .def("getNodeNumSiblings", [](TreeLike& self, NodeId nodeId) {
            return topology(self).getNodeNumSiblings(nodeId);
        }, "nodeId"_a, "Return the number of siblings of nodeId.")
        .def("getNumProperParts", [](TreeLike& self, NodeId nodeId) {
            return topology(self).getNumProperParts(nodeId);
        }, "nodeId"_a, "Return the number of direct proper parts owned by nodeId.")
        .def("getNodeTimePreOrder", [](TreeLike& self, NodeId nodeId) {
            return topology(self).getNodeTimePreOrder(nodeId);
        }, "nodeId"_a, "Return the preorder timestamp of nodeId.")
        .def("getNodeTimePostOrder", [](TreeLike& self, NodeId nodeId) {
            return topology(self).getNodeTimePostOrder(nodeId);
        }, "nodeId"_a, "Return the postorder timestamp of nodeId.")
        .def("getProperParts", [](TreeLike& self, NodeId nodeId) {
            return collectNodeIds(topology(self).getProperParts(nodeId));
        }, "nodeId"_a, "Return the proper parts owned directly by a node.")
        .def("getConnectedComponent", [](TreeLike& self, NodeId nodeId) {
            auto range = topology(self).getConnectedComponent(nodeId);
            return py::make_iterator(range.begin(), range.end());
        }, py::keep_alive<0, 1>(), "nodeId"_a, "Iterate over all proper parts in the connected component represented by nodeId.")
        .def("reconstructNode", [](TreeLike& self, NodeId nodeId) {
            return MorphologicalTreePybind::reconstructNode(topology(self), nodeId);
        }, "nodeId"_a, "Reconstruct a binary mask for the connected component represented by nodeId.")
        .def("getPostOrderNodes", [](TreeLike& self, std::optional<NodeId> rootNodeId) {
            return rootNodeId.has_value()
                ? collectNodeIds(topology(self).getPostOrderNodes(*rootNodeId))
                : collectNodeIds(topology(self).getPostOrderNodes());
        }, "rootNodeId"_a = std::nullopt,
            "Return post-order traversal node ids under `rootNodeId`, or under the tree root.")
        .def("getIteratorBreadthFirstTraversal", [](TreeLike& self, std::optional<NodeId> rootNodeId) {
            return rootNodeId.has_value()
                ? collectNodeIds(topology(self).getIteratorBreadthFirstTraversal(*rootNodeId))
                : collectNodeIds(topology(self).getIteratorBreadthFirstTraversal());
        }, "rootNodeId"_a = std::nullopt,
            "Return breadth-first traversal node ids under `rootNodeId`, or under the tree root.")
        .def("getPathToRootNodes", [](TreeLike& self, NodeId nodeId) {
            return collectNodeIds(topology(self).getPathToRootNodes(nodeId));
        }, "nodeId"_a,
            "Return the path from `nodeId` to the current root.")
        .def("getPathBetweenNodes", [](TreeLike& self, NodeId sourceNodeId, NodeId targetNodeId) {
            return collectNodeIds(topology(self).getPathBetweenNodes(sourceNodeId, targetNodeId));
        }, "sourceNodeId"_a, "targetNodeId"_a,
            "Return the upward path from `sourceNodeId` toward `targetNodeId`.")
        .def("getNodeSubtree", [](TreeLike& self, NodeId nodeId) {
            return collectNodeIds(topology(self).getNodeSubtree(nodeId));
        }, "nodeId"_a,
            "Return all alive node ids in the subtree rooted at `nodeId`.")
        .def("getDescendants", [](TreeLike& self, NodeId nodeId) {
            return collectNodeIds(topology(self).getDescendants(nodeId));
        }, "nodeId"_a,
            "Return all strict descendant node ids of `nodeId`.")
        .def("getNodeParent", [](TreeLike& self, NodeId nodeId) {
            return topology(self).getNodeParent(nodeId);
        }, "nodeId"_a,
            "Return the parent node id of `nodeId`, or `InvalidNode` when absent.")
        .def("getProperPartOwner", [](TreeLike& self, NodeId properPartId) {
            return topology(self).getProperPartOwner(properPartId);
        }, "properPartId"_a,
            "Return the internal node that directly owns `properPartId`.")
        .def("getHigraNodeId", [](TreeLike& self, NodeId nodeId) {
            return topology(self).getHigraNodeId(nodeId);
        }, "nodeId"_a, "Return the preserved imported Higra node id for a live internal NodeId, or InvalidNode.")
        .def("getNumChildren", [](TreeLike& self, NodeId nodeId) {
            return topology(self).getNumChildren(nodeId);
        }, "nodeId"_a,
            "Return the number of direct children of `nodeId`.")
        .def("getFirstChild", [](TreeLike& self, NodeId nodeId) {
            return topology(self).getFirstChild(nodeId);
        }, "nodeId"_a,
            "Return the first child of `nodeId`, or `InvalidNode`.")
        .def("getNextSibling", [](TreeLike& self, NodeId nodeId) {
            return topology(self).getNextSibling(nodeId);
        }, "nodeId"_a,
            "Return the next sibling of `nodeId`, or `InvalidNode`.")
        .def("isNode", [](TreeLike& self, NodeId nodeId) {
            return topology(self).isNode(nodeId);
        }, "nodeId"_a,
            "Return true when `nodeId` is in the internal-node slot domain.")
        .def("isProperPart", [](TreeLike& self, NodeId nodeId) {
            return topology(self).isProperPart(nodeId);
        }, "nodeId"_a,
            "Return true when the id is in the proper-part domain.")
        .def("isAlive", [](TreeLike& self, NodeId nodeId) {
            return topology(self).isAlive(nodeId);
        }, "nodeId"_a,
            "Return true when `nodeId` is an alive internal node.")
        .def("isRoot", [](TreeLike& self, NodeId nodeId) {
            return topology(self).isRoot(nodeId);
        }, "nodeId"_a,
            "Return true when `nodeId` is the current root.")
        .def("isLeaf", [](TreeLike& self, NodeId nodeId) {
            return topology(self).isLeaf(nodeId);
        }, "nodeId"_a,
            "Return true when `nodeId` is an alive leaf node.")
        .def("hasChild", [](TreeLike& self, NodeId parentId, NodeId childId) {
            return topology(self).hasChild(parentId, childId);
        }, "parentId"_a, "childId"_a,
            "Return true when `childId` is a direct child of `parentId`.")
        .def("pruneNode", [](TreeLike& self, NodeId nodeId) {
            self.pruneNode(nodeId);
        }, "nodeId"_a,
            "Prune one node from the topology, preserving a valid rooted tree.")
        .def("mergeNodeIntoParent", [](TreeLike& self, NodeId nodeId) {
            self.mergeNodeIntoParent(nodeId);
        }, "nodeId"_a,
            "Merge `nodeId` into its parent and reattach descendants/proper parts.")
        .def_property_readonly("treeType", [](TreeLike& self) { return topology(self).getTreeType(); },
            "Semantic kind of the morphological hierarchy.")
        .def_property_readonly("hasAdjacencyRelation", [](TreeLike& self) { return topology(self).hasAdjacencyRelation(); },
            "Whether the tree stores an image adjacency relation.")
        .def_property_readonly("hasTreeOfShapesAdjacencyPolicy", [](TreeLike& self) { return topology(self).hasTreeOfShapesAdjacencyPolicy(); },
            "Whether a tree-of-shapes min/max adjacency policy is available.")
        .def("getTreeOfShapesMinTreeAdjacencyRadius", [](TreeLike& self) {
            return topology(self).getTreeOfShapesMinTreeAdjacencyRadius();
        }, "Return the min-tree adjacency radius used by a tree-of-shapes policy.")
        .def("getTreeOfShapesMaxTreeAdjacencyRadius", [](TreeLike& self) {
            return topology(self).getTreeOfShapesMaxTreeAdjacencyRadius();
        }, "Return the max-tree adjacency radius used by a tree-of-shapes policy.")
        .def("getTreeOfShapesMinTreeAdjacencyRelation", [](TreeLike& self) -> const AdjacencyRelation& {
            const AdjacencyRelation* adjacency = topology(self).getTreeOfShapesMinTreeAdjacencyRelation();
            if (!adjacency) {
                throw std::runtime_error("Tree-of-shapes adjacency policy is not available.");
            }
            return *adjacency;
        }, py::return_value_policy::reference_internal,
            "Return the borrowed min-tree adjacency relation used by a tree of shapes.")
        .def("getTreeOfShapesMaxTreeAdjacencyRelation", [](TreeLike& self) -> const AdjacencyRelation& {
            const AdjacencyRelation* adjacency = topology(self).getTreeOfShapesMaxTreeAdjacencyRelation();
            if (!adjacency) {
                throw std::runtime_error("Tree-of-shapes adjacency policy is not available.");
            }
            return *adjacency;
        }, py::return_value_policy::reference_internal,
            "Return the borrowed max-tree adjacency relation used by a tree of shapes.")
        .def_property_readonly("numRows", [](TreeLike& self) { return topology(self).getNumRowsOfImage(); },
            "Number of rows in the original image domain.")
        .def_property_readonly("numCols", [](TreeLike& self) { return topology(self).getNumColsOfImage(); },
            "Number of columns in the original image domain.")
        .def_property_readonly("numNodes", [](TreeLike& self) { return topology(self).getNumNodes(); },
            "Number of currently alive internal nodes.");

    if constexpr (std::is_same_v<TreeLike, WeightedMorphologicalTree<std::uint8_t>>) {
        cls.def("getAltitude", [](TreeLike& self, NodeId nodeId) {
                return getAltitudeValue(self, nodeId);
            }, "nodeId"_a, "Return the altitude associated with nodeId.")
            .def("getNodeResidue", [](TreeLike& self, NodeId nodeId) {
                return getResidueValue(self, nodeId);
            }, "nodeId"_a, "Return the residue between nodeId and its parent.")
            .def("reconstructionImage", [](TreeLike& self) {
                return reconstructionImageOf(self);
            }, "Reconstruct the current tree into a 2D image using the attached image domain.")
            .def("exportHigraHierarchy", [](TreeLike& self) {
                return exportHigraHierarchyOf(self);
            }, "Export the current rooted tree to a new compact Higra (parent, altitude) representation.");
    }
}

} // namespace

void initMorphologicalTree(py::module_& m) {
    py::enum_<ToSInterpolation>(m, "ToSInterpolation", py::module_local(false))
        .value("SelfDual", ToSInterpolation::SelfDual)
        .value("Min4cMax8c", ToSInterpolation::Min4cMax8c)
        .value("Min8cMax4c", ToSInterpolation::Min8cMax4c)
        .export_values();

    py::enum_<NodeIdSpace>(m, "NodeIdSpace", py::module_local(false))
        .value("MORPHOLOGICAL_TREE", NodeIdSpace::MORPHOLOGICAL_TREE)
        .value("HIGRA", NodeIdSpace::HIGRA)
        .export_values();

    py::enum_<MorphologicalTreeKind>(m, "MorphologicalTreeKind", py::module_local(false))
        .value("MAX_TREE", MorphologicalTreeKind::MAX_TREE)
        .value("MIN_TREE", MorphologicalTreeKind::MIN_TREE)
        .value("TREE_OF_SHAPES", MorphologicalTreeKind::TREE_OF_SHAPES)
        .value("SELF_DUAL_RESIDUAL_TREE", MorphologicalTreeKind::SELF_DUAL_RESIDUAL_TREE)
        .export_values();

    py::class_<MorphologicalTree, std::shared_ptr<MorphologicalTree>>(m, "MorphologicalTreeBase", py::module_local(false));
    auto treeCls = py::class_<MorphologicalTreePybind, std::shared_ptr<MorphologicalTreePybind>>(m, "MorphologicalTree", py::module_local(false),
        "Morphological tree topology type with a NodeId-first query API. Python construction is centralized in MorphologicalTreeFactory, "
        "which returns WeightedMorphologicalTree<std::uint8_t> instances exposing the same topology queries plus weighted operations.");
    bindTreeQueryApi<MorphologicalTreePybind>(treeCls);

    auto weightedCls = py::class_<WeightedMorphologicalTree<std::uint8_t>, std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>>>(m, "WeightedMorphologicalTree", py::module_local(false),
        "Wrapper pairing MorphologicalTree topology with an external dense altitude buffer. "
        "Imports can preserve an original Higra node-id domain; exports always create a new compact Higra domain.");

    py::class_<MorphologicalTreeFactory>(m, "MorphologicalTreeFactory", py::module_local(false))
        .def_static("createMaxTree", [](UInt8InputArray input, double radius) {
            return std::make_shared<WeightedMorphologicalTree<std::uint8_t>>(
                MorphologicalTreeFactory::createMaxTree(imageFromArray(input), radius));
        }, "input"_a, "radius"_a = 1.5,
            "Create a max-tree from a 2D C-contiguous `np.uint8` image.")
        .def_static("createMinTree", [](UInt8InputArray input, double radius) {
            return std::make_shared<WeightedMorphologicalTree<std::uint8_t>>(
                MorphologicalTreeFactory::createMinTree(imageFromArray(input), radius));
        }, "input"_a, "radius"_a = 1.5,
            "Create a min-tree from a 2D C-contiguous `np.uint8` image.")
        .def_static("createTreeOfShapes", [](UInt8InputArray input, ToSInterpolation interpolation, int infinitySeedRow, int infinitySeedCol) {
            return std::make_shared<WeightedMorphologicalTree<std::uint8_t>>(
                MorphologicalTreeFactory::createTreeOfShapes(imageFromArray(input), interpolation, infinitySeedRow, infinitySeedCol));
        },
            "input"_a,
            "interpolation"_a = ToSInterpolation::SelfDual,
            "infinitySeedRow"_a = ToSDefaultInfinityRow,
            "infinitySeedCol"_a = ToSDefaultInfinityCol,
            "Create a tree of shapes from a 2D C-contiguous `np.uint8` image.")
        .def_static("createFromHigraParent", [](const std::vector<NodeId>& parent, py::object altitudeInput, int rows, int cols, MorphologicalTreeKind kind, std::optional<double> radius) {
            const std::vector<std::uint8_t> altitude = pythonUInt8AltitudeVector(
                altitudeInput,
                "MorphologicalTreeFactory.createFromHigraParent altitude");
            if (parent.size() != altitude.size()) {
                throw std::invalid_argument("parent and altitude must have the same size");
            }

            return std::make_shared<WeightedMorphologicalTree<std::uint8_t>>(
                MorphologicalTreeFactory::createFromHigraParent(
                    std::span<const NodeId>(parent),
                    std::span<const std::uint8_t>(altitude),
                    rows,
                    cols,
                    kind,
                    radius ? std::optional<AdjacencyRelation>(std::in_place, rows, cols, *radius) : std::nullopt));
        },
            "parent"_a,
            "altitude"_a,
            "rows"_a,
            "cols"_a,
            "kind"_a,
            "radius"_a = py::none(),
            "Create a weighted tree from an imported static Higra parent/altitude representation [leaves | internal nodes]. "
            "The imported Higra node-id domain is preserved until the tree is edited.");

    weightedCls
        .def("setAltitude", [](WeightedMorphologicalTree<std::uint8_t>& tree, NodeId nodeId, py::object altitudeInput) {
            if (!tree.topology().isNode(nodeId)) {
                throw std::invalid_argument("invalid NodeId for altitude update");
            }
            const std::uint8_t altitude = pythonUInt8AltitudeValue(
                altitudeInput,
                nodeId,
                "WeightedMorphologicalTree.setAltitude");
            tree.setAltitude(nodeId, altitude);
        }, "nodeId"_a, "altitude"_a, "Set one node altitude while preserving max-tree/min-tree altitude order.")
        .def("setAltitudeBuffer", [](WeightedMorphologicalTree<std::uint8_t>& tree, py::object altitudeInput) {
            const std::vector<std::uint8_t> altitude = pythonUInt8AltitudeVector(
                altitudeInput,
                "WeightedMorphologicalTree.setAltitudeBuffer");
            PybindUtils::requireVectorSize(altitude, static_cast<std::size_t>(tree.topology().getNumInternalNodeSlots()), "altitude");
            tree.setAltitudeBuffer(altitude);
        }, "altitude"_a, "Replace the dense altitude buffer indexed by internal NodeId after validating altitude order.")
        .def_property(
            "altitude",
            [](const WeightedMorphologicalTree<std::uint8_t>& tree) {
                return AltitudeBuffer<std::uint8_t>(tree.getAltitudeBuffer());
            },
            [](WeightedMorphologicalTree<std::uint8_t>& tree, py::object altitudeInput) {
                const std::vector<std::uint8_t> altitude = pythonUInt8AltitudeVector(
                    altitudeInput,
                    "WeightedMorphologicalTree.altitude");
                PybindUtils::requireVectorSize(altitude, static_cast<std::size_t>(tree.topology().getNumInternalNodeSlots()), "altitude");
                tree.setAltitudeBuffer(altitude);
            },
            "Dense altitude buffer indexed by internal NodeId. Assignments validate max-tree/min-tree altitude order.")
        .def("validateAltitudeBufferShape", static_cast<void (WeightedMorphologicalTree<std::uint8_t>::*)() const>(&WeightedMorphologicalTree<std::uint8_t>::validateAltitudeBufferShape),
            "Validate that the dense altitude buffer covers every internal node slot.")
        .def("validateMonotoneAltitude", static_cast<void (WeightedMorphologicalTree<std::uint8_t>::*)() const>(&WeightedMorphologicalTree<std::uint8_t>::validateMonotoneAltitude),
            "Validate the max-tree/min-tree altitude order where applicable.")
        .def("projectNodeValuesToExportedHigra", &projectNodeValuesToExportedHigraOf,
            "nodeValues"_a,
            "attributes"_a,
            "Project a node-indexed scalar or 2D attribute buffer to the compact Higra layout produced by exportHigraHierarchy().")
        .def("project_node_values_to_exported_higra", &projectNodeValuesToExportedHigraOf,
            "nodeValues"_a,
            "attributes"_a,
            "Project a node-indexed scalar or 2D attribute buffer to the compact Higra layout produced by exportHigraHierarchy().");

    bindTreeQueryApi<WeightedMorphologicalTree<std::uint8_t>>(weightedCls);
}

} // namespace mmcfilters::pybindings
