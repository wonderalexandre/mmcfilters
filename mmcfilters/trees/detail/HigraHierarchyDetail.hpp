#pragma once

/**
 * @file HigraHierarchyDetail.hpp
 * @brief Higra hierarchy import conversion and compact export layout.
 */

#include "../MorphologicalTree.hpp"
#include "../MorphologicalTreeSemantics.hpp"
#include "TreeKernelAccess.hpp"
#include "../NativeHierarchyView.hpp"
#include "../NativeHierarchyOwned.hpp"
#include "../../utils/Altitude.hpp"
#include "../../utils/Contract.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <optional>
#include <sstream>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace mmcfilters::detail {

// Validation and import from the Higra leaf-first representation.

/**
 * @brief Validates the native buffers produced by Higra import conversion.
 */
class NativeHierarchyValidation {
    /**
     * @brief Validates descriptor.
     *
     * @param hierarchy Native hierarchy data to validate or convert.
     */
    template <AltitudeValue T> static void validateDescriptor(const NativeHierarchyView<T>& hierarchy) {
        if (hierarchy.parent.size() > static_cast<std::size_t>(std::numeric_limits<NodeId>::max()) ||
            hierarchy.smallestNodeMap.size() > static_cast<std::size_t>(std::numeric_limits<NodeId>::max())) {
            throw std::invalid_argument("Native hierarchy domain exceeds the NodeId range.");
        }
        if (hierarchy.parent.empty()) {
            throw std::invalid_argument("Native hierarchy requires at least one internal node.");
        }
        if (hierarchy.smallestNodeMap.empty()) {
            throw std::invalid_argument("Native hierarchy requires at least one proper part.");
        }
        if (hierarchy.nodeAltitudes.size() != hierarchy.parent.size()) {
            throw std::invalid_argument("Native hierarchy altitude size must match the internal-node domain.");
        }

        const int numNodes = static_cast<int>(hierarchy.parent.size());
        if (hierarchy.root < 0 || hierarchy.root >= numNodes) {
            throw std::invalid_argument("Native hierarchy requires a valid root node id.");
        }
        if (hierarchy.gridDomain2D && hierarchy.smallestNodeMap.size() != hierarchy.gridDomain2D->size("Native hierarchy 2D pixel domain")) {
            throw std::invalid_argument("Native hierarchy pixel domain must match the attached 2D grid.");
        }

        switch (hierarchy.semantics.kind) {
        case MorphologicalTreeKind::Generic:
        case MorphologicalTreeKind::MaxTree:
        case MorphologicalTreeKind::MinTree:
        case MorphologicalTreeKind::TreeOfShapes:
        case MorphologicalTreeKind::ExtremumResidualTree:
            break;
        default:
            throw std::invalid_argument("Native hierarchy kind is not supported.");
        }

        validateMorphologicalTreeSemantics(hierarchy.semantics);
        const auto requireMatchingGrid = [&hierarchy](const RegularGridAdjacency2D& adjacency, const char* context) {
            if (!hierarchy.gridDomain2D) {
                throw std::invalid_argument(std::string(context) + " requires an attached 2D grid domain.");
            }
            if (adjacency.getNumRows() != hierarchy.gridDomain2D->rows || adjacency.getNumColumns() != hierarchy.gridDomain2D->columns) {
                throw std::invalid_argument(std::string(context) + " must match the attached 2D grid.");
            }
        };
        if (const auto* context = std::get_if<SharedAdjacencyContext>(&hierarchy.semantics.constructionContext)) {
            requireMatchingGrid(context->adjacency, "SharedAdjacencyContext adjacency");
        } else if (const auto* convention = std::get_if<TopographicConvention>(&hierarchy.semantics.constructionContext)) {
            if (convention->infinityPixel < 0) {
                throw std::invalid_argument("TopographicConvention infinity pixel must be non-negative.");
            }
            if (!hierarchy.gridDomain2D) {
                throw std::invalid_argument("TopographicConvention requires an attached 2D grid domain.");
            }
            const std::int64_t extension = convention->domainExtension == TopographicDomainExtension::ExteriorRing ? 1 : -1;
            const std::int64_t activeRows = 2 * static_cast<std::int64_t>(hierarchy.gridDomain2D->rows) + extension;
            const std::int64_t activeColumns = 2 * static_cast<std::int64_t>(hierarchy.gridDomain2D->columns) + extension;
            if (activeRows <= 0 || activeColumns <= 0 || static_cast<std::int64_t>(convention->infinityPixel) >= activeRows * activeColumns) {
                throw std::invalid_argument("TopographicConvention infinity pixel must belong to the active topographic domain.");
            }
            if (const auto* immersion = std::get_if<ComplementaryGridImmersion>(&convention->immersion)) {
                requireMatchingGrid(immersion->complementaryAdjacencies.minAdjacency, "Topographic minimum adjacency");
                requireMatchingGrid(immersion->complementaryAdjacencies.maxAdjacency, "Topographic maximum adjacency");
            }
        }
    }

  public:
    /**
     * @brief Validates structure, finite altitudes, and the declared altitude order
     * before the converted Higra input can enter trusted materialization.
     *
     * @param hierarchy Hierarchy data.
     */
    template <AltitudeValue T> static void validateComplete(const NativeHierarchyView<T>& hierarchy) {
        validateDescriptor(hierarchy);

        const int numNodes = static_cast<int>(hierarchy.parent.size());
        std::vector<int> remainingChildren(hierarchy.parent.size(), 0);

        for (NodeId node = 0; node < numNodes; ++node) {
            const auto nodeIndex = static_cast<std::size_t>(node);
            const T altitude = hierarchy.nodeAltitudes[nodeIndex];
            if constexpr (std::is_floating_point_v<T>) {
                if (!std::isfinite(altitude)) {
                    throw std::invalid_argument("Native hierarchy requires finite node altitudes.");
                }
            }
            const NodeId parent = hierarchy.parent[nodeIndex];
            if (node == hierarchy.root) {
                if (parent != node) {
                    throw std::invalid_argument("Native hierarchy root must be self-parented.");
                }
                continue;
            }
            if (parent < 0 || parent >= numNodes) {
                throw std::invalid_argument("Native hierarchy parent is outside the internal-node domain.");
            }
            if (parent == node) {
                throw std::invalid_argument("Native hierarchy contains a detached self-parented non-root node.");
            }
            const T parentAltitude = hierarchy.nodeAltitudes[static_cast<std::size_t>(parent)];
            if (hierarchy.semantics.nodeAltitudeOrder == NodeAltitudeOrder::Increasing && parentAltitude >= altitude) {
                throw std::runtime_error("Native hierarchy altitudes must be strictly increasing from parent to child.");
            }
            if (hierarchy.semantics.nodeAltitudeOrder == NodeAltitudeOrder::Decreasing && parentAltitude <= altitude) {
                throw std::runtime_error("Native hierarchy altitudes must be strictly decreasing from parent to child.");
            }
            ++remainingChildren[static_cast<std::size_t>(parent)];
        }

        std::vector<std::uint8_t> subtreeHasSupport(hierarchy.parent.size(), std::uint8_t{0});
        for (PixelId pixel = 0; pixel < static_cast<PixelId>(hierarchy.smallestNodeMap.size()); ++pixel) {
            const NodeId smallestNode = hierarchy.smallestNodeMap[static_cast<std::size_t>(pixel)];
            if (smallestNode < 0 || smallestNode >= numNodes) {
                throw std::invalid_argument("Native hierarchy smallest node is outside the internal-node domain.");
            }
            subtreeHasSupport[static_cast<std::size_t>(smallestNode)] = 1;
        }

        std::vector<NodeId> ready;
        ready.reserve(hierarchy.parent.size());
        for (NodeId node = 0; node < numNodes; ++node) {
            if (remainingChildren[static_cast<std::size_t>(node)] == 0) {
                ready.push_back(node);
            }
        }

        std::size_t processedNodes = 0;
        for (std::size_t index = 0; index < ready.size(); ++index) {
            const NodeId node = ready[index];
            if (subtreeHasSupport[static_cast<std::size_t>(node)] == 0) {
                throw std::invalid_argument("Native hierarchy contains a node whose subtree support is empty.");
            }
            ++processedNodes;
            if (node == hierarchy.root) {
                continue;
            }

            const NodeId parent = hierarchy.parent[static_cast<std::size_t>(node)];
            subtreeHasSupport[static_cast<std::size_t>(parent)] = 1;
            int& parentRemaining = remainingChildren[static_cast<std::size_t>(parent)];
            --parentRemaining;
            if (parentRemaining == 0) {
                ready.push_back(parent);
            }
        }
        if (processedNodes != hierarchy.parent.size()) {
            throw std::invalid_argument("Native hierarchy contains a cycle or nodes disconnected from the root.");
        }
    }
};



/**
 * @brief Result of adapting a compact leaf-first Higra hierarchy.
 *
 * `internalNodeOffset` records the affine mapping from dense internal node ids
 * back to the imported external id domain. It is interoperability provenance,
 * not morphological-tree semantics.
 */
template <AltitudeValue T> struct HigraHierarchyImport {
    /** @brief Hierarchy. */
    NativeHierarchyOwned<T> hierarchy;
    /** @brief Dense node identifier of the internal node offset. */
    NodeId internalNodeOffset = InvalidNode;
};

/**
 * @brief Converts Higra parent/altitude arrays into the generic native format.
 *
 * The adapter owns all leaf-first id policy. `MorphologicalTree` subsequently
 * materializes only independent node-parent and smallest-node-map domains.
 * Higra leaves occupy `[0, numPixels)`; internal node ids are converted to the
 * native domain by subtracting `numPixels`.
 *
 * @param parent Parent node.
 * @param higraAltitude Altitude or level.
 * @param rows Number of rows in the domain.
 * @param columns Number of columns in the domain.
 * @param kind Morphological-tree family.
 * @param adjacency Adjacency relation.
 * @return The converted Higra parent/altitude arrays into the generic native format.
 */
template <AltitudeValue T>
[[nodiscard]] HigraHierarchyImport<T> adaptHigraHierarchy(std::span<const NodeId> parent, std::span<const T> higraAltitude, int rows, int columns,
                                                          MorphologicalTreeKind kind, std::optional<RegularGridAdjacency2D> adjacency) {
    if (kind == MorphologicalTreeKind::ExtremumResidualTree) {
        throw std::invalid_argument(
            "Higra import does not accept residual-tree kinds because residual-tree proper parts are not represented by the Higra leaf layout.");
    }
    if (!adjacency && (kind == MorphologicalTreeKind::MaxTree || kind == MorphologicalTreeKind::MinTree)) {
        throw std::invalid_argument("Higra import of max/min trees requires an explicit adjacency relation.");
    }
    if (adjacency && kind == MorphologicalTreeKind::TreeOfShapes) {
        throw std::invalid_argument(
            "Higra tree-of-shapes import cannot derive a complete TopographicConvention from one adjacency relation; omit it to record NoConstructionContext.");
    }

    const GridDomain2D gridDomain{rows, columns};
    MorphologicalTreeConstructionContext constructionContext = NoConstructionContext{};
    if (adjacency) {
        constructionContext = SharedAdjacencyContext{std::move(*adjacency)};
    }
    MorphologicalTreeSemantics semantics = makeMorphologicalTreeSemantics(kind, std::move(constructionContext));
    const int numPixels = static_cast<int>(gridDomain.size("Higra import 2D pixel domain"));
    if (const auto* sharedContext = std::get_if<SharedAdjacencyContext>(&semantics.constructionContext)) {
        if (sharedContext->adjacency.getNumRows() * sharedContext->adjacency.getNumColumns() != numPixels) {
            throw std::invalid_argument("Higra leaf count must match image domain size.");
        }
    }

    if (parent.size() > static_cast<std::size_t>(std::numeric_limits<NodeId>::max())) {
        throw std::invalid_argument("Higra parent array exceeds the supported node-id range.");
    }
    const int numVertices = static_cast<int>(parent.size());
    if (numPixels <= 0 || numPixels >= numVertices) {
        throw std::invalid_argument("Higra parent array must contain leaves followed by at least one internal node.");
    }
    const int numNodes = numVertices - numPixels;
    if (higraAltitude.size() != static_cast<std::size_t>(numVertices)) {
        throw std::invalid_argument("Higra altitude buffer size must match the preserved imported Higra hierarchy.");
    }
    const auto internalNodeId = [numPixels, numVertices](NodeId higraNodeId, const char* invalidIdMessage) -> NodeId {
        if (higraNodeId < numPixels || higraNodeId >= numVertices) {
            throw std::invalid_argument(invalidIdMessage);
        }
        return higraNodeId - numPixels;
    };
    const auto validateFiniteAltitude = [&](NodeId externalNode) {
        if constexpr (std::is_floating_point_v<T>) {
            const long double level = static_cast<long double>(higraAltitude[static_cast<std::size_t>(externalNode)]);
            if (!std::isfinite(level)) {
                std::ostringstream message;
                message << "ValuedMorphologicalTree Higra altitude input requires finite floating-point altitudes; value at index " << externalNode << " is "
                        << level << ".";
                throw std::invalid_argument(message.str());
            }
        }
    };

    std::vector<NodeId> smallestNodeMap(static_cast<std::size_t>(numPixels), InvalidNode);
    for (PixelId pixel = 0; pixel < numPixels; ++pixel) {
        validateFiniteAltitude(pixel);
        smallestNodeMap[static_cast<std::size_t>(pixel)] =
            internalNodeId(parent[static_cast<std::size_t>(pixel)], "Each Higra leaf must point to an internal node.");
    }

    std::vector<NodeId> nodeParent(static_cast<std::size_t>(numNodes), InvalidNode);
    std::vector<T> altitude(static_cast<std::size_t>(numNodes));
    NodeId root = InvalidNode;
    for (NodeId node = 0; node < numNodes; ++node) {
        const NodeId higraNode = numPixels + node;
        const NodeId parentHigraNode = parent[static_cast<std::size_t>(higraNode)];
        validateFiniteAltitude(higraNode);
        altitude[static_cast<std::size_t>(node)] = higraAltitude[static_cast<std::size_t>(higraNode)];

        if (higraNode == parentHigraNode) {
            if (root != InvalidNode) {
                throw std::invalid_argument("A Higra hierarchy must encode exactly one root.");
            }
            nodeParent[static_cast<std::size_t>(node)] = node;
            root = node;
            continue;
        }

        nodeParent[static_cast<std::size_t>(node)] =
            internalNodeId(parentHigraNode, "Each Higra internal node must point to another internal node or to itself.");
    }
    if (root == InvalidNode) {
        throw std::invalid_argument("A Higra hierarchy must encode exactly one root.");
    }

    const NativeHierarchyView<T> view{nodeParent, smallestNodeMap, altitude, root, gridDomain, semantics};
    NativeHierarchyValidation::validateComplete(view);
    return {makeNativeHierarchyOwned<T>(std::move(nodeParent), std::move(smallestNodeMap), std::move(altitude), root, gridDomain, std::move(semantics)),
            numPixels};
}

// Compact layout shared by hierarchy export and attribute projection.

/**
 * @brief Compact Higra id mapping shared by hierarchy export and attribute projection.
 *
 * @details
 * The public APIs expose two related operations:
 * - `TreeAltitudeAlgorithms::exportHigraHierarchy()` exports parent and altitude
 *   arrays;
 * - `AttributeComputation::projectNodeValuesToExportedHigra()`
 *   projects node-indexed attribute buffers to the same compact id domain.
 *
 * This internal layout keeps both operations aligned without making the
 * attribute pipeline a friend of `ValuedMorphologicalTree<std::uint8_t>`.
 */
struct ExportedHigraLayout {
    /// Live internal nodes in the order assigned to compact Higra internal ids.
    std::vector<NodeId> sortedNodes;
    /// Dense internal `NodeId` slot -> compact Higra node id. Dead slots map to `InvalidNode`.
    std::vector<NodeId> nodeToHigra;
    /// Pixel identifiers in compact Higra leaf order (row-major).
    std::vector<PixelId> properParts;
    /// Number of leaf/pixel ids at the beginning of the compact Higra domain.
    int numLeaves = 0;
    /// Total compact Higra domain size: `numLeaves + number of live internal nodes`.
    int numVertices = 0;
};

/**
 * @brief Computes the compact Higra node-id layout for a live rooted tree.
 *
 * @details
 * The layout follows the Higra convention used by the public export API: leaves
 * occupy `[0, numPixels)` and live internal nodes occupy the remaining ids.
 * Globally monotone hierarchies are sorted by altitude direction, with
 * post-order time as the deterministic tie-breaker. The output mapping stores
 * those ranks temporarily, avoiding a separate rank buffer. For an unconstrained
 * hierarchy, no single altitude direction can order every mixed-polarity
 * branch, so nodes are ordered directly by post-order. Both policies guarantee
 * that every non-root internal node is emitted before its parent, as required
 * by the compact Higra tree layout. The function only computes the id mapping;
 * callers remain responsible for filling parent, altitude, or attribute
 * buffers.
 *
 * @param tree Tree topology.
 * @param altitude Altitude data indexed by node identifier.
 * @return The computed compact Higra node-id layout for a live rooted tree.
 *
 * @throws std::logic_error when an edit session is active.
 * @throws std::runtime_error when the topology is not one connected rooted live
 * component or when the altitude buffer does not cover the dense internal-node
 * domain.
 *
 */
template <AltitudeValue T> inline ExportedHigraLayout computeExportedHigraLayout(const MorphologicalTree& tree, std::span<const T> altitude) {
    tree.requireNotEditing("Higra hierarchy export");
    MMCFILTERS_CONTRACT_REQUIRE(altitude.size() == static_cast<size_t>(tree.numInternalNodeSlots()),
                                throw std::runtime_error("Altitude buffer size must match the dense internal-node domain."));

    if (tree.root() == InvalidNode || !tree.isAlive(tree.root())) {
        throw std::runtime_error("Cannot export a tree without a valid rooted component.");
    }

    std::vector<NodeId> sortedNodes;
    sortedNodes.reserve(static_cast<size_t>(tree.numNodes()));
    const NodeAltitudeOrder nodeAltitudeOrder = tree.nodeAltitudeOrder();
    if (nodeAltitudeOrder == NodeAltitudeOrder::Unconstrained) {
        for (NodeId nodeId : tree.postOrder()) {
            sortedNodes.push_back(nodeId);
        }
    } else {
        for (NodeId nodeId : TreeKernelAccess::subtree(tree, tree.root())) {
            sortedNodes.push_back(nodeId);
        }
    }

    if (static_cast<int>(sortedNodes.size()) != tree.numNodes()) {
        throw std::runtime_error("Cannot export a forest or a tree with detached alive nodes to a compact Higra representation.");
    }

    const int numLeaves = tree.numPixels();
    const int numAliveNodes = static_cast<int>(sortedNodes.size());
    const int numVertices = numLeaves + numAliveNodes;
    std::vector<NodeId> nodeToHigra(static_cast<size_t>(tree.numInternalNodeSlots()), InvalidNode);

    if (nodeAltitudeOrder != NodeAltitudeOrder::Unconstrained) {
        // Reuse the output buffer for tie-breaking ranks until ids are assigned.
        NodeId nextPostOrderRank = 0;
        for (NodeId nodeId : tree.postOrder()) {
            nodeToHigra[static_cast<size_t>(nodeId)] = nextPostOrderRank++;
        }
        const bool sortAscendingAltitude = nodeAltitudeOrder == NodeAltitudeOrder::Decreasing;
        std::stable_sort(sortedNodes.begin(), sortedNodes.end(), [&](NodeId lhs, NodeId rhs) {
            const T altL = altitude[static_cast<size_t>(lhs)];
            const T altR = altitude[static_cast<size_t>(rhs)];
            if (altL != altR) {
                return sortAscendingAltitude ? altL < altR : altL > altR;
            }
            return nodeToHigra[static_cast<size_t>(lhs)] < nodeToHigra[static_cast<size_t>(rhs)];
        });
    }

    for (NodeId i = 0; i < numAliveNodes; ++i) {
        const NodeId oldNodeId = sortedNodes[static_cast<size_t>(i)];
        nodeToHigra[static_cast<size_t>(oldNodeId)] = numLeaves + i;
    }

    std::vector<PixelId> properParts(static_cast<size_t>(numLeaves), InvalidPixel);
    std::iota(properParts.begin(), properParts.end(), PixelId{0});

    return {std::move(sortedNodes), std::move(nodeToHigra), std::move(properParts), numLeaves, numVertices};
}

} // namespace mmcfilters::detail
