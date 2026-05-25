#pragma once

#include "../MorphologicalTree.hpp"
#include "../../utils/Altitude.hpp"

#include <algorithm>
#include <numeric>
#include <span>
#include <stdexcept>
#include <vector>

namespace mmcfilters::detail {

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
 * attribute pipeline a friend of `WeightedMorphologicalTree<std::uint8_t>`.
 */
struct ExportedHigraLayout {
    /// Live internal nodes in the order assigned to compact Higra internal ids.
    std::vector<NodeId> sortedNodes;
    /// Dense internal `NodeId` slot -> compact Higra node id. Dead slots map to `InvalidNode`.
    std::vector<NodeId> nodeToHigra;
    /// Proper parts in compact Higra leaf order. This is currently row-major proper-part id order.
    std::vector<NodeId> properParts;
    /// Number of leaf/proper-part ids at the beginning of the compact Higra domain.
    NodeId numLeaves = 0;
    /// Total compact Higra domain size: `numLeaves + number of live internal nodes`.
    NodeId numVertices = 0;
};

/**
 * @brief Reads an altitude from a dense internal-node buffer with bounds checks.
 */
template<AltitudeValue T>
inline T exportedHigraAltitudeAt(std::span<const T> altitude, NodeId nodeId) {
    if (nodeId < 0 || static_cast<size_t>(nodeId) >= altitude.size()) {
        throw std::invalid_argument("Altitude access requires a valid internal NodeId.");
    }
    return altitude[static_cast<size_t>(nodeId)];
}

/**
 * @brief Computes the compact Higra node-id layout for a live rooted tree.
 *
 * @details
 * The layout follows the Higra convention used by the public export API: leaves
 * occupy `[0, numProperParts)` and live internal nodes occupy the remaining ids.
 * Internal nodes are sorted by altitude direction, with post-order time as the
 * deterministic tie-breaker. The function only computes the id mapping; callers
 * remain responsible for filling parent, altitude, or attribute buffers.
 *
 * @throws std::runtime_error when the topology is not one connected rooted live
 * component or when the altitude buffer does not cover the dense internal-node
 * domain.
 */
template<AltitudeValue T>
inline ExportedHigraLayout computeExportedHigraLayout(const MorphologicalTree& tree, std::span<const T> altitude) {
    if (altitude.size() != static_cast<size_t>(tree.getNumInternalNodeSlots())) {
        throw std::runtime_error("Altitude buffer size must match the dense internal-node domain.");
    }

    if (tree.getRoot() == InvalidNode || !tree.isAlive(tree.getRoot())) {
        throw std::runtime_error("Cannot export a tree without a valid rooted component.");
    }

    std::vector<NodeId> exportedNodes;
    exportedNodes.reserve(static_cast<size_t>(tree.getNumNodes()));
    for (NodeId nodeId : tree.getNodeSubtree(tree.getRoot())) {
        exportedNodes.push_back(nodeId);
    }

    if (static_cast<int>(exportedNodes.size()) != tree.getNumNodes()) {
        throw std::runtime_error("Cannot export a forest or a tree with detached alive nodes to a compact Higra representation.");
    }

    const NodeId numLeaves = tree.getNumTotalProperParts();
    const NodeId numAliveNodes = static_cast<NodeId>(exportedNodes.size());
    const NodeId numVertices = numLeaves + numAliveNodes;
    auto sortedNodes = exportedNodes;

    bool sortAscendingAltitude = true;
    switch (tree.getTreeType()) {
        case MorphologicalTreeKind::MAX_TREE:
            sortAscendingAltitude = false;
            break;
        case MorphologicalTreeKind::MIN_TREE:
            sortAscendingAltitude = true;
            break;
        case MorphologicalTreeKind::TREE_OF_SHAPES:
        case MorphologicalTreeKind::SELF_DUAL_RESIDUAL_TREE:
            for (NodeId nodeId : sortedNodes) {
                if (tree.isRoot(nodeId)) {
                    continue;
                }

                const NodeId parentNodeId = tree.getNodeParent(nodeId);
                if (parentNodeId == InvalidNode || !tree.isAlive(parentNodeId)) {
                    throw std::runtime_error("Cannot export a node whose parent is not part of the rooted alive component.");
                }
                if (exportedHigraAltitudeAt(altitude, nodeId) > exportedHigraAltitudeAt(altitude, parentNodeId)) {
                    sortAscendingAltitude = false;
                }
            }
            break;
        default:
            throw std::invalid_argument("Unsupported tree type for compact Higra export.");
    }

    std::stable_sort(
        sortedNodes.begin(),
        sortedNodes.end(),
        [&](NodeId lhs, NodeId rhs) {
            const T altL = exportedHigraAltitudeAt(altitude, lhs);
            const T altR = exportedHigraAltitudeAt(altitude, rhs);
            if (altL != altR) {
                return sortAscendingAltitude ? altL < altR : altL > altR;
            }
            return tree.getNodeTimePostOrder(lhs) < tree.getNodeTimePostOrder(rhs);
        });

    std::vector<NodeId> nodeToHigra(static_cast<size_t>(tree.getNumInternalNodeSlots()), InvalidNode);
    for (NodeId i = 0; i < numAliveNodes; ++i) {
        const NodeId oldNodeId = sortedNodes[static_cast<size_t>(i)];
        nodeToHigra[static_cast<size_t>(oldNodeId)] = numLeaves + i;
    }

    std::vector<NodeId> properParts(static_cast<size_t>(numLeaves), InvalidNode);
    std::iota(properParts.begin(), properParts.end(), NodeId{0});

    return {std::move(sortedNodes), std::move(nodeToHigra), std::move(properParts), numLeaves, numVertices};
}

} // namespace mmcfilters::detail
