#pragma once

#include "../MorphologicalTree.hpp"
#include "CommittedTreeAccess.hpp"
#include "../../utils/Altitude.hpp"
#include "../../utils/Contract.hpp"

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
 *
 * @param altitude Altitude data indexed by node identifier.
 * @param nodeId Identifier of the node used by the operation.
 * @return The requested altitude from a dense internal-node buffer with bounds checks.
 */
template <AltitudeValue T> inline T exportedHigraAltitudeAt(std::span<const T> altitude, NodeId nodeId) {
    MMCFILTERS_CONTRACT_REQUIRE(nodeId >= 0 && static_cast<size_t>(nodeId) < altitude.size(),
                                throw std::invalid_argument("Altitude access requires a valid internal NodeId."));
    return altitude[static_cast<size_t>(nodeId)];
}

/**
 * @brief Computes the compact Higra node-id layout for a live rooted tree.
 *
 * @details
 * The layout follows the Higra convention used by the public export API: leaves
 * occupy `[0, numProperParts)` and live internal nodes occupy the remaining ids.
 * Globally monotone hierarchies are sorted by altitude direction, with
 * post-order time as the deterministic tie-breaker. For an unconstrained
 * hierarchy, no single altitude direction can order every mixed-polarity
 * branch, so nodes are ordered directly by post-order. Both policies guarantee
 * that every non-root internal node is emitted before its parent, as required
 * by the compact Higra tree layout. The function only computes the id mapping;
 * callers remain responsible for filling parent, altitude, or attribute
 * buffers.
 *
 * @param tree Tree topology used by the operation.
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
    MMCFILTERS_CONTRACT_REQUIRE(altitude.size() == static_cast<size_t>(tree.getNumInternalNodeSlots()),
                                throw std::runtime_error("Altitude buffer size must match the dense internal-node domain."));

    if (tree.getRoot() == InvalidNode || !tree.isAlive(tree.getRoot())) {
        throw std::runtime_error("Cannot export a tree without a valid rooted component.");
    }

    std::vector<NodeId> exportedNodes;
    exportedNodes.reserve(static_cast<size_t>(tree.getNumNodes()));
    for (NodeId nodeId : CommittedTreeAccess::subtree(tree, tree.getRoot())) {
        exportedNodes.push_back(nodeId);
    }

    if (static_cast<int>(exportedNodes.size()) != tree.getNumNodes()) {
        throw std::runtime_error("Cannot export a forest or a tree with detached alive nodes to a compact Higra representation.");
    }

    const NodeId numLeaves = tree.getNumTotalProperParts();
    const NodeId numAliveNodes = static_cast<NodeId>(exportedNodes.size());
    const NodeId numVertices = numLeaves + numAliveNodes;
    auto sortedNodes = exportedNodes;
    std::vector<int> postOrderRank(static_cast<size_t>(tree.getNumInternalNodeSlots()), 0);
    int nextPostOrderRank = 0;
    for (NodeId nodeId : tree.getPostOrderNodes()) {
        postOrderRank[static_cast<size_t>(nodeId)] = nextPostOrderRank++;
    }

    const AltitudeOrder altitudeOrder = tree.getAltitudeOrder();
    if (altitudeOrder == AltitudeOrder::UNCONSTRAINED) {
        std::stable_sort(sortedNodes.begin(), sortedNodes.end(), [&](NodeId lhs, NodeId rhs) {
            return postOrderRank[static_cast<size_t>(lhs)] < postOrderRank[static_cast<size_t>(rhs)];
        });
    } else {
        const bool sortAscendingAltitude = altitudeOrder == AltitudeOrder::DECREASING_FROM_ROOT;
        std::stable_sort(sortedNodes.begin(), sortedNodes.end(), [&](NodeId lhs, NodeId rhs) {
            const T altL = altitude[static_cast<size_t>(lhs)];
            const T altR = altitude[static_cast<size_t>(rhs)];
            if (altL != altR) {
                return sortAscendingAltitude ? altL < altR : altL > altR;
            }
            return postOrderRank[static_cast<size_t>(lhs)] < postOrderRank[static_cast<size_t>(rhs)];
        });
    }

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
