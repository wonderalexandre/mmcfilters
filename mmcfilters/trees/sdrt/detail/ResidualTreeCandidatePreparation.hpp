#pragma once

/**
 * @file ResidualTreeCandidatePreparation.hpp
 * @brief Common support and boundary preparation for residual-tree candidates.
 */

#include "FlatZonePartition.hpp"
#include "ResidualTreeCandidateContext.hpp"
#include "../../ValuedMorphologicalTree.hpp"

#include <cassert>
#include <cstddef>
#include <stdexcept>

namespace mmcfilters::sdrt::detail {

/**
 * @brief Prepares the support and dual-tree boundary data required to contract one extremum.
 *
 * The operation is independent of saturated eligibility. It derives all
 * candidate support metadata from the current flat-zone partition and fills
 * reusable scratch used by both construction modes.
 *
 * @param partition Current union-only flat-zone partition.
 * @param context Reusable candidate scratch populated by the operation.
 * @param candidateNode Current leaf in the primal component tree.
 * @param primal Component tree that owns `candidateNode`.
 * @param dual Opposite-polarity component tree.
 * @return Whether the candidate support contains the declared infinity pixel.
 */
template <AltitudeValue T>
[[nodiscard]] bool prepareResidualTreeCandidate(FlatZonePartition<T>& partition, ResidualTreeCandidateContext& context, NodeId candidateNode,
                                                const ValuedMorphologicalTree<T>& primal, const ValuedMorphologicalTree<T>& dual) {
    const MorphologicalTree& dualTopology = dual.topology();
    const bool dualIsMaxTree = dualTopology.kind() == MorphologicalTreeKind::MaxTree;

    context.boundaryPixelMarks.resetAll();
    context.boundarySmallestNodeMarks.resetAll();
    context.supportSmallestNodeMarks.resetAll();
    context.boundarySmallestNodes.clear();
    context.boundaryPixels.clear();
    context.supportPixels = {};
    context.supportSmallestNodes.clear();
    context.dualExtremalSmallestNode = InvalidNode;
    context.wholeSupportSmallestNode = InvalidNode;

    const auto view = partition.viewForNode(candidateNode, primal);
    context.flatZoneRepresentative = view.representative;
    context.supportPixels = view.supportPixels;
    const auto boundaryPixelsByIncidence = view.externalPixelsByIncidence;

    T dualExtremalAltitude{};
    bool supportSmallestNodeCertifiedFromFlatZone = false;
    if (!context.supportPixels.empty()) {
        const NodeId smallestNodeId = dualTopology.smallestNode(context.supportPixels.front());
        if (smallestNodeId != InvalidNode && dualTopology.isAlive(smallestNodeId) &&
            dualTopology.properPartCardinality(smallestNodeId) == static_cast<int>(context.supportPixels.size())) {
#ifndef NDEBUG
            for (PixelId pixel : context.supportPixels) {
                assert(dualTopology.smallestNode(pixel) == smallestNodeId);
            }
#endif
            context.supportSmallestNodes.push_back(smallestNodeId);
            context.dualExtremalSmallestNode = smallestNodeId;
            context.wholeSupportSmallestNode = smallestNodeId;
            supportSmallestNodeCertifiedFromFlatZone = true;
        }
    }

    if (!supportSmallestNodeCertifiedFromFlatZone) {
        for (PixelId pixel : context.supportPixels) {
            const NodeId smallestNodeId = dualTopology.smallestNode(pixel);
            if (smallestNodeId == InvalidNode || !dualTopology.isAlive(smallestNodeId)) {
                throw std::runtime_error("Min/max residual support has an invalid opposite-tree smallest node.");
            }
            if (!context.supportSmallestNodeMarks.isMarked(static_cast<std::size_t>(smallestNodeId))) {
                context.supportSmallestNodeMarks.mark(static_cast<std::size_t>(smallestNodeId));
                context.supportSmallestNodes.push_back(smallestNodeId);
                const T smallestNodeAltitude = dual.nodeAltitude(smallestNodeId);
                if (context.dualExtremalSmallestNode == InvalidNode || (dualIsMaxTree && smallestNodeAltitude < dualExtremalAltitude) ||
                    (!dualIsMaxTree && smallestNodeAltitude > dualExtremalAltitude)) {
                    context.dualExtremalSmallestNode = smallestNodeId;
                    dualExtremalAltitude = smallestNodeAltitude;
                }
            }
        }
    }
    if (!supportSmallestNodeCertifiedFromFlatZone && context.supportSmallestNodes.size() == 1) {
        const NodeId smallestNodeId = context.supportSmallestNodes.front();
        if (dualTopology.properPartCardinality(smallestNodeId) == static_cast<int>(context.supportPixels.size())) {
            context.wholeSupportSmallestNode = smallestNodeId;
        }
    }

    for (PixelId neighbor : boundaryPixelsByIncidence) {
        const std::size_t neighborIndex = static_cast<std::size_t>(neighbor);
        if (!context.boundaryPixelMarks.isMarked(neighborIndex)) {
            context.boundaryPixelMarks.mark(neighborIndex);
            context.boundaryPixels.push_back(neighbor);
        }
        const NodeId smallestNodeId = dualTopology.smallestNode(neighbor);
        if (smallestNodeId == InvalidNode || !dualTopology.isAlive(smallestNodeId)) {
            throw std::runtime_error("Min/max residual boundary has an invalid opposite-tree smallest node.");
        }
        if (!context.boundarySmallestNodeMarks.isMarked(static_cast<std::size_t>(smallestNodeId))) {
            context.boundarySmallestNodeMarks.mark(static_cast<std::size_t>(smallestNodeId));
            context.boundarySmallestNodes.push_back(smallestNodeId);
        }
    }
    if (context.boundarySmallestNodes.empty()) {
        throw std::runtime_error("A non-root min/max leaf has no external boundary.");
    }
    return view.containsInfinityPixel;
}

} // namespace mmcfilters::sdrt::detail
