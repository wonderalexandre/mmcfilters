#pragma once

/**
 * @file ResidualTreeCandidatePreparation.hpp
 * @brief Common support and boundary preparation for residual-tree candidates.
 */

#include "FlatZonePartition.hpp"
#include "ResidualTreeCandidateContext.hpp"
#include "../../WeightedMorphologicalTree.hpp"

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
 * @return Whether the candidate support contains the configured exterior seed.
 */
template <AltitudeValue T>
[[nodiscard]] bool prepareResidualTreeCandidate(FlatZonePartition<T>& partition, ResidualTreeCandidateContext& context, NodeId candidateNode,
                                                const WeightedMorphologicalTree<T>& primal, const WeightedMorphologicalTree<T>& dual) {
    const MorphologicalTree& dualTopology = dual.topology();
    const bool dualIsMaxTree = dualTopology.getDescriptiveKind() == MorphologicalTreeKind::MAX_TREE;

    context.boundaryPixelMarks.resetAll();
    context.boundaryOwnerMarks.resetAll();
    context.supportOwnerMarks.resetAll();
    context.boundaryOwners.clear();
    context.boundaryPixels.clear();
    context.supportPixels = {};
    context.supportOwners.clear();
    context.dualExtremalOwner = InvalidNode;
    context.wholeSupportOwner = InvalidNode;

    const auto view = partition.viewForNode(candidateNode, primal);
    context.flatZoneRepresentative = view.representative;
    context.supportPixels = view.supportPixels;
    const auto boundaryPixelsByIncidence = view.externalPixelsByIncidence;

    T dualExtremalAltitude{};
    bool supportOwnerCertifiedFromFlatZone = false;
    if (!context.supportPixels.empty()) {
        const NodeId owner = dualTopology.getProperPartOwner(context.supportPixels.front());
        if (owner != InvalidNode && dualTopology.isAlive(owner) &&
            dualTopology.getNumProperParts(owner) == static_cast<int>(context.supportPixels.size())) {
#ifndef NDEBUG
            for (NodeId pixel : context.supportPixels) {
                assert(dualTopology.getProperPartOwner(pixel) == owner);
            }
#endif
            context.supportOwners.push_back(owner);
            context.dualExtremalOwner = owner;
            context.wholeSupportOwner = owner;
            supportOwnerCertifiedFromFlatZone = true;
        }
    }

    if (!supportOwnerCertifiedFromFlatZone) {
        for (NodeId pixel : context.supportPixels) {
            const NodeId owner = dualTopology.getProperPartOwner(pixel);
            if (owner == InvalidNode || !dualTopology.isAlive(owner)) {
                throw std::runtime_error("Min/max residual support has an invalid opposite-tree owner.");
            }
            if (!context.supportOwnerMarks.isMarked(static_cast<std::size_t>(owner))) {
                context.supportOwnerMarks.mark(static_cast<std::size_t>(owner));
                context.supportOwners.push_back(owner);
                const T ownerAltitude = dual.getAltitude(owner);
                if (context.dualExtremalOwner == InvalidNode || (dualIsMaxTree && ownerAltitude < dualExtremalAltitude) ||
                    (!dualIsMaxTree && ownerAltitude > dualExtremalAltitude)) {
                    context.dualExtremalOwner = owner;
                    dualExtremalAltitude = ownerAltitude;
                }
            }
        }
    }
    if (!supportOwnerCertifiedFromFlatZone && context.supportOwners.size() == 1) {
        const NodeId owner = context.supportOwners.front();
        if (dualTopology.getNumProperParts(owner) == static_cast<int>(context.supportPixels.size())) {
            context.wholeSupportOwner = owner;
        }
    }

    for (NodeId neighbor : boundaryPixelsByIncidence) {
        const std::size_t neighborIndex = static_cast<std::size_t>(neighbor);
        if (!context.boundaryPixelMarks.isMarked(neighborIndex)) {
            context.boundaryPixelMarks.mark(neighborIndex);
            context.boundaryPixels.push_back(neighbor);
        }
        const NodeId owner = dualTopology.getProperPartOwner(neighbor);
        if (owner == InvalidNode || !dualTopology.isAlive(owner)) {
            throw std::runtime_error("Min/max residual boundary has an invalid opposite-tree owner.");
        }
        if (!context.boundaryOwnerMarks.isMarked(static_cast<std::size_t>(owner))) {
            context.boundaryOwnerMarks.mark(static_cast<std::size_t>(owner));
            context.boundaryOwners.push_back(owner);
        }
    }
    if (context.boundaryOwners.empty()) {
        throw std::runtime_error("A non-root min/max leaf has no external boundary.");
    }
    return view.containsExteriorSeed;
}

} // namespace mmcfilters::sdrt::detail
