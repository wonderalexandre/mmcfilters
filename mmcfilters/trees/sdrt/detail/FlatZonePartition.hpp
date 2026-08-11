#pragma once

/**
 * @file FlatZonePartition.hpp
 * @brief Union-find partition and boundary storage for residual-tree flat zones.
 */

#include "ResidualTreeCandidateTypes.hpp"
#include "ResidualTreeRegionTypes.hpp"
#include "../../WeightedMorphologicalTree.hpp"
#include "../../../utils/GenerationStampSet.hpp"
#include "../../../utils/Image.hpp"
#include "../../../utils/RegularGridAdjacency2D.hpp"

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcfilters::sdrt::detail {

/**
 * @brief Union-only boundary index for the flat zones of a residual-tree construction.
 *
 * Elementary levelings merge flat zones and never split them. This index
 * therefore remains independent of the internal rewiring of the mutable
 * max-tree and min-tree. Boundary bags keep one outside-pixel entry per
 * original crossing incidence and discard entries that become internal
 * lazily after a union.
 */
template <AltitudeValue T> class FlatZonePartition {
  public:
    /** @brief Residual and image altitude type. */
    using altitude_t = T;
    /** @brief Shared image type supplying pixel altitudes. */
    using image_ptr_t = ImagePtr<altitude_t>;
    /** @brief Weighted component-tree type tracked by boundary storage. */
    using tree_t = WeightedMorphologicalTree<altitude_t>;
    /** @brief Candidate properties derived from the current flat-zone partition. */
    using CandidateDescriptor = ResidualTreeCandidateDescriptor;

    /** @brief Non-owning view of a current merged flat zone. */
    struct View {
        NodeId representative = InvalidNode;               ///< Dense flat-zone representative.
        std::span<const NodeId> supportPixels;             ///< Current zone support.
        std::span<const NodeId> externalPixelsByIncidence; ///< Outside incidences.
        NodeId stableSpatialKey = InvalidNode;             ///< Tie-breaking key.
        bool containsExteriorSeed = false;                 ///< Whether the zone owns the configured exterior seed.
    };

  private:
    /** @brief Stores one current flat zone and its boundary incidences. */
    struct Zone {
        /** @brief Stores the level. */
        altitude_t level{};
        /** @brief Stores the members. */
        std::vector<NodeId> members;
        /** @brief Stores the external pixels by incidence. */
        std::vector<NodeId> externalPixelsByIncidence;
        /** @brief Stores the stable spatial key. */
        NodeId stableSpatialKey = InvalidNode;
        /** @brief Indicates whether the flat zone contains the infinity pixel. */
        bool containsExteriorSeed = false;
    };

    /** @brief Stores the parent. */
    std::vector<NodeId> parent_;
    /** @brief Stores the zones. */
    std::vector<Zone> zones_;

    /**
     * @brief Finds the current disjoint-set representative of a pixel.
     *
     * @param pixel Row-major pixel identifier.
     * @return The matching node identifier, or the operation-specific sentinel when absent.
     */
    [[nodiscard]] NodeId findRepresentative(NodeId pixel) {
        NodeId representative = pixel;
        while (parent_[static_cast<std::size_t>(representative)] != representative) {
            representative = parent_[static_cast<std::size_t>(representative)];
        }
        while (pixel != representative) {
            const NodeId next = parent_[static_cast<std::size_t>(pixel)];
            parent_[static_cast<std::size_t>(pixel)] = representative;
            pixel = next;
        }
        return representative;
    }

    /**
     * @brief Unites two equal-altitude pixels during partition initialization.
     *
     * @param first First value processed by the operation.
     * @param second Second value processed by the operation.
     * @param componentSize Initialization-only component sizes used for balanced unions.
     * @return Representative of the initialized disjoint set.
     */
    NodeId uniteEqualLevelPixels(NodeId first, NodeId second, std::vector<std::size_t>& componentSize) {
        first = findRepresentative(first);
        second = findRepresentative(second);
        if (first == second) {
            return first;
        }
        if (componentSize[static_cast<std::size_t>(first)] < componentSize[static_cast<std::size_t>(second)]) {
            std::swap(first, second);
        }
        parent_[static_cast<std::size_t>(second)] = first;
        componentSize[static_cast<std::size_t>(first)] += componentSize[static_cast<std::size_t>(second)];
        return first;
    }

    /**
     * @brief Moves entries from the smaller list into the larger list.
     *
     * @param target List that receives all entries.
     * @param source List consumed and cleared by the merge.
     */
    static void mergeListsSmallToLarge(std::vector<NodeId>& target, std::vector<NodeId>& source) {
        if (target.size() < source.size()) {
            target.swap(source);
        }
        target.insert(target.end(), source.begin(), source.end());
        source.clear();
    }

    /**
     * @brief Unites two current flat zones and combines their metadata.
     *
     * @param first First value processed by the operation.
     * @param second Second value processed by the operation.
     * @param level Altitude level processed by the operation.
     * @return Representative of the merged current flat zone.
     */
    NodeId uniteFlatZonesAtLevel(NodeId first, NodeId second, altitude_t level) {
        first = findRepresentative(first);
        second = findRepresentative(second);
        if (first == second) {
            zones_[static_cast<std::size_t>(first)].level = level;
            return first;
        }
        if (zones_[static_cast<std::size_t>(first)].members.size() < zones_[static_cast<std::size_t>(second)].members.size()) {
            std::swap(first, second);
        }
        Zone& target = zones_[static_cast<std::size_t>(first)];
        Zone& source = zones_[static_cast<std::size_t>(second)];
        mergeListsSmallToLarge(target.members, source.members);
        mergeListsSmallToLarge(target.externalPixelsByIncidence, source.externalPixelsByIncidence);
        target.stableSpatialKey = std::min(target.stableSpatialKey, source.stableSpatialKey);
        target.containsExteriorSeed = target.containsExteriorSeed || source.containsExteriorSeed;
        target.level = level;
        parent_[static_cast<std::size_t>(second)] = first;
        return first;
    }

    /**
     * @brief Resolves and validates the flat zone represented by one component-tree proper part.
     *
     * @param tree Component tree that owns the proper part.
     * @param nodeId Node whose direct proper part identifies the flat zone.
     * @return Current flat-zone representative, or `InvalidNode` for an empty proper part.
     */
    [[nodiscard]] NodeId representativeForNode(const tree_t& tree, NodeId nodeId) {
        const MorphologicalTree& topology = tree.topology();
        const auto properParts = topology.getProperParts(nodeId);
        const auto first = properParts.begin();
        if (first == properParts.end()) {
            return InvalidNode;
        }
        const NodeId representative = findRepresentative(*first);
        const Zone& zone = zones_[static_cast<std::size_t>(representative)];
        if (zone.members.size() != static_cast<std::size_t>(topology.getNumProperParts(nodeId))) {
            throw std::runtime_error("The flat-zone partition diverged from the component-tree proper part.");
        }
        return representative;
    }

    /**
     * @brief Removes lazily retained internal incidences and returns one current zone view.
     *
     * @param representative Any pixel or representative belonging to the requested flat zone.
     * @return Non-owning view of the current support, boundary incidences, and metadata.
     */
    [[nodiscard]] View currentView(NodeId representative) {
        representative = findRepresentative(representative);
        Zone& zone = zones_[static_cast<std::size_t>(representative)];
        auto output = zone.externalPixelsByIncidence.begin();
        for (NodeId neighbor : zone.externalPixelsByIncidence) {
            if (findRepresentative(neighbor) != representative) {
                *output++ = neighbor;
            }
        }
        zone.externalPixelsByIncidence.erase(output, zone.externalPixelsByIncidence.end());
        return View{representative, zone.members, zone.externalPixelsByIncidence, zone.stableSpatialKey, zone.containsExteriorSeed};
    }

  public:
    /**
     * @brief Returns the current flat-zone representative of every pixel.
     *
     * @return Region representative indexed by row-major pixel.
     */
    [[nodiscard]] std::vector<RegionId> representativesByPixel() {
        std::vector<RegionId> representativeByPixel(parent_.size(), InvalidRegion);
        for (NodeId pixel = 0; pixel < static_cast<NodeId>(parent_.size()); ++pixel) {
            representativeByPixel[static_cast<std::size_t>(pixel)] = static_cast<RegionId>(findRepresentative(pixel));
        }
        return representativeByPixel;
    }

    /**
     * @brief Returns the current flat-zone representative of a pixel.
     *
     * @param pixel Row-major pixel identifier.
     * @return Current disjoint-set representative of `pixel`.
     */
    [[nodiscard]] NodeId representativeOf(NodeId pixel) {
        if (pixel < 0 || pixel >= static_cast<NodeId>(parent_.size())) {
            throw std::out_of_range("Flat-zone representative query lies outside the active domain.");
        }
        return findRepresentative(pixel);
    }

    /**
     * @brief Returns the current flat zone containing one pixel.
     * @param pixel Pixel belonging to the requested flat zone.
     * @return Current support, external incidences, and stable metadata.
     */
    [[nodiscard]] View viewForPixel(NodeId pixel) {
        if (pixel < 0 || pixel >= static_cast<NodeId>(parent_.size())) {
            throw std::out_of_range("Flat-zone view query lies outside the active domain.");
        }
        return currentView(pixel);
    }

    /**
     * @brief Initializes the partition without an exterior seed.
     *
     * This overload is used by unrestricted residual-tree construction.
     *
     * @param image Image processed by the operation.
     * @param adjacency Adjacency relation used by the operation.
     */
    void initialize(const image_ptr_t& image, const RegularGridAdjacency2D& adjacency) { initialize(image, adjacency, InvalidNode); }

    /**
     * @brief Initializes the partition and tracks one configured exterior seed.
     *
     * @param image Image processed by the operation.
     * @param adjacency Adjacency relation used by the operation.
     * @param exteriorSeed Row-major pixel used as the exterior seed.
     */
    void initialize(const image_ptr_t& image, const RegularGridAdjacency2D& adjacency, NodeId exteriorSeed) {
        if (exteriorSeed != InvalidNode && (exteriorSeed < 0 || exteriorSeed >= image->getSize())) {
            throw std::out_of_range("Flat-zone partition exterior seed lies outside the image domain.");
        }
        parent_.clear();
        zones_.clear();
        const std::size_t numPixels = static_cast<std::size_t>(image->getSize());
        parent_.resize(numPixels);
        zones_.resize(numPixels);
        std::iota(parent_.begin(), parent_.end(), NodeId{0});
        std::vector<std::size_t> componentSize(numPixels, 1);

        for (NodeId pixel = 0; pixel < image->getSize(); ++pixel) {
            for (NodeId neighbor : adjacency.getForwardNeighborIndices(pixel)) {
                if ((*image)[pixel] == (*image)[neighbor]) {
                    uniteEqualLevelPixels(pixel, neighbor, componentSize);
                }
            }
        }
        for (NodeId pixel = 0; pixel < image->getSize(); ++pixel) {
            const NodeId representative = findRepresentative(pixel);
            Zone& zone = zones_[static_cast<std::size_t>(representative)];
            zone.level = (*image)[pixel];
            zone.members.push_back(pixel);
            if (zone.stableSpatialKey == InvalidNode) {
                zone.stableSpatialKey = pixel;
            } else {
                zone.stableSpatialKey = std::min(zone.stableSpatialKey, pixel);
            }
            zone.containsExteriorSeed = zone.containsExteriorSeed || pixel == exteriorSeed;
        }
        for (NodeId pixel = 0; pixel < image->getSize(); ++pixel) {
            const NodeId pixelRepresentative = findRepresentative(pixel);
            for (NodeId neighbor : adjacency.getForwardNeighborIndices(pixel)) {
                const NodeId neighborRepresentative = findRepresentative(neighbor);
                if (pixelRepresentative == neighborRepresentative) {
                    continue;
                }
                zones_[static_cast<std::size_t>(pixelRepresentative)].externalPixelsByIncidence.push_back(neighbor);
                zones_[static_cast<std::size_t>(neighborRepresentative)].externalPixelsByIncidence.push_back(pixel);
            }
        }
    }

    /**
     * @brief Returns the proper-part boundary required for the current candidate.
     *
     * @param candidateNode Component-tree node whose boundary is requested.
     * @param primal Component tree that owns `candidateNode`.
     * @return View of the current flat-zone members and external boundary incidences.
     */
    [[nodiscard]] View viewForNode(NodeId candidateNode, const tree_t& primal) {
        const MorphologicalTree& topology = primal.topology();
        const NodeId representative = representativeForNode(primal, candidateNode);
        if (representative == InvalidNode) {
            throw std::runtime_error("A min/max residual candidate has an empty proper part.");
        }
#ifndef NDEBUG
        const Zone& zone = zones_[static_cast<std::size_t>(representative)];
        for (NodeId pixel : zone.members) {
            if (topology.getProperPartOwner(pixel) != candidateNode) {
                throw std::runtime_error("The flat-zone partition has a foreign proper part.");
            }
        }
#endif
        return currentView(representative);
    }

    /**
     * @brief Computes the partition descriptor of a current component-tree node.
     *
     * @param tree Tree processed by the operation.
     * @param nodeId Identifier of the node processed by the operation.
     * @return Descriptor containing area, stable spatial key, and infinity membership.
     */
    [[nodiscard]] CandidateDescriptor describeNode(const tree_t& tree, NodeId nodeId) {
        const NodeId representative = representativeForNode(tree, nodeId);
        if (representative == InvalidNode) {
            return CandidateDescriptor{};
        }
        const Zone& zone = zones_[static_cast<std::size_t>(representative)];
        return CandidateDescriptor{static_cast<int>(zone.members.size()), zone.stableSpatialKey, zone.containsExteriorSeed};
    }

    /**
     * @brief Collects adjacent flat-zone representatives at the merging level.
     *
     * @param selectedRepresentative Current flat-zone representative selected for contraction.
     * @param targetLevel Target altitude level of the merge.
     * @param boundaryPixelsByIncidence External boundary pixels, repeated by incidence.
     * @param marks Generation-stamped set used to deduplicate collected representatives.
     * @param mergeRepresentatives Representatives collected for the merge.
     */
    void collectAdjacentRepresentativesAtLevel(NodeId selectedRepresentative, altitude_t targetLevel, std::span<const NodeId> boundaryPixelsByIncidence,
                                               GenerationStampSet& marks, std::vector<NodeId>& mergeRepresentatives) {
        marks.resetAll();
        mergeRepresentatives.clear();
        selectedRepresentative = findRepresentative(selectedRepresentative);
        for (NodeId pixel : boundaryPixelsByIncidence) {
            const NodeId representative = findRepresentative(pixel);
            if (representative == selectedRepresentative || zones_[static_cast<std::size_t>(representative)].level != targetLevel ||
                marks.isMarked(static_cast<std::size_t>(representative))) {
                continue;
            }
            marks.mark(static_cast<std::size_t>(representative));
            mergeRepresentatives.push_back(representative);
        }
        if (mergeRepresentatives.empty()) {
            throw std::runtime_error("An elementary leveling found no adjacent flat zone at its merging level.");
        }
    }

    /**
     * @brief Merges current flat zones at a common altitude level.
     *
     * @param selectedRepresentative Current flat-zone representative selected for contraction.
     * @param mergeRepresentatives Representatives collected for the merge.
     * @param targetLevel Target altitude level of the merge.
     * @return Representative of the flat zone after all requested merges.
     */
    [[nodiscard]] NodeId mergeFlatZonesAtLevel(NodeId selectedRepresentative, std::span<const NodeId> mergeRepresentatives, altitude_t targetLevel) {
        NodeId representative = findRepresentative(selectedRepresentative);
        for (NodeId neighborRepresentative : mergeRepresentatives) {
            representative = uniteFlatZonesAtLevel(representative, neighborRepresentative, targetLevel);
        }
        zones_[static_cast<std::size_t>(representative)].level = targetLevel;
        return representative;
    }
};

} // namespace mmcfilters::sdrt::detail
