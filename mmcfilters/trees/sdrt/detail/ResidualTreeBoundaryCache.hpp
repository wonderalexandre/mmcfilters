#pragma once

/**
 * @file ResidualTreeBoundaryCache.hpp
 * @brief Recomputed and incremental boundary storage for residual-tree events.
 */

#include "ResidualTreeCandidateAgenda.hpp"
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

/** @brief Boundary-storage types shared by residual-tree construction for one altitude type. */
template <AltitudeValue T> struct ResidualTreeBoundaryTypes {
    /** @brief Residual and image altitude type. */
    using altitude_t = T;
    /** @brief Shared image type supplying pixel altitudes. */
    using image_ptr_t = ImagePtr<altitude_t>;
    /** @brief Weighted component-tree type tracked by boundary storage. */
    using tree_t = WeightedMorphologicalTree<altitude_t>;
    /** @brief Candidate properties required by incremental flat-zone boundaries. */
    using CandidateDescriptor = typename ResidualTreeAgendaTypes<altitude_t>::CandidateDescriptor;

    /** @brief Recomputes one proper part and its external incidences on demand. */
    class RecomputedProperPartBoundary {
      public:
        /** @brief Scratch storage valid until the next acquisition. */
        struct Bag {
            std::vector<NodeId> supportPixels;             ///< Current direct proper part.
            std::vector<NodeId> externalPixelsByIncidence; ///< Outside incidences.
        };

      private:
        /** @brief Stores the bag. */
        Bag bag_;

      public:
        /**
         * @brief Returns the proper-part boundary required for the current candidate.
         *
         * @param nodeId Identifier of the node processed by the operation.
         * @param tree Tree processed by the operation.
         * @param adjacency Adjacency relation used by the operation.
         * @return A reference to the requested object.
         */
        [[nodiscard]] const Bag& acquire(NodeId nodeId, const tree_t& tree, const RegularGridAdjacency2D& adjacency) {
            const MorphologicalTree& topology = tree.topology();
            bag_.supportPixels.clear();
            bag_.externalPixelsByIncidence.clear();
            for (NodeId pixel : topology.getProperParts(nodeId)) {
                bag_.supportPixels.push_back(pixel);
                for (NodeId neighbor : adjacency.getNeighborIndices(pixel)) {
                    if (topology.getProperPartOwner(neighbor) != nodeId) {
                        bag_.externalPixelsByIncidence.push_back(neighbor);
                    }
                }
            }
            return bag_;
        }
    };

    /**
     * @brief Union-only boundary index for the flat zones of the current image.
     *
     * Elementary levelings merge flat zones and never split them. This index
     * therefore remains independent of the internal rewiring of the mutable
     * max-tree and min-tree. Boundary bags keep one outside-pixel entry per
     * original crossing incidence and discard entries that become internal
     * lazily after a union.
     */
    class IncrementalFlatZoneBoundaryCache {
      public:
        /** @brief Non-owning view of a current merged flat zone. */
        struct View {
            NodeId root = InvalidNode;                         ///< Dense flat-zone root.
            std::span<const NodeId> supportPixels;             ///< Current zone support.
            std::span<const NodeId> externalPixelsByIncidence; ///< Outside incidences.
            NodeId stableSpatialKey = InvalidNode;             ///< Tie-breaking key.
            bool containsInfinity = false;                     ///< Whether the zone owns the exterior seed.
        };

      private:
        /** @brief Stores one incremental flat-zone boundary set. */
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
            bool containsInfinity = false;
        };

        /** @brief Indicates whether incremental flat-zone caching is enabled. */
        bool enabled_ = true;
        /** @brief Stores the parent. */
        std::vector<NodeId> parent_;
        /** @brief Stores the size. */
        std::vector<std::size_t> size_;
        /** @brief Stores the zones. */
        std::vector<Zone> zones_;

        /**
         * @brief Finds the current disjoint-set representative of a pixel.
         *
         * @param pixel Row-major pixel identifier.
         * @return The matching node identifier, or the operation-specific sentinel when absent.
         */
        [[nodiscard]] NodeId findRoot(NodeId pixel) {
            NodeId root = pixel;
            while (parent_[static_cast<std::size_t>(root)] != root) {
                root = parent_[static_cast<std::size_t>(root)];
            }
            while (pixel != root) {
                const NodeId next = parent_[static_cast<std::size_t>(pixel)];
                parent_[static_cast<std::size_t>(pixel)] = root;
                pixel = next;
            }
            return root;
        }

        /**
         * @brief Unites two equal-altitude pixels during cache initialization.
         *
         * @param first First value processed by the operation.
         * @param second Second value processed by the operation.
         * @return Representative of the initialized disjoint set.
         */
        NodeId uniteInitial(NodeId first, NodeId second) {
            first = findRoot(first);
            second = findRoot(second);
            if (first == second) {
                return first;
            }
            if (size_[static_cast<std::size_t>(first)] < size_[static_cast<std::size_t>(second)]) {
                std::swap(first, second);
            }
            parent_[static_cast<std::size_t>(second)] = first;
            size_[static_cast<std::size_t>(first)] += size_[static_cast<std::size_t>(second)];
            return first;
        }

        /**
         * @brief Moves boundary entries from the smaller set into the larger set.
         *
         * @param target Boundary list that receives all entries.
         * @param source Boundary list consumed and cleared by the merge.
         */
        static void mergeSmallToLarge(std::vector<NodeId>& target, std::vector<NodeId>& source) {
            if (target.size() < source.size()) {
                target.swap(source);
            }
            target.insert(target.end(), source.begin(), source.end());
            source.clear();
        }

        /**
         * @brief Unites two current flat zones and preserves their boundary cache.
         *
         * @param first First value processed by the operation.
         * @param second Second value processed by the operation.
         * @param level Altitude level processed by the operation.
         * @return Representative of the merged current flat zone.
         */
        NodeId uniteCurrent(NodeId first, NodeId second, altitude_t level) {
            first = findRoot(first);
            second = findRoot(second);
            if (first == second) {
                zones_[static_cast<std::size_t>(first)].level = level;
                return first;
            }
            if (zones_[static_cast<std::size_t>(first)].members.size() < zones_[static_cast<std::size_t>(second)].members.size()) {
                std::swap(first, second);
            }
            Zone& target = zones_[static_cast<std::size_t>(first)];
            Zone& source = zones_[static_cast<std::size_t>(second)];
            mergeSmallToLarge(target.members, source.members);
            mergeSmallToLarge(target.externalPixelsByIncidence, source.externalPixelsByIncidence);
            target.stableSpatialKey = std::min(target.stableSpatialKey, source.stableSpatialKey);
            target.containsInfinity = target.containsInfinity || source.containsInfinity;
            target.level = level;
            parent_[static_cast<std::size_t>(second)] = first;
            size_[static_cast<std::size_t>(first)] += size_[static_cast<std::size_t>(second)];
            return first;
        }

      public:
        /**
         * @brief Constructs a `IncrementalFlatZoneBoundaryCache` instance.
         *
         * @param enabled Whether the optional strategy is enabled.
         */
        explicit IncrementalFlatZoneBoundaryCache(bool enabled = true) : enabled_(enabled) {}

        /**
         * @brief Returns the current flat-zone representative of every pixel.
         *
         * @return Region representative indexed by row-major pixel.
         */
        [[nodiscard]] std::vector<detail::RegionId> currentRegionByPixel() {
            std::vector<detail::RegionId> regionByPixel(parent_.size(), detail::InvalidRegion);
            for (NodeId pixel = 0; pixel < static_cast<NodeId>(parent_.size()); ++pixel) {
                regionByPixel[static_cast<std::size_t>(pixel)] = static_cast<detail::RegionId>(findRoot(pixel));
            }
            return regionByPixel;
        }

        /**
         * @brief Returns the current flat-zone root of a pixel.
         *
         * @param pixel Row-major pixel identifier.
         * @return Current disjoint-set representative of `pixel`.
         */
        [[nodiscard]] NodeId currentRoot(NodeId pixel) {
            if (!enabled_ || pixel < 0 || pixel >= static_cast<NodeId>(parent_.size())) {
                throw std::out_of_range("Flat-zone root query lies outside the active domain.");
            }
            return findRoot(pixel);
        }

        /**
         * @brief Initializes the data structure from the current component trees.
         *
         * @param image Image processed by the operation.
         * @param adjacency Adjacency relation used by the operation.
         * @param infinityPixel Row-major pixel used as the exterior seed.
         */
        void initialize(const image_ptr_t& image, const RegularGridAdjacency2D& adjacency, NodeId infinityPixel) {
            parent_.clear();
            size_.clear();
            zones_.clear();
            if (!enabled_) {
                return;
            }
            const std::size_t numPixels = static_cast<std::size_t>(image->getSize());
            parent_.resize(numPixels);
            size_.assign(numPixels, 1);
            zones_.resize(numPixels);
            std::iota(parent_.begin(), parent_.end(), NodeId{0});

            for (NodeId pixel = 0; pixel < image->getSize(); ++pixel) {
                for (NodeId neighbor : adjacency.getForwardNeighborIndices(pixel)) {
                    if ((*image)[pixel] == (*image)[neighbor]) {
                        uniteInitial(pixel, neighbor);
                    }
                }
            }
            for (NodeId pixel = 0; pixel < image->getSize(); ++pixel) {
                const NodeId root = findRoot(pixel);
                Zone& zone = zones_[static_cast<std::size_t>(root)];
                zone.level = (*image)[pixel];
                zone.members.push_back(pixel);
                if (zone.stableSpatialKey == InvalidNode) {
                    zone.stableSpatialKey = pixel;
                } else {
                    zone.stableSpatialKey = std::min(zone.stableSpatialKey, pixel);
                }
                zone.containsInfinity = zone.containsInfinity || pixel == infinityPixel;
            }
            for (NodeId pixel = 0; pixel < image->getSize(); ++pixel) {
                const NodeId pixelRoot = findRoot(pixel);
                for (NodeId neighbor : adjacency.getForwardNeighborIndices(pixel)) {
                    const NodeId neighborRoot = findRoot(neighbor);
                    if (pixelRoot == neighborRoot) {
                        continue;
                    }
                    zones_[static_cast<std::size_t>(pixelRoot)].externalPixelsByIncidence.push_back(neighbor);
                    zones_[static_cast<std::size_t>(neighborRoot)].externalPixelsByIncidence.push_back(pixel);
                }
            }
        }

        /**
         * @brief Returns the proper-part boundary required for the current candidate.
         *
         * @param candidateNode Component-tree node whose boundary is requested.
         * @param primal Component tree that owns `candidateNode`.
         * @return View of the cached flat-zone members and external boundary incidences.
         */
        [[nodiscard]] View acquire(NodeId candidateNode, const tree_t& primal) {
            if (!enabled_) {
                throw std::logic_error("The incremental flat-zone boundary cache is disabled.");
            }
            const MorphologicalTree& topology = primal.topology();
            const auto properParts = topology.getProperParts(candidateNode);
            const auto first = properParts.begin();
            if (first == properParts.end()) {
                throw std::runtime_error("A min/max residual candidate has an empty proper part.");
            }
            const NodeId root = findRoot(*first);
            Zone& zone = zones_[static_cast<std::size_t>(root)];
            if (zone.members.size() != static_cast<std::size_t>(topology.getNumProperParts(candidateNode))) {
                throw std::runtime_error("The incremental flat-zone boundary cache diverged from the component-tree proper part.");
            }
#ifndef NDEBUG
            for (NodeId pixel : zone.members) {
                if (topology.getProperPartOwner(pixel) != candidateNode) {
                    throw std::runtime_error("The incremental flat-zone boundary cache has a foreign proper part.");
                }
            }
#endif
            auto output = zone.externalPixelsByIncidence.begin();
            for (NodeId neighbor : zone.externalPixelsByIncidence) {
                if (findRoot(neighbor) == root) {
                    continue;
                }
                *output++ = neighbor;
            }
            zone.externalPixelsByIncidence.erase(output, zone.externalPixelsByIncidence.end());
            return View{root, zone.members, zone.externalPixelsByIncidence, zone.stableSpatialKey, zone.containsInfinity};
        }

        /**
         * @brief Computes the cached descriptor of a current component-tree node.
         *
         * @param tree Tree processed by the operation.
         * @param nodeId Identifier of the node processed by the operation.
         * @return Descriptor containing area, stable spatial key, and infinity membership.
         */
        [[nodiscard]] CandidateDescriptor describeNode(const tree_t& tree, NodeId nodeId) {
            const MorphologicalTree& topology = tree.topology();
            const auto properParts = topology.getProperParts(nodeId);
            const auto first = properParts.begin();
            if (first == properParts.end()) {
                return CandidateDescriptor{};
            }
            const NodeId root = findRoot(*first);
            const Zone& zone = zones_[static_cast<std::size_t>(root)];
            if (zone.members.size() != static_cast<std::size_t>(topology.getNumProperParts(nodeId))) {
                throw std::runtime_error("Flat-zone metadata diverged from a component-tree proper part.");
            }
            return CandidateDescriptor{static_cast<int>(zone.members.size()), zone.stableSpatialKey, zone.containsInfinity};
        }

        /**
         * @brief Collects the current flat-zone roots that must be merged.
         *
         * @param selectedRoot Current flat-zone root selected for contraction.
         * @param targetLevel Target altitude level of the merge.
         * @param boundaryPixelsByIncidence External boundary pixels, repeated by incidence.
         * @param marks Generation-stamped set used to deduplicate collected roots.
         * @param mergeRoots Roots collected or updated by the merge.
         */
        void collectMergeRoots(NodeId selectedRoot, altitude_t targetLevel, std::span<const NodeId> boundaryPixelsByIncidence, GenerationStampSet& marks,
                               std::vector<NodeId>& mergeRoots) {
            marks.resetAll();
            mergeRoots.clear();
            selectedRoot = findRoot(selectedRoot);
            for (NodeId pixel : boundaryPixelsByIncidence) {
                const NodeId root = findRoot(pixel);
                if (root == selectedRoot || zones_[static_cast<std::size_t>(root)].level != targetLevel || marks.isMarked(static_cast<std::size_t>(root))) {
                    continue;
                }
                marks.mark(static_cast<std::size_t>(root));
                mergeRoots.push_back(root);
            }
            if (mergeRoots.empty()) {
                throw std::runtime_error("An elementary leveling found no adjacent flat zone at its merging level.");
            }
        }

        /**
         * @brief Merges current flat zones at a common altitude level.
         *
         * @param selectedRoot Current flat-zone root selected for contraction.
         * @param mergeRoots Roots collected or updated by the merge.
         * @param targetLevel Target altitude level of the merge.
         * @return Representative of the flat zone after all requested merges.
         */
        [[nodiscard]] NodeId mergeAtLevel(NodeId selectedRoot, std::span<const NodeId> mergeRoots, altitude_t targetLevel) {
            NodeId root = findRoot(selectedRoot);
            for (NodeId neighborRoot : mergeRoots) {
                root = uniteCurrent(root, neighborRoot, targetLevel);
            }
            zones_[static_cast<std::size_t>(root)].level = targetLevel;
            return root;
        }
    };
};

} // namespace mmcfilters::sdrt::detail
