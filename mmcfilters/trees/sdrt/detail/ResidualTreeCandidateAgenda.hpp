#pragma once

/**
 * @file ResidualTreeCandidateAgenda.hpp
 * @brief Deterministic candidate agenda shared by residual-tree modes.
 */

#include "../ResidualTreePolicies.hpp"
#include "../../WeightedMorphologicalTree.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <set>
#include <stdexcept>
#include <vector>

namespace mmcfilters::sdrt::detail {

/** @brief Mode-independent candidate, ordering, and agenda types for one residual altitude type. */
template <AltitudeValue T> struct ResidualTreeAgendaTypes {
    /** @brief Weighted component-tree type processed by the agenda. */
    using tree_t = WeightedMorphologicalTree<T>;

    /** @brief Enumerates the supported polarity values. */
    enum class Polarity : std::uint8_t { Max, Min };

    /** @brief Describes one residual-tree candidate and its deterministic ordering keys. */
    struct Candidate {
        /** @brief Stores the area. */
        int area = 0;
        /** @brief Stores the stable spatial key. */
        NodeId stableSpatialKey = InvalidNode;
        /** @brief Stores the polarity. */
        Polarity polarity = Polarity::Max;
        /** @brief Stores the node identifier. */
        NodeId nodeId = InvalidNode;
    };

    /** @brief Caches the derived properties of one residual-tree candidate. */
    struct CandidateDescriptor {
        /** @brief Stores the area. */
        int area = 0;
        /** @brief Stores the stable spatial key. */
        NodeId stableSpatialKey = InvalidNode;
        /** @brief Indicates whether the candidate support contains the infinity pixel. */
        bool containsInfinity = false;
    };

    /** @brief Orders residual-tree candidates according to the configured tie policy. */
    struct CandidateLess {
        /** @brief Stores the policy. */
        SdrtTiePolicy policy = SdrtTiePolicy::ContrastInvariantSpatial;

        /**
         * @brief Returns the deterministic ordering rank of a component-tree polarity.
         *
         * @param polarity Component-tree polarity to process.
         * @return `0` for max-tree candidates and `1` for min-tree candidates.
         */
        [[nodiscard]] static int polarityRank(Polarity polarity) noexcept { return polarity == Polarity::Max ? 0 : 1; }

        /**
         * @brief Compares two candidates according to the configured deterministic order.
         *
         * @param lhs Left-hand value of the comparison.
         * @param rhs Right-hand value of the comparison.
         * @return `true` when the left-hand candidate precedes the right-hand candidate.
         */
        [[nodiscard]] bool operator()(const Candidate& lhs, const Candidate& rhs) const noexcept {
            if (lhs.area != rhs.area) {
                return lhs.area < rhs.area;
            }
            if (policy == SdrtTiePolicy::MaxBeforeMinThenSpatial && lhs.polarity != rhs.polarity) {
                return polarityRank(lhs.polarity) < polarityRank(rhs.polarity);
            }
            if (lhs.stableSpatialKey != rhs.stableSpatialKey) {
                return lhs.stableSpatialKey < rhs.stableSpatialKey;
            }
            if (lhs.polarity != rhs.polarity) {
                return polarityRank(lhs.polarity) < polarityRank(rhs.polarity);
            }
            return lhs.nodeId < rhs.nodeId;
        }
    };

    /** @brief Stores the active agenda entry associated with one component-tree node. */
    struct AgendaSlot {
        /** @brief Indicates whether the slot currently owns an agenda candidate. */
        bool active = false;
        /** @brief Stores the candidate. */
        Candidate candidate;
    };

    /** @brief Maintains the deterministic agenda of eligible component-tree leaves. */
    class CandidateAgenda {
      private:
        /** @brief Stores the policy. */
        SdrtTiePolicy policy_;
        /** @brief Stores the infinity pixel. */
        NodeId infinityPixel_ = InvalidNode;
        /** @brief Indicates whether the infinity pixel must be excluded from candidates. */
        bool excludeInfinity_ = true;
        /** @brief Stores the candidates. */
        std::set<Candidate, CandidateLess> candidates_;
        /** @brief Stores the max slots. */
        std::vector<AgendaSlot> maxSlots_;
        /** @brief Stores the min slots. */
        std::vector<AgendaSlot> minSlots_;

        /**
         * @brief Returns the agenda-slot array for one component-tree polarity.
         *
         * @param polarity Component-tree polarity to process.
         * @return A reference to the requested object.
         */
        [[nodiscard]] std::vector<AgendaSlot>& slots(Polarity polarity) { return polarity == Polarity::Max ? maxSlots_ : minSlots_; }

        /**
         * @brief Returns the agenda-slot array for one component-tree polarity.
         *
         * @param polarity Component-tree polarity to process.
         * @return A reference to the requested object.
         */
        [[nodiscard]] const std::vector<AgendaSlot>& slots(Polarity polarity) const { return polarity == Polarity::Max ? maxSlots_ : minSlots_; }

        /**
         * @brief Checks whether a node identifier indexes an agenda-slot array.
         *
         * @param side Agenda-slot array to validate.
         * @param nodeId Identifier of the node processed by the operation.
         * @return `true` when the node identifier is within the slot array; otherwise `false`.
         */
        [[nodiscard]] static bool validSlot(const std::vector<AgendaSlot>& side, NodeId nodeId) noexcept {
            return nodeId >= 0 && nodeId < static_cast<NodeId>(side.size());
        }

        /**
         * @brief Removes one node from the active candidate agenda.
         *
         * @param polarity Component-tree polarity to process.
         * @param nodeId Identifier of the node processed by the operation.
         */
        void remove(Polarity polarity, NodeId nodeId) {
            auto& side = slots(polarity);
            if (!validSlot(side, nodeId)) {
                return;
            }
            auto& slot = side[static_cast<std::size_t>(nodeId)];
            if (!slot.active) {
                return;
            }
            if (candidates_.erase(slot.candidate) != 1) {
                throw std::runtime_error("Min/max residual agenda violates its slot/set invariant.");
            }
            slot = AgendaSlot{};
        }

        /**
         * @brief Indexes all current regional extrema of one component tree.
         *
         * @param tree Tree processed by the operation.
         * @param polarity Component-tree polarity to process.
         */
        void indexTree(const tree_t& tree, Polarity polarity) {
            for (NodeId nodeId : tree.topology().getAliveNodeIds()) {
                update(tree, polarity, nodeId);
            }
        }

      public:
        /**
         * @brief Constructs a `CandidateAgenda` instance.
         *
         * @param policy Policy controlling the operation.
         * @param infinityPixel Row-major pixel used as the exterior seed.
         * @param excludeInfinity Whether the infinity pixel must be excluded from the agenda.
         */
        CandidateAgenda(SdrtTiePolicy policy, NodeId infinityPixel, bool excludeInfinity)
            : policy_(policy), infinityPixel_(infinityPixel), excludeInfinity_(excludeInfinity), candidates_(CandidateLess{policy}) {}

        /**
         * @brief Initializes the data structure from the current component trees.
         *
         * @param maxTree Max-tree consumed by the operation.
         * @param minTree Min-tree consumed by the operation.
         */
        void initialize(const tree_t& maxTree, const tree_t& minTree) {
            candidates_ = std::set<Candidate, CandidateLess>(CandidateLess{policy_});
            maxSlots_.assign(static_cast<std::size_t>(maxTree.topology().getNumInternalNodeSlots()), AgendaSlot{});
            minSlots_.assign(static_cast<std::size_t>(minTree.topology().getNumInternalNodeSlots()), AgendaSlot{});
            indexTree(maxTree, Polarity::Max);
            indexTree(minTree, Polarity::Min);
        }

        /**
         * @brief Refreshes the active agenda entry of one component-tree node.
         *
         * @param tree Tree processed by the operation.
         * @param polarity Component-tree polarity to process.
         * @param nodeId Identifier of the node processed by the operation.
         */
        void update(const tree_t& tree, Polarity polarity, NodeId nodeId) {
            auto& side = slots(polarity);
            if (!validSlot(side, nodeId)) {
                return;
            }
            remove(polarity, nodeId);

            const MorphologicalTree& topology = tree.topology();
            if (!topology.isNode(nodeId) || !topology.isAlive(nodeId) || nodeId == topology.getRoot() || !topology.isLeaf(nodeId)) {
                return;
            }
            const NodeId parent = topology.getNodeParent(nodeId);
            if (parent == InvalidNode || parent == nodeId) {
                throw std::runtime_error("Min/max residual leaf violates the parent invariant.");
            }
            const int area = topology.getNumProperParts(nodeId);
            if (area <= 0) {
                return;
            }
            NodeId spatial = std::numeric_limits<NodeId>::max();
            bool containsInfinity = false;
            for (NodeId pixel : topology.getProperParts(nodeId)) {
                spatial = std::min(spatial, pixel);
                containsInfinity = containsInfinity || pixel == infinityPixel_;
            }
            if (spatial == std::numeric_limits<NodeId>::max()) {
                throw std::runtime_error("Min/max residual candidate has an empty proper part.");
            }
            if (excludeInfinity_ && containsInfinity) {
                return;
            }

            Candidate candidate{area, spatial, polarity, nodeId};
            const auto [_, inserted] = candidates_.insert(candidate);
            if (!inserted) {
                throw std::runtime_error("Min/max residual agenda received a duplicate candidate.");
            }
            auto& slot = side[static_cast<std::size_t>(nodeId)];
            slot.active = true;
            slot.candidate = candidate;
        }

        /**
         * @brief Updates an agenda entry from a precomputed candidate descriptor.
         *
         * @param tree Tree processed by the operation.
         * @param polarity Component-tree polarity to process.
         * @param nodeId Identifier of the node processed by the operation.
         * @param descriptor Descriptor associated with the candidate.
         */
        void updateKnown(const tree_t& tree, Polarity polarity, NodeId nodeId, const CandidateDescriptor& descriptor) {
            auto& side = slots(polarity);
            if (!validSlot(side, nodeId)) {
                return;
            }
            remove(polarity, nodeId);

            const MorphologicalTree& topology = tree.topology();
            if (!topology.isNode(nodeId) || !topology.isAlive(nodeId) || nodeId == topology.getRoot() || !topology.isLeaf(nodeId)) {
                return;
            }
            if (descriptor.area <= 0 || descriptor.stableSpatialKey == InvalidNode || topology.getNumProperParts(nodeId) != descriptor.area) {
                throw std::runtime_error("Min/max residual flat-zone metadata violates the agenda invariant.");
            }
            if (excludeInfinity_ && descriptor.containsInfinity) {
                return;
            }

            Candidate candidate{descriptor.area, descriptor.stableSpatialKey, polarity, nodeId};
            const auto [_, inserted] = candidates_.insert(candidate);
            if (!inserted) {
                throw std::runtime_error("Min/max residual agenda received a duplicate metadata candidate.");
            }
            auto& slot = side[static_cast<std::size_t>(nodeId)];
            slot.active = true;
            slot.candidate = candidate;
        }

        /**
         * @brief Removes a rejected node from the active candidate agenda.
         *
         * @param polarity Component-tree polarity to process.
         * @param nodeId Identifier of the node processed by the operation.
         */
        void reject(Polarity polarity, NodeId nodeId) { remove(polarity, nodeId); }

        /**
         * @brief Selects the first candidate in deterministic agenda order.
         *
         * @return The next candidate, or `std::nullopt` when the agenda is empty.
         */
        [[nodiscard]] std::optional<Candidate> select() const {
            if (candidates_.empty()) {
                return std::nullopt;
            }
            return *candidates_.begin();
        }

        /**
         * @brief Checks whether a node currently has an active agenda entry.
         *
         * @param polarity Component-tree polarity to process.
         * @param nodeId Identifier of the node processed by the operation.
         * @return `true` when the node has an active agenda entry; otherwise `false`.
         */
        [[nodiscard]] bool contains(Polarity polarity, NodeId nodeId) const {
            const auto& side = slots(polarity);
            return validSlot(side, nodeId) && side[static_cast<std::size_t>(nodeId)].active;
        }
    };
};

} // namespace mmcfilters::sdrt::detail
