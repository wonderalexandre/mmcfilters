#pragma once

/**
 * @file ResidualTreeCandidateAgenda.hpp
 * @brief Deterministic candidate ordering shared by residual-tree modes.
 */

#include "ResidualTreeCandidateTypes.hpp"
#include "../ResidualTreePolicies.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <set>
#include <stdexcept>
#include <vector>

namespace mmcfilters::sdrt::detail {

/** @brief Candidate, ordering, and agenda types shared by all residual altitude types. */
struct ResidualTreeAgendaTypes {
    /** @brief Enumerates the supported component-tree polarities. */
    enum class Polarity : std::uint8_t { Max, Min };

    /** @brief Candidate entry stored in deterministic agenda order. */
    struct Candidate {
        int area = 0;                                  ///< Current flat-zone area.
        NodeId stableSpatialKey = InvalidNode;         ///< Stable spatial tie-breaking key.
        Polarity polarity = Polarity::Max;             ///< Component-tree polarity.
        NodeId nodeId = InvalidNode;                   ///< Current component-tree node.
    };

    /** @brief Neutral candidate metadata supplied by the flat-zone partition. */
    using CandidateDescriptor = ResidualTreeCandidateDescriptor;

    /** @brief Orders candidates according to the configured deterministic tie policy. */
    struct CandidateLess {
        SdrtTiePolicy policy = SdrtTiePolicy::ContrastInvariantSpatial; ///< Configured tie policy.

        /**
         * @brief Maps a polarity to its deterministic secondary ordering rank.
         * @param polarity Polarity to rank.
         * @return Zero for max-tree candidates and one for min-tree candidates.
         */
        [[nodiscard]] static int polarityRank(Polarity polarity) noexcept { return polarity == Polarity::Max ? 0 : 1; }

        /**
         * @brief Compares two candidates according to the configured tie policy.
         * @param lhs Left candidate.
         * @param rhs Right candidate.
         * @return `true` when `lhs` precedes `rhs` in agenda order.
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

    /** @brief Active agenda entry associated with one component-tree node slot. */
    struct AgendaSlot {
        bool active = false; ///< Whether this slot currently owns a candidate.
        Candidate candidate; ///< Candidate stored in the ordered set.
    };

    /** @brief Maintains deterministic ordering and membership of current extrema. */
    class CandidateAgenda {
      private:
        SdrtTiePolicy policy_;                           ///< Configured deterministic tie policy.
        bool excludeExteriorSeed_ = true;               ///< Whether exterior-containing candidates are excluded.
        std::set<Candidate, CandidateLess> candidates_; ///< Ordered active candidates.
        std::vector<AgendaSlot> maxSlots_;               ///< Candidate slots for the max-tree.
        std::vector<AgendaSlot> minSlots_;               ///< Candidate slots for the min-tree.

        /**
         * @brief Selects the mutable slot array for one polarity.
         * @param polarity Component-tree polarity.
         * @return Mutable slot array associated with `polarity`.
         */
        [[nodiscard]] std::vector<AgendaSlot>& slots(Polarity polarity) { return polarity == Polarity::Max ? maxSlots_ : minSlots_; }

        /**
         * @brief Selects the immutable slot array for one polarity.
         * @param polarity Component-tree polarity.
         * @return Immutable slot array associated with `polarity`.
         */
        [[nodiscard]] const std::vector<AgendaSlot>& slots(Polarity polarity) const {
            return polarity == Polarity::Max ? maxSlots_ : minSlots_;
        }

        /**
         * @brief Checks whether a node identifier addresses an existing slot.
         * @param side Slot array to inspect.
         * @param nodeId Candidate node identifier.
         * @return `true` when `nodeId` lies inside `side`.
         */
        [[nodiscard]] static bool validSlot(const std::vector<AgendaSlot>& side, NodeId nodeId) noexcept {
            return nodeId >= 0 && nodeId < static_cast<NodeId>(side.size());
        }

        /**
         * @brief Removes an active candidate while preserving the slot/set invariant.
         * @param polarity Component-tree polarity.
         * @param nodeId Candidate node identifier.
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

      public:
        /**
         * @brief Creates an empty agenda.
         * @param policy Deterministic equal-area ordering policy.
         * @param excludeExteriorSeed Whether exterior-containing candidates must be excluded.
         */
        CandidateAgenda(SdrtTiePolicy policy, bool excludeExteriorSeed)
            : policy_(policy), excludeExteriorSeed_(excludeExteriorSeed), candidates_(CandidateLess{policy}) {}

        /**
         * @brief Resets ordered entries and sizes the per-polarity node slots.
         * @param maxNodeSlots Number of max-tree node slots.
         * @param minNodeSlots Number of min-tree node slots.
         */
        void reset(std::size_t maxNodeSlots, std::size_t minNodeSlots) {
            candidates_ = std::set<Candidate, CandidateLess>(CandidateLess{policy_});
            maxSlots_.assign(maxNodeSlots, AgendaSlot{});
            minSlots_.assign(minNodeSlots, AgendaSlot{});
        }

        /**
         * @brief Inserts or replaces one current candidate using partition metadata.
         * @param polarity Component-tree polarity.
         * @param nodeId Current component-tree node identifier.
         * @param descriptor Metadata supplied by the current flat-zone partition.
         */
        void upsert(Polarity polarity, NodeId nodeId, const CandidateDescriptor& descriptor) {
            auto& side = slots(polarity);
            if (!validSlot(side, nodeId)) {
                throw std::out_of_range("Min/max residual candidate lies outside its agenda slot domain.");
            }
            remove(polarity, nodeId);
            if (descriptor.area <= 0 || descriptor.stableSpatialKey == InvalidNode) {
                throw std::runtime_error("Min/max residual agenda received invalid flat-zone metadata.");
            }
            if (excludeExteriorSeed_ && descriptor.containsExteriorSeed) {
                return;
            }

            const Candidate candidate{descriptor.area, descriptor.stableSpatialKey, polarity, nodeId};
            const auto [_, inserted] = candidates_.insert(candidate);
            if (!inserted) {
                throw std::runtime_error("Min/max residual agenda received a duplicate candidate.");
            }
            auto& slot = side[static_cast<std::size_t>(nodeId)];
            slot.active = true;
            slot.candidate = candidate;
        }

        /**
         * @brief Removes a node from the active agenda, if present.
         * @param polarity Component-tree polarity.
         * @param nodeId Candidate node identifier.
         */
        void reject(Polarity polarity, NodeId nodeId) { remove(polarity, nodeId); }

        /** @return First candidate in deterministic order, or `std::nullopt` when empty. */
        [[nodiscard]] std::optional<Candidate> select() const {
            if (candidates_.empty()) {
                return std::nullopt;
            }
            return *candidates_.begin();
        }

        /**
         * @brief Checks whether a node currently owns an active agenda entry.
         * @param polarity Component-tree polarity.
         * @param nodeId Candidate node identifier.
         * @return `true` when the corresponding slot is active.
         */
        [[nodiscard]] bool contains(Polarity polarity, NodeId nodeId) const {
            const auto& side = slots(polarity);
            return validSlot(side, nodeId) && side[static_cast<std::size_t>(nodeId)].active;
        }
    };
};

} // namespace mmcfilters::sdrt::detail
