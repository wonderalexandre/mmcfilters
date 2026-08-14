#pragma once

/**
 * @file ResidualTreeCandidateAgenda.hpp
 * @brief Deterministic candidate ordering shared by residual-tree modes.
 */

#include "ResidualTreeCandidateTypes.hpp"
#include "../ResidualEvolution.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <set>
#include <stdexcept>
#include <vector>

namespace mmcfilters::sdrt::detail {

/** @brief Candidate, ordering, and agenda types shared by all residual altitude types. */
struct ResidualTreeAgendaTypes {
    /** @brief Public polarity vocabulary used by the synchronized evolution. */
    using Polarity = sdrt::Polarity;

    /** @brief Candidate entry stored in deterministic agenda order. */
    struct Candidate {
        SelfDualResidualKey residualKey;                 ///< Canonical support-only scheduling key.
        Polarity polarity = Polarity::Maximum;           ///< Component-tree polarity, excluded from ordering.
        NodeId nodeId = InvalidNode;                     ///< Current component-tree node, excluded from ordering.
    };

    /** @brief Neutral candidate metadata supplied by the flat-zone partition. */
    using CandidateDescriptor = ResidualTreeCandidateDescriptor;

    /** @brief Orders candidates only by the canonical self-dual residual key. */
    struct CandidateLess {
        SelfDualResidualOrder residualOrder; ///< Canonical contrast-invariant key order.

        /**
         * @brief Compares two candidates without consulting polarity, altitude, or node id.
         * @param lhs Left candidate.
         * @param rhs Right candidate.
         * @return `true` when `lhs` precedes `rhs` in agenda order.
         */
        [[nodiscard]] bool operator()(const Candidate& lhs, const Candidate& rhs) const {
            return residualOrder.compareResidualCandidates(lhs.residualKey, rhs.residualKey);
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
        SelfDualResidualSchedule schedule_;              ///< Unique canonical residual schedule.
        bool excludeInfinityPixel_ = true;               ///< Whether infinity-pixel candidates are excluded.
        std::set<Candidate, CandidateLess> candidates_; ///< Ordered active candidates.
        std::vector<AgendaSlot> maxSlots_;               ///< Candidate slots for the max-tree.
        std::vector<AgendaSlot> minSlots_;               ///< Candidate slots for the min-tree.

        /**
         * @brief Selects the mutable slot array for one polarity.
         * @param polarity Component-tree polarity.
         * @return Mutable slot array associated with `polarity`.
         */
        [[nodiscard]] std::vector<AgendaSlot>& slots(Polarity polarity) {
            return polarity == Polarity::Maximum ? maxSlots_ : minSlots_;
        }

        /**
         * @brief Selects the immutable slot array for one polarity.
         * @param polarity Component-tree polarity.
         * @return Immutable slot array associated with `polarity`.
         */
        [[nodiscard]] const std::vector<AgendaSlot>& slots(Polarity polarity) const {
            return polarity == Polarity::Maximum ? maxSlots_ : minSlots_;
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
         * @param spatialOrder Total order defining the spatial-minimum key coordinate.
         * @param excludeInfinityPixel Whether candidates containing the infinity pixel must be excluded.
         */
        CandidateAgenda(SpatialOrder spatialOrder, bool excludeInfinityPixel)
            : schedule_(std::move(spatialOrder)), excludeInfinityPixel_(excludeInfinityPixel),
              candidates_(CandidateLess{schedule_.residualOrder()}) {}

        /**
         * @brief Resets ordered entries and sizes the per-polarity node slots.
         * @param maxNodeSlots Number of max-tree node slots.
         * @param minNodeSlots Number of min-tree node slots.
         */
        void reset(std::size_t maxNodeSlots, std::size_t minNodeSlots) {
            candidates_ = std::set<Candidate, CandidateLess>(CandidateLess{schedule_.residualOrder()});
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
            if (descriptor.supportCardinality == 0 || descriptor.spatialMinimum == InvalidPixel) {
                throw std::runtime_error("Min/max residual agenda received invalid flat-zone metadata.");
            }
            if (excludeInfinityPixel_ && descriptor.containsInfinityPixel) {
                return;
            }

            const Candidate candidate{SelfDualResidualKey{descriptor.supportCardinality, descriptor.spatialMinimum}, polarity, nodeId};
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
