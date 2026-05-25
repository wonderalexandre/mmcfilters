#pragma once

#include "../../WeightedMorphologicalTree.hpp"
#include "../../../utils/Common.hpp"

#include <algorithm>
#include <optional>
#include <set>
#include <stdexcept>
#include <vector>

namespace mmcfilters::sdrt::detail {

/**
 * @brief Identifies which component-tree side owns an SDRT pruning candidate.
 */
enum class ResidualPolarity {
    Max,
    Min
};

/**
 * @brief Snapshot of one currently prunable leaf candidate.
 *
 * `area` is copied from `MorphologicalTree::getNumProperParts(nodeId)` at the
 * moment the candidate is inserted into the agenda. The builder revalidates the
 * snapshot before pruning because dual-tree edits can make an old candidate
 * stale before it reaches the front of the priority queue.
 */
struct CandidateLeaf {
    ResidualPolarity polarity = ResidualPolarity::Max;
    NodeId nodeId = InvalidNode;
    int area = 0;
};

/**
 * Incremental priority agenda used by the SDRT builder to choose the next
 * component-tree leaf to prune.
 *
 * The self-dual residual-tree construction evolves two mutable component trees
 * at the same time: a max-tree and a min-tree. At every iteration the builder
 * must pick the smallest currently prunable leaf from either side. Equal-area
 * candidates are resolved deterministically by selecting the max-tree side
 * first. This class keeps that frontier explicitly instead of rebuilding it
 * from the whole pair of trees after each pruning operation.
 *
 * A node is present in the agenda only when it is:
 * - inside the node-id domain reserved for that tree side;
 * - alive in the current topology;
 * - not the root;
 * - a leaf;
 * - backed by at least one proper part.
 *
 * Each side stores the same information in two complementary forms. The slot
 * vector is indexed by NodeId and supports O(1) membership/state checks. The
 * ordered set stores the active candidates sorted by (area, node id), so the
 * next leaf can be read from begin(). The fundamental invariant is:
 *
 *     side.slots[n].active == true
 *         iff side.candidates contains CandidateKey{side.slots[n].area, n}
 *
 * update() is the only mutating operation that should be used after
 * initialization. It first removes any previous version of the node, then
 * re-inserts it only if the node is still a valid leaf candidate in the
 * supplied tree. This makes affected-node updates cheap and also handles nodes
 * that stopped being candidates after a dual-tree synchronization step.
 */
template<AltitudeValue T>
class LeafPrunePriorityAgenda {
private:
    using tree_t = WeightedMorphologicalTree<T>;

    /**
     * @brief Ordering key for the set-backed priority queue.
     *
     * Node id is part of the key to make equal-area candidates deterministic
     * and uniquely addressable within one polarity.
     */
    struct CandidateKey {
        int area = 0;
        NodeId nodeId = InvalidNode;
    };

    /**
     * @brief Sorts candidates by increasing area, then increasing node id.
     */
    struct CandidateKeyLess {
        bool operator()(const CandidateKey& lhs, const CandidateKey& rhs) const {
            if (lhs.area != rhs.area) {
                return lhs.area < rhs.area;
            }
            return lhs.nodeId < rhs.nodeId;
        }
    };

    /**
     * @brief O(1) side-table entry for a node currently tracked by the agenda.
     */
    struct CandidateSlot {
        bool active = false;
        int area = 0;
    };

    /**
     * @brief Candidate state for one evolving component tree.
     *
     * `slots` is indexed by internal `NodeId`; `candidates` is the ordered
     * priority queue. `maxArea` is kept as a defensive upper bound for area
     * values that should be possible in the tree domain.
     */
    struct AgendaSide {
        std::vector<CandidateSlot> slots;
        std::set<CandidateKey, CandidateKeyLess> candidates;
        int maxArea = 0;

        /**
         * @brief Clears a side and sizes it to the current internal node domain.
         */
        void reset(std::size_t nodeCount, int totalArea) {
            slots.assign(nodeCount, CandidateSlot{});
            candidates.clear();
            maxArea = std::max(0, totalArea);
        }
    };

    AgendaSide maxSide_;
    AgendaSide minSide_;

    /**
     * @brief Returns whether a node id can index the side's slot table.
     *
     * Nodes outside this range can appear in affected-node lists after edits
     * and are intentionally ignored by agenda synchronization.
     */
    static bool isValidSlot(const AgendaSide& side, NodeId nodeId) {
        return nodeId != InvalidNode &&
               nodeId >= 0 &&
               nodeId < static_cast<NodeId>(side.slots.size());
    }

    /**
     * @brief Checks that an area can belong to the represented image domain.
     */
    static bool isValidArea(const AgendaSide& side, int area) {
        return area > 0 && area <= side.maxArea;
    }

    /**
     * @brief Removes the current version of a node from one side.
     *
     * The method preserves the slot/set invariant. Invariant failures are hard
     * errors because they mean a previous update left the agenda internally
     * inconsistent.
     */
    static void removeNode(AgendaSide& side, NodeId nodeId) {
        if (!isValidSlot(side, nodeId)) {
            return;
        }
        auto& slot = side.slots[static_cast<std::size_t>(nodeId)];
        if (!slot.active) {
            return;
        }
        if (!isValidArea(side, slot.area)) {
            throw std::runtime_error("SDRT candidate agenda found an invalid active area.");
        }
        const auto erased = side.candidates.erase(CandidateKey{slot.area, nodeId});
        if (erased != 1) {
            throw std::runtime_error("SDRT candidate agenda violates the slot/set invariant.");
        }
        slot = CandidateSlot{};
    }

    /**
     * @brief Reads the best candidate in one side without mutating the agenda.
     *
     * The candidate is revalidated against the slot vector before being
     * exposed, so any mismatch between the two internal representations is
     * detected close to the point where it would affect pruning order.
     */
    static std::optional<CandidateLeaf> topInSide(const AgendaSide& side, ResidualPolarity polarity) {
        if (side.candidates.empty()) {
            return std::nullopt;
        }
        const CandidateKey& key = *side.candidates.begin();
        if (!isValidSlot(side, key.nodeId)) {
            throw std::runtime_error("SDRT candidate agenda top is outside the slot domain.");
        }
        const auto& slot = side.slots[static_cast<std::size_t>(key.nodeId)];
        if (!slot.active || slot.area != key.area) {
            throw std::runtime_error("SDRT candidate agenda top violates the slot/set invariant.");
        }
        return CandidateLeaf{polarity, key.nodeId, key.area};
    }

    /**
     * @brief Indexes all currently alive candidate leaves from one tree side.
     */
    static void indexTreeSide(const tree_t& tree,
                              ResidualPolarity polarity,
                              LeafPrunePriorityAgenda& agenda) {
        const MorphologicalTree& topology = tree.topology();
        for (NodeId nodeId : topology.getAliveNodeIds()) {
            // update() performs all candidate checks, so initialization and
            // incremental maintenance share the exact same admission rules.
            agenda.update(tree, polarity, nodeId);
        }
    }

public:
    /**
     * @brief Builds both agenda sides from the initial max/min tree state.
     *
     * This is intended for the beginning of SDRT construction, before any
     * pruning has happened.
     */
    void initialize(const tree_t& maxTree,
                    const tree_t& minTree) {
        maxSide_.reset(static_cast<std::size_t>(maxTree.topology().getNumInternalNodeSlots()),
                       maxTree.topology().getNumTotalProperParts());
        minSide_.reset(static_cast<std::size_t>(minTree.topology().getNumInternalNodeSlots()),
                       minTree.topology().getNumTotalProperParts());
        indexTreeSide(maxTree, ResidualPolarity::Max, *this);
        indexTreeSide(minTree, ResidualPolarity::Min, *this);
    }

    /**
     * @brief Reports whether a node is currently active in one agenda side.
     *
     * The builder uses this to decide whether an affected dual-tree node needs
     * to be revalidated even when it may no longer be a valid leaf after the
     * edit.
     */
    bool contains(ResidualPolarity polarity, NodeId nodeId) const {
        const AgendaSide& side = polarity == ResidualPolarity::Max ? maxSide_ : minSide_;
        return isValidSlot(side, nodeId) && side.slots[static_cast<std::size_t>(nodeId)].active;
    }

    /**
     * @brief Synchronizes one node with the supplied tree side.
     *
     * This method is safe to call for nodes that are no longer alive, no longer
     * leaves, or outside the side's slot domain. Valid candidates are inserted
     * with their current area; invalid candidates are simply absent after the
     * call.
     */
    void update(const tree_t& tree, ResidualPolarity polarity, NodeId nodeId) {
        AgendaSide& side = polarity == ResidualPolarity::Max ? maxSide_ : minSide_;
        if (!isValidSlot(side, nodeId)) {
            return;
        }
        removeNode(side, nodeId);

        const MorphologicalTree& topology = tree.topology();
        if (nodeId == topology.getRoot() ||
            !topology.isNode(nodeId) ||
            !topology.isAlive(nodeId) ||
            !topology.isLeaf(nodeId)) {
            return;
        }

        const NodeId parentId = topology.getNodeParent(nodeId);
        if (parentId == InvalidNode || parentId == nodeId) {
            throw std::runtime_error("SDRT candidate leaf violates the parent invariant.");
        }

        const int area = topology.getNumProperParts(nodeId);
        if (area <= 0) {
            return;
        }
        if (area > side.maxArea) {
            throw std::runtime_error("SDRT candidate agenda received an invalid candidate area.");
        }

        auto& slot = side.slots[static_cast<std::size_t>(nodeId)];
        const auto [_, inserted] = side.candidates.insert(CandidateKey{area, nodeId});
        if (!inserted) {
            throw std::runtime_error("SDRT candidate agenda received a duplicate candidate.");
        }
        slot.active = true;
        slot.area = area;
    }

    /**
     * @brief Returns the globally smallest leaf across max and min sides.
     *
     * Equal-area leaves are resolved by selecting the max-tree side first.
     */
    std::optional<CandidateLeaf> select() const {
        auto maxCandidate = topInSide(maxSide_, ResidualPolarity::Max);
        auto minCandidate = topInSide(minSide_, ResidualPolarity::Min);
        if (!maxCandidate.has_value()) {
            return minCandidate;
        }
        if (!minCandidate.has_value()) {
            return maxCandidate;
        }

        if (maxCandidate->area == minCandidate->area) {
            return maxCandidate;
        }
        return maxCandidate->area < minCandidate->area ? maxCandidate : minCandidate;
    }
};


} // namespace mmcfilters::sdrt::detail
