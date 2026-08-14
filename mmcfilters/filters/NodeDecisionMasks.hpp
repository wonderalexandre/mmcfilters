#pragma once

#include <cstddef>
#include <utility>
#include <vector>

namespace mmcfilters {

/**
 * @brief Dense Boolean decisions where true means preserve the node.
 *
 * The mask is indexed by the dense internal `NodeId` domain. It is deliberately
 * not implicitly interchangeable with `NodePruningMask`.
 */
class NodePreservationMask {
  private:
    std::vector<bool> decisions_; ///< Dense decisions indexed by the internal node-slot domain.

  public:
    /** @brief Creates an empty preservation mask. */
    NodePreservationMask() = default;

    /**
     * @brief Creates a preservation mask from explicit decisions.
     * @param decisions Dense decisions where `true` preserves the corresponding node.
     */
    explicit NodePreservationMask(std::vector<bool> decisions) : decisions_(std::move(decisions)) {}

    /**
     * @brief Creates a preservation mask whose entries share one decision.
     * @param size Number of dense node slots represented by the mask.
     * @param preserveNode Decision assigned to every slot.
     */
    NodePreservationMask(std::size_t size, bool preserveNode) : decisions_(size, preserveNode) {}

    /** @brief Returns the number of dense node slots. @return Number of decisions. */
    [[nodiscard]] std::size_t size() const noexcept { return decisions_.size(); }

    /** @brief Tests whether the mask contains no decisions. @return True when the mask is empty. */
    [[nodiscard]] bool empty() const noexcept { return decisions_.empty(); }

    /**
     * @brief Returns one preservation decision.
     * @param index Dense node-slot index.
     * @return True when the indexed node must be preserved.
     */
    [[nodiscard]] bool operator[](std::size_t index) const { return decisions_[index]; }

    /** @brief Returns all preservation decisions. @return Read-only dense decision buffer. */
    [[nodiscard]] const std::vector<bool>& decisions() const noexcept { return decisions_; }
};

/**
 * @brief Dense Boolean decisions where true means prune the node.
 *
 * The polarity is the Boolean complement of `NodePreservationMask`, but the two
 * types require an explicit conversion.
 */
class NodePruningMask {
  private:
    std::vector<bool> decisions_; ///< Dense decisions indexed by the internal node-slot domain.

  public:
    /** @brief Creates an empty pruning mask. */
    NodePruningMask() = default;

    /**
     * @brief Creates a pruning mask from explicit decisions.
     * @param decisions Dense decisions where `true` prunes the corresponding node.
     */
    explicit NodePruningMask(std::vector<bool> decisions) : decisions_(std::move(decisions)) {}

    /**
     * @brief Creates a pruning mask whose entries share one decision.
     * @param size Number of dense node slots represented by the mask.
     * @param pruneNode Decision assigned to every slot.
     */
    NodePruningMask(std::size_t size, bool pruneNode) : decisions_(size, pruneNode) {}

    /** @brief Returns the number of dense node slots. @return Number of decisions. */
    [[nodiscard]] std::size_t size() const noexcept { return decisions_.size(); }

    /** @brief Tests whether the mask contains no decisions. @return True when the mask is empty. */
    [[nodiscard]] bool empty() const noexcept { return decisions_.empty(); }

    /**
     * @brief Returns one pruning decision.
     * @param index Dense node-slot index.
     * @return True when the indexed node must be pruned.
     */
    [[nodiscard]] bool operator[](std::size_t index) const { return decisions_[index]; }

    /** @brief Returns all pruning decisions. @return Read-only dense decision buffer. */
    [[nodiscard]] const std::vector<bool>& decisions() const noexcept { return decisions_; }
};

/** @brief Explicitly complements preservation decisions into pruning decisions. */
[[nodiscard]] inline NodePruningMask toNodePruningMask(const NodePreservationMask& nodePreservationMask) {
    std::vector<bool> nodePruningDecisions(nodePreservationMask.size(), false);
    for (std::size_t index = 0; index < nodePreservationMask.size(); ++index) {
        nodePruningDecisions[index] = !nodePreservationMask[index];
    }
    return NodePruningMask(std::move(nodePruningDecisions));
}

/** @brief Explicitly complements pruning decisions into preservation decisions. */
[[nodiscard]] inline NodePreservationMask toNodePreservationMask(const NodePruningMask& nodePruningMask) {
    std::vector<bool> nodePreservationDecisions(nodePruningMask.size(), false);
    for (std::size_t index = 0; index < nodePruningMask.size(); ++index) {
        nodePreservationDecisions[index] = !nodePruningMask[index];
    }
    return NodePreservationMask(std::move(nodePreservationDecisions));
}

} // namespace mmcfilters
