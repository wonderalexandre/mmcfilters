#pragma once

#include "../../../../trees/MorphologicalTree.hpp"
#include "../../../../trees/detail/CommittedTreeAccess.hpp"

#include <algorithm>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

namespace mmcfilters::attributes::computers::detail::distance_transform_approx {

/**
 * @brief Groups live nodes into topology-derived bottom-up level sets.
 *
 * Leaves occupy level zero and every other node occupies one level above its
 * highest child. Nodes in one level are therefore pairwise incomparable and
 * may be updated as one independent DIFT batch without using altitude or tree
 * kind.
 */
class MorphologicalTreeTopologicalLevelIndex {
  public:
    MorphologicalTreeTopologicalLevelIndex(const MorphologicalTreeTopologicalLevelIndex&) = delete;
    MorphologicalTreeTopologicalLevelIndex& operator=(const MorphologicalTreeTopologicalLevelIndex&) = delete;
    MorphologicalTreeTopologicalLevelIndex(MorphologicalTreeTopologicalLevelIndex&&) = delete;
    MorphologicalTreeTopologicalLevelIndex& operator=(MorphologicalTreeTopologicalLevelIndex&&) = delete;

    explicit MorphologicalTreeTopologicalLevelIndex(const MorphologicalTree& tree)
        : tree_(tree), mutationVersion_(tree.getMutationVersion()), levelByNode_(static_cast<std::size_t>(tree.numInternalNodeSlots()), -1) {
        tree_.requireNotEditing("MorphologicalTreeTopologicalLevelIndex");
        int maximumLevel = -1;
        for (NodeId node : tree_.postOrder()) {
            int nodeLevel = 0;
            for (NodeId child : ::mmcfilters::detail::CommittedTreeAccess::children(tree_, node)) {
                const int childLevel = levelByNode_[static_cast<std::size_t>(child)];
                nodeLevel = std::max(nodeLevel, childLevel + 1);
            }
            levelByNode_[static_cast<std::size_t>(node)] = nodeLevel;
            maximumLevel = std::max(maximumLevel, nodeLevel);
        }
        levelOffsets_.assign(static_cast<std::size_t>(maximumLevel) + 2, 0);
        for (NodeId node : tree_.postOrder()) {
            ++levelOffsets_[static_cast<std::size_t>(levelByNode_[static_cast<std::size_t>(node)]) + 1];
        }
        for (std::size_t index = 1; index < levelOffsets_.size(); ++index) {
            levelOffsets_[index] += levelOffsets_[index - 1];
        }
        nodes_.assign(static_cast<std::size_t>(tree_.numNodes()), InvalidNode);
        std::vector<std::size_t> cursors = levelOffsets_;
        for (NodeId node : tree_.postOrder()) {
            nodes_[cursors[static_cast<std::size_t>(levelByNode_[static_cast<std::size_t>(node)])]++] = node;
        }
    }

    [[nodiscard]] int numLevels() const {
        requireStableTree();
        return static_cast<int>(levelOffsets_.size()) - 1;
    }

    [[nodiscard]] int level(NodeId node) const {
        requireLiveNode(node);
        return levelByNode_[static_cast<std::size_t>(node)];
    }

    [[nodiscard]] std::span<const NodeId> nodesAtLevel(int levelValue) const {
        requireStableTree();
        if (levelValue < 0 || levelValue >= static_cast<int>(levelOffsets_.size()) - 1) {
            throw std::out_of_range("Topological level index received an invalid level.");
        }
        const std::size_t levelIndex = static_cast<std::size_t>(levelValue);
        return std::span<const NodeId>(nodes_).subspan(levelOffsets_[levelIndex], levelOffsets_[levelIndex + 1] - levelOffsets_[levelIndex]);
    }

    /**
     * @brief Returns the number of levels after the caller established tree stability.
     */
    [[nodiscard]] int establishedNumLevels() const noexcept { return static_cast<int>(levelOffsets_.size()) - 1; }

    /**
     * @brief Returns nodes at a caller-established topological level.
     */
    [[nodiscard]] std::span<const NodeId> establishedNodesAtLevel(int levelValue) const noexcept {
        const std::size_t levelIndex = static_cast<std::size_t>(levelValue);
        return std::span<const NodeId>(nodes_).subspan(levelOffsets_[levelIndex], levelOffsets_[levelIndex + 1] - levelOffsets_[levelIndex]);
    }

  private:
    void requireStableTree() const { tree_.requireMutationVersion(mutationVersion_, "MorphologicalTreeTopologicalLevelIndex"); }

    void requireLiveNode(NodeId node) const {
        requireStableTree();
        if (!tree_.isAlive(node)) {
            throw std::out_of_range("Topological level index received a non-live node id.");
        }
    }

    const MorphologicalTree& tree_;
    std::size_t mutationVersion_ = 0;
    std::vector<int> levelByNode_;
    std::vector<std::size_t> levelOffsets_;
    std::vector<NodeId> nodes_;
};

} // namespace mmcfilters::attributes::computers::detail::distance_transform_approx
