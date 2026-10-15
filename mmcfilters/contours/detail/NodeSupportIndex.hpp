#pragma once

#include "../../trees/MorphologicalTree.hpp"
#include "../../trees/detail/CommittedTreeAccess.hpp"
#include "../../utils/Common.hpp"

#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace mmcfilters::contours::detail {

/**
 * @brief Half-open interval in the shared pre-order support storage.
 */
struct SupportInterval {
    /// Inclusive offset of the first support pixel.
    std::size_t begin = 0;
    /// Exclusive offset past the last support pixel.
    std::size_t end = 0;

    /**
     * @brief Returns the support cardinality represented by the interval.
     * @return Number of pixels in the interval.
     */
    [[nodiscard]] std::size_t size() const noexcept { return end - begin; }
};

/**
 * @brief Dense immutable support index for one committed morphological tree.
 *
 * Proper parts are appended once in pre-order. Consequently, all proper parts
 * in a node subtree occupy one contiguous interval. Support cardinalities are
 * accumulated bottom-up to establish interval ends. Dead node slots retain
 * empty intervals.
 *
 * The referenced tree must outlive the index. Returned spans borrow the
 * index-owned pixel storage and remain valid only while the index exists and
 * the tree stays at the captured mutation version.
 */
class NodeSupportIndex {
  public:
    NodeSupportIndex(const NodeSupportIndex&) = delete;
    NodeSupportIndex& operator=(const NodeSupportIndex&) = delete;
    NodeSupportIndex(NodeSupportIndex&&) = delete;
    NodeSupportIndex& operator=(NodeSupportIndex&&) = delete;

    /**
     * @brief Builds the O(P + N) support index for a committed 2D tree.
     * @param tree Stable source tree with a regular 2D domain.
     */
    explicit NodeSupportIndex(const MorphologicalTree& tree)
        : tree_(tree), mutationVersion_(tree.getMutationVersion()), domain_(requireDomain(tree)),
          supportIntervalByNode_(static_cast<std::size_t>(tree.numInternalNodeSlots())) {
        supportPixels_.reserve(static_cast<std::size_t>(tree_.numPixels()));
        std::vector<NodeId> nodesInPreorder;
        nodesInPreorder.reserve(static_cast<std::size_t>(tree_.numNodes()));

        for (NodeId node : ::mmcfilters::detail::CommittedTreeAccess::subtree(tree_, tree_.root())) {
            const std::size_t nodeIndex = static_cast<std::size_t>(node);
            nodesInPreorder.push_back(node);
            supportIntervalByNode_[nodeIndex].begin = supportPixels_.size();

            for (PixelId pixel : ::mmcfilters::detail::CommittedTreeAccess::properParts(tree_, node)) {
                supportPixels_.push_back(pixel);
                ++supportIntervalByNode_[nodeIndex].end; // Initially the proper-part cardinality.
            }
        }

        for (auto iterator = nodesInPreorder.rbegin(); iterator != nodesInPreorder.rend(); ++iterator) {
            const NodeId node = *iterator;
            const std::size_t nodeIndex = static_cast<std::size_t>(node);
            const std::size_t supportCardinality = supportIntervalByNode_[nodeIndex].end;
            supportIntervalByNode_[nodeIndex].end = supportIntervalByNode_[nodeIndex].begin + supportCardinality;
            if (node != tree_.root()) {
                const NodeId parent = ::mmcfilters::detail::CommittedTreeAccess::nodeParent(tree_, node);
                const std::size_t parentIndex = static_cast<std::size_t>(parent);
                supportIntervalByNode_[parentIndex].end += supportCardinality;
            }
        }
    }

    /**
     * @brief Returns the contiguous support pixels of one live node.
     * @param node Live node identifier.
     * @return Borrowed span containing the node support.
     */
    [[nodiscard]] std::span<const PixelId> support(NodeId node) const {
        const SupportInterval interval = supportInterval(node);
        return std::span<const PixelId>(supportPixels_).subspan(interval.begin, interval.size());
    }

    /**
     * @brief Returns the half-open support interval of one live node.
     * @param node Live node identifier.
     * @return Half-open interval in the shared support buffer.
     */
    [[nodiscard]] SupportInterval supportInterval(NodeId node) const {
        requireLiveNode(node);
        return supportIntervalByNode_[static_cast<std::size_t>(node)];
    }

    /**
     * @brief Returns support for a caller-established live node.
     * @param node Node already established as live and compatible.
     * @return Borrowed span containing the node support.
     */
    [[nodiscard]] std::span<const PixelId> establishedSupport(NodeId node) const noexcept {
        const SupportInterval interval = establishedSupportInterval(node);
        return std::span<const PixelId>(supportPixels_).subspan(interval.begin, interval.size());
    }

    /**
     * @brief Returns the support interval of a caller-established live node.
     * @param node Node already established as live and compatible.
     * @return Half-open interval in the shared support buffer.
     */
    [[nodiscard]] SupportInterval establishedSupportInterval(NodeId node) const noexcept {
        return supportIntervalByNode_[static_cast<std::size_t>(node)];
    }

    /**
     * @brief Returns the number of proper-part pixels stored exactly once.
     * @return Number of pixels in the shared support buffer.
     */
    [[nodiscard]] std::size_t numIndexedPixels() const {
        tree_.requireMutationVersion(mutationVersion_, "NodeSupportIndex");
        return supportPixels_.size();
    }

    /**
     * @brief Rejects use with a different tree, mutation version, or domain.
     * @param tree Candidate tree to compare with the indexed tree.
     */
    void requireCompatibleTree(const MorphologicalTree& tree) const {
        if (&tree != &tree_) {
            throw std::invalid_argument("Node support index belongs to a different tree instance.");
        }
        if (tree_.getMutationVersion() != mutationVersion_) {
            throw std::logic_error("Node support index belongs to an older tree mutation version.");
        }
        const GridDomain2D domain = tree_.requireGridDomain2D("NodeSupportIndex");
        if (domain.rows != domain_.rows || domain.columns != domain_.columns) {
            throw std::logic_error("Node support index no longer matches the captured 2D domain.");
        }
    }

  private:
    /**
     * @brief Validates and returns the regular 2D domain of a tree.
     * @param tree Tree whose domain is required.
     * @return Valid non-empty 2D grid domain.
     */
    [[nodiscard]] static GridDomain2D requireDomain(const MorphologicalTree& tree) {
        tree.requireNotEditing("NodeSupportIndex");
        if (!tree.hasGridDomain2D()) {
            throw std::invalid_argument("Node support indexing requires a regular 2D pixel domain.");
        }
        const GridDomain2D domain = *tree.gridDomain2D();
        if (domain.rows <= 0 || domain.columns <= 0 || domain.columns > std::numeric_limits<int>::max() / domain.rows ||
            domain.rows * domain.columns != tree.numPixels()) {
            throw std::invalid_argument("Node support indexing requires a consistent non-empty 2D pixel domain.");
        }
        return domain;
    }

    /**
     * @brief Rejects a node that is not live in the indexed tree.
     * @param node Candidate node identifier.
     */
    void requireLiveNode(NodeId node) const {
        tree_.requireMutationVersion(mutationVersion_, "NodeSupportIndex");
        if (!tree_.isNode(node) || !::mmcfilters::detail::CommittedTreeAccess::isAlive(tree_, node)) {
            throw std::out_of_range("Node support index received a non-live node id.");
        }
    }

    const MorphologicalTree& tree_;                 ///< Indexed source tree.
    std::size_t mutationVersion_ = 0;               ///< Mutation version captured at construction.
    GridDomain2D domain_{};                         ///< Validated image domain.
    std::vector<PixelId> supportPixels_;            ///< Proper-part pixels in tree pre-order.
    std::vector<SupportInterval> supportIntervalByNode_; ///< Support slice for each internal node slot.
};

} // namespace mmcfilters::contours::detail
