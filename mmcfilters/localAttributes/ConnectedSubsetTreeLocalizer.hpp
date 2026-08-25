#pragma once

#include "../trees/MorphologicalTree.hpp"
#include "../trees/detail/CommittedTreeAccess.hpp"
#include "../utils/Common.hpp"

#include <cstddef>
#include <optional>
#include <stdexcept>

namespace mmcfilters::local_attributes {

/** @brief Relative row-column offset of one translated local sample. */
struct WindowOffset {
    int rowOffset = 0;    ///< Positive values move down.
    int columnOffset = 0; ///< Positive values move right.

    /** @brief Compares both coordinate displacements. */
    friend bool operator==(const WindowOffset& lhs, const WindowOffset& rhs) = default;
};

namespace detail::kernel {

/**
 * @brief Returns the inclusion join of two smallest nodes in a connected-subset tree.
 *
 * In a rooted laminar hierarchy the join is the lowest common ancestor.  This
 * single primitive covers equal, comparable, and cross-branch smallest-node pairs.
 */
inline NodeId connectedSubsetJoin(const MorphologicalTree& tree, NodeId lhsSmallestNode, NodeId rhsSmallestNode) {
    return ::mmcfilters::detail::CommittedTreeAccess::lowestCommonAncestor(tree, lhsSmallestNode, rhsSmallestNode);
}

/** @brief Compatibility spelling for the anchored-entry join of two smallest nodes. */
inline NodeId anchoredEntryFromSmallestNodes(const MorphologicalTree& tree, NodeId anchorSmallestNode, NodeId sampleSmallestNode) {
    return connectedSubsetJoin(tree, anchorSmallestNode, sampleSmallestNode);
}

/** @brief Locates an absolute sample on the branch of its anchor. */
inline NodeId anchoredEntry(const MorphologicalTree& tree, PixelId anchorPixel, PixelId samplePixel) {
    const NodeId anchorSmallestNode = ::mmcfilters::detail::CommittedTreeAccess::smallestNodeMap(tree, anchorPixel);
    const NodeId sampleSmallestNode = ::mmcfilters::detail::CommittedTreeAccess::smallestNodeMap(tree, samplePixel);
    return connectedSubsetJoin(tree, anchorSmallestNode, sampleSmallestNode);
}

/** @brief Locates a translated sample, returning `InvalidNode` outside the grid. */
inline NodeId anchoredEntry(const MorphologicalTree& tree, PixelId anchorPixel, WindowOffset offset) {
    const GridDomain2D& domain = ::mmcfilters::detail::CommittedTreeAccess::gridDomain2D(tree);
    const int anchorRow = anchorPixel / domain.columns;
    const int anchorColumn = anchorPixel % domain.columns;
    const int sampleRow = anchorRow + offset.rowOffset;
    const int sampleColumn = anchorColumn + offset.columnOffset;
    if (sampleRow < 0 || sampleRow >= domain.rows || sampleColumn < 0 || sampleColumn >= domain.columns) {
        return InvalidNode;
    }
    return anchoredEntry(tree, anchorPixel, static_cast<PixelId>(sampleRow * domain.columns + sampleColumn));
}

} // namespace detail::kernel

/**
 * @brief Validated public view of connected-subset-tree localization.
 *
 * The hot finite-window kernels use the equivalent functions in
 * `detail::kernel` after validating once at their public boundary.  This view
 * is intended for explicit theorem-facing queries and therefore validates
 * pixels, nodes, the grid domain, and captured-tree stability per call.
 */
class ConnectedSubsetTreeLocalizer {
  public:
    /** @brief Captures one established tree without building another index. @param tree Stable connected-subset tree. */
    explicit ConnectedSubsetTreeLocalizer(const MorphologicalTree& tree) : tree_(tree), mutationVersion_(tree.getMutationVersion()) {
        tree_.requireNotEditing("ConnectedSubsetTreeLocalizer");
    }

    /** @brief Returns the inclusion-smallest node of one valid pixel. @param pixel Valid pixel. @return Pixel's smallest node. */
    [[nodiscard]] NodeId smallestNode(PixelId pixel) const {
        requireStableTree();
        if (!tree_.isPixel(pixel)) {
            throw std::out_of_range("ConnectedSubsetTreeLocalizer::smallestNode received an invalid pixel.");
        }
        const NodeId value = tree_.smallestNode(pixel);
        if (!tree_.isAlive(value)) {
            throw std::logic_error("ConnectedSubsetTreeLocalizer::smallestNode found a pixel without a live smallest node.");
        }
        return value;
    }

    /**
     * @brief Returns the inclusion join, equivalently the LCA, of two live smallest nodes.
     * @param lhsSmallestNode First live node.
     * @param rhsSmallestNode Second live node.
     * @return Inclusion-smallest common ancestor.
     */
    [[nodiscard]] NodeId join(NodeId lhsSmallestNode, NodeId rhsSmallestNode) const {
        requireStableTree();
        if (!tree_.isAlive(lhsSmallestNode) || !tree_.isAlive(rhsSmallestNode)) {
            throw std::invalid_argument("ConnectedSubsetTreeLocalizer::join requires two live nodes.");
        }
        return detail::kernel::connectedSubsetJoin(tree_, lhsSmallestNode, rhsSmallestNode);
    }

    /**
     * @brief Returns the first node on the anchor branch containing an absolute sample.
     * @param anchorPixel Valid anchor pixel.
     * @param samplePixel Valid translated-sample pixel.
     * @return Anchored entry, or no value for an invalid pixel.
     */
    [[nodiscard]] std::optional<NodeId> anchoredEntry(PixelId anchorPixel, PixelId samplePixel) const {
        requireStableTree();
        if (!tree_.isPixel(anchorPixel) || !tree_.isPixel(samplePixel)) {
            return std::nullopt;
        }
        return detail::kernel::anchoredEntry(tree_, anchorPixel, samplePixel);
    }

    /**
     * @brief Returns the first node containing a translated sample.
     * @param anchorPixel Valid anchor pixel.
     * @param offset Relative row-column sample offset.
     * @return Anchored entry, or no value outside the grid.
     */
    [[nodiscard]] std::optional<NodeId> anchoredEntry(PixelId anchorPixel, WindowOffset offset) const {
        requireStableTree();
        const GridDomain2D& domain = tree_.requireGridDomain2D("ConnectedSubsetTreeLocalizer::anchoredEntry");
        if (anchorPixel < 0 || anchorPixel >= domain.rows * domain.columns) {
            return std::nullopt;
        }
        const NodeId entry = detail::kernel::anchoredEntry(tree_, anchorPixel, offset);
        return entry == InvalidNode ? std::nullopt : std::optional<NodeId>{entry};
    }

  private:
    /** @brief Rejects queries after the captured tree has mutated. */
    void requireStableTree() const { tree_.requireMutationVersion(mutationVersion_, "ConnectedSubsetTreeLocalizer"); }

    const MorphologicalTree& tree_;       ///< Captured connected-subset tree.
    std::size_t mutationVersion_ = 0;     ///< Mutation version captured at construction.
};

} // namespace mmcfilters::local_attributes
