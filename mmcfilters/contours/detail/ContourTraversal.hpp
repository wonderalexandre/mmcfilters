#pragma once

#include "NodeSupportIndex.hpp"
#include "ContourLifetimeIndex.hpp"

#include <algorithm>
#include <cstddef>
#include <memory>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcfilters::contours::detail {

/**
 * @brief Resumable post-order contour traversal with O(P + N) working storage.
 *
 * The child with the largest support is the heavy child. The parent reuses its
 * vector and releases the other child vectors after merging them. The traversal
 * retains no emitted contours. The source tree and supplied indexes must outlive
 * the traversal.
 */
class ContourTraversal {
  public:
    /**
     * @brief Builds the traversal schedule without computing a contour.
     * @param tree Stable source tree.
     * @param supportIndex Support index built from `tree`.
     * @param contourLifetimes Contour lifetimes built from `tree`.
     */
    ContourTraversal(const MorphologicalTree& tree, const NodeSupportIndex& supportIndex,
                     const ContourLifetimeIndex& contourLifetimes)
        : tree_(tree), contourLifetimes_(contourLifetimes), mutationVersion_(tree.getMutationVersion()),
          contourByNode_(static_cast<std::size_t>(tree.numInternalNodeSlots())),
          positionByPixel_(static_cast<std::size_t>(tree.numPixels()), InvalidPixel) {
        supportIndex.requireCompatibleTree(tree);
        const std::size_t numSlots = static_cast<std::size_t>(tree.numInternalNodeSlots());
        heavyChildByNode_.assign(numSlots, InvalidNode);
        for (NodeId node : tree.postOrder()) {
            std::size_t heavySupportSize = 0;
            for (NodeId child : ::mmcfilters::detail::CommittedTreeAccess::children(tree, node)) {
                const std::size_t childSupportSize = supportIndex.establishedSupportInterval(child).size();
                if (heavyChildByNode_[static_cast<std::size_t>(node)] == InvalidNode || childSupportSize > heavySupportSize) {
                    heavyChildByNode_[static_cast<std::size_t>(node)] = child;
                    heavySupportSize = childSupportSize;
                }
            }
        }

        postOrder_.reserve(static_cast<std::size_t>(tree.numNodes()));
        std::vector<std::pair<NodeId, bool>> traversalStack;
        traversalStack.emplace_back(tree.root(), false);
        while (!traversalStack.empty()) {
            const auto [node, expanded] = traversalStack.back();
            traversalStack.pop_back();
            if (expanded) {
                postOrder_.push_back(node);
                continue;
            }
            traversalStack.emplace_back(node, true);
            const NodeId heavyChild = heavyChildByNode_[static_cast<std::size_t>(node)];
            if (heavyChild != InvalidNode) {
                traversalStack.emplace_back(heavyChild, false);
            }
            for (NodeId child : ::mmcfilters::detail::CommittedTreeAccess::children(tree, node)) {
                if (child != heavyChild) {
                    traversalStack.emplace_back(child, false);
                }
            }
        }
    }

    /**
     * @brief Advances to the next live node contour.
     * @return True when a new current contour is available.
     */
    bool advance() { return advance([](PixelId) {}, [](PixelId) {}); }

    /**
     * @brief Advances and reports changes from the selected heavy child.
     * @param onAdded Called for every pixel added while building the contour.
     * @param onRemoved Called for every pixel removed while building the contour.
     * @return True when a new current contour is available.
     */
    template <class AdditionConsumer, class RemovalConsumer>
    bool advance(AdditionConsumer&& onAdded, RemovalConsumer&& onRemoved) {
        requireStableTree();
        currentNode_ = InvalidNode;
        if (nextNodeIndex_ == postOrder_.size()) {
            // Release the root's final vector when the cursor reaches its end.
            contourByNode_.clear();
            return false;
        }
        const NodeId node = postOrder_[nextNodeIndex_++];
        const NodeId heavyChild = heavyChildByNode_[static_cast<std::size_t>(node)];
        std::size_t requiredCapacity = contourLifetimes_.establishedAdditions(node).size();
        for (NodeId child : ::mmcfilters::detail::CommittedTreeAccess::children(tree_, node)) {
            const auto& childContour = contourByNode_[static_cast<std::size_t>(child)];
            if (!childContour) {
                throw std::logic_error("Contour traversal found a missing child contour.");
            }
            requiredCapacity += childContour->size();
        }

        std::unique_ptr<std::vector<PixelId>> nodeContour;
        if (heavyChild == InvalidNode) {
            nodeContour = std::make_unique<std::vector<PixelId>>();
        } else {
            nodeContour = std::move(contourByNode_[static_cast<std::size_t>(heavyChild)]);
        }
        // Grow geometrically when the reused child needs more capacity.
        if (requiredCapacity > nodeContour->capacity()) {
            nodeContour->reserve(std::max(requiredCapacity, nodeContour->capacity() * 2));
        }

        for (NodeId child : ::mmcfilters::detail::CommittedTreeAccess::children(tree_, node)) {
            if (child == heavyChild) {
                continue;
            }
            auto& childContour = contourByNode_[static_cast<std::size_t>(child)];
            for (PixelId pixel : *childContour) {
                positionByPixel_[static_cast<std::size_t>(pixel)] = static_cast<PixelId>(nodeContour->size());
                nodeContour->push_back(pixel);
                onAdded(pixel);
            }
            childContour.reset();
        }

        for (PixelId pixel : contourLifetimes_.establishedAdditions(node)) {
            if (positionByPixel_[static_cast<std::size_t>(pixel)] != InvalidPixel) {
                throw std::logic_error("Contour traversal found an invalid or duplicate node addition.");
            }
            positionByPixel_[static_cast<std::size_t>(pixel)] = static_cast<PixelId>(nodeContour->size());
            nodeContour->push_back(pixel);
            onAdded(pixel);
        }

        for (PixelId pixel : contourLifetimes_.establishedRemovals(node)) {
            const PixelId position = positionByPixel_[static_cast<std::size_t>(pixel)];
            if (position == InvalidPixel || static_cast<std::size_t>(position) >= nodeContour->size() ||
                (*nodeContour)[static_cast<std::size_t>(position)] != pixel) {
                throw std::logic_error("Contour traversal removed a pixel outside the active node contour.");
            }
            const PixelId movedPixel = nodeContour->back();
            (*nodeContour)[static_cast<std::size_t>(position)] = movedPixel;
            positionByPixel_[static_cast<std::size_t>(movedPixel)] = position;
            nodeContour->pop_back();
            positionByPixel_[static_cast<std::size_t>(pixel)] = InvalidPixel;
            onRemoved(pixel);
        }

        if (nodeContour->empty()) {
            throw std::logic_error("Contour traversal produced an empty live-node contour.");
        }
        contourByNode_[static_cast<std::size_t>(node)] = std::move(nodeContour);
        currentNode_ = node;
        return true;
    }

    /**
     * @brief Borrows the current node and contour.
     * @return Current node and its borrowed contour pixels.
     */
    [[nodiscard]] std::pair<NodeId, std::span<const PixelId>> current() const {
        requireStableTree();
        if (currentNode_ == InvalidNode) {
            throw std::out_of_range("Contour traversal has no current node.");
        }
        return {currentNode_, *contourByNode_[static_cast<std::size_t>(currentNode_)]};
    }

    /**
     * @brief Returns the selected heavy child of an established live node.
     * @param node Node already established as live.
     * @return Child with the largest support, or `InvalidNode` for a leaf.
     */
    [[nodiscard]] NodeId heavyChild(NodeId node) const noexcept { return heavyChildByNode_[static_cast<std::size_t>(node)]; }

    /**
     * @brief Returns the child-before-parent traversal schedule.
     * @return Borrowed post-order node sequence.
     */
    [[nodiscard]] std::span<const NodeId> postOrder() const noexcept { return postOrder_; }

    /** @brief Rejects access after a topology mutation. */
    void requireStableTree() const { tree_.requireMutationVersion(mutationVersion_, "ContourTraversal"); }

  private:
    const MorphologicalTree& tree_;                 ///< Source tree.
    const ContourLifetimeIndex& contourLifetimes_;  ///< Contour membership intervals and events.
    std::size_t mutationVersion_;                   ///< Tree mutation version captured at construction.
    std::vector<std::unique_ptr<std::vector<PixelId>>> contourByNode_; ///< Live temporary contour for each node slot.
    std::vector<PixelId> positionByPixel_;          ///< Position of each pixel in the active reused contour.
    std::vector<NodeId> heavyChildByNode_;          ///< Largest-support child selected for each node.
    std::vector<NodeId> postOrder_;                 ///< Child-before-parent traversal schedule.
    std::size_t nextNodeIndex_ = 0;                 ///< Next position in `postOrder_`.
    NodeId currentNode_ = InvalidNode;              ///< Node yielded by the most recent successful advance.
};

} // namespace mmcfilters::contours::detail
