#pragma once

#include "NodeSupportBoxIndex.hpp"
#include "../../../../contours/detail/ContourTraversal.hpp"

namespace mmcfilters::attributes::computers::detail::distance_transform {

/**
 * @brief One contour and its changes from the preceding heavy child.
 *
 * The update and its pixel spans are valid until the traversal advances or is
 * destroyed. The parent may reuse the backing vector after advancement.
 */
class NodeContourUpdate {
  public:
    NodeContourUpdate(const NodeContourUpdate&) = delete;
    NodeContourUpdate& operator=(const NodeContourUpdate&) = delete;
    NodeContourUpdate(NodeContourUpdate&&) = delete;
    NodeContourUpdate& operator=(NodeContourUpdate&&) = delete;

    /**
     * @brief Returns the represented live node.
     */
    [[nodiscard]] NodeId node() const noexcept { return node_; }

    /**
     * @brief Returns the node's foreground A4 contour.
     */
    [[nodiscard]] std::span<const PixelId> contour() const noexcept { return contourPixels_; }

    /**
     * @brief Tests whether the update starts a heavy path.
     */
    [[nodiscard]] bool startsHeavyPath() const noexcept { return startsHeavyPath_; }

    /**
     * @brief Returns the direct heavy child reused by this transition, if any.
     */
    [[nodiscard]] NodeId heavyChild() const noexcept { return heavyChild_; }

    /**
     * @brief Returns the topmost node of this node's heavy path.
     */
    [[nodiscard]] NodeId heavyPathTop() const noexcept { return heavyPathTop_; }

    /**
     * @brief Returns contour pixels added to the preceding heavy child's contour.
     */
    [[nodiscard]] std::span<const PixelId> addedContourPixels() const noexcept { return addedContourPixels_; }

    /**
     * @brief Returns contour pixels removed from the preceding heavy child's contour.
     */
    [[nodiscard]] std::span<const PixelId> removedContourPixels() const noexcept { return removedContourPixels_; }

  private:
    template <class Consumer>
    friend void forEachEstablishedContourUpdate(const MorphologicalTree&, const NodeSupportBoxIndex&, Consumer&&);

    NodeContourUpdate(NodeId node, std::span<const PixelId> contourPixels, NodeId heavyChild, NodeId heavyPathTop,
                     std::span<const PixelId> addedContourPixels,
                     std::span<const PixelId> removedContourPixels)
        : node_(node), contourPixels_(contourPixels), startsHeavyPath_(heavyChild == InvalidNode), heavyChild_(heavyChild), heavyPathTop_(heavyPathTop),
          addedContourPixels_(addedContourPixels), removedContourPixels_(removedContourPixels) {}

    NodeId node_ = InvalidNode;
    std::span<const PixelId> contourPixels_;
    bool startsHeavyPath_ = true;
    NodeId heavyChild_ = InvalidNode;
    NodeId heavyPathTop_ = InvalidNode;
    std::span<const PixelId> addedContourPixels_;
    std::span<const PixelId> removedContourPixels_;
};

/** @brief Emits contour updates after the caller establishes tree and index compatibility. */
template <class Consumer>
void forEachEstablishedContourUpdate(const MorphologicalTree& tree, const NodeSupportBoxIndex& supportIndex, Consumer&& consumer) {
    const ::mmcfilters::contours::detail::ContourLifetimeIndex contourLifetimes(tree);
    ::mmcfilters::contours::detail::ContourTraversal traversal(tree, supportIndex, contourLifetimes);
    std::vector<NodeId> heavyPathTopByNode(static_cast<std::size_t>(tree.numInternalNodeSlots()), InvalidNode);
    const auto postOrder = traversal.postOrder();
    for (auto nodeIterator = postOrder.rbegin(); nodeIterator != postOrder.rend(); ++nodeIterator) {
        const NodeId node = *nodeIterator;
        const NodeId parent = ::mmcfilters::detail::CommittedTreeAccess::nodeParent(tree, node);
        heavyPathTopByNode[static_cast<std::size_t>(node)] =
            node != tree.root() && traversal.heavyChild(parent) == node ? heavyPathTopByNode[static_cast<std::size_t>(parent)] : node;
    }
    std::vector<PixelId> addedPixels;
    std::vector<PixelId> removedPixels;
    for (;;) {
        addedPixels.clear();
        removedPixels.clear();
        if (!traversal.advance([&](PixelId pixel) { addedPixels.push_back(pixel); }, [&](PixelId pixel) { removedPixels.push_back(pixel); })) {
            break;
        }
        const auto [node, contourPixels] = traversal.current();
        const NodeContourUpdate update(node, contourPixels, traversal.heavyChild(node), heavyPathTopByNode[static_cast<std::size_t>(node)], addedPixels,
                                       removedPixels);
        consumer(update);
        traversal.requireStableTree();
    }
}

/** @brief Emits the contour update required for each distance transform. */
template <class Consumer>
void forEachContourUpdate(const MorphologicalTree& tree, const NodeSupportBoxIndex& supportIndex, Consumer&& consumer) {
    tree.requireNotEditing("forEachContourUpdate");
    supportIndex.requireCompatibleTree(tree);
    forEachEstablishedContourUpdate(tree, supportIndex, std::forward<Consumer>(consumer));
}

} // namespace mmcfilters::attributes::computers::detail::distance_transform
