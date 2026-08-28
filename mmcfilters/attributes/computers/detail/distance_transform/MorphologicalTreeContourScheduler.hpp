#pragma once

#include "MorphologicalTreeRegionIndex.hpp"
#include "../../../../contours/detail/MorphologicalTreeBoundaryLifetimeIndex.hpp"
#include "../../../../trees/MorphologicalTree.hpp"
#include "../../../../utils/Common.hpp"

#include <cstddef>
#include <memory>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcfilters::attributes::computers::detail::distance_transform {

/**
 * @brief Borrowed contour of one node during a bottom-up hierarchy schedule.
 *
 * The frame and its pixel span are valid only for the synchronous consumer
 * callback. The backing vector may be reused by the parent immediately after
 * the callback returns.
 */
class NodeContourFrame {
  public:
    NodeContourFrame(const NodeContourFrame&) = delete;
    NodeContourFrame& operator=(const NodeContourFrame&) = delete;
    NodeContourFrame(NodeContourFrame&&) = delete;
    NodeContourFrame& operator=(NodeContourFrame&&) = delete;

    /**
     * @brief Returns the represented live node.
     */
    [[nodiscard]] NodeId node() const noexcept { return node_; }

    /**
     * @brief Returns its exact foreground 4-neighbour contour pixels.
     */
    [[nodiscard]] std::span<const PixelId> pixels() const noexcept { return pixels_; }

    /**
     * @brief Tests whether this is the first bottom-up frame emitted on its heavy path.
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
     * @brief Returns sites inserted relative to the immediately preceding heavy child.
     */
    [[nodiscard]] std::span<const PixelId> additionsFromHeavyChild() const noexcept { return additionsFromHeavyChild_; }

    /**
     * @brief Returns sites removed relative to the immediately preceding heavy child.
     */
    [[nodiscard]] std::span<const PixelId> removalsFromHeavyChild() const noexcept { return removalsFromHeavyChild_; }

  private:
    friend class MorphologicalTreeContourScheduler;

    NodeContourFrame(NodeId node, std::span<const PixelId> pixels, NodeId heavyChild, NodeId heavyPathTop, std::span<const PixelId> additionsFromHeavyChild,
                     std::span<const PixelId> removalsFromHeavyChild)
        : node_(node), pixels_(pixels), startsHeavyPath_(heavyChild == InvalidNode), heavyChild_(heavyChild), heavyPathTop_(heavyPathTop),
          additionsFromHeavyChild_(additionsFromHeavyChild), removalsFromHeavyChild_(removalsFromHeavyChild) {}

    NodeId node_ = InvalidNode;
    std::span<const PixelId> pixels_;
    bool startsHeavyPath_ = true;
    NodeId heavyChild_ = InvalidNode;
    NodeId heavyPathTop_ = InvalidNode;
    std::span<const PixelId> additionsFromHeavyChild_;
    std::span<const PixelId> removalsFromHeavyChild_;
};

/**
 * @brief Exact bottom-up contour scheduler shared by every 2D tree kind.
 *
 * Compact local additions/removals are aggregated destructively in post-order.
 * A parent reuses the contour storage of its support-largest child and transfers
 * only the other child contours. Choosing the heavy child by support size makes
 * every transferred pixel enter a support at least twice as large, independently
 * of later contour removals.
 */
class MorphologicalTreeContourScheduler {
  public:
    /**
     * @brief Emits every live-node contour once in post-order.
     *
     * @param tree Stable committed tree on a regular 2D domain.
     * @param regions Support intervals used to select the heavy child.
     * @param consumer Synchronous callback receiving one borrowed contour.
     */
    template <class Consumer> static void forEachNode(const MorphologicalTree& tree, const MorphologicalTreeRegionIndex& regions, Consumer&& consumer) {
        tree.requireNotEditing("MorphologicalTreeContourScheduler::forEachNode");
        regions.requireCompatibleTree(tree);
        forEachEstablishedNode(tree, regions, std::forward<Consumer>(consumer));
    }

    /**
     * @brief Emits contours after the caller established tree/index compatibility.
     * @param tree Established committed tree.
     * @param regions Established support index for `tree`.
     * @param consumer Synchronous callback receiving one borrowed contour.
     */
    template <class Consumer>
    static void forEachEstablishedNode(const MorphologicalTree& tree, const MorphologicalTreeRegionIndex& regions, Consumer&& consumer) {
        const std::size_t mutationVersion = tree.getMutationVersion();
        const ::mmcfilters::contours::detail::MorphologicalTreeBoundaryLifetimeIndex boundaries(tree);
        const std::size_t numSlots = static_cast<std::size_t>(tree.numInternalNodeSlots());
        std::vector<std::unique_ptr<std::vector<PixelId>>> contours(numSlots);
        std::vector<PixelId> positions(static_cast<std::size_t>(tree.numPixels()), InvalidPixel);
        std::vector<NodeId> heavyChildren(numSlots, InvalidNode);
        for (NodeId node : tree.postOrder()) {
            std::size_t heavySupportSize = 0;
            for (NodeId child : ::mmcfilters::detail::CommittedTreeAccess::children(tree, node)) {
                const std::size_t childSupportSize = regions.establishedSupportInterval(child).size();
                if (heavyChildren[static_cast<std::size_t>(node)] == InvalidNode || childSupportSize > heavySupportSize) {
                    heavyChildren[static_cast<std::size_t>(node)] = child;
                    heavySupportSize = childSupportSize;
                }
            }
        }

        std::vector<NodeId> heavyPathTops(numSlots, InvalidNode);
        heavyPathTops[static_cast<std::size_t>(tree.root())] = tree.root();
        for (NodeId node : ::mmcfilters::detail::CommittedTreeAccess::subtree(tree, tree.root())) {
            for (NodeId child : ::mmcfilters::detail::CommittedTreeAccess::children(tree, node)) {
                heavyPathTops[static_cast<std::size_t>(child)] =
                    heavyChildren[static_cast<std::size_t>(node)] == child ? heavyPathTops[static_cast<std::size_t>(node)] : child;
            }
        }

        std::vector<NodeId> scheduledPostOrder;
        scheduledPostOrder.reserve(static_cast<std::size_t>(tree.numNodes()));
        std::vector<std::pair<NodeId, bool>> traversalStack;
        traversalStack.emplace_back(tree.root(), false);
        while (!traversalStack.empty()) {
            const auto [node, expanded] = traversalStack.back();
            traversalStack.pop_back();
            if (expanded) {
                scheduledPostOrder.push_back(node);
                continue;
            }
            traversalStack.emplace_back(node, true);
            const NodeId heavyChild = heavyChildren[static_cast<std::size_t>(node)];
            if (heavyChild != InvalidNode) {
                traversalStack.emplace_back(heavyChild, false);
            }
            for (NodeId child : ::mmcfilters::detail::CommittedTreeAccess::children(tree, node)) {
                if (child != heavyChild) {
                    traversalStack.emplace_back(child, false);
                }
            }
        }

        std::vector<PixelId> transitionAdditions;
        std::vector<PixelId> transitionRemovals;

        for (NodeId node : scheduledPostOrder) {
            const NodeId heavyChild = heavyChildren[static_cast<std::size_t>(node)];
            std::size_t combinedContourSize = boundaries.establishedAdditions(node).size();
            for (NodeId child : ::mmcfilters::detail::CommittedTreeAccess::children(tree, node)) {
                const auto& childContour = contours[static_cast<std::size_t>(child)];
                if (!childContour) {
                    throw std::logic_error("Hierarchical contour scheduling found an unmaterialized child contour.");
                }
                combinedContourSize += childContour->size();
            }

            std::unique_ptr<std::vector<PixelId>> nodeContour;
            if (heavyChild == InvalidNode) {
                nodeContour = std::make_unique<std::vector<PixelId>>();
            } else {
                nodeContour = std::move(contours[static_cast<std::size_t>(heavyChild)]);
            }
            nodeContour->reserve(combinedContourSize);
            transitionAdditions.clear();
            transitionRemovals.clear();

            for (NodeId child : ::mmcfilters::detail::CommittedTreeAccess::children(tree, node)) {
                if (child == heavyChild) {
                    continue;
                }
                auto& childContour = contours[static_cast<std::size_t>(child)];
                for (PixelId pixel : *childContour) {
                    positions[static_cast<std::size_t>(pixel)] = static_cast<PixelId>(nodeContour->size());
                    nodeContour->push_back(pixel);
                    transitionAdditions.push_back(pixel);
                }
                childContour.reset();
            }

            for (PixelId pixel : boundaries.establishedAdditions(node)) {
                if (positions[static_cast<std::size_t>(pixel)] != InvalidPixel) {
                    throw std::logic_error("Hierarchical contour scheduling found an invalid or duplicate local addition.");
                }
                positions[static_cast<std::size_t>(pixel)] = static_cast<PixelId>(nodeContour->size());
                nodeContour->push_back(pixel);
                transitionAdditions.push_back(pixel);
            }

            for (PixelId pixel : boundaries.establishedRemovals(node)) {
                const PixelId position = positions[static_cast<std::size_t>(pixel)];
                if (position == InvalidPixel || static_cast<std::size_t>(position) >= nodeContour->size() ||
                    (*nodeContour)[static_cast<std::size_t>(position)] != pixel) {
                    throw std::logic_error("Hierarchical contour scheduling removed a pixel outside the active node contour.");
                }
                const PixelId movedPixel = nodeContour->back();
                (*nodeContour)[static_cast<std::size_t>(position)] = movedPixel;
                positions[static_cast<std::size_t>(movedPixel)] = position;
                nodeContour->pop_back();
                positions[static_cast<std::size_t>(pixel)] = InvalidPixel;
                transitionRemovals.push_back(pixel);
            }

            if (nodeContour->empty()) {
                throw std::logic_error("Hierarchical contour scheduling produced an empty live-node contour.");
            }
            const NodeContourFrame frame(node, *nodeContour, heavyChild, heavyPathTops[static_cast<std::size_t>(node)], transitionAdditions,
                                         transitionRemovals);
            consumer(frame);
            tree.requireMutationVersion(mutationVersion, "MorphologicalTreeContourScheduler::forEachNode callback");
            contours[static_cast<std::size_t>(node)] = std::move(nodeContour);
        }
    }
};

} // namespace mmcfilters::attributes::computers::detail::distance_transform
