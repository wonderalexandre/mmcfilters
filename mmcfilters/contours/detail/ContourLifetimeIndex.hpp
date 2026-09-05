#pragma once

#include "../../trees/MorphologicalTree.hpp"
#include "../../trees/detail/CommittedTreeAccess.hpp"
#include "../../utils/Common.hpp"

#include <array>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcfilters::contours::detail {

/**
 * @brief One pixel's exact contour lifetime in the node hierarchy.
 *
 * A lifetime starts at the pixel's smallest node and follows its ancestor path
 * up to, but excluding, `endNodeExclusive`. `InvalidNode` means that the
 * lifetime reaches the root. Equal start and end nodes mean that the pixel is
 * interior in every support that contains it.
 */
struct ContourLifetime {
    NodeId startNode = InvalidNode;         ///< Smallest node where the lifetime starts.
    NodeId endNodeExclusive = InvalidNode; ///< First ancestor support where the pixel is no longer on the boundary.

    /**
     * @brief Tests whether the pixel is a contour site in at least one node.
     * @return True when the lifetime contains at least one node.
     */
    [[nodiscard]] bool nonEmpty() const noexcept { return startNode != endNodeExclusive; }
};

/**
 * @brief Exact A4 contour-lifetime index derived only from tree topology.
 *
 * Let o(p) be the smallest node of p. For a pixel not touching the global
 * domain boundary, its first interior node is the LCA of o(p) and the smallest
 * nodes of all four side neighbors. Thus its contour activity is one connected
 * path interval. A pixel on the image boundary has an exterior neighbor that never
 * enters a support, so its interval reaches the root.
 */
class ContourLifetimeIndex {
  public:
    ContourLifetimeIndex(const ContourLifetimeIndex&) = delete;
    ContourLifetimeIndex& operator=(const ContourLifetimeIndex&) = delete;
    ContourLifetimeIndex(ContourLifetimeIndex&&) = delete;
    ContourLifetimeIndex& operator=(ContourLifetimeIndex&&) = delete;

    /**
     * @brief Builds all lifetimes and compact start/stop event lists.
     * @param tree Stable tree with a two-dimensional grid domain.
     */
    explicit ContourLifetimeIndex(const MorphologicalTree& tree)
        : tree_(tree), mutationVersion_(tree.getMutationVersion()), domain_(requireDomain(tree)), lifetimeByPixel_(static_cast<std::size_t>(tree.numPixels())),
          preorderIndexByNode_(static_cast<std::size_t>(tree.numInternalNodeSlots()), -1),
          subtreeEndIndexByNode_(static_cast<std::size_t>(tree.numInternalNodeSlots()), -1),
          additionOffsets_(static_cast<std::size_t>(tree.numInternalNodeSlots()) + 1, 0),
          removalOffsets_(static_cast<std::size_t>(tree.numInternalNodeSlots()) + 1, 0) {
        const std::span<const NodeId> smallestNodes = tree_.smallestNodeMap();
        for (PixelId pixel = 0; pixel < tree_.numPixels(); ++pixel) {
            lifetimeByPixel_[static_cast<std::size_t>(pixel)].startNode = smallestNodes[static_cast<std::size_t>(pixel)];
        }

        buildPreOrderIntervals(tree_, preorderIndexByNode_, subtreeEndIndexByNode_);
        const std::size_t interiorRows = domain_.rows > 2 ? static_cast<std::size_t>(domain_.rows - 2) : 0;
        const std::size_t interiorColumns = domain_.columns > 2 ? static_cast<std::size_t>(domain_.columns - 2) : 0;
        const std::size_t numLcaQueries = interiorRows * interiorColumns;
        const auto queryPixel = [&](std::size_t queryIndex) {
            const int row = static_cast<int>(queryIndex / interiorColumns) + 1;
            const int column = static_cast<int>(queryIndex % interiorColumns) + 1;
            return static_cast<PixelId>(row * domain_.columns + column);
        };
        const auto lcaQuery = [&](std::size_t queryIndex) {
            const PixelId pixel = queryPixel(queryIndex);
            const NodeId startNode = smallestNodes[static_cast<std::size_t>(pixel)];
            const std::array<NodeId, 5> neighborhoodSmallestNodes{
                startNode,
                smallestNodes[static_cast<std::size_t>(pixel - domain_.columns)],
                smallestNodes[static_cast<std::size_t>(pixel + domain_.columns)],
                smallestNodes[static_cast<std::size_t>(pixel - 1)],
                smallestNodes[static_cast<std::size_t>(pixel + 1)],
            };
            NodeId firstInPreOrder = startNode;
            NodeId lastInPreOrder = startNode;
            for (NodeId smallestNode : neighborhoodSmallestNodes) {
                if (preorderIndexByNode_[static_cast<std::size_t>(smallestNode)] < preorderIndexByNode_[static_cast<std::size_t>(firstInPreOrder)]) {
                    firstInPreOrder = smallestNode;
                }
                if (preorderIndexByNode_[static_cast<std::size_t>(smallestNode)] > preorderIndexByNode_[static_cast<std::size_t>(lastInPreOrder)]) {
                    lastInPreOrder = smallestNode;
                }
            }
            return std::pair<NodeId, NodeId>{firstInPreOrder, lastInPreOrder};
        };
        ::mmcfilters::detail::CommittedTreeAccess::forEachLowestCommonAncestor(
            tree_, numLcaQueries, lcaQuery,
            [&](std::size_t queryIndex, NodeId endNodeExclusive) {
                lifetimeByPixel_[static_cast<std::size_t>(queryPixel(queryIndex))].endNodeExclusive = endNodeExclusive;
            });

        for (PixelId pixel = 0; pixel < tree_.numPixels(); ++pixel) {
            const ContourLifetime& lifetime = lifetimeByPixel_[static_cast<std::size_t>(pixel)];
            const NodeId startNode = lifetime.startNode;
            const NodeId endNodeExclusive = lifetime.endNodeExclusive;
            if (startNode != endNodeExclusive) {
                ++additionOffsets_[static_cast<std::size_t>(startNode) + 1];
                if (endNodeExclusive != InvalidNode) {
                    ++removalOffsets_[static_cast<std::size_t>(endNodeExclusive) + 1];
                }
            }
        }

        prefixSum(additionOffsets_);
        prefixSum(removalOffsets_);
        additionPixels_.resize(additionOffsets_.back());
        removalPixels_.resize(removalOffsets_.back());
        std::vector<std::size_t> nextAdditionOffset = additionOffsets_;
        std::vector<std::size_t> nextRemovalOffset = removalOffsets_;
        for (PixelId pixel = 0; pixel < tree_.numPixels(); ++pixel) {
            const ContourLifetime& lifetime = lifetimeByPixel_[static_cast<std::size_t>(pixel)];
            if (!lifetime.nonEmpty()) {
                continue;
            }
            additionPixels_[nextAdditionOffset[static_cast<std::size_t>(lifetime.startNode)]++] = pixel;
            if (lifetime.endNodeExclusive != InvalidNode) {
                removalPixels_[nextRemovalOffset[static_cast<std::size_t>(lifetime.endNodeExclusive)]++] = pixel;
            }
        }
    }

    /**
     * @brief Returns the captured lifetime of a valid domain pixel.
     * @param pixel Valid row-major pixel identifier.
     * @return Captured contour lifetime.
     */
    [[nodiscard]] const ContourLifetime& lifetime(PixelId pixel) const {
        requireStableTree();
        if (!tree_.isPixel(pixel)) {
            throw std::out_of_range("Contour lifetime index received an invalid pixel id.");
        }
        return lifetimeByPixel_[static_cast<std::size_t>(pixel)];
    }

    /**
     * @brief Returns pixels whose contour lifetime starts at one live node.
     * @param node Live start-event node.
     * @return Borrowed span of pixel identifiers.
     */
    [[nodiscard]] std::span<const PixelId> additions(NodeId node) const {
        requireLiveNode(node);
        return eventSpan(additionPixels_, additionOffsets_, node);
    }

    /**
     * @brief Returns pixels whose contour lifetime stops before one live node.
     * @param node Live stop-event node.
     * @return Borrowed span of pixel identifiers.
     */
    [[nodiscard]] std::span<const PixelId> removals(NodeId node) const {
        requireLiveNode(node);
        return eventSpan(removalPixels_, removalOffsets_, node);
    }

    /**
     * @brief Returns start events for a caller-established live node.
     * @param node Established live node identifier.
     * @return Borrowed span of pixel identifiers.
     */
    [[nodiscard]] std::span<const PixelId> establishedAdditions(NodeId node) const noexcept {
        return eventSpan(additionPixels_, additionOffsets_, node);
    }

    /**
     * @brief Returns stop events for a caller-established live node.
     * @param node Established live node identifier.
     * @return Borrowed span of pixel identifiers.
     */
    [[nodiscard]] std::span<const PixelId> establishedRemovals(NodeId node) const noexcept {
        return eventSpan(removalPixels_, removalOffsets_, node);
    }

    /**
     * @brief Tests whether an established pixel is on an established node contour.
     * @param pixel Established row-major pixel identifier.
     * @param node Established live support node.
     * @return True when the pixel belongs to the node contour.
     */
    [[nodiscard]] bool establishedIsContourPixel(PixelId pixel, NodeId node) const noexcept {
        const ContourLifetime& value = lifetimeByPixel_[static_cast<std::size_t>(pixel)];
        if (!isAncestorInCapturedHierarchy(node, value.startNode)) {
            return false;
        }
        if (value.endNodeExclusive == InvalidNode) {
            return true;
        }
        return value.endNodeExclusive != node && isAncestorInCapturedHierarchy(value.endNodeExclusive, node);
    }

    /**
     * @brief Tests whether one pixel is on a live node contour.
     * @param pixel Valid row-major pixel identifier.
     * @param node Live support node.
     * @return True when the pixel belongs to the node contour.
     */
    [[nodiscard]] bool isContourPixel(PixelId pixel, NodeId node) const {
        requireLiveNode(node);
        if (!tree_.isPixel(pixel)) {
            throw std::out_of_range("Contour lifetime index received an invalid pixel id.");
        }
        return establishedIsContourPixel(pixel, node);
    }

    /**
     * @brief Returns the number of contour lifetime start events.
     * @return Number of stored additions.
     */
    [[nodiscard]] std::size_t numAdditions() const {
        requireStableTree();
        return additionPixels_.size();
    }

    /**
     * @brief Returns the number of finite lifetime-stop events.
     * @return Number of stored removals.
     */
    [[nodiscard]] std::size_t numRemovals() const {
        requireStableTree();
        return removalPixels_.size();
    }

  private:
    /**
     * @brief Validates the source tree and returns its grid domain.
     * @param tree Source hierarchy.
     * @return Validated two-dimensional domain.
     */
    [[nodiscard]] static GridDomain2D requireDomain(const MorphologicalTree& tree) {
        tree.requireNotEditing("ContourLifetimeIndex");
        return tree.requireGridDomain2D("ContourLifetimeIndex");
    }

    /**
     * @brief Assigns a dense preorder rank to every live node without recursion.
     *
     * For any non-empty node set, the LCA of the complete set equals the LCA of
     * its first and last nodes in depth-first preorder. This reduces the query
     * over a pixel and its four side neighbors to one `MorphologicalTree` LCA
     * query.
     * @param tree Stable hierarchy to traverse.
     * @param entryIndexByNode Dense output array of preorder entry ranks.
     * @param subtreeEndIndexByNode Dense output array of exclusive subtree-end ranks.
     */
    static void buildPreOrderIntervals(const MorphologicalTree& tree, std::vector<int>& entryIndexByNode, std::vector<int>& subtreeEndIndexByNode) {
        int nextPreorderIndex = 0;
        std::vector<std::pair<NodeId, bool>> traversalStack;
        traversalStack.reserve(static_cast<std::size_t>(tree.numNodes()) + 1);
        traversalStack.emplace_back(tree.root(), false);
        while (!traversalStack.empty()) {
            const auto [node, exiting] = traversalStack.back();
            traversalStack.pop_back();
            const std::size_t nodeIndex = static_cast<std::size_t>(node);
            if (exiting) {
                subtreeEndIndexByNode[nodeIndex] = nextPreorderIndex;
                continue;
            }
            entryIndexByNode[nodeIndex] = nextPreorderIndex++;
            traversalStack.emplace_back(node, true);
            for (NodeId child : ::mmcfilters::detail::CommittedTreeAccess::children(tree, node)) {
                traversalStack.emplace_back(child, false);
            }
        }
    }

    /**
     * @brief Tests ancestry using the index-owned iterative DFS intervals.
     * @param ancestor Candidate ancestor.
     * @param descendant Candidate descendant.
     * @return True when `ancestor` contains `descendant`.
     */
    [[nodiscard]] bool isAncestorInCapturedHierarchy(NodeId ancestor, NodeId descendant) const noexcept {
        const int ancestorEntry = preorderIndexByNode_[static_cast<std::size_t>(ancestor)];
        const int descendantEntry = preorderIndexByNode_[static_cast<std::size_t>(descendant)];
        return ancestorEntry <= descendantEntry && descendantEntry < subtreeEndIndexByNode_[static_cast<std::size_t>(ancestor)];
    }

    /**
     * @brief Converts per-slot counts to CSR offsets in place.
     * @param offsets Count array with one sentinel slot.
     */
    static void prefixSum(std::vector<std::size_t>& offsets) {
        for (std::size_t index = 1; index < offsets.size(); ++index) {
            offsets[index] += offsets[index - 1];
        }
    }

    /**
     * @brief Returns one node's borrowed range from compact event storage.
     * @param values Compact event values.
     * @param offsets CSR offsets.
     * @param node Dense node slot.
     * @return Borrowed event span.
     */
    [[nodiscard]] static std::span<const PixelId> eventSpan(const std::vector<PixelId>& values, const std::vector<std::size_t>& offsets,
                                                            NodeId node) noexcept {
        const std::size_t nodeIndex = static_cast<std::size_t>(node);
        return std::span<const PixelId>(values).subspan(offsets[nodeIndex], offsets[nodeIndex + 1] - offsets[nodeIndex]);
    }

    /**
     * @brief Rejects queries after the captured tree has mutated.
     */
    void requireStableTree() const { tree_.requireMutationVersion(mutationVersion_, "ContourLifetimeIndex"); }

    /**
     * @brief Validates captured-tree stability and a node identifier.
     * @param node Node that must still be live.
     */
    void requireLiveNode(NodeId node) const {
        requireStableTree();
        if (!tree_.isAlive(node)) {
            throw std::out_of_range("Contour lifetime index received a non-live node id.");
        }
    }

    const MorphologicalTree& tree_;                 ///< Captured source hierarchy.
    std::size_t mutationVersion_ = 0;               ///< Mutation version captured at construction.
    GridDomain2D domain_{};                         ///< Validated two-dimensional image domain.
    std::vector<ContourLifetime> lifetimeByPixel_;       ///< Dense lifetime table indexed by pixel.
    std::vector<int> preorderIndexByNode_;              ///< Dense preorder entry ranks indexed by node.
    std::vector<int> subtreeEndIndexByNode_;            ///< Dense exclusive subtree-end ranks indexed by node.
    std::vector<std::size_t> additionOffsets_;      ///< CSR offsets for lifetime-start events.
    std::vector<std::size_t> removalOffsets_;       ///< CSR offsets for lifetime-stop events.
    std::vector<PixelId> additionPixels_;           ///< Compact lifetime-start pixel identifiers.
    std::vector<PixelId> removalPixels_;            ///< Compact lifetime-stop pixel identifiers.
};

} // namespace mmcfilters::contours::detail
