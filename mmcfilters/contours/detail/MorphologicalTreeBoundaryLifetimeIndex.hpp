#pragma once

#include "../../trees/MorphologicalTree.hpp"
#include "../../utils/Common.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcfilters::contours::detail {

/**
 * @brief One pixel's exact foreground-contour lifetime in the node hierarchy.
 *
 * An active lifetime starts at the pixel owner and follows its ancestor path up
 * to, but excluding, `stopExclusive`. `InvalidNode` means that the lifetime
 * reaches the root. When owner equals stopExclusive, the pixel is interior in
 * every support that contains it.
 */
struct BoundaryLifetime {
    NodeId owner = InvalidNode;
    NodeId stopExclusive = InvalidNode;

    /** @brief Tests whether the pixel is a contour site in at least one node. */
    [[nodiscard]] bool active() const noexcept { return owner != stopExclusive; }
};

/**
 * @brief Exact A4 contour-lifetime index derived only from tree topology.
 *
 * Let o(p) be the inclusion-smallest owner of p. For a pixel not touching the
 * global domain boundary, its first interior node is the LCA of o(p) and the
 * owners of all four side-neighbours. Thus its contour activity is one connected
 * path interval. A global-boundary pixel has an exterior neighbour that never
 * enters a support, so its interval reaches the root.
 */
class MorphologicalTreeBoundaryLifetimeIndex {
  public:
    MorphologicalTreeBoundaryLifetimeIndex(const MorphologicalTreeBoundaryLifetimeIndex&) = delete;
    MorphologicalTreeBoundaryLifetimeIndex& operator=(const MorphologicalTreeBoundaryLifetimeIndex&) = delete;
    MorphologicalTreeBoundaryLifetimeIndex(MorphologicalTreeBoundaryLifetimeIndex&&) = delete;
    MorphologicalTreeBoundaryLifetimeIndex& operator=(MorphologicalTreeBoundaryLifetimeIndex&&) = delete;

    /** @brief Builds all lifetimes and compact start/stop event lists. */
    explicit MorphologicalTreeBoundaryLifetimeIndex(const MorphologicalTree& tree)
        : tree_(tree), mutationVersion_(tree.getMutationVersion()), domain_(requireDomain(tree)), lifetimes_(static_cast<std::size_t>(tree.numPixels())),
          preOrderIndices_(static_cast<std::size_t>(tree.numInternalNodeSlots()), -1),
          subtreeEndIndices_(static_cast<std::size_t>(tree.numInternalNodeSlots()), -1),
          additionOffsets_(static_cast<std::size_t>(tree.numInternalNodeSlots()) + 1, 0),
          removalOffsets_(static_cast<std::size_t>(tree.numInternalNodeSlots()) + 1, 0) {
        const std::span<const NodeId> owners = tree_.smallestNodeMap();
        for (PixelId pixel = 0; pixel < tree_.numPixels(); ++pixel) {
            const NodeId owner = owners[static_cast<std::size_t>(pixel)];
            if (!tree_.isAlive(owner)) {
                throw std::logic_error("Boundary-lifetime indexing found a pixel without a live owner.");
            }

            lifetimes_[static_cast<std::size_t>(pixel)].owner = owner;
        }

        buildPreOrderIntervals(tree_, preOrderIndices_, subtreeEndIndices_);
        std::vector<OfflineLcaQuery> offlineQueries;
        offlineQueries.reserve(static_cast<std::size_t>(tree_.numPixels()));
        for (PixelId pixel = 0; pixel < tree_.numPixels(); ++pixel) {
            const NodeId owner = owners[static_cast<std::size_t>(pixel)];

            const int row = pixel / domain_.columns;
            const int column = pixel % domain_.columns;
            NodeId stopExclusive = InvalidNode;
            if (row > 0 && row + 1 < domain_.rows && column > 0 && column + 1 < domain_.columns) {
                const std::array<NodeId, 5> neighborhoodOwners{
                    owner,
                    owners[static_cast<std::size_t>(pixel - domain_.columns)],
                    owners[static_cast<std::size_t>(pixel + domain_.columns)],
                    owners[static_cast<std::size_t>(pixel - 1)],
                    owners[static_cast<std::size_t>(pixel + 1)],
                };
                NodeId firstInPreOrder = owner;
                NodeId lastInPreOrder = owner;
                for (NodeId neighborhoodOwner : neighborhoodOwners) {
                    if (!tree_.isAlive(neighborhoodOwner) || preOrderIndices_[static_cast<std::size_t>(neighborhoodOwner)] < 0) {
                        throw std::logic_error("Boundary-lifetime indexing found a neighbour without a live connected owner.");
                    }
                    if (preOrderIndices_[static_cast<std::size_t>(neighborhoodOwner)] < preOrderIndices_[static_cast<std::size_t>(firstInPreOrder)]) {
                        firstInPreOrder = neighborhoodOwner;
                    }
                    if (preOrderIndices_[static_cast<std::size_t>(neighborhoodOwner)] > preOrderIndices_[static_cast<std::size_t>(lastInPreOrder)]) {
                        lastInPreOrder = neighborhoodOwner;
                    }
                }

                if (firstInPreOrder == lastInPreOrder) {
                    stopExclusive = firstInPreOrder;
                } else {
                    offlineQueries.push_back({firstInPreOrder, lastInPreOrder, pixel});
                }
            }

            lifetimes_[static_cast<std::size_t>(pixel)].stopExclusive = stopExclusive;
        }
        resolveOfflineLcas(tree_, offlineQueries, lifetimes_);

        for (PixelId pixel = 0; pixel < tree_.numPixels(); ++pixel) {
            const BoundaryLifetime& lifetime = lifetimes_[static_cast<std::size_t>(pixel)];
            const NodeId owner = lifetime.owner;
            const NodeId stopExclusive = lifetime.stopExclusive;
            if (owner != stopExclusive) {
                ++additionOffsets_[static_cast<std::size_t>(owner) + 1];
                if (stopExclusive != InvalidNode) {
                    ++removalOffsets_[static_cast<std::size_t>(stopExclusive) + 1];
                }
            }
        }

        prefixSum(additionOffsets_);
        prefixSum(removalOffsets_);
        additionPixels_.resize(additionOffsets_.back());
        removalPixels_.resize(removalOffsets_.back());
        std::vector<std::size_t> additionCursor = additionOffsets_;
        std::vector<std::size_t> removalCursor = removalOffsets_;
        for (PixelId pixel = 0; pixel < tree_.numPixels(); ++pixel) {
            const BoundaryLifetime& lifetime = lifetimes_[static_cast<std::size_t>(pixel)];
            if (!lifetime.active()) {
                continue;
            }
            additionPixels_[additionCursor[static_cast<std::size_t>(lifetime.owner)]++] = pixel;
            if (lifetime.stopExclusive != InvalidNode) {
                removalPixels_[removalCursor[static_cast<std::size_t>(lifetime.stopExclusive)]++] = pixel;
            }
        }
        tree_.requireMutationVersion(mutationVersion_, "MorphologicalTreeBoundaryLifetimeIndex construction");
    }

    /** @brief Returns the captured lifetime of a valid domain pixel. */
    [[nodiscard]] const BoundaryLifetime& lifetime(PixelId pixel) const {
        requireStableTree();
        if (!tree_.isPixel(pixel)) {
            throw std::out_of_range("Boundary-lifetime index received an invalid pixel id.");
        }
        return lifetimes_[static_cast<std::size_t>(pixel)];
    }

    /** @brief Returns pixels whose contour lifetime starts at one live node. */
    [[nodiscard]] std::span<const PixelId> additions(NodeId node) const {
        requireLiveNode(node);
        return eventSpan(additionPixels_, additionOffsets_, node);
    }

    /** @brief Returns pixels whose contour lifetime stops before one live node. */
    [[nodiscard]] std::span<const PixelId> removals(NodeId node) const {
        requireLiveNode(node);
        return eventSpan(removalPixels_, removalOffsets_, node);
    }

    /** @brief Tests the exact contour predicate for one pixel and live node. */
    [[nodiscard]] bool isBoundaryAt(PixelId pixel, NodeId node) const {
        requireLiveNode(node);
        const BoundaryLifetime& value = lifetime(pixel);
        if (!isAncestorInCapturedHierarchy(node, value.owner)) {
            return false;
        }
        if (value.stopExclusive == InvalidNode) {
            return true;
        }
        return value.stopExclusive != node && isAncestorInCapturedHierarchy(value.stopExclusive, node);
    }

    /** @brief Returns the number of active-lifetime start events. */
    [[nodiscard]] std::size_t numAdditions() const {
        requireStableTree();
        return additionPixels_.size();
    }

    /** @brief Returns the number of finite lifetime-stop events. */
    [[nodiscard]] std::size_t numRemovals() const {
        requireStableTree();
        return removalPixels_.size();
    }

  private:
    /** @brief One pairwise LCA request whose result completes a pixel lifetime. */
    struct OfflineLcaQuery {
        NodeId first = InvalidNode;
        NodeId second = InvalidNode;
        PixelId pixel = InvalidPixel;
    };

    enum class OfflineTraversalEventKind : std::uint8_t { Enter, MergeChild, Exit };

    /** @brief Explicit event used by the iterative Tarjan traversal. */
    struct OfflineTraversalEvent {
        OfflineTraversalEventKind kind = OfflineTraversalEventKind::Enter;
        NodeId node = InvalidNode;
        NodeId child = InvalidNode;
    };

    [[nodiscard]] static GridDomain2D requireDomain(const MorphologicalTree& tree) {
        tree.requireNotEditing("MorphologicalTreeBoundaryLifetimeIndex");
        return tree.requireGridDomain2D("MorphologicalTreeBoundaryLifetimeIndex");
    }

    /**
     * @brief Assigns a dense preorder rank to every live node without recursion.
     *
     * For any non-empty node set, the LCA of the complete set equals the LCA of
     * its first and last nodes in depth-first preorder. This reduces each
     * five-owner contour query to one offline pairwise query.
     */
    static void buildPreOrderIntervals(const MorphologicalTree& tree, std::vector<int>& entries, std::vector<int>& subtreeEnds) {
        int nextIndex = 0;
        std::vector<std::pair<NodeId, bool>> traversalStack;
        traversalStack.reserve(static_cast<std::size_t>(tree.numNodes()) + 1);
        traversalStack.emplace_back(tree.root(), false);
        while (!traversalStack.empty()) {
            const auto [node, exiting] = traversalStack.back();
            traversalStack.pop_back();
            if (!tree.isAlive(node)) {
                throw std::logic_error("Boundary-lifetime indexing found an invalid live-node hierarchy.");
            }
            const std::size_t nodeIndex = static_cast<std::size_t>(node);
            if (exiting) {
                subtreeEnds[nodeIndex] = nextIndex;
                continue;
            }
            if (entries[nodeIndex] >= 0) {
                throw std::logic_error("Boundary-lifetime indexing found a cyclic live-node hierarchy.");
            }
            entries[nodeIndex] = nextIndex++;
            traversalStack.emplace_back(node, true);
            for (NodeId child : tree.children(node)) {
                traversalStack.emplace_back(child, false);
            }
        }
        if (nextIndex != tree.numNodes()) {
            throw std::logic_error("Boundary-lifetime indexing found a disconnected live-node hierarchy.");
        }
    }

    /** @brief Tests ancestry using the index-owned iterative DFS intervals. */
    [[nodiscard]] bool isAncestorInCapturedHierarchy(NodeId ancestor, NodeId descendant) const noexcept {
        const int ancestorEntry = preOrderIndices_[static_cast<std::size_t>(ancestor)];
        const int descendantEntry = preOrderIndices_[static_cast<std::size_t>(descendant)];
        return ancestorEntry <= descendantEntry && descendantEntry < subtreeEndIndices_[static_cast<std::size_t>(ancestor)];
    }

    /** @brief Returns the representative of one initialized disjoint set. */
    [[nodiscard]] static NodeId findSet(std::vector<NodeId>& setParents, NodeId node) {
        NodeId root = node;
        while (setParents[static_cast<std::size_t>(root)] != root) {
            root = setParents[static_cast<std::size_t>(root)];
        }
        while (setParents[static_cast<std::size_t>(node)] != node) {
            const NodeId next = setParents[static_cast<std::size_t>(node)];
            setParents[static_cast<std::size_t>(node)] = root;
            node = next;
        }
        return root;
    }

    /** @brief Unites two initialized sets by rank and returns their representative. */
    [[nodiscard]] static NodeId uniteSets(std::vector<NodeId>& setParents, std::vector<std::uint8_t>& setRanks, NodeId lhs, NodeId rhs) {
        NodeId lhsRoot = findSet(setParents, lhs);
        NodeId rhsRoot = findSet(setParents, rhs);
        if (lhsRoot == rhsRoot) {
            return lhsRoot;
        }
        if (setRanks[static_cast<std::size_t>(lhsRoot)] < setRanks[static_cast<std::size_t>(rhsRoot)]) {
            std::swap(lhsRoot, rhsRoot);
        }
        setParents[static_cast<std::size_t>(rhsRoot)] = lhsRoot;
        if (setRanks[static_cast<std::size_t>(lhsRoot)] == setRanks[static_cast<std::size_t>(rhsRoot)]) {
            ++setRanks[static_cast<std::size_t>(lhsRoot)];
        }
        return lhsRoot;
    }

    /**
     * @brief Resolves all requested LCAs by iterative offline Tarjan traversal.
     *
     * Query references use CSR-like storage. The traversal never calls the
     * tree's lazy LCA service, so construction uses O(P + N) temporary memory,
     * performs no recursive descent, and leaves persistent tree caches untouched.
     */
    static void resolveOfflineLcas(const MorphologicalTree& tree, std::span<const OfflineLcaQuery> queries, std::vector<BoundaryLifetime>& lifetimes) {
        if (queries.empty()) {
            return;
        }

        const std::size_t numSlots = static_cast<std::size_t>(tree.numInternalNodeSlots());
        std::vector<std::size_t> queryOffsets(numSlots + 1, 0);
        for (const OfflineLcaQuery& query : queries) {
            if (!tree.isAlive(query.first) || !tree.isAlive(query.second) || !tree.isPixel(query.pixel)) {
                throw std::logic_error("Boundary-lifetime indexing produced an invalid offline LCA query.");
            }
            ++queryOffsets[static_cast<std::size_t>(query.first) + 1];
            ++queryOffsets[static_cast<std::size_t>(query.second) + 1];
        }
        prefixSum(queryOffsets);

        std::vector<NodeId> otherEndpoints(queryOffsets.back(), InvalidNode);
        std::vector<std::size_t> queryIndices(queryOffsets.back(), 0);
        std::vector<std::size_t> queryCursors = queryOffsets;
        for (std::size_t queryIndex = 0; queryIndex < queries.size(); ++queryIndex) {
            const OfflineLcaQuery& query = queries[queryIndex];
            const std::size_t firstPosition = queryCursors[static_cast<std::size_t>(query.first)]++;
            otherEndpoints[firstPosition] = query.second;
            queryIndices[firstPosition] = queryIndex;
            const std::size_t secondPosition = queryCursors[static_cast<std::size_t>(query.second)]++;
            otherEndpoints[secondPosition] = query.first;
            queryIndices[secondPosition] = queryIndex;
        }

        std::vector<NodeId> setParents(numSlots, InvalidNode);
        std::vector<std::uint8_t> setRanks(numSlots, std::uint8_t{0});
        std::vector<NodeId> ancestors(numSlots, InvalidNode);
        std::vector<std::uint8_t> finished(numSlots, std::uint8_t{0});
        std::vector<NodeId> resolvedQueries(queries.size(), InvalidNode);
        std::vector<OfflineTraversalEvent> traversalStack;
        traversalStack.reserve(static_cast<std::size_t>(tree.numNodes()) * 2 + 1);
        traversalStack.push_back({OfflineTraversalEventKind::Enter, tree.root(), InvalidNode});

        while (!traversalStack.empty()) {
            const OfflineTraversalEvent event = traversalStack.back();
            traversalStack.pop_back();
            const std::size_t nodeIndex = static_cast<std::size_t>(event.node);
            if (event.kind == OfflineTraversalEventKind::Enter) {
                setParents[nodeIndex] = event.node;
                ancestors[nodeIndex] = event.node;
                traversalStack.push_back({OfflineTraversalEventKind::Exit, event.node, InvalidNode});
                for (NodeId child : tree.children(event.node)) {
                    traversalStack.push_back({OfflineTraversalEventKind::MergeChild, event.node, child});
                    traversalStack.push_back({OfflineTraversalEventKind::Enter, child, InvalidNode});
                }
            } else if (event.kind == OfflineTraversalEventKind::MergeChild) {
                const NodeId representative = uniteSets(setParents, setRanks, event.node, event.child);
                ancestors[static_cast<std::size_t>(representative)] = event.node;
            } else {
                finished[nodeIndex] = 1;
                for (std::size_t reference = queryOffsets[nodeIndex]; reference < queryOffsets[nodeIndex + 1]; ++reference) {
                    const NodeId other = otherEndpoints[reference];
                    if (finished[static_cast<std::size_t>(other)] == 0) {
                        continue;
                    }
                    const NodeId representative = findSet(setParents, other);
                    resolvedQueries[queryIndices[reference]] = ancestors[static_cast<std::size_t>(representative)];
                }
            }
        }

        for (std::size_t queryIndex = 0; queryIndex < queries.size(); ++queryIndex) {
            const NodeId lca = resolvedQueries[queryIndex];
            if (!tree.isAlive(lca)) {
                throw std::logic_error("Boundary-lifetime indexing failed to resolve an offline LCA query.");
            }
            lifetimes[static_cast<std::size_t>(queries[queryIndex].pixel)].stopExclusive = lca;
        }
    }

    static void prefixSum(std::vector<std::size_t>& offsets) {
        for (std::size_t index = 1; index < offsets.size(); ++index) {
            offsets[index] += offsets[index - 1];
        }
    }

    [[nodiscard]] static std::span<const PixelId> eventSpan(const std::vector<PixelId>& values, const std::vector<std::size_t>& offsets, NodeId node) {
        const std::size_t nodeIndex = static_cast<std::size_t>(node);
        return std::span<const PixelId>(values).subspan(offsets[nodeIndex], offsets[nodeIndex + 1] - offsets[nodeIndex]);
    }

    void requireStableTree() const { tree_.requireMutationVersion(mutationVersion_, "MorphologicalTreeBoundaryLifetimeIndex"); }

    void requireLiveNode(NodeId node) const {
        requireStableTree();
        if (!tree_.isAlive(node)) {
            throw std::out_of_range("Boundary-lifetime index received a non-live node id.");
        }
    }

    const MorphologicalTree& tree_;
    std::size_t mutationVersion_ = 0;
    GridDomain2D domain_{};
    std::vector<BoundaryLifetime> lifetimes_;
    std::vector<int> preOrderIndices_;
    std::vector<int> subtreeEndIndices_;
    std::vector<std::size_t> additionOffsets_;
    std::vector<std::size_t> removalOffsets_;
    std::vector<PixelId> additionPixels_;
    std::vector<PixelId> removalPixels_;
};

} // namespace mmcfilters::contours::detail
