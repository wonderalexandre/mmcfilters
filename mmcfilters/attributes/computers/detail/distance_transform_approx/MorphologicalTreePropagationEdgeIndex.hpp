#pragma once

#include "../../../../trees/MorphologicalTree.hpp"
#include "../../../../trees/detail/CommittedTreeAccess.hpp"
#include "../../../../utils/Common.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcfilters::attributes::computers::detail::distance_transform_approx {

/**
 * @brief One undirected A8 pixel edge activated when its owner LCA is processed.
 */
struct PropagationEdge2D {
    PixelId first = InvalidPixel;
    PixelId second = InvalidPixel;
    std::uint8_t firstDirectionBit = 0;
    std::uint8_t secondDirectionBit = 0;
};

/**
 * @brief Topology-derived activation events for all undirected A8 pixel edges.
 *
 * An edge `(p,q)` becomes valid exactly at the LCA of the inclusion-smallest
 * nodes of its endpoints. This
 * lets several incomparable regions coexist in one DIFT workspace without
 * allowing propagation between siblings before their parent is processed.
 * LCA ownership is resolved through the tree's shared Euler/RMQ cache. The
 * cache is topology-generic, is reused by subsequent attribute computations,
 * and is invalidated automatically whenever the tree topology changes.
 */
class MorphologicalTreePropagationEdgeIndex {
  public:
    MorphologicalTreePropagationEdgeIndex(const MorphologicalTreePropagationEdgeIndex&) = delete;
    MorphologicalTreePropagationEdgeIndex& operator=(const MorphologicalTreePropagationEdgeIndex&) = delete;
    MorphologicalTreePropagationEdgeIndex(MorphologicalTreePropagationEdgeIndex&&) = delete;
    MorphologicalTreePropagationEdgeIndex& operator=(MorphologicalTreePropagationEdgeIndex&&) = delete;

    explicit MorphologicalTreePropagationEdgeIndex(const MorphologicalTree& tree)
        : tree_(tree), mutationVersion_(tree.getMutationVersion()), domain_(requireDomain(tree)) {
        buildEvents();
    }

    [[nodiscard]] std::span<const PropagationEdge2D> activations(NodeId node) const {
        requireLiveNode(node);
        const std::size_t index = static_cast<std::size_t>(node);
        return std::span<const PropagationEdge2D>(edges_).subspan(offsets_[index], offsets_[index + 1] - offsets_[index]);
    }

    [[nodiscard]] std::size_t numEdges() const {
        requireStableTree();
        return edges_.size();
    }

    /**
     * @brief Returns activations for a caller-established live node.
     */
    [[nodiscard]] std::span<const PropagationEdge2D> establishedActivations(NodeId node) const noexcept {
        const std::size_t index = static_cast<std::size_t>(node);
        return std::span<const PropagationEdge2D>(edges_).subspan(offsets_[index], offsets_[index + 1] - offsets_[index]);
    }

    /**
     * @brief Returns the edge count after the caller established tree stability.
     */
    [[nodiscard]] std::size_t establishedNumEdges() const noexcept { return edges_.size(); }

  private:
    struct PendingEdge {
        PropagationEdge2D edge;
        NodeId activationNode = InvalidNode;
    };

    [[nodiscard]] static GridDomain2D requireDomain(const MorphologicalTree& tree) {
        tree.requireNotEditing("MorphologicalTreePropagationEdgeIndex");
        return tree.requireGridDomain2D("MorphologicalTreePropagationEdgeIndex");
    }

    void buildEvents() {
        const std::span<const NodeId> owners = tree_.smallestNodeMap();
        constexpr std::array<std::pair<int, int>, 4> forwardOffsets{{{0, 1}, {1, -1}, {1, 0}, {1, 1}}};
        std::vector<PendingEdge> pending;
        pending.reserve(static_cast<std::size_t>(tree_.numPixels()) * 4);
        for (PixelId pixel = 0; pixel < tree_.numPixels(); ++pixel) {
            const int row = pixel / domain_.columns;
            const int column = pixel % domain_.columns;
            for (const auto [rowOffset, columnOffset] : forwardOffsets) {
                const int neighbourRow = row + rowOffset;
                const int neighbourColumn = column + columnOffset;
                if (neighbourRow < 0 || neighbourRow >= domain_.rows || neighbourColumn < 0 || neighbourColumn >= domain_.columns) {
                    continue;
                }
                const PixelId neighbour = neighbourRow * domain_.columns + neighbourColumn;
                const NodeId firstOwner = owners[static_cast<std::size_t>(pixel)];
                const NodeId secondOwner = owners[static_cast<std::size_t>(neighbour)];
                const NodeId activationNode =
                    ::mmcfilters::detail::CommittedTreeAccess::lowestCommonAncestor(tree_, firstOwner, secondOwner);
                pending.push_back({{pixel, neighbour, directionBit(rowOffset, columnOffset),
                                    directionBit(-rowOffset, -columnOffset)},
                                   activationNode});
            }
        }

        offsets_.assign(static_cast<std::size_t>(tree_.numInternalNodeSlots()) + 1, 0);
        for (const PendingEdge& value : pending) {
            ++offsets_[static_cast<std::size_t>(value.activationNode) + 1];
        }
        for (std::size_t index = 1; index < offsets_.size(); ++index) {
            offsets_[index] += offsets_[index - 1];
        }
        edges_.resize(pending.size());
        std::vector<std::size_t> cursors = offsets_;
        for (const PendingEdge& value : pending) {
            edges_[cursors[static_cast<std::size_t>(value.activationNode)]++] = value.edge;
        }
    }

    [[nodiscard]] static constexpr std::uint8_t directionBit(int rowOffset, int columnOffset) noexcept {
        constexpr std::array<std::uint8_t, 9> bits{
            std::uint8_t{1}, std::uint8_t{2}, std::uint8_t{4}, std::uint8_t{8}, std::uint8_t{0},
            std::uint8_t{16}, std::uint8_t{32}, std::uint8_t{64}, std::uint8_t{128}};
        return bits[static_cast<std::size_t>((rowOffset + 1) * 3 + (columnOffset + 1))];
    }

    void requireStableTree() const { tree_.requireMutationVersion(mutationVersion_, "MorphologicalTreePropagationEdgeIndex"); }

    void requireLiveNode(NodeId node) const {
        requireStableTree();
        if (!tree_.isAlive(node)) {
            throw std::out_of_range("Propagation-edge index received a non-live node id.");
        }
    }

    const MorphologicalTree& tree_;
    std::size_t mutationVersion_ = 0;
    GridDomain2D domain_;
    std::vector<std::size_t> offsets_;
    std::vector<PropagationEdge2D> edges_;
};

} // namespace mmcfilters::attributes::computers::detail::distance_transform_approx
