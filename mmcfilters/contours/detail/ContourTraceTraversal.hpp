#pragma once

#include "ContourBoundaryTracer.hpp"
#include "ContourEdgeDeltaStore.hpp"
#include "../../trees/MorphologicalTree.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcfilters::contours::detail {

/**
 * @brief Resumable post-order traversal of ordered node contour traces.
 *
 * Child edge sets are retained only until their parent is built. The parent
 * takes the largest child buffer and releases the others after merging them.
 * The current trace remains valid until the next call to `advance()`.
 */
class ContourTraceTraversal {
  public:
    /**
     * @brief Prepares an independent traversal over shared immutable indexes.
     * @param tree Stable source tree.
     * @param edgeDeltas Compact edge additions and removals by node.
     * @param connectivityByNode Foreground connectivity for each node slot.
     */
    ContourTraceTraversal(const MorphologicalTree& tree, const ContourEdgeDeltaStore& edgeDeltas,
                          std::span<const ForegroundConnectivity> connectivityByNode)
        : tree_(tree), edgeDeltas_(edgeDeltas), connectivityByNode_(connectivityByNode), mutationVersion_(tree.getMutationVersion()),
          edgesByNode_(static_cast<std::size_t>(tree.numInternalNodeSlots())),
          edgeMarks_(static_cast<std::size_t>(4 * tree.numPixels()), 0), tracer_(tree.numRows(), tree.numColumns()) {
        postOrder_.reserve(static_cast<std::size_t>(tree.numNodes()));
        for (NodeId node : tree.postOrder()) {
            postOrder_.push_back(node);
        }
    }

    /**
     * @brief Advances to the next ordered node trace.
     * @return True when a new current trace is available.
     */
    bool advance() {
        requireStableTree();
        currentNode_ = InvalidNode;
        currentBoundaries_.clear();
        if (nextNodeIndex_ == postOrder_.size()) {
            edgesByNode_.clear();
            return false;
        }

        const NodeId node = postOrder_[nextNodeIndex_++];
        const NodeId reusableChild = childWithLargestEdgeSet(node);
        std::size_t requiredCapacity = edgeDeltas_.additions(node).size();
        for (NodeId child : tree_.children(node)) {
            const auto& childEdges = edgesByNode_[static_cast<std::size_t>(child)];
            if (!childEdges) {
                throw std::logic_error("Contour trace traversal found a missing child edge set.");
            }
            requiredCapacity += childEdges->size();
        }

        std::unique_ptr<std::vector<int>> nodeEdges;
        if (reusableChild == InvalidNode) {
            nodeEdges = std::make_unique<std::vector<int>>();
        } else {
            nodeEdges = std::move(edgesByNode_[static_cast<std::size_t>(reusableChild)]);
        }
        if (requiredCapacity > nodeEdges->capacity()) {
            nodeEdges->reserve(std::max(requiredCapacity, nodeEdges->capacity() * 2));
        }

        advanceMarkGeneration();
        for (int packedEdge : *nodeEdges) {
            markPresent(packedEdge);
        }
        for (NodeId child : tree_.children(node)) {
            if (child == reusableChild) {
                continue;
            }
            auto& childEdges = edgesByNode_[static_cast<std::size_t>(child)];
            for (int packedEdge : *childEdges) {
                appendIfAbsent(*nodeEdges, packedEdge);
            }
            childEdges.reset();
        }
        for (int packedEdge : edgeDeltas_.additions(node)) {
            appendIfAbsent(*nodeEdges, packedEdge);
        }
        for (int packedEdge : edgeDeltas_.removals(node)) {
            markAbsent(packedEdge);
        }

        std::size_t nextKeptEdge = 0;
        for (int packedEdge : *nodeEdges) {
            if (isPresent(packedEdge)) {
                (*nodeEdges)[nextKeptEdge++] = packedEdge;
            }
        }
        nodeEdges->resize(nextKeptEdge);
        if (nodeEdges->empty()) {
            throw std::logic_error("Contour trace traversal produced an empty edge set for a live node.");
        }

        tracer_.trace(*nodeEdges, currentBoundaries_, connectivityByNode_[static_cast<std::size_t>(node)]);
        edgesByNode_[static_cast<std::size_t>(node)] = std::move(nodeEdges);
        currentNode_ = node;
        return true;
    }

    /**
     * @brief Borrows the current node identifier and ordered trace.
     * @return Current node and trace view.
     */
    [[nodiscard]] std::pair<NodeId, ContourTraceView> current() const {
        requireStableTree();
        if (currentNode_ == InvalidNode) {
            throw std::out_of_range("Contour trace traversal has no current node.");
        }
        const auto& edges = edgesByNode_[static_cast<std::size_t>(currentNode_)];
        return {currentNode_, ContourTraceView(*edges, currentBoundaries_)};
    }

    /** @brief Rejects access after a topology mutation. */
    void requireStableTree() const { tree_.requireMutationVersion(mutationVersion_, "ContourTraceTraversal"); }

  private:
    /**
     * @brief Selects the child whose edge-vector capacity should be reused.
     * @param node Parent node.
     * @return Child with the largest edge set, or `InvalidNode` for a leaf.
     */
    [[nodiscard]] NodeId childWithLargestEdgeSet(NodeId node) const {
        NodeId largestChild = InvalidNode;
        std::size_t largestSize = 0;
        for (NodeId child : tree_.children(node)) {
            const auto& childEdges = edgesByNode_[static_cast<std::size_t>(child)];
            if (childEdges && (largestChild == InvalidNode || childEdges->size() > largestSize)) {
                largestChild = child;
                largestSize = childEdges->size();
            }
        }
        return largestChild;
    }

    /** @brief Starts a fresh membership-mark generation. */
    void advanceMarkGeneration() {
        ++edgeMarkGeneration_;
        if (edgeMarkGeneration_ == 0) {
            std::fill(edgeMarks_.begin(), edgeMarks_.end(), 0);
            edgeMarkGeneration_ = 1;
        }
    }

    /** @brief Marks one valid packed edge as present. @param packedEdge Packed contour edge. */
    void markPresent(int packedEdge) {
        if (packedEdge >= 0 && packedEdge < static_cast<int>(edgeMarks_.size())) {
            edgeMarks_[static_cast<std::size_t>(packedEdge)] = edgeMarkGeneration_;
        }
    }

    /** @brief Marks one valid packed edge as absent. @param packedEdge Packed contour edge. */
    void markAbsent(int packedEdge) {
        if (packedEdge >= 0 && packedEdge < static_cast<int>(edgeMarks_.size())) {
            edgeMarks_[static_cast<std::size_t>(packedEdge)] = 0;
        }
    }

    /**
     * @brief Tests whether one packed edge is present in the active generation.
     * @param packedEdge Packed contour edge.
     * @return True when the edge is present.
     */
    [[nodiscard]] bool isPresent(int packedEdge) const {
        return packedEdge >= 0 && packedEdge < static_cast<int>(edgeMarks_.size()) &&
               edgeMarks_[static_cast<std::size_t>(packedEdge)] == edgeMarkGeneration_;
    }

    /** @brief Appends one packed edge when it is absent. @param edges Destination edge set. @param packedEdge Packed contour edge. */
    void appendIfAbsent(std::vector<int>& edges, int packedEdge) {
        if (!isPresent(packedEdge)) {
            markPresent(packedEdge);
            edges.push_back(packedEdge);
        }
    }

    const MorphologicalTree& tree_;                         ///< Stable source topology.
    const ContourEdgeDeltaStore& edgeDeltas_;               ///< Shared immutable edge changes.
    std::span<const ForegroundConnectivity> connectivityByNode_; ///< Connectivity for each node slot.
    std::size_t mutationVersion_ = 0;                       ///< Captured tree mutation version.
    std::vector<std::unique_ptr<std::vector<int>>> edgesByNode_; ///< Active child and current edge sets.
    std::vector<uint16_t> edgeMarks_;                       ///< Packed-edge membership generations.
    uint16_t edgeMarkGeneration_ = 1;                       ///< Active membership generation.
    ContourBoundaryTracer tracer_;                          ///< Reusable ordered-boundary tracer.
    std::vector<ContourBoundary> currentBoundaries_;        ///< Boundaries of the current node.
    std::vector<NodeId> postOrder_;                         ///< Child-before-parent traversal schedule.
    std::size_t nextNodeIndex_ = 0;                         ///< Next schedule position.
    NodeId currentNode_ = InvalidNode;                      ///< Current node identifier.
};

} // namespace mmcfilters::contours::detail
