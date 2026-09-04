#pragma once

#include "../localAttributes/FiniteWindowLocalAttributeComputer.hpp"
#include "../trees/MorphologicalTree.hpp"
#include "../trees/ValuedMorphologicalTreeView.hpp"
#include "../trees/detail/MorphologicalTreeConstructionContextQueries.hpp"
#include "../trees/detail/TreeTraversalDetail.hpp"
#include "../utils/Common.hpp"
#include "../utils/Image.hpp"
#include "detail/ContourEdgeDeltaStore.hpp"
#include "detail/PendingEdgeLists.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mmcfilters {

/**
 * @brief Side of a support pixel occupied by a contour edge.
 */
enum class ContourSide : uint8_t { North = 0, West = 1, East = 2, South = 3 };

/**
 * @brief Whether a contour boundary is external or surrounds a hole.
 */
enum class ContourBoundaryKind : uint8_t { External, Internal };

/**
 * @brief One contour edge represented by its support pixel and side.
 */
struct ContourEdge {
    /// Row-major support-pixel index incident to the boundary edge.
    PixelId pixel = InvalidPixel;
    /// Side of the support pixel occupied by the boundary edge.
    ContourSide side = ContourSide::North;

    /// Compares the support pixel and side.
    friend bool operator==(const ContourEdge&, const ContourEdge&) = default;
};

/**
 * @brief One ordered external or internal contour boundary.
 */
struct ContourBoundary {
    /// Whether this boundary is an external boundary or an internal hole.
    ContourBoundaryKind kind = ContourBoundaryKind::External;
    /// First edge in the shared ordered-edge buffer.
    uint32_t edgeOffset = 0;
    /// Number of consecutive edges in the boundary.
    uint32_t edgeCount = 0;
    /// Doubled signed area under the contour orientation convention.
    int doubledSignedArea = 0;
};

/**
 * @brief Lazy contour-edge construction and ordered boundary tracing.
 *
 * Boundary primitives are oriented grid edges. Edges are cached lazily for each
 * node, then traced into ordered boundaries on demand. The orientation convention
 * keeps the support pixel on the right side of each directed edge in image
 * coordinates, where rows grow downward and columns grow rightward. With this
 * convention, external boundaries have positive doubled signed area and internal
 * boundaries have negative doubled signed area.
 */
class ContourTraceComputation {
  private:
    /** @brief Foreground connectivity, or Unknown when construction information is insufficient. */
    enum class ForegroundConnectivity : uint8_t { Unknown, Four, Eight };

    /// Compact edge additions and removals indexed by node.
    using EdgeDeltas = contours::detail::ContourEdgeDeltaStore;

    /** @brief Buffers prepared before lazy edge caching begins. */
    struct ConstructionData {
        /// Packed edge additions and removals for every node.
        EdgeDeltas edgeDeltas;
        /// Capacity estimate for the shared cached-edge buffer.
        int estimatedCachedEdgeCount = 0;
        /// Foreground connectivity used to resolve diagonal contacts at each node.
        std::vector<ForegroundConnectivity> connectivityByNode;
    };

  public:
    /**
     * @brief Packs one pixel-side edge into a compact integer id.
     *
     * @param pixel Pixel identifier.
     * @param side Side selected by the operation.
     * @return Packed contour edge identifier.
     */
    [[nodiscard]] static int packEdge(PixelId pixel, ContourSide side) { return (4 * pixel) + static_cast<int>(side); }

    /**
     * @brief Unpacks one compact edge id.
     *
     * @param packedEdge Packed boundary-edge identifier.
     * @return Support pixel and side represented by `packedEdge`.
     */
    [[nodiscard]] static ContourEdge unpackEdge(int packedEdge) {
        if (packedEdge < 0) {
            return {};
        }
        const int sideIndex = packedEdge & 3;
        return ContourEdge{packedEdge / 4, static_cast<ContourSide>(sideIndex)};
    }

  private:
    /** @brief Cardinal direction of an oriented contour edge. */
    enum class Direction : uint8_t { North = 0, East = 1, South = 2, West = 3 };

    /** @brief Storage strategy for the mapping from vertices to outgoing edges. */
    enum class VertexIndexMode : uint8_t { Dense, Sparse };

    /** @brief Stores one oriented contour edge and its accumulated geometry. */
    struct DirectedEdge {
        /// Compact identifier of the source pixel and side.
        int packedEdge = -1;
        /// Grid vertex where the directed edge begins.
        int startVertex = -1;
        /// Grid vertex where the directed edge ends.
        int endVertex = -1;
        /// Direction from `startVertex` to `endVertex`.
        Direction direction = Direction::North;
        /// Contribution of this edge to twice the signed boundary area.
        int doubledSignedArea = 0;

        /**
         * @brief Constructs a default `DirectedEdge`.
         */
        DirectedEdge() = default;

        /**
         * @brief Constructs `DirectedEdge` from the supplied inputs.
         *
         * @param packedEdge Packed boundary-edge identifier.
         * @param startVertex Starting vertex identifier.
         * @param endVertex Terminal vertex of the traced directed edge.
         * @param direction Direction code of the traced edge.
         * @param doubledSignedArea Twice the signed area accumulated for the boundary.
         */
        DirectedEdge(int packedEdge, int startVertex, int endVertex, Direction direction, int doubledSignedArea)
            : packedEdge(packedEdge), startVertex(startVertex), endVertex(endVertex), direction(direction), doubledSignedArea(doubledSignedArea) {}
    };

    const MorphologicalTree& tree_; ///< Stable source tree.
    std::vector<ForegroundConnectivity> connectivityByNode_; ///< Foreground connectivity for each internal node slot.
    std::size_t treeMutationVersion_ = 0; ///< Tree mutation version captured at construction.
    mutable EdgeDeltas edgeDeltas_; ///< Edge changes retained until all node edge sets are cached.

    mutable std::vector<int> cachedPackedEdges_; ///< Shared storage for cached contour edges.
    mutable std::vector<uint32_t> edgeOffsets_; ///< First cached edge for each internal node slot.
    mutable std::vector<uint32_t> edgeCounts_; ///< Number of cached edges for each internal node slot.
    mutable std::vector<uint8_t> edgesCached_; ///< Whether each internal node slot has a cached edge set.
    mutable std::size_t numNodesWithCachedEdges_ = 0; ///< Number of live nodes with cached edge sets.

    mutable std::vector<ContourBoundary> cachedBoundaries_; ///< Shared storage for cached boundary descriptors.
    mutable std::vector<uint32_t> boundaryOffsets_; ///< First cached boundary for each internal node slot.
    mutable std::vector<uint32_t> boundaryCounts_; ///< Number of cached boundaries for each internal node slot.
    mutable std::vector<uint8_t> boundariesCached_; ///< Whether each internal node slot has cached boundaries.
    mutable std::size_t numNodesWithCachedBoundaries_ = 0; ///< Number of live nodes with cached boundaries.

    mutable std::vector<uint16_t> edgeMarks_; ///< Generation marks used while composing node edge sets.
    mutable uint16_t edgeMarkGeneration_ = 1; ///< Active generation in `edgeMarks_`.
    mutable bool edgeScratchReleased_ = false; ///< Whether edge construction scratch storage was released.

    mutable std::vector<DirectedEdge> directedEdges_; ///< Geometry for the node edges currently being traced.
    mutable std::vector<int> outgoingEdgeHeads_; ///< Dense first-outgoing-edge index for each grid vertex.
    mutable std::vector<int> nextOutgoingEdges_; ///< Linked-list successor for each directed edge.
    mutable std::vector<int> touchedVertices_; ///< Dense vertex entries that must be reset after tracing.
    mutable std::vector<int> sparseVertexKeys_; ///< Grid-vertex keys in the sparse outgoing-edge table.
    mutable std::vector<int> sparseOutgoingEdgeHeads_; ///< First outgoing edge stored in each sparse table slot.
    mutable std::vector<uint32_t> sparseSlotGenerations_; ///< Active generation of each sparse table slot.
    mutable std::vector<int> touchedSparseSlots_; ///< Sparse slots populated by the current node trace.
    mutable uint32_t sparseGeneration_ = 1; ///< Active generation in the sparse vertex table.
    mutable std::size_t numIndividualNodeTraces_ = 0; ///< Individual traces performed before selecting dense indexing.
    mutable std::vector<uint32_t> visitedGenerations_; ///< Visit generation for each edge of the current node.
    mutable uint32_t visitGeneration_ = 1; ///< Active generation in `visitedGenerations_`.
    mutable std::vector<int> orderedNodeEdges_; ///< Ordered packed edges produced for the current node.
    mutable std::vector<ContourBoundary> nodeBoundaries_; ///< Boundary descriptors produced for the current node.
    mutable bool traceScratchReleased_ = false; ///< Whether boundary-tracing scratch storage was released.

    /**
     * @brief Takes ownership of the prepared buffers without copying them.
     * @param tree Stable source topology.
     * @param constructionData Edge changes and foreground connectivity.
     */
    ContourTraceComputation(const MorphologicalTree& tree, ConstructionData constructionData)
        : tree_(tree), connectivityByNode_(std::move(constructionData.connectivityByNode)), treeMutationVersion_(tree.getMutationVersion()),
          edgeDeltas_(std::move(constructionData.edgeDeltas)),
          edgeOffsets_(static_cast<std::size_t>(tree.numInternalNodeSlots()), 0),
          edgeCounts_(static_cast<std::size_t>(tree.numInternalNodeSlots()), 0),
          edgesCached_(static_cast<std::size_t>(tree.numInternalNodeSlots()), 0),
          boundaryOffsets_(static_cast<std::size_t>(tree.numInternalNodeSlots()), 0),
          boundaryCounts_(static_cast<std::size_t>(tree.numInternalNodeSlots()), 0),
          boundariesCached_(static_cast<std::size_t>(tree.numInternalNodeSlots()), 0),
          edgeMarks_(static_cast<std::size_t>(4 * tree.numRows() * tree.numColumns()), 0) {
        if (constructionData.estimatedCachedEdgeCount > 0) {
            cachedPackedEdges_.reserve(static_cast<std::size_t>(constructionData.estimatedCachedEdgeCount));
        }
    }

  public:
    /**
     * @brief Immutable range over unpacked boundary edges.
     */
    class EdgeRange {
      public:
        /**
         * @brief Forward iterator that unpacks boundary edges on dereference.
         */
        class iterator {
          public:
            /// Standard category for a multi-pass forward iterator.
            using iterator_category = std::forward_iterator_tag;
            /// Unpacked edge value yielded by dereference.
            using value_type = ContourEdge;
            /// Signed iterator-distance type.
            using difference_type = std::ptrdiff_t;
            /// No pointer type is exposed because dereference returns a value.
            using pointer = void;
            /// Value-returning reference type.
            using reference = value_type;

            /**
             * @brief Creates an empty iterator.
             */
            iterator() = default;

            /**
             * @brief Creates an iterator over a packed-edge buffer position.
             *
             * @param packedEdges Shared packed edge buffer.
             * @param index Zero-based index.
             */
            iterator(const std::vector<int>* packedEdges, std::size_t index) : packedEdges_(packedEdges), index_(index) {}

            /**
             * @brief Returns the unpacked edge at the current position.
             *
             * @return The unpacked edge at the current position.
             */
            value_type operator*() const { return ContourTraceComputation::unpackEdge((*packedEdges_)[index_]); }

            /**
             * @brief Advances to the next packed edge.
             *
             * @return Mutable reference to the updated object.
             */
            iterator& operator++() {
                ++index_;
                return *this;
            }

            /**
             * @brief Advances and returns the previous iterator position.
             *
             * @return Iterator position before the advancement.
             */
            iterator operator++(int) {
                iterator previous(*this);
                ++(*this);
                return previous;
            }

            /// Compares the backing buffer and position.
            friend bool operator==(const iterator& lhs, const iterator& rhs) {
                return lhs.packedEdges_ == rhs.packedEdges_ && lhs.index_ == rhs.index_;
            }

            /// Returns true when two iterator positions differ.
            friend bool operator!=(const iterator& lhs, const iterator& rhs) { return !(lhs == rhs); }

          private:
            /// Shared packed-edge buffer traversed by this iterator.
            const std::vector<int>* packedEdges_ = nullptr;
            /// Current zero-based position in `packedEdges_`.
            std::size_t index_ = 0;
        };

        /**
         * @brief Creates an empty edge range.
         */
        EdgeRange() = default;

        /**
         * @brief Creates a view over `size` packed edges starting at `offset`.
         *
         * @param packedEdges Shared packed edge buffer.
         * @param offset Offset into the underlying storage.
         * @param size Number of edges in the range.
         */
        EdgeRange(const std::vector<int>* packedEdges, std::size_t offset, std::size_t size)
            : packedEdges_(packedEdges), offset_(offset), size_(size) {}

        /**
         * @brief Returns an iterator to the first edge.
         *
         * @return An iterator to the first edge.
         */
        iterator begin() const { return iterator(packedEdges_, offset_); }

        /**
         * @brief Returns the exclusive end iterator.
         *
         * @return The exclusive end iterator.
         */
        iterator end() const { return iterator(packedEdges_, offset_ + size_); }

        /**
         * @brief Returns whether the range contains no edges.
         *
         * @return Whether the range contains no edges.
         */
        [[nodiscard]] bool empty() const noexcept { return size_ == 0; }

        /**
         * @brief Returns the number of edges in the range.
         *
         * @return The number of edges in the range.
         */
        [[nodiscard]] std::size_t size() const noexcept { return size_; }

      private:
        /// Shared packed-edge buffer viewed by this range.
        const std::vector<int>* packedEdges_ = nullptr;
        /// Position of the first edge in `packedEdges_`.
        std::size_t offset_ = 0;
        /// Number of edges in the range.
        std::size_t size_ = 0;
    };

    /**
     * @brief Returns the unordered contour edges of one node.
     *
     * @param node Node identifier.
     * @return Borrowed range over the node's contour edges.
     */
    [[nodiscard]] EdgeRange edges(NodeId node) const {
        requireStableTree("ContourTraceComputation::edges");
        requireLiveTraceNode(node, "ContourTraceComputation::edges");
        cacheSubtreeEdges(node);
        return EdgeRange(&cachedPackedEdges_, edgeOffsets_[static_cast<std::size_t>(node)], edgeCounts_[static_cast<std::size_t>(node)]);
    }

    /**
     * @brief Returns an owning copy of the boundary metadata for one node.
     *
     * The returned vector remains valid when later lazy queries cache
     * boundaries for other nodes.
     *
     * @param node Node identifier.
     * @return An owning copy of the boundary metadata for one node.
     */
    [[nodiscard]] std::vector<ContourBoundary> boundaries(NodeId node) const {
        requireStableTree("ContourTraceComputation::boundaries");
        requireLiveTraceNode(node, "ContourTraceComputation::boundaries");
        ensureNodeBoundariesTraced(node);
        const auto nodeBoundaryOffset = static_cast<std::size_t>(boundaryOffsets_[static_cast<std::size_t>(node)]);
        const auto numNodeBoundaries = static_cast<std::size_t>(boundaryCounts_[static_cast<std::size_t>(node)]);
        if (numNodeBoundaries == 0) {
            return {};
        }
        return std::vector<ContourBoundary>(cachedBoundaries_.begin() + static_cast<std::ptrdiff_t>(nodeBoundaryOffset),
                                             cachedBoundaries_.begin() + static_cast<std::ptrdiff_t>(nodeBoundaryOffset + numNodeBoundaries));
    }

    /**
     * @brief Returns the ordered edges belonging to one boundary.
     *
     * @param boundary Contour boundary descriptor.
     * @return The ordered edges belonging to one boundary.
     */
    [[nodiscard]] EdgeRange boundaryEdges(const ContourBoundary& boundary) const {
        requireStableTree("ContourTraceComputation::boundaryEdges");
        const std::size_t edgeOffset = boundary.edgeOffset;
        const std::size_t numEdges = boundary.edgeCount;
        if (edgeOffset > cachedPackedEdges_.size() || numEdges > cachedPackedEdges_.size() - edgeOffset) {
            throw std::invalid_argument("ContourBoundary does not belong to this ContourTraceComputation.");
        }
        return EdgeRange(&cachedPackedEdges_, edgeOffset, numEdges);
    }

    /**
     * @brief Caches edges and traces ordered boundaries for every live node.
     */
    void traceAll() const {
        requireStableTree("ContourTraceComputation::traceAll");
        traceSubtree(tree_.root());
    }

    /**
     * @brief Tests whether every live node has cached ordered boundaries.
     *
     * @return True after every live node has been traced.
     */
    [[nodiscard]] bool hasTracedAllBoundaries() const {
        requireStableTree("ContourTraceComputation::hasTracedAllBoundaries");
        for (NodeId node : tree_.aliveNodeIds()) {
            if (!boundariesCached_[static_cast<std::size_t>(node)]) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Tests whether one node has cached contour edges.
     * @param node Live node identifier.
     * @return True when `edges(node)` has cached the edge set.
     */
    [[nodiscard]] bool hasCachedEdges(NodeId node) const {
        requireStableTree("ContourTraceComputation::hasCachedEdges");
        requireLiveTraceNode(node, "ContourTraceComputation::hasCachedEdges");
        return static_cast<bool>(edgesCached_[static_cast<std::size_t>(node)]);
    }

    /**
     * @brief Tests whether one node has cached ordered boundaries.
     * @param node Live node identifier.
     * @return True when the node's boundaries have been traced.
     */
    [[nodiscard]] bool hasTracedBoundaries(NodeId node) const {
        requireStableTree("ContourTraceComputation::hasTracedBoundaries");
        requireLiveTraceNode(node, "ContourTraceComputation::hasTracedBoundaries");
        return static_cast<bool>(boundariesCached_[static_cast<std::size_t>(node)]);
    }

  private:
    /**
     * @brief Rejects access after a topology mutation.
     * @param context Operation name used in diagnostics.
     */
    void requireStableTree(const char* context) const { tree_.requireMutationVersion(treeMutationVersion_, context); }

    /**
     * @brief Rejects a node outside the live internal node domain.
     * @param node Candidate node identifier.
     * @param context Operation name used in diagnostics.
     */
    void requireLiveTraceNode(NodeId node, const char* context) const {
        if (!tree_.isAlive(node)) {
            throw std::invalid_argument(std::string(context) + " requires a live internal NodeId.");
        }
    }

    /**
     * @brief Converts a checked buffer index to `uint32_t`.
     * @param value Buffer index or count.
     * @param context Operation name used in diagnostics.
     * @return `value` represented as `uint32_t`.
     */
    static uint32_t checkedUint32(std::size_t value, const char* context) {
        if (value > static_cast<std::size_t>(std::numeric_limits<uint32_t>::max())) {
            throw std::overflow_error(std::string(context) + " exceeds uint32_t.");
        }
        return static_cast<uint32_t>(value);
    }

    /**
     * @brief Clears a vector and releases its capacity.
     * @param values Vector whose storage is released.
     */
    template <class T> static void releaseVector(std::vector<T>& values) { std::vector<T>().swap(values); }

    /**
     * @brief Tests whether every live node has a cached edge set.
     * @return True when edge construction scratch storage can be released.
     */
    [[nodiscard]] bool allEdgesCached() const { return numNodesWithCachedEdges_ == static_cast<std::size_t>(tree_.numNodes()); }

    /**
     * @brief Tests whether every live node has cached ordered boundaries.
     * @return True when boundary tracing scratch storage can be released.
     */
    [[nodiscard]] bool allBoundariesCached() const { return numNodesWithCachedBoundaries_ == static_cast<std::size_t>(tree_.numNodes()); }

    /**
     * @brief Releases edge construction buffers after all edge sets are cached.
     */
    void releaseEdgeScratchIfComplete() const {
        if (edgeScratchReleased_ || !allEdgesCached()) {
            return;
        }

        edgeDeltas_ = EdgeDeltas{};
        releaseVector(edgeMarks_);
        edgeMarkGeneration_ = 1;
        edgeScratchReleased_ = true;
    }

    /**
     * @brief Releases the sparse vertex index.
     */
    void releaseSparseTraceOutgoingScratch() const {
        releaseVector(sparseVertexKeys_);
        releaseVector(sparseOutgoingEdgeHeads_);
        releaseVector(sparseSlotGenerations_);
        releaseVector(touchedSparseSlots_);
        sparseGeneration_ = 1;
    }

    /**
     * @brief Releases boundary tracing buffers after all nodes are traced.
     */
    void releaseTraceScratchIfComplete() const {
        if (traceScratchReleased_ || !allBoundariesCached()) {
            return;
        }

        resetDenseVertexIndex();
        releaseVector(directedEdges_);
        releaseVector(outgoingEdgeHeads_);
        releaseVector(nextOutgoingEdges_);
        releaseVector(touchedVertices_);
        releaseSparseTraceOutgoingScratch();
        releaseVector(visitedGenerations_);
        releaseVector(orderedNodeEdges_);
        releaseVector(nodeBoundaries_);
        numIndividualNodeTraces_ = 0;
        visitGeneration_ = 1;
        traceScratchReleased_ = true;
    }

    /**
     * @brief Starts a new edge mark generation.
     */
    void advanceEdgeMarkGeneration() const {
        ++edgeMarkGeneration_;
        if (edgeMarkGeneration_ == 0) {
            std::fill(edgeMarks_.begin(), edgeMarks_.end(), 0);
            edgeMarkGeneration_ = 1;
        }
    }

    /**
     * @brief Appends an edge that is absent from the current set.
     * @param edges Current node edge set.
     * @param packedEdge Packed boundary-edge identifier.
     */
    void appendIfAbsent(std::vector<int>& edges, int packedEdge) const {
        if (packedEdge < 0 || packedEdge >= static_cast<int>(edgeMarks_.size())) {
            return;
        }
        if (edgeMarks_[static_cast<std::size_t>(packedEdge)] != edgeMarkGeneration_) {
            edgeMarks_[static_cast<std::size_t>(packedEdge)] = edgeMarkGeneration_;
            edges.push_back(packedEdge);
        }
    }

    /**
     * @brief Removes an edge from the current set when present.
     * @param packedEdge Packed boundary-edge identifier.
     */
    void removeIfPresent(int packedEdge) const {
        if (packedEdge >= 0 && packedEdge < static_cast<int>(edgeMarks_.size())) {
            edgeMarks_[static_cast<std::size_t>(packedEdge)] = 0;
        }
    }

    /**
     * @brief Returns an iterator to the first cached edge of a node.
     *
     * @param node Node identifier.
     * @return Iterator to the first packed edge.
     */
    std::vector<int>::const_iterator cachedEdgeBegin(NodeId node) const {
        return cachedPackedEdges_.begin() + static_cast<std::ptrdiff_t>(edgeOffsets_[static_cast<std::size_t>(node)]);
    }

    /**
     * @brief Returns an iterator past the last cached edge of a node.
     *
     * @param node Node identifier.
     * @return Iterator past the last packed edge.
     */
    std::vector<int>::const_iterator cachedEdgeEnd(NodeId node) const {
        return cachedEdgeBegin(node) + static_cast<std::ptrdiff_t>(edgeCounts_[static_cast<std::size_t>(node)]);
    }

    /**
     * @brief Appends one node's contour edges to the shared cache.
     * @param node Live node identifier.
     * @param edges Packed contour edges for `node`.
     */
    void cacheNodeEdges(NodeId node, const std::vector<int>& edges) const {
        edgeOffsets_[static_cast<std::size_t>(node)] = checkedUint32(cachedPackedEdges_.size(), "cached trace edge offset");
        edgeCounts_[static_cast<std::size_t>(node)] = checkedUint32(edges.size(), "cached trace edge count");
        cachedPackedEdges_.insert(cachedPackedEdges_.end(), edges.begin(), edges.end());
        const auto nodeIndex = static_cast<std::size_t>(node);
        if (!edgesCached_[nodeIndex]) {
            edgesCached_[nodeIndex] = 1;
            ++numNodesWithCachedEdges_;
        }
    }

    /**
     * @brief Materializes every uncached edge set in a subtree.
     * @param subtreeRoot Root of the requested subtree.
     */
    void cacheSubtreeEdges(NodeId subtreeRoot) const {
        requireStableTree("ContourTraceComputation::cacheSubtreeEdges");
        requireLiveTraceNode(subtreeRoot, "ContourTraceComputation::cacheSubtreeEdges");
        if (edgesCached_[static_cast<std::size_t>(subtreeRoot)]) {
            return;
        }

        std::vector<std::pair<NodeId, bool>> traversalStack;
        traversalStack.emplace_back(subtreeRoot, false);
        std::vector<int> nodeEdges;

        while (!traversalStack.empty()) {
            const auto [node, expanded] = traversalStack.back();
            traversalStack.pop_back();
            if (edgesCached_[static_cast<std::size_t>(node)]) {
                continue;
            }
            if (!expanded) {
                traversalStack.emplace_back(node, true);
                for (NodeId child : tree_.children(node)) {
                    if (!edgesCached_[static_cast<std::size_t>(child)]) {
                        traversalStack.emplace_back(child, false);
                    }
                }
                continue;
            }

            nodeEdges.clear();
            const auto additions = edgeDeltas_.additions(node);
            std::size_t requiredCapacity = additions.size();
            for (NodeId child : tree_.children(node)) {
                requiredCapacity += static_cast<std::size_t>(edgeCounts_[static_cast<std::size_t>(child)]);
            }
            nodeEdges.reserve(requiredCapacity);
            advanceEdgeMarkGeneration();

            for (NodeId child : tree_.children(node)) {
                for (auto edge = cachedEdgeBegin(child); edge != cachedEdgeEnd(child); ++edge) {
                    appendIfAbsent(nodeEdges, *edge);
                }
            }

            for (int packedEdge : additions) {
                appendIfAbsent(nodeEdges, packedEdge);
            }

            for (int packedEdge : edgeDeltas_.removals(node)) {
                removeIfPresent(packedEdge);
            }

            std::size_t nextKeptEdge = 0;
            for (int packedEdge : nodeEdges) {
                if (edgeMarks_[static_cast<std::size_t>(packedEdge)] == edgeMarkGeneration_) {
                    nodeEdges[nextKeptEdge++] = packedEdge;
                }
            }
            nodeEdges.resize(nextKeptEdge);
            cacheNodeEdges(node, nodeEdges);
        }
        releaseEdgeScratchIfComplete();
    }

    /**
     * @brief Computes the row-major identifier of a grid vertex.
     *
     * @param row Zero-based row coordinate.
     * @param column Zero-based column coordinate.
     * @param numVertexColumns Number of grid vertices per row.
     * @return Row-major vertex identifier.
     */
    static int vertexId(int row, int column, int numVertexColumns) { return (row * numVertexColumns) + column; }

    /** @brief Geometry and doubled signed-area contribution of one directed edge. */
    struct EdgeGeometry {
        /// Grid vertex where the directed edge begins.
        int startVertex = -1;
        /// Grid vertex where the directed edge ends.
        int endVertex = -1;
        /// Cardinal direction from the start vertex to the end vertex.
        Direction direction = Direction::North;
        /// Contribution of the edge to twice the signed boundary area.
        int doubledSignedArea = 0;
    };

    /**
     * @brief Converts a packed edge to its directed grid geometry.
     *
     * @param packedEdge Packed boundary-edge identifier.
     * @return Oriented endpoints, direction, and signed-area contribution.
     */
    EdgeGeometry edgeGeometry(int packedEdge) const {
        const ContourEdge edge = ContourTraceComputation::unpackEdge(packedEdge);
        const int imageColumns = tree_.numColumns();
        const int numVertexColumns = imageColumns + 1;
        const auto [row, column] = ImageUtils::to2D(edge.pixel, imageColumns);

        switch (edge.side) {
        case ContourSide::North:
            return {vertexId(row, column, numVertexColumns), vertexId(row, column + 1, numVertexColumns), Direction::East, -row};
        case ContourSide::East:
            return {vertexId(row, column + 1, numVertexColumns), vertexId(row + 1, column + 1, numVertexColumns), Direction::South, column + 1};
        case ContourSide::South:
            return {vertexId(row + 1, column + 1, numVertexColumns), vertexId(row + 1, column, numVertexColumns), Direction::West, row + 1};
        case ContourSide::West:
            return {vertexId(row + 1, column, numVertexColumns), vertexId(row, column, numVertexColumns), Direction::North, -column};
        }
        throw std::runtime_error("Invalid contour side.");
    }

    /**
     * @brief Ranks one successor direction for the selected connectivity.
     * @param incoming Incoming trace direction.
     * @param outgoing Outgoing trace direction.
     * @param connectivity Digital foreground connectivity.
     * @return Smaller values for preferred successors.
     */
    static int successorPriority(Direction incoming, Direction outgoing, ForegroundConnectivity connectivity) {
        const int quarterTurns = (static_cast<int>(outgoing) - static_cast<int>(incoming) + 4) & 3;
        // Support stays on the right: turning right pairs around the same
        // foreground pixel (4/8); turning left pairs around the background
        // pixel (8/4). The choice is independent of traversal history.
        const std::array<int, 4> priorities =
            connectivity == ForegroundConnectivity::Eight ? std::array<int, 4>{1, 2, 3, 0} : std::array<int, 4>{1, 0, 3, 2};
        return priorities[static_cast<std::size_t>(quarterTurns)];
    }

    /**
     * @brief Selects the successor of one directed contour edge.
     *
     * @param current Current directed edge.
     * @param outgoingHead First candidate edge leaving the terminal vertex.
     * @param directedEdges Directed edges that form the traced contour graph.
     * @param outgoingNext Next-edge links for the outgoing adjacency lists.
     * @param connectivity Foreground connectivity determining diagonal pairing.
     * @return Selected next edge.
     */
    static int selectSuccessor(const DirectedEdge& current, int outgoingHead, const std::vector<DirectedEdge>& directedEdges,
                               const std::vector<int>& outgoingNext, ForegroundConnectivity connectivity) {
        if (connectivity == ForegroundConnectivity::Unknown) {
            throw std::invalid_argument("Diagonal contour tracing requires canonical 4/8 construction adjacency; "
                                        "for complementary shapes, pass a valued-tree view with distinct node/parent altitudes.");
        }
        int bestEdgeIndex = -1;
        int bestPriority = std::numeric_limits<int>::max();
        for (int candidateIndex = outgoingHead; candidateIndex != -1; candidateIndex = outgoingNext[static_cast<std::size_t>(candidateIndex)]) {
            const int priority = successorPriority(current.direction, directedEdges[static_cast<std::size_t>(candidateIndex)].direction, connectivity);
            if (priority < bestPriority || (priority == bestPriority && directedEdges[static_cast<std::size_t>(candidateIndex)].packedEdge <
                                                                            directedEdges[static_cast<std::size_t>(bestEdgeIndex)].packedEdge)) {
                bestEdgeIndex = candidateIndex;
                bestPriority = priority;
            }
        }
        return bestEdgeIndex;
    }

    /**
     * @brief Clears the used entries of the dense vertex index.
     */
    void resetDenseVertexIndex() const {
        for (int vertex : touchedVertices_) {
            const auto vertexIndex = static_cast<std::size_t>(vertex);
            outgoingEdgeHeads_[vertexIndex] = -1;
        }
        touchedVertices_.clear();
    }

    /**
     * @brief Returns the number of vertices in the image grid.
     * @return Number of vertices in the regular grid underlying the image.
     */
    [[nodiscard]] std::size_t numImageVertices() const {
        return static_cast<std::size_t>((tree_.numRows() + 1) * (tree_.numColumns() + 1));
    }

    /**
     * @brief Allocates the dense vertex index when first needed.
     */
    void ensureDenseVertexIndex() const {
        if (!outgoingEdgeHeads_.empty()) {
            return;
        }
        outgoingEdgeHeads_.assign(numImageVertices(), -1);
    }

    /**
     * @brief Rounds a positive value up to a power of two.
     * @param value Positive lower bound.
     * @return Smallest power of two greater than or equal to `value`.
     */
    static std::size_t nextPowerOfTwoAtLeast(std::size_t value) {
        std::size_t result = 1;
        while (result < value) {
            result <<= 1;
        }
        return result;
    }

    /**
     * @brief Chooses a sparse vertex table size for one contour.
     * @param edgeCount Number of boundary edges.
     * @return Power-of-two table size with a load factor below one half.
     */
    static std::size_t sparseVertexTableSize(std::size_t edgeCount) {
        const std::size_t minimumTableSize = std::max<std::size_t>(2, (2 * edgeCount) + 1);
        return nextPowerOfTwoAtLeast(minimumTableSize);
    }

    /**
     * @brief Returns the number of individual queries allowed before using the dense index.
     * @return Maximum number of sparse individual-node traces.
     */
    static constexpr std::size_t sparseIndexQueryLimit() {
        // Use sparse adjacency for a few queries on individual nodes, then use
        // dense storage before repeated access accumulates hash-table overhead.
        return 8;
    }

    /**
     * @brief Chooses the sparse vertex index when it uses less temporary storage.
     * @param edgeCount Number of boundary edges.
     * @param appendToSharedCache Whether descriptors are appended directly to the shared cache.
     * @return True when sparse indexing is selected.
     */
    [[nodiscard]] bool shouldUseSparseVertexIndex(std::size_t edgeCount, bool appendToSharedCache) const {
        if (appendToSharedCache || !outgoingEdgeHeads_.empty()) {
            return false;
        }
        if (numIndividualNodeTraces_ >= sparseIndexQueryLimit()) {
            return false;
        }
        const std::size_t sparseBytes = sparseVertexTableSize(edgeCount) * ((2 * sizeof(int)) + sizeof(uint32_t));
        const std::size_t denseBytes = numImageVertices() * sizeof(int);
        return sparseBytes < denseBytes;
    }

    /**
     * @brief Prepares the sparse vertex index for one contour.
     *
     * @param edgeCount Number of boundary edges.
     */
    void prepareSparseVertexIndex(std::size_t edgeCount) const {
        const std::size_t tableSize = sparseVertexTableSize(edgeCount);
        if (sparseVertexKeys_.size() != tableSize) {
            sparseVertexKeys_.assign(tableSize, 0);
            sparseOutgoingEdgeHeads_.assign(tableSize, -1);
            sparseSlotGenerations_.assign(tableSize, 0);
        }
        touchedSparseSlots_.clear();

        ++sparseGeneration_;
        if (sparseGeneration_ == 0) {
            std::fill(sparseSlotGenerations_.begin(), sparseSlotGenerations_.end(), 0);
            sparseGeneration_ = 1;
        }
    }

    /**
     * @brief Computes the initial sparse table slot for a grid vertex.
     *
     * @param vertex Image-grid vertex identifier.
     * @return Sparse-table slot selected for the vertex.
     */
    [[nodiscard]] std::size_t sparseVertexSlot(int vertex) const {
        const auto hashValue = static_cast<uint32_t>(vertex) * uint32_t{2654435761u};
        return static_cast<std::size_t>(hashValue) & (sparseVertexKeys_.size() - 1);
    }

    /**
     * @brief Returns the first outgoing edge stored for a vertex.
     *
     * @param vertex Image-grid vertex identifier.
     * @return Directed edge index or `-1` when the vertex is absent.
     */
    [[nodiscard]] int sparseOutgoingEdgeHead(int vertex) const {
        std::size_t slot = sparseVertexSlot(vertex);
        while (sparseSlotGenerations_[slot] == sparseGeneration_) {
            if (sparseVertexKeys_[slot] == vertex) {
                return sparseOutgoingEdgeHeads_[slot];
            }
            slot = (slot + 1) & (sparseVertexKeys_.size() - 1);
        }
        return -1;
    }

    /**
     * @brief Finds or inserts a vertex in the sparse index.
     *
     * @param vertex Image-grid vertex identifier.
     * @return Slot containing the existing or newly inserted vertex.
     */
    [[nodiscard]] std::size_t findOrInsertSparseVertex(int vertex) const {
        std::size_t slot = sparseVertexSlot(vertex);
        while (sparseSlotGenerations_[slot] == sparseGeneration_) {
            if (sparseVertexKeys_[slot] == vertex) {
                return slot;
            }
            slot = (slot + 1) & (sparseVertexKeys_.size() - 1);
        }

        sparseSlotGenerations_[slot] = sparseGeneration_;
        sparseVertexKeys_[slot] = vertex;
        sparseOutgoingEdgeHeads_[slot] = -1;
        touchedSparseSlots_.push_back(static_cast<int>(slot));
        return slot;
    }

    /**
     * @brief Returns the first outgoing directed edge of a trace vertex.
     *
     * @param vertex Image-grid vertex identifier.
     * @return Directed-edge index at the head of the vertex list, or the empty sentinel.
     */
    template <VertexIndexMode Mode> [[nodiscard]] int outgoingEdgeHead(int vertex) const {
        if constexpr (Mode == VertexIndexMode::Sparse) {
            return sparseOutgoingEdgeHead(vertex);
        } else {
            return outgoingEdgeHeads_[static_cast<std::size_t>(vertex)];
        }
    }

    /**
     * @brief Prepares the selected vertex index.
     *
     * @param edgeCount Number of boundary edges.
     */
    template <VertexIndexMode Mode> void prepareVertexIndex(std::size_t edgeCount) const {
        if constexpr (Mode == VertexIndexMode::Sparse) {
            prepareSparseVertexIndex(edgeCount);
        } else {
            if (!sparseVertexKeys_.empty()) {
                releaseSparseTraceOutgoingScratch();
            }
            ensureDenseVertexIndex();
            resetDenseVertexIndex();
        }
    }

    /**
     * @brief Clears the selected vertex index after tracing.
     */
    template <VertexIndexMode Mode> void resetVertexIndex() const {
        if constexpr (Mode == VertexIndexMode::Dense) {
            resetDenseVertexIndex();
        }
    }

    /**
     * @brief Adds a directed edge to its start vertex.
     *
     * @param startVertex Starting vertex identifier.
     * @param edgeIndex Boundary-edge index.
     */
    template <VertexIndexMode Mode> void addOutgoingEdge(int startVertex, int edgeIndex) const {
        if constexpr (Mode == VertexIndexMode::Sparse) {
            const std::size_t slot = findOrInsertSparseVertex(startVertex);
            nextOutgoingEdges_.push_back(sparseOutgoingEdgeHeads_[slot]);
            sparseOutgoingEdgeHeads_[slot] = edgeIndex;
        } else {
            const auto vertexIndex = static_cast<std::size_t>(startVertex);
            if (outgoingEdgeHeads_[vertexIndex] == -1) {
                touchedVertices_.push_back(startVertex);
            }
            nextOutgoingEdges_.push_back(outgoingEdgeHeads_[vertexIndex]);
            outgoingEdgeHeads_[vertexIndex] = edgeIndex;
        }
    }

    /**
     * @brief Starts a new visited-edge generation.
     *
     * @param edgeCount Number of boundary edges.
     */
    void advanceVisitGeneration(std::size_t edgeCount) const {
        if (visitedGenerations_.size() < edgeCount) {
            visitedGenerations_.resize(edgeCount, 0);
        }
        ++visitGeneration_;
        if (visitGeneration_ == 0) {
            std::fill(visitedGenerations_.begin(), visitedGenerations_.end(), 0);
            visitGeneration_ = 1;
        }
    }

    /**
     * @brief Tests whether an edge belongs to the current boundary traversal.
     * @param edgeIndex Boundary-edge index.
     * @return True when the edge was visited in the active generation.
     */
    [[nodiscard]] bool isVisited(int edgeIndex) const { return visitedGenerations_[static_cast<std::size_t>(edgeIndex)] == visitGeneration_; }

    /**
     * @brief Marks an edge as visited in the current generation.
     *
     * @param edgeIndex Boundary-edge index.
     */
    void markVisited(int edgeIndex) const { visitedGenerations_[static_cast<std::size_t>(edgeIndex)] = visitGeneration_; }

    /**
     * @brief Reserves boundary descriptors for a subtree trace.
     * @param subtreeRoot Root of the requested subtree.
     */
    void reserveBoundaryCapacityForSubtree(NodeId subtreeRoot) const {
        std::size_t pendingNonEmptyNodes = 0;
        std::size_t pendingEdgeCount = 0;
        for (NodeId node : tree_.subtreeNodes(subtreeRoot)) {
            if (boundariesCached_[static_cast<std::size_t>(node)]) {
                continue;
            }
            const auto edgeCount = static_cast<std::size_t>(edgeCounts_[static_cast<std::size_t>(node)]);
            if (edgeCount == 0) {
                continue;
            }
            ++pendingNonEmptyNodes;
            pendingEdgeCount += edgeCount;
        }
        if (pendingNonEmptyNodes == 0) {
            return;
        }

        const std::size_t edgeBasedEstimate = std::max<std::size_t>(1, pendingEdgeCount / 16);
        const std::size_t additionalCapacity = std::max(pendingNonEmptyNodes, edgeBasedEstimate);
        const std::size_t targetCapacity = cachedBoundaries_.size() + additionalCapacity;
        if (cachedBoundaries_.capacity() < targetCapacity) {
            cachedBoundaries_.reserve(targetCapacity);
        }
    }

    /**
     * @brief Tests whether tracing can reorder the node's cached edge segment.
     * @param node Node identifier.
     * @return True when no uncached parent still depends on the current order.
     */
    bool canReuseEdgeStorage(NodeId node) const {
        if (tree_.isRoot(node)) {
            return true;
        }
        const NodeId parent = tree_.parent(node);
        return parent == InvalidNode || parent == node || edgesCached_[static_cast<std::size_t>(parent)] != 0;
    }

    /**
     * @brief Traces one node when its boundaries are not cached.
     * @param node Live node identifier.
     */
    void ensureNodeBoundariesTraced(NodeId node) const {
        requireStableTree("ContourTraceComputation::ensureNodeBoundariesTraced");
        requireLiveTraceNode(node, "ContourTraceComputation::ensureNodeBoundariesTraced");
        if (boundariesCached_[static_cast<std::size_t>(node)]) {
            return;
        }
        cacheSubtreeEdges(node);
        traceNodeBoundaries(node, false);
    }

    /**
     * @brief Traces every boundary missing from a subtree.
     * @param subtreeRoot Root of the requested subtree.
     */
    void traceSubtree(NodeId subtreeRoot) const {
        requireStableTree("ContourTraceComputation::traceSubtree");
        requireLiveTraceNode(subtreeRoot, "ContourTraceComputation::traceSubtree");
        cacheSubtreeEdges(subtreeRoot);
        reserveBoundaryCapacityForSubtree(subtreeRoot);

        std::vector<std::pair<NodeId, bool>> traversalStack;
        traversalStack.emplace_back(subtreeRoot, false);

        while (!traversalStack.empty()) {
            const auto [node, expanded] = traversalStack.back();
            traversalStack.pop_back();
            if (!expanded) {
                traversalStack.emplace_back(node, true);
                for (NodeId child : tree_.children(node)) {
                    if (!boundariesCached_[static_cast<std::size_t>(child)]) {
                        traversalStack.emplace_back(child, false);
                    }
                }
                continue;
            }
            if (!boundariesCached_[static_cast<std::size_t>(node)]) {
                traceNodeBoundaries(node, true);
            }
        }
    }

    /**
     * @brief Traces and caches all contour boundaries associated with a tree node.
     *
     * @param node Node identifier.
     * @param appendToSharedCache Whether descriptors are appended directly to the shared cache.
     */
    void traceNodeBoundaries(NodeId node, bool appendToSharedCache) const {
        const std::size_t edgeCount = static_cast<std::size_t>(edgeCounts_[static_cast<std::size_t>(node)]);
        const bool useSparseIndex = shouldUseSparseVertexIndex(edgeCount, appendToSharedCache);
        if (!appendToSharedCache) {
            ++numIndividualNodeTraces_;
        }
        if (useSparseIndex) {
            traceNodeBoundariesWithVertexIndex<VertexIndexMode::Sparse>(node, appendToSharedCache, edgeCount);
        } else {
            traceNodeBoundariesWithVertexIndex<VertexIndexMode::Dense>(node, appendToSharedCache, edgeCount);
        }
    }

    /**
     * @brief Traces a node contour with the selected vertex index.
     *
     * @param node Node identifier.
     * @param appendToSharedCache Whether descriptors are appended directly to the shared cache.
     * @param edgeCount Number of boundary edges.
     */
    template <VertexIndexMode Mode> void traceNodeBoundariesWithVertexIndex(NodeId node, bool appendToSharedCache, std::size_t edgeCount) const {
        prepareVertexIndex<Mode>(edgeCount);
        directedEdges_.clear();
        directedEdges_.reserve(edgeCount);
        nextOutgoingEdges_.clear();
        nextOutgoingEdges_.reserve(edgeCount);

        for (auto edge = cachedEdgeBegin(node); edge != cachedEdgeEnd(node); ++edge) {
            const EdgeGeometry geometry = edgeGeometry(*edge);
            const int edgeIndex = static_cast<int>(directedEdges_.size());
            directedEdges_.emplace_back(*edge, geometry.startVertex, geometry.endVertex, geometry.direction, geometry.doubledSignedArea);
            addOutgoingEdge<Mode>(geometry.startVertex, edgeIndex);
        }

        advanceVisitGeneration(directedEdges_.size());
        orderedNodeEdges_.clear();
        orderedNodeEdges_.reserve(directedEdges_.size());
        nodeBoundaries_.clear();
        if (!appendToSharedCache) {
            nodeBoundaries_.reserve((directedEdges_.size() / 4) + 1);
        }

        std::size_t outputEdgeOffset = static_cast<std::size_t>(edgeOffsets_[static_cast<std::size_t>(node)]);
        const bool reuseCachedEdgeStorage = canReuseEdgeStorage(node);
        if (!reuseCachedEdgeStorage) {
            outputEdgeOffset = cachedPackedEdges_.size();
        }
        const std::size_t nodeBoundaryOffset = cachedBoundaries_.size();

        try {
            for (int startEdgeIndex = 0; startEdgeIndex < static_cast<int>(directedEdges_.size()); ++startEdgeIndex) {
                if (isVisited(startEdgeIndex)) {
                    continue;
                }

                const std::size_t firstBoundaryEdgeIndex = orderedNodeEdges_.size();
                int doubledSignedArea = 0;
                int currentEdgeIndex = startEdgeIndex;

                // A cycle may revisit a grid vertex. Only returning to the
                // initial directed edge closes the successor permutation.
                do {
                    if (isVisited(currentEdgeIndex)) {
                        throw std::runtime_error("Contour successor revisited an edge before closing its cycle.");
                    }
                    const DirectedEdge& edge = directedEdges_[static_cast<std::size_t>(currentEdgeIndex)];
                    markVisited(currentEdgeIndex);
                    orderedNodeEdges_.push_back(edge.packedEdge);
                    doubledSignedArea += edge.doubledSignedArea;

                    const int firstSuccessorIndex = outgoingEdgeHead<Mode>(edge.endVertex);
                    if (firstSuccessorIndex == -1) {
                        throw std::runtime_error("Contour boundary contains an edge without a successor.");
                    }
                    if (nextOutgoingEdges_[static_cast<std::size_t>(firstSuccessorIndex)] == -1) {
                        currentEdgeIndex = firstSuccessorIndex;
                    } else {
                        currentEdgeIndex = selectSuccessor(edge, firstSuccessorIndex, directedEdges_, nextOutgoingEdges_,
                                                      connectivityByNode_[static_cast<std::size_t>(node)]);
                        if (currentEdgeIndex == -1) {
                            throw std::runtime_error("Contour boundary contains an unresolved successor.");
                        }
                    }
                } while (currentEdgeIndex != startEdgeIndex);

                const std::size_t boundaryEdgeCount = orderedNodeEdges_.size() - firstBoundaryEdgeIndex;
                if (boundaryEdgeCount == 0) {
                    continue;
                }

                const ContourBoundaryKind kind = doubledSignedArea >= 0 ? ContourBoundaryKind::External : ContourBoundaryKind::Internal;
                if (appendToSharedCache) {
                    cachedBoundaries_.push_back(ContourBoundary{kind, checkedUint32(outputEdgeOffset + firstBoundaryEdgeIndex, "global boundary edge offset"),
                                                                checkedUint32(boundaryEdgeCount, "boundary edge count"), doubledSignedArea});
                } else {
                    nodeBoundaries_.push_back(ContourBoundary{kind, checkedUint32(firstBoundaryEdgeIndex, "node boundary edge offset"),
                                                               checkedUint32(boundaryEdgeCount, "boundary edge count"), doubledSignedArea});
                }
            }
        } catch (...) {
            if (appendToSharedCache) {
                cachedBoundaries_.resize(nodeBoundaryOffset);
            }
            resetVertexIndex<Mode>();
            throw;
        }
        resetVertexIndex<Mode>();

        if (orderedNodeEdges_.size() != edgeCount) {
            if (appendToSharedCache) {
                cachedBoundaries_.resize(nodeBoundaryOffset);
            }
            throw std::runtime_error("Contour trace boundary traversal did not cover every cached edge.");
        }
        if (reuseCachedEdgeStorage) {
            std::copy(orderedNodeEdges_.begin(), orderedNodeEdges_.end(),
                      cachedPackedEdges_.begin() + static_cast<std::ptrdiff_t>(outputEdgeOffset));
        } else {
            cachedPackedEdges_.insert(cachedPackedEdges_.end(), orderedNodeEdges_.begin(), orderedNodeEdges_.end());
        }

        if (!appendToSharedCache) {
            for (ContourBoundary& boundary : nodeBoundaries_) {
                boundary.edgeOffset = checkedUint32(outputEdgeOffset + boundary.edgeOffset, "global boundary edge offset");
            }
            cachedBoundaries_.insert(cachedBoundaries_.end(), nodeBoundaries_.begin(), nodeBoundaries_.end());
        }

        boundaryOffsets_[static_cast<std::size_t>(node)] = checkedUint32(nodeBoundaryOffset, "boundary offset");
        boundaryCounts_[static_cast<std::size_t>(node)] = checkedUint32(cachedBoundaries_.size() - nodeBoundaryOffset, "boundary count");
        const auto nodeIndex = static_cast<std::size_t>(node);
        if (!boundariesCached_[nodeIndex]) {
            boundariesCached_[nodeIndex] = 1;
            ++numNodesWithCachedBoundaries_;
        }
        releaseTraceScratchIfComplete();
    }

  private:
    /**
     * @brief Converts canonical grid adjacency to digital foreground connectivity.
     * @param adjacency Regular-grid adjacency retained by tree construction.
     * @return Four, Eight, or Unknown when the adjacency is noncanonical.
     */
    static ForegroundConnectivity foregroundConnectivity(const RegularGridAdjacency2D& adjacency) {
        if (adjacency.is4connectivity()) {
            return ForegroundConnectivity::Four;
        }
        if (adjacency.is8connectivity()) {
            return ForegroundConnectivity::Eight;
        }
        return ForegroundConnectivity::Unknown;
    }

    /**
     * @brief Returns lower and upper connectivity from retained construction semantics.
     * @param tree Tree whose construction metadata defines shape connectivity.
     * @return Foreground connectivity for lower and upper shapes, respectively.
     */
    static std::array<ForegroundConnectivity, 2> shapeForegroundConnectivities(const MorphologicalTree& tree) {
        if (const auto* adjacency = detail::constructionAdjacency(tree)) {
            const auto connectivity = foregroundConnectivity(*adjacency);
            return {connectivity, connectivity};
        }
        if (const auto* adjacencies = detail::complementaryAdjacencies(tree)) {
            return {foregroundConnectivity(adjacencies->minAdjacency), foregroundConnectivity(adjacencies->maxAdjacency)};
        }
        if (const auto* convention = tree.topographicConvention();
            convention && std::holds_alternative<SelfDualSpanImmersion>(convention->immersion)) {
            // The projected self-dual span supports use the existing 4/4
            // convention, also used by the scalar bitquad projection.
            return {ForegroundConnectivity::Four, ForegroundConnectivity::Four};
        }
        return {ForegroundConnectivity::Unknown, ForegroundConnectivity::Unknown};
    }

    /**
     * @brief Captures topology-only choices; mixed lower/upper choices remain unresolved.
     * @param tree Stable source tree.
     * @return Foreground connectivity for each internal node slot.
     */
    static std::vector<ForegroundConnectivity> foregroundConnectivityByNode(const MorphologicalTree& tree) {
        const auto [lowerShape, upperShape] = shapeForegroundConnectivities(tree);
        return std::vector<ForegroundConnectivity>(static_cast<std::size_t>(tree.numInternalNodeSlots()),
                                                   lowerShape == upperShape ? lowerShape : ForegroundConnectivity::Unknown);
    }

    /**
     * @brief Captures shape polarity so lazy tracing never borrows the altitude buffer.
     * @param view Current valued tree view.
     * @return Foreground connectivity for each internal node slot.
     */
    template <AltitudeValue T>
    static std::vector<ForegroundConnectivity> foregroundConnectivityByNode(const ValuedMorphologicalTreeView<T>& view) {
        const auto& tree = view.topology();
        auto connectivityByNode = foregroundConnectivityByNode(tree);
        const auto [lowerShape, upperShape] = shapeForegroundConnectivities(tree);
        if (lowerShape != upperShape) {
            for (NodeId node : tree.aliveNodeIds()) {
                if (tree.isRoot(node)) {
                    continue; // The rectangular domain boundary has no diagonal ambiguity.
                }
                const auto altitude = view.nodeAltitude(node);
                const auto parentAltitude = view.nodeAltitude(tree.parent(node));
                if (altitude < parentAltitude) {
                    connectivityByNode[static_cast<std::size_t>(node)] = lowerShape;
                } else if (altitude > parentAltitude) {
                    connectivityByNode[static_cast<std::size_t>(node)] = upperShape;
                }
            }
        }
        return connectivityByNode;
    }

    /// Temporary packed-edge lists used before compaction by node.
    using PendingEdgeLists = contours::detail::PendingEdgeLists;

    /**
     * @brief Returns the pixel adjacent to one side of a source pixel.
     *
     * @param tree Tree topology.
     * @param pixel Pixel identifier.
     * @param side Contour or boundary side.
     * @return Identifier of the neighboring pixel.
     */
    static PixelId adjacentPixel(const MorphologicalTree& tree, PixelId pixel, ContourSide side) {
        const int rows = tree.numRows();
        const int columns = tree.numColumns();
        const auto [row, column] = ImageUtils::to2D(pixel, columns);

        switch (side) {
        case ContourSide::North:
            return row == 0 ? -1 : ImageUtils::to1D(row - 1, column, columns);
        case ContourSide::West:
            return column == 0 ? -1 : ImageUtils::to1D(row, column - 1, columns);
        case ContourSide::East:
            return column == columns - 1 ? -1 : ImageUtils::to1D(row, column + 1, columns);
        case ContourSide::South:
            return row == rows - 1 ? -1 : ImageUtils::to1D(row + 1, column, columns);
        }
        return -1;
    }

    /**
     * @brief Builds compact edge changes for every tree node.
     * @param tree Source topology.
     * @return Edge changes and an output size estimate.
     */
    [[nodiscard]] static ConstructionData prepareEdgeDeltas(const MorphologicalTree& tree) {
        if (tree.numRows() <= 0 || tree.numColumns() <= 0) {
            throw std::invalid_argument("Contour tracing requires a non-empty image domain.");
        }
        if (!tree.isAlive(tree.root())) {
            throw std::invalid_argument("Contour tracing requires a live tree root.");
        }

        const int numNodes = tree.numInternalNodeSlots();
        const int totalPixels = tree.numRows() * tree.numColumns();
        const int numPossibleEdges = 4 * totalPixels;
        const int estimatedCachedEdgeCount = std::max(totalPixels, 1);

        PendingEdgeLists pendingAdditions(numNodes, estimatedCachedEdgeCount);
        PendingEdgeLists pendingRemovals(numNodes, estimatedCachedEdgeCount);

        static constexpr std::array<ContourSide, 4> sides{ContourSide::North, ContourSide::West, ContourSide::East,
                                                               ContourSide::South};

        detail::traversePostOrder(
            tree, tree.root(), [](NodeId) -> void {}, [](NodeId, NodeId) -> void {},
            [&](NodeId nodeId) {
                for (PixelId pixel : tree.properPart(nodeId)) {
                    for (ContourSide side : sides) {
                        const int packedEdge = packEdge(pixel, side);
                        const PixelId neighbor = adjacentPixel(tree, pixel, side);
                        if (neighbor < 0) {
                            pendingAdditions.add(nodeId, packedEdge);
                            continue;
                        }

                        const NodeId entry = local_attributes::detail::kernel::anchoredEntry(tree, pixel, neighbor);
                        if (entry == InvalidNode || entry == nodeId) {
                            continue;
                        }

                        pendingAdditions.add(nodeId, packedEdge);
                        pendingRemovals.add(entry, packedEdge);
                    }
                }
            });

        return {EdgeDeltas::compact(pendingAdditions, pendingRemovals, numPossibleEdges), estimatedCachedEdgeCount, {}};
    }

    /**
     * @brief Validates the topology and prepares edge changes and connectivity.
     * @param tree Stable source tree.
     * @return Construction buffers for topology-only tracing.
     */
    [[nodiscard]] static ConstructionData prepareConstructionData(const MorphologicalTree& tree) {
        tree.requireNotEditing("ContourTraceComputation");
        ConstructionData constructionData = prepareEdgeDeltas(tree);
        constructionData.connectivityByNode = foregroundConnectivityByNode(tree);
        return constructionData;
    }

    /**
     * @brief Captures lower and upper shape connectivity from a current valued view.
     * @param view Current valued tree view.
     * @return Construction buffers with connectivity resolved from node polarity.
     */
    template <AltitudeValue T> [[nodiscard]] static ConstructionData prepareConstructionData(const ValuedMorphologicalTreeView<T>& view) {
        view.requireTopologyUnchanged("ContourTraceComputation");
        view.topology().requireNotEditing("ContourTraceComputation");
        ConstructionData constructionData = prepareEdgeDeltas(view.topology());
        constructionData.connectivityByNode = foregroundConnectivityByNode(view);
        return constructionData;
    }

  public:
    /**
     * @brief Prepares compact edge changes for lazy edge and boundary queries.
     *
     * The committed tree must outlive this object and remain unchanged. Edges
     * and boundaries are cached only when requested; traceAll() prepares all.
     * @param tree Source topology with a non-empty regular 2D domain.
     */
    explicit ContourTraceComputation(const MorphologicalTree& tree) : ContourTraceComputation(tree, prepareConstructionData(tree)) {}

    /**
     * @brief Constructs traces with the valued view's lower/upper shape connectivity.
     *
     * Connectivity is captured during construction. The altitude span is not
     * retained; the source topology must outlive this object and stay unchanged.
     * @param view Current valued view supplying topology and shape polarity.
     */
    template <AltitudeValue T>
    explicit ContourTraceComputation(const ValuedMorphologicalTreeView<T>& view)
        : ContourTraceComputation(view.topology(), prepareConstructionData(view)) {}
};

} // namespace mmcfilters
