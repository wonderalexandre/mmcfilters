#pragma once

#include "../trees/MorphologicalTree.hpp"
#include "../trees/WeightedTreeView.hpp"
#include "../trees/detail/ProperPartEntryNode.hpp"
#include "../trees/detail/TreeTraversalDetail.hpp"
#include "../utils/Common.hpp"
#include "detail/ContourTraceDeltaStore.hpp"
#include "detail/PendingPixelLists.hpp"

#include <algorithm>
#include <array>
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
#include <chrono>
#endif
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mmcfilters {

/**
 * @brief Oriented side of one support pixel used as a boundary trace edge.
 */
enum class ContourTraceSide : uint8_t {
    North = 0,
    West = 1,
    East = 2,
    South = 3
};

/**
 * @brief Geometric class of a traced boundary loop.
 */
enum class ContourLoopKind : uint8_t {
    External,
    Internal
};

/**
 * @brief One unpacked boundary edge attached to a support pixel.
 */
struct ContourTraceEdge {
    /// Row-major support-pixel index incident to the boundary edge.
    int pixel = -1;
    /// Side of the support pixel occupied by the boundary edge.
    ContourTraceSide side = ContourTraceSide::North;

    /// Compares the support pixel and side.
    friend bool operator==(const ContourTraceEdge&, const ContourTraceEdge&) = default;
};

/**
 * @brief Metadata for one ordered boundary loop.
 */
struct ContourTraceLoop {
    /// Whether this loop is an external boundary or an internal hole.
    ContourLoopKind kind = ContourLoopKind::External;
    /// First edge in the shared ordered-edge buffer.
    uint32_t edgeOffset = 0;
    /// Number of consecutive edges in the loop.
    uint32_t edgeCount = 0;
    /// Doubled signed area under the trace-orientation convention.
    int signedArea2 = 0;
};

/**
 * @brief Incremental side-level contour extraction and boundary-loop tracing.
 *
 * Boundary primitives are oriented grid edges. Edges are materialized lazily per
 * node, then traced into ordered loops on demand. The orientation convention
 * keeps the support pixel on the right side of each directed edge in image
 * coordinates, where rows grow downward and columns grow rightward. With this
 * convention, external loops have positive doubled signed area and internal
 * loops have negative doubled signed area.
 */
class ContourTraceComputation {
public:
    /// Compact storage for local contour additions and removals.
    using LocalTraceDeltas = detail::ContourTraceDeltaStore;

    /**
     * @brief Packs one pixel-side edge into a compact integer id.
     */
    [[nodiscard]] static int packEdge(int pixel, ContourTraceSide side) {
        return (4 * pixel) + static_cast<int>(side);
    }

    /**
     * @brief Unpacks one compact edge id.
     */
    [[nodiscard]] static ContourTraceEdge unpackEdge(int packedEdge) {
        if (packedEdge < 0) {
            return {};
        }
        const int side = packedEdge & 3;
        return ContourTraceEdge{
            packedEdge / 4,
            static_cast<ContourTraceSide>(side)};
    }

    /**
     * @brief Lazy result for side-level contour traces.
     */
    struct IncrementalContourTraces {
    private:
        friend class ContourTraceComputation;

        enum class Direction : uint8_t {
            North = 0,
            East = 1,
            South = 2,
            West = 3
        };

        enum class TraceAdjacencyMode : uint8_t {
            Dense,
            Sparse
        };

        struct DirectedEdge {
            int packedEdge = -1;
            int startVertex = -1;
            int endVertex = -1;
            Direction direction = Direction::North;
            int signedArea2 = 0;

            DirectedEdge() = default;

            DirectedEdge(
                int packedEdge,
                int startVertex,
                int endVertex,
                Direction direction,
                int signedArea2)
                : packedEdge(packedEdge),
                  startVertex(startVertex),
                  endVertex(endVertex),
                  direction(direction),
                  signedArea2(signedArea2) {}
        };

#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
    public:
        struct TraceProfileStats {
            std::size_t nodesTraced = 0;
            std::size_t edgesTraced = 0;
            std::size_t loopsTraced = 0;
            std::size_t outgoingVertices = 0;
            std::size_t singleOutgoingVertices = 0;
            std::size_t multiOutgoingVertices = 0;
            std::size_t maxOutgoingDegree = 0;
            std::size_t closedLoopStops = 0;
            std::size_t missingOutgoingStops = 0;
            std::size_t singleSuccessorSteps = 0;
            std::size_t singleSuccessorVisitedStops = 0;
            std::size_t ambiguousSuccessorSteps = 0;
            std::size_t ambiguousSuccessorDeadEnds = 0;
            std::size_t successorCandidateScans = 0;
            std::size_t successorVisitedSkips = 0;
            std::size_t successorUnvisitedCandidates = 0;
            std::int64_t profileCountersNs = 0;
            std::int64_t buildAdjacencyNs = 0;
            std::int64_t walkLoopsNs = 0;
            std::int64_t resetOutgoingNs = 0;
            std::int64_t commitEdgesNs = 0;
            std::int64_t commitLoopsNs = 0;
            std::int64_t releaseScratchNs = 0;
        };

    private:
        mutable TraceProfileStats* activeTraceProfile_ = nullptr;
#endif

        const MorphologicalTree& tree;
        std::size_t treeMutationVersion_ = 0;
        mutable LocalTraceDeltas localDeltas_;

        mutable std::vector<int> cachedEdgeValues_;
        mutable std::vector<uint32_t> cachedEdgeOffset_;
        mutable std::vector<uint32_t> cachedEdgeSize_;
        mutable std::vector<uint8_t> cachedEdgeReady_;
        mutable std::size_t cachedEdgeReadyCount_ = 0;

        mutable std::vector<ContourTraceLoop> cachedLoopInfos_;
        mutable std::vector<uint32_t> cachedLoopInfoOffset_;
        mutable std::vector<uint32_t> cachedLoopInfoSize_;
        mutable std::vector<uint8_t> cachedLoopReady_;
        mutable std::size_t cachedLoopReadyCount_ = 0;
        mutable std::size_t cachedLoopEdgeCount_ = 0;

        mutable std::vector<uint16_t> edgeMark_;
        mutable uint16_t markGeneration_ = 1;
        mutable bool edgeMaterializationScratchReleased_ = false;

        mutable std::vector<DirectedEdge> traceDirectedEdges_;
        mutable std::vector<int> traceOutgoingHead_;
        mutable std::vector<int> traceOutgoingNext_;
        mutable std::vector<int> traceTouchedVertices_;
        mutable std::vector<int> traceSparseVertexKeys_;
        mutable std::vector<int> traceSparseOutgoingHead_;
        mutable std::vector<uint32_t> traceSparseSlotGeneration_;
        mutable std::vector<int> traceSparseTouchedSlots_;
        mutable uint32_t traceSparseGeneration_ = 1;
        mutable std::size_t nodeLocalLoopTraceCount_ = 0;
        mutable std::vector<uint32_t> traceVisitedGeneration_;
        mutable uint32_t traceVisitGeneration_ = 1;
        mutable std::vector<int> traceNodeLoopEdges_;
        mutable std::vector<ContourTraceLoop> traceNodeLoops_;
        mutable bool traceScratchReleased_ = false;

        IncrementalContourTraces(
            const MorphologicalTree& tree,
            LocalTraceDeltas localDeltas,
            int capacityHint)
            : tree(tree),
              treeMutationVersion_(tree.getMutationVersion()),
              localDeltas_(std::move(localDeltas)),
              cachedEdgeOffset_(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), 0),
              cachedEdgeSize_(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), 0),
              cachedEdgeReady_(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), 0),
              cachedLoopInfoOffset_(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), 0),
              cachedLoopInfoSize_(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), 0),
              cachedLoopReady_(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), 0),
              edgeMark_(static_cast<std::size_t>(4 * tree.getNumRowsOfImage() * tree.getNumColsOfImage()), 0) {
            if (capacityHint > 0) {
                cachedEdgeValues_.reserve(static_cast<std::size_t>(capacityHint));
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
                using value_type = ContourTraceEdge;
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
                 * @param values Packed-edge buffer traversed by the iterator.
                 * @param index Zero-based position in the buffer.
                 */
                iterator(const std::vector<int>* values, std::size_t index)
                    : values_(values), index_(index) {}

                /**
                 * @brief Returns the unpacked edge at the current position.
                 *
                 * @return Unpacked edge at the current position.
                 */
                value_type operator*() const {
                    return ContourTraceComputation::unpackEdge((*values_)[index_]);
                }

                /**
                 * @brief Advances to the next packed edge.
                 *
                 * @return Reference to the advanced iterator.
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
                    iterator tmp(*this);
                    ++(*this);
                    return tmp;
                }

                /// Compares the backing buffer and position.
                friend bool operator==(const iterator& lhs, const iterator& rhs) {
                    return lhs.values_ == rhs.values_ && lhs.index_ == rhs.index_;
                }

                /// Returns true when two iterator positions differ.
                friend bool operator!=(const iterator& lhs, const iterator& rhs) {
                    return !(lhs == rhs);
                }

            private:
                const std::vector<int>* values_ = nullptr;
                std::size_t index_ = 0;
            };

            /**
             * @brief Creates an empty edge range.
             */
            EdgeRange() = default;

            /**
             * @brief Creates a view over `size` packed edges starting at `offset`.
             *
             * @param values Packed-edge buffer viewed by the range.
             * @param offset First position in the buffer.
             * @param size Number of edges in the range.
             */
            EdgeRange(const std::vector<int>* values, std::size_t offset, std::size_t size)
                : values_(values), offset_(offset), size_(size) {}

            /**
             * @brief Returns an iterator to the first edge.
             *
             * @return Iterator to the first edge.
             */
            iterator begin() const {
                return iterator(values_, offset_);
            }

            /**
             * @brief Returns the exclusive end iterator.
             *
             * @return Exclusive end iterator.
             */
            iterator end() const {
                return iterator(values_, offset_ + size_);
            }

            /**
             * @brief Returns whether the range contains no edges.
             *
             * @return `true` when the range is empty.
             */
            [[nodiscard]] bool empty() const noexcept {
                return size_ == 0;
            }

            /**
             * @brief Returns the number of edges in the range.
             *
             * @return Number of edges in the range.
             */
            [[nodiscard]] std::size_t size() const noexcept {
                return size_;
            }

        private:
            const std::vector<int>* values_ = nullptr;
            std::size_t offset_ = 0;
            std::size_t size_ = 0;
        };

#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
        /**
         * @brief Cache and allocation diagnostics for internal profiling builds.
         */
        struct StorageStats {
            std::size_t addDeltaValues = 0;
            std::size_t removeDeltaValues = 0;
            std::size_t cachedEdgeValues = 0;
            std::size_t cachedLoopEdges = 0;
            std::size_t cachedLoops = 0;
            std::size_t cachedEdgeReadyNodes = 0;
            std::size_t cachedLoopReadyNodes = 0;
            std::size_t traceDenseOutgoingSlots = 0;
            std::size_t traceSparseOutgoingSlots = 0;
            std::size_t approxAllocatedBytes = 0;
        };

        [[nodiscard]] StorageStats storageStats() const {
            requireStableTree("IncrementalContourTraces::storageStats");
            StorageStats stats;
            stats.addDeltaValues = localDeltas_.addValues.size();
            stats.removeDeltaValues = localDeltas_.removeValues.size();
            stats.cachedEdgeValues = cachedEdgeValues_.size();
            stats.cachedLoopEdges = cachedLoopEdgeCount_;
            stats.cachedLoops = cachedLoopInfos_.size();
            stats.cachedEdgeReadyNodes = cachedEdgeReadyCount_;
            stats.cachedLoopReadyNodes = cachedLoopReadyCount_;
            stats.traceDenseOutgoingSlots = traceOutgoingHead_.capacity();
            stats.traceSparseOutgoingSlots = traceSparseVertexKeys_.capacity();
            stats.approxAllocatedBytes =
                localDeltas_.addValues.capacity() * sizeof(int) +
                localDeltas_.removeValues.capacity() * sizeof(int) +
                localDeltas_.addSpans.capacity() * sizeof(LocalTraceDeltas::Span) +
                localDeltas_.removeSpans.capacity() * sizeof(LocalTraceDeltas::Span) +
                cachedEdgeValues_.capacity() * sizeof(int) +
                cachedEdgeOffset_.capacity() * sizeof(uint32_t) +
                cachedEdgeSize_.capacity() * sizeof(uint32_t) +
                cachedEdgeReady_.capacity() * sizeof(uint8_t) +
                cachedLoopInfos_.capacity() * sizeof(ContourTraceLoop) +
                cachedLoopInfoOffset_.capacity() * sizeof(uint32_t) +
                cachedLoopInfoSize_.capacity() * sizeof(uint32_t) +
                cachedLoopReady_.capacity() * sizeof(uint8_t) +
                edgeMark_.capacity() * sizeof(uint16_t) +
                traceDirectedEdges_.capacity() * sizeof(DirectedEdge) +
                traceOutgoingHead_.capacity() * sizeof(int) +
                traceOutgoingNext_.capacity() * sizeof(int) +
                traceTouchedVertices_.capacity() * sizeof(int) +
                traceSparseVertexKeys_.capacity() * sizeof(int) +
                traceSparseOutgoingHead_.capacity() * sizeof(int) +
                traceSparseSlotGeneration_.capacity() * sizeof(uint32_t) +
                traceSparseTouchedSlots_.capacity() * sizeof(int) +
                traceVisitedGeneration_.capacity() * sizeof(uint32_t) +
                traceNodeLoopEdges_.capacity() * sizeof(int) +
                traceNodeLoops_.capacity() * sizeof(ContourTraceLoop);
            return stats;
        }
#endif

        /**
         * @brief Returns unordered materialized boundary edges for one node.
         */
        [[nodiscard]] EdgeRange getEdges(NodeId node) const {
            requireStableTree("IncrementalContourTraces::getEdges");
            requireLiveTraceNode(node, "IncrementalContourTraces::getEdges");
            ensureEdgesMaterialized(node);
            return EdgeRange(
                &cachedEdgeValues_,
                cachedEdgeOffset_[static_cast<std::size_t>(node)],
                cachedEdgeSize_[static_cast<std::size_t>(node)]);
        }

        /**
         * @brief Returns loop metadata for one node.
         */
        [[nodiscard]] std::span<const ContourTraceLoop> getLoops(NodeId node) const {
            requireStableTree("IncrementalContourTraces::getLoops");
            requireLiveTraceNode(node, "IncrementalContourTraces::getLoops");
            ensureNodeLoopsMaterialized(node);
            const auto offset = static_cast<std::size_t>(cachedLoopInfoOffset_[static_cast<std::size_t>(node)]);
            const auto size = static_cast<std::size_t>(cachedLoopInfoSize_[static_cast<std::size_t>(node)]);
            if (size == 0) {
                return {};
            }
            return std::span<const ContourTraceLoop>(cachedLoopInfos_.data() + offset, size);
        }

        /**
         * @brief Returns the ordered edges belonging to one loop.
         */
        [[nodiscard]] EdgeRange getLoopEdges(const ContourTraceLoop& loop) const {
            requireStableTree("IncrementalContourTraces::getLoopEdges");
            const std::size_t offset = loop.edgeOffset;
            const std::size_t size = loop.edgeCount;
            if (offset > cachedEdgeValues_.size() || size > cachedEdgeValues_.size() - offset) {
                throw std::invalid_argument("ContourTraceLoop does not belong to this contour-trace result.");
            }
            return EdgeRange(&cachedEdgeValues_, offset, size);
        }

        /**
         * @brief Materializes and traces every live node.
         */
        void materializeAll() const {
            requireStableTree("IncrementalContourTraces::materializeAll");
            ensureLoopsMaterialized(tree.getRoot());
        }

        /**
         * @brief Returns whether every live node has a traced-loop cache.
         *
         * @return `true` when all live nodes have been traced.
         */
        [[nodiscard]] bool isMaterialized() const {
            requireStableTree("IncrementalContourTraces::isMaterialized");
            for (NodeId node : tree.getAliveNodeIds()) {
                if (!cachedLoopReady_[static_cast<std::size_t>(node)]) {
                    return false;
                }
            }
            return true;
        }

        /**
         * @brief Returns whether the boundary edges of one node are cached.
         *
         * @param node Node whose edge cache is queried.
         * @return `true` when the node's boundary edges are cached.
         */
        [[nodiscard]] bool isEdgeMaterialized(NodeId node) const {
            requireStableTree("IncrementalContourTraces::isEdgeMaterialized");
            requireLiveTraceNode(node, "IncrementalContourTraces::isEdgeMaterialized");
            return static_cast<bool>(cachedEdgeReady_[static_cast<std::size_t>(node)]);
        }

        /**
         * @brief Returns whether the boundary loops of one node are cached.
         *
         * @param node Node whose loop cache is queried.
         * @return `true` when the node's boundary loops are cached.
         */
        [[nodiscard]] bool isNodeTraced(NodeId node) const {
            requireStableTree("IncrementalContourTraces::isNodeTraced");
            requireLiveTraceNode(node, "IncrementalContourTraces::isNodeTraced");
            return static_cast<bool>(cachedLoopReady_[static_cast<std::size_t>(node)]);
        }

#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
        [[nodiscard]] TraceProfileStats profileMaterializeAllLoops() const {
            requireStableTree("IncrementalContourTraces::profileMaterializeAllLoops");
            TraceProfileStats stats;
            TraceProfileStats* previousProfile = activeTraceProfile_;
            activeTraceProfile_ = &stats;
            try {
                ensureLoopsMaterialized(tree.getRoot());
            } catch (...) {
                activeTraceProfile_ = previousProfile;
                throw;
            }
            activeTraceProfile_ = previousProfile;
            return stats;
        }
#endif

    private:
        void requireStableTree(const char* context) const {
            tree.requireMutationVersion(treeMutationVersion_, context);
        }

        void requireLiveTraceNode(NodeId node, const char* context) const {
            if (!tree.isAlive(node)) {
                throw std::invalid_argument(std::string(context) + " requires a live internal NodeId.");
            }
        }

        static uint32_t checkedU32(std::size_t value, const char* context) {
            if (value > static_cast<std::size_t>(std::numeric_limits<uint32_t>::max())) {
                throw std::overflow_error(std::string(context) + " exceeds uint32_t.");
            }
            return static_cast<uint32_t>(value);
        }

        template <class T>
        static void releaseVector(std::vector<T>& values) {
            std::vector<T>().swap(values);
        }

#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
        using TraceProfileClock = std::chrono::steady_clock;

        [[nodiscard]] static TraceProfileClock::time_point traceProfileNow() {
            return TraceProfileClock::now();
        }

        static std::int64_t traceProfileElapsedNs(
            TraceProfileClock::time_point start,
            TraceProfileClock::time_point end) {
            return std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        }
#endif

        [[nodiscard]] bool allEdgesMaterialized() const {
            return cachedEdgeReadyCount_ == static_cast<std::size_t>(tree.getNumNodes());
        }

        [[nodiscard]] bool allLoopsMaterialized() const {
            return cachedLoopReadyCount_ == static_cast<std::size_t>(tree.getNumNodes());
        }

        void releaseEdgeMaterializationScratchIfComplete() const {
            if (edgeMaterializationScratchReleased_ || !allEdgesMaterialized()) {
                return;
            }

            localDeltas_ = LocalTraceDeltas{};
            releaseVector(edgeMark_);
            markGeneration_ = 1;
            edgeMaterializationScratchReleased_ = true;
        }

        void releaseSparseTraceOutgoingScratch() const {
            releaseVector(traceSparseVertexKeys_);
            releaseVector(traceSparseOutgoingHead_);
            releaseVector(traceSparseSlotGeneration_);
            releaseVector(traceSparseTouchedSlots_);
            traceSparseGeneration_ = 1;
        }

        void releaseTraceScratchIfComplete() const {
            if (traceScratchReleased_ || !allLoopsMaterialized()) {
                return;
            }

            resetTraceOutgoingHeads();
            releaseVector(traceDirectedEdges_);
            releaseVector(traceOutgoingHead_);
            releaseVector(traceOutgoingNext_);
            releaseVector(traceTouchedVertices_);
            releaseSparseTraceOutgoingScratch();
            releaseVector(traceVisitedGeneration_);
            releaseVector(traceNodeLoopEdges_);
            releaseVector(traceNodeLoops_);
            nodeLocalLoopTraceCount_ = 0;
            traceVisitGeneration_ = 1;
            traceScratchReleased_ = true;
        }

        void nextMarkGeneration() const {
            ++markGeneration_;
            if (markGeneration_ == 0) {
                std::fill(edgeMark_.begin(), edgeMark_.end(), 0);
                markGeneration_ = 1;
            }
        }

        void addIfUnmarked(std::vector<int>& values, int packedEdge) const {
            if (packedEdge < 0 || packedEdge >= static_cast<int>(edgeMark_.size())) {
                return;
            }
            if (edgeMark_[static_cast<std::size_t>(packedEdge)] != markGeneration_) {
                edgeMark_[static_cast<std::size_t>(packedEdge)] = markGeneration_;
                values.push_back(packedEdge);
            }
        }

        void removeIfMarked(int packedEdge) const {
            if (packedEdge >= 0 && packedEdge < static_cast<int>(edgeMark_.size())) {
                edgeMark_[static_cast<std::size_t>(packedEdge)] = 0;
            }
        }

        std::vector<int>::const_iterator cachedEdgeBegin(NodeId node) const {
            return cachedEdgeValues_.begin() + static_cast<std::ptrdiff_t>(cachedEdgeOffset_[static_cast<std::size_t>(node)]);
        }

        std::vector<int>::const_iterator cachedEdgeEnd(NodeId node) const {
            return cachedEdgeBegin(node) + static_cast<std::ptrdiff_t>(cachedEdgeSize_[static_cast<std::size_t>(node)]);
        }

        void commitMaterializedEdges(NodeId node, const std::vector<int>& values) const {
            cachedEdgeOffset_[static_cast<std::size_t>(node)] = checkedU32(cachedEdgeValues_.size(), "cached trace edge offset");
            cachedEdgeSize_[static_cast<std::size_t>(node)] = checkedU32(values.size(), "cached trace edge size");
            cachedEdgeValues_.insert(cachedEdgeValues_.end(), values.begin(), values.end());
            const auto index = static_cast<std::size_t>(node);
            if (!cachedEdgeReady_[index]) {
                cachedEdgeReady_[index] = 1;
                ++cachedEdgeReadyCount_;
            }
        }

        void ensureEdgesMaterialized(NodeId root) const {
            requireStableTree("IncrementalContourTraces::ensureEdgesMaterialized");
            requireLiveTraceNode(root, "IncrementalContourTraces::ensureEdgesMaterialized");
            if (cachedEdgeReady_[static_cast<std::size_t>(root)]) {
                return;
            }

            std::vector<std::pair<NodeId, bool>> stack;
            stack.emplace_back(root, false);
            std::vector<int> values;

            while (!stack.empty()) {
                const auto [node, expanded] = stack.back();
                stack.pop_back();
                if (cachedEdgeReady_[static_cast<std::size_t>(node)]) {
                    continue;
                }
                if (!expanded) {
                    stack.emplace_back(node, true);
                    for (NodeId child : tree.getChildren(node)) {
                        if (!cachedEdgeReady_[static_cast<std::size_t>(child)]) {
                            stack.emplace_back(child, false);
                        }
                    }
                    continue;
                }

                values.clear();
                const auto additions = localDeltas_.additions(node);
                std::size_t reserveSize = additions.size();
                for (NodeId child : tree.getChildren(node)) {
                    reserveSize += static_cast<std::size_t>(cachedEdgeSize_[static_cast<std::size_t>(child)]);
                }
                values.reserve(reserveSize);
                nextMarkGeneration();

                for (NodeId child : tree.getChildren(node)) {
                    for (auto it = cachedEdgeBegin(child); it != cachedEdgeEnd(child); ++it) {
                        addIfUnmarked(values, *it);
                    }
                }

                for (int packedEdge : additions) {
                    addIfUnmarked(values, packedEdge);
                }

                for (int packedEdge : localDeltas_.removals(node)) {
                    removeIfMarked(packedEdge);
                }

                std::size_t writeIndex = 0;
                for (int packedEdge : values) {
                    if (edgeMark_[static_cast<std::size_t>(packedEdge)] == markGeneration_) {
                        values[writeIndex++] = packedEdge;
                    }
                }
                values.resize(writeIndex);
                commitMaterializedEdges(node, values);
            }
            releaseEdgeMaterializationScratchIfComplete();
        }

        static int vertexId(int row, int col, int numVertexCols) {
            return (row * numVertexCols) + col;
        }

        struct OrientedGeometry {
            int startVertex = -1;
            int endVertex = -1;
            Direction direction = Direction::North;
            int signedArea2 = 0;
        };

        OrientedGeometry orientedGeometry(int packedEdge) const {
            const ContourTraceEdge edge = ContourTraceComputation::unpackEdge(packedEdge);
            const int cols = tree.getNumColsOfImage();
            const int numVertexCols = cols + 1;
            const auto [row, col] = ImageUtils::to2D(edge.pixel, cols);

            switch (edge.side) {
                case ContourTraceSide::North:
                    return {
                        vertexId(row, col, numVertexCols),
                        vertexId(row, col + 1, numVertexCols),
                        Direction::East,
                        -row};
                case ContourTraceSide::East:
                    return {
                        vertexId(row, col + 1, numVertexCols),
                        vertexId(row + 1, col + 1, numVertexCols),
                        Direction::South,
                        col + 1};
                case ContourTraceSide::South:
                    return {
                        vertexId(row + 1, col + 1, numVertexCols),
                        vertexId(row + 1, col, numVertexCols),
                        Direction::West,
                        row + 1};
                case ContourTraceSide::West:
                    return {
                        vertexId(row + 1, col, numVertexCols),
                        vertexId(row, col, numVertexCols),
                        Direction::North,
                        -col};
            }
            throw std::runtime_error("Invalid contour trace side.");
        }

        static int turnPriority(Direction incoming, Direction outgoing) {
            static constexpr std::array<int, 4> priorities{
                1, // straight
                0, // right
                3, // back
                2  // left
            };
            const int turn = (static_cast<int>(outgoing) - static_cast<int>(incoming) + 4) & 3;
            return priorities[static_cast<std::size_t>(turn)];
        }

        static int chooseNextEdge(
            const DirectedEdge& current,
            int outgoingHead,
            const std::vector<DirectedEdge>& directedEdges,
            const std::vector<int>& outgoingNext,
            const std::vector<uint32_t>& visitedGeneration,
            uint32_t visitGeneration
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
            ,
            std::size_t& candidateScans,
            std::size_t& visitedSkips,
            std::size_t& unvisitedCandidates
#endif
        ) {
            int best = -1;
            int bestPriority = std::numeric_limits<int>::max();
            for (int candidate = outgoingHead; candidate != -1; candidate = outgoingNext[static_cast<std::size_t>(candidate)]) {
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
                ++candidateScans;
#endif
                if (visitedGeneration[static_cast<std::size_t>(candidate)] == visitGeneration) {
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
                    ++visitedSkips;
#endif
                    continue;
                }
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
                ++unvisitedCandidates;
#endif
                const int priority = turnPriority(
                    current.direction,
                    directedEdges[static_cast<std::size_t>(candidate)].direction);
                if (priority < bestPriority ||
                    (priority == bestPriority && directedEdges[static_cast<std::size_t>(candidate)].packedEdge < directedEdges[static_cast<std::size_t>(best)].packedEdge)) {
                    best = candidate;
                    bestPriority = priority;
                }
            }
            return best;
        }

        void resetTraceOutgoingHeads() const {
            for (int vertex : traceTouchedVertices_) {
                const auto vertexIndex = static_cast<std::size_t>(vertex);
                traceOutgoingHead_[vertexIndex] = -1;
            }
            traceTouchedVertices_.clear();
        }

        [[nodiscard]] std::size_t imageVertexCount() const {
            return static_cast<std::size_t>((tree.getNumRowsOfImage() + 1) * (tree.getNumColsOfImage() + 1));
        }

        void ensureTraceOutgoingHeadStorage() const {
            if (!traceOutgoingHead_.empty()) {
                return;
            }
            traceOutgoingHead_.assign(imageVertexCount(), -1);
        }

        static std::size_t nextPowerOfTwoAtLeast(std::size_t value) {
            std::size_t result = 1;
            while (result < value) {
                result <<= 1;
            }
            return result;
        }

        static std::size_t sparseOutgoingTableSize(std::size_t edgeCount) {
            const std::size_t target = std::max<std::size_t>(2, (2 * edgeCount) + 1);
            return nextPowerOfTwoAtLeast(target);
        }

        static constexpr std::size_t sparseNodeLocalTraceLimit() {
            // Keep sparse adjacency for interactive point queries, then switch
            // to dense storage before random all-node access pays hash overhead.
            return 8;
        }

        [[nodiscard]] bool shouldUseSparseTraceAdjacency(
            std::size_t edgeCount,
            bool streamLoopInfosDirectly) const {
            if (streamLoopInfosDirectly || !traceOutgoingHead_.empty()) {
                return false;
            }
            if (nodeLocalLoopTraceCount_ >= sparseNodeLocalTraceLimit()) {
                return false;
            }
            const std::size_t sparseBytes =
                sparseOutgoingTableSize(edgeCount) *
                ((2 * sizeof(int)) + sizeof(uint32_t));
            const std::size_t denseBytes = imageVertexCount() * sizeof(int);
            return sparseBytes < denseBytes;
        }

        void prepareSparseTraceOutgoingHeads(std::size_t edgeCount) const {
            const std::size_t tableSize = sparseOutgoingTableSize(edgeCount);
            if (traceSparseVertexKeys_.size() != tableSize) {
                traceSparseVertexKeys_.assign(tableSize, 0);
                traceSparseOutgoingHead_.assign(tableSize, -1);
                traceSparseSlotGeneration_.assign(tableSize, 0);
            }
            traceSparseTouchedSlots_.clear();

            ++traceSparseGeneration_;
            if (traceSparseGeneration_ == 0) {
                std::fill(traceSparseSlotGeneration_.begin(), traceSparseSlotGeneration_.end(), 0);
                traceSparseGeneration_ = 1;
            }
        }

        [[nodiscard]] std::size_t sparseTraceSlot(int vertex) const {
            const auto hash = static_cast<uint32_t>(vertex) * uint32_t{2654435761u};
            return static_cast<std::size_t>(hash) & (traceSparseVertexKeys_.size() - 1);
        }

        [[nodiscard]] int sparseTraceOutgoingHead(int vertex) const {
            std::size_t slot = sparseTraceSlot(vertex);
            while (traceSparseSlotGeneration_[slot] == traceSparseGeneration_) {
                if (traceSparseVertexKeys_[slot] == vertex) {
                    return traceSparseOutgoingHead_[slot];
                }
                slot = (slot + 1) & (traceSparseVertexKeys_.size() - 1);
            }
            return -1;
        }

        [[nodiscard]] std::size_t findOrInsertSparseTraceVertex(int vertex) const {
            std::size_t slot = sparseTraceSlot(vertex);
            while (traceSparseSlotGeneration_[slot] == traceSparseGeneration_) {
                if (traceSparseVertexKeys_[slot] == vertex) {
                    return slot;
                }
                slot = (slot + 1) & (traceSparseVertexKeys_.size() - 1);
            }

            traceSparseSlotGeneration_[slot] = traceSparseGeneration_;
            traceSparseVertexKeys_[slot] = vertex;
            traceSparseOutgoingHead_[slot] = -1;
            traceSparseTouchedSlots_.push_back(static_cast<int>(slot));
            return slot;
        }

        template <TraceAdjacencyMode Mode>
        [[nodiscard]] int traceOutgoingHead(int vertex) const {
            if constexpr (Mode == TraceAdjacencyMode::Sparse) {
                return sparseTraceOutgoingHead(vertex);
            } else {
                return traceOutgoingHead_[static_cast<std::size_t>(vertex)];
            }
        }

        template <TraceAdjacencyMode Mode>
        void prepareTraceAdjacency(std::size_t edgeCount) const {
            if constexpr (Mode == TraceAdjacencyMode::Sparse) {
                prepareSparseTraceOutgoingHeads(edgeCount);
            } else {
                if (!traceSparseVertexKeys_.empty()) {
                    releaseSparseTraceOutgoingScratch();
                }
                ensureTraceOutgoingHeadStorage();
                resetTraceOutgoingHeads();
            }
        }

        template <TraceAdjacencyMode Mode>
        void resetTraceAdjacency() const {
            if constexpr (Mode == TraceAdjacencyMode::Dense) {
                resetTraceOutgoingHeads();
            }
        }

        template <TraceAdjacencyMode Mode>
        void pushTraceOutgoingEdge(
            int startVertex,
            int edgeIndex) const {
            if constexpr (Mode == TraceAdjacencyMode::Sparse) {
                const std::size_t slot = findOrInsertSparseTraceVertex(startVertex);
                traceOutgoingNext_.push_back(traceSparseOutgoingHead_[slot]);
                traceSparseOutgoingHead_[slot] = edgeIndex;
            } else {
                const auto vertexIndex = static_cast<std::size_t>(startVertex);
                if (traceOutgoingHead_[vertexIndex] == -1) {
                    traceTouchedVertices_.push_back(startVertex);
                }
                traceOutgoingNext_.push_back(traceOutgoingHead_[vertexIndex]);
                traceOutgoingHead_[vertexIndex] = edgeIndex;
            }
        }

        void nextTraceVisitGeneration(std::size_t edgeCount) const {
            if (traceVisitedGeneration_.size() < edgeCount) {
                traceVisitedGeneration_.resize(edgeCount, 0);
            }
            ++traceVisitGeneration_;
            if (traceVisitGeneration_ == 0) {
                std::fill(traceVisitedGeneration_.begin(), traceVisitedGeneration_.end(), 0);
                traceVisitGeneration_ = 1;
            }
        }

        [[nodiscard]] bool isTraceVisited(int edgeIndex) const {
            return traceVisitedGeneration_[static_cast<std::size_t>(edgeIndex)] == traceVisitGeneration_;
        }

        void markTraceVisited(int edgeIndex) const {
            traceVisitedGeneration_[static_cast<std::size_t>(edgeIndex)] = traceVisitGeneration_;
        }

        void reserveLoopInfoCapacityForGlobalTrace(NodeId root) const {
            std::size_t pendingNonEmptyNodes = 0;
            std::size_t pendingEdgeCount = 0;
            for (NodeId node : tree.getNodeSubtree(root)) {
                if (cachedLoopReady_[static_cast<std::size_t>(node)]) {
                    continue;
                }
                const auto edgeCount = static_cast<std::size_t>(cachedEdgeSize_[static_cast<std::size_t>(node)]);
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
            const std::size_t targetCapacity = cachedLoopInfos_.size() + additionalCapacity;
            if (cachedLoopInfos_.capacity() < targetCapacity) {
                cachedLoopInfos_.reserve(targetCapacity);
            }
        }

        bool canReuseEdgeSegmentForLoops(NodeId node) const {
            if (tree.isRoot(node)) {
                return true;
            }
            const NodeId parent = tree.getNodeParent(node);
            return parent == InvalidNode ||
                   parent == node ||
                   cachedEdgeReady_[static_cast<std::size_t>(parent)] != 0;
        }

        void ensureNodeLoopsMaterialized(NodeId node) const {
            requireStableTree("IncrementalContourTraces::ensureNodeLoopsMaterialized");
            requireLiveTraceNode(node, "IncrementalContourTraces::ensureNodeLoopsMaterialized");
            if (cachedLoopReady_[static_cast<std::size_t>(node)]) {
                return;
            }
            ensureEdgesMaterialized(node);
            traceNodeLoops(node, false);
        }

        void ensureLoopsMaterialized(NodeId root) const {
            requireStableTree("IncrementalContourTraces::ensureLoopsMaterialized");
            requireLiveTraceNode(root, "IncrementalContourTraces::ensureLoopsMaterialized");
            ensureEdgesMaterialized(root);
            reserveLoopInfoCapacityForGlobalTrace(root);

            std::vector<std::pair<NodeId, bool>> stack;
            stack.emplace_back(root, false);

            while (!stack.empty()) {
                const auto [node, expanded] = stack.back();
                stack.pop_back();
                if (!expanded) {
                    stack.emplace_back(node, true);
                    for (NodeId child : tree.getChildren(node)) {
                        if (!cachedLoopReady_[static_cast<std::size_t>(child)]) {
                            stack.emplace_back(child, false);
                        }
                    }
                    continue;
                }
                if (!cachedLoopReady_[static_cast<std::size_t>(node)]) {
                    traceNodeLoops(node, true);
                }
            }
        }

        void traceNodeLoops(NodeId node, bool streamLoopInfosDirectly) const {
            const std::size_t edgeCount = static_cast<std::size_t>(cachedEdgeSize_[static_cast<std::size_t>(node)]);
            const bool useSparseAdjacency = shouldUseSparseTraceAdjacency(edgeCount, streamLoopInfosDirectly);
            if (!streamLoopInfosDirectly) {
                ++nodeLocalLoopTraceCount_;
            }
            if (useSparseAdjacency) {
                traceNodeLoopsWithAdjacency<TraceAdjacencyMode::Sparse>(
                    node,
                    streamLoopInfosDirectly,
                    edgeCount);
            } else {
                traceNodeLoopsWithAdjacency<TraceAdjacencyMode::Dense>(
                    node,
                    streamLoopInfosDirectly,
                    edgeCount);
            }
        }

        template <TraceAdjacencyMode Mode>
        void traceNodeLoopsWithAdjacency(
            NodeId node,
            bool streamLoopInfosDirectly,
            std::size_t edgeCount) const {
            prepareTraceAdjacency<Mode>(edgeCount);
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
            const auto buildAdjacencyStart = traceProfileNow();
            std::size_t nodeOutgoingVertices = 0;
            std::size_t nodeSingleOutgoingVertices = 0;
            std::size_t nodeMultiOutgoingVertices = 0;
            std::size_t nodeMaxOutgoingDegree = 0;
            std::size_t nodeClosedLoopStops = 0;
            std::size_t nodeMissingOutgoingStops = 0;
            std::size_t nodeSingleSuccessorSteps = 0;
            std::size_t nodeSingleSuccessorVisitedStops = 0;
            std::size_t nodeAmbiguousSuccessorSteps = 0;
            std::size_t nodeAmbiguousSuccessorDeadEnds = 0;
            std::size_t nodeSuccessorCandidateScans = 0;
            std::size_t nodeSuccessorVisitedSkips = 0;
            std::size_t nodeSuccessorUnvisitedCandidates = 0;
#endif
            traceDirectedEdges_.clear();
            traceDirectedEdges_.reserve(edgeCount);
            traceOutgoingNext_.clear();
            traceOutgoingNext_.reserve(edgeCount);

            for (auto it = cachedEdgeBegin(node); it != cachedEdgeEnd(node); ++it) {
                const OrientedGeometry geometry = orientedGeometry(*it);
                const int edgeIndex = static_cast<int>(traceDirectedEdges_.size());
                traceDirectedEdges_.emplace_back(
                    *it,
                    geometry.startVertex,
                    geometry.endVertex,
                    geometry.direction,
                    geometry.signedArea2);
                pushTraceOutgoingEdge<Mode>(geometry.startVertex, edgeIndex);
            }
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
            const auto buildAdjacencyEnd = traceProfileNow();
            const auto profileCountersStart = traceProfileNow();
            if constexpr (Mode == TraceAdjacencyMode::Dense) {
                for (int vertex : traceTouchedVertices_) {
                    std::size_t degree = 0;
                    for (int candidate = traceOutgoingHead_[static_cast<std::size_t>(vertex)];
                         candidate != -1;
                         candidate = traceOutgoingNext_[static_cast<std::size_t>(candidate)]) {
                        ++degree;
                    }
                    if (degree == 0) {
                        continue;
                    }
                    ++nodeOutgoingVertices;
                    if (degree == 1) {
                        ++nodeSingleOutgoingVertices;
                    } else {
                        ++nodeMultiOutgoingVertices;
                    }
                    nodeMaxOutgoingDegree = std::max(nodeMaxOutgoingDegree, degree);
                }
            } else {
                for (int touchedSlot : traceSparseTouchedSlots_) {
                    std::size_t degree = 0;
                    const auto slot = static_cast<std::size_t>(touchedSlot);
                    for (int candidate = traceSparseOutgoingHead_[slot];
                         candidate != -1;
                         candidate = traceOutgoingNext_[static_cast<std::size_t>(candidate)]) {
                        ++degree;
                    }
                    if (degree == 0) {
                        continue;
                    }
                    ++nodeOutgoingVertices;
                    if (degree == 1) {
                        ++nodeSingleOutgoingVertices;
                    } else {
                        ++nodeMultiOutgoingVertices;
                    }
                    nodeMaxOutgoingDegree = std::max(nodeMaxOutgoingDegree, degree);
                }
            }
            const auto profileCountersEnd = traceProfileNow();
#endif

            nextTraceVisitGeneration(traceDirectedEdges_.size());
            traceNodeLoopEdges_.clear();
            traceNodeLoopEdges_.reserve(traceDirectedEdges_.size());
            traceNodeLoops_.clear();
            if (!streamLoopInfosDirectly) {
                traceNodeLoops_.reserve((traceDirectedEdges_.size() / 4) + 1);
            }

            std::size_t globalEdgeOffset = static_cast<std::size_t>(cachedEdgeOffset_[static_cast<std::size_t>(node)]);
            const bool reuseEdgeSegment = canReuseEdgeSegmentForLoops(node);
            if (!reuseEdgeSegment) {
                globalEdgeOffset = cachedEdgeValues_.size();
            }
            const std::size_t loopInfoOffset = cachedLoopInfos_.size();

#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
            const auto walkLoopsStart = traceProfileNow();
#endif
            try {
                for (int startEdge = 0; startEdge < static_cast<int>(traceDirectedEdges_.size()); ++startEdge) {
                    if (isTraceVisited(startEdge)) {
                        continue;
                    }

                    const int loopStartVertex = traceDirectedEdges_[static_cast<std::size_t>(startEdge)].startVertex;
                    const std::size_t loopEdgeOffset = traceNodeLoopEdges_.size();
                    int signedArea2 = 0;
                    int current = startEdge;

                    // Successor selection only returns unvisited edges, so the
                    // hot loop does not need a second visited check here.
                    while (current != -1) {
                        const DirectedEdge& edge = traceDirectedEdges_[static_cast<std::size_t>(current)];
                        markTraceVisited(current);
                        traceNodeLoopEdges_.push_back(edge.packedEdge);
                        signedArea2 += edge.signedArea2;

                        if (edge.endVertex == loopStartVertex) {
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
                            ++nodeClosedLoopStops;
#endif
                            break;
                        }

                        const int outgoingHead = traceOutgoingHead<Mode>(edge.endVertex);
                        if (outgoingHead == -1) {
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
                            ++nodeMissingOutgoingStops;
#endif
                            break;
                        }
                        if (traceOutgoingNext_[static_cast<std::size_t>(outgoingHead)] == -1) {
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
                            ++nodeSingleSuccessorSteps;
#endif
                            if (isTraceVisited(outgoingHead)) {
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
                                ++nodeSingleSuccessorVisitedStops;
#endif
                                current = -1;
                            } else {
                                current = outgoingHead;
                            }
                        } else {
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
                            ++nodeAmbiguousSuccessorSteps;
#endif
                            current = chooseNextEdge(
                                edge,
                                outgoingHead,
                                traceDirectedEdges_,
                                traceOutgoingNext_,
                                traceVisitedGeneration_,
                                traceVisitGeneration_
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
                                ,
                                nodeSuccessorCandidateScans,
                                nodeSuccessorVisitedSkips,
                                nodeSuccessorUnvisitedCandidates
#endif
                            );
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
                            if (current == -1) {
                                ++nodeAmbiguousSuccessorDeadEnds;
                            }
#endif
                        }
                    }

                    const std::size_t loopEdgeCount = traceNodeLoopEdges_.size() - loopEdgeOffset;
                    if (loopEdgeCount == 0) {
                        continue;
                    }

                    const ContourLoopKind kind =
                        signedArea2 >= 0 ? ContourLoopKind::External : ContourLoopKind::Internal;
                    if (streamLoopInfosDirectly) {
                        cachedLoopInfos_.push_back(ContourTraceLoop{
                            kind,
                            checkedU32(globalEdgeOffset + loopEdgeOffset, "global loop edge offset"),
                            checkedU32(loopEdgeCount, "loop edge count"),
                            signedArea2});
                    } else {
                        traceNodeLoops_.push_back(ContourTraceLoop{
                            kind,
                            checkedU32(loopEdgeOffset, "local loop edge offset"),
                            checkedU32(loopEdgeCount, "loop edge count"),
                            signedArea2});
                    }
                }
            } catch (...) {
                if (streamLoopInfosDirectly) {
                    cachedLoopInfos_.resize(loopInfoOffset);
                }
                resetTraceAdjacency<Mode>();
                throw;
            }
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
            const auto walkLoopsEnd = traceProfileNow();
            const auto resetOutgoingStart = traceProfileNow();
#endif
            resetTraceAdjacency<Mode>();
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
            const auto resetOutgoingEnd = traceProfileNow();
#endif

            if (traceNodeLoopEdges_.size() != edgeCount) {
                if (streamLoopInfosDirectly) {
                    cachedLoopInfos_.resize(loopInfoOffset);
                }
                throw std::runtime_error("Contour trace loop traversal did not cover every materialized edge.");
            }
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
            const auto commitEdgesStart = traceProfileNow();
#endif
            if (reuseEdgeSegment) {
                std::copy(
                    traceNodeLoopEdges_.begin(),
                    traceNodeLoopEdges_.end(),
                    cachedEdgeValues_.begin() + static_cast<std::ptrdiff_t>(globalEdgeOffset));
            } else {
                cachedEdgeValues_.insert(cachedEdgeValues_.end(), traceNodeLoopEdges_.begin(), traceNodeLoopEdges_.end());
            }
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
            const auto commitEdgesEnd = traceProfileNow();
            const auto commitLoopsStart = traceProfileNow();
#endif
            cachedLoopEdgeCount_ += traceNodeLoopEdges_.size();

            if (!streamLoopInfosDirectly) {
                for (ContourTraceLoop& loop : traceNodeLoops_) {
                    loop.edgeOffset = checkedU32(globalEdgeOffset + loop.edgeOffset, "global loop edge offset");
                }
                cachedLoopInfos_.insert(cachedLoopInfos_.end(), traceNodeLoops_.begin(), traceNodeLoops_.end());
            }

            cachedLoopInfoOffset_[static_cast<std::size_t>(node)] = checkedU32(loopInfoOffset, "loop info offset");
            cachedLoopInfoSize_[static_cast<std::size_t>(node)] =
                checkedU32(cachedLoopInfos_.size() - loopInfoOffset, "loop info size");
            const auto index = static_cast<std::size_t>(node);
            if (!cachedLoopReady_[index]) {
                cachedLoopReady_[index] = 1;
                ++cachedLoopReadyCount_;
            }
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
            const auto commitLoopsEnd = traceProfileNow();
            const auto releaseScratchStart = traceProfileNow();
#endif
            releaseTraceScratchIfComplete();
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
            const auto releaseScratchEnd = traceProfileNow();
            if (activeTraceProfile_ != nullptr) {
                activeTraceProfile_->nodesTraced += 1;
                activeTraceProfile_->edgesTraced += edgeCount;
                activeTraceProfile_->loopsTraced += cachedLoopInfos_.size() - loopInfoOffset;
                activeTraceProfile_->outgoingVertices += nodeOutgoingVertices;
                activeTraceProfile_->singleOutgoingVertices += nodeSingleOutgoingVertices;
                activeTraceProfile_->multiOutgoingVertices += nodeMultiOutgoingVertices;
                activeTraceProfile_->maxOutgoingDegree =
                    std::max(activeTraceProfile_->maxOutgoingDegree, nodeMaxOutgoingDegree);
                activeTraceProfile_->closedLoopStops += nodeClosedLoopStops;
                activeTraceProfile_->missingOutgoingStops += nodeMissingOutgoingStops;
                activeTraceProfile_->singleSuccessorSteps += nodeSingleSuccessorSteps;
                activeTraceProfile_->singleSuccessorVisitedStops += nodeSingleSuccessorVisitedStops;
                activeTraceProfile_->ambiguousSuccessorSteps += nodeAmbiguousSuccessorSteps;
                activeTraceProfile_->ambiguousSuccessorDeadEnds += nodeAmbiguousSuccessorDeadEnds;
                activeTraceProfile_->successorCandidateScans += nodeSuccessorCandidateScans;
                activeTraceProfile_->successorVisitedSkips += nodeSuccessorVisitedSkips;
                activeTraceProfile_->successorUnvisitedCandidates += nodeSuccessorUnvisitedCandidates;
                activeTraceProfile_->profileCountersNs +=
                    traceProfileElapsedNs(profileCountersStart, profileCountersEnd);
                activeTraceProfile_->buildAdjacencyNs += traceProfileElapsedNs(buildAdjacencyStart, buildAdjacencyEnd);
                activeTraceProfile_->walkLoopsNs += traceProfileElapsedNs(walkLoopsStart, walkLoopsEnd);
                activeTraceProfile_->resetOutgoingNs += traceProfileElapsedNs(resetOutgoingStart, resetOutgoingEnd);
                activeTraceProfile_->commitEdgesNs += traceProfileElapsedNs(commitEdgesStart, commitEdgesEnd);
                activeTraceProfile_->commitLoopsNs += traceProfileElapsedNs(commitLoopsStart, commitLoopsEnd);
                activeTraceProfile_->releaseScratchNs += traceProfileElapsedNs(releaseScratchStart, releaseScratchEnd);
            }
#endif
        }
    };

private:
    using PendingEdgeLists = detail::PendingPixelLists;

    struct ExtractedTraceDeltas {
        LocalTraceDeltas deltas;
        int capacityHint = 0;
    };

    static int neighborPixel(
        const MorphologicalTree& tree,
        int pixel,
        ContourTraceSide side) {
        const int rows = tree.getNumRowsOfImage();
        const int cols = tree.getNumColsOfImage();
        const auto [row, col] = ImageUtils::to2D(pixel, cols);

        switch (side) {
            case ContourTraceSide::North:
                return row == 0 ? -1 : ImageUtils::to1D(row - 1, col, cols);
            case ContourTraceSide::West:
                return col == 0 ? -1 : ImageUtils::to1D(row, col - 1, cols);
            case ContourTraceSide::East:
                return col == cols - 1 ? -1 : ImageUtils::to1D(row, col + 1, cols);
            case ContourTraceSide::South:
                return row == rows - 1 ? -1 : ImageUtils::to1D(row + 1, col, cols);
        }
        return -1;
    }

    [[nodiscard]] static ExtractedTraceDeltas extractTraceDeltasImpl(const MorphologicalTree& tree) {
        if (tree.getNumRowsOfImage() <= 0 || tree.getNumColsOfImage() <= 0) {
            throw std::invalid_argument("Contour trace extraction requires a non-empty image domain.");
        }
        if (!tree.isAlive(tree.getRoot())) {
            throw std::invalid_argument("Contour trace extraction requires a live tree root.");
        }

        const int numNodes = tree.getNumInternalNodeSlots();
        const int totalPixels = tree.getNumRowsOfImage() * tree.getNumColsOfImage();
        const int totalPackedEdges = 4 * totalPixels;
        const int capacityHint = std::max(totalPixels, 1);

        PendingEdgeLists localEdgeAdditions(numNodes, capacityHint);
        PendingEdgeLists localEdgeRemovals(numNodes, capacityHint);

        static constexpr std::array<ContourTraceSide, 4> sides{
            ContourTraceSide::North,
            ContourTraceSide::West,
            ContourTraceSide::East,
            ContourTraceSide::South};

        detail::traversePostOrder(
            tree,
            tree.getRoot(),
            [](NodeId) -> void {},
            [](NodeId, NodeId) -> void {},
            [&](NodeId nodeId) {
                for (int pixel : tree.getProperParts(nodeId)) {
                    for (ContourTraceSide side : sides) {
                        const int packedEdge = packEdge(pixel, side);
                        const int neighbor = neighborPixel(tree, pixel, side);
                        if (neighbor < 0) {
                            localEdgeAdditions.add(nodeId, packedEdge);
                            continue;
                        }

                        const NodeId entry = detail::properPartEntryNode(tree, pixel, neighbor);
                        if (entry == InvalidNode || entry == nodeId) {
                            continue;
                        }

                        localEdgeAdditions.add(nodeId, packedEdge);
                        localEdgeRemovals.add(entry, packedEdge);
                    }
                }
            });

        return {
            LocalTraceDeltas::fromPendingEdgeLists(localEdgeAdditions, localEdgeRemovals, totalPackedEdges),
            capacityHint
        };
    }

public:
    /**
     * @brief Extracts compact local boundary-edge additions/removals.
     */
    [[nodiscard]] static LocalTraceDeltas extractTraceDeltas(const MorphologicalTree& tree) {
        ExtractedTraceDeltas extracted = extractTraceDeltasImpl(tree);
        return std::move(extracted.deltas);
    }

    template<AltitudeValue T>
    /**
     * @brief Extracts compact local boundary-edge additions/removals from a weighted view.
     */
    [[nodiscard]] static LocalTraceDeltas extractTraceDeltas(const WeightedTreeView<T>& tree) {
        tree.requireTopologyUnchanged("ContourTraceComputation::extractTraceDeltas");
        return extractTraceDeltas(tree.topology());
    }

    /**
     * @brief Runs incremental side-level contour extraction and returns lazy traces.
     */
    [[nodiscard]] static IncrementalContourTraces extract(const MorphologicalTree& tree) {
        ExtractedTraceDeltas extracted = extractTraceDeltasImpl(tree);
        return IncrementalContourTraces(tree, std::move(extracted.deltas), extracted.capacityHint);
    }

    template<AltitudeValue T>
    /**
     * @brief Runs incremental side-level contour extraction on a weighted view.
     */
    [[nodiscard]] static IncrementalContourTraces extract(const WeightedTreeView<T>& tree) {
        tree.requireTopologyUnchanged("ContourTraceComputation::extract");
        return extract(tree.topology());
    }
};

} // namespace mmcfilters
