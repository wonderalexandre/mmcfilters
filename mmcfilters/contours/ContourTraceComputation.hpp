#pragma once

#include "../localAttributes/FiniteWindowLocalAttributeComputer.hpp"
#include "../trees/MorphologicalTree.hpp"
#include "../trees/ValuedMorphologicalTreeView.hpp"
#include "../trees/detail/TreeTraversalDetail.hpp"
#include "../utils/Common.hpp"
#include "../utils/Image.hpp"
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
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mmcfilters {

/**
 * @brief Oriented side of one support pixel used as a boundary trace edge.
 */
enum class ContourTraceSide : uint8_t { North = 0, West = 1, East = 2, South = 3 };

/**
 * @brief Geometric class of a traced boundary loop.
 */
enum class ContourLoopKind : uint8_t { External, Internal };

/**
 * @brief One unpacked boundary edge attached to a support pixel.
 */
struct ContourTraceEdge {
    /// Row-major support-pixel index incident to the boundary edge.
    PixelId pixel = InvalidPixel;
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
  private:
    /** @brief Defines the `LocalTraceDeltas` alias used by the component. */
    using LocalTraceDeltas = detail::ContourTraceDeltaStore;

  public:
    /**
     * @brief Packs one pixel-side edge into a compact integer id.
     *
     * @param pixel Pixel identifier.
     * @param side Side selected by the operation.
     * @return The packed pixel-side edge into a compact integer id.
     */
    [[nodiscard]] static int packEdge(PixelId pixel, ContourTraceSide side) { return (4 * pixel) + static_cast<int>(side); }

    /**
     * @brief Unpacks one compact edge id.
     *
     * @param packedEdge Packed boundary-edge identifier.
     * @return The unpacked compact edge id.
     */
    [[nodiscard]] static ContourTraceEdge unpackEdge(int packedEdge) {
        if (packedEdge < 0) {
            return {};
        }
        const int side = packedEdge & 3;
        return ContourTraceEdge{packedEdge / 4, static_cast<ContourTraceSide>(side)};
    }

    /**
     * @brief Lazy result for side-level contour traces.
     */
    struct IncrementalContourTraces {
      private:
        friend class ContourTraceComputation;

        /** @brief Enumerates the supported direction values. */
        enum class Direction : uint8_t { North = 0, East = 1, South = 2, West = 3 };

        /** @brief Enumerates the supported trace adjacency mode values. */
        enum class TraceAdjacencyMode : uint8_t { Dense, Sparse };

        /** @brief Stores one oriented contour edge and its accumulated geometry. */
        struct DirectedEdge {
            /** @brief Packed edge. */
            int packedEdge = -1;
            /** @brief Start vertex. */
            int startVertex = -1;
            /** @brief End vertex. */
            int endVertex = -1;
            /** @brief Direction. */
            Direction direction = Direction::North;
            /** @brief Signed area2. */
            int signedArea2 = 0;

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
             * @param signedArea2 Twice the signed area accumulated for the loop.
             */
            DirectedEdge(int packedEdge, int startVertex, int endVertex, Direction direction, int signedArea2)
                : packedEdge(packedEdge), startVertex(startVertex), endVertex(endVertex), direction(direction), signedArea2(signedArea2) {}
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

        /** @brief References the tree used by the component. */
        const MorphologicalTree& tree;
        /** @brief Tree mutation version used to detect stale derived state. */
        std::size_t treeMutationVersion_ = 0;
        /** @brief Local deltas. */
        mutable LocalTraceDeltas localDeltas_;

        /** @brief Cached edge values buffer. */
        mutable std::vector<int> cachedEdgeValues_;
        /** @brief Cached edge offset buffer. */
        mutable std::vector<uint32_t> cachedEdgeOffset_;
        /** @brief Cached edge size buffer. */
        mutable std::vector<uint32_t> cachedEdgeSize_;
        /** @brief Cached edge ready buffer. */
        mutable std::vector<uint8_t> cachedEdgeReady_;
        /** @brief Cached edge ready count. */
        mutable std::size_t cachedEdgeReadyCount_ = 0;

        /** @brief Cached loop infos buffer. */
        mutable std::vector<ContourTraceLoop> cachedLoopInfos_;
        /** @brief Cached loop info offset buffer. */
        mutable std::vector<uint32_t> cachedLoopInfoOffset_;
        /** @brief Cached loop info size buffer. */
        mutable std::vector<uint32_t> cachedLoopInfoSize_;
        /** @brief Cached loop ready buffer. */
        mutable std::vector<uint8_t> cachedLoopReady_;
        /** @brief Cached loop ready count. */
        mutable std::size_t cachedLoopReadyCount_ = 0;
        /** @brief Cached loop edge count. */
        mutable std::size_t cachedLoopEdgeCount_ = 0;

        /** @brief Edge mark buffer. */
        mutable std::vector<uint16_t> edgeMark_;
        /** @brief Mark generation. */
        mutable uint16_t markGeneration_ = 1;
        /** @brief Indicates whether edge-materialization scratch storage was released. */
        mutable bool edgeMaterializationScratchReleased_ = false;

        /** @brief Trace directed edges buffer. */
        mutable std::vector<DirectedEdge> traceDirectedEdges_;
        /** @brief Trace outgoing head buffer. */
        mutable std::vector<int> traceOutgoingHead_;
        /** @brief Trace outgoing next buffer. */
        mutable std::vector<int> traceOutgoingNext_;
        /** @brief Trace touched vertices buffer. */
        mutable std::vector<int> traceTouchedVertices_;
        /** @brief Trace sparse vertex keys buffer. */
        mutable std::vector<int> traceSparseVertexKeys_;
        /** @brief Trace sparse outgoing head buffer. */
        mutable std::vector<int> traceSparseOutgoingHead_;
        /** @brief Trace sparse slot generation buffer. */
        mutable std::vector<uint32_t> traceSparseSlotGeneration_;
        /** @brief Trace sparse touched slots buffer. */
        mutable std::vector<int> traceSparseTouchedSlots_;
        /** @brief Trace sparse generation. */
        mutable uint32_t traceSparseGeneration_ = 1;
        /** @brief Node local loop trace count. */
        mutable std::size_t nodeLocalLoopTraceCount_ = 0;
        /** @brief Trace visited generation buffer. */
        mutable std::vector<uint32_t> traceVisitedGeneration_;
        /** @brief Trace visit generation. */
        mutable uint32_t traceVisitGeneration_ = 1;
        /** @brief Trace node loop edges buffer. */
        mutable std::vector<int> traceNodeLoopEdges_;
        /** @brief Trace node loops buffer. */
        mutable std::vector<ContourTraceLoop> traceNodeLoops_;
        /** @brief Indicates whether trace-construction scratch storage was released. */
        mutable bool traceScratchReleased_ = false;

        /**
         * @brief Constructs `IncrementalContourTraces` from the supplied inputs.
         *
         * @param tree Tree topology.
         * @param localDeltas Local contour deltas to append to the trace store.
         * @param capacityHint Estimated number of entries used to reserve storage.
         */
        IncrementalContourTraces(const MorphologicalTree& tree, LocalTraceDeltas localDeltas, int capacityHint)
            : tree(tree), treeMutationVersion_(tree.getMutationVersion()), localDeltas_(std::move(localDeltas)),
              cachedEdgeOffset_(static_cast<std::size_t>(tree.numInternalNodeSlots()), 0),
              cachedEdgeSize_(static_cast<std::size_t>(tree.numInternalNodeSlots()), 0),
              cachedEdgeReady_(static_cast<std::size_t>(tree.numInternalNodeSlots()), 0),
              cachedLoopInfoOffset_(static_cast<std::size_t>(tree.numInternalNodeSlots()), 0),
              cachedLoopInfoSize_(static_cast<std::size_t>(tree.numInternalNodeSlots()), 0),
              cachedLoopReady_(static_cast<std::size_t>(tree.numInternalNodeSlots()), 0),
              edgeMark_(static_cast<std::size_t>(4 * tree.numRows() * tree.numColumns()), 0) {
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
                 * @param values Values read or written by the operation.
                 * @param index Zero-based index.
                 */
                iterator(const std::vector<int>* values, std::size_t index) : values_(values), index_(index) {}

                /**
                 * @brief Returns the unpacked edge at the current position.
                 *
                 * @return The unpacked edge at the current position.
                 */
                value_type operator*() const { return ContourTraceComputation::unpackEdge((*values_)[index_]); }

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
                    iterator tmp(*this);
                    ++(*this);
                    return tmp;
                }

                /// Compares the backing buffer and position.
                friend bool operator==(const iterator& lhs, const iterator& rhs) { return lhs.values_ == rhs.values_ && lhs.index_ == rhs.index_; }

                /// Returns true when two iterator positions differ.
                friend bool operator!=(const iterator& lhs, const iterator& rhs) { return !(lhs == rhs); }

              private:
                /** @brief Values buffer. */
                const std::vector<int>* values_ = nullptr;
                /** @brief Index. */
                std::size_t index_ = 0;
            };

            /**
             * @brief Creates an empty edge range.
             */
            EdgeRange() = default;

            /**
             * @brief Creates a view over `size` packed edges starting at `offset`.
             *
             * @param values Values read or written by the operation.
             * @param offset Offset into the underlying storage.
             * @param size Number.
             */
            EdgeRange(const std::vector<int>* values, std::size_t offset, std::size_t size) : values_(values), offset_(offset), size_(size) {}

            /**
             * @brief Returns an iterator to the first edge.
             *
             * @return An iterator to the first edge.
             */
            iterator begin() const { return iterator(values_, offset_); }

            /**
             * @brief Returns the exclusive end iterator.
             *
             * @return The exclusive end iterator.
             */
            iterator end() const { return iterator(values_, offset_ + size_); }

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
            /** @brief Values buffer. */
            const std::vector<int>* values_ = nullptr;
            /** @brief Offset. */
            std::size_t offset_ = 0;
            /** @brief Size. */
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
                localDeltas_.addValues.capacity() * sizeof(int) + localDeltas_.removeValues.capacity() * sizeof(int) +
                localDeltas_.addSpans.capacity() * sizeof(LocalTraceDeltas::Span) + localDeltas_.removeSpans.capacity() * sizeof(LocalTraceDeltas::Span) +
                cachedEdgeValues_.capacity() * sizeof(int) + cachedEdgeOffset_.capacity() * sizeof(uint32_t) + cachedEdgeSize_.capacity() * sizeof(uint32_t) +
                cachedEdgeReady_.capacity() * sizeof(uint8_t) + cachedLoopInfos_.capacity() * sizeof(ContourTraceLoop) +
                cachedLoopInfoOffset_.capacity() * sizeof(uint32_t) + cachedLoopInfoSize_.capacity() * sizeof(uint32_t) +
                cachedLoopReady_.capacity() * sizeof(uint8_t) + edgeMark_.capacity() * sizeof(uint16_t) +
                traceDirectedEdges_.capacity() * sizeof(DirectedEdge) + traceOutgoingHead_.capacity() * sizeof(int) +
                traceOutgoingNext_.capacity() * sizeof(int) + traceTouchedVertices_.capacity() * sizeof(int) + traceSparseVertexKeys_.capacity() * sizeof(int) +
                traceSparseOutgoingHead_.capacity() * sizeof(int) + traceSparseSlotGeneration_.capacity() * sizeof(uint32_t) +
                traceSparseTouchedSlots_.capacity() * sizeof(int) + traceVisitedGeneration_.capacity() * sizeof(uint32_t) +
                traceNodeLoopEdges_.capacity() * sizeof(int) + traceNodeLoops_.capacity() * sizeof(ContourTraceLoop);
            return stats;
        }
#endif

        /**
         * @brief Returns unordered materialized boundary edges for one node.
         *
         * @param node Node identifier.
         * @return Unordered materialized boundary edges for one node.
         */
        [[nodiscard]] EdgeRange getEdges(NodeId node) const {
            requireStableTree("IncrementalContourTraces::getEdges");
            requireLiveTraceNode(node, "IncrementalContourTraces::getEdges");
            ensureEdgesMaterialized(node);
            return EdgeRange(&cachedEdgeValues_, cachedEdgeOffset_[static_cast<std::size_t>(node)], cachedEdgeSize_[static_cast<std::size_t>(node)]);
        }

        /**
         * @brief Returns an owning copy of the loop metadata for one node.
         *
         * The returned vector remains valid when later lazy queries materialize
         * loops for other nodes.
         *
         * @param node Node identifier.
         * @return An owning copy of the loop metadata for one node.
         */
        [[nodiscard]] std::vector<ContourTraceLoop> getLoops(NodeId node) const {
            requireStableTree("IncrementalContourTraces::getLoops");
            requireLiveTraceNode(node, "IncrementalContourTraces::getLoops");
            ensureNodeLoopsMaterialized(node);
            const auto offset = static_cast<std::size_t>(cachedLoopInfoOffset_[static_cast<std::size_t>(node)]);
            const auto size = static_cast<std::size_t>(cachedLoopInfoSize_[static_cast<std::size_t>(node)]);
            if (size == 0) {
                return {};
            }
            return std::vector<ContourTraceLoop>(cachedLoopInfos_.begin() + static_cast<std::ptrdiff_t>(offset),
                                                 cachedLoopInfos_.begin() + static_cast<std::ptrdiff_t>(offset + size));
        }

        /**
         * @brief Returns the ordered edges belonging to one loop.
         *
         * @param loop Contour loop descriptor.
         * @return The ordered edges belonging to one loop.
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
            ensureLoopsMaterialized(tree.root());
        }

        /**
         * @brief Returns whether loop traces are materialized for every live node.
         *
         * @return Whether loop traces are materialized for every live node.
         */
        [[nodiscard]] bool isMaterialized() const {
            requireStableTree("IncrementalContourTraces::isMaterialized");
            for (NodeId node : tree.aliveNodeIds()) {
                if (!cachedLoopReady_[static_cast<std::size_t>(node)]) {
                    return false;
                }
            }
            return true;
        }

        /**
         * @brief Returns whether packed boundary edges are materialized for `node`.
         *
         * @param node Node identifier.
         * @return Whether packed boundary edges are materialized for node.
         */
        [[nodiscard]] bool isEdgeMaterialized(NodeId node) const {
            requireStableTree("IncrementalContourTraces::isEdgeMaterialized");
            requireLiveTraceNode(node, "IncrementalContourTraces::isEdgeMaterialized");
            return static_cast<bool>(cachedEdgeReady_[static_cast<std::size_t>(node)]);
        }

        /**
         * @brief Returns whether ordered loops are materialized for `node`.
         *
         * @param node Node identifier.
         * @return Whether ordered loops are materialized for node.
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
                ensureLoopsMaterialized(tree.root());
            } catch (...) {
                activeTraceProfile_ = previousProfile;
                throw;
            }
            activeTraceProfile_ = previousProfile;
            return stats;
        }
#endif

      private:
        /**
         * @brief Validates stable tree.
         *
         * @param context Operation name used in diagnostics.
         */
        void requireStableTree(const char* context) const { tree.requireMutationVersion(treeMutationVersion_, context); }

        /**
         * @brief Validates live trace node.
         *
         * @param node Node identifier.
         * @param context Operation name used in diagnostics.
         */
        void requireLiveTraceNode(NodeId node, const char* context) const {
            if (!tree.isAlive(node)) {
                throw std::invalid_argument(std::string(context) + " requires a live internal NodeId.");
            }
        }

        /**
         * @brief Checks and converts u32.
         *
         * @param value Value.
         * @param context Operation name used in diagnostics.
         * @return Owned native hierarchy storage.
         */
        static uint32_t checkedU32(std::size_t value, const char* context) {
            if (value > static_cast<std::size_t>(std::numeric_limits<uint32_t>::max())) {
                throw std::overflow_error(std::string(context) + " exceeds uint32_t.");
            }
            return static_cast<uint32_t>(value);
        }

        /**
         * @brief Releases vector.
         *
         * @param values Values read or written by the operation.
         */
        template <class T> static void releaseVector(std::vector<T>& values) { std::vector<T>().swap(values); }

#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
        using TraceProfileClock = std::chrono::steady_clock;

        [[nodiscard]] static TraceProfileClock::time_point traceProfileNow() { return TraceProfileClock::now(); }

        static std::int64_t traceProfileElapsedNs(TraceProfileClock::time_point start, TraceProfileClock::time_point end) {
            return std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        }
#endif

        /**
         * @brief Checks whether the edge representation is cached for every tree node.
         *
         * @return True when the documented condition holds; otherwise false.
         */
        [[nodiscard]] bool allEdgesMaterialized() const { return cachedEdgeReadyCount_ == static_cast<std::size_t>(tree.numNodes()); }

        /**
         * @brief Checks whether the loop representation is cached for every tree node.
         *
         * @return True when the documented condition holds; otherwise false.
         */
        [[nodiscard]] bool allLoopsMaterialized() const { return cachedLoopReadyCount_ == static_cast<std::size_t>(tree.numNodes()); }

        /**
         * @brief Releases edge materialization scratch if complete.
         */
        void releaseEdgeMaterializationScratchIfComplete() const {
            if (edgeMaterializationScratchReleased_ || !allEdgesMaterialized()) {
                return;
            }

            localDeltas_ = LocalTraceDeltas{};
            releaseVector(edgeMark_);
            markGeneration_ = 1;
            edgeMaterializationScratchReleased_ = true;
        }

        /**
         * @brief Releases sparse trace outgoing scratch.
         */
        void releaseSparseTraceOutgoingScratch() const {
            releaseVector(traceSparseVertexKeys_);
            releaseVector(traceSparseOutgoingHead_);
            releaseVector(traceSparseSlotGeneration_);
            releaseVector(traceSparseTouchedSlots_);
            traceSparseGeneration_ = 1;
        }

        /**
         * @brief Releases trace scratch if complete.
         */
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

        /**
         * @brief Advances mark generation.
         */
        void nextMarkGeneration() const {
            ++markGeneration_;
            if (markGeneration_ == 0) {
                std::fill(edgeMark_.begin(), edgeMark_.end(), 0);
                markGeneration_ = 1;
            }
        }

        /**
         * @brief Adds if unmarked.
         *
         * @param values Values read or written by the operation.
         * @param packedEdge Packed boundary-edge identifier.
         */
        void addIfUnmarked(std::vector<int>& values, int packedEdge) const {
            if (packedEdge < 0 || packedEdge >= static_cast<int>(edgeMark_.size())) {
                return;
            }
            if (edgeMark_[static_cast<std::size_t>(packedEdge)] != markGeneration_) {
                edgeMark_[static_cast<std::size_t>(packedEdge)] = markGeneration_;
                values.push_back(packedEdge);
            }
        }

        /**
         * @brief Removes if marked.
         *
         * @param packedEdge Packed boundary-edge identifier.
         */
        void removeIfMarked(int packedEdge) const {
            if (packedEdge >= 0 && packedEdge < static_cast<int>(edgeMark_.size())) {
                edgeMark_[static_cast<std::size_t>(packedEdge)] = 0;
            }
        }

        /**
         * @brief Returns an iterator to the first cached edge of a node.
         *
         * @param node Node identifier.
         * @return Values produced by the operation.
         */
        std::vector<int>::const_iterator cachedEdgeBegin(NodeId node) const {
            return cachedEdgeValues_.begin() + static_cast<std::ptrdiff_t>(cachedEdgeOffset_[static_cast<std::size_t>(node)]);
        }

        /**
         * @brief Returns an iterator past the last cached edge of a node.
         *
         * @param node Node identifier.
         * @return Values produced by the operation.
         */
        std::vector<int>::const_iterator cachedEdgeEnd(NodeId node) const {
            return cachedEdgeBegin(node) + static_cast<std::ptrdiff_t>(cachedEdgeSize_[static_cast<std::size_t>(node)]);
        }

        /**
         * @brief Materializes one node's boundary-edge sequence in the shared cache.
         *
         * @param node Node identifier.
         * @param values Values read or written by the operation.
         */
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

        /**
         * @brief Ensures edges materialized.
         *
         * @param root Root node of the operation.
         */
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
                    for (NodeId child : tree.children(node)) {
                        if (!cachedEdgeReady_[static_cast<std::size_t>(child)]) {
                            stack.emplace_back(child, false);
                        }
                    }
                    continue;
                }

                values.clear();
                const auto additions = localDeltas_.additions(node);
                std::size_t reserveSize = additions.size();
                for (NodeId child : tree.children(node)) {
                    reserveSize += static_cast<std::size_t>(cachedEdgeSize_[static_cast<std::size_t>(child)]);
                }
                values.reserve(reserveSize);
                nextMarkGeneration();

                for (NodeId child : tree.children(node)) {
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

        /**
         * @brief Computes the row-major identifier of a grid vertex.
         *
         * @param row Zero-based row coordinate.
         * @param column Zero-based column coordinate.
         * @param numVertexColumns Count.
         * @return Row-major vertex identifier.
         */
        static int vertexId(int row, int column, int numVertexColumns) { return (row * numVertexColumns) + column; }

        /** @brief Stores signed geometric measurements for one oriented contour. */
        struct OrientedGeometry {
            /** @brief Start vertex. */
            int startVertex = -1;
            /** @brief End vertex. */
            int endVertex = -1;
            /** @brief Direction. */
            Direction direction = Direction::North;
            /** @brief Signed area2. */
            int signedArea2 = 0;
        };

        /**
         * @brief Computes the oriented geometry accumulated by a traced loop.
         *
         * @param packedEdge Packed boundary-edge identifier.
         * @return Oriented endpoints, direction, and signed-area contribution.
         */
        OrientedGeometry orientedGeometry(int packedEdge) const {
            const ContourTraceEdge edge = ContourTraceComputation::unpackEdge(packedEdge);
            const int columns = tree.numColumns();
            const int numVertexColumns = columns + 1;
            const auto [row, column] = ImageUtils::to2D(edge.pixel, columns);

            switch (edge.side) {
            case ContourTraceSide::North:
                return {vertexId(row, column, numVertexColumns), vertexId(row, column + 1, numVertexColumns), Direction::East, -row};
            case ContourTraceSide::East:
                return {vertexId(row, column + 1, numVertexColumns), vertexId(row + 1, column + 1, numVertexColumns), Direction::South, column + 1};
            case ContourTraceSide::South:
                return {vertexId(row + 1, column + 1, numVertexColumns), vertexId(row + 1, column, numVertexColumns), Direction::West, row + 1};
            case ContourTraceSide::West:
                return {vertexId(row + 1, column, numVertexColumns), vertexId(row, column, numVertexColumns), Direction::North, -column};
            }
            throw std::runtime_error("Invalid contour trace side.");
        }

        /**
         * @brief Returns priority.
         *
         * @param incoming Incoming trace direction.
         * @param outgoing Outgoing trace direction.
         * @return Priority.
         */
        static int turnPriority(Direction incoming, Direction outgoing) {
            const int turn = (static_cast<int>(outgoing) - static_cast<int>(incoming) + 4) & 3;
            static constexpr std::array<int, 4> priorities{
                1, // straight
                0, // right
                3, // back
                2  // left
            };
            return priorities[static_cast<std::size_t>(turn)];
        }

        /**
         * @brief Chooses next edge.
         *
         * @param current Current item in the traversal.
         * @param outgoingHead Head indices of each vertex outgoing-edge list.
         * @param directedEdges Directed edges that form the traced contour graph.
         * @param outgoingNext Next-edge links for the outgoing adjacency lists.
         * @param visitedGeneration Generation marks used to identify visited directed edges.
         * @param visitGeneration Active visit-mark generation.
         * @return Selected next edge.
         */
        static int chooseNextEdge(const DirectedEdge& current, int outgoingHead, const std::vector<DirectedEdge>& directedEdges,
                                  const std::vector<int>& outgoingNext, const std::vector<uint32_t>& visitedGeneration, uint32_t visitGeneration
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
                                  ,
                                  std::size_t& candidateScans, std::size_t& visitedSkips, std::size_t& unvisitedCandidates
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
                const int priority = turnPriority(current.direction, directedEdges[static_cast<std::size_t>(candidate)].direction);
                if (priority < bestPriority || (priority == bestPriority && directedEdges[static_cast<std::size_t>(candidate)].packedEdge <
                                                                                directedEdges[static_cast<std::size_t>(best)].packedEdge)) {
                    best = candidate;
                    bestPriority = priority;
                }
            }
            return best;
        }

        /**
         * @brief Resets trace outgoing heads.
         */
        void resetTraceOutgoingHeads() const {
            for (int vertex : traceTouchedVertices_) {
                const auto vertexIndex = static_cast<std::size_t>(vertex);
                traceOutgoingHead_[vertexIndex] = -1;
            }
            traceTouchedVertices_.clear();
        }

        /**
         * @brief Returns vertex count.
         *
         * @return Vertex count.
         */
        [[nodiscard]] std::size_t imageVertexCount() const {
            return static_cast<std::size_t>((tree.numRows() + 1) * (tree.numColumns() + 1));
        }

        /**
         * @brief Ensures trace outgoing head storage.
         */
        void ensureTraceOutgoingHeadStorage() const {
            if (!traceOutgoingHead_.empty()) {
                return;
            }
            traceOutgoingHead_.assign(imageVertexCount(), -1);
        }

        /**
         * @brief Advances power of two at least.
         *
         * @param value Value.
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
         * @brief Returns outgoing table size.
         *
         * @param edgeCount Number of boundary edges.
         * @return Outgoing table size.
         */
        static std::size_t sparseOutgoingTableSize(std::size_t edgeCount) {
            const std::size_t target = std::max<std::size_t>(2, (2 * edgeCount) + 1);
            return nextPowerOfTwoAtLeast(target);
        }

        /**
         * @brief Returns node local trace limit.
         *
         * @return Node local trace limit.
         */
        static constexpr std::size_t sparseNodeLocalTraceLimit() {
            // Keep sparse adjacency for interactive point queries, then switch
            // to dense storage before random all-node access pays hash overhead.
            return 8;
        }

        /**
         * @brief Tests whether use sparse trace adjacency holds.
         *
         * @param edgeCount Number of boundary edges.
         * @param streamLoopInfosDirectly Whether loop descriptors are streamed directly.
         * @return True when use sparse trace adjacency; otherwise false.
         */
        [[nodiscard]] bool shouldUseSparseTraceAdjacency(std::size_t edgeCount, bool streamLoopInfosDirectly) const {
            if (streamLoopInfosDirectly || !traceOutgoingHead_.empty()) {
                return false;
            }
            if (nodeLocalLoopTraceCount_ >= sparseNodeLocalTraceLimit()) {
                return false;
            }
            const std::size_t sparseBytes = sparseOutgoingTableSize(edgeCount) * ((2 * sizeof(int)) + sizeof(uint32_t));
            const std::size_t denseBytes = imageVertexCount() * sizeof(int);
            return sparseBytes < denseBytes;
        }

        /**
         * @brief Prepares sparse trace outgoing heads.
         *
         * @param edgeCount Number of boundary edges.
         */
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

        /**
         * @brief Computes the sparse-table slot for a trace vertex.
         *
         * @param vertex Image-grid vertex identifier.
         * @return Sparse-table slot selected for the vertex.
         */
        [[nodiscard]] std::size_t sparseTraceSlot(int vertex) const {
            const auto hash = static_cast<uint32_t>(vertex) * uint32_t{2654435761u};
            return static_cast<std::size_t>(hash) & (traceSparseVertexKeys_.size() - 1);
        }

        /**
         * @brief Returns trace outgoing head.
         *
         * @param vertex Image-grid vertex identifier.
         * @return Trace outgoing head.
         */
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

        /**
         * @brief Finds or inserts a vertex in the sparse trace table.
         *
         * @param vertex Image-grid vertex identifier.
         * @return Slot containing the existing or newly inserted vertex.
         */
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

        /**
         * @brief Returns the first outgoing directed edge of a trace vertex.
         *
         * @param vertex Image-grid vertex identifier.
         * @return Directed-edge index at the head of the vertex list, or the empty sentinel.
         */
        template <TraceAdjacencyMode Mode> [[nodiscard]] int traceOutgoingHead(int vertex) const {
            if constexpr (Mode == TraceAdjacencyMode::Sparse) {
                return sparseTraceOutgoingHead(vertex);
            } else {
                return traceOutgoingHead_[static_cast<std::size_t>(vertex)];
            }
        }

        /**
         * @brief Prepares trace adjacency.
         *
         * @param edgeCount Number of boundary edges.
         */
        template <TraceAdjacencyMode Mode> void prepareTraceAdjacency(std::size_t edgeCount) const {
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

        /**
         * @brief Resets trace adjacency.
         */
        template <TraceAdjacencyMode Mode> void resetTraceAdjacency() const {
            if constexpr (Mode == TraceAdjacencyMode::Dense) {
                resetTraceOutgoingHeads();
            }
        }

        /**
         * @brief Pushes trace outgoing edge.
         *
         * @param startVertex Starting vertex identifier.
         * @param edgeIndex Boundary-edge index.
         */
        template <TraceAdjacencyMode Mode> void pushTraceOutgoingEdge(int startVertex, int edgeIndex) const {
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

        /**
         * @brief Advances trace visit generation.
         *
         * @param edgeCount Number of boundary edges.
         */
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

        /**
         * @brief Tests whether trace visited holds.
         *
         * @param edgeIndex Boundary-edge index.
         * @return True when trace visited; otherwise false.
         */
        [[nodiscard]] bool isTraceVisited(int edgeIndex) const { return traceVisitedGeneration_[static_cast<std::size_t>(edgeIndex)] == traceVisitGeneration_; }

        /**
         * @brief Marks trace visited.
         *
         * @param edgeIndex Boundary-edge index.
         */
        void markTraceVisited(int edgeIndex) const { traceVisitedGeneration_[static_cast<std::size_t>(edgeIndex)] = traceVisitGeneration_; }

        /**
         * @brief Reserves loop-metadata storage for a complete-tree trace.
         *
         * @param root Root node of the operation.
         */
        void reserveLoopInfoCapacityForGlobalTrace(NodeId root) const {
            std::size_t pendingNonEmptyNodes = 0;
            std::size_t pendingEdgeCount = 0;
            for (NodeId node : tree.subtreeNodes(root)) {
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

        /**
         * @brief Tests whether reuse edge segment for loops holds.
         *
         * @param node Node identifier.
         * @return True when reuse edge segment for loops; otherwise false.
         */
        bool canReuseEdgeSegmentForLoops(NodeId node) const {
            if (tree.isRoot(node)) {
                return true;
            }
            const NodeId parent = tree.parent(node);
            return parent == InvalidNode || parent == node || cachedEdgeReady_[static_cast<std::size_t>(parent)] != 0;
        }

        /**
         * @brief Ensures node loops materialized.
         *
         * @param node Node identifier.
         */
        void ensureNodeLoopsMaterialized(NodeId node) const {
            requireStableTree("IncrementalContourTraces::ensureNodeLoopsMaterialized");
            requireLiveTraceNode(node, "IncrementalContourTraces::ensureNodeLoopsMaterialized");
            if (cachedLoopReady_[static_cast<std::size_t>(node)]) {
                return;
            }
            ensureEdgesMaterialized(node);
            traceNodeLoops(node, false);
        }

        /**
         * @brief Ensures loops materialized.
         *
         * @param root Root node of the operation.
         */
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
                    for (NodeId child : tree.children(node)) {
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

        /**
         * @brief Traces and caches all contour loops associated with a tree node.
         *
         * @param node Node identifier.
         * @param streamLoopInfosDirectly Whether loop descriptors are streamed directly.
         */
        void traceNodeLoops(NodeId node, bool streamLoopInfosDirectly) const {
            const std::size_t edgeCount = static_cast<std::size_t>(cachedEdgeSize_[static_cast<std::size_t>(node)]);
            const bool useSparseAdjacency = shouldUseSparseTraceAdjacency(edgeCount, streamLoopInfosDirectly);
            if (!streamLoopInfosDirectly) {
                ++nodeLocalLoopTraceCount_;
            }
            if (useSparseAdjacency) {
                traceNodeLoopsWithAdjacency<TraceAdjacencyMode::Sparse>(node, streamLoopInfosDirectly, edgeCount);
            } else {
                traceNodeLoopsWithAdjacency<TraceAdjacencyMode::Dense>(node, streamLoopInfosDirectly, edgeCount);
            }
        }

        /**
         * @brief Traces a node contour using the supplied adjacency representation.
         *
         * @param node Node identifier.
         * @param streamLoopInfosDirectly Whether loop descriptors are streamed directly.
         * @param edgeCount Number of boundary edges.
         */
        template <TraceAdjacencyMode Mode> void traceNodeLoopsWithAdjacency(NodeId node, bool streamLoopInfosDirectly, std::size_t edgeCount) const {
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
                traceDirectedEdges_.emplace_back(*it, geometry.startVertex, geometry.endVertex, geometry.direction, geometry.signedArea2);
                pushTraceOutgoingEdge<Mode>(geometry.startVertex, edgeIndex);
            }
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
            const auto buildAdjacencyEnd = traceProfileNow();
            const auto profileCountersStart = traceProfileNow();
            if constexpr (Mode == TraceAdjacencyMode::Dense) {
                for (int vertex : traceTouchedVertices_) {
                    std::size_t degree = 0;
                    for (int candidate = traceOutgoingHead_[static_cast<std::size_t>(vertex)]; candidate != -1;
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
                    for (int candidate = traceSparseOutgoingHead_[slot]; candidate != -1; candidate = traceOutgoingNext_[static_cast<std::size_t>(candidate)]) {
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
                            current = chooseNextEdge(edge, outgoingHead, traceDirectedEdges_, traceOutgoingNext_, traceVisitedGeneration_, traceVisitGeneration_
#ifdef MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
                                                     ,
                                                     nodeSuccessorCandidateScans, nodeSuccessorVisitedSkips, nodeSuccessorUnvisitedCandidates
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

                    const ContourLoopKind kind = signedArea2 >= 0 ? ContourLoopKind::External : ContourLoopKind::Internal;
                    if (streamLoopInfosDirectly) {
                        cachedLoopInfos_.push_back(ContourTraceLoop{kind, checkedU32(globalEdgeOffset + loopEdgeOffset, "global loop edge offset"),
                                                                    checkedU32(loopEdgeCount, "loop edge count"), signedArea2});
                    } else {
                        traceNodeLoops_.push_back(ContourTraceLoop{kind, checkedU32(loopEdgeOffset, "local loop edge offset"),
                                                                   checkedU32(loopEdgeCount, "loop edge count"), signedArea2});
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
                std::copy(traceNodeLoopEdges_.begin(), traceNodeLoopEdges_.end(), cachedEdgeValues_.begin() + static_cast<std::ptrdiff_t>(globalEdgeOffset));
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
            cachedLoopInfoSize_[static_cast<std::size_t>(node)] = checkedU32(cachedLoopInfos_.size() - loopInfoOffset, "loop info size");
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
                activeTraceProfile_->maxOutgoingDegree = std::max(activeTraceProfile_->maxOutgoingDegree, nodeMaxOutgoingDegree);
                activeTraceProfile_->closedLoopStops += nodeClosedLoopStops;
                activeTraceProfile_->missingOutgoingStops += nodeMissingOutgoingStops;
                activeTraceProfile_->singleSuccessorSteps += nodeSingleSuccessorSteps;
                activeTraceProfile_->singleSuccessorVisitedStops += nodeSingleSuccessorVisitedStops;
                activeTraceProfile_->ambiguousSuccessorSteps += nodeAmbiguousSuccessorSteps;
                activeTraceProfile_->ambiguousSuccessorDeadEnds += nodeAmbiguousSuccessorDeadEnds;
                activeTraceProfile_->successorCandidateScans += nodeSuccessorCandidateScans;
                activeTraceProfile_->successorVisitedSkips += nodeSuccessorVisitedSkips;
                activeTraceProfile_->successorUnvisitedCandidates += nodeSuccessorUnvisitedCandidates;
                activeTraceProfile_->profileCountersNs += traceProfileElapsedNs(profileCountersStart, profileCountersEnd);
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
    /** @brief Defines the `PendingEdgeLists` alias used by the component. */
    using PendingEdgeLists = detail::PendingPixelLists;

    /** @brief Owns the trace deltas extracted from an incremental contour computation. */
    struct ExtractedTraceDeltas {
        /** @brief Deltas. */
        LocalTraceDeltas deltas;
        /** @brief Capacity hint. */
        int capacityHint = 0;
    };

    /**
     * @brief Returns the neighboring pixel reached in the requested direction.
     *
     * @param tree Tree topology.
     * @param pixel Pixel identifier.
     * @param side Contour or boundary side.
     * @return Identifier of the neighboring pixel.
     */
    static PixelId neighborPixel(const MorphologicalTree& tree, PixelId pixel, ContourTraceSide side) {
        const int rows = tree.numRows();
        const int columns = tree.numColumns();
        const auto [row, column] = ImageUtils::to2D(pixel, columns);

        switch (side) {
        case ContourTraceSide::North:
            return row == 0 ? -1 : ImageUtils::to1D(row - 1, column, columns);
        case ContourTraceSide::West:
            return column == 0 ? -1 : ImageUtils::to1D(row, column - 1, columns);
        case ContourTraceSide::East:
            return column == columns - 1 ? -1 : ImageUtils::to1D(row, column + 1, columns);
        case ContourTraceSide::South:
            return row == rows - 1 ? -1 : ImageUtils::to1D(row + 1, column, columns);
        }
        return -1;
    }

    /**
     * @brief Extracts trace deltas impl.
     *
     * @param tree Tree topology.
     * @return Extracted trace deltas impl.
     */
    [[nodiscard]] static ExtractedTraceDeltas extractTraceDeltasImpl(const MorphologicalTree& tree) {
        if (tree.numRows() <= 0 || tree.numColumns() <= 0) {
            throw std::invalid_argument("Contour trace extraction requires a non-empty image domain.");
        }
        if (!tree.isAlive(tree.root())) {
            throw std::invalid_argument("Contour trace extraction requires a live tree root.");
        }

        const int numNodes = tree.numInternalNodeSlots();
        const int totalPixels = tree.numRows() * tree.numColumns();
        const int totalPackedEdges = 4 * totalPixels;
        const int capacityHint = std::max(totalPixels, 1);

        PendingEdgeLists localEdgeAdditions(numNodes, capacityHint);
        PendingEdgeLists localEdgeRemovals(numNodes, capacityHint);

        static constexpr std::array<ContourTraceSide, 4> sides{ContourTraceSide::North, ContourTraceSide::West, ContourTraceSide::East,
                                                               ContourTraceSide::South};

        detail::traversePostOrder(
            tree, tree.root(), [](NodeId) -> void {}, [](NodeId, NodeId) -> void {},
            [&](NodeId nodeId) {
                for (PixelId pixel : tree.properPart(nodeId)) {
                    for (ContourTraceSide side : sides) {
                        const int packedEdge = packEdge(pixel, side);
                        const PixelId neighbor = neighborPixel(tree, pixel, side);
                        if (neighbor < 0) {
                            localEdgeAdditions.add(nodeId, packedEdge);
                            continue;
                        }

                        const NodeId entry = local_attributes::detail::kernel::anchoredEntry(tree, pixel, neighbor);
                        if (entry == InvalidNode || entry == nodeId) {
                            continue;
                        }

                        localEdgeAdditions.add(nodeId, packedEdge);
                        localEdgeRemovals.add(entry, packedEdge);
                    }
                }
            });

        return {LocalTraceDeltas::fromPendingEdgeLists(localEdgeAdditions, localEdgeRemovals, totalPackedEdges), capacityHint};
    }

  public:
    /**
     * @brief Runs incremental side-level contour extraction and returns lazy traces.
     *
     * @param tree Tree topology.
     * @return Result produced by running incremental side-level contour extraction and returns lazy traces.
     */
    [[nodiscard]] static IncrementalContourTraces extract(const MorphologicalTree& tree) {
        ExtractedTraceDeltas extracted = extractTraceDeltasImpl(tree);
        return IncrementalContourTraces(tree, std::move(extracted.deltas), extracted.capacityHint);
    }

    template <AltitudeValue T>
    /**
     * @brief Runs incremental side-level contour extraction on a valued-tree view.
     *
     * @param tree Tree topology.
     * @return Result produced by running incremental side-level contour extraction on a valued-tree view.
     */
    [[nodiscard]] static IncrementalContourTraces extract(const ValuedMorphologicalTreeView<T>& tree) {
        tree.requireTopologyUnchanged("ContourTraceComputation::extract");
        return extract(tree.topology());
    }
};

} // namespace mmcfilters
