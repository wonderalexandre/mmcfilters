#pragma once

/*
 * Overview
 * --------
 * This file implements the arena-based incremental contour computation
 * described for component trees in:
 *
 *   D. J. Da Silva et al., "Incremental component tree contour computation",
 *   Pattern Recognition Letters, 2025.
 *
 * The paper computes contours in the original image domain by counting exposed
 * pixel sides; therefore the contour adjacency is always the 4-neighbourhood,
 * independently from the adjacency used to build the tree. For a tree of
 * shapes, this implementation applies the same 4-connected side-contour
 * definition to the node supports projected onto the original image domain.
 *
 * The design focuses on memory reuse and locality during the post-order passes
 * used to extract and aggregate contour pixels.
 *
 * 1. Boundary-lifetime extraction
 *    Every pixel contributes one contour-birth event at its smallest owner and,
 *    when it becomes interior, one stop event at the LCA of the owners in its
 *    four-neighbourhood. The same semantic event source is consumed by the
 *    incremental contour and max-distance implementations.
 *
 * 2. Compact event storage
 *    Birth/stop events are copied into a CSR-like delta store
 *    (`ContourDeltaStore`) so materialization remains cache-friendly.
 *
 * 3. Aggregation phase
 *    `ensureSubtreeMaterialized()` performs a second post-order pass on
 *    demand. It accumulates local contour pixels, applies deferred removals,
 *    and removes duplicates using generation-marked scratch storage.
 *    Materialized contours are cached per subtree in contiguous vectors, so
 *    subsequent reads of an already materialised node pay only the iteration
 *    cost. Regular `getContour(node)` reads materialize and cache the requested
 *    subtree when needed, so repeated or broad iteration is incremental.
 *
 * See docs/contours.md for the public API, invariants, complexity, memory
 * notes, and benchmark interpretation.
 */

#include "../utils/Common.hpp"
#include "../trees/MorphologicalTree.hpp"
#include "../trees/ValuedMorphologicalTreeView.hpp"
#include "../trees/detail/TreeTraversalDetail.hpp"
#include "../utils/Contract.hpp"
#include "detail/ContourDeltaStore.hpp"
#include "detail/MorphologicalTreeBoundaryLifetimeIndex.hpp"

namespace mmcfilters {

namespace detail::kernel {

/** @brief Radius selecting four-neighbour side adjacency for contour kernels. */
inline constexpr double ContourSideAdjacencyRadius = 1.0;

/** @brief Compact contour deltas together with their expected materialization capacity. */
struct ExtractedContourDeltas {
    ContourDeltaStore deltas; ///< Per-node local contour additions and removals.
    int capacityHint = 0;     ///< Initial capacity suggested for materialized contours.
};

/**
 * @brief Extracts local contour deltas from an established tree and grid domain.
 * @param tree Established tree topology and grid domain.
 * @return Extracted deltas and capacity hint.
 */
inline ExtractedContourDeltas extractContourDeltas(const MorphologicalTree& tree) {
    const int totalPixels = tree.numPixels();
    const int capacityHint = std::max(totalPixels / 4, 1);
    const contours::detail::MorphologicalTreeBoundaryLifetimeIndex events(tree);
    return {ContourDeltaStore::fromEventSource(tree, events), capacityHint};
}

/**
 * @brief Extracts only the compact local contour-delta store.
 * @param tree Established tree topology and grid domain.
 * @return Per-node local contour additions and removals.
 */
inline ContourDeltaStore extractLocalContourDeltas(const MorphologicalTree& tree) {
    ExtractedContourDeltas extracted = extractContourDeltas(tree);
    return std::move(extracted.deltas);
}

} // namespace detail::kernel

/**
 * @brief Arena-based incremental contour extraction and aggregation for `MorphologicalTree`.
 */
class ContoursComputedIncrementally {
  public:
    /// Radius of the 4-neighbour side-contour adjacency used by the algorithm.
    static constexpr double ContourSideAdjacencyRadius = detail::kernel::ContourSideAdjacencyRadius;

    /// Compact local contour additions/removals indexed by internal node id.
    using LocalContourDeltas = detail::ContourDeltaStore;

  private:
    /** @brief Defines the `ContourDeltaStore` alias used by the component. */
    using ContourDeltaStore = LocalContourDeltas;

  public:
    /**
     * @brief Incremental contour result stored as compact local deltas.
     *
     * It keeps local contour additions/removals as compact spans and the materialized
     * representation exposed to callers through range-based iteration helpers.
     */
    struct IncrementalContours {
      private:
        friend class ContoursComputedIncrementally;

        /// Tree used to interpret node ids and child/parent relations. Not owned.
        const MorphologicalTree& tree;
        /// Tree mutation version captured when local contour deltas were computed.
        std::size_t treeMutationVersion_ = 0;
        /// Immutable compact local additions/removals produced by extraction.
        ContourDeltaStore localDeltas_;
        /// Concatenated storage for all materialized contour slices.
        mutable std::vector<PixelId> cachedContourValues_;
        /// Per-node offset into `cachedContourValues_`; valid only when ready.
        mutable std::vector<uint32_t> cachedContourOffset_;
        /// Per-node slice length in `cachedContourValues_`; valid only when ready.
        mutable std::vector<uint32_t> cachedContourSize_;
        /// Per-node materialization flag. A ready node has valid offset/size.
        mutable std::vector<uint8_t> cachedContourReady_;
        /// Scratch generation marks used while building one materialized contour.
        mutable std::vector<uint16_t> pixelMark_;
        /// Current non-zero generation in `pixelMark_`.
        mutable uint16_t markGeneration_ = 1;

        /**
         * @brief Builds a lazy incremental-contour result from compact local deltas.
         *
         * @param tree Tree used by the contour computation.
         * @param localDeltas Compacted local deltas per node.
         * @param capacityHint Initial reservation hint for aggregate contours.
         *
         */
        IncrementalContours(const MorphologicalTree& tree, ContourDeltaStore localDeltas, int capacityHint)
            : tree(tree), treeMutationVersion_(tree.getMutationVersion()), localDeltas_(std::move(localDeltas)),
              cachedContourOffset_(static_cast<std::size_t>(tree.numInternalNodeSlots()), 0),
              cachedContourSize_(static_cast<std::size_t>(tree.numInternalNodeSlots()), 0),
              cachedContourReady_(static_cast<std::size_t>(tree.numInternalNodeSlots()), 0),
              pixelMark_(tree.numRows() * tree.numColumns(), 0) {
            if (capacityHint > 0) {
                cachedContourValues_.reserve(static_cast<size_t>(capacityHint));
            }
        }

      public:
        /**
         * @brief Allocation and cache diagnostics for an incremental contour result.
         */
        struct StorageStats {
            /// Number of compacted local contour-addition pixels.
            std::size_t addDeltaValues = 0;
            /// Number of compacted local contour-removal pixels.
            std::size_t removeDeltaValues = 0;
            /// Number of pixels already committed to materialized contour cache.
            std::size_t cachedContourValues = 0;
            /// Reserved capacity of the materialized contour cache.
            std::size_t cachedContourCapacity = 0;
            /// Number of live-node slots with ready materialized contours.
            std::size_t cachedContourReadyNodes = 0;
            /// Approximate bytes reserved by the owned vectors in this result.
            std::size_t approxAllocatedBytes = 0;
        };

        /**
         * @brief Cache-aware range over the contour of one node.
         *
         * The range is cheap to copy. Its first `begin()` or `end()` call may
         * materialize the requested subtree; later reads over ready nodes only
         * iterate cached contiguous storage.
         */
        class ContourRange {
          public:
            /// Immutable iterator over materialized contour pixel indices.
            using iterator = std::vector<PixelId>::const_iterator;

            /**
             * @brief Creates a range for `node` backed by the contour source.
             *
             * @param source Incremental contour computation serving the range.
             * @param node Node identifier.
             */
            ContourRange(const IncrementalContours* source, NodeId node) : source_(source), node_(node) {}

            /**
             * @brief Returns the first contour pixel, materializing on demand.
             *
             * @return The first contour pixel, materializing on demand.
             */
            iterator begin() const {
                ensureReadable_();
                return source_->cachedContourBegin(node_);
            }

            /**
             * @brief Returns the end sentinel for the contour range.
             *
             * @return The end sentinel for the contour range.
             */
            iterator end() const {
                ensureReadable_();
                return source_->cachedContourEnd(node_);
            }

            /**
             * @brief Returns true when the node contour is empty.
             *
             * @return True when the node contour is empty.
             */
            bool empty() const { return begin() == end(); }

          private:
            /**
             * @brief Ensures readable.
             */
            void ensureReadable_() const { source_->ensureSubtreeMaterialized(node_); }

            /** @brief References the incremental contour computation serving this range. */
            const IncrementalContours* source_ = nullptr;
            /** @brief Dense node identifier held by this record. */
            NodeId node_ = InvalidNode;
        };

        /**
         * @brief Lazy range over all live nodes as `(nodeId, contourRange)` pairs.
         *
         * This is only an all-node traversal adapter. It does not use a second
         * contour representation; each returned `ContourRange` still materializes
         * through `getContour`'s cache-aware path.
         */
        class ContoursByNodeRange {
          public:
            /**
             * @brief Forward iterator that yields `(nodeId, ContourRange)` pairs.
             */
            class iterator {
              public:
                /// Standard iterator category exposed for STL compatibility.
                using iterator_category = std::forward_iterator_tag;

                /// Pair of live node id and lazy contour range.
                using value_type = std::pair<NodeId, ContourRange>;

                /// Signed distance type exposed for STL compatibility.
                using difference_type = std::ptrdiff_t;

                /// This iterator returns values by value, so no pointer is exposed.
                using pointer = void;

                /// Value type returned by dereference.
                using reference = value_type;

                /**
                 * @brief Constructs a default `iterator`.
                 */
                iterator() = default;

                /**
                 * @brief Creates an iterator over `source` starting at `node`.
                 *
                 * @param source Incremental contour computation serving the iterator.
                 * @param node Node identifier.
                 */
                iterator(const IncrementalContours* source, NodeId node) : source_(source), node_(node) { settle_(); }

                /**
                 * @brief Returns the current `(nodeId, contourRange)` pair.
                 *
                 * @return The current (nodeId, contourRange) pair.
                 */
                value_type operator*() const { return {node_, ContourRange(source_, node_)}; }

                /**
                 * @brief Advances to the next live node.
                 *
                 * @return Mutable reference to the updated object.
                 */
                iterator& operator++() {
                    ++node_;
                    settle_();
                    return *this;
                }

                /**
                 * @brief Advances to the next live node and returns the previous iterator.
                 *
                 * @return Iterator position before the advancement.
                 */
                iterator operator++(int) {
                    iterator tmp(*this);
                    ++(*this);
                    return tmp;
                }

                /**
                 * @brief Returns true when both iterators refer to the same node position.
                 */
                friend bool operator==(const iterator& lhs, const iterator& rhs) { return lhs.node_ == rhs.node_; }

                /**
                 * @brief Returns true when iterators refer to different node positions.
                 */
                friend bool operator!=(const iterator& lhs, const iterator& rhs) { return !(lhs == rhs); }

              private:
                /**
                 * @brief Advances to the next node with a materialized contour.
                 */
                void settle_() {
                    if (!source_) {
                        node_ = InvalidNode;
                        return;
                    }
                    const int numNodeSlots = source_->tree.numInternalNodeSlots();
                    while (node_ >= 0 && node_ < numNodeSlots && !source_->tree.isAlive(node_)) {
                        ++node_;
                    }
                    if (node_ >= numNodeSlots) {
                        node_ = numNodeSlots;
                    }
                }

                /** @brief References the incremental contour computation serving this iterator. */
                const IncrementalContours* source_ = nullptr;
                /** @brief Dense node identifier held by this record. */
                NodeId node_ = 0;
            };

            /**
             * @brief Creates a lazy all-node contour range backed by the contour source.
             *
             * @param source Incremental contour computation serving the range.
             */
            explicit ContoursByNodeRange(const IncrementalContours* source) : source_(source) {}

            /**
             * @brief Returns an iterator positioned at the first live node.
             *
             * @return An iterator positioned at the first live node.
             */
            iterator begin() const {
                source_->requireStableTree("IncrementalContours::contoursByNode");
                return iterator(source_, 0);
            }

            /**
             * @brief Returns the all-node range sentinel.
             *
             * @return The all-node range sentinel.
             */
            iterator end() const {
                source_->requireStableTree("IncrementalContours::contoursByNode");
                return iterator(source_, source_->tree.numInternalNodeSlots());
            }

          private:
            /** @brief References the incremental contour computation serving this range. */
            const IncrementalContours* source_ = nullptr;
        };

        /**
         * @brief Returns allocation-oriented diagnostics for benchmarks.
         *
         * The byte count is approximate: it reports vector capacities owned by
         * this object and does not include allocator metadata or referenced tree
         * storage.
         *
         * @return Allocation-oriented diagnostics for benchmarks.
         */
        [[nodiscard]] StorageStats storageStats() const {
            requireStableTree("IncrementalContours::storageStats");
            StorageStats stats;
            stats.addDeltaValues = localDeltas_.addValues.size();
            stats.removeDeltaValues = localDeltas_.removeValues.size();
            stats.cachedContourValues = cachedContourValues_.size();
            stats.cachedContourCapacity = cachedContourValues_.capacity();
            stats.cachedContourReadyNodes = static_cast<std::size_t>(std::count(cachedContourReady_.begin(), cachedContourReady_.end(), uint8_t{1}));
            stats.approxAllocatedBytes = localDeltas_.addValues.capacity() * sizeof(int) + localDeltas_.removeValues.capacity() * sizeof(int) +
                                         localDeltas_.addSpans.capacity() * sizeof(ContourDeltaStore::Span) +
                                         localDeltas_.removeSpans.capacity() * sizeof(ContourDeltaStore::Span) + cachedContourValues_.capacity() * sizeof(int) +
                                         cachedContourOffset_.capacity() * sizeof(uint32_t) + cachedContourSize_.capacity() * sizeof(uint32_t) +
                                         cachedContourReady_.capacity() * sizeof(uint8_t) + pixelMark_.capacity() * sizeof(uint16_t);
            return stats;
        }

        /**
         * @brief Returns cache-aware range that iterates the contour of nodeId.
         *
         * @param node Node identifier.
         * @return Cache-aware range that iterates the contour of `nodeId`.
         *
         */
        [[nodiscard]] ContourRange getContour(NodeId node) const {
            requireStableTree("IncrementalContours::getContour");
            requireLiveContourNode(node, "IncrementalContours::getContour");
            return ContourRange(this, node);
        }

        /**
         * @brief Returns lazy range over all live nodes as (nodeId, contourRange) pairs.
         *
         * @return Lazy range over all live nodes as `(nodeId, contourRange)` pairs.
         */
        [[nodiscard]] ContoursByNodeRange contoursByNode() const {
            requireStableTree("IncrementalContours::contoursByNode");
            return ContoursByNodeRange(this);
        }

        /**
         * @brief Materializes and caches every live-node contour.
         *
         * This is a prefetch operation for broad workloads. It uses the same
         * materialization path as `getContour(root)`.
         */
        void materializeAll() const {
            requireStableTree("IncrementalContours::materializeAll");
            ensureSubtreeMaterialized(tree.root());
        }

        /**
         * @brief Tests whether every live-node contour has already been materialized.
         *
         * @return True if every live-node contour has already been materialized; otherwise false.
         */
        [[nodiscard]] bool isMaterialized() const {
            requireStableTree("IncrementalContours::isMaterialized");
            for (NodeId node : tree.aliveNodeIds()) {
                if (!cachedContourReady_[node]) {
                    return false;
                }
            }
            return true;
        }

        /**
         * @brief Tests whether one live node contour is already materialized.
         *
         * @param node Node identifier.
         * @return True if one live node contour is already materialized; otherwise false.
         */
        [[nodiscard]] bool isContourMaterialized(NodeId node) const {
            requireStableTree("IncrementalContours::isContourMaterialized");
            requireLiveContourNode(node, "IncrementalContours::isContourMaterialized");
            return static_cast<bool>(cachedContourReady_[node]);
        }

      private:
        /**
         * @brief Validates stable tree.
         *
         * @param context Operation name used in diagnostics.
         */
        void requireStableTree(const char* context) const { tree.requireMutationVersion(treeMutationVersion_, context); }

        /**
         * @brief Rejects invalid or dead nodes before exposing contour ranges.
         *
         * @param node Node identifier.
         * @param context Operation context or diagnostic label.
         */
        void requireLiveContourNode(NodeId node, const char* context) const {
            if (!tree.isAlive(node)) {
                throw std::invalid_argument(std::string(context) + " requires a live internal NodeId.");
            }
        }

        /**
         * @brief Returns an iterator to the first cached contour element of a node.
         *
         * @param node Node identifier.
         * @return Values produced by the operation.
         */
        std::vector<PixelId>::const_iterator cachedContourBegin(NodeId node) const {
            return cachedContourValues_.begin() + static_cast<std::ptrdiff_t>(cachedContourOffset_[node]);
        }

        /**
         * @brief Returns an iterator past the last cached contour element of a node.
         *
         * @param node Node identifier.
         * @return Values produced by the operation.
         */
        std::vector<PixelId>::const_iterator cachedContourEnd(NodeId node) const {
            return cachedContourBegin(node) + static_cast<std::ptrdiff_t>(cachedContourSize_[node]);
        }

        /**
         * @brief Advances mark generation.
         */
        void nextMarkGeneration() const {
            ++markGeneration_;
            if (markGeneration_ == 0) {
                std::fill(pixelMark_.begin(), pixelMark_.end(), 0);
                markGeneration_ = 1;
            }
        }

        /**
         * @brief Adds if unmarked.
         *
         * @param values Values read or written by the operation.
         * @param pixel Pixel identifier.
         */
        void addIfUnmarked(std::vector<PixelId>& values, PixelId pixel) const {
            if (pixel < 0 || pixel >= static_cast<int>(pixelMark_.size())) {
                return;
            }
            if (pixelMark_[pixel] != markGeneration_) {
                pixelMark_[pixel] = markGeneration_;
                values.push_back(pixel);
            }
        }

        /**
         * @brief Removes if marked.
         *
         * @param pixel Pixel identifier.
         */
        void removeIfMarked(PixelId pixel) const {
            if (pixel >= 0 && pixel < static_cast<int>(pixelMark_.size())) {
                pixelMark_[pixel] = 0;
            }
        }

        /**
         * @brief Stores a materialized node contour in the shared cache.
         *
         * @param node Node identifier.
         * @param values Values read or written by the operation.
         */
        void commitMaterializedContour(NodeId node, const std::vector<PixelId>& values) const {
            cachedContourOffset_[node] = checkedU32(cachedContourValues_.size(), "cached contour offset");
            cachedContourSize_[node] = checkedU32(values.size(), "cached contour size");
            cachedContourValues_.insert(cachedContourValues_.end(), values.begin(), values.end());
            cachedContourReady_[node] = 1;
        }

        /**
         * @brief Checks and converts u32.
         *
         * @param value Value.
         * @param context Operation name used in diagnostics.
         * @return Value converted to `uint32_t` after range validation.
         */
        static uint32_t checkedU32(std::size_t value, const char* context) {
            if (value > static_cast<std::size_t>(std::numeric_limits<uint32_t>::max())) {
                throw std::overflow_error(std::string(context) + " exceeds uint32_t.");
            }
            return static_cast<uint32_t>(value);
        }

        /**
         * @brief Materializes missing contour caches in the subtree rooted at `root`.
         *
         * The traversal is post-order. Each newly materialized node reuses
         * already materialized child contours, adds its local contour pixels, applies
         * deferred removals, and commits a contiguous cached slice.
         *
         * Already-ready children are not expanded again. This is the mechanism
         * that makes repeated and broad iteration incremental.
         *
         * @param root Root node of the traversal.
         */
        void ensureSubtreeMaterialized(NodeId root) const {
            requireStableTree("IncrementalContours::ensureSubtreeMaterialized");
            requireLiveContourNode(root, "IncrementalContours::ensureSubtreeMaterialized");
            if (cachedContourReady_[root]) {
                return;
            }

            std::vector<std::pair<NodeId, bool>> stack;
            stack.emplace_back(root, false);
            std::vector<PixelId> values;

            while (!stack.empty()) {
                const auto [node, expanded] = stack.back();
                stack.pop_back();
                if (cachedContourReady_[node]) {
                    continue;
                }
                if (!expanded) {
                    stack.emplace_back(node, true);
                    for (NodeId child : tree.children(node)) {
                        if (!cachedContourReady_[child]) {
                            stack.emplace_back(child, false);
                        }
                    }
                    continue;
                }

                values.clear();
                const auto additions = localDeltas_.additions(node);
                std::size_t reserveSize = additions.size();
                for (NodeId child : tree.children(node)) {
                    reserveSize += static_cast<std::size_t>(cachedContourSize_[child]);
                }
                values.reserve(reserveSize);
                nextMarkGeneration();

                for (NodeId child : tree.children(node)) {
                    for (auto it = cachedContourBegin(child); it != cachedContourEnd(child); ++it) {
                        addIfUnmarked(values, *it);
                    }
                }

                for (int value : additions) {
                    addIfUnmarked(values, value);
                }

                for (int rem : localDeltas_.removals(node)) {
                    removeIfMarked(rem);
                }

                std::size_t writeIndex = 0;
                for (int value : values) {
                    if (pixelMark_[static_cast<std::size_t>(value)] == markGeneration_) {
                        values[writeIndex++] = value;
                    }
                }
                values.resize(writeIndex);
                commitMaterializedContour(node, values);
            }
        }
    };

  public:
    /**
     * @brief Extracts compact local contour additions/removals without materializing aggregate contours.
     *
     * This exposes the CSR-like delta store used internally by
     * `IncrementalContours`. It is intended for algorithms that consume local
     * contour events directly, such as incremental distance-transform updates.
     *
     * @param tree Tree topology.
     * @return The extracted compact local contour additions/removals without materializing aggregate contours.
     */
    [[nodiscard]] static LocalContourDeltas extractContourDeltas(const MorphologicalTree& tree) {
        MMCFILTERS_CONTRACT_REQUIRE(tree.numRows() > 0 && tree.numColumns() > 0,
                                    throw std::invalid_argument("Contour extraction requires a non-empty image domain."));
        MMCFILTERS_CONTRACT_REQUIRE(tree.isAlive(tree.root()), throw std::invalid_argument("Contour extraction requires a live tree root."));
        return detail::kernel::extractLocalContourDeltas(tree);
    }

  public:
    template <AltitudeValue T>
    /**
     * @brief Extracts compact local contour additions/removals from a valued-tree view.
     *
     * @param tree Tree topology.
     * @return The extracted compact local contour additions/removals from a valued-tree view.
     */
    [[nodiscard]] static LocalContourDeltas extractContourDeltas(const ValuedMorphologicalTreeView<T>& tree) {
        MMCFILTERS_CONTRACT_CHECKED_ONLY(tree.requireTopologyUnchanged("ContoursComputedIncrementally::extractContourDeltas"));
        MMCFILTERS_CONTRACT_REQUIRE(tree.topology().numRows() > 0 && tree.topology().numColumns() > 0,
                                    throw std::invalid_argument("Contour extraction requires a non-empty image domain."));
        MMCFILTERS_CONTRACT_REQUIRE(tree.topology().isAlive(tree.topology().root()),
                                    throw std::invalid_argument("Contour extraction requires a live tree root."));
        return detail::kernel::extractLocalContourDeltas(tree.topology());
    }

    /**
     * @brief Runs incremental contour computation and returns compact contours.
     *
     * @details
     * This is the component-tree contour algorithm from Da Silva et al. (PRL
     * 2025). For max-trees and min-trees it follows the paper directly. For
     * tree-of-shapes inputs, it computes 4-connected side contours of each
     * projected node support in the original image domain; this ToS use is an
     * implementation extension rather than a claim from the paper.
     *
     * The contour neighbourhood is intentionally fixed to 4-connectivity
     * because the contour is defined through exposed pixel sides, not through
     * the adjacency relation used to construct the input tree.
     *
     * @param tree Morphological tree (max-tree, min-tree, or ToS) on which the
     * contour computation is performed.
     * @return `IncrementalContours` containing compact deltas for contour access.
     *
     * Example:
     * @code
     * auto contours = ContoursComputedIncrementally::extractCompactContours(tree);
     * auto contour = contours.getContour(nodeId);
     * std::vector<PixelId> pixels(contour.begin(), contour.end());
     * @endcode
     */
    [[nodiscard]] static IncrementalContours extractCompactContours(const MorphologicalTree& tree) {
        MMCFILTERS_CONTRACT_REQUIRE(tree.numRows() > 0 && tree.numColumns() > 0,
                                    throw std::invalid_argument("Contour extraction requires a non-empty image domain."));
        MMCFILTERS_CONTRACT_REQUIRE(tree.isAlive(tree.root()), throw std::invalid_argument("Contour extraction requires a live tree root."));
        detail::kernel::ExtractedContourDeltas extracted = detail::kernel::extractContourDeltas(tree);
        return IncrementalContours(tree, std::move(extracted.deltas), extracted.capacityHint);
    }

    template <AltitudeValue T>
    /**
     * @brief Runs incremental contour computation on a valued-tree view.
     *
     * @param tree Tree topology.
     * @return Result produced by running incremental contour computation on a valued-tree view.
     */
    [[nodiscard]] static IncrementalContours extractCompactContours(const ValuedMorphologicalTreeView<T>& tree) {
        MMCFILTERS_CONTRACT_CHECKED_ONLY(tree.requireTopologyUnchanged("ContoursComputedIncrementally::extractCompactContours"));
        MMCFILTERS_CONTRACT_REQUIRE(tree.topology().numRows() > 0 && tree.topology().numColumns() > 0,
                                    throw std::invalid_argument("Contour extraction requires a non-empty image domain."));
        MMCFILTERS_CONTRACT_REQUIRE(tree.topology().isAlive(tree.topology().root()),
                                    throw std::invalid_argument("Contour extraction requires a live tree root."));
        detail::kernel::ExtractedContourDeltas extracted = detail::kernel::extractContourDeltas(tree.topology());
        return IncrementalContours(tree.topology(), std::move(extracted.deltas), extracted.capacityHint);
    }
};

} // namespace mmcfilters
