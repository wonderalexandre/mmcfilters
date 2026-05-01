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
 * 1. PendingPixelLists
 *    A lightweight transient store keeps per-node pixel lists inside one
 *    contiguous buffer. It avoids repeated heap allocation while contour pixels
 *    are inserted, forwarded to ancestors, and removed.
 *
 * 2. recycle() and consumeInto()
 *    Removed contour pixels do not normally reappear during the same pass, yet
 *    other nodes still need fresh slots. Without recycling, every `add()` call
 *    would keep extending `entries` and memory usage would scale with the
 *    total number of list operations rather than with the simultaneous peak.
 *    `consumeInto()` drains a list, forwards the values to a temporary
 *    container, and returns the slots to the free list through `recycle()`.
 *
 * 3. extractCompactContours()
 *    The incremental traversal first records local contour additions and
 *    removals in transient pixel lists. These lists are compacted into a CSR-like
 *    delta store (`ContourDeltaStore`) before the result is returned. A
 *    temporary list store, `pendingContourRemovals`, stores pixels that must be
 *    removed at some ancestor, typically the lowest common ancestor of two
 *    incomparable nodes.
 *
 * 4. Aggregation phase
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
#include "../attributes/AttributeComputedIncrementally.hpp"
#include "../utils/AdjacencyRelation.hpp"
#include "detail/PendingPixelLists.hpp"
#include "detail/ContourDeltaStore.hpp"

namespace mmcfilters {

/**
 * @brief Arena-based incremental contour extraction and aggregation for `MorphologicalTree`.
 */
class ContoursComputedIncrementally {
public:
    static constexpr double ContourSideAdjacencyRadius = 1.0;

private:
    using PendingPixelLists = detail::PendingPixelLists;
    using ContourDeltaStore = detail::ContourDeltaStore;

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
        /// Immutable compact local additions/removals produced by extraction.
        ContourDeltaStore localDeltas_;
        /// Concatenated storage for all materialized contour slices.
        mutable std::vector<int> cachedContourValues_;
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
         * @param tree ponteiro para a árvore utilizada na computação.
         * @param localDeltas deltas locais compactados por nó.
         * @param capacityHint sugestão para reserva inicial dos contornos agregados.
         */
        IncrementalContours(const MorphologicalTree& tree, ContourDeltaStore localDeltas, int capacityHint)
            : tree(tree),
              localDeltas_(std::move(localDeltas)),
              cachedContourOffset_(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), 0),
              cachedContourSize_(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), 0),
              cachedContourReady_(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), 0),
              pixelMark_(tree.getNumRowsOfImage() * tree.getNumColsOfImage(), 0) {
            if (capacityHint > 0) {
                cachedContourValues_.reserve(static_cast<size_t>(capacityHint));
            }
        }

    public:
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
         * @brief Range cache-aware para iterar o contorno de um nó específico.
         *
         * The range is cheap to copy. Its first `begin()` or `end()` call may
         * materialize the requested subtree; later reads over ready nodes only
         * iterate cached contiguous storage.
         */
        class ContourRange {
        public:
            using iterator = std::vector<int>::const_iterator;

            ContourRange(const IncrementalContours* owner, NodeId node): owner_(owner), node_(node) {}

            iterator begin() const {
                ensureReadable_();
                return owner_->cachedContourBegin(node_);
            }
            iterator end() const {
                ensureReadable_();
                return owner_->cachedContourEnd(node_);
            }

            bool empty() const { return begin() == end(); }

        private:
            void ensureReadable_() const {
                owner_->ensureSubtreeMaterialized(node_);
            }

            const IncrementalContours* owner_ = nullptr;
            NodeId node_ = InvalidNode;
        };

        /**
         * @brief Range lazy que percorre todos os nós, devolvendo pares `(nodeId, contourRange)`.
         *
         * This is only an all-node traversal adapter. It does not use a second
         * contour representation; each returned `ContourRange` still materializes
         * through `getContour`'s cache-aware path.
         */
        class ContoursByNodeRange {
        public:
            /**
             * @brief Iterador forward que devolve pares `(nodeId, ContourRange)`.
             */
            class iterator {
            public:
                using iterator_category = std::forward_iterator_tag;
                using value_type = std::pair<NodeId, ContourRange>;
                using difference_type = std::ptrdiff_t;
                using pointer = void;
                using reference = value_type;

                iterator() = default;
                iterator(const IncrementalContours* owner, NodeId node): owner_(owner), node_(node) {
                    settle_();
                }

                value_type operator*() const { return {node_, ContourRange(owner_, node_)}; }

                iterator& operator++() {
                    ++node_;
                    settle_();
                    return *this;
                }

                iterator operator++(int) {
                    iterator tmp(*this);
                    ++(*this);
                    return tmp;
                }

                friend bool operator==(const iterator& lhs, const iterator& rhs) {
                    return lhs.node_ == rhs.node_;
                }

                friend bool operator!=(const iterator& lhs, const iterator& rhs) {
                    return !(lhs == rhs);
                }

            private:
                void settle_() {
                    if (!owner_) {
                        node_ = InvalidNode;
                        return;
                    }
                    const NodeId numNodeSlots = owner_->tree.getNumInternalNodeSlots();
                    while (node_ >= 0 && node_ < numNodeSlots &&
                           !owner_->tree.isAlive(node_)) {
                        ++node_;
                    }
                    if (node_ >= numNodeSlots) {
                        node_ = numNodeSlots;
                    }
                }

                const IncrementalContours* owner_ = nullptr;
                NodeId node_ = 0;
            };

            explicit ContoursByNodeRange(const IncrementalContours* owner)
                : owner_(owner) {}

            iterator begin() const { return iterator(owner_, 0); }
            iterator end() const { return iterator(owner_, owner_->tree.getNumInternalNodeSlots()); }

        private:
            const IncrementalContours* owner_ = nullptr;
        };

        /**
         * @brief Returns allocation-oriented diagnostics for benchmarks.
         *
         * The byte count is approximate: it reports vector capacities owned by
         * this object and does not include allocator metadata or referenced tree
         * storage.
         */
        StorageStats storageStats() const noexcept {
            StorageStats stats;
            stats.addDeltaValues = localDeltas_.addValues.size();
            stats.removeDeltaValues = localDeltas_.removeValues.size();
            stats.cachedContourValues = cachedContourValues_.size();
            stats.cachedContourCapacity = cachedContourValues_.capacity();
            stats.cachedContourReadyNodes = static_cast<std::size_t>(
                std::count(cachedContourReady_.begin(), cachedContourReady_.end(), uint8_t{1}));
            stats.approxAllocatedBytes =
                localDeltas_.addValues.capacity() * sizeof(int) +
                localDeltas_.removeValues.capacity() * sizeof(int) +
                localDeltas_.addSpans.capacity() * sizeof(ContourDeltaStore::Span) +
                localDeltas_.removeSpans.capacity() * sizeof(ContourDeltaStore::Span) +
                cachedContourValues_.capacity() * sizeof(int) +
                cachedContourOffset_.capacity() * sizeof(uint32_t) +
                cachedContourSize_.capacity() * sizeof(uint32_t) +
                cachedContourReady_.capacity() * sizeof(uint8_t) +
                pixelMark_.capacity() * sizeof(uint16_t);
            return stats;
        }

        /// @return Cache-aware range that iterates the contour of `nodeId`.
        ContourRange getContour(NodeId node) const {
            requireLiveContourNode(node, "IncrementalContours::getContour");
            return ContourRange(this, node);
        }

        /// @return Lazy range over all live nodes as `(nodeId, contourRange)` pairs.
        ContoursByNodeRange contoursByNode() const {
            return ContoursByNodeRange(this);
        }

        /**
         * @brief Materializes and caches every live-node contour.
         *
         * This is a prefetch operation for broad workloads. It uses the same
         * materialization path as `getContour(root)`.
         */
        void materializeAll() const {
            ensureSubtreeMaterialized(tree.getRoot());
        }

        /**
         * @brief Tests whether every live-node contour has already been materialized.
         */
        bool isMaterialized() const noexcept {
            for (NodeId node : tree.getAliveNodeIds()) {
                if (!cachedContourReady_[node]) {
                    return false;
                }
            }
            return true;
        }

        /**
         * @brief Tests whether one live node contour is already materialized.
         */
        bool isContourMaterialized(NodeId node) const {
            requireLiveContourNode(node, "IncrementalContours::isContourMaterialized");
            return static_cast<bool>(cachedContourReady_[node]);
        }

    private:
        /**
         * @brief Rejects invalid or dead nodes before exposing contour ranges.
         */
        void requireLiveContourNode(NodeId node, const char* context) const {
            if (!tree.isAlive(node)) {
                throw std::invalid_argument(std::string(context) + " requires a live internal NodeId.");
            }
        }

        std::vector<int>::const_iterator cachedContourBegin(NodeId node) const {
            return cachedContourValues_.begin() + static_cast<std::ptrdiff_t>(cachedContourOffset_[node]);
        }

        std::vector<int>::const_iterator cachedContourEnd(NodeId node) const {
            return cachedContourBegin(node) + static_cast<std::ptrdiff_t>(cachedContourSize_[node]);
        }

        void nextMarkGeneration() const {
            ++markGeneration_;
            if (markGeneration_ == 0) {
                std::fill(pixelMark_.begin(), pixelMark_.end(), 0);
                markGeneration_ = 1;
            }
        }

        void addIfUnmarked(std::vector<int>& values, int pixel) const {
            if (pixel < 0 || pixel >= static_cast<int>(pixelMark_.size())) {
                return;
            }
            if (pixelMark_[pixel] != markGeneration_) {
                pixelMark_[pixel] = markGeneration_;
                values.push_back(pixel);
            }
        }

        void removeIfMarked(int pixel) const {
            if (pixel >= 0 && pixel < static_cast<int>(pixelMark_.size())) {
                pixelMark_[pixel] = 0;
            }
        }

        void commitMaterializedContour(NodeId node, const std::vector<int>& values) const {
            cachedContourOffset_[node] = checkedU32(cachedContourValues_.size(), "cached contour offset");
            cachedContourSize_[node] = checkedU32(values.size(), "cached contour size");
            cachedContourValues_.insert(cachedContourValues_.end(), values.begin(), values.end());
            cachedContourReady_[node] = 1;
        }

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
         */
        void ensureSubtreeMaterialized(NodeId root) const {
            requireLiveContourNode(root, "IncrementalContours::ensureSubtreeMaterialized");
            if (cachedContourReady_[root]) {
                return;
            }

            std::vector<std::pair<NodeId, bool>> stack;
            stack.emplace_back(root, false);
            std::vector<int> values;

            while (!stack.empty()) {
                const auto [node, expanded] = stack.back();
                stack.pop_back();
                if (cachedContourReady_[node]) {
                    continue;
                }
                if (!expanded) {
                    stack.emplace_back(node, true);
                    for (NodeId child : tree.getChildren(node)) {
                        if (!cachedContourReady_[child]) {
                            stack.emplace_back(child, false);
                        }
                    }
                    continue;
                }

                values.clear();
                const auto additions = localDeltas_.additions(node);
                std::size_t reserveSize = additions.size();
                for (NodeId child : tree.getChildren(node)) {
                    reserveSize += static_cast<std::size_t>(cachedContourSize_[child]);
                }
                values.reserve(reserveSize);
                nextMarkGeneration();

                for (NodeId child : tree.getChildren(node)) {
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

    /**
     * @brief Executa a computação incremental e devolve contornos compactados.
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
     * @param tree árvore morfológica (máx-tree, mín-tree ou ToS) sobre a qual o cálculo será realizado.
     * @return Estrutura `IncrementalContours` contendo deltas compactos para acesso aos contornos.
     *
     * Exemplo:
     * @code
     * auto contours = ContoursComputedIncrementally::extractCompactContours(tree);
     * auto contour = contours.getContour(nodeId);
     * std::vector<int> pixels(contour.begin(), contour.end());
     * @endcode
     */
    static IncrementalContours extractCompactContours(const MorphologicalTree& tree) {
        if (tree.getNumRowsOfImage() <= 0 || tree.getNumColsOfImage() <= 0) {
            throw std::invalid_argument("Contour extraction requires a non-empty image domain.");
        }
        if (!tree.isAlive(tree.getRoot())) {
            throw std::invalid_argument("Contour extraction requires a live tree root.");
        }

        const int numNodes = tree.getNumInternalNodeSlots();
        const int totalPixels = tree.getNumRowsOfImage() * tree.getNumColsOfImage();

        const int capacityHint = std::max(totalPixels / 4, 1);
        // listas temporárias para deltas locais antes da compactação final
        PendingPixelLists localContourPixels(numNodes, capacityHint);
        PendingPixelLists localRemovalPixels(numNodes, capacityHint);
        // lista temporária que carrega pixels até o ancestral correto antes da remoção
        PendingPixelLists pendingContourRemovals(numNodes, capacityHint);

        // contador auxiliar para saber quando um pixel deixa de ser contorno (cf. artigo original)
        std::vector<int> ncount(totalPixels, 0);
        // reuso de armazenamento para pixels que devem ser removidos neste nó
        std::vector<int> removalBuffer;
        removalBuffer.reserve(64);
        AdjacencyRelation adj4(tree.getNumRowsOfImage(), tree.getNumColsOfImage(), ContourSideAdjacencyRadius);
        AttributeComputedIncrementally::traversePostOrder(
            tree,
            tree.getRoot(),
            [](NodeId) -> void {},
            [](NodeId, NodeId) -> void {},
            [&](NodeId nodePId) {
                const NodeId nodeP = nodePId;
                // remove e processa todas as remoções pendentes deste nó
                removalBuffer.clear();
                pendingContourRemovals.consumeInto(nodeP, removalBuffer);
                for (int p : removalBuffer) {
                    bool isPixelToBeRemoved = true;
                    for (int r : adj4.getNeighborPixels(p)) {
                        NodeId nodeRId = tree.getSmallestComponent(r);
                        if (tree.isStrictAncestor(nodeRId, nodePId)) {
                            pendingContourRemovals.add(nodeRId, p);
                            isPixelToBeRemoved = false;
                        } else if (!tree.isComparable(nodePId, nodeRId)) {
                            NodeId otherNodeLCA = tree.getLowestCommonAncestor(nodePId, nodeRId);
                            pendingContourRemovals.add(otherNodeLCA, p);
                            isPixelToBeRemoved = false;
                        }
                    }
                    if (!adj4.isBorderDomainImage(p) && isPixelToBeRemoved) {
                        localRemovalPixels.add(nodeP, p);
                    }
                }

                // percorre os proper parts pertencentes ao nó atual
                for (int p : tree.getProperParts(nodePId)) {
                    if (adj4.isBorderDomainImage(p)) {
                        ncount[p]++;
                    }

                    for (int q : adj4.getNeighborPixels(p)) {
                        NodeId nodeQId = tree.getSmallestComponent(q);
                        if (!tree.isComparable(nodePId, nodeQId)) { // contorno será tratado pelo LCA
                            NodeId nodeLCA = tree.getLowestCommonAncestor(nodePId, nodeQId);
                            pendingContourRemovals.add(nodeLCA, p);
                            ncount[p]++;
                        } else if (tree.isStrictDescendant(nodePId, nodeQId)) { // pixel ainda é fronteira
                            ncount[p]++;
                        } else if (tree.isStrictAncestor(nodePId, nodeQId)) {
                            ncount[q]--;
                            if (ncount[q] == 0) {
                                localRemovalPixels.add(nodeP, q);
                            }
                        }
                    }

                    if (ncount[p] > 0) {
                        localContourPixels.add(nodeP, p);
                    }
                }
            });

        return IncrementalContours(
            tree,
            ContourDeltaStore::fromPendingPixelLists(localContourPixels, localRemovalPixels, totalPixels),
            capacityHint);
    }
};

} // namespace mmcfilters
