#pragma once

/*
 * Overview
 * --------
 * This file implements an arena-based incremental contour computation for
 * morphological trees. The design focuses on memory reuse and locality during
 * the post-order passes used to extract and aggregate contour pixels.
 *
 * 1. ListArena
 *    A lightweight arena stores many singly linked lists inside one contiguous
 *    buffer (`entries`). Each node owns one list through `head`, which avoids
 *    repeated heap allocation while contour pixels are inserted and removed.
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
 *    The incremental traversal builds an `IncrementalContours` object
 *    (`result`) that stores the exposed contour and removal arenas. A
 *    temporary arena, `contoursToRemoveLCA`, stores pixels that must be
 *    removed at some ancestor, typically the lowest common ancestor of two
 *    incomparable nodes.
 *
 * 4. Aggregation phase
 *    `ensureAggregated()` performs a second post-order pass on demand. It
 *    accumulates local contour pixels, applies deferred removals, and removes
 *    duplicates using a compact bitmap. The final result is stored in the
 *    `aggregated` arena, so subsequent reads pay only the already materialised
 *    access cost.
 */

#include "../utils/Common.hpp"
#include "../trees/MorphologicalTree.hpp"
#include "../attributes/AttributeComputedIncrementally.hpp"
#include "../utils/AdjacencyRelation.hpp"

namespace mmcfilters {

/**
 * @brief Arena-based incremental contour extraction and aggregation for `MorphologicalTree`.
 */
class ContoursComputedIncrementally {
public:
    static inline NodeId nodeIdOf(const MorphologicalTree&, NodeId nodeId) noexcept {
        return nodeId;
    }

    /**
     * @brief Arena storing many integer lists inside a single contiguous buffer.
     */
    struct ListArena {
        /**
         * @brief Single linked-list entry stored in the arena.
         */
        struct Entry {
            int value;
            int next;
        };

        /**
         * @brief Forward iterator over one arena-managed list.
         */
        class const_iterator {
        public:
            using iterator_category = std::forward_iterator_tag;
            using value_type = int;
            using difference_type = std::ptrdiff_t;
            using pointer = const int*;
            using reference = const int&;

            const_iterator() = default;
            const_iterator(const ListArena* arena, int index): arena_(arena), index_(index) {}

            int operator*() const { return arena_->entries[index_].value; }

            const_iterator& operator++() {
                index_ = (index_ == -1) ? -1 : arena_->entries[index_].next;
                return *this;
            }

            const_iterator operator++(int) {
                const_iterator tmp(*this);
                ++(*this);
                return tmp;
            }

            friend bool operator==(const const_iterator& lhs, const const_iterator& rhs) {
                return lhs.index_ == rhs.index_;
            }

            friend bool operator!=(const const_iterator& lhs, const const_iterator& rhs) {
                return !(lhs == rhs);
            }

        private:
            const ListArena* arena_ = nullptr;
            int index_ = -1;
        };

        /**
         * @brief Lightweight range wrapper for one arena-managed list.
         */
        class Range {
        public:
            Range(const ListArena* arena, int head) : arena_(arena), head_(head) {}

            const_iterator begin() const { return const_iterator(arena_, head_); }
            const_iterator end() const { return const_iterator(arena_, -1); }
            bool empty() const { return head_ == -1; }

        private:
            const ListArena* arena_;
            int head_;
        }; 

        /// Builds an empty arena without reserving storage.
        ListArena() = default;

        std::vector<Entry> entries;
        std::vector<int> head;
        std::vector<int> size;
        int freeHead = -1;

        /// Builds an arena for `numNodes` lists with an optional capacity hint.
        explicit ListArena(int numNodes, int capacityHint = 0): head(numNodes, -1), size(numNodes, 0){
            if (capacityHint > 0) {
                entries.reserve(capacityHint);
            }
        }

        /// Removes every element associated with `node` and returns the slots to the free list.
        void clearNode(int node) {
            int idx = head[node];
            while (idx != -1) {
                const int next = entries[idx].next;
                recycle(idx);
                idx = next;
            }
            head[node] = -1;
            size[node] = 0;
        }

        /// Inserts `value` into the list owned by `node`.
        void add(int node, int value) {
            const int slot = allocate();
            entries[slot] = Entry{value, head[node]};
            head[node] = slot;
            ++size[node];
        }

        Range range(int node) const { return Range(this, head[node]); }

        /**
         * @brief Moves every element of `node` into `out` while recycling the used slots.
         *
         * @tparam OutputContainer container que expõe `push_back(int)`.
         * @param node índice da lista a ser drenada.
         * @param out  recipiente para onde os valores serão copiados.
         */
        template <typename OutputContainer>
        void consumeInto(int node, OutputContainer& out) {
            int idx = head[node];
            while (idx != -1) {
                const int next = entries[idx].next;
                out.push_back(entries[idx].value);
                recycle(idx);
                idx = next;
            }
            head[node] = -1;
            size[node] = 0;
        }

        /// @return Number of elements currently associated with `node`.
        int sizeOf(int node) const { return size[node]; }

        const_iterator begin(int node) const { return const_iterator(this, head[node]); }
        const_iterator end() const { return const_iterator(this, -1); }

    private:
        /// Returns a free slot, either from the free list or via `push_back`.
        int allocate() {
            if (freeHead == -1) {
                entries.push_back(Entry{0, -1});
                return static_cast<int>(entries.size() - 1);
            }
            const int idx = freeHead;
            freeHead = entries[idx].next;
            return idx;
        }

        /// Returns an index to the free list.
        void recycle(int idx) {
            entries[idx].next = freeHead;
            freeHead = idx;
        }
    };

    /**
     * @brief Incremental contour result stored in arenas.
     *
     * It keeps raw contour lists, deferred removals, and the lazily aggregated
     * representation exposed to callers through range-based iteration helpers.
     */
    struct IncrementalContours {
        const MorphologicalTree& tree;
        ListArena contours;
        ListArena removals;
        mutable ListArena aggregated;
        mutable bool aggregatedReady_ = false;

        /**
         * @param tree ponteiro para a árvore utilizada na computação.
         * @param numNodes número de nós (define o tamanho das arenas).
         * @param capacityHint sugestão para reserva inicial de entradas.
         */
        IncrementalContours(const MorphologicalTree& tree, int numNodes, int capacityHint)
            : tree(tree), contours(numNodes, capacityHint), removals(numNodes, capacityHint) {}

        /**
         * @brief Proxy para iterar o contorno agregado de um nó específico
         *        usando a arena consolidada.
         */
        class ContourProxy {
        public:
            /**
             * @brief Iterador forward sobre os pixels do contorno agregado.
             */
            class iterator {
            public:
                using iterator_category = std::forward_iterator_tag;
                using value_type = int;
                using difference_type = std::ptrdiff_t;
                using pointer = const int*;
                using reference = const int&;

                iterator() = default;
                iterator(const ContourProxy*, ListArena::const_iterator current): current_(current) {}

                int operator*() const { return *current_; }

                iterator& operator++() {
                    ++current_;
                    return *this;
                }

                iterator operator++(int) {
                    iterator tmp(*this);
                    ++(*this);
                    return tmp;
                }

                friend bool operator==(const iterator& lhs, const iterator& rhs) {
                    return lhs.current_ == rhs.current_;
                }

                friend bool operator!=(const iterator& lhs, const iterator& rhs) {
                    return !(lhs == rhs);
                }

            private:
                ListArena::const_iterator current_;
            };

            ContourProxy(const IncrementalContours* owner, NodeId node): owner_(owner), node_(node) {}

            iterator begin() const {
                owner_->ensureAggregated();
                return iterator(this, owner_->aggregated.begin(node_));
            }
            iterator end() const {
                owner_->ensureAggregated();
                return iterator(this, owner_->aggregated.end());
            }

            bool empty() const { return begin() == end(); }

        private:
            const IncrementalContours* owner_ = nullptr;
            NodeId node_ = InvalidNode;
        };

        /**
         * @brief Range lazy que percorre todos os nós, devolvendo pares `(nodeId, proxy)`.
         */
        class ContoursLazyRange {
        public:
            /**
             * @brief Iterador forward que devolve pares `(nodeId, ContourProxy)`.
             */
            class iterator {
            public:
                using iterator_category = std::forward_iterator_tag;
                using value_type = std::pair<NodeId, ContourProxy>;
                using difference_type = std::ptrdiff_t;
                using pointer = void;
                using reference = value_type;

                iterator() = default;
                iterator(const IncrementalContours* owner, NodeId node): owner_(owner), node_(node) {
                    settle_();
                }

                value_type operator*() const { return {node_, ContourProxy(owner_, node_)}; }

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
                           !owner_->tree.isAlive(nodeIdOf(owner_->tree, node_))) {
                        ++node_;
                    }
                    if (node_ >= numNodeSlots) {
                        node_ = numNodeSlots;
                    }
                }

                const IncrementalContours* owner_ = nullptr;
                NodeId node_ = 0;
            };

            explicit ContoursLazyRange(const IncrementalContours* owner)
                : owner_(owner) {}

            iterator begin() const { return iterator(owner_, 0); }
            iterator end() const { return iterator(owner_, owner_->tree.getNumInternalNodeSlots()); }

        private:
            const IncrementalContours* owner_ = nullptr;
        };

        /// @return Proxy that iterates the aggregated contour of `nodeId`.
        ContourProxy contour(NodeId node) const { return ContourProxy(this, node); }

        /// @return Lazy range over all live nodes as `(nodeId, contourProxy)` pairs.
        ContoursLazyRange contoursLazy() const { return ContoursLazyRange(this); }

        /**
         * @brief Executa `visitor(value)` para cada pixel do contorno agregado de `node`.
         */
        template <typename Visitor>
        void forEachContourPixel(NodeId node, Visitor&& visitor) const {
            for (int value : contour(node)) {
                visitor(value);
            }
        }

        /**
         * @brief Constrói um vetor com o contorno agregado do nó.
         */
        std::vector<int> buildContourVector(NodeId node) const {
            std::vector<int> values;
            values.reserve(static_cast<std::size_t>(contours.sizeOf(node)));
            forEachContourPixel(node, [&](int px) { values.push_back(px); });
            return values;
        }

        /**
         * @brief Copia o contorno agregado para um iterador de saída (ex.: `back_inserter`).
         */
        template <typename OutputIterator>
        void copyContour(NodeId node, OutputIterator out) const {
            for (int value : contour(node)) {
                *out++ = value;
            }
        }

    private:
        /**
         * @brief Constrói, quando necessário, a arena agregada com contornos definitivos.
         *
         * A rotina percorre a árvore em pós-ordem, acumula os pixels herdados dos filhos,
         * adiciona os pixels do próprio nó e remove aqueles marcados em `removals`. Um
         * bitmap temporário evita duplicidades e permite descartar remoções em O(1).
         * O resultado final é gravado em `aggregated`, deixando a iteração subsequente
         * livre de novas alocações.
         */
        void ensureAggregated() const {
            if (aggregatedReady_) {
                return;
            }

            aggregated = ListArena(tree.getNumInternalNodeSlots(), static_cast<int>(contours.entries.size()));
            std::vector<std::vector<int>> accumulator(tree.getNumInternalNodeSlots());
            // bitmap linear (1 byte por pixel) utilizado para deduplicar/remover em O(1)
            std::vector<uint8_t> bitmap(tree.getNumRowsOfImage() * tree.getNumColsOfImage(), 0);

            // percorre em pós-ordem para propagar primeiro os filhos
            auto traversal = tree.getPostOrderNodes();
            for (NodeId nodeId : traversal) {
                const NodeId node = nodeId;
                auto& values = accumulator[node];

                //acumula contornos dos filhos
                for (NodeId childNodeId : tree.getChildren(nodeId)) {
                    const NodeId child = childNodeId;
                    auto& childValues = accumulator[child];
                    values.insert(values.end(), childValues.begin(), childValues.end());
                    childValues.clear();
                }
                //adiciona contornos do próprio nó
                for (int value : contours.range(node)) {
                    values.push_back(value);
                }

                //remove duplicatas locais de values: “compactar” o vetor no próprio lugar
                std::size_t writeIndex = 0;
                for (int value : values) {
                    if (!bitmap[value]) {
                        bitmap[value] = 1;
                        values[writeIndex++] = value; 
                    }
                }
                values.resize(writeIndex);

                for (int rem : removals.range(node)) {
                    if (rem >= 0 && rem < static_cast<int>(bitmap.size())) {
                        bitmap[rem] = 0;
                    }
                }

                writeIndex = 0;
                for (int value : values) {
                    if (bitmap[value]) {
                        aggregated.add(node, value);
                        values[writeIndex++] = value;
                        bitmap[value] = 0;
                    }
                }
                values.resize(writeIndex);
            }

            aggregatedReady_ = true;
        }
    };

    /**
     * @brief Executa a computação incremental e devolve contornos compactados.
     *
     * @param tree árvore morfológica (máx-tree, mín-tree ou ToS) sobre a qual o cálculo será realizado.
     * @return Estrutura `IncrementalContours` contendo arenas para acesso aos contornos.
     *
     * Exemplo:
     * @code
     * auto contours = ContoursComputedIncrementally::extractCompactContours(tree);
     * auto proxy = contours.contour(nodeId);
     * std::vector<int> pixels(proxy.begin(), proxy.end());
     * @endcode
     */
    static IncrementalContours extractCompactContours(MorphologicalTree& tree) {
        const int numNodes = tree.getNumInternalNodeSlots();
        const int totalPixels = tree.getNumRowsOfImage() * tree.getNumColsOfImage();
        const AdjacencyRelation* adjacencyContext = tree.getAdjacencyRelation();
        if (adjacencyContext == nullptr) {
            throw std::invalid_argument("Contour extraction requires an adjacency relation.");
        }

        // estrutura final que será exposta ao chamador (contornos + remoções)
        IncrementalContours result(tree, numNodes, std::max(totalPixels / 4, 1));
        // arena temporária que carrega pixels até o LCA correto antes da remoção
        ListArena contoursToRemoveLCA(numNodes, std::max(totalPixels / 4, 1));

        // contador auxiliar para saber quando um pixel deixa de ser contorno (cf. artigo original)
        std::vector<int> ncount(totalPixels, 0);
        // reuso de armazenamento para pixels que devem ser removidos neste nó
        std::vector<int> removalBuffer;
        removalBuffer.reserve(64);
        AdjacencyRelation adj4(tree.getNumRowsOfImage(), tree.getNumColsOfImage(), 1);
        AttributeComputedIncrementally::traversePostOrder(
            tree,
            tree.getRoot(),
            [](NodeId) -> void {},
            [](NodeId, NodeId) -> void {},
            [&](NodeId nodePId) {
                const NodeId nodeP = nodePId;
                // remove e processa todas as remoções pendentes deste nó
                removalBuffer.clear();
                contoursToRemoveLCA.consumeInto(nodeP, removalBuffer);
                for (int p : removalBuffer) {
                    bool isPixelToBeRemoved = true;
                    for (int r : adj4.getNeighborPixels(p)) {
                        NodeId nodeRId = tree.getSmallestComponent(r);
                        if (tree.isStrictAncestor(nodeRId, nodePId)) {
                            contoursToRemoveLCA.add(nodeRId, p);
                            isPixelToBeRemoved = false;
                        } else if (!tree.isComparable(nodePId, nodeRId)) {
                            NodeId otherNodeLCA = tree.getLowestCommonAncestor(nodePId, nodeRId);
                            contoursToRemoveLCA.add(otherNodeLCA, p);
                            isPixelToBeRemoved = false;
                        }
                    }
                    if (!adj4.isBorderDomainImage(p) && isPixelToBeRemoved) {
                        result.removals.add(nodeP, p);
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
                            contoursToRemoveLCA.add(nodeLCA, p);
                            ncount[p]++;
                        } else if (tree.isStrictDescendant(nodePId, nodeQId)) { // pixel ainda é fronteira
                            ncount[p]++;
                        } else if (tree.isStrictAncestor(nodePId, nodeQId)) {
                            ncount[q]--;
                            if (ncount[q] == 0) {
                                result.removals.add(nodeP, q);
                            }
                        }
                    }

                    if (ncount[p] > 0) {
                        result.contours.add(nodeP, p);
                    }
                }
            });

        return result;
    }
};

} // namespace mmcfilters
