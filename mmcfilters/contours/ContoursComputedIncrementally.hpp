#pragma once

/*
 * Visão Geral
 * -----------
 * Aqui implementamos uma versão baseada em arena para o cálculo
 * incremental de contornos em árvores morfológicas. O desenho nasceu das
 * discussões sobre performance e reaproveitamento de memória e captura alguns
 * pontos principais:
 *
 * 1. ListArena
 *    Arena leve que armazena várias listas simplesmente encadeadas dentro de um
 *    único buffer contíguo (`entries`). Cada nó da árvore possui a sua lista
 *    (controlada por `head`). Com isso, evitamos alocações dinâmicas
 *    frequentes enquanto adicionamos ou removemos pixels durante a travessia em
 *    pós-ordem.
 *
 * 2. recycle() e consumeInto()
 *    Mesmo que um pixel removido jamais volte a ser contorno na mesma execução,
 *    outros pixels ou nós ainda precisarão de novos slots. Sem reciclagem, cada
 *    chamada a `add()` faria um `push_back` e deixaria entradas antigas sem uso,
 *    fazendo o vetor crescer com o número total de operações e não apenas com o
 *    pico simultâneo. `consumeInto()` drena uma lista, entrega seus valores a um
 *    recipiente temporário e devolve os nós ao freelist via `recycle()`,
 *    permitindo reutilização barata (O(1) e amigável ao cache) em inserções
 *    futuras.
 *
 * 3. extractCompactContours()
 *    Durante a travessia incremental criamos uma instância de
 *    `IncrementalContours` (referida como `result`). Ela guarda as arenas que
 *    serão expostas ao usuário (`contours`, `removals`) e que poderão ser
 *    iteradas depois. O auxiliar `contoursToRemoveLCA` é outra `ListArena`
 *    temporária: quando descobrimos que um pixel precisa ser removido em um
 *    ancestral (normalmente o LCA de dois nós incomparáveis), enfileiramos esse
 *    pixel até que o ancestral seja visitado.
 *
 *    Em cada nó drenamos (`consumeInto`) as remoções pendentes daquele nó,
 *    propagamos para ancestrais quando necessário ou registramos a remoção
 *    definitiva em `result.removals`. Os contornos válidos são adicionados a
 *    `result.contours`. Mais tarde, quando o chamador pede o contorno
 *    consolidado (`contour(node)`), a agregação lazy combina os pixels locais
 *    com os dos filhos, subtraindo as remoções acumuladas.
 *
 * 4. Etapa de agregação
 *    O método `ensureAggregated()` realiza uma segunda passagem em pós-ordem
 *    assim que algum iterador público é requisitado. Ele acumula os pixels por
 *    nó, aplica as remoções e elimina duplicidades usando um bitmap temporário
 *    bem leve (a ordem final corresponde à ordem de descoberta, não à raster).
 *    O resultado é gravado em outra arena (`aggregated`). Dessa forma, leituras
 *    subsequentes ficam rápidas e pagamos o custo de consolidação apenas uma
 *    vez.
 *
 * Em resumo: `ListArena` evita alocações durante a fase incremental, `recycle`
 * impede o crescimento desnecessário da arena, `result` guarda o resultado
 * visível pelo usuário e `contoursToRemoveLCA` carrega remoções adiadas até que
 * o ancestral correto seja processado.
 */

#include "../utils/Common.hpp"
#include "../trees/MorphologicalTree.hpp"
#include "../attributes/AttributeComputedIncrementally.hpp"
#include "../utils/AdjacencyRelation.hpp"

namespace mmcfilters {

/**
 * @brief Estruturas auxiliares baseadas em arena para armazenar contornos incrementais.
 */
class ContoursComputedIncrementally {
public:
    /**
     * @brief Arena que guarda diversas listas encadeadas de inteiros em memória contígua.
     */
    struct ListArena {
        /**
         * @brief Nó simples da lista: valor + ponteiro (índice) para o próximo elemento.
         */
        struct Entry {
            int value;
            int next;
        };

        /**
         * @brief Iterador forward constante sobre uma lista armazenada no arena.
         */
        /**
         * @brief Iterador forward que percorre uma lista pertencente ao arena.
         *
         * Implementa interface de forward iterator padrão permitindo uso em
         * range-based for (`for (int v : arena.range(node))`).
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
         * @brief View helper que expõe begin/end para um dado identificador de lista.
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

        /// Cria arena vazia sem reservar memória.
        ListArena() = default;

        std::vector<Entry> entries;
        std::vector<int> head;
        std::vector<int> size;
        int freeHead = -1;

        /// Cria arena com `numNodes` listas potenciais e hint opcional de capacidade.
        explicit ListArena(int numNodes, int capacityHint = 0): head(numNodes, -1), size(numNodes, 0){
            if (capacityHint > 0) {
                entries.reserve(capacityHint);
            }
        }

        /// Remove todos os elementos associados a `node`, devolvendo-os ao freelist.
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

        /// Insere `value` na lista pertencente a `node`.
        void add(int node, int value) {
            const int slot = allocate();
            entries[slot] = Entry{value, head[node]};
            head[node] = slot;
            ++size[node];
        }

        Range range(int node) const { return Range(this, head[node]); }

        /**
         * @brief Move todos os elementos da lista `node` para `out` reutilizando os slots.
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

        /// @return número de elementos atualmente associados a `node`.
        int sizeOf(int node) const { return size[node]; }

        const_iterator begin(int node) const { return const_iterator(this, head[node]); }
        const_iterator end() const { return const_iterator(this, -1); }

    private:
        /// Obtém um slot livre (do freelist ou via push_back).
        int allocate() {
            if (freeHead == -1) {
                entries.push_back(Entry{0, -1});
                return static_cast<int>(entries.size() - 1);
            }
            const int idx = freeHead;
            freeHead = entries[idx].next;
            return idx;
        }

        /// Devolve um índice para o freelist.
        void recycle(int idx) {
            entries[idx].next = freeHead;
            freeHead = idx;
        }
    };

    /**
     * @brief Resultado da computação incremental armazenado em arenas.
     *
     * Contém os contornos crus, as informações de remoção e a arena agregada construída
     * sob demanda. As funções de utilidade expõem iteração conveniente sobre os pixels
     * sem obrigar o chamador a conhecer detalhes internos de armazenamento.
     */
    struct IncrementalContours {
        const MorphologicalTree* tree;
        ListArena contours;
        ListArena removals;
        mutable ListArena aggregated;
        mutable bool aggregatedReady_ = false;

        /**
         * @param tree ponteiro para a árvore utilizada na computação.
         * @param numNodes número de nós (define o tamanho das arenas).
         * @param capacityHint sugestão para reserva inicial de entradas.
         */
        IncrementalContours(const MorphologicalTree* tree, int numNodes, int capacityHint)
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
                iterator(const IncrementalContours* owner, NodeId node): owner_(owner), node_(node) {}

                value_type operator*() const { return {node_, ContourProxy(owner_, node_)}; }

                iterator& operator++() {
                    ++node_;
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
                const IncrementalContours* owner_ = nullptr;
                NodeId node_ = 0;
            };

            explicit ContoursLazyRange(const IncrementalContours* owner)
                : owner_(owner) {}

            iterator begin() const { return iterator(owner_, 0); }
            iterator end() const { return iterator(owner_, owner_->tree->getNumNodes()); }

        private:
            const IncrementalContours* owner_ = nullptr;
        };

        /// @return proxy que permite iterar o contorno agregado do nó informado.
        ContourProxy contour(NodeId node) const { return ContourProxy(this, node); }

        /// @return range lazy sobre todos os nós da árvore, produzindo pares `(nodeId, proxy)`.
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

            aggregated = ListArena(tree->getNumNodes(), static_cast<int>(contours.entries.size()));
            std::vector<std::vector<int>> accumulator(tree->getNumNodes());
            // bitmap linear (1 byte por pixel) utilizado para deduplicar/remover em O(1)
            std::vector<uint8_t> bitmap(tree->getNumRowsOfImage() * tree->getNumColsOfImage(), 0);

            // percorre em pós-ordem para propagar primeiro os filhos
            auto traversal = const_cast<MorphologicalTree*>(tree)->getIteratorPostOrderTraversalById();
            for (NodeId node : traversal) {
                auto& values = accumulator[node];

                //acumula contornos dos filhos
                for (NodeId child : tree->getChildrenById(node)) {
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
     * auto contours = ContoursComputedIncrementally::extractCompactContours(tree.get());
     * auto proxy = contours.contour(nodeId);
     * std::vector<int> pixels(proxy.begin(), proxy.end());
     * @endcode
     */
    static IncrementalContours extractCompactContours(MorphologicalTree* tree) {
        const int numNodes = tree->getNumNodes();
        const int totalPixels = tree->getNumRowsOfImage() * tree->getNumColsOfImage();

        // estrutura final que será exposta ao chamador (contornos + remoções)
        IncrementalContours result(tree, numNodes, std::max(totalPixels / 4, 1));
        // arena temporária que carrega pixels até o LCA correto antes da remoção
        ListArena contoursToRemoveLCA(numNodes, std::max(totalPixels / 4, 1));

        // contador auxiliar para saber quando um pixel deixa de ser contorno (cf. artigo original)
        std::vector<int> ncount(totalPixels, 0);
        // reuso de armazenamento para pixels que devem ser removidos neste nó
        std::vector<int> removalBuffer;
        removalBuffer.reserve(64);
        auto adj4 = std::make_shared<AdjacencyRelation>(tree->getNumRowsOfImage(), tree->getNumColsOfImage(), 1);
        LCAEulerRMQ lca(tree);

        AttributeComputedIncrementally::computerAttribute(
            tree,
            tree->getRootById(),
            [](NodeId) -> void {},
            [](NodeId, NodeId) -> void {},
            [&](NodeId nodeP) {
                // remove e processa todas as remoções pendentes deste nó
                removalBuffer.clear();
                contoursToRemoveLCA.consumeInto(nodeP, removalBuffer);
                for (int p : removalBuffer) {
                    bool isPixelToBeRemoved = true;
                    for (int r : adj4->getNeighborPixels(p)) {
                        NodeId nodeR = tree->getSCById(r);
                        if (tree->isStrictAncestor(nodeR, nodeP)) {
                            contoursToRemoveLCA.add(nodeR, p);
                            isPixelToBeRemoved = false;
                        } else if (!tree->isComparable(nodeP, nodeR)) {
                            NodeId otherNodeLCA = lca.findLowestCommonAncestor(nodeP, nodeR);
                            contoursToRemoveLCA.add(otherNodeLCA, p);
                            isPixelToBeRemoved = false;
                        }
                    }
                    if (!adj4->isBorderDomainImage(p) && isPixelToBeRemoved) {
                        result.removals.add(nodeP, p);
                    }
                }

                // percorre os CNPs pertencentes ao nó atual
                for (int p : tree->getCNPsById(nodeP)) {
                    if (adj4->isBorderDomainImage(p)) {
                        ncount[p]++;
                    }

                    for (int q : adj4->getNeighborPixels(p)) {
                        NodeId nodeQ = tree->getSCById(q);
                    if (!tree->isComparable(nodeP, nodeQ)) { // contorno será tratado pelo LCA
                        NodeId nodeLCA = lca.findLowestCommonAncestor(nodeP, nodeQ);
                        contoursToRemoveLCA.add(nodeLCA, p);
                        ncount[p]++;
                    } else if (tree->isStrictDescendant(nodeP, nodeQ)) { // pixel ainda é fronteira
                        ncount[p]++;
                    } else if (tree->isStrictAncestor(nodeP, nodeQ)) {
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
