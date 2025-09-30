#pragma once

#include "../utils/AdjacencyRelation.hpp"
#include "../utils/Common.hpp"

namespace mmcfilters {

// Forward declaration do builder externo
class NodeMT;
class MorphologicalTree;

/**
 * @brief Arena de nós para Morphological trees com armazenamento contíguo e acesso rápido.
 *
 * A `NodeMTArena<CNPsType>` gerencia os atributos estruturais e dados de cada nó da árvore
 * de componentes (Morphological tree) em vetores paralelos, proporcionando eficiência de cache
 * e operações O(1) para acesso a campos.
 *
 * ## Estrutura de dados
 * Cada vetor armazena um atributo específico dos nós:
 * - `repNode`: representante do conjunto de pixels (Union-Find).
 * - `threshold2`: nível (gray-level ou threshold máximo).
 * - `areaCC`: área acumulada do componente conexo.
 * - `repCNPs`: representantes de CNPs (pixels ou flat zones, dependendo do tipo template).
 * - Ponteiros estruturais: `parentId`, `firstChildId`, `nextSiblingId`, `prevSiblingId`, `lastChildId`.
 * - `childCount`: número de filhos diretos (cache).
 *
 * ## Operações principais
 * - `allocate(rep, thr1, thr2)`: cria um novo nó e inicializa seus campos com valores padrão.
 * - `reserve(n)`: reserva espaço para `n` nós, evitando realocações.
 * - `size()`: retorna o número de nós já alocados.
 * - `children(id)`: retorna um range leve (`ChildRange`) para iterar sobre os filhos de um nó.
 * - `getRepsOfCC(id)`: retorna um range BFS (`RepsOfCCRangeById`) que percorre todos os representantes
 *   da subárvore enraizada em `id`.
 * - **Reutilização de IDs**: `releaseNode(id)` libera o slot e `allocate(...)` reaproveita IDs disponíveis sem crescer os vetores.
 *
 * ## Iteradores
 * - `ChildRange`: itera os filhos diretos de um nó via `range-for`.
 * - `RepsOfCCIteratorById`: percorre em BFS os representantes armazenados nos nós de uma CC.
 *
 * ## Exemplo
 * @code
 * NodeMTArena arena;
 * arena.reserve(1000);
 * NodeId root = arena.allocate(rep, thr2);
 * for (NodeId c : arena.children(root)) {
 *     // processa cada filho
 * }
 * for (int rep : arena.getRepsOfCC(root)) {
 *     // processa reps da subárvore
 * }
 * @endcode
 *
 */
class NodeMTArena {
    private:
    friend class MorphologicalTree; 
    friend class NodeMT; 

    std::vector<NodeId>   repNode;       // representante do UF
    std::vector<int>      threshold2;    // level (max threshold)
    std::vector<int32_t>  areaCC;
    
    std::vector<int> timePostOrder;  //tempo de entrada durante o percurso pos-ordem (2 incremento) 
    std::vector<int> timePreOrder;  //tempo de saia durante o percurso pos-ordem (2 incremento) 
	
    
    std::vector<NodeId>   parentId;      // InvalidNode = raiz
    std::vector<NodeId>   firstChildId;  // InvalidNode = sem filhos
    std::vector<NodeId>   nextSiblingId; // InvalidNode = sem próximo
    std::vector<NodeId>   prevSiblingId; // InvalidNode = sem anterior
    std::vector<NodeId>   lastChildId;   // InvalidNode = sem filhos
    std::vector<int>      childCount;    // cache para ter acesso a quantidade de filhos diretos
    
    // Lista de IDs livres para reutilização (LIFO)
    std::vector<NodeId>   freeIds;

    public:
   
    // --- gerenciamento ---
    inline NodeId allocate(NodeId rep, int thr2) {
        // 1) Reutiliza ID livre, se houver
        if (!freeIds.empty()) {
            NodeId id = freeIds.back();
            freeIds.pop_back();

            // Reinicializa o slot existente
            repNode[id]      = rep;              // representante UF
            threshold2[id]   = thr2;             // level
            areaCC[id]       = 0;
            
            parentId[id]     = InvalidNode;
            firstChildId[id] = InvalidNode;
            nextSiblingId[id]= InvalidNode;
            prevSiblingId[id]= InvalidNode;
            lastChildId[id]  = InvalidNode;
            childCount[id]   = 0;
            timePreOrder[id] = -1;
            timePostOrder[id] = -1;
            return id;
        }

        // 2) Caso contrário, aloca um novo slot no final (comportamento antigo)
        NodeId id = static_cast<NodeId>(repNode.size());
        repNode.push_back(rep);
        threshold2.push_back(thr2);
        areaCC.push_back(0);
        parentId.push_back(InvalidNode);
        firstChildId.push_back(InvalidNode);
        nextSiblingId.push_back(InvalidNode);
        prevSiblingId.push_back(InvalidNode);
        lastChildId.push_back(InvalidNode);
        childCount.push_back(0);
        timePreOrder.push_back(-1);
        timePostOrder.push_back(-1);
        return id;
    }

    //aloca espaço para n nós (sem inicialização)
    inline void reserve(size_t n) {
        repNode.reserve(n); threshold2.reserve(n); areaCC.reserve(n);
        parentId.reserve(n); firstChildId.reserve(n); nextSiblingId.reserve(n); 
        prevSiblingId.reserve(n); lastChildId.reserve(n); childCount.reserve(n);
        timePreOrder.reserve(n); timePostOrder.reserve(n);
    }

    // Marca o ID como livre para reutilização. Pré-condição: nó já está desconectado (sem pai e sem filhos).
    inline void releaseNode(NodeId id) noexcept {
        // Zera/normaliza campos observáveis
        repNode[id]      = InvalidNode;      // marcador de slot livre 
        threshold2[id]   = 0;
        areaCC[id]       = 0;
        parentId[id]     = InvalidNode;
        firstChildId[id] = InvalidNode;
        nextSiblingId[id]= InvalidNode;
        prevSiblingId[id]= InvalidNode;
        lastChildId[id]  = InvalidNode;
        childCount[id]   = 0;
        timePreOrder[id] = -1;
        timePostOrder[id] = -1;

        freeIds.push_back(id);
    }

    // Consulta simples: retorna true se o slot parece livre
    inline bool isFree(NodeId id) const noexcept {
        return id >= 0 && id < static_cast<NodeId>(repNode.size()) && repNode[id] == InvalidNode;
    }

    size_t size() const { return repNode.size(); }


    // ============================================
    // ADIÇÃO: Range leve para filhos (range-for)
    // ============================================
    /**
     * @brief Range leve para iterar sobre os filhos diretos de um nó.
     */
    class ChildRange {
    public:
        // iterador 'input' mínimo para range-for
        /**
         * @brief Iterador que percorre sequencialmente os filhos de um nó.
         */
        class iterator {
        public:
            iterator(NodeId cur, const NodeMTArena* arena) noexcept
            : cur_(cur), arena_(arena) {}

            NodeId operator*() const noexcept { return cur_; }

            iterator& operator++() noexcept {
                cur_ = (cur_ == InvalidNode) ? InvalidNode : arena_->nextSiblingId[cur_];
                return *this;
            }

            bool operator!=(const iterator& other) const noexcept {
                return cur_ != other.cur_;
            }

        private:
            NodeId cur_;
            const NodeMTArena* arena_;
        };

        ChildRange(NodeId first, const NodeMTArena* arena) noexcept
        : first_(first), arena_(arena) {}

        iterator begin() const noexcept { return iterator(first_, arena_); }
        iterator end()   const noexcept { return iterator(InvalidNode,    arena_); }

        // açúcares úteis:
        bool empty() const noexcept { return first_ == InvalidNode; }
        NodeId front() const noexcept { return first_; }

    private:
        NodeId first_;
        const NodeMTArena* arena_;
    };

    // uso: for (NodeId c : arena.children(parentId)) { ... }
    inline ChildRange children(NodeId id) const noexcept {
        return ChildRange(firstChildId[id], this);
    }



    // Itera os representantes (ints) em BFS sobre a CC (subárvore) do nó.
    /**
     * @brief Iterador BFS que percorre representantes de componentes conectados.
     */
    class RepsOfCCIteratorById {
    private:
        const NodeMTArena* arena_ = nullptr;
        FastQueue<int> q_;               // BFS por IDs
        const int* curPtr_  = nullptr;    // ponteiro p/ bloco atual
        const int* curEnd_  = nullptr;
        int singleBuf_      = InvalidNode;         // buffer p/ CNPsType==Pixels

        void enqueueChildrenOf(int nid) {
            for(int c: arena_->children(nid)) 
                q_.push(c);
        }
        void loadBlockFromId(int nid) {
            singleBuf_ = arena_->repNode[nid];
            curPtr_ = &singleBuf_;
            curEnd_ = &singleBuf_ + 1;
        }
        void advanceToNextNodeWithReps() {
            curPtr_ = curEnd_ = nullptr;
            while (!q_.empty()) {
                int nid = q_.pop();
                // enfileira filhos primeiro (BFS)
                enqueueChildrenOf(nid);
                // carrega reps do nó corrente
                loadBlockFromId(nid);
                if (curPtr_ != curEnd_) return; // achou bloco não-vazio
            }
            // fim
            curPtr_ = curEnd_ = nullptr;
        }

    public:
        using iterator_category = std::input_iterator_tag;
        using value_type        = int;
        using difference_type   = std::ptrdiff_t;
        using pointer           = const int*;
        using reference         = const int&;

        RepsOfCCIteratorById(const NodeMTArena* arena, NodeId root, bool isEnd) {
            if (!isEnd && root >= 0) {
                arena_ = arena;
                // inicia fila com a RAIZ da CC (inclui reps do próprio nó)
                q_.push(root);
                advanceToNextNodeWithReps();
            }
        }

        reference operator*()  const { return *curPtr_; }
        pointer   operator->() const { return  curPtr_; }

        RepsOfCCIteratorById& operator++() {
            if (!curPtr_) return *this;         // já no fim
            ++curPtr_;
            if (curPtr_ == curEnd_) {           // acabou bloco do nó → próximo nó com reps
                advanceToNextNodeWithReps();
            }
            return *this;
        }

        RepsOfCCIteratorById operator++(int) { auto tmp = *this; ++(*this); return tmp; }

        bool operator==(const RepsOfCCIteratorById& other) const {
            const bool endA = (curPtr_ == nullptr && curEnd_ == nullptr && q_.empty());
            const bool endB = (other.curPtr_ == nullptr && other.curEnd_ == nullptr && other.q_.empty());
            if (endA || endB) return endA == endB;
            // estado interno diferente => não igual (suficiente p/ uso como input iterator)
            return false;
        }
        bool operator!=(const RepsOfCCIteratorById& other) const { return !(*this == other); }
    };

    /**
     * @brief Range para percorrer representantes de uma subárvore por BFS.
     */
    class RepsOfCCRangeById {
    private:
        const NodeMTArena* arena_;
        NodeId root_;
    public:
        explicit RepsOfCCRangeById(const NodeMTArena* arena, NodeId root): arena_(arena), root_(root) {}

        RepsOfCCIteratorById begin() const { return RepsOfCCIteratorById(arena_, root_, false); }
        RepsOfCCIteratorById end()   const { return RepsOfCCIteratorById(arena_, root_, true ); }
    };

    // Exponha um método público para uso direto em range-for:
    inline RepsOfCCRangeById getRepsOfCC(NodeId id) const {
        return RepsOfCCRangeById(this, id);
    }

};

} // namespace mmcfilters