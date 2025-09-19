#include <list>
#include <vector>
#include <stack>
#include <queue>
#include <iterator>
#include <utility>
#include <functional> 
#include <optional>
#include <stdexcept>
#include <type_traits>

#ifndef NODECT_H
#define NODECT_H

#include "../include/Common.hpp"
#include "../include/ComponentTree.hpp"


/**
 * @brief Handle (proxy leve) para nós da Component Tree com acesso O(1) via NodeArena.
 *
 * `NodeCT` é um *wrapper* leve que referencia um nó da árvore de componentes
 * (identificado por `id`) e delega leituras/escritas para a `NodeArena` da
 * árvore dona (`ComponentTree<CNPsType>* tree`). Não possui propriedade do nó nem
 * alocação própria — serve como *view* handle seguro e barato de copiar.
 *
 * ## Principais capacidades
 * - **Acesso O(1)** aos atributos do nó (nível, área, representante, CNPs).
 * - **Relações estruturais**: `getParent()`, `addChild()`, `removeChild()`,
 *   `spliceChildren()` e as variantes por `Id` otimizadas (operam só na arena).
 * - **Ranges/Iteradores** para percursos:
 *   - `getChildren()` → range de filhos (por *proxy* `NodeCT`).
 *   - `getNodesOfPathToRoot()` → caminho até a raiz.
 *   - `getIteratorPostOrderTraversal()` → pós-ordem.
 *   - `getIteratorBreadthFirstTraversal()` → largura (BFS).
 *   - `getRepsOfCC()` → itera representantes (CNPs) em BFS na subárvore.
 * - **Métricas**: `computerNumDescendants()` e `computerNumFlatzoneDescendants()`.
 *
 * ## Semântica de tipo (CNPsType)
 * - `Pixels` (alias para `int`): cada nó guarda um único representante; métodos
 *   como `addRepCNPs(int)` substituem o valor.
 * - `FlatZones` (ex.: `std::vector<int>`): cada nó pode armazenar vários reps;
 *   há métodos especializados (habilitados por SFINAE) para adicionar/remover
 *   e consultar contagens de flat-zones.
 *
 * ## Invariantes
 * - Este handle é válido se e somente se `tree != nullptr` e `id >= 0`.
 * - As operações que alteram a estrutura (ex.: `addChild`, `setParent`) são
 *   encaminhadas à `ComponentTree`, que mantém a coerência dos vetores da arena.
 *
 * ## Complexidade (típica)
 * - Acesso/atribuição de campo: O(1).
 * - Inserção/remoção de relação pai↔filho (arena-based): O(1).
 * - Iteração por ranges/iteradores: O(k) no número de nós/itens visitados.
 *
 * ## Exemplo mínimo
 * @code
 * ComponentTree<Pixels> T = ... ;
 * NodeCT<Pixels> n = T.proxy(rootId);
 * if (n) {
 *     int lvl = n.getLevel();
 *     for (auto c : n.getChildren()) {
 *         // processa filhos por proxy
 *     }
 *     for (int rep : n.getRepsOfCC()) {
 *         // reps em BFS da subárvore de n
 *     }
 * }
 * @endcode
 *
 */
class NodeCT{
private:
    friend class ComponentTree;

    NodeId id = -1;
    ComponentTree* tree = nullptr;  // árvore dona
    
public:
    
    // Construtor padrão
    NodeCT() = default;
    NodeCT(ComponentTree* owner, NodeId i) : id(i), tree(owner) {}

    //id >= 0 é o operador booleno
    explicit operator bool() const noexcept { return tree != nullptr && id >= 0; }

    //id é o operador == e !=
    bool operator==(const NodeCT& other) const noexcept { return tree == other.tree && id == other.id; }
    bool operator!=(const NodeCT& other) const noexcept { return !(*this == other); }

    //id é o operador int
    operator int() const noexcept { return id; }

    //acesso e escrita de dados na arena
    inline int getRepNode() const noexcept{ return tree->arena.repNode[id]; }
    inline int getIndex() const noexcept{ return id; }
    inline int getLevel() const noexcept{ return tree->arena.threshold2[id]; }
    inline int32_t getArea() const noexcept{ return tree->arena.areaCC[id]; }
    inline int getNumChildren() const noexcept{ return tree->arena.childCount[id]; }
    inline int getNumSiblings() const noexcept{ return tree->getNumSiblingsById(id); }
    inline void setArea(int32_t a) noexcept { tree->arena.areaCC[id] = a; }
    inline void setLevel(int lv) noexcept { tree->arena.threshold2[id] = lv; }
    
    inline int getTimePostOrder() const noexcept{ return tree->arena.timePostOrder[id]; }
    inline void setTimePostOrder(int time) noexcept{ tree->arena.timePostOrder[id] = time; } 
    inline int getTimePreOrder() const noexcept{ return tree->arena.timePreOrder[id]; }
    inline void setTimePreOrder(int time)noexcept { tree->arena.timePreOrder[id] = time; }

    //propaga para versão por id 
    inline auto getPixelsOfCC() const { return tree->getPixelsOfCCById(id); } //iterador
    inline auto getCNPs() const { return tree->getCNPsById(id); } //iterador 
    inline auto getNodesDescendants() const { return tree->getNodesDescendantsById(id); } //iterador
    
    inline int getNumCNPs() const { return tree->getNumCNPsById(id);} 
    inline bool isChild(NodeCT node) const noexcept { return tree && node && tree->hasChildById(id, node.getIndex()); }
    inline NodeCT getParent() { if (!tree || tree->arena.parentId[id] < 0) return {}; return NodeCT(tree, tree->arena.parentId[id]);}
    inline void setParent(NodeCT node) { if (!tree) return;  tree->setParentById(id, node.getIndex()); }
    inline void addChild(NodeCT child) { if (!tree || !child) return; tree->addChildById(id, child.getIndex()); }
    inline void removeChild(NodeCT child, bool releaseNode) { if (!tree || !child) return; tree->removeChildById(id, child.getIndex(), releaseNode); }
    inline void spliceChildren(NodeCT from) { if (!tree || !from || from.getIndex() == id) return; tree->spliceChildrenById(id, from.getIndex()); }
    inline bool isLeaf() const { return tree->isLeafById(id); }
    inline int getResidue() const noexcept{ return tree->getResidueById(id); }
    inline bool isMaxtreeNode() const noexcept{ return tree->isMaxtreeNodeById(id); }

          
    // Conta descendentes (exclui o próprio nó)
    int getNumDescendants() {
        return tree->getNumDescendantsById(id);
    }


    // Ranges existentes que devolvem filhos (por ponteiro lógico) continuam,
    // mas internamente usam `tree->proxy(id)` (que agora devolve handle)
    class ChildIdRange {
        int cur; ComponentTree* T;
    public:
        struct It {
            int id; ComponentTree* T;
            bool operator!=(const It& o) const { return id != o.id; }
            void operator++() { id = (id == -1) ? -1 : T->arena.nextSiblingId[id]; }
            NodeCT operator*() const { return T->proxy(id); }
        };
        It begin() const { return {cur, T}; }
        It end()   const { return {-1,  T}; }
        ChildIdRange(int first, ComponentTree* t) : cur(first), T(t) {}
    };
    auto getChildren() const {
        return ChildIdRange(tree ? tree->arena.firstChildId[id] : -1, tree);
    }

    //============= Iterator para iterar os nodes do caminho até o root==============//
    class InternalIteratorNodesOfPathToRoot {
    private:
        ComponentTree* T_ = nullptr;
        NodeId curId_ = -1;
    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = NodeCT;
        using difference_type = std::ptrdiff_t;
        using pointer = NodeCT*;
        using reference = NodeCT; // retornamos por valor (handle leve)

        InternalIteratorNodesOfPathToRoot(NodeCT obj) : T_(obj.tree), curId_(obj ? obj.getIndex() : -1) {}

        InternalIteratorNodesOfPathToRoot& operator++() {
            if (T_ && curId_ != -1) {
                curId_ = T_->getParentById(curId_);
            }
            return *this;
        }

        bool operator==(const InternalIteratorNodesOfPathToRoot& other) const {
            return curId_ == other.curId_;
        }
        bool operator!=(const InternalIteratorNodesOfPathToRoot& other) const { return !(*this == other); }

        value_type operator*() const { return (curId_ == -1) ? NodeCT() : NodeCT(T_, curId_); }
    };
    
    class IteratorNodesOfPathToRoot {
    private:
        ComponentTree* T_ = nullptr;
        NodeId startId_ = -1;
    public:
        explicit IteratorNodesOfPathToRoot(NodeCT obj) : T_(obj.tree), startId_(obj ? obj.getIndex() : -1) {}
        InternalIteratorNodesOfPathToRoot begin() const { return InternalIteratorNodesOfPathToRoot(NodeCT(T_, startId_)); }
        InternalIteratorNodesOfPathToRoot end() const { return InternalIteratorNodesOfPathToRoot(NodeCT()); }
    };
    
    // Chamador usa this como shared_ptr:
    IteratorNodesOfPathToRoot getNodesOfPathToRoot() { return IteratorNodesOfPathToRoot(NodeCT(tree, id)); }
    



    //============= Iterator para iterar os nodes de um percuso em pos-ordem ==============//
    class InternalIteratorPostOrderTraversal {
    private:
        struct Item { int id; bool expanded; };

        ComponentTree* T_ = nullptr;
        FastStack<Item> st_;
        NodeId current_ = -1;

        inline void settle_() noexcept {
            while (!st_.empty()) {
                Item &top = st_.top();
                if (!top.expanded) {
                    top.expanded = true;
                    // empilha filhos (direita→esquerda na prática; inverta se quiser L→R)
                    for (int c = T_->arena.firstChildId[top.id]; c != -1; c = T_->arena.nextSiblingId[c]) {
                        st_.push(Item{c, false});
                    }
                } else {
                    current_ = st_.top().id;  // todos os filhos já emitidos
                    return;
                }
            }
            current_ = -1; // fim
        }

    public:
        using iterator_category = std::input_iterator_tag;
        using value_type        = NodeCT;
        using difference_type   = std::ptrdiff_t;
        using pointer           = NodeCT*;     // apenas para conformidade
        using reference         = NodeCT;      // retornaremos por valor (handle leve)

        InternalIteratorPostOrderTraversal(NodeCT root) noexcept {
            if (root) {
                T_ = root.tree;
                st_.push({root.getIndex(), false});
                settle_();
            } else {
                current_ = -1;
            }
        }

        inline InternalIteratorPostOrderTraversal& operator++() noexcept {
            if (!st_.empty()) st_.pop();  // consome o atual
            settle_();                    // posiciona no próximo
            return *this;
        }

        // devolve o PROXY (handle leve) para o nó atual
        inline value_type operator*() const noexcept {
            return (current_ >= 0) ? NodeCT(T_, current_) : NodeCT();
        }

        inline bool operator==(const InternalIteratorPostOrderTraversal& other) const noexcept {
            const bool aEnd = (current_ == -1);
            const bool bEnd = (other.current_ == -1);
            if (aEnd || bEnd) return aEnd == bEnd;
            // fora do fim, basta comparar o id atual
            return current_ == other.current_;
        }
        inline bool operator!=(const InternalIteratorPostOrderTraversal& other) const noexcept {
            return !(*this == other);
        }
    };

    class IteratorPostOrderTraversal {
        ComponentTree* T_ = nullptr; NodeId rootId_ = -1;
    public:
        explicit IteratorPostOrderTraversal(NodeCT root) : T_(root.tree), rootId_(root ? root.getIndex() : -1) {}
        InternalIteratorPostOrderTraversal begin() const { return InternalIteratorPostOrderTraversal(NodeCT(T_, rootId_)); }
        InternalIteratorPostOrderTraversal end()   const { return InternalIteratorPostOrderTraversal(NodeCT()); }
    };

    IteratorPostOrderTraversal getIteratorPostOrderTraversal() {
        return IteratorPostOrderTraversal(NodeCT(tree, id));
    }



    //============= Iterator para iterar os nodes de um percuso em largura ==============//
    class InternalIteratorBreadthFirstTraversal {
    private:
        ComponentTree* T_ = nullptr; // apenas leitura/encaminhamento
        FastQueue<int> q_;                     // guarda somente NodeIds

    public:
        using iterator_category = std::input_iterator_tag;
        using value_type        = NodeCT;    // devolvemos PROXY
        using difference_type   = std::ptrdiff_t;
        using pointer           = void;               // não expomos ponteiro real
        using reference         = NodeCT;   // retornamos por valor (handle leve)

        // Constrói a partir de um proxy raiz. Enfileira apenas o id da raiz.
        explicit InternalIteratorBreadthFirstTraversal(NodeCT root) noexcept {
            if (root) {
                T_ = root.tree;
                q_.push(root.getIndex());
            }
        }

        // Pré-incremento: consome o nó da frente e enfileira seus filhos por ID
        InternalIteratorBreadthFirstTraversal& operator++() noexcept {
            if (!q_.empty()) {
                int u = q_.front();
                q_.pop();
                for(int c: T_->arena.children(u))
                    q_.push(c);
            }
            return *this;
        }

        // Desreferencia: devolve um PROXY (NodeCT) para o nó atual (frente da fila)
        value_type operator*() const noexcept {
            return q_.empty() ? NodeCT() : NodeCT(T_, q_.front());
        }

        bool operator==(const InternalIteratorBreadthFirstTraversal& other) const noexcept {
            return q_.empty() == other.q_.empty();
        }
        bool operator!=(const InternalIteratorBreadthFirstTraversal& other) const noexcept {
            return !(*this == other);
        }
    };

    class IteratorBreadthFirstTraversal {
    private:
        ComponentTree* T_ = nullptr; NodeId rootId_ = -1;
    public:
        explicit IteratorBreadthFirstTraversal(NodeCT root) : T_(root.tree), rootId_(root ? root.getIndex() : -1) {}
        InternalIteratorBreadthFirstTraversal begin() const { return InternalIteratorBreadthFirstTraversal(NodeCT(T_, rootId_)); }
        InternalIteratorBreadthFirstTraversal end()   const { return InternalIteratorBreadthFirstTraversal(NodeCT()); }
    };

    // Método para expor o iterador na classe NodeCT
    IteratorBreadthFirstTraversal getIteratorBreadthFirstTraversal() {
        return IteratorBreadthFirstTraversal(NodeCT(tree, id));
    }
    
};



#endif
