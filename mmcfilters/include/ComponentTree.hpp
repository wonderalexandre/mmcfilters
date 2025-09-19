#include <list>
#include <vector>
#include <array>
#include <unordered_set>
#include <utility>
#include <optional>
#include <functional>
#include <variant>
#include <span>

#include "../include/AdjacencyRelation.hpp"
#include "../include/Common.hpp"


#ifndef COMPONENT_TREE_H
#define COMPONENT_TREE_H

// Forward declaration do builder externo
class IMorphologicalTreeBuilder; 
class NodeCT;

/**
 * @brief Arena de nós para Component Trees com armazenamento contíguo e acesso rápido.
 *
 * A `NodeArena<CNPsType>` gerencia os atributos estruturais e dados de cada nó da árvore
 * de componentes (component tree) em vetores paralelos, proporcionando eficiência de cache
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
 * NodeArena<std::vector<int>> arena;
 * arena.reserve(1000);
 * NodeId root = arena.allocate(rep, thr1, thr2);
 * for (NodeId c : arena.children(root)) {
 *     // processa cada filho
 * }
 * for (int rep : arena.getRepsOfCC(root)) {
 *     // processa reps da subárvore
 * }
 * @endcode
 *
 */
class NodeArena {
    private:
    friend class ComponentTree; 
    friend class NodeCT; 

    std::vector<int>      repNode;       // representante do UF
    std::vector<int>      threshold2;    // level (max threshold)
    //std::vector<int>    threshold1;    // min threshold
    std::vector<int32_t>  areaCC;
    
    std::vector<int> timePostOrder;  //tempo de entrada durante o percurso pos-ordem (2 incremento) 
    std::vector<int> timePreOrder;  //tempo de saia durante o percurso pos-ordem (2 incremento) 
	
    
    std::vector<NodeId>   parentId;      // -1 = raiz
    std::vector<NodeId>   firstChildId;  // -1 = sem filhos
    std::vector<NodeId>   nextSiblingId; // -1 = sem próximo
    std::vector<NodeId>   prevSiblingId; // -1 = sem anterior
    std::vector<NodeId>   lastChildId;   // -1 = sem filhos
    std::vector<int>      childCount;    // cache para ter acesso a quantidade de filhos diretos
    
    // Lista de IDs livres para reutilização (LIFO)
    std::vector<NodeId>   freeIds;

    public:
   
    // --- gerenciamento ---
    inline NodeId allocate(int rep, int thr2) {
        // 1) Reutiliza ID livre, se houver
        if (!freeIds.empty()) {
            NodeId id = freeIds.back();
            freeIds.pop_back();

            // Reinicializa o slot existente
            repNode[id]      = rep;              // representante UF
            threshold2[id]   = thr2;             // level
            areaCC[id]       = 0;
            
            parentId[id]     = -1;
            firstChildId[id] = -1;
            nextSiblingId[id]= -1;
            prevSiblingId[id]= -1;
            lastChildId[id]  = -1;
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
        parentId.push_back(-1);
        firstChildId.push_back(-1);
        nextSiblingId.push_back(-1);
        prevSiblingId.push_back(-1);
        lastChildId.push_back(-1);
        childCount.push_back(0);
        timePreOrder.push_back(-1);
        timePostOrder.push_back(-1);
        return id;
    }

    inline void reserve(size_t n) {
        repNode.reserve(n); threshold2.reserve(n); areaCC.reserve(n);
        parentId.reserve(n); firstChildId.reserve(n); nextSiblingId.reserve(n); 
        prevSiblingId.reserve(n); lastChildId.reserve(n); childCount.reserve(n);
        timePreOrder.reserve(n); timePostOrder.reserve(n);
    }

    // Marca o ID como livre para reutilização. Pré-condição: nó já está desconectado (sem pai e sem filhos).
    inline void releaseNode(NodeId id) noexcept {
        // Zera/normaliza campos observáveis
        repNode[id]      = -1;      // marcador de slot livre (não conflita com IDs válidos >= 0)
        threshold2[id]   = 0;
        areaCC[id]       = 0;
        parentId[id]     = -1;
        firstChildId[id] = -1;
        nextSiblingId[id]= -1;
        prevSiblingId[id]= -1;
        lastChildId[id]  = -1;
        childCount[id]   = 0;
        timePreOrder[id] = -1;
        timePostOrder[id] = -1;

        freeIds.push_back(id);
    }

    // Consulta simples: retorna true se o slot parece livre
    inline bool isFree(NodeId id) const noexcept {
        return id >= 0 && id < static_cast<NodeId>(repNode.size()) && repNode[id] == -1;
    }

    size_t size() const { return repNode.size(); }


    // ============================================
    // ADIÇÃO: Range leve para filhos (range-for)
    // ============================================
    class ChildRange {
    public:
        // iterador 'input' mínimo para range-for
        class iterator {
        public:
            iterator(NodeId cur, const NodeArena* arena) noexcept
            : cur_(cur), arena_(arena) {}

            NodeId operator*() const noexcept { return cur_; }

            iterator& operator++() noexcept {
                cur_ = (cur_ == -1) ? -1 : arena_->nextSiblingId[cur_];
                return *this;
            }

            bool operator!=(const iterator& other) const noexcept {
                return cur_ != other.cur_;
            }

        private:
            NodeId cur_;
            const NodeArena* arena_;
        };

        ChildRange(NodeId first, const NodeArena* arena) noexcept
        : first_(first), arena_(arena) {}

        iterator begin() const noexcept { return iterator(first_, arena_); }
        iterator end()   const noexcept { return iterator(-1,    arena_); }

        // açúcares úteis:
        bool empty() const noexcept { return first_ == -1; }
        NodeId front() const noexcept { return first_; }

    private:
        NodeId first_;
        const NodeArena* arena_;
    };

    // uso: for (NodeId c : arena.children(parentId)) { ... }
    inline ChildRange children(NodeId id) const noexcept {
        return ChildRange(firstChildId[id], this);
    }



    // Itera os representantes (ints) em BFS sobre a CC (subárvore) do nó.
    class RepsOfCCIteratorById {
    private:
        const NodeArena* arena_ = nullptr;
        FastQueue<int> q_;               // BFS por IDs
        const int* curPtr_  = nullptr;    // ponteiro p/ bloco atual
        const int* curEnd_  = nullptr;
        int singleBuf_      = -1;         // buffer p/ CNPsType==Pixels

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

        RepsOfCCIteratorById(const NodeArena* arena, NodeId root, bool isEnd) {
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

    class RepsOfCCRangeById {
    private:
        const NodeArena* arena_;
        NodeId root_;
    public:
        explicit RepsOfCCRangeById(const NodeArena* arena, NodeId root): arena_(arena), root_(root) {}

        RepsOfCCIteratorById begin() const { return RepsOfCCIteratorById(arena_, root_, false); }
        RepsOfCCIteratorById end()   const { return RepsOfCCIteratorById(arena_, root_, true ); }
    };

    // Exponha um método público para uso direto em range-for:
    inline RepsOfCCRangeById getRepsOfCC(NodeId id) const {
        return RepsOfCCRangeById(this, id);
    }

};


/**
 * @brief Estrutura de árvore de componentes (Component Tree) para imagens, com suporte a Pixels ou FlatZones.
 *
 * `ComponentTree<CNPsType>` organiza hierarquicamente regiões conexas de uma imagem
 * (componentes) em uma estrutura de árvore, permitindo análise multiescala e operações
 * de filtragem baseadas em atributos. O parâmetro de template `CNPsType` define se a
 * árvore é construída diretamente por pixels (`Pixels`) ou por flat-zones (`FlatZones`).
 *
 * ## Características principais
 * - **Construção**: via Union-Find (otimizado) usando pixels ou flat-zones.
 * - **Estrutura interna**: dados de nós armazenados em `NodeArena`, com acesso rápido O(1).
 * - **Proxy**: interface de nó exposta via `NodeCT<CNPsType>`, que encapsula acesso e
 *   relações pai/filho.
 *
 * ## Exemplo mínimo
 * @code
 * ImageUInt8Ptr img = ...;
 * auto adj = std::make_shared<AdjacencyRelation>(rows, cols, 1.5);
 * ComponentTree<Pixels> T(img, true, adj); // max-tree por pixels
 *
 * NodeCT<Pixels> root = T.getRoot();
 * for (auto c : root.getChildren()) {
 *     int area = c.getArea();
 * }
 * auto recon = T.reconstructionImage();
 * @endcode
 *
 * @tparam CNPsType Define o tipo de construção da árvore: `Pixels` ou `FlatZones`.
 */
class ComponentTree : public std::enable_shared_from_this<ComponentTree> {
    friend class BuilderComponentTreeByUnionFind; 
    friend class BuilderTreeOfShapeByUnionFind; 
    friend class NodeCT; 

protected:
    NodeId root;
    int numRows;
    int numCols;
    AdjacencyRelationPtr adj; //disk of a given ratio: ratio(1) for 4-connect and ratio(1.5) for 8-connect 
    int numNodes;
    int treeType; //0-mintree, 1-maxtree, 2-tree of shapes


    std::vector<NodeId> pixelToNodeId; //Mapeamento dos pixels representantes para NodeID. Para adquirir todos os representantes valido utilize o método getRepCNPs
    std::shared_ptr<PixelSetManager> pixelBuffer; PixelSetManager::View pixelView; //gerenciamento de pixels da arvore
    NodeArena arena; // armazenamento indexado dos dados de todos os nós 
    
    NodeId makeNode(int repNode, NodeId parentId, int threshold2);
    void reserveNodes(int expected) { arena.reserve(expected);}
    void reconstruction(NodeId node, uint8_t* data);
    void build(const ImageUInt8Ptr& imgPtr, IMorphologicalTreeBuilder& builderUF);

	void computerIncrementalAttributes(NodeId root, 
										std::function<void(NodeId)> preProcessing,
										std::function<void(NodeId, NodeId)> mergeChildren,
										std::function<void(NodeId)> postProcessing ){
		
		preProcessing(root);
		for(NodeId child: getChildrenById(root)){
			computerIncrementalAttributes(child, preProcessing, mergeChildren, postProcessing);
			mergeChildren(root, child);
		}
		postProcessing(root);
	}

    void computerTreeAttributes();
    ComponentTree() = default;
public:

   	static const int MAX_TREE = 0;
	static const int MIN_TREE = 1;
	static const int TREE_OF_SHAPES = 2;

	ComponentTree(ImageUInt8Ptr img, std::string ToSInperpolation="self-dual");
	explicit ComponentTree(ImageUInt8Ptr img, bool isMaxtree, double radius = 1.5);
    ComponentTree(ImageUInt8Ptr img, const char* ToSInterpolation) : ComponentTree(img, std::string(ToSInterpolation)) {}

    virtual ~ComponentTree() = default;

    static ComponentTreePtr create(int rows, int cols, bool isMaxtree, AdjacencyRelationPtr adj) {
        struct Enabler : public ComponentTree {
            Enabler(int r, int c, bool m, AdjacencyRelationPtr a) {
                this->numRows = r;
                this->numCols = c;
                this->treeType = m ? MAX_TREE : MIN_TREE;
                this->adj = a;
            }
        };
        return std::make_shared<Enabler>(rows, cols, isMaxtree, adj);
    }
    

    template <typename PixelType>
    static ComponentTreePtr createFromAttributeMapping(ImagePtr<PixelType> attrMappingPtr, ImageUInt8Ptr imgPtr, bool isMaxtree, double radius);


    

    NodeCT proxy(NodeId id) const;
    NodeCT getSC(int p) const noexcept;
    void setSC(int p, NodeCT node);
    NodeCT getRoot();
    void setRoot(NodeCT n);

    int getTreeType() const noexcept{return treeType;}
    inline bool isMaxtree()const noexcept{ return treeType == MAX_TREE;}
    inline int getNumNodes()const noexcept{ return numNodes; }
    inline int getNumRowsOfImage()const noexcept{ return numRows;}
    inline int getNumColsOfImage()const noexcept{ return numCols;}
    inline AdjacencyRelationPtr getAdjacencyRelation() noexcept {return adj;}
    inline bool isAncestor(NodeId u, NodeId v) const noexcept {return arena.timePreOrder[u] <= arena.timePreOrder[v] && arena.timePostOrder[u] >= arena.timePostOrder[v];}
    inline bool isDescendant(NodeId u, NodeId v) const noexcept {return arena.timePreOrder[v] <= arena.timePreOrder[u] && arena.timePostOrder[v] >= arena.timePostOrder[u];}
    inline bool isComparable(NodeId u, NodeId v) const noexcept {return isAncestor(u, v) || isAncestor(v, u);}
    inline bool isStrictAncestor(NodeId u, NodeId v) const noexcept {return u != v && arena.timePreOrder[u] <= arena.timePreOrder[v] && arena.timePostOrder[u] >= arena.timePostOrder[v];}
    inline bool isStrictDescendant(NodeId u, NodeId v) const noexcept { return u != v && arena.timePreOrder[v] <= arena.timePreOrder[u] && arena.timePostOrder[v] >= arena.timePostOrder[u];}
    inline bool isStrictComparable(NodeId u, NodeId v) const noexcept { return isStrictAncestor(u, v) || isStrictAncestor(v, u);}
	//inline int getDepth() const noexcept {return depth;}

    inline NodeArena::RepsOfCCRangeById getRepCNPs() const { return arena.getRepsOfCC(root); } //iterador
    inline NodeArena::RepsOfCCRangeById getRepCNPsOfCCById(NodeId id) const { return arena.getRepsOfCC(id); } //iterador
    inline auto getPixelsOfCCById(NodeId id) const{ return pixelBuffer->getPixelsBySet(arena.getRepsOfCC(id));  } //iterador
    inline auto getPixelsOfFlatzone(int repFZ) const{ return pixelBuffer->getPixelsBySet(repFZ); } //iterador
    inline int getLevelById(NodeId id) const noexcept{return arena.threshold2[id];}
    inline void setLevelById(NodeId id, int level){ arena.threshold2[id] = level;}
    inline int getResidueById(NodeId id) const noexcept{  if(arena.threshold2[id] == -1) return arena.threshold2[id]; else return arena.threshold2[id] - arena.threshold2[arena.parentId[id]]; }
    
    void releaseNode(NodeId id) noexcept { arena.releaseNode(id); numNodes--; }
    inline int32_t getAreaById(NodeId id) const noexcept{ return arena.areaCC[id];}
    inline void setAreaById(NodeId id, int32_t area) noexcept { arena.areaCC[id] = area;}
    inline int getTimePostOrderById(NodeId id) const noexcept{ return arena.timePostOrder[id]; }
    inline void setTimePostOrderById(NodeId id, int time) noexcept{ arena.timePostOrder[id] = time; } 
    inline int getTimePreOrderById(NodeId id) const noexcept{ return arena.timePreOrder[id]; }
    inline void setTimePreOrderById(NodeId id, int time) noexcept { arena.timePreOrder[id] = time; }

    inline NodeId getParentById(NodeId id) const noexcept {return arena.parentId[id];}
    inline int getNumChildrenById(NodeId id) const noexcept{ return arena.childCount[id]; }
    inline int getNumSiblingsById(NodeId id) const noexcept{ return arena.parentId[id] == -1 ? 0 : arena.childCount[ arena.parentId[id] ] - 1; }
    inline bool hasChildById(NodeId nodeId, NodeId childId) const { return arena.parentId[childId] == nodeId;}
    inline bool isLeaf(NodeId id) const noexcept{ return arena.childCount[id] == 0; }
    inline bool isMaxtreeNodeById(NodeId id) const noexcept{ return arena.parentId[id] != -1 && arena.threshold2[id] > arena.threshold2[arena.parentId[id]];}
    inline NodeArena::ChildRange getChildrenById(NodeId id) const {return arena.children(id);}    
    inline NodeId getRootById()const noexcept { return this->root;}
    inline void setRootById(NodeId n){ setParentById(n, -1); this->root = n; }
    inline NodeId getSCById(int p) const noexcept{ return this->pixelToNodeId[p];}
    inline void setSCById(int p, NodeId id) noexcept { this->pixelToNodeId[p] = id;}    
    inline auto getCNPsById(NodeId id) const { return pixelBuffer->getPixelsBySet(arena.repNode[id]);}
    inline int getNumCNPsById(NodeId id) const { return this->pixelBuffer->numPixelsInSet(arena.repNode[id]);}


    inline void removeChildById(int parentId, int childId, bool release);
    inline void addChildById(int parentId, int childId);
    inline void spliceChildrenById(int toId, int fromId);
    inline void setParentById(NodeId nodeId, NodeId parentId);

    void computerArea(NodeId node);
    ImageUInt8Ptr reconstructionImage();
    std::vector<NodeId> getLeaves();

    void prunning(NodeId nodeId);
    void mergeWithParent(NodeCT node);
    void mergeWithParent(std::vector<int>& flatzone);
    

    static bool validateStructure(ComponentTreePtr tree)  {
        return validateStructure(tree.get());
    }
    static bool validateStructure(const ComponentTree* tree){
        const auto erroMsg = [&](std::string msg){ 
            std::cerr << "❌ Erro " << msg << "\n";
            return false;
        };
        const auto infoMsg = [&](std::string msg){ 
            std::cout << "✅ " << msg << "\n";
        };

        if (tree->root < 0 || tree->root >= (int)tree->arena.size())
            return erroMsg("root inválido");

        // 1) Exigir exatamente 1 raiz (DESCONSIDERANDO slots liberados pela free-list)
        int roots = 0;
        for (int id = 0; id < (int)tree->arena.size(); ++id) {
            if (!tree->arena.isFree(id) && tree->arena.parentId[id] == -1) {
                ++roots;
            }
        }
        if (roots != 1)
            return erroMsg("1: A árvore NÃO possui exatamente uma raiz; a soma de parentId == -1 (desconsiderando slots livres) é "+ std::to_string(roots) +" mas deveria ser 1");
        else
            infoMsg("A árvore contém exatamente 1 raiz (excluindo slots livres)");
    

        // 2) Pai consistente
        for (int id = 0; id < (int)tree->arena.size(); ++id) {
            if (tree->arena.parentId[id] != -1) {
                if (tree->arena.parentId[id] < 0 || tree->arena.parentId[id] >= (int)tree->arena.size()) 
                    return erroMsg("2: O parentId="+ std::to_string(id)+" está fora do range [0, "+ std::to_string(tree->arena.size()) + "]");
                if (id == tree->arena.parentId[id]) 
                    return erroMsg("3: O parentId="+std::to_string(id)+" está apontando para si mesmo");
            }
        }
        infoMsg("A estrutura de parentesco arena.parentId está consistente");

        // 3) Encadeamento filhos/irmãos + lastChildId + childCount
        for (int u = 0; u < (int)tree->arena.size(); ++u) {
            int cnt = 0, last = -1;
            if (tree->arena.firstChildId[u] == -1) {
                if (tree->arena.lastChildId[u] != -1 || tree->arena.childCount[u] != 0) 
                    return erroMsg("4: Nó sem filhos mas last/childCount incoerentes");
            } else {
                // Caminha pela lista de filhos de u
                for (int c = tree->arena.firstChildId[u]; c != -1; c = tree->arena.nextSiblingId[c]) {
                    if (c < 0 || c >= (int)tree->arena.size()) 
                        return erroMsg("5: Estrutura de filhos e irmãos (firstChildId/nextSiblingId) fora do range [0, "+ std::to_string(tree->arena.size()) + "]");
                    if (tree->arena.parentId[c] != u) 
                        return erroMsg("6: Filho com parentId diferente do pai");
                    
                        // Checa simetria prev/next
                    if (tree->arena.prevSiblingId[c] != last) 
                        return erroMsg("7: Estrutura de irmãos prevSiblingId está inconsistente");
                    if (last != -1 && tree->arena.nextSiblingId[last] != c) 
                        return erroMsg("8: Estrutura de irmãos nextSiblingId está inconsistente");
                    last = c; 
                    ++cnt;
                }
                if (last != tree->arena.lastChildId[u]) 
                    return erroMsg("8: Estrutura de filhos lastChildId não bate com último encadeado");
                if (cnt != tree->arena.childCount[u]) 
                    return erroMsg("9: Estrutura que armazena a quantidade de filhos childCount não bate com encadeamento");
            }
        }
        infoMsg("A estrutura de filhos/irmãos está consistente");
        return true;
    }
    
    static std::vector<NodeId> getNodesThreshold(ComponentTreePtr tree, int areaThreshold){
        std::vector<NodeId> lista;
        FastQueue<NodeId> queue;
        queue.push(tree->root);

        while (!queue.empty()) {
            NodeId id = queue.pop();
            if (tree->arena.areaCC[id] > areaThreshold) {
                for(NodeId c: tree->arena.children(id)){
                    queue.push(c);
                }
            } else {
                lista.push_back(id);
            }
        }
        
        return lista;
    }
    
    // ====================== Iteradores por ID (sem proxy) ====================== //

    // ================== Iterador de NodeIds VÁLIDOS — versão otimizada ================== //
    class InternalIteratorValidNodeIds {
    private:
        const int* rep_;        // ponteiro p/ arena.repNode[0]
        NodeId cur_;            // posição atual
        NodeId end_;            // N = arena.repNode.size()

        // avança cur_ até um id válido ou coloca cur_ = end_
        inline void settle_() noexcept {
            while (cur_ < end_ && rep_[cur_] == -1) ++cur_;
        }

    public:
        using iterator_category = std::input_iterator_tag;
        using value_type        = NodeId;
        using difference_type   = std::ptrdiff_t;
        using pointer           = const NodeId*;
        using reference         = const NodeId&;

        inline InternalIteratorValidNodeIds(ComponentTree* T, NodeId start) noexcept
        : rep_(T ? T->arena.repNode.data() : nullptr),
        cur_(T ? start : 0),
        end_(T ? static_cast<NodeId>(T->arena.repNode.size()) : 0) {
            settle_();
        }

        inline InternalIteratorValidNodeIds& operator++() noexcept {
            ++cur_;
            settle_();
            return *this;
        }

        inline NodeId operator*() const noexcept { return cur_; }

        // iterador input: comparar só posição é suficiente
        inline bool operator==(const InternalIteratorValidNodeIds& other) const noexcept { return cur_ == other.cur_; }
        inline bool operator!=(const InternalIteratorValidNodeIds& other) const noexcept { return cur_ != other.cur_; }
    };

    class IteratorValidNodeIds {
    private:
        ComponentTree* T_ = nullptr;
    public:
        inline explicit IteratorValidNodeIds(ComponentTree* T) noexcept : T_(T) {}

        inline InternalIteratorValidNodeIds begin() const noexcept { return InternalIteratorValidNodeIds(T_, 0); }
        // sentinela: construtor com T=nullptr dá end_=0, cur_=0 → iguala quando begin() também chega ao fim
        inline InternalIteratorValidNodeIds end()   const noexcept { return InternalIteratorValidNodeIds(nullptr, 0); }
    };

    /** Range para iterar NodeId válidos (exclui slots livres) */
    inline IteratorValidNodeIds getNodeIds() noexcept { return IteratorValidNodeIds(this); }
    inline IteratorValidNodeIds getNodeIds() const noexcept { return IteratorValidNodeIds(const_cast<ComponentTree*>(this)); }


    // Pós-ordem (retorna NodeId)
    class InternalIteratorPostOrderTraversalId {
    private:
        struct Item { int id; bool expanded; };

        ComponentTree* T_ = nullptr;
        FastStack<Item> st_;
        NodeId current_ = -1;

        // Avança até o próximo nó a ser emitido (ou deixa current_ = -1 se acabou)
        void settle_() noexcept {
            while (!st_.empty()) {
                Item &top = st_.top();
                if (!top.expanded) {
                    top.expanded = true;

                    // A ordem resultante não é garantida (costuma ser direita->esquerda).
                    for (int c = T_->arena.firstChildId[top.id]; c != -1; c = T_->arena.nextSiblingId[c]) {
                        st_.push(Item{c, false});
                    }
                    // volta ao loop: agora o topo será algum filho
                } else {
                    current_ = top.id;      // todos os filhos já emitidos
                    return;
                }
            }
            current_ = -1; // fim
        }

    public:
        using iterator_category = std::input_iterator_tag;
        using value_type        = NodeId;
        using difference_type   = std::ptrdiff_t;
        using pointer           = const NodeId*;
        using reference         = const NodeId&;

        InternalIteratorPostOrderTraversalId(ComponentTree* T, NodeId rootId) noexcept : T_(T) {
            if (T_ && rootId >= 0) {
                st_.push(Item{rootId, false});
                settle_(); // posiciona no primeiro elemento
            } else {
                current_ = -1;
            }
        }

        // pré-incremento
        InternalIteratorPostOrderTraversalId& operator++() noexcept {
            if (!st_.empty()) st_.pop();  // consome o atual
            settle_();                    // acha o próximo
            return *this;
        }

        // desreferência → NodeId atual
        NodeId operator*() const noexcept { return current_; }

        bool operator==(const InternalIteratorPostOrderTraversalId& other) const noexcept {
            return current_ == other.current_;
        }
        bool operator!=(const InternalIteratorPostOrderTraversalId& other) const noexcept {
            return !(*this == other);
        }
    };

    class IteratorPostOrderTraversalId {
    private:
        ComponentTree* T_ = nullptr;
        int rootId_ = -1;
    public:
        explicit IteratorPostOrderTraversalId(ComponentTree* T, int rootId) noexcept : T_(T), rootId_(rootId) {}

        InternalIteratorPostOrderTraversalId begin() const noexcept {
            return InternalIteratorPostOrderTraversalId(T_, rootId_);
        }
        InternalIteratorPostOrderTraversalId end() const noexcept {
            return InternalIteratorPostOrderTraversalId(nullptr, -1);
        }
    };

    auto getIteratorPostOrderTraversalById(int id) {        
        return IteratorPostOrderTraversalId(this, id);
    }
    auto getIteratorPostOrderTraversalById() {        
        return IteratorPostOrderTraversalId(this, root);
    }

    // Largura (BFS) — retorna NodeId
    class InternalIteratorBreadthFirstTraversalId {
    private:
        ComponentTree* T_ = nullptr;
        FastQueue<int> q_;

    public:
        using iterator_category = std::input_iterator_tag;
        using value_type        = int;
        using difference_type   = std::ptrdiff_t;
        using pointer           = const int*;
        using reference         = const int&;

        InternalIteratorBreadthFirstTraversalId(ComponentTree* T, int rootId) noexcept : T_(T) {
            if (T_ && rootId >= 0) q_.push(rootId);
        }

        InternalIteratorBreadthFirstTraversalId& operator++() noexcept {
            if (!q_.empty()) {
                int u = q_.pop();
                for (int c: T_->arena.children(u)) {
                    q_.push(c);
                }
            }
            return *this;
        }

        int operator*() const noexcept { return q_.front(); }

        bool operator==(const InternalIteratorBreadthFirstTraversalId& other) const noexcept {
            return q_.empty() == other.q_.empty();
        }
        bool operator!=(const InternalIteratorBreadthFirstTraversalId& other) const noexcept {
            return !(*this == other);
        }
    };

    class IteratorBreadthFirstTraversalId {
    private:
        ComponentTree* T_ = nullptr;
        int rootId_ = -1;
    public:
        explicit IteratorBreadthFirstTraversalId(ComponentTree* T, int rootId) noexcept : T_(T), rootId_(rootId) {}

        InternalIteratorBreadthFirstTraversalId begin() const noexcept {
            return InternalIteratorBreadthFirstTraversalId(T_, rootId_);
        }
        InternalIteratorBreadthFirstTraversalId end() const noexcept {
            return InternalIteratorBreadthFirstTraversalId(nullptr, -1);
        }
    };

    IteratorBreadthFirstTraversalId getIteratorBreadthFirstTraversalById(NodeId id) noexcept {
        return IteratorBreadthFirstTraversalId(this, id);
    }
    IteratorBreadthFirstTraversalId getIteratorBreadthFirstTraversalById() noexcept {
        return IteratorBreadthFirstTraversalId(this, root);
    }
    // ================== Fim dos iteradores por ID (sem proxy) ================== //

    // ================== Iterador para caminho até a raiz por NodeId (sem proxy) ================== //
    class InternalIteratorNodesOfPathToRootId {
    private:
        ComponentTree* T_ = nullptr;
        NodeId currentId_ = -1;

    public:
        using iterator_category = std::input_iterator_tag;
        using value_type        = NodeId;
        using difference_type   = std::ptrdiff_t;
        using pointer           = const NodeId*;
        using reference         = const NodeId&;

        InternalIteratorNodesOfPathToRootId(ComponentTree* T, NodeId startId) noexcept : T_(T), currentId_(startId) {}

        InternalIteratorNodesOfPathToRootId& operator++() noexcept {
            if (T_ && currentId_ != -1) {
                currentId_ = T_->arena.parentId[currentId_];
            }
            return *this;
        }

        NodeId operator*() const noexcept { return currentId_; }

        bool operator==(const InternalIteratorNodesOfPathToRootId& other) const noexcept {
            return currentId_ == other.currentId_;
        }
        bool operator!=(const InternalIteratorNodesOfPathToRootId& other) const noexcept {
            return !(*this == other);
        }
    };

    class IteratorNodesOfPathToRootId {
    private:
        ComponentTree* T_ = nullptr;
        NodeId startId_ = -1;

    public:
        IteratorNodesOfPathToRootId(ComponentTree* T, NodeId startId) noexcept : T_(T), startId_(startId) {}

        InternalIteratorNodesOfPathToRootId begin() const noexcept {
            return InternalIteratorNodesOfPathToRootId(T_, startId_);
        }
        InternalIteratorNodesOfPathToRootId end() const noexcept {
            return InternalIteratorNodesOfPathToRootId(nullptr, -1);
        }
    };

    IteratorNodesOfPathToRootId getNodesOfPathToRootById(NodeId id) noexcept {
        return IteratorNodesOfPathToRootId(this, id);
    }


    // ================== Iterador de descendentes por ID (sem proxy) ================== //
    /**
     * @brief Iterador em largura (BFS) sobre os descendentes de um nó.
     *
     * Percorre todos os nós descendentes do nó informado, excluindo o
     * próprio nó raiz do percurso. A ordem é em largura (BFS).
     * Uso típico:
     *   for (NodeId u : tree->getNodesDescendantsById(id)) { ... }
     */
    class InternalIteratorNodesDescendantsId {
    private:
        ComponentTree* T_ = nullptr;
        FastQueue<int> q_;

    public:
        using iterator_category = std::input_iterator_tag;
        using value_type        = NodeId;
        using difference_type   = std::ptrdiff_t;
        using pointer           = const NodeId*;
        using reference         = const NodeId&;

        explicit InternalIteratorNodesDescendantsId(ComponentTree* T, NodeId rootId) noexcept : T_(T) {
            if (T_ && rootId >= 0) {
                // Descendentes EXCLUEM o próprio nó: inicia com os filhos diretos
                for (int c : T_->arena.children(rootId)) q_.push(c);
            }
        }

        InternalIteratorNodesDescendantsId& operator++() noexcept {
            if (!q_.empty()) {
                int u = q_.pop();
                for (int v : T_->arena.children(u)) q_.push(v);
            }
            return *this;
        }

        NodeId operator*() const noexcept { return q_.front(); }

        bool operator==(const InternalIteratorNodesDescendantsId& other) const noexcept {
            return q_.empty() == other.q_.empty();
        }
        bool operator!=(const InternalIteratorNodesDescendantsId& other) const noexcept {
            return !(*this == other);
        }
    };

    /**
     * @brief Range leve para iterar os descendentes de `id` via range-for.
     */
    class IteratorNodesDescendantsId {
    private:
        ComponentTree* T_ = nullptr;
        NodeId rootId_ = -1;

    public:
        explicit IteratorNodesDescendantsId(ComponentTree* T, NodeId rootId) noexcept : T_(T), rootId_(rootId) {}

        InternalIteratorNodesDescendantsId begin() const noexcept { return InternalIteratorNodesDescendantsId(T_, rootId_); }
        InternalIteratorNodesDescendantsId end()   const noexcept { return InternalIteratorNodesDescendantsId(nullptr, -1); }
    };

    /**
     * @brief Retorna um range para iterar os descendentes (exclui o próprio nó).
     * @param id Nó base do qual se deseja percorrer os descendentes.
     * @return Range compatível com range-for de `NodeId`.
     */
    IteratorNodesDescendantsId getNodesDescendantsById(NodeId id) noexcept {
        return IteratorNodesDescendantsId(this, id);
    }

};


/**
 * Método Euler Tour + RMQ
  
 Etapa 1: Euler Tour
    Realiza um DFS na árvore e registra:
	1.	A ordem dos nós visitados → euler[]
	2.	A profundidade de cada nó na árvore durante o percurso → depth[]
	3.	A índice da primeira ocorrência de cada nó no vetor euler → firstOccurrence[]

  Etapa 2: RMQ na profundidade
    Para responder LCA(u, v):
	1.	Pegue i = firstOccurrence[u], 
              j = firstOccurrence[v]
	2.	Realize um RMQ (Range Minimum Query) sobre o vetor depth[] entre as posições min(i, j) e max(i, j) no vetor euler[].
	3.	O resultado do RMQ será o índice do nó com menor profundidade entre u e v no caminho — ou seja, o LCA!

    Exemplo:
      0
     / \
    1   2
   /
  3
  Índices:         0  1  2  3  4  5  6
  euler =         [0, 1, 3, 1, 0, 2, 0]
  depth =         [0, 1, 2, 1, 0, 1, 0]
  firstOccurrence=[0, 1, 5, 2         ]
    
  LCA(3, 2) = 0
    i = firstOccurrence[3] = 2
    j = firstOccurrence[2] = 5
    RMQ: 
      1. Descobrir o intervalo no vetor depth: depth[2..5] = [2, 1, 0, 1]
      2. Encontrar a posição do menor valor: O mínimo é 0, que ocorre em depth[4]
      3. O correspondente em no vetor euler: euler[4] = 0 que é o indice do LCA
	
 */
class LCAEulerRMQ {
private:
    std::vector<int> euler;            // timePreOrder dos nós na ordem de visita
    std::vector<int> depth;            // profundidade associada a cada posição em euler
    std::vector<int> firstOccurrence;  // [timePreOrder] = posição no vetor euler
    std::vector<std::vector<int>> st;  // Sparse Table para RMQ
    ComponentTree* tree;
public:
    LCAEulerRMQ(ComponentTree* tree): tree(tree) {
        
        int n = tree->getNumNodes();
        euler.reserve(n);
        depth.reserve(n);
        firstOccurrence.resize(n, -1);

        depthFirstTraversal(tree->getRootById(), 0);
        buildSparseTable();
    }

	NodeId findLowestCommonAncestor(NodeId u, NodeId v) {
        int i = firstOccurrence[u];
        int j = firstOccurrence[v];
        if (i > j) std::swap(i, j);
        int idx = rmq(i, j);
        return euler[idx];
    }


private:
    void depthFirstTraversal(NodeId timeNode, int d) {
        if (firstOccurrence[timeNode] == -1)
            firstOccurrence[timeNode] = euler.size();

        euler.push_back(timeNode);
        depth.push_back(d);

        for (NodeId child : tree->getChildrenById(timeNode)) {
            depthFirstTraversal(child, d + 1);
            euler.push_back(timeNode);
            depth.push_back(d);
        }
    }

    void buildSparseTable() {
        int n = depth.size();
        int logn = std::log2(n) + 1;
        st.assign(n, std::vector<int>(logn));

        for (int i = 0; i < n; ++i)
            st[i][0] = i;

        for (int j = 1; (1 << j) <= n; ++j) {
            for (int i = 0; i + (1 << j) <= n; ++i) {
                int l = st[i][j - 1];
                int r = st[i + (1 << (j - 1))][j - 1];
                st[i][j] = (depth[l] < depth[r]) ? l : r;
            }
        }
    }

    int rmq(int l, int r) {
        int len = r - l + 1;
        int k = std::log2(len);
        int a = st[l][k];
        int b = st[r - (1 << k) + 1][k];
        return (depth[a] < depth[b]) ? a : b;
    }
};


#include "../include/ComponentTree.tpp"

#endif
