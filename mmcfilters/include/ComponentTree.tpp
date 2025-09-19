#include <list>
#include <array>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <iomanip>
#include <utility>
#include <algorithm> 

#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"
#include "../include/AdjacencyRelation.hpp"
#include "../include/BuilderMorphologicalTreeByUnionFind.hpp"


//ComponentTree(ImageUInt8Ptr img, bool isMaxtree, double radius = 1.5);
//ComponentTree(ImageUInt8Ptr img, std::string ToSInperpolation="self-dual");
    

ComponentTree::ComponentTree(ImageUInt8Ptr img, bool isMaxtree, double radius) : 
    numRows(img->getNumRows()), numCols(img->getNumCols()), treeType( (isMaxtree? ComponentTree::MAX_TREE : ComponentTree::MIN_TREE) ), 
    adj(std::make_shared<AdjacencyRelation>(numRows, numCols, radius)), numNodes(0){   

        this->pixelToNodeId.resize(numRows * numCols, -1);
        BuilderComponentTreeByUnionFind builderUF(adj.get(), isMaxtree);
        build(img, builderUF);
}

ComponentTree::ComponentTree(ImageUInt8Ptr img, std::string ToSInperpolation):
    numRows(img->getNumRows()), numCols(img->getNumCols()), treeType(ComponentTree::TREE_OF_SHAPES), 
    adj(nullptr), numNodes(0){   
        
        this->pixelToNodeId.resize(numRows * numCols, -1);
        BuilderTreeOfShapeByUnionFind builderUF( ToSInperpolation == "4c8c" );
        build(img, builderUF);
}


void ComponentTree::build(const ImageUInt8Ptr& imgPtr, IMorphologicalTreeBuilder& builderUF){ 
    
    auto [parent, orderedPixels, numNodes] = builderUF.createTreeByUnionFind(imgPtr);

    int numPixels = imgPtr->getSize();
    auto img = imgPtr->rawData();

    this->reserveNodes(numNodes);
    this->pixelBuffer = std::make_shared<PixelSetManager>(numPixels, numNodes);
    this->pixelView = this->pixelBuffer->view();
    int indice = 0;
    for (int i = 0; i < numPixels; i++) {
        int p = orderedPixels[i];

        //Construção da árvore e arena
        if (p == parent[p]) {
            pixelToNodeId[p] = this->root = this->makeNode(p, -1, img[p]);
        } else if (img[p] != img[parent[p]]) {
            pixelToNodeId[p] = this->makeNode(p, pixelToNodeId[parent[p]], img[p]);
        } else {
            pixelToNodeId[p] = pixelToNodeId[parent[p]];
        }

        //Construção de PixelSetManager
        if (p == parent[p] || img[p] != img[parent[p]]) {
            pixelView.indexToPixel[indice] = p;
            pixelView.pixelToIndex[p] = indice;
            pixelView.sizeSets[indice] = 1;
            pixelView.pixelsNext[p] = p;
            indice++;
        } else {
            pixelView.pixelsNext[p] = pixelView.pixelsNext[parent[p]];
            pixelView.pixelsNext[parent[p]] = p;
            int idx = pixelView.pixelToIndex[parent[p]];
            pixelView.sizeSets[idx]++;
        }
    }


	computerTreeAttributes();
}


template <typename PixelType>
ComponentTreePtr ComponentTree::createFromAttributeMapping(ImagePtr<PixelType> attrMappingPtr, ImageUInt8Ptr imgPtr, bool isMaxtree, double radius) {
    AdjacencyRelationPtr adj = std::make_shared<AdjacencyRelation>(imgPtr->getNumRows(), imgPtr->getNumCols(), radius);	
    ComponentTreePtr tree = ComponentTree::create(imgPtr->getNumRows(), imgPtr->getNumCols(), isMaxtree, adj);
    BuilderComponentTreeByUnionFind builderUF(adj.get(), isMaxtree);
    auto [parent, orderedPixels, numNodes] = builderUF.createTreeByUnionFind(attrMappingPtr);

    int numPixels = imgPtr->getSize();
    auto img = imgPtr->rawData();
    auto attrMapping = attrMappingPtr->rawData();

    tree->reserveNodes(numNodes);
    tree->pixelBuffer = std::make_shared<PixelSetManager>(numPixels, numNodes);
    tree->pixelView = tree->pixelBuffer->view();
    int indice = 0;
    float epsilon = 1e-5f; // Tolerância para comparação de pixels flutuantes
    
    auto sameLevel = [&](PixelType a, PixelType b){ 
        if constexpr (std::is_floating_point_v<PixelType>)
            return std::fabs(attrMapping[a] - attrMapping[parent[b]]) < epsilon; 
        else
            return attrMapping[a] == attrMapping[b]; 
    };

    for (int i = 0; i < numPixels; i++) {
        int p = orderedPixels[i];
        //Construção da árvore e arena
        if (p == parent[p]) {
            tree->pixelToNodeId[p] = tree->root = tree->makeNode(p, -1, img[p]);
        } 
        else if (!sameLevel(p, parent[p])) {
            tree->pixelToNodeId[p] = tree->makeNode(p, tree->pixelToNodeId[parent[p]], img[p]);
        } 
        else {
            tree->pixelToNodeId[p] = tree->pixelToNodeId[parent[p]];
        }
        
        //Construção de PixelSetManager
        if (p == parent[p] || !sameLevel(p, parent[p])) {
            tree->pixelView.indexToPixel[indice] = p;
            tree->pixelView.pixelToIndex[p] = indice;
            tree->pixelView.sizeSets[indice] = 1;
            tree->pixelView.pixelsNext[p] = p;
            indice++;
        } else {
            tree->pixelView.pixelsNext[p] = tree->pixelView.pixelsNext[parent[p]];
            tree->pixelView.pixelsNext[parent[p]] = p;
            int idx = tree->pixelView.pixelToIndex[parent[p]];
            tree->pixelView.sizeSets[idx]++;
        }
    }
    tree->computerTreeAttributes();

    //ajustar o level de cada nó
    for(NodeId id: tree->getNodeIds()){
        int media = 0;
        for(int p: tree->getCNPsById(id)){
            media += img[p];
        }
        media /=  tree->getNumCNPsById(id); 
        tree->setLevelById(id, media);
    }
    return tree;
}



/*inline NodeCT ComponentTree::proxy(NodeId id) const {
    return NodeCT(const_cast<ComponentTree*>(this), id);
}*/

inline NodeCT ComponentTree::proxy(NodeId id) const {
    if (id < 0 || id >= static_cast<NodeId>(arena.size())) {
        throw std::out_of_range("Node ID is out of range in arena!");
    }
    if (arena.isFree(id)) { 
        throw std::runtime_error("Node ID refers to a freed node in arena!");
    }
    return NodeCT(const_cast<ComponentTree*>(this), id);
}

inline NodeCT ComponentTree::getSC(int p) const noexcept { return proxy(this->pixelToNodeId[p]); }
inline void ComponentTree::setSC(int p, NodeCT node){ setSCById(p, node); }
inline NodeCT ComponentTree::getRoot() { return proxy(this->root); }
inline void ComponentTree::setRoot(NodeCT n){ setRootById(n); }


NodeId ComponentTree::makeNode(int repNode, NodeId parentId, int threshold2){
    // Aloca ID contíguo
    NodeId id = this->arena.allocate(repNode, threshold2);

    // Encadeia no pai (por ID) se houver
    if (parentId >= 0) {
        addChildById(parentId, id);
    }
    else{
        std::cout<<"Root:" << repNode << std::endl;
    }
    //contador de nós
    this->numNodes++;    
    return id;
}

inline void ComponentTree::setParentById(NodeId nodeId, NodeId parentId) {
    if (parentId == arena.parentId[nodeId]) return;
    if (parentId == -1) {
        if (arena.parentId[nodeId] != -1) 
            removeChildById(arena.parentId[nodeId], parentId, false);
    } else {
        addChildById(parentId, nodeId);
    }        
}

void ComponentTree::addChildById(int parentId, int childId) {
    if (parentId < 0 || childId < 0) return;

    // Se o filho já tem pai, desconecta antes de ler P/C
    if (arena.parentId[childId] != -1) {
        removeChildById(arena.parentId[childId], childId, false);
    }

    arena.parentId[childId]      = parentId;
    arena.prevSiblingId[childId] = arena.lastChildId[parentId];
    arena.nextSiblingId[childId] = -1;

    if (arena.firstChildId[parentId] == -1) {
        arena.firstChildId[parentId] = arena.lastChildId[parentId] = childId;
    } else {
        arena.nextSiblingId[arena.lastChildId[parentId]] = childId;
        arena.lastChildId[parentId] = childId;
    }
    ++arena.childCount[parentId];

    
    //assert(ComponentTree::validateStructure(this) && "ComponentTree topology invariant failed after addChildById");
    
}

// Remove um filho 'childId' da lista encadeada de filhos do pai 'parentId'.
inline void ComponentTree::removeChildById(int parentId, int childId, bool release) {
    if (parentId < 0 || childId < 0) return;
    if (arena.parentId[childId] != parentId) return;

    const int prev = arena.prevSiblingId[childId];
    const int next = arena.nextSiblingId[childId];

    if (prev == -1) arena.firstChildId[parentId] = next;
    else            arena.nextSiblingId[prev] = next;

    if (next == -1) arena.lastChildId[parentId] = prev;
    else            arena.prevSiblingId[next] = prev;

    if (arena.childCount[parentId] > 0) 
        --arena.childCount[parentId];

    arena.parentId[childId] = -1;
    arena.prevSiblingId[childId] = -1;
    arena.nextSiblingId[childId] = -1;
    if(release){
        releaseNode(childId);
        //assert(ComponentTree::validateStructure(this) && "ComponentTree topology invariant failed after removeChildById");
    }
}


// Move todos os filhos de 'fromId' para o fim da lista de filhos de 'toId'.
inline void ComponentTree::spliceChildrenById(int toId, int fromId) {
    if (toId < 0 || fromId < 0 || toId == fromId) return;

    NodeId firstFrom = arena.firstChildId[fromId];
    if (firstFrom == -1) return; // nada para mover

    // 1) todos os filhos de 'fromId' passam a ter pai 'toId'
    for (int c = arena.firstChildId[fromId]; c != -1; c = arena.nextSiblingId[c]) {
        arena.parentId[c] = toId;
    }

    // 2) concatena a lista de filhos de 'fromId' no fim da lista de 'toId'
    if (arena.firstChildId[toId] == -1) {
        // 'toId' não tinha filhos — vira exatamente a lista de 'fromId'
        arena.firstChildId[toId] = arena.firstChildId[fromId];
        arena.lastChildId[toId]  = arena.lastChildId[fromId];
        // o primeiro filho já tem prevSiblingId == -1 porque era o primeiro de 'fromId'
    } else {
        // 'toId' já tinha filhos — encadeia no final
        arena.nextSiblingId[ arena.lastChildId[toId] ] = arena.firstChildId[fromId];
        arena.prevSiblingId[ arena.firstChildId[fromId] ] = arena.lastChildId[toId];
        arena.lastChildId[toId] = arena.lastChildId[fromId];
    }

    // 3) atualiza contadores
    arena.childCount[toId] += arena.childCount[fromId];

    // 4) zera a lista de 'fromId'
    arena.firstChildId[fromId] = -1;
    arena.lastChildId[fromId]  = -1;
    arena.childCount[fromId]   = 0;

   // assert(ComponentTree::validateStructure(this) && "ComponentTree topology invariant failed after spliceChildrenById");
}




void ComponentTree::computerTreeAttributes(){

	int timer = 0;
	//int maxDepth = 0;
	//std::vector<int> depth(this->numNodes);
	computerIncrementalAttributes(this->root,
		[&](NodeId nodeId) -> void { //pre-processing
            arena.areaCC[nodeId] = getNumCNPsById(nodeId);;
            arena.timePreOrder[nodeId] = timer++;
			//depth[nodeId] = arena.parentId[nodeId] == -1 ? 0 : depth[arena.parentId[nodeId]] + 1;
		},
		[&](NodeId parentId, NodeId childId) -> void { //merge-processing
            arena.areaCC[parentId] += arena.areaCC[childId];

		},
		[&](NodeId nodeId) -> void { // post-processing
			arena.timePostOrder[nodeId] = timer++;
            //maxDepth = std::max(maxDepth, depth[nodeId]);
		}
	);
	//this->depth = maxDepth;
}

/*
template <typename CNPsType>
void ComponentTree<CNPsType>::computerArea(NodeId id){
    int32_t area = getNumCNPsById(id);
    for (NodeId c : arena.children(id)) {
        computerArea(c);
        area += arena.areaCC[c]; // já calculado na recursão
    }
    arena.areaCC[id] = area;
}*/

void ComponentTree::prunning(NodeId nodeId){
    assert(nodeId && "node is invalid");
    assert(getParentById(nodeId) != -1 && "node is root");

    const int parentId = arena.parentId[nodeId]; 
    if(parentId >= 0){ 
        // 1) desconecta 'node' do pai
        removeChildById(parentId, nodeId, false);

        // 2) BFS ID-first na subárvore para redirecionar UF/SC e contabilizar remoções
        FastQueue<NodeId> q;
        q.push(nodeId);

        const int parentRep = arena.repNode[parentId];
        while (!q.empty()) {
            NodeId curId = q.pop();

            // enfileira filhos por IDs
            for (NodeId c : arena.children(curId)) {
                q.push(c);
            }

            // une representantes no UF e atualiza pixel->node para o representante
            int repCur = arena.repNode[curId];
            setSCById(repCur, parentId);
            pixelBuffer->mergeSetsByRep(parentRep, repCur);

            // release no node removido
            releaseNode(curId);
        }
    }
}



void ComponentTree::mergeWithParent(NodeCT node)
{
    if (!node || !node.getParent()) return;

    const int nodeId   = node.getIndex();
    const int parentId = arena.parentId[nodeId];

    // 1) tira 'node' da lista de filhos do pai
    removeChildById(parentId, nodeId, false);

    // 2) move os filhos de 'node' para o pai (preserva ordem) — tudo por IDs
    spliceChildrenById(parentId, nodeId);

    // 3) une representantes e atualiza pixel->node para o representante do nó colapsado
    const int parentRep = arena.repNode[parentId];
    const int nodeRep   = arena.repNode[nodeId];
    setSCById(nodeRep, parentId);
    pixelBuffer->mergeSetsByRep(parentRep, nodeRep);

    // 4) marca o nó como desconectado
    releaseNode(nodeId);
}


void ComponentTree::mergeWithParent(std::vector<int>& flatzone){
    int idFlatzone = flatzone.front();
    NodeCT node = proxy(this->pixelToNodeId[idFlatzone]);
    if(getNumCNPsById(node) == static_cast<int>(flatzone.size())) {
        this->mergeWithParent(node);
    }
    else{
        //TODO: pensar em como otimizar esse caso
        //winner ganha os pixels de flatzone e o loser perde esses pixels
/*
        NodeP parent = node.getParent();
        int repWinner = parent.getRepNode(); //representante do pai
        int repLoser  = node.getRepNode();   //representante do filho

        //1. Recupera índices dos representantes
        int idxRootWinner = pixelView.pixelToIndex[repWinner];
        int idxRootLoser  = pixelView.pixelToIndex[repLoser];

        //2. Atualiza a quantidade de cnps
        pixelView.sizeSets[idxRootWinner] += flatzone->size();
        pixelView.sizeSets[idxRootLoser] -= flatzone->size();

        for( int p: *flatzone) {				
            parent.addRepCNPs(p);
            this->pixelToNode[p] = parent;	


            // 3. Splice O(1) das listas circulares (pixels)
            int nextWinner = pixelView.pixelsNext[repWinner];
            int nextLoser  = pixelView.pixelsNext[repLoser];
            pixelView.pixelsNext[repWinner] = nextLoser;
            pixelView.pixelsNext[repLoser]  = nextWinner;

            // 4. Invalida slot perdedor
            pixelView.sizeSets[idxRootLoser]  = 0;
            pixelView.indexToPixel[idxRootLoser] = -1;

            // 5. Redireciona lookups pelo antigo rep pixel
            pixelView.pixelToIndex[repLoser] = idxRootWinner;
        }
        */
    }
}




std::vector<NodeId> ComponentTree::getLeaves(){
    std::vector<NodeId> leaves;
    FastQueue<NodeId> s;
    s.push(this->root);

    while (!s.empty()) {
        NodeId id = s.pop();
        if (arena.childCount[id] == 0) {
            leaves.push_back(id);
        } else {
            for(NodeId c: arena.children(id)){
                s.push(c);
            }
        }
    }
    return leaves;
}




ImageUInt8Ptr ComponentTree::reconstructionImage(){
    ImageUInt8Ptr imgPtr = ImageUInt8::create(this->numRows, this->numCols);
    this->reconstruction(this->root, imgPtr->rawData());
    return imgPtr;
}


// Pixels
inline void ComponentTree::reconstruction(NodeId id, uint8_t* imgOut) {
    assert(imgOut && "Erro: Ponteiro de saída da imagem é nulo!");
    for (int p : pixelBuffer->getPixelsBySet(arena.repNode[id])) {
        imgOut[p] = static_cast<uint8_t>(arena.threshold2[id]);
    }
    for (int c : arena.children(id)) {
        reconstruction(c, imgOut);
    }
}


