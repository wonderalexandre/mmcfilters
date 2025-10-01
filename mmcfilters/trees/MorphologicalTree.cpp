

#include "../trees/MorphologicalTree.hpp"
#include "../trees/NodeMT.hpp"
#include "../trees/BuilderMorphologicalTreeByUnionFind.hpp"

namespace mmcfilters {


MorphologicalTree::MorphologicalTree(ImageUInt8Ptr img, bool isMaxtree, double radius) : 
    numRows(img->getNumRows()), numCols(img->getNumCols()), treeType( (isMaxtree? MorphologicalTree::MAX_TREE : MorphologicalTree::MIN_TREE) ), 
    adj(std::make_shared<AdjacencyRelation>(numRows, numCols, radius)), numNodes(0){   

        this->pixelToNodeId.resize(numRows * numCols, InvalidNode);
        BuilderComponentTree builderUF(adj.get(), isMaxtree);
        build(img, builderUF);
}

MorphologicalTree::MorphologicalTree(ImageUInt8Ptr img, std::string ToSInperpolation):
    numRows(img->getNumRows()), numCols(img->getNumCols()), treeType(MorphologicalTree::TREE_OF_SHAPES), 
    adj(nullptr), numNodes(0){   
        
        this->pixelToNodeId.resize(numRows * numCols, InvalidNode);
        BuilderTreeOfShape builderUF( ToSInperpolation == "4c8c" );
        build(img, builderUF);
}

NodeMT MorphologicalTree::getSC(int p) const { return proxy(this->pixelToNodeId[p]); }
NodeMT MorphologicalTree::getRoot() const { return proxy(this->root); }
void MorphologicalTree::setSC(int p, NodeMT node){ setSCById(p, node); }
void MorphologicalTree::setRoot(NodeMT n){ setRootById(n); }

void MorphologicalTree::build(const ImageUInt8Ptr& imgPtr, IMorphologicalTreeBuilder& builderUF){ 
    
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
            pixelToNodeId[p] = this->root = this->makeNode(p, InvalidNode, img[p]);
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


NodeMT MorphologicalTree::proxy(NodeId id) const {
    if (id < 0 || id >= static_cast<NodeId>(arena.size())) {
        throw std::out_of_range("Node ID is out of range in arena!");
    }
    if (arena.isFree(id)) { 
        throw std::runtime_error("Node ID refers to a freed node in arena!");
    }
    return NodeMT(const_cast<MorphologicalTree*>(this), id);
}

NodeId MorphologicalTree::makeNode(int repNode, NodeId parentId, int threshold2){
    // Aloca ID contíguo
    NodeId id = this->arena.allocate(repNode, threshold2);

    // Encadeia no pai (por ID) se houver
    if (parentId >= 0) {
        addChildById(parentId, id);
    }
    //contador de nós
    this->numNodes++;    
    return id;
}


void MorphologicalTree::computerTreeAttributes(){

	int timer = 0;
	//int maxDepth = 0;
	//std::vector<int> depth(this->numNodes);
	computerIncrementalAttributes(this, this->root,
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

void MorphologicalTree::prunning(NodeId nodeId){
    assert(nodeId && "node is invalid");
    assert(getParentById(nodeId) != InvalidNode && "node is root");

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



void MorphologicalTree::mergeWithParent(NodeId nodeId)
{
    if (getParentById(nodeId) != InvalidNode) return;

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


void MorphologicalTree::mergeWithParent(std::vector<int>& flatzone){
    int idFlatzone = flatzone.front();
    NodeMT node = proxy(this->pixelToNodeId[idFlatzone]);
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




std::vector<NodeId> MorphologicalTree::getLeaves(){
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

ImageUInt8Ptr MorphologicalTree::reconstructionImage(){
    ImageUInt8Ptr imgPtr = ImageUInt8::create(this->numRows, this->numCols);
    this->reconstruction(this->root, imgPtr->rawData());
    return imgPtr;
}


// Pixels
inline void MorphologicalTree::reconstruction(NodeId id, uint8_t* imgOut) {
    assert(imgOut && "Erro: Ponteiro de saída da imagem é nulo!");
    for (int p : pixelBuffer->getPixelsBySet(arena.repNode[id])) {
        imgOut[p] = static_cast<uint8_t>(arena.threshold2[id]);
    }
    for (int c : arena.children(id)) {
        reconstruction(c, imgOut);
    }
}

} // namespace mmcfilters
