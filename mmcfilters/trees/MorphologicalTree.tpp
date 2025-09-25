#include "../trees/MorphologicalTree.hpp"
#include "../trees/BuilderMorphologicalTreeByUnionFind.hpp"

namespace mmcfilters {

template <typename PixelType>
MorphologicalTreePtr MorphologicalTree::createFromAttributeMapping(ImagePtr<PixelType> attrMappingPtr, ImageUInt8Ptr imgPtr, bool isMaxtree, double radius) {
    AdjacencyRelationPtr adj = std::make_shared<AdjacencyRelation>(imgPtr->getNumRows(), imgPtr->getNumCols(), radius);
    MorphologicalTreePtr tree = MorphologicalTree::create(imgPtr->getNumRows(), imgPtr->getNumCols(), isMaxtree, adj);
    return createFromAttributeMapping(tree, attrMappingPtr, imgPtr, isMaxtree, radius);
}


template <typename PixelType>
MorphologicalTreePtr MorphologicalTree::createFromAttributeMapping(MorphologicalTreePtr tree, ImagePtr<PixelType> attrMappingPtr, ImageUInt8Ptr imgPtr, bool isMaxtree, double radius) {
    AdjacencyRelationPtr adj = std::make_shared<AdjacencyRelation>(imgPtr->getNumRows(), imgPtr->getNumCols(), radius);
    
    BuilderComponentTree builderUF(adj.get(), isMaxtree);
    auto [parent, orderedPixels, numNodes] = builderUF.createTreeByUnionFind(attrMappingPtr);

    int numPixels = imgPtr->getSize();
    auto img = imgPtr->rawData();
    auto attrMapping = attrMappingPtr->rawData();

    tree->reserveNodes(numNodes);
    tree->pixelBuffer = std::make_shared<PixelSetManager>(numPixels, numNodes);
    tree->pixelView = tree->pixelBuffer->view();
    int indice = 0;
    float epsilon = 1e-5f; // Tolerância para comparação de pixels flutuantes
    
    auto sameLevel = [&](int pIdx, int qIdx){
        if constexpr (std::is_floating_point_v<PixelType>) {
            return std::fabs(static_cast<double>(attrMapping[pIdx]) - static_cast<double>(attrMapping[qIdx])) < epsilon;
        } else {
            return attrMapping[pIdx] == attrMapping[qIdx];
        }
    };

    for (int i = 0; i < numPixels; i++) {
        int p = orderedPixels[i];
        //Construção da árvore e arena
        if (p == parent[p]) {
            tree->pixelToNodeId[p] = tree->root = tree->makeNode(p, InvalidNode, img[p]);
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

} // namespace mmcfilters

