#pragma once

#include "../trees/MorphologicalTree.hpp"
#include "../utils/AdjacencyRelation.hpp"
#include "../utils/Common.hpp"
#include "../attributes/AttributeComputedIncrementally.hpp"
#include "../contours/ContoursComputedIncrementally.hpp"


namespace mmcfilters {



/**
 * @brief Estrutura auxiliar que registra informações de um extremo regional.
 *
 * Mantém o nó folha original, o nó de corte onde o extremo deixa de ser
 * dominante e o valor de extinção associado calculado durante o processo de
 * ordenação dos extremos.
 */
struct RegionalExtremaNode{
    NodeId leaf;
    NodeId cutoffNode;
    float extinction;
    
    RegionalExtremaNode(NodeId leaf, NodeId cutoffNode, float extinction) : leaf(leaf), cutoffNode(cutoffNode), extinction(extinction) { }
};

/**
 * @brief Calcula e armazena valores de extinção para extremos regionais.
 *
 * A classe percorre uma árvore morfológica e identifica as folhas que
 * sobrevivem sob diferentes critérios de supressão. Ela permite gerar mapas de
 * saliência ponderados ou discretos e aplicar filtragens mantendo apenas os
 * extremos mais persistentes.
 */
class ExtinctionValues{

    protected:
        std::vector<RegionalExtremaNode> regionalExtremaNodes;
        MorphologicalTree* tree;
        std::shared_ptr<float[]> attribute;

    public:
	    ExtinctionValues(MorphologicalTreePtr tree, std::shared_ptr<float[]> attr): ExtinctionValues(tree.get(), attr) {}    
        ExtinctionValues(MorphologicalTree* tree, std::shared_ptr<float[]> attr): tree(tree), attribute(attr) { 
            std::vector<NodeId> leaves = tree->getLeaves();
            regionalExtremaNodes.reserve(leaves.size());
            std::vector<uint8_t> visited(tree->getNumNodes(), false); //inicializa com false
            for(NodeId leaf: leaves){
                float extinction = 0;
                NodeId cutoffNode = leaf;
                NodeId parent = tree->getParentById(cutoffNode);
                bool flag = true;
                while (flag  &&  parent != InvalidNode) {
                    if (tree->getNumChildrenById(parent) > 1) {
                        for(NodeId son: tree->getChildrenById(parent) ){  // verifica se possui irmao com atributo maior
                            if(flag){
                                if (visited[son]  &&  son != cutoffNode  &&  attr[son] == attr[cutoffNode]) { //EMPATE Grimaud,92
                                    flag = false;
                                }
                                else if (son != cutoffNode  &&  attr[son] > attr[cutoffNode]) {
                                    flag = false;
                                }
                                visited[son] = true;
                            }
                        }
                    }
                    if (flag) {
                        cutoffNode = parent;
                        parent = tree->getParentById(cutoffNode);
                    }
                }
                if(parent != InvalidNode)
                    extinction = attr[cutoffNode];
                regionalExtremaNodes.emplace_back( leaf, cutoffNode, extinction );
                
            }

            // Ordena pelas extremas mais persistentes
            std::sort(regionalExtremaNodes.begin(), regionalExtremaNodes.end(), [](const auto& a, const auto& b) {
                return a.extinction > b.extinction;
            });
        }

        ImageFloatPtr saliencyMap(int extremaToKeep, bool unweighted = true){
            std::vector<uint8_t> keep(tree->getNumNodes(), false); //inicializa com false
            std::vector<float> extinctionByNode(tree->getNumNodes(), 0.0f);
            int leafToKeep = std::min(extremaToKeep, static_cast<int>(regionalExtremaNodes.size()));
            for (int i = 0; i < leafToKeep; ++i) {
                NodeId cutoffNode = this->regionalExtremaNodes[i].cutoffNode;
                keep[cutoffNode] = true;
                if(unweighted)
                    extinctionByNode[cutoffNode] = leafToKeep - i; // importance (discrete)    
                else
                    extinctionByNode[cutoffNode] = this->regionalExtremaNodes[i].extinction; //importance    
            }

            ImageFloatPtr imgOutputPtr = ImageFloat::create(tree->getNumRowsOfImage(), tree->getNumColsOfImage(), 0);
            auto saliencyOutput = imgOutputPtr->rawData();

            auto contours = ContoursComputedIncrementally::extractCompactContours(tree);
            for (auto&& [node, contour] : contours.contoursLazy()) {
                if (keep[node]) {
                    for (int p : contour) {
                        saliencyOutput[p] = extinctionByNode[node];
                    }
                }
            }

            return imgOutputPtr;
        }


        ImageUInt8Ptr filtering(int extremaToKeep){
            std::vector<uint8_t> criterion(tree->getNumNodes(), false); //inicializa com false
            int leafToKeep = std::min(extremaToKeep, static_cast<int>(regionalExtremaNodes.size()));
            for(int i=0; i < leafToKeep; i++){
                criterion[regionalExtremaNodes[i].leaf] = true;
            }
            for(NodeId node: tree->getIteratorPostOrderTraversalById()){
                NodeId parent = tree->getParentById(node);
                if (parent != InvalidNode && criterion[node]) {
                    criterion[parent] = true;
                }
            }
            ImageUInt8Ptr imgOutputPtr = ImageUInt8::create(tree->getNumRowsOfImage(), tree->getNumColsOfImage(), 0);
            auto imgOutput = imgOutputPtr->rawData();
            std::stack<NodeId> s;
            s.push(tree->getRootById());
            while(!s.empty()){
                NodeId node = s.top(); s.pop();
                int level = tree->getLevelById(node);
                for (int pixel : tree->getCNPsById(node)){
                    imgOutput[pixel] = level;
                }
                for (NodeId child: tree->getChildrenById(node)){
                    if(criterion[child]){
                        s.push(child);
                    }else{
                        for(int pixel: tree->getPixelsOfCCById(child)){
                            imgOutput[pixel] = level;
                        }
                    }
                }
            }
            return imgOutputPtr;
        }

        std::vector<RegionalExtremaNode>& getExtinctionValues() { return regionalExtremaNodes; }

};

} // namespace mmcfilters
