#

#include "../include/ComponentTree.hpp"
#include "../include/AdjacencyRelation.hpp"
#include "../include/Common.hpp"
#include "../include/AttributeComputedIncrementally.hpp"
#include <algorithm>

#define PI 3.14159265358979323846


struct RegionalExtremaNodeCT{
    NodeId leaf;
    NodeId cutoffNode;
    float extinction;
    
    RegionalExtremaNodeCT(NodeId leaf, NodeId cutoffNode, float extinction) : leaf(leaf), cutoffNode(cutoffNode), extinction(extinction) { }
};

class ExtinctionValuesCT{

    private:
        std::vector<RegionalExtremaNodeCT> regionalExtremaNodes;
        ComponentTree* tree;
        std::shared_ptr<float[]> attribute;

    public:
	    ExtinctionValuesCT(ComponentTreePtr tree, std::shared_ptr<float[]> attr): ExtinctionValuesCT(tree.get(), attr) {}    
        ExtinctionValuesCT(ComponentTree* tree, std::shared_ptr<float[]> attr): tree(tree), attribute(attr) { 
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
                    extinction = attr[cutoffNode->getIndex()];
                regionalExtremaNodes.push_back( leaf, cutoffNode, extinction );
                
            }

            // Ordena pelas extremas mais persistentes
            std::sort(regionalExtremaNodes.begin(), regionalExtremaNodes.end(), [](const auto& a, const auto& b) {
                return a->extinction > b->extinction;
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

            auto contoursCT = AttributeComputedIncrementallyCT::extractCompactContours(tree);
            for (auto&& [node, contourNode] : contoursCT->contoursLazy()) {
                if (keep[node]) {
                    for (int p : contourNode) {
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
            for(NodeId node: tree->getNodeIds()){
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

        std::vector<RegionalExtremaNodeCT>& getExtinctionValues() { return regionalExtremaNodes; }

};

#endif // EXTINCTION_VALUES_H


