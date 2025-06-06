#include "../include/ExtinctionValues.hpp"
#include "../include/AttributeComputedIncrementally.hpp"
#include <algorithm>

ExtinctionValues::ExtinctionValues(MorphologicalTreePtr tree, std::shared_ptr<float[]> attr): tree(tree), attribute(attr) { 
    std::list<NodeMTPtr> leaves = tree->getLeaves();
    regionalExtremaNodes.reserve(leaves.size());
    std::unique_ptr<bool[]> visited(new bool[tree->getNumNodes()]()); //inicializa com false
    for(NodeMTPtr leaf: leaves){
        float extinction = 0;
        NodeMTPtr cutoffNode = leaf;
        NodeMTPtr parent = cutoffNode->getParent();
        bool flag = true;
        while (flag  &&  parent != nullptr) {
            if (parent->getChildren().size() > 1) {
                for(NodeMTPtr son: parent->getChildren()){  // verifica se possui irmao com atributo maior
                    if(flag){
                        if (visited[son->getIndex()]  &&  son != cutoffNode  &&  attr[son->getIndex()] == attr[cutoffNode->getIndex()]) { //EMPATE Grimaud,92
                            flag = false;
                        }
                        else if (son != cutoffNode  &&  attr[son->getIndex()] > attr[cutoffNode->getIndex()]) {
                            flag = false;
                        }
                        visited[son->getIndex()] = true;
                    }
                }
            }
            if (flag) {
                cutoffNode = parent;
                parent = cutoffNode->getParent();
            }
        }
        if(parent != nullptr)
            extinction = attr[cutoffNode->getIndex()];
        regionalExtremaNodes.push_back( std::make_shared<RegionalExtremaNode>(leaf, cutoffNode, extinction) );
        
    }

    // Ordena pelas extremas mais persistentes
    std::sort(regionalExtremaNodes.begin(), regionalExtremaNodes.end(), [](const auto& a, const auto& b) {
        return a->extinction > b->extinction;
    });
}

ImageFloatPtr ExtinctionValues::saliencyMap(int extremaToKeep, bool unweighted){ 
    std::unique_ptr<bool[]> keep(new bool[tree->getNumNodes()]());
    std::vector<float> extinctionByNode(tree->getNumNodes(), 0.0f);

    int leafToKeep = std::min(extremaToKeep, static_cast<int>(regionalExtremaNodes.size()));
    for (int i = 0; i < leafToKeep; ++i) {
        auto cutoffNode = this->regionalExtremaNodes[i]->cutoffNode;
        keep[cutoffNode->getIndex()] = true;
        if(unweighted)
            extinctionByNode[cutoffNode->getIndex()] = leafToKeep - i; // importance (discrete)    
        else
            extinctionByNode[cutoffNode->getIndex()] = this->regionalExtremaNodes[i]->extinction; //importance    
    }

    ImageFloatPtr imgOutputPtr = ImageFloat::create(tree->getNumRowsOfImage(), tree->getNumColsOfImage(), 0);
    auto saliencyOutput = imgOutputPtr->rawData();

    auto contoursMT = AttributeComputedIncrementally::extractCompactContours(tree);
    for (auto&& [node, contourNode] : contoursMT->contoursLazy()) {

        if (keep[node->getIndex()]) {
            for (int p : contourNode) {
                saliencyOutput[p] = extinctionByNode[node->getIndex()];
            }
        }
    }

    return imgOutputPtr;
}


ImageUInt8Ptr ExtinctionValues::filtering(int extremaToKeep){
    
    std::unique_ptr<bool[]> criterion(new bool[tree->getNumNodes()]());
    int leafToKeep = std::min(extremaToKeep, static_cast<int>(regionalExtremaNodes.size()));
    for(int i=0; i < leafToKeep; i++){
        criterion[regionalExtremaNodes[i]->leaf->getIndex()] = true;
    }
    for(NodeMTPtr node: tree->getRoot()->getIteratorPostOrderTraversal()){
        NodeMTPtr parent = node->getParent();
        if (parent && criterion[node->getIndex()]) {
            criterion[parent->getIndex()] = true;
        }
    }
    ImageUInt8Ptr imgOutputPtr = ImageUInt8::create(tree->getNumRowsOfImage(), tree->getNumColsOfImage(), 0);
    auto imgOutput = imgOutputPtr->rawData();
    std::stack<NodeMTPtr> s;
    s.push(tree->getRoot());
    while(!s.empty()){
        NodeMTPtr node = s.top(); s.pop();
        int level = node->getLevel();
        for (int pixel : node->getCNPs()){
            imgOutput[pixel] = level;
        }
        for (NodeMTPtr child: node->getChildren()){
            if(criterion[child->getIndex()]){
                s.push(child);
            }else{
                for(int pixel: child->getPixelsOfCC()){
                    imgOutput[pixel] = level;
                }
            }
        }
    }
    return imgOutputPtr;
}
