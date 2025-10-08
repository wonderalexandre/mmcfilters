#pragma once

#include "../utils/Common.hpp"
#include "../trees/MorphologicalTree.hpp"
#include "../attributes/AttributeComputedIncrementally.hpp"
#include "../attributes/ComputerMSER.hpp"
#include "../trees/NodeRes.hpp"
#include "../trees/ResidualTree.hpp"



namespace mmcfilters {


class AttributeFilters;
using AttributeFiltersPtr = std::shared_ptr<AttributeFilters>;

/**
 * @brief Conjunto de operações de filtragem baseadas em atributos sobre árvores morfológicas.
 *
 * Implementa diferentes regras (direta, subtrativa, poda, extinção) para
 * reconstruir imagens a partir de critérios definidos sobre atributos dos nós.
 */
class AttributeFilters{
    protected:
        MorphologicalTree* tree;

    public:

    AttributeFilters(MorphologicalTree* tree);
    AttributeFilters(MorphologicalTreePtr tree);

    ~AttributeFilters();

    std::vector<bool> getAdaptativeCriterion(std::vector<bool>& criterion, int delta);

    ImageUInt8Ptr filteringByPruningMin(std::shared_ptr<float[]> attr, float threshold);

    ImageUInt8Ptr filteringByPruningMax(std::shared_ptr<float[]> attr, float threshold);

    ImageUInt8Ptr filteringByPruningMin(std::vector<bool>& criterion);

    ImageUInt8Ptr filteringByPruningMax(std::vector<bool>& criterion);

    ImageUInt8Ptr filteringByDirectRule(std::vector<bool>& criterion);

    ImageUInt8Ptr filteringBySubtractiveRule(std::vector<bool>& criterion);

    ImageFloatPtr filteringBySubtractiveScoreRule(std::vector<float>& prob);


    static void filteringBySubtractiveScoreRule(MorphologicalTreePtr tree, std::vector<float>& prob, ImageFloatPtr imgOutputPtr){ return filteringBySubtractiveScoreRule(tree.get(), prob, imgOutputPtr);}
    static void filteringBySubtractiveScoreRule(MorphologicalTree* tree, std::vector<float>& prob, ImageFloatPtr imgOutputPtr){
        std::unique_ptr<float[]> mapLevel(new float[tree->getNumNodes()]);
        
        //the root is always kept
        mapLevel[0] = tree->getLevelById( tree->getRootById() );

        for(NodeId node: tree->getNodeIds()){
            if(tree->getParentById(node) != InvalidNode){ 
                int residue = tree->getResidueById(node);
                mapLevel[node] =  (float)mapLevel[tree->getParentById(node)] + (residue * prob[node]);
            }
        }
        auto imgOutput = imgOutputPtr->rawData();
        for(NodeId node: tree->getNodeIds()){
            for (int pixel : tree->getCNPsById(node)){
                imgOutput[pixel] = mapLevel[node];
            }
        }
    }

    
    static void filteringByResidualRule(ResidualTree* rtree, std::shared_ptr<float[]> attribute, float threshold, ImageUInt8Ptr imgOutputPtr){
        std::stack<NodeResPtr> s;
        for (NodeResPtr node : rtree->getRoot()->getChildren()){
            s.push(node);
        }
        MorphologicalTree* ctree = rtree->getCTree();
        std::unique_ptr<int[]> mapLevel(new int[ctree->getNumNodes()]);
        for(NodeId nodeCT: ctree->getNodeIds()){
            mapLevel[nodeCT] = 0;
        } 

        while (!s.empty()){
            NodeResPtr node = s.top(); s.pop();
            for (NodeId nodeCT : node->getNodeInNr()){
                if(ctree->getParentById(nodeCT) != InvalidNode){
                    if(attribute[node->getRootNr()] > threshold)
                        mapLevel[nodeCT] =  mapLevel[ctree->getParentById(nodeCT)] + ctree->getResidueById(nodeCT);
                    else
                        mapLevel[nodeCT] =  mapLevel[ctree->getParentById(nodeCT)];
                }
            }            
            for (NodeResPtr child : node->getChildren()){
                s.push(child);
            }
        }

        auto imgOutput = imgOutputPtr->rawData();
        auto restOfImage = rtree->getRestOfImage()->rawData();
        for(NodeId node:  ctree->getNodeIds()){
            for (int pixel : ctree->getCNPsById(node)){
                if(ctree->isMaxtree())
                    imgOutput[pixel] = restOfImage[pixel] + mapLevel[node];
                else
                    imgOutput[pixel] = restOfImage[pixel] - mapLevel[node];
            }
        }

    }


    static void filteringBySubtractiveRule(MorphologicalTreePtr tree, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr){ return filteringBySubtractiveRule(tree.get(), criterion, imgOutputPtr); }
    static void filteringBySubtractiveRule(MorphologicalTree* tree, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr){
        std::unique_ptr<int[]> mapLevel(new int[tree->getNumNodes()]);
        //the root is always kept
        mapLevel[0] = tree->getLevelById( tree->getRootById() );

        for(NodeId node: tree->getNodeIds()){
            if(tree->getParentById(node) != InvalidNode){ 
                if(criterion[node]){
                    mapLevel[node] = mapLevel[tree->getParentById(node)] + tree->getResidueById(node);
                }
                else
                    mapLevel[node] = mapLevel[tree->getParentById(node)];
            }

        }

        auto imgOutput = imgOutputPtr->rawData();
        for(NodeId node: tree->getNodeIds()){
            for (int pixel : tree->getCNPsById(node)){
                imgOutput[pixel] = mapLevel[node];
            }
        }
    }

    static void filteringByDirectRule(MorphologicalTreePtr tree, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr){ return filteringByDirectRule(tree.get(), criterion, imgOutputPtr); }
    static void filteringByDirectRule(MorphologicalTree* tree, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr){
        std::unique_ptr<int[]> mapLevel(new int[tree->getNumNodes()]);

        //the root is always kept
        mapLevel[0] = tree->getLevelById( tree->getRootById() );

        for(NodeId node: tree->getNodeIds()){
            if(tree->getParentById(node) != InvalidNode){ 
                if(criterion[node])
                    mapLevel[node] = tree->getLevelById(node);
                else
                    mapLevel[node] = mapLevel[tree->getParentById(node)];
            }

        }
        auto imgOutput = imgOutputPtr->rawData();
        for(NodeId node: tree->getNodeIds()){
            for (int pixel : tree->getCNPsById(node)){
                imgOutput[pixel] = mapLevel[node];
            }
        }
    }

    static void filteringByPruningMin(MorphologicalTreePtr tree, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr){ return filteringByPruningMin(tree.get(), criterion, imgOutputPtr); }
    static void filteringByPruningMin(MorphologicalTree* tree, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr){
        std::stack<NodeId> s;
        s.push(tree->getRootById());
        auto imgOutput = imgOutputPtr->rawData();
        while(!s.empty()){
            NodeId node = s.top(); s.pop();
            for (int pixel : tree->getCNPsById(node)){
                imgOutput[pixel] = tree->getLevelById(node);;
            }
            for (NodeId child: tree->getChildrenById(node)){
                if(criterion[child]){
                    s.push(child);
                }else{
                    for(int pixel:  tree->getPixelsOfCCById(child)){
                        imgOutput[pixel] = tree->getLevelById(child);
                    }
                }
            }
        }
    }

    static void filteringByPruningMax(MorphologicalTreePtr tree, std::vector<bool>& _criterion, ImageUInt8Ptr imgOutputPtr){ return filteringByPruningMax(tree.get(), _criterion, imgOutputPtr); }
    static void filteringByPruningMax(MorphologicalTree* tree, std::vector<bool>& _criterion, ImageUInt8Ptr imgOutputPtr){
        std::vector<uint8_t> criterion(tree->getNumNodes(), false);
        AttributeComputedIncrementally::computerAttribute(tree, tree->getRootById(),
            [&criterion, _criterion](NodeId node) -> void { //pre-processing
                if(!_criterion[node])
                    criterion[node] = true;
                else
                    criterion[node] = false;
            },
            [&criterion](NodeId parent, NodeId child) -> void { 
                criterion[parent] = (criterion[parent] & criterion[child]);
            },
            [](NodeId) -> void { //post-processing
                                        
            }
        );
        auto imgOutput = imgOutputPtr->rawData();
        std::stack<NodeId> s;
        s.push(tree->getRootById());
        while(!s.empty()){
            NodeId node = s.top(); s.pop();
            for (int pixel : tree->getCNPsById(node)){
                imgOutput[pixel] = tree->getLevelById(node);
            }
            for (NodeId child: tree->getChildrenById(node)){
                if(!criterion[child]){
                    s.push(child);
                }else{
                    for(int pixel: tree->getPixelsOfCCById(child)){
                        imgOutput[pixel] = tree->getLevelById(child);
                    }
                }
            }
        }
    }

    static void filteringByPruningMin(MorphologicalTreePtr tree, std::shared_ptr<float[]> attribute, float threshold, ImageUInt8Ptr imgOutputPtr){ return filteringByPruningMin(tree.get(), attribute, threshold, imgOutputPtr); }
    static void filteringByPruningMin(MorphologicalTree* tree, std::shared_ptr<float[]> attribute, float threshold, ImageUInt8Ptr imgOutputPtr){
        auto imgOutput = imgOutputPtr->rawData();
        std::stack<NodeId> s;
        s.push(tree->getRootById());
        while(!s.empty()){
            NodeId node = s.top(); s.pop();
            for (int pixel : tree->getCNPsById(node)){
                imgOutput[pixel] = tree->getLevelById(node);
            }
            for (NodeId child: tree->getChildrenById(node)){
                if(attribute[child] > threshold){
                    s.push(child);
                }else{
                    for(int pixel: tree->getPixelsOfCCById(child)){
                        imgOutput[pixel] =  tree->getLevelById(node);
                    }
                }
                
            }
        }
    }

    static void filteringByPruningMax(MorphologicalTreePtr tree, std::shared_ptr<float[]> attribute, float threshold, ImageUInt8Ptr imgOutputPtr){ return filteringByPruningMax(tree.get(), attribute, threshold, imgOutputPtr); }
    static void filteringByPruningMax(MorphologicalTree* tree, std::shared_ptr<float[]> attribute, float threshold, ImageUInt8Ptr imgOutputPtr){
        std::vector<uint8_t> criterion(tree->getNumNodes(), false);
        AttributeComputedIncrementally::computerAttribute(tree, tree->getRootById(),
            [&criterion, attribute, threshold](NodeId node) -> void { //pre-processing
                if(attribute[node] <= threshold)
                    criterion[node] = true;
            },
            [&criterion](NodeId parent, NodeId child) -> void { 
                criterion[parent] = (criterion[parent] & criterion[child]);
            },
            [](NodeId) -> void { //post-processing
                                        
            }
        );
        auto imgOutput = imgOutputPtr->rawData();
        std::stack<NodeId> s;
        s.push(tree->getRootById());
        while(!s.empty()){
            NodeId node = s.top(); s.pop();
            for (int pixel : tree->getCNPsById(node)){
                imgOutput[pixel] = tree->getLevelById(node);
            }
            for (NodeId child: tree->getChildrenById(node)){
                if(!criterion[child]){
                    s.push(child);
                }else{
                    for(int pixel: tree->getPixelsOfCCById(child)){
                        imgOutput[pixel] =  tree->getLevelById(node);
                    }
                }
            }
        }
    }

    static std::vector<bool> getAdaptativeCriterion(MorphologicalTreePtr tree, std::shared_ptr<float[]> attribute, float threshold, int delta){ return getAdaptativeCriterion(tree.get(), attribute, threshold, delta); }
    static std::vector<bool> getAdaptativeCriterion(MorphologicalTree* tree, std::shared_ptr<float[]> attribute, float threshold, int delta){
		
        ComputerMSER mser(tree);
		std::vector<uint8_t> isMSER = mser.computerMSER(delta);

		std::vector<float> stability = mser.getStabilities();
		std::vector<bool> isPruned(tree->getNumNodes(), false);
		for(NodeId node: tree->getNodeIds()){
            if(attribute[node] < threshold){ //node pruned

                if(std::isnan(stability[node])){
                    isPruned[node] = true;
                }else{
                    
                    //NodeId nodeMax = mser.getNodeInPathWithMaxStability(node, isMSER);
                    //isPruned[nodeMax] = true;
                    
                    float max = stability[node];
                    NodeId indexDescMaxStability = mser.descendantWithMaxStability(node);
                    NodeId indexAscMaxStability = mser.ascendantWithMaxStability(node);
                    float maxDesc = stability[indexDescMaxStability];
                    float maxAnc = stability[indexAscMaxStability];
                    
                    if(max >= maxDesc && max >= maxAnc) {
                        isPruned[node] = true;
                    }else if (maxDesc >= max && maxDesc >= maxAnc) {
                        isPruned[indexDescMaxStability] = true;
                    }else {
                        isPruned[indexAscMaxStability] = true;
                    }
                    
                }
			}
			
		}
        return isPruned;
    }

    static std::vector<bool> getAdaptativeCriterion(MorphologicalTreePtr tree, std::vector<bool>& criterion, int delta){ return getAdaptativeCriterion(tree.get(), criterion, delta); }
    static std::vector<bool> getAdaptativeCriterion(MorphologicalTree* tree, std::vector<bool>& criterion, int delta){
		
        ComputerMSER mser(tree);
		std::vector<uint8_t> isMSER = mser.computerMSER(delta);

		std::vector<float> stability = mser.getStabilities();
		std::vector<bool> isPruned(tree->getNumNodes(), false);
		for(NodeId node: tree->getNodeIds()){
            if(!criterion[node]){ //node pruned

                if(std::isnan(stability[node])){
                    isPruned[node] = true;
                }else{
                    
                    //NodeId nodeMax = mser.getNodeInPathWithMaxStability(node, isMSER);
                    //isPruned[nodeMax] = true;
                    
                    float max = stability[node];
                    NodeId indexDescMaxStability = mser.descendantWithMaxStability(node);
                    NodeId indexAscMaxStability = mser.ascendantWithMaxStability(node);
                    float maxDesc = stability[indexDescMaxStability];
                    float maxAnc = stability[indexAscMaxStability];
                    
                    if(max >= maxDesc && max >= maxAnc) {
                        isPruned[node] = true;
                    }else if (maxDesc >= max && maxDesc >= maxAnc) {
                        isPruned[indexDescMaxStability] = true;
                    }else {
                        isPruned[indexAscMaxStability] = true;
                    }
                    
                }
			}
			
		}
        return isPruned;
    }

	
};

} // namespace mmcfilters

