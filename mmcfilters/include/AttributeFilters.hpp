#include "../include/Common.hpp"
#include "../include/NodeMT.hpp"
#include "../include/MorphologicalTree.hpp"
#include "../include/AttributeComputedIncrementally.hpp"
#include "../include/NodeRes.hpp"
#include "../include/AttributeOpeningPrimitivesFamily.hpp"
#include "../include/ResidualTree.hpp"
#include "../include/ComputerMSER.hpp"

#include <stack>
#include <vector>
#include <limits.h>



#ifndef ATTRIBUTE_FILTERS_H
#define ATTRIBUTE_FILTERS_H

#define UNDEF -999999999999

class AttrbuteFilters;
using AttributeFiltersPtr = std::shared_ptr<AttrbuteFilters>;

class AttributeFilters{
    protected:
        MorphologicalTreePtr tree;

    public:

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

    ImageUInt8Ptr filteringByExtinctionValue(MorphologicalTreePtr tree, std::shared_ptr<float[]> attribute, int numLeaf);

    static void filteringBySubtractiveScoreRule(MorphologicalTreePtr tree, std::vector<float>& prob, ImageFloatPtr imgOutputPtr){
        std::unique_ptr<float[]> mapLevel(new float[tree->getNumNodes()]);
        
        //the root is always kept
        mapLevel[0] = tree->getRoot()->getLevel();

        for(NodeMTPtr node: tree->getIndexNode()){
            if(node->getParent() != nullptr){ 
                int residue = node->getResidue();
                if(node->isMaxtreeNode())
                    mapLevel[node->getIndex()] =  (float)mapLevel[node->getParent()->getIndex()] + (residue * prob[node->getIndex()]);
                else
                    mapLevel[node->getIndex()] = (float) mapLevel[node->getParent()->getIndex()] - (residue * prob[node->getIndex()]);
            }
        }
        auto imgOutput = imgOutputPtr->rawData();
        for(NodeMTPtr node: tree->getIndexNode()){
            for (int pixel : node->getCNPs()){
                imgOutput[pixel] = mapLevel[node->getIndex()];
            }
        }
    }


    static void filteringByResidualRule(ResidualTree* rtree, std::shared_ptr<float[]> attribute, float threshold, ImageUInt8Ptr imgOutputPtr){
        std::stack<NodeResPtr> s;
        for (NodeResPtr node : rtree->getRoot()->getChildren()){
            s.push(node);
        }
        MorphologicalTreePtr ctree = rtree->getCTree();
        std::unique_ptr<int[]> mapLevel(new int[ctree->getNumNodes()]);
        for(NodeMTPtr nodeCT: ctree->getIndexNode()){
            mapLevel[nodeCT->getIndex()] = 0;
        } 

        while (!s.empty()){
            NodeResPtr node = s.top(); s.pop();
            for (NodeMTPtr nodeCT : node->getNodeInNr()){
                if(nodeCT->getParent() != nullptr){
                    if(attribute[node->getRootNr()->getIndex()] > threshold)
                        mapLevel[nodeCT->getIndex()] =  mapLevel[nodeCT->getParent()->getIndex()] + nodeCT->getResidue();
                    else
                        mapLevel[nodeCT->getIndex()] =  mapLevel[nodeCT->getParent()->getIndex()];
                }
            }            
            for (NodeResPtr child : node->getChildren()){
                s.push(child);
            }
        }

        auto imgOutput = imgOutputPtr->rawData();
        auto restOfImage = rtree->getRestOfImage()->rawData();
        for(NodeMTPtr node:  ctree->getIndexNode()){
            for (int pixel : node->getCNPs()){
                if(ctree->isMaxtree())
                    imgOutput[pixel] = restOfImage[pixel] + mapLevel[node->getIndex()];
                else
                    imgOutput[pixel] = restOfImage[pixel] - mapLevel[node->getIndex()];
            }
        }

    }

    static void filteringBySubtractiveRule(MorphologicalTreePtr tree, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr){
        std::unique_ptr<int[]> mapLevel(new int[tree->getNumNodes()]);
        //the root is always kept
        mapLevel[0] = tree->getRoot()->getLevel();

        for(NodeMTPtr node: tree->getIndexNode()){
            if(node->getParent() != nullptr){ 
                if(criterion[node->getIndex()]){
                    if(node->isMaxtreeNode())
                        mapLevel[node->getIndex()] = mapLevel[node->getParent()->getIndex()] + node->getResidue();
                    else
                        mapLevel[node->getIndex()] = mapLevel[node->getParent()->getIndex()] - node->getResidue();
                }
                else
                    mapLevel[node->getIndex()] = mapLevel[node->getParent()->getIndex()];
            }

        }

        auto imgOutput = imgOutputPtr->rawData();
        for(NodeMTPtr node: tree->getIndexNode()){
            for (int pixel : node->getCNPs()){
                imgOutput[pixel] = mapLevel[node->getIndex()];
            }
        }
    }

    static void filteringByDirectRule(MorphologicalTreePtr tree, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr){
        std::unique_ptr<int[]> mapLevel(new int[tree->getNumNodes()]);

        //the root is always kept
        mapLevel[0] = tree->getRoot()->getLevel();

        for(NodeMTPtr node: tree->getIndexNode()){
            if(node->getParent() != nullptr){ 
                if(criterion[node->getIndex()])
                    mapLevel[node->getIndex()] = node->getLevel();
                else
                    mapLevel[node->getIndex()] = mapLevel[node->getParent()->getIndex()];
            }

        }
        auto imgOutput = imgOutputPtr->rawData();
        for(NodeMTPtr node: tree->getIndexNode()){
            for (int pixel : node->getCNPs()){
                imgOutput[pixel] = mapLevel[node->getIndex()];
            }
        }
    }

    static void filteringByPruningMin(MorphologicalTreePtr tree, std::vector<bool>& criterion, ImageUInt8Ptr imgOutputPtr){
        std::stack<NodeMTPtr> s;
        s.push(tree->getRoot());
        auto imgOutput = imgOutputPtr->rawData();
        while(!s.empty()){
            NodeMTPtr node = s.top(); s.pop();
            for (int pixel : node->getCNPs()){
                imgOutput[pixel] = node->getLevel();;
            }
            for (NodeMTPtr child: node->getChildren()){
                if(criterion[child->getIndex()]){
                    s.push(child);
                }else{
                    for(int pixel: child->getPixelsOfCC()){
                        imgOutput[pixel] = child->getLevel();
                    }
                }
            }
        }
    }

    static void filteringByPruningMax(MorphologicalTreePtr tree, std::vector<bool>& _criterion, ImageUInt8Ptr imgOutputPtr){
        std::unique_ptr<bool[]> criterion(new bool[tree->getNumNodes()]);
        AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
            [&criterion, _criterion](NodeMTPtr node) -> void { //pre-processing
                if(!_criterion[node->getIndex()])
                    criterion[node->getIndex()] = true;
                else
                    criterion[node->getIndex()] = false;
            },
            [&criterion](NodeMTPtr parent, NodeMTPtr child) -> void { 
                criterion[parent->getIndex()] = (criterion[parent->getIndex()] & criterion[child->getIndex()]);
            },
            [](NodeMTPtr node) -> void { //post-processing
                                        
            }
        );
        auto imgOutput = imgOutputPtr->rawData();
        std::stack<NodeMTPtr> s;
        s.push(tree->getRoot());
        while(!s.empty()){
            NodeMTPtr node = s.top(); s.pop();
            for (int pixel : node->getCNPs()){
                imgOutput[pixel] = node->getLevel();
            }
            for (NodeMTPtr child: node->getChildren()){
                if(!criterion[child->getIndex()]){
                    s.push(child);
                }else{
                    for(int pixel: child->getPixelsOfCC()){
                        imgOutput[pixel] = child->getLevel();
                    }
                }
            }
        }
    }

    
    static void filteringByPruningMin(MorphologicalTreePtr tree, std::shared_ptr<float[]> attribute, float threshold, ImageUInt8Ptr imgOutputPtr){
        auto imgOutput = imgOutputPtr->rawData();
        std::stack<NodeMTPtr> s;
        s.push(tree->getRoot());
        while(!s.empty()){
            NodeMTPtr node = s.top(); s.pop();
            for (int pixel : node->getCNPs()){
                imgOutput[pixel] = node->getLevel();
            }
            for (NodeMTPtr child: node->getChildren()){
                if(attribute[child->getIndex()] > threshold){
                    s.push(child);
                }else{
                    for(int pixel: child->getPixelsOfCC()){
                        imgOutput[pixel] =  node->getLevel();
                    }
                }
                
            }
        }
    }

    
    static void filteringByPruningMax(MorphologicalTreePtr tree, std::shared_ptr<float[]> attribute, float threshold, ImageUInt8Ptr imgOutputPtr){
        
        std::unique_ptr<bool[]> criterion(new bool[tree->getNumNodes()]);
        AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
            [&criterion, attribute, threshold](NodeMTPtr node) -> void { //pre-processing
                if(attribute[node->getIndex()] <= threshold)
                    criterion[node->getIndex()] = true;
                else
                    criterion[node->getIndex()] = false;
            },
            [&criterion, attribute, threshold](NodeMTPtr parent, NodeMTPtr child) -> void { 
                criterion[parent->getIndex()] = (criterion[parent->getIndex()] & criterion[child->getIndex()]);
            },
            [&criterion, attribute, threshold](NodeMTPtr node) -> void { //post-processing
                                        
            }
        );
        auto imgOutput = imgOutputPtr->rawData();
        std::stack<NodeMTPtr> s;
        s.push(tree->getRoot());
        while(!s.empty()){
            NodeMTPtr node = s.top(); s.pop();
            for (int pixel : node->getCNPs()){
                imgOutput[pixel] = node->getLevel();
            }
            for (NodeMTPtr child: node->getChildren()){
                if(!criterion[child->getIndex()]){
                    s.push(child);
                }else{
                    for(int pixel: child->getPixelsOfCC()){
                        imgOutput[pixel] =  node->getLevel();
                    }
                }
            }
        }
    }


    static std::vector<bool> getAdaptativeCriterion(MorphologicalTreePtr tree, std::shared_ptr<float[]> attribute, float threshold, int delta){
		
        ComputerMSER mser(tree);
		std::vector<bool> isMSER = mser.computerMSER(delta);

		std::vector<double> stability = mser.getStabilities();
		std::vector<bool> isPruned(tree->getNumNodes(), false);
		for(NodeMTPtr node: tree->getIndexNode()){
            if(attribute[node->getIndex()] < threshold){ //node pruned

                if(stability[node->getIndex()] == UNDEF){
                    isPruned[node->getIndex()] = true;
                }else{
                    
                    //NodeMTPtr nodeMax = mser.getNodeInPathWithMaxStability(node, isMSER);
                    //isPruned[nodeMax->getIndex()] = true;
                    
                    double max = stability[node->getIndex()];
                    int indexDescMaxStability = mser.descendantWithMaxStability(node)->getIndex();
                    int indexAscMaxStability = mser.ascendantWithMaxStability(node)->getIndex();
                    double maxDesc = stability[indexDescMaxStability];
                    double maxAnc = stability[indexAscMaxStability];
                    
                    if(max >= maxDesc && max >= maxAnc) {
                        isPruned[node->getIndex()] = true;
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

    static std::vector<bool> getAdaptativeCriterion(MorphologicalTreePtr tree, std::vector<bool>& criterion, int delta){
		
        ComputerMSER mser(tree);
		std::vector<bool> isMSER = mser.computerMSER(delta);

		std::vector<double> stability = mser.getStabilities();
		std::vector<bool> isPruned(tree->getNumNodes(), false);
		for(NodeMTPtr node: tree->getIndexNode()){
            if(!criterion[node->getIndex()]){ //node pruned

                if(stability[node->getIndex()] == UNDEF){
                    isPruned[node->getIndex()] = true;
                }else{
                    
                    //NodeMTPtr nodeMax = mser.getNodeInPathWithMaxStability(node, isMSER);
                    //isPruned[nodeMax->getIndex()] = true;
                    
                    double max = stability[node->getIndex()];
                    int indexDescMaxStability = mser.descendantWithMaxStability(node)->getIndex();
                    int indexAscMaxStability = mser.ascendantWithMaxStability(node)->getIndex();
                    double maxDesc = stability[indexDescMaxStability];
                    double maxAnc = stability[indexAscMaxStability];
                    
                    if(max >= maxDesc && max >= maxAnc) {
                        isPruned[node->getIndex()] = true;
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


#endif