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

class AttributeFilters{
    protected:
        MorphologicalTreePtr tree;

    public:

    AttributeFilters(MorphologicalTreePtr tree);

    ~AttributeFilters();

    std::vector<bool> getAdaptativeCriterion(std::vector<bool>& criterion, int delta);

    int* filteringByPruningMin(float* attr, float threshold);

    int* filteringByPruningMax(float* attr, float threshold);

    int* filteringByPruningMin(std::vector<bool>& criterion);

    int* filteringByPruningMax(std::vector<bool>& criterion);

    int* filteringByDirectRule(std::vector<bool>& criterion);

    int* filteringBySubtractiveRule(std::vector<bool>& criterion);

    float* filteringBySubtractiveScoreRule(std::vector<float>& prob);

    static void filteringBySubtractiveScoreRule(MorphologicalTreePtr tree, std::vector<float>& prob, float *imgOutput){
        std::unique_ptr<float[]> mapLevel(new float[tree->getNumNodes()]);
        
        //the root is always kept
        mapLevel[0] = tree->getRoot()->getLevel();

        for(NodeMTPtr node: tree->getIndexNode()){
            if(node->getParent() != nullptr){ 
                int h = (int)std::abs(node->getLevel() - node->getParent()->getLevel());
                mapLevel[node->getIndex()] = (float) mapLevel[node->getParent()->getIndex()] + (h * prob[node->getIndex()]);
            }

        }
        for(NodeMTPtr node: tree->getIndexNode()){
            for (int pixel : node->getCNPs()){
                imgOutput[pixel] = mapLevel[node->getIndex()];
            }
        }
    }



    static void filteringByResidualRule(ResidualTree* rtree, float* attribute, float threshold, int* imgOutput){
        std::stack<NodeRes*> s;
        for (NodeRes *node : rtree->getRoot()->getChildren()){
            s.push(node);
        }
        MorphologicalTreePtr ctree = rtree->getCTree();
        std::unique_ptr<int[]> mapLevel(new int[ctree->getNumNodes()]);
        for(NodeMTPtr nodeCT: ctree->getIndexNode()){
            mapLevel[nodeCT->getIndex()] = 0;
        } 

        while (!s.empty()){
            NodeRes *node = s.top(); s.pop();
            for (NodeMTPtr nodeCT : node->getNodeInNr()){
                if(nodeCT->getParent() != nullptr){
                    if(attribute[node->getRootNr()->getIndex()] > threshold)
                        mapLevel[nodeCT->getIndex()] =  mapLevel[nodeCT->getParent()->getIndex()] + nodeCT->getResidue();
                    else
                        mapLevel[nodeCT->getIndex()] =  mapLevel[nodeCT->getParent()->getIndex()];
                }
            }            
            for (NodeRes *child : node->getChildren()){
                s.push(child);
            }
        }

        int* restOfImage = rtree->getRestOfImage();
        for(NodeMTPtr node:  ctree->getIndexNode()){
            for (int pixel : node->getCNPs()){
                if(ctree->isMaxtree())
                    imgOutput[pixel] = restOfImage[pixel] + mapLevel[node->getIndex()];
                else
                    imgOutput[pixel] = restOfImage[pixel] - mapLevel[node->getIndex()];
            }
        }

    }

    static void filteringBySubtractiveRule(MorphologicalTreePtr tree, std::vector<bool>& criterion, int *imgOutput){
        std::unique_ptr<int[]> mapLevel(new int[tree->getNumNodes()]);
        //the root is always kept
        mapLevel[0] = tree->getRoot()->getLevel();

        for(NodeMTPtr node: tree->getIndexNode()){
            if(node->getParent() != nullptr){ 
                if(criterion[node->getIndex()]){
                    int h = (int)std::abs(node->getLevel() - node->getParent()->getLevel());
                    if(!node->isMaxtreeNode())
                        h = -h;
                    mapLevel[node->getIndex()] = mapLevel[node->getParent()->getIndex()] + h;
                }
                else
                    mapLevel[node->getIndex()] = mapLevel[node->getParent()->getIndex()];
            }

        }
        for(NodeMTPtr node: tree->getIndexNode()){
            for (int pixel : node->getCNPs()){
                imgOutput[pixel] = mapLevel[node->getIndex()];
            }
        }
    }

    static void filteringByDirectRule(MorphologicalTreePtr tree, std::vector<bool>& criterion, int *imgOutput){
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
        for(NodeMTPtr node: tree->getIndexNode()){
            for (int pixel : node->getCNPs()){
                imgOutput[pixel] = mapLevel[node->getIndex()];
            }
        }
        /*std::stack<NodeMTPtr> s;
        s.push(tree->getRoot());
        std::stack<int> sLevel;
        sLevel.push(tree->getRoot()->getLevel());
        criterion[0] = true; //the root is always kept
        
        while(!s.empty()){
            NodeMTPtr node = s.top(); s.pop();
            int level = sLevel.top(); sLevel.pop();
            for (int pixel : node->getCNPs()){
                imgOutput[pixel] = level;
            }

            for (NodeMTPtr child: node->getChildren()){
                s.push(child);
                if(criterion[child->getIndex()]){
                    sLevel.push(child->getLevel());
                }else{
                    sLevel.push(level);
                }
            }
        }*/
    }

    static void filteringByPruningMin(MorphologicalTreePtr tree, std::vector<bool>& criterion, int *imgOutput){
        std::stack<NodeMTPtr> s;
        s.push(tree->getRoot());
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

    static void filteringByPruningMax(MorphologicalTreePtr tree, std::vector<bool>& _criterion, int *imgOutput){
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


    static void filteringByPruningMin(MorphologicalTreePtr tree, float *attribute, float threshold, int *imgOutput){
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

    static void filteringByPruningMax(MorphologicalTreePtr tree, float *attribute, float threshold, int *imgOutput){
        
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


    static std::vector<bool> getAdaptativeCriterion(MorphologicalTreePtr tree, double *attribute, float threshold, int delta){
		
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