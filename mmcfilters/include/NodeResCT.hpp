#include <list>
#include <algorithm>

#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"
#include "../include/Common.hpp"

#ifndef NODE_RES_H
#define NODE_RES_H

class NodeResCT;
using NodeResCTPtr = std::shared_ptr<NodeResCT>;

class NodeResCT : public std::enable_shared_from_this<NodeResCT> {

    private:
        NodeResCTPtr parent;
        std::list<NodeResCTPtr> children;
        int associeatedIndex;
        bool desirableResidue;
        int levelNodeNotInNR;
        ComponentTree* tree; //ponteiro para a árvore original
        //Nr(i) is subtree of the component tree where the variable root is the root of the subtree
        NodeId rootNr; //first node in Nr(i)
        std::list<NodeId> nodes; //nodes belongs to Nr(i)
    public:
        
        NodeResCT(ComponentTree* tree, NodeId rootNr, int associeatedIndex, bool desirableResidue): 
        tree(tree), rootNr(rootNr), associeatedIndex(associeatedIndex), desirableResidue(desirableResidue) {}

        void addNodeInNr(NodeId node){
            this->nodes.push_back(node);
        }

        void addChild(NodeResCTPtr child){
            this->children.push_back(child);
        }

        void setParent(NodeResCTPtr parent){
            this->parent = parent;
        }

        int getAssocieatedIndex(){
            return this->associeatedIndex;
        }

        bool isDesirableResidue(){
            return this->desirableResidue;
        }

        std::list<NodeId> getNodeInNr(){
            return this->nodes;
        }

        std::list<NodeResCTPtr> getChildren(){
            return this->children;
        }

        NodeId getRootNr(){
            return this->rootNr;
        }

        NodeResCTPtr getParent(){
            return this->parent;
        }

        int getLevelNodeNotInNR(){
            return this->levelNodeNotInNR;
        }

        void setLevelNodeNotInNR(int level){
            this->levelNodeNotInNR = level;
        }
        
        bool belongsToNr(NodeId node){
            return std::find(this->nodes.begin(), this->nodes.end(), node) != this->nodes.end();
        }

        auto getPixelsOfCNPs() const{ //iterador
            std::vector<int> reps(this->nodes.size());
            for(NodeId node: this->nodes)
                reps.push_back( tree->getRepNodeById(node) ); 
            return tree->getPixelsOfFlatzones(reps);  
        } 



        
};


#endif