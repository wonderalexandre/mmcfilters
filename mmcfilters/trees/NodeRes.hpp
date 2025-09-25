#pragma once

#include "../trees/MorphologicalTree.hpp"
#include "../utils/Common.hpp"

namespace mmcfilters {

class NodeRes;
using NodeResPtr = std::shared_ptr<NodeRes>;

/**
 * @brief Representa um nó da árvore residual utilizada em filtragens por
 * resíduos.
 *
 * Cada instância armazena relações de parentesco, o conjunto de nós do
 * subárvore residual e informações auxiliares para reconstruir pixels
 * associados às regiões preservadas durante a filtragem.
 */
class NodeRes : public std::enable_shared_from_this<NodeRes> {

    private:
        NodeResPtr parent;
        std::list<NodeResPtr> children;
        int associeatedIndex;
        bool desirableResidue;
        int levelNodeNotInNR;
        MorphologicalTree* tree; //ponteiro para a árvore original
        //Nr(i) is subtree of the component tree where the variable root is the root of the subtree
        NodeId rootNr; //first node in Nr(i)
        std::list<NodeId> nodes; //nodes belongs to Nr(i)
    public:
        
        NodeRes(MorphologicalTree* tree, NodeId rootNr, int associeatedIndex, bool desirableResidue): 
         associeatedIndex(associeatedIndex), desirableResidue(desirableResidue), tree(tree), rootNr(rootNr) {}

        void addNodeInNr(NodeId node){
            this->nodes.push_back(node);
        }

        void addChild(NodeResPtr child){
            this->children.push_back(child);
        }

        void setParent(NodeResPtr parent){
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

        std::list<NodeResPtr> getChildren(){
            return this->children;
        }

        NodeId getRootNr(){
            return this->rootNr;
        }

        NodeResPtr getParent(){
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

} // namespace mmcfilters

