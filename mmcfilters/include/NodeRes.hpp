#include <list>
#include <algorithm>

#include "../include/NodeMT.hpp"
#include "../include/MorphologicalTree.hpp"
#include "../include/Common.hpp"

#ifndef NODE_RES_H
#define NODE_RES_H

class NodeRes;
using NodeResPtr = std::shared_ptr<NodeRes>;

class NodeRes : public std::enable_shared_from_this<NodeRes> {

    private:
        NodeResPtr parent;
        std::list<NodeResPtr> children;
        int associeatedIndex;
        bool desirableResidue;
        int levelNodeNotInNR;

        //Nr(i) is subtree of the component tree where the variable root is the root of the subtree
        NodeMTPtr rootNr; //first node in Nr(i)
        std::list<NodeMTPtr> nodes; //nodes belongs to Nr(i)
    public:
        
        NodeRes(NodeMTPtr rootNr, int associeatedIndex, bool desirableResidue);

        void addNodeInNr(NodeMTPtr node);
        void addChild(NodeResPtr child);
        void setParent(NodeResPtr parent);
        NodeResPtr getParent();
        int getAssocieatedIndex();
        bool isDesirableResidue();
        std::list<NodeResPtr> getChildren();
        std::list<NodeMTPtr> getNodeInNr();
        NodeMTPtr getRootNr();
        int getLevelNodeNotInNR();
        void setLevelNodeNotInNR(int level);
        bool belongsToNr(NodeMTPtr node){
            return std::find(this->nodes.begin(), this->nodes.end(), node) != this->nodes.end();
        }

        
        
   class InternalIteratorPixelsOfCNPs{
		private:
			NodeMTPtr currentNode;
			std::stack<NodeMTPtr> s;
			std::list<int>::iterator iter;
			int countCNPs;
			using iterator_category = std::input_iterator_tag;
            using value_type = int; 
		public:
			InternalIteratorPixelsOfCNPs(NodeResPtr instance, int numCNPs)  {
				this->countCNPs = numCNPs;
				for (NodeMTPtr child: instance->nodes){
					s.push(child);
				}
                this->currentNode = s.top(); s.pop();
                this->iter = this->currentNode->getCNPs().begin();
			}
			InternalIteratorPixelsOfCNPs& operator++() { 
			    this->iter++; 
				if(this->iter == this->currentNode->getCNPs().end()){
					if(!s.empty()){
            			this->currentNode = s.top(); s.pop();
						this->iter = this->currentNode->getCNPs().begin();
					}
				}
				this->countCNPs++;
				return *this; 
            }
            bool operator==(InternalIteratorPixelsOfCNPs other) const { 
                return this->countCNPs == other.countCNPs; 
            }
            bool operator!=(InternalIteratorPixelsOfCNPs other) const { 
                return !(*this == other);
            }
            int operator*() const { 
                return (*this->iter); 
            }  
    };
	class IteratorPixelsOfCNPs{
		private:
            NodeResPtr instance;
			int sumCPNs;
		public:
			IteratorPixelsOfCNPs(NodeResPtr obj, int _sumCNPs): instance(obj), sumCPNs(_sumCNPs) {}
			InternalIteratorPixelsOfCNPs begin(){ return InternalIteratorPixelsOfCNPs(instance, 0); }
            InternalIteratorPixelsOfCNPs end(){ return InternalIteratorPixelsOfCNPs(instance, sumCPNs); }
	};	
	IteratorPixelsOfCNPs& getPixelsOfCNPs(){
        int sumCPNs = 0;
        for(NodeMTPtr node: this->nodes){
            sumCPNs += node->getCNPs().size();
        }
	    IteratorPixelsOfCNPs *iter = new IteratorPixelsOfCNPs(this->shared_from_this(), sumCPNs);
    	return *iter;
	}


        
};


#endif