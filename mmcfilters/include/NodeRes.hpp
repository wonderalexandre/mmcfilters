#include <list>
#include <algorithm>

#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"

#ifndef NODE_RES_H
#define NODE_RES_H

class NodeRes {

    private:
        NodeRes* parent;
        std::list<NodeRes*> children;
        int associeatedIndex;
        bool desirableResidue;
        int levelNodeNotInNR;

        //Nr(i) is subtree of the component tree where the variable root is the root of the subtree
        NodeCT* rootNr; //first node in Nr(i)
        std::list<NodeCT*> nodes; //nodes belongs to Nr(i)
    public:
        
        NodeRes(NodeCT* rootNr, int associeatedIndex, bool desirableResidue);

        void addNodeInNr(NodeCT* node);
        void addChild(NodeRes* child);
        void setParent(NodeRes* parent);
        NodeRes* getParent();
        int getAssocieatedIndex();
        bool isDesirableResidue();
        std::list<NodeRes*> getChildren();
        std::list<NodeCT*> getNodeInNr();
        NodeCT* getRootNr();
        int getLevelNodeNotInNR();
        void setLevelNodeNotInNR(int level);
        bool belongsToNr(NodeCT* node){
            return std::find(this->nodes.begin(), this->nodes.end(), node) != this->nodes.end();
        }

        
        
   class InternalIteratorPixelsOfCNPs{
		private:
			NodeCT *currentNode;
			std::stack<NodeCT*> s;
			std::list<int>::iterator iter;
			int countCNPs;
			using iterator_category = std::input_iterator_tag;
            using value_type = int; 
		public:
			InternalIteratorPixelsOfCNPs(NodeRes* instance, int numCNPs)  {
				this->countCNPs = numCNPs;
				for (NodeCT *child: instance->nodes){
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
			NodeRes *instance;
			int sumCPNs;
		public:
			IteratorPixelsOfCNPs(NodeRes *obj, int _sumCNPs): instance(obj), sumCPNs(_sumCNPs) {}
			InternalIteratorPixelsOfCNPs begin(){ return InternalIteratorPixelsOfCNPs(instance, 0); }
            InternalIteratorPixelsOfCNPs end(){ return InternalIteratorPixelsOfCNPs(instance, sumCPNs); }
	};	
	IteratorPixelsOfCNPs& getPixelsOfCNPs(){
        int sumCPNs = 0;
        for(NodeCT* node: this->nodes){
            sumCPNs += node->getCNPs().size();
        }
	    IteratorPixelsOfCNPs *iter = new IteratorPixelsOfCNPs(this, sumCPNs);
    	return *iter;
	}


        
};


#endif