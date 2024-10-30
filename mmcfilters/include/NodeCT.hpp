#include <list>
#include <stack>
#include <iterator>
#include <utility>

#ifndef NODECT_H
#define NODECT_H

class NodeCT {
private:
	int index; 
    int rep;
    int level;
	int residue;
	NodeCT* parent;
	std::list<int> cnps;
    std::list<NodeCT*> children;
    bool isMaxtree;
	int areaCC;
	int numDescendants;

public:
	
    NodeCT();
    NodeCT(int index, int rep, NodeCT* parent, int level);
    void addCNPs(int p);
    void addChild(NodeCT* child);
	int getRep();
	int getIndex();
	void setIndex(int index);
	int getResidue();
	void setResidue(int residue);
	int getLevel();
	int getAreaCC();
	bool isMaxtreeNode();
	void setAreaCC(int area);
	int getNumDescendants();
	void setNumDescendants(int num);
	NodeCT* getParent();
	void setParent(NodeCT* parent);
	std::list<int> getCNPs();
	std::list<NodeCT*> getChildren();
	int getNumSiblings();

	
///////////////////////////////////////////////////
    class InternalIteratorPixelsOfCC{
		private:
			NodeCT *currentNode;
			std::stack<NodeCT*> s;
			std::list<int>::iterator iter;
			int countArea;
			using iterator_category = std::input_iterator_tag;
            using value_type = int; 
		public:
			InternalIteratorPixelsOfCC(NodeCT *obj, int area)  {
				this->currentNode = obj;
				this->countArea =area;
				this->iter = this->currentNode->cnps.begin();
				for (NodeCT *child: this->currentNode->getChildren()){
					s.push(child);
				}	
			}
			InternalIteratorPixelsOfCC& operator++() { 
			    this->iter++; 
				if(this->iter == this->currentNode->cnps.end()){
					if(!s.empty()){
            			this->currentNode = s.top(); s.pop();
						this->iter = this->currentNode->cnps.begin();
						for (NodeCT *child: currentNode->getChildren()){
                		    s.push(child);
						}
					}
				}
				this->countArea++;
				return *this; 
            }
            bool operator==(InternalIteratorPixelsOfCC other) const { 
                return this->countArea == other.countArea; 
            }
            bool operator!=(InternalIteratorPixelsOfCC other) const { 
                return !(*this == other);
            }
            int operator*() const { 
                return (*this->iter); 
            }  
    };
	class IteratorPixelsOfCC{
		private:
			NodeCT *instance;
			int area;
		public:
			IteratorPixelsOfCC(NodeCT *obj, int _area): instance(obj), area(_area) {}
			InternalIteratorPixelsOfCC begin(){ return InternalIteratorPixelsOfCC(instance, 0); }
            InternalIteratorPixelsOfCC end(){ return InternalIteratorPixelsOfCC(instance, area); }
	};	
	IteratorPixelsOfCC& getPixelsOfCC(){
	    IteratorPixelsOfCC *iter = new IteratorPixelsOfCC(this, this->areaCC);
    	return *iter;
	}



///////////////////////////////////////////////////
	class InternalIteratorNodesOfPathToRoot{
		private:
			NodeCT *currentNode;
			int index;
			using iterator_category = std::input_iterator_tag;
            using value_type = NodeCT; 
		public:
			InternalIteratorNodesOfPathToRoot(NodeCT *obj, int index)  {
				this->currentNode = obj;
				this->index = index;	
			}
			InternalIteratorNodesOfPathToRoot& operator++() { 
				this->index = this->currentNode->index;
				if(this->currentNode != nullptr){
					this->currentNode = this->currentNode->parent;
				}
				return *this; 
			}
			bool operator==(InternalIteratorNodesOfPathToRoot other) const { 
                return this->index == other.index; 
            }
            bool operator!=(InternalIteratorNodesOfPathToRoot other) const { 
                return !(*this == other);
            }
            NodeCT* operator*()  { 
                return (this->currentNode); 
            }  
	};
	class IteratorNodesOfPathToRoot{
		private:
			NodeCT *instance;
		public:
			IteratorNodesOfPathToRoot(NodeCT *obj): instance(obj){}
			InternalIteratorNodesOfPathToRoot begin(){ return InternalIteratorNodesOfPathToRoot(instance, instance->index); }
            InternalIteratorNodesOfPathToRoot end(){ return InternalIteratorNodesOfPathToRoot(instance, 0); }
	};
	IteratorNodesOfPathToRoot& getNodesOfPathToRoot(){
	    IteratorNodesOfPathToRoot *iter = new IteratorNodesOfPathToRoot(this);
    	return *iter;
	}

////////////////////////////////////////////////
    class InternalIteratorNodesDescendants{
		private:
			NodeCT *currentNode;
			std::stack<NodeCT*> s;
			int numDescendants;
			using iterator_category = std::input_iterator_tag;
            using value_type = NodeCT; 
		public:
			InternalIteratorNodesDescendants(NodeCT *obj, int numDescendants)  {
				this->numDescendants = numDescendants;
				this->currentNode = obj;
				for (NodeCT *child: obj->getChildren()){
					s.push(child);
				}
					
			}
			InternalIteratorNodesDescendants& operator++() { 
			    if(!s.empty()){
            		this->currentNode = s.top(); s.pop();
					for (NodeCT *child: currentNode->getChildren()){
            		    s.push(child);
					}
				}
				this->numDescendants += 1;
				return *this; 
            }
            bool operator==(InternalIteratorNodesDescendants other) const { 
                return this->numDescendants == other.numDescendants; 
            }
            bool operator!=(InternalIteratorNodesDescendants other) const { 
                return !(*this == other);
            }
            NodeCT* operator*() { 
                return (this->currentNode); 
            }  
    };
	class IteratorNodesDescendants{
		private:
			NodeCT *instance;
			int numDescendants;
		public:
			IteratorNodesDescendants(NodeCT *obj, int _numDescendants): instance(obj), numDescendants(_numDescendants) {}
			InternalIteratorNodesDescendants begin(){ return InternalIteratorNodesDescendants(instance, 0); }
            InternalIteratorNodesDescendants end(){ return InternalIteratorNodesDescendants(instance, numDescendants+1); }
	};	
	IteratorNodesDescendants& getNodesDescendants(){
	    IteratorNodesDescendants *iter = new IteratorNodesDescendants(this, this->numDescendants);
    	return *iter;
	}
	
};

#endif