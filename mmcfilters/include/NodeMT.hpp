#include <list>
#include <stack>
#include <iterator>
#include <utility>
#include <memory>
#include "../include/Common.hpp"

#ifndef NODECT_H
#define NODECT_H


class NodeMT : public std::enable_shared_from_this<NodeMT> {
private:
	int index; //tempo de entrada durante o percurso pos-ordem (1 incremento) 
    int level;
	int areaCC;
    int timePostOrder;  //tempo de entrada durante o percurso pos-ordem (2 incremento) 
    int timePreOrder;  //tempo de saia durante o percurso pos-ordem (2 incremento) 
	NodeMTPtr parent;
    std::list<NodeMTPtr> children;
	std::list<int> cnps;
public:
	
    NodeMT();
    NodeMT(int index, NodeMTPtr parent, int level);
    void addCNPs(int p);
    void setLevel(int level);
    void addChild(NodeMTPtr child);
	int getIndex();
	void setIndex(int index);
	int getResidue();
	int getLevel();
	int getAreaCC();
	bool isMaxtreeNode();
	void setAreaCC(int area);
	int getNumDescendants();
	NodeMTPtr getParent();
	void setParent(NodeMTPtr parent);
	std::list<int>& getCNPs();
	std::list<NodeMTPtr>& getChildren();
	int getNumSiblings();
    int getTimePostOrder();
    int getTimePreOrder();
    void setTimePostOrder(int time);
    void setTimePreOrder(int time);
    bool isLeaf();

	
//============= Iterator para iterar os nodes do caminho até o root==============//
class InternalIteratorNodesOfPathToRoot {
    private:
        NodeMTPtr currentNode;
    
    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = NodeMTPtr;
        using difference_type = std::ptrdiff_t;
        using pointer = NodeMTPtr*;
        using reference = NodeMTPtr&;
    
        InternalIteratorNodesOfPathToRoot(NodeMTPtr obj) : currentNode(obj) {}
    
        InternalIteratorNodesOfPathToRoot& operator++() {
            if (currentNode) {
                currentNode = currentNode->getParent();  // Retorna outro shared_ptr
            }
            return *this;
        }
    
        bool operator==(const InternalIteratorNodesOfPathToRoot& other) const {
            return currentNode == other.currentNode;
        }
    
        bool operator!=(const InternalIteratorNodesOfPathToRoot& other) const {
            return !(*this == other);
        }
    
        reference operator*() {
            return currentNode;
        }
    };
    
    class IteratorNodesOfPathToRoot {
    private:
        NodeMTPtr instance;
    
    public:
        explicit IteratorNodesOfPathToRoot(NodeMTPtr obj) : instance(obj) {}
    
        InternalIteratorNodesOfPathToRoot begin() const { return InternalIteratorNodesOfPathToRoot(instance); }
        InternalIteratorNodesOfPathToRoot end() const { return InternalIteratorNodesOfPathToRoot(nullptr); }
    };
    
    // Chamador usa this como shared_ptr:
    IteratorNodesOfPathToRoot getNodesOfPathToRoot() {
        return IteratorNodesOfPathToRoot(this->shared_from_this());
    }
    

/////////
// **Iterador para coletar ramos em pós-ordem**
// Classe do Iterador Pós-Ordem por Ramos
class InternalIteratorBranchPostOrderTraversal {
    private:
        std::stack<NodeMTPtr> processingStack;
        std::stack<NodeMTPtr> postOrderStack;
        std::list<std::list<NodeMTPtr>> branches;
        typename std::list<std::list<NodeMTPtr>>::iterator branchIterator;
    
    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = std::list<NodeMTPtr>;
        using pointer = std::list<NodeMTPtr>*;
        using reference = std::list<NodeMTPtr>&;
    
        InternalIteratorBranchPostOrderTraversal(NodeMTPtr root) {
            if (!root) return;
    
            std::stack<NodeMTPtr> tempStack;
            tempStack.push(root);
    
            while (!tempStack.empty()) {
                auto current = tempStack.top();
                tempStack.pop();
                postOrderStack.push(current);
    
                for (auto& child : current->getChildren()) {
                    tempStack.push(child);
                }
            }
    
            std::list<NodeMTPtr> currentBranch;
            while (!postOrderStack.empty()) {
                auto node = postOrderStack.top();
                postOrderStack.pop();
    
                if (!currentBranch.empty()) {
                    auto lastNode = currentBranch.back();
                    if (lastNode->getParent() && lastNode->getParent()->getChildren().back() != lastNode) {
                        branches.push_back(currentBranch);
                        currentBranch.clear();
                    }
                }
    
                currentBranch.push_back(node);
            }
    
            if (!currentBranch.empty()) {
                branches.push_back(currentBranch);
            }
    
            branchIterator = branches.begin();
        }
    
        InternalIteratorBranchPostOrderTraversal& operator++() {
            if (branchIterator != branches.end()) {
                ++branchIterator;
            }
            return *this;
        }
    
        reference operator*() {
            return *branchIterator;
        }
    
        bool operator==(const InternalIteratorBranchPostOrderTraversal& other) const {
            return branchIterator == other.branchIterator;
        }
    
        bool operator!=(const InternalIteratorBranchPostOrderTraversal& other) const {
            return !(*this == other);
        }
    };
    
    // Classe externa
    class IteratorBranchPostOrderTraversal {
    private:
        NodeMTPtr root;
    public:
        explicit IteratorBranchPostOrderTraversal(NodeMTPtr root) : root(root) {}
    
        InternalIteratorBranchPostOrderTraversal begin() { return InternalIteratorBranchPostOrderTraversal(root); }
        InternalIteratorBranchPostOrderTraversal end() { return InternalIteratorBranchPostOrderTraversal(nullptr); }
    };
    
    // Método da classe
    IteratorBranchPostOrderTraversal getIteratorBranchPostOrderTraversal() {
        return IteratorBranchPostOrderTraversal(this->shared_from_this());
    }
    


//============= Iterator para iterar os nodes de um percuso em pos-ordem ==============//
	class InternalIteratorPostOrderTraversal {
    private:
        std::stack<NodeMTPtr> nodeStack;
        std::stack<NodeMTPtr> outputStack;
    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = NodeMTPtr;
        using difference_type = std::ptrdiff_t;
        using pointer = NodeMTPtr*;
        using reference = NodeMTPtr&; 

        InternalIteratorPostOrderTraversal(NodeMTPtr root) {
            if (root) {
                nodeStack.push(root);
                while (!nodeStack.empty()) {
                    NodeMTPtr current = nodeStack.top();nodeStack.pop();
                    outputStack.push(current);
                    for (NodeMTPtr child : current->getChildren()) {
                        nodeStack.push(child);
                    }
                }
            }
        }

        InternalIteratorPostOrderTraversal& operator++() {
            if (!outputStack.empty()) {
                outputStack.pop();
            }
            return *this;
        }

        reference operator*() {  
            return outputStack.top();  // Retorna ponteiro!
        }

        bool operator==(const InternalIteratorPostOrderTraversal& other) const {
            return (outputStack.empty() == other.outputStack.empty());
        }

        bool operator!=(const InternalIteratorPostOrderTraversal& other) const {
            return !(*this == other);
        }
    };

	class IteratorPostOrderTraversal {
    private:
        NodeMTPtr root;
    public:
        explicit IteratorPostOrderTraversal(NodeMTPtr root) : root(root) {}

        InternalIteratorPostOrderTraversal begin() { return InternalIteratorPostOrderTraversal(root); }
        InternalIteratorPostOrderTraversal end() { return InternalIteratorPostOrderTraversal(nullptr); }
    };

    IteratorPostOrderTraversal getIteratorPostOrderTraversal() { return IteratorPostOrderTraversal(this->shared_from_this()); }



//============= Iterator para iterar os nodes de um percuso em largura ==============//
    class InternalIteratorBreadthFirstTraversal {
    private:
        std::queue<NodeMTPtr> nodeQueue;

    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = NodeMTPtr;
        using difference_type = std::ptrdiff_t;
        using pointer = NodeMTPtr*;
        using reference = NodeMTPtr&; // Retorna ponteiro!

        InternalIteratorBreadthFirstTraversal(NodeMTPtr root) {
            if (root) {
                nodeQueue.push(root);
            }
        }

        InternalIteratorBreadthFirstTraversal& operator++() {
            if (!nodeQueue.empty()) {
                NodeMTPtr current = nodeQueue.front();
                nodeQueue.pop();
                for (NodeMTPtr child : current->getChildren()) {
                    nodeQueue.push(child);
                }
            }
            return *this;
        }

        reference operator*() {
            return nodeQueue.front();
        }

        bool operator==(const InternalIteratorBreadthFirstTraversal& other) const {
            return nodeQueue.empty() == other.nodeQueue.empty();
        }

        bool operator!=(const InternalIteratorBreadthFirstTraversal& other) const {
            return !(*this == other);
        }
    };

    class IteratorBreadthFirstTraversal {
    private:
    NodeMTPtr root;

    public:
        explicit IteratorBreadthFirstTraversal(NodeMTPtr root) : root(root) {}

        InternalIteratorBreadthFirstTraversal begin() { return InternalIteratorBreadthFirstTraversal(root); }
        InternalIteratorBreadthFirstTraversal end() { return InternalIteratorBreadthFirstTraversal(nullptr); }
    };

    // Método para expor o iterador na classe NodeCT
    IteratorBreadthFirstTraversal getIteratorBreadthFirstTraversal() { 
        return IteratorBreadthFirstTraversal(this->shared_from_this()); 
    }




//============= Iterator para iterar os pixels de um CC==============//
class InternalIteratorPixelsOfCC{
    private:
        NodeMTPtr currentNode;
        std::stack<NodeMTPtr> s;

        std::list<int>::iterator iter;
        int countArea;
        using iterator_category = std::input_iterator_tag;
        using value_type = int; 
    public:
        InternalIteratorPixelsOfCC(NodeMTPtr obj, int area)  {
            this->currentNode = obj;
            this->countArea =area;
            this->iter = this->currentNode->cnps.begin();
            for (NodeMTPtr child: this->currentNode->getChildren()){
                s.push(child);
            }	
        }
        InternalIteratorPixelsOfCC& operator++() { 
            this->iter++; 
            if(this->iter == this->currentNode->cnps.end()){
                if(!s.empty()){
                    this->currentNode = s.top(); s.pop();
                    this->iter = this->currentNode->cnps.begin();
                    for (NodeMTPtr child: currentNode->getChildren()){
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

    class IteratorPixelsOfCC {
    private:
        NodeMTPtr instance;
        int area;
    public:
        explicit IteratorPixelsOfCC(NodeMTPtr obj, int _area) : instance(obj), area(_area) {}

        auto begin() {
            return InternalIteratorPixelsOfCC(instance, 0);
        }

        auto end() {
            return InternalIteratorPixelsOfCC(instance, area);
            
        }

    };
    IteratorPixelsOfCC getPixelsOfCC() {
        return IteratorPixelsOfCC(this->shared_from_this(), this->areaCC);
    }




////////////////////////////////////////////////
    class InternalIteratorNodesDescendants{
		private:
			NodeMTPtr currentNode;
			std::stack<NodeMTPtr> s;
			int numDescendants;
			using iterator_category = std::input_iterator_tag;
            using value_type = NodeMT; 
		public:
			InternalIteratorNodesDescendants(NodeMTPtr obj, int numDescendants)  {
				this->numDescendants = numDescendants;
				this->currentNode = obj;
				for (NodeMTPtr child: obj->getChildren()){
					s.push(child);
				}
					
			}
			InternalIteratorNodesDescendants& operator++() { 
			    if(!s.empty()){
            		this->currentNode = s.top(); s.pop();
					for (NodeMTPtr child: currentNode->getChildren()){
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
            NodeMTPtr operator*() { 
                return (this->currentNode); 
            }  
    };
	class IteratorNodesDescendants{
		private:
			NodeMTPtr instance;
			int numDescendants;
		public:
			IteratorNodesDescendants(NodeMTPtr obj, int _numDescendants): instance(obj), numDescendants(_numDescendants) {}
			InternalIteratorNodesDescendants begin(){ return InternalIteratorNodesDescendants(instance, 0); }
            InternalIteratorNodesDescendants end(){ return InternalIteratorNodesDescendants(instance, numDescendants+1); }
	};	
	IteratorNodesDescendants& getNodesDescendants(){
	    IteratorNodesDescendants *iter = new IteratorNodesDescendants(this->shared_from_this(), this->getNumDescendants());
    	return *iter;
	}
	
};

#endif