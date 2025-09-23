#pragma once

#include "../include/Common.hpp"
#include "../include/AttributeComputedIncrementally.hpp"
#include "../include/MorphologicalTree.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>



/**
 * @brief Estrutura auxiliar para construir contornos compactos a partir da árvore.
 *
 * Mantém listas de pixels que devem ser adicionados ou removidos por nó e
 * provê percursos convenientes para inspecionar contornos, componentes conexos
 * e pares ordenados durante a execução incremental.
 */
class Contours{
private:
	std::vector<std::list<int>> contours;
	std::vector<std::list<int>> contoursToRemove;
	MorphologicalTree* tree;
	
public:
	Contours(MorphologicalTree* tree): tree(tree), contours(tree->getNumNodes()), contoursToRemove(tree->getNumNodes()){}
    Contours(MorphologicalTreePtr tree): Contours(tree.get()) {}

	void add(NodeId node, int pixel){
		contours[node].push_back(pixel);
	}
	void remove(NodeId node, int pixel){
		contoursToRemove[node].push_back(pixel);
	}

	std::unordered_set<int> getContour(NodeId nodeSubtree) {
		std::unordered_set<int> contour;
		AttributeComputedIncrementally::computerAttribute(tree, nodeSubtree,
			[](NodeId node) -> void {},  // pre-processing
			[](NodeId parent, NodeId child) -> void { }, // merge-processing
			[&contour, this](NodeId node) -> void { //post-processing
				for(int p: this->contours[node]){
					contour.insert(p);
				}
				for(int p: this->contoursToRemove[node]){
					contour.erase(p);
				}
			}
		);
		return contour;
	}

	void visitContours(std::function<void(NodeId, const std::unordered_set<int>&)> visitor) {
		const int numNodes = tree->getNumNodes();
	
		std::vector<std::unique_ptr<std::unordered_set<int>>> contoursByNodes(numNodes);
	
		AttributeComputedIncrementally::computerAttribute(tree, tree->getRootById(),
			[](NodeId) -> void {},
	
			// merge: funde filhos no pai, usando o maior conjunto como base
			[&contoursByNodes](NodeId parent, NodeId child) -> void {
				auto& parentContour = contoursByNodes[parent];
				auto& childContour = contoursByNodes[child];
	
				if (!parentContour) {
					parentContour = std::move(childContour);
				} else {
					if (childContour->size() > parentContour->size()) {
						std::swap(parentContour, childContour);
					}
					parentContour->insert(childContour->begin(), childContour->end());
					childContour.reset(); 
				}
			},
	
			// pós-processamento: consolida e chama visitor
			[this, &contoursByNodes, &visitor](NodeId node) -> void {
				auto& contour = contoursByNodes[node];
				if (!contour) {
					contour = std::make_unique<std::unordered_set<int>>();
				}
	
				for (int p : this->contours[node]) {
					contour->insert(p);
				}
	
				for (int p : this->contoursToRemove[node]) {
					contour->erase(p);
				}
	
				visitor(node, *contour);
			}
		);
	}

	void visitContoursAndCCs(std::function<void(NodeId, const std::list<int>&, const std::unordered_set<int>&)> visitor) {
		const int numNodes = tree->getNumNodes();
	
		std::vector<std::unique_ptr<std::unordered_set<int>>> contoursByNodes(numNodes);
		std::vector<std::unique_ptr<std::list<int>>> CCsByNodes(numNodes);
	
		AttributeComputedIncrementally::computerAttribute(tree, tree->getRootById(),
			[](NodeId) -> void {},
	
			[&CCsByNodes, &contoursByNodes](NodeId parent, NodeId child) -> void {
				// --- Contornos ---
				auto& parentContour = contoursByNodes[parent];
				auto& childContour = contoursByNodes[child];
				if (!parentContour) {
					parentContour = std::move(childContour);
				} else {
					if (childContour->size() > parentContour->size()) {
						std::swap(parentContour, childContour);
					}
					parentContour->insert(childContour->begin(), childContour->end());
					//childContour.reset();
				}

				// --- Componentes Conexos ---
				auto& parentCC = CCsByNodes[parent];
				auto& childCC = CCsByNodes[child];
				if (!parentCC) {
					parentCC = std::move(childCC);
				} else {
					if (childCC->size() > parentCC->size()) {
						std::swap(parentCC, childCC);
					}
					parentCC->insert(parentCC->end(), childCC->begin(), childCC->end());
					//childCC.reset();
				}
			},
	
			// post-processing
			[this, &contoursByNodes, &CCsByNodes, &visitor](NodeId node) -> void {
				// --- Contornos ---
				auto& contour = contoursByNodes[node];
				if (!contour) {
					contour = std::make_unique<std::unordered_set<int>>();
				}
				for (int p : this->contours[node]) {
					contour->insert(p);
				}
				for (int p : this->contoursToRemove[node]) {
					contour->erase(p);
				}
				
				// --- Componentes Conexos ---
				auto& cc = CCsByNodes[node];
				if (!cc) {
					cc = std::make_unique<std::list<int>>();
				}
                
                //TODO: Melhorar a transferencia dos CNPs
                for(int p: this->tree->getCNPsById(node))
                    cc->push_back(p);
				//cc->insert(cc->end(), cnpsNode.begin(), cnpsNode.end());
	
				
				visitor(node, *cc, *contour);
			}
		);
	}


        /**
         * @brief Iterador pós-ordem que fornece pares (nó, contorno) consolidados.
         */
    class ContourPostOrderIterator {
	private:
		using value_type = std::pair<NodeId, std::unordered_set<int>>;
		using reference = value_type&;
		using pointer = value_type*;
		using iterator_category = std::input_iterator_tag;

		Contours* contoursCT;
        MorphologicalTree* tree;
		std::stack<NodeId> outputStack;
		std::vector<std::unique_ptr<std::unordered_set<int>>> contoursByNodes;
		value_type currentValue;

		void advance() {
			if (!outputStack.empty()) {
				NodeId node = outputStack.top(); outputStack.pop();

				// Merge dos filhos (igual antes)
				for (NodeId child : tree->getChildrenById(node)) {
					auto& parentContour = contoursByNodes[node];
					auto& childContour = contoursByNodes[child];
					if (!parentContour) {
						parentContour = std::move(childContour);
					} else if (childContour) {
						if (childContour->size() > parentContour->size()) {
							std::swap(parentContour, childContour);
						}
						parentContour->insert(childContour->begin(), childContour->end());
					}
					if (childContour) childContour.reset();
				}

				// Pós-processamento: aplica inserção/remoção
				auto& contour = contoursByNodes[node];
				if (!contour) contour = std::make_unique<std::unordered_set<int>>();
				for (int p : contoursCT->contours[node])
					contour->insert(p);
				for (int p : contoursCT->contoursToRemove[node])
					contour->erase(p);

				currentValue = std::make_pair(node, *contour);
			}
		}

	public:
		ContourPostOrderIterator(MorphologicalTree* tree, NodeId root, Contours* contoursCT) : tree(tree), contoursCT(contoursCT){
			// Travessia prévia para montar outputStack em pós-ordem
			if(root != InvalidNode){
				std::stack<NodeId> tempStack;
				tempStack.push(root);
				while (!tempStack.empty()) {
					NodeId current = tempStack.top(); tempStack.pop();
					outputStack.push(current);
					for (NodeId child : tree->getChildrenById(current)) {
						tempStack.push(child);
					}
				}
				int numNodes = tree->getNumDescendantsById(root) +1;
				contoursByNodes.resize(numNodes);
				advance();
			}
			
		}

		// Pré-incremento
		ContourPostOrderIterator& operator++() {
			advance();
			return *this;
		}

		reference operator*() {
			return currentValue;
		}

		bool operator==(const ContourPostOrderIterator& other) const {
			return outputStack.empty() && other.outputStack.empty();
		}
		bool operator!=(const ContourPostOrderIterator& other) const {
			return !(*this == other);
		}
		ContourPostOrderIterator(const ContourPostOrderIterator&) = delete;
		ContourPostOrderIterator& operator=(const ContourPostOrderIterator&) = delete;
		ContourPostOrderIterator(ContourPostOrderIterator&&) = default;
		ContourPostOrderIterator& operator=(ContourPostOrderIterator&&) = default;
	};

        /**
         * @brief Range que produz iteradores pós-ordem sobre contornos compactados.
         */
    class ContourPostOrderRange {
	private:
		NodeId root;
        MorphologicalTree* tree;
		Contours* contoursCT;

	public:
		ContourPostOrderRange(MorphologicalTree* tree, NodeId root, Contours* contoursCT) : root(root), tree(tree), contoursCT(contoursCT) {}

		ContourPostOrderIterator begin() { return ContourPostOrderIterator(tree, root, contoursCT); }
		ContourPostOrderIterator end() { return ContourPostOrderIterator(tree, InvalidNode, contoursCT); }
	};

	
	ContourPostOrderRange contoursLazy() {
		return ContourPostOrderRange(tree, tree->getRootById(), this);
	}

};




/**
 * @brief Funções utilitárias para computar contornos em árvores morfologicas de forma incremental.
 *
 * A classe expõe algoritmos genéricos de travessia pós-ordem permitindo
 * compor etapas de pré-processamento, mesclagem e pós-processamento sem criar
 * estruturas auxiliares temporárias. 
 */
class ContoursComputedIncrementally {
public:
    
    
    static std::shared_ptr<Contours> extractCompactContours(MorphologicalTreePtr tree){ return extractCompactContours(tree.get()); }
    static std::shared_ptr<Contours> extractCompactContours(MorphologicalTree* tree){
        std::shared_ptr<Contours> contoursMT = std::make_shared<Contours>(tree);
        
        std::vector<std::vector<int>> contoursToRemoveLCA(tree->getNumNodes());
        std::vector<std::int8_t> ncount(tree->getNumRowsOfImage() * tree->getNumColsOfImage(), 0);
        AdjacencyRelationPtr adj4 = std::make_shared<AdjacencyRelation>(tree->getNumRowsOfImage(), tree->getNumColsOfImage(), 1);
        LCAEulerRMQ lca(tree);	

        AttributeComputedIncrementally::computerAttribute(tree, tree->getRootById(),
            [](NodeId node) -> void { // pre-processing

            },
            [](NodeId parent, NodeId child) -> void { // merge-processing
                
            },
            [&contoursMT, &contoursToRemoveLCA, &lca, &ncount, tree, adj4](NodeId nodeP) -> void { // post-processing
                std::vector<int> &NcontourToRemoveLCA = contoursToRemoveLCA[nodeP];

                NodeId nodeLCA = nodeP;
                for(int p: NcontourToRemoveLCA){ //pixels que sao contornos de nodes descendentes ao NodeAtual
                    bool isPixelToBeRemoved = true;
                    for (int r : adj4->getNeighborPixels(p)) { //Existe um nodeQ ascendente de NodeAtual contendo p como contorno? (p, q) in A
                        NodeId nodeR = tree->getSCById(r); 
                        if (tree->isStrictAncestor(nodeR, nodeLCA)){
                            contoursToRemoveLCA[nodeR].push_back(p); 
                            isPixelToBeRemoved = false;	
                        }else if(!tree->isComparable(nodeLCA, nodeR)) {
                            NodeId otherNodeLCA = lca.findLowestCommonAncestor(nodeLCA, nodeR);
                            contoursToRemoveLCA[otherNodeLCA].push_back(p);
                            isPixelToBeRemoved = false;
                        }
                    }
                    if(!adj4->isBorderDomainImage(p) && isPixelToBeRemoved){
                        contoursMT->remove(nodeLCA, p);
                    }
                }
            
                for (int p : tree->getCNPsById(nodeP)) {
                    if (adj4->isBorderDomainImage(p)){
                        ncount[p]++;
                    }

                    for (int q : adj4->getNeighborPixels(p)) {
                        NodeId nodeQ = tree->getSCById(q); 
                        if(!tree->isComparable(nodeP, nodeQ)){ //se os nodeP e nodeQ não sao comparaveis, então p pode ser removido pelo LCA de nodeP e nodeQ 
                            NodeId nodeLCA = lca.findLowestCommonAncestor(nodeP, nodeQ);
                            contoursToRemoveLCA[nodeLCA].push_back(p);
                            ncount[p]++;
                        }
                        else if(tree->isStrictDescendant(nodeP, nodeQ)){  //maxtree:  SC(p) \subset SC(q) <=> f(p) > f(q)
                            ncount[p]++;
                        }else if (tree->isStrictAncestor(nodeP, nodeQ)) { ////maxtree:  SC(q) \subset SC(p) <=> f(p) < f(q)
                            ncount[q]--;
                            if (ncount[q] == 0) {
                                contoursMT->remove(nodeP, q);
                            }
                        }
                    }
                    if (ncount[p] > 0){
                        contoursMT->add(nodeP, p);
                    }
                }

            }
        );
                    
        return contoursMT;
    }

    static std::vector<std::unordered_set<int>> extractNonCompactContours(MorphologicalTreePtr tree){ return extractNonCompactContours(tree.get()); }
    static std::vector<std::unordered_set<int>> extractNonCompactContours(MorphologicalTree* tree){
        std::vector<std::unordered_set<int>> contours(tree->getNumNodes());
        std::vector<std::vector<int>> contoursToRemoveLCA(tree->getNumNodes());
        std::vector<std::int8_t> ncount(tree->getNumRowsOfImage() * tree->getNumColsOfImage(), 0);
        AdjacencyRelationPtr adj4 = std::make_shared<AdjacencyRelation>(tree->getNumRowsOfImage(), tree->getNumColsOfImage(), 1);
        LCAEulerRMQ lca(tree);	

        AttributeComputedIncrementally::computerAttribute(tree, tree->getRootById(),
            [](NodeId node) -> void { // pre-processing

            },
            [&contours, &ncount, tree, adj4](NodeId parent, NodeId child) -> void { // merge-processing
                std::unordered_set<int> &Ncontour = contours[parent];
                for (int p : contours[child]){
                    Ncontour.insert(p);
                }
            },
            [&contours, &contoursToRemoveLCA, &lca, &ncount, tree, adj4](NodeId nodeP) -> void { // post-processing
                // Initialise contours of node "N"
                std::unordered_set<int> &Ncontour = contours[nodeP];
                std::vector<int> &NcontourToRemoveLCA = contoursToRemoveLCA[nodeP];
                NodeId nodeLCA = nodeP;
                for(int p: NcontourToRemoveLCA){ //pixels que sao contornos de nodes descendentes ao NodeAtual
                    bool isPixelToBeRemoved = true;
                    
                    for (int r : adj4->getNeighborPixels(p)) { //Existe um nodeQ ascendente de NodeAtual contendo p como contorno? (p, q) in A
                        NodeId nodeR = tree->getSCById(r); 
                        if (tree->isStrictAncestor(nodeR, nodeLCA)){
                            contoursToRemoveLCA[nodeR].push_back(p); 
                            isPixelToBeRemoved = false;	
                        }else if(!tree->isComparable(nodeLCA, nodeR)) {
                            NodeId otherNodeLCA = lca.findLowestCommonAncestor(nodeLCA, nodeR);
                            contoursToRemoveLCA[otherNodeLCA].push_back(p);
                            isPixelToBeRemoved = false;
                        }
                    }
                    if(!adj4->isBorderDomainImage(p) && isPixelToBeRemoved){
                        Ncontour.erase(p);
                    }
                }
            
                for (int p : tree->getCNPsById(nodeP)) {
                    if (adj4->isBorderDomainImage(p)){
                        ncount[p]++;
                    }
                    for (int q : adj4->getNeighborPixels(p)) {
                        NodeId nodeQ = tree->getSCById(q); 
                        if(!tree->isComparable(nodeP, nodeQ)){ //se os nodeP e nodeQ não sao comparaveis, então p pode ser removido pelo LCA de nodeP e nodeQ 
                            NodeId nodeLCA = lca.findLowestCommonAncestor(nodeP, nodeQ);
                            contoursToRemoveLCA[nodeLCA].push_back(p);
                            ncount[p]++;
                        }
                        else if(tree->isStrictDescendant(nodeP, nodeQ)){  //maxtree:  SC(p) \subset SC(q) <=> f(p) > f(q)
                            ncount[p]++;
                        }else if (tree->isStrictAncestor(nodeP, nodeQ)) { ////maxtree:  SC(q) \subset SC(p) <=> f(p) < f(q)
                            ncount[q]--;
                            if (ncount[q] == 0) {
                                Ncontour.erase(q);
                            }
                        }
                    }

                    if (ncount[p] > 0){
                        Ncontour.insert(p);
                    }
                }

            }
        );
                    
        return contours;
    }


};