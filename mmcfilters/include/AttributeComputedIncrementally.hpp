#pragma once

#include "../include/Common.hpp"
#include "../include/AttributeNames.hpp"
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



//Forward declaration
class Contours; 
class AttributeComputer; 

/**
 * @brief Funções utilitárias para computar atributos em árvores de forma incremental.
 *
 * A classe expõe algoritmos genéricos de travessia pós-ordem permitindo
 * compor etapas de pré-processamento, mesclagem e pós-processamento sem criar
 * estruturas auxiliares temporárias. Também oferece helpers para extração de
 * contornos e avaliação de atributos disponíveis.
 */
class AttributeComputedIncrementally {
public:
    
 /*
    virtual void preProcessing(NodeId v);

    virtual void mergeChildren(NodeId parent, NodeId child);

    virtual void postProcessing(NodeId parent);

    void computerAttribute(MorphologicalTree* tree, NodeId root);

    static void computerAttribute(MorphologicalTree* tree, NodeId root, 
										std::function<void(NodeId)> preProcessing,
										std::function<void(NodeId, NodeId)> mergeChildren,
										std::function<void(NodeId)> postProcessing ){
		
		preProcessing(root);
		for(NodeId child: tree->getChildrenById(root)){
			AttributeComputedIncrementally::computerAttribute(tree, child, preProcessing, mergeChildren, postProcessing);
			mergeChildren(root, child);
		}
		postProcessing(root);
	}*/

    
    template<class PreProcessing, class MergeProcessing, class PostProcessing>
    static void computerAttribute(MorphologicalTree* tree, NodeId root, 
                                    PreProcessing&& preProcessing, 
                                    MergeProcessing&& mergeProcessing, 
                                    PostProcessing&& postProcessing) {
        preProcessing(root);
        for (NodeId child : tree->getChildrenById(root)) {
            AttributeComputedIncrementally::computerAttribute(tree, child, preProcessing, mergeProcessing, postProcessing); // passar por ref (sem cópia)
            mergeProcessing(root, child);
        }
        postProcessing(root);
    }


    
	static std::shared_ptr<Contours> extractCompactContours(MorphologicalTree* tree);
    static std::shared_ptr<Contours> extractCompactContours(MorphologicalTreePtr tree){ return extractCompactContours(tree.get()); }

	static std::vector<std::unordered_set<int>> extractNonCompactContours(MorphologicalTree* tree);
    static std::vector<std::unordered_set<int>> extractNonCompactContours(MorphologicalTreePtr tree){ return extractNonCompactContours(tree.get()); }

    static std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>> computeAttributesByComputer(MorphologicalTree* tree, std::shared_ptr<AttributeComputer> comp, const DependencyMap& available = {});
    static std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>> computeAttributesByComputer(MorphologicalTreePtr tree, std::shared_ptr<AttributeComputer> comp, const DependencyMap& available = {}){ return computeAttributesByComputer(tree.get(), comp, available); }
	
	static std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>> computeSingleAttribute(MorphologicalTree* tree, AttributeOrGroup attr, const DependencyMap& availableDeps = {});
    static std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>> computeSingleAttribute(MorphologicalTreePtr tree, AttributeOrGroup attr, const DependencyMap& availableDeps = {}){ return computeSingleAttribute(tree.get(), attr, availableDeps); }
	
	static std::pair<std::shared_ptr<AttributeNamesWithDelta>, std::shared_ptr<float[]>> computeSingleAttributeWithDelta(MorphologicalTree* tree, Attribute attribute, int delta, std::string padding="last-padding", const DependencyMap& availableDeps={});
    static std::pair<std::shared_ptr<AttributeNamesWithDelta>, std::shared_ptr<float[]>> computeSingleAttributeWithDelta(MorphologicalTreePtr tree, Attribute attribute, int delta, std::string padding="last-padding", const DependencyMap& availableDeps={}){ return computeSingleAttributeWithDelta(tree.get(), attribute, delta, padding, availableDeps); }

	static std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>> computeAttributes(MorphologicalTree* tree, const std::vector<AttributeOrGroup>& attributes,const DependencyMap& providedDependencies={});
	static std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>> computeAttributes(MorphologicalTreePtr tree, const std::vector<AttributeOrGroup>& attributes,const DependencyMap& providedDependencies={}){ return computeAttributes(tree.get(), attributes, providedDependencies); }

	static ImageFloatPtr computerAttributeMapping(MorphologicalTree* tree, Attribute attribute);
    static ImageFloatPtr computerAttributeMapping(MorphologicalTreePtr tree, Attribute attribute){ return computerAttributeMapping(tree.get(), attribute); }
    
};









/*
Computacao incremental de countours
*/
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


