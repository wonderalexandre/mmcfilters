

#ifndef ATTRIBUTE_COMPUTED_INCREMENTALLY_H
#define ATTRIBUTE_COMPUTED_INCREMENTALLY_H

#include "../include/NodeMT.hpp"
#include "../include/ImageUtils.hpp"
#include "../include/MorphologicalTree.hpp"
#include "../include/Common.hpp"
#include <iterator>  
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits> // Para usar std::numeric_limits<float>::epsilon()
#include <unordered_map>
#include <utility>
#include <array>
#include <variant>

#define PI 3.14159265358979323846


enum class InvalidAttribute{
	INVALID
}; 

enum class GeometricAttribute {
    AREA, 
    VOLUME,
    LEVEL,
	DYNAMICS,
    MEAN_LEVEL,
    VARIANCE_LEVEL,
	BOUNDING_BOX,	
    BOX_WIDTH,
    BOX_HEIGHT,
    RECTANGULARITY,
    RATIO_WH,
	CENTRAL_MOMENT,
    CENTRAL_MOMENT_00,
	CENTRAL_MOMENT_20,
    CENTRAL_MOMENT_02,
    CENTRAL_MOMENT_11,
    CENTRAL_MOMENT_30,
    CENTRAL_MOMENT_03,
    CENTRAL_MOMENT_21,
    CENTRAL_MOMENT_12,
    AXIS_ORIENTATION, 
    LENGTH_MAJOR_AXIS,
    LENGTH_MINOR_AXIS,
    ECCENTRICITY,
    COMPACTNESS,
	HU_MOMENT,
	INERTIA,
    HU_MOMENT_1,
    HU_MOMENT_2,
    HU_MOMENT_3,
    HU_MOMENT_4,
    HU_MOMENT_5,
    HU_MOMENT_6,
    HU_MOMENT_7
};


enum class StructuralAttribute {
    HEIGHT,
    DEPTH,
    IS_LEAF,
    IS_ROOT,
    NUM_CHILDREN,
    NUM_SIBLINGS,
    NUM_DESCENDANTS,
    NUM_LEAF_DESCENDANTS,
    LEAF_RATIO,
    BALANCE,
    AVG_CHILD_HEIGHT
};



using Attribute = std::variant<GeometricAttribute, StructuralAttribute, InvalidAttribute>;

namespace std {
    template <>
    struct hash<StructuralAttribute> {
        size_t operator()(const StructuralAttribute& attr) const {
            return static_cast<size_t>(attr);
        }
    };

    template <>
    struct hash<GeometricAttribute> {
        size_t operator()(const GeometricAttribute& attr) const {
            return static_cast<size_t>(attr);
        }
    };

    template <>
    struct hash<Attribute> {
        size_t operator()(const Attribute& attr) const {
            return std::visit([](auto&& a) -> size_t {
                return std::hash<std::decay_t<decltype(a)>>{}(a);
            }, attr);
        }
    };
}


class AttributeNames {
	public:
		std::unordered_map<Attribute, int> indexMap;
		const int NUM_ATTRIBUTES;
	
		AttributeNames(std::unordered_map<Attribute, int>&& map) : indexMap(std::move(map)), NUM_ATTRIBUTES(indexMap.size()) {}
	
		static AttributeNames geometric(int n) {
			std::unordered_map<Attribute, int> map;
			int i = 0;
			for (int attr = static_cast<int>(GeometricAttribute::AREA);
				 attr <= static_cast<int>(GeometricAttribute::HU_MOMENT_7); ++attr) {
				map[static_cast<GeometricAttribute>(attr)] = i++ * n;
			}
			return AttributeNames(std::move(map));
		}
	
		static AttributeNames structural(int n) {
			std::unordered_map<Attribute, int> map;
			int i = 0;
			for (int attr = static_cast<int>(StructuralAttribute::HEIGHT);
				 attr <= static_cast<int>(StructuralAttribute::AVG_CHILD_HEIGHT); ++attr) {
				map[static_cast<StructuralAttribute>(attr)] = i++ * n;
			}
			return AttributeNames(std::move(map));
		}
	
		int getIndex(const Attribute& attr) const {
			return indexMap.at(attr);
		}
	
		int linearIndex(int nodeIndex, const Attribute& attr) const {
			return nodeIndex * NUM_ATTRIBUTES + getIndex(attr);
		}

		static std::string toString(const Attribute& attr) {
			return std::visit([](auto&& a) -> std::string {
				return toString(a);
			}, attr);
		}

		static std::string toString(GeometricAttribute attr) {
			switch (attr) {
				case GeometricAttribute::AREA: return "AREA";
				case GeometricAttribute::VOLUME: return "VOLUME";
				case GeometricAttribute::LEVEL: return "LEVEL";
				case GeometricAttribute::DYNAMICS: return "DYNAMICS";
				case GeometricAttribute::MEAN_LEVEL: return "MEAN_LEVEL";
				case GeometricAttribute::VARIANCE_LEVEL: return "VARIANCE_LEVEL";
				//case GeometricAttribute::STANDARD_DEVIATION: return "STANDARD_DEVIATION";
				case GeometricAttribute::BOX_WIDTH: return "BOX_WIDTH";
				case GeometricAttribute::BOX_HEIGHT: return "BOX_HEIGHT";
				case GeometricAttribute::RECTANGULARITY: return "RECTANGULARITY";
				case GeometricAttribute::RATIO_WH: return "RATIO_WH";
				case GeometricAttribute::CENTRAL_MOMENT_20: return "CENTRAL_MOMENT_20";
				case GeometricAttribute::CENTRAL_MOMENT_02: return "CENTRAL_MOMENT_02";
				case GeometricAttribute::CENTRAL_MOMENT_11: return "CENTRAL_MOMENT_11";
				case GeometricAttribute::CENTRAL_MOMENT_30: return "CENTRAL_MOMENT_30";
				case GeometricAttribute::CENTRAL_MOMENT_03: return "CENTRAL_MOMENT_03";
				case GeometricAttribute::CENTRAL_MOMENT_21: return "CENTRAL_MOMENT_21";
				case GeometricAttribute::CENTRAL_MOMENT_12: return "CENTRAL_MOMENT_12";
				case GeometricAttribute::AXIS_ORIENTATION: return "AXIS_ORIENTATION";
				case GeometricAttribute::LENGTH_MAJOR_AXIS: return "LENGTH_MAJOR_AXIS";
				case GeometricAttribute::LENGTH_MINOR_AXIS: return "LENGTH_MINOR_AXIS";
				case GeometricAttribute::ECCENTRICITY: return "ECCENTRICITY";
				case GeometricAttribute::COMPACTNESS: return "COMPACTNESS";
				case GeometricAttribute::INERTIA: return "INERTIA";
				case GeometricAttribute::HU_MOMENT_1: return "HU_MOMENT_1";
				case GeometricAttribute::HU_MOMENT_2: return "HU_MOMENT_2";
				case GeometricAttribute::HU_MOMENT_3: return "HU_MOMENT_3";
				case GeometricAttribute::HU_MOMENT_4: return "HU_MOMENT_4";
				case GeometricAttribute::HU_MOMENT_5: return "HU_MOMENT_5";
				case GeometricAttribute::HU_MOMENT_6: return "HU_MOMENT_6";
				case GeometricAttribute::HU_MOMENT_7: return "HU_MOMENT_7";
				default: return "UNKNOWN";
			}
		}


		static std::string toString(StructuralAttribute attr) {
			switch (attr) {
				case StructuralAttribute::HEIGHT: return "HEIGHT";
				case StructuralAttribute::DEPTH: return "DEPTH";
				case StructuralAttribute::IS_LEAF: return "IS_LEAF";
				case StructuralAttribute::IS_ROOT: return "IS_ROOT";
				case StructuralAttribute::NUM_CHILDREN: return "NUM_CHILDREN";
				case StructuralAttribute::NUM_SIBLINGS: return "NUM_SIBLINGS";
				case StructuralAttribute::NUM_DESCENDANTS: return "NUM_DESCENDANTS";
				case StructuralAttribute::NUM_LEAF_DESCENDANTS: return "NUM_LEAF_DESCENDANTS";
				case StructuralAttribute::LEAF_RATIO: return "LEAF_RATIO";
				case StructuralAttribute::BALANCE: return "BALANCE";
				case StructuralAttribute::AVG_CHILD_HEIGHT: return "AVG_CHILD_HEIGHT";
				default: return "UNKNOWN";
			}
		}

			
		static Attribute parseAttribute(const std::string& str) {
			if (auto g = parseGeometricAttribute(str)) {
				return *g;
			}
			if (auto s = parseStructuralAttribute(str)) {
				return *s;
			}
			return InvalidAttribute::INVALID;

		}
			

		static std::optional<GeometricAttribute> parseGeometricAttribute(const std::string& str) {
			if (str == "AREA") return GeometricAttribute::AREA;
			else if (str == "VOLUME") return GeometricAttribute::VOLUME;
			else if (str == "LEVEL") return GeometricAttribute::LEVEL;
			else if (str == "DYNAMICS") return GeometricAttribute::DYNAMICS;
			else if (str == "MEAN_LEVEL") return GeometricAttribute::MEAN_LEVEL;
			else if (str == "VARIANCE_LEVEL") return GeometricAttribute::VARIANCE_LEVEL;
			//else if (str == "STANDARD_DEVIATION") return GeometricAttribute::STANDARD_DEVIATION;
			else if (str == "BOX_WIDTH") return GeometricAttribute::BOX_WIDTH;
			else if (str == "BOX_HEIGHT") return GeometricAttribute::BOX_HEIGHT;
			else if (str == "RECTANGULARITY") return GeometricAttribute::RECTANGULARITY;
			else if (str == "RATIO_WH") return GeometricAttribute::RATIO_WH;
			else if (str == "CENTRAL_MOMENT_20") return GeometricAttribute::CENTRAL_MOMENT_20;
			else if (str == "CENTRAL_MOMENT_02") return GeometricAttribute::CENTRAL_MOMENT_02;
			else if (str == "CENTRAL_MOMENT_11") return GeometricAttribute::CENTRAL_MOMENT_11;
			else if (str == "CENTRAL_MOMENT_30") return GeometricAttribute::CENTRAL_MOMENT_30;
			else if (str == "CENTRAL_MOMENT_03") return GeometricAttribute::CENTRAL_MOMENT_03;
			else if (str == "CENTRAL_MOMENT_21") return GeometricAttribute::CENTRAL_MOMENT_21;
			else if (str == "CENTRAL_MOMENT_12") return GeometricAttribute::CENTRAL_MOMENT_12;
			else if (str == "AXIS_ORIENTATION") return GeometricAttribute::AXIS_ORIENTATION;
			else if (str == "LENGTH_MAJOR_AXIS") return GeometricAttribute::LENGTH_MAJOR_AXIS;
			else if (str == "LENGTH_MINOR_AXIS") return GeometricAttribute::LENGTH_MINOR_AXIS;
			else if (str == "ECCENTRICITY") return GeometricAttribute::ECCENTRICITY;
			else if (str == "COMPACTNESS") return GeometricAttribute::COMPACTNESS;
			else if (str == "INERTIA") return GeometricAttribute::INERTIA;
			else if (str == "HU_MOMENT_1") return GeometricAttribute::HU_MOMENT_1;
			else if (str == "HU_MOMENT_2") return GeometricAttribute::HU_MOMENT_2;
			else if (str == "HU_MOMENT_3") return GeometricAttribute::HU_MOMENT_3;
			else if (str == "HU_MOMENT_4") return GeometricAttribute::HU_MOMENT_4;
			else if (str == "HU_MOMENT_5") return GeometricAttribute::HU_MOMENT_5;
			else if (str == "HU_MOMENT_6") return GeometricAttribute::HU_MOMENT_6;
			else if (str == "HU_MOMENT_7") return GeometricAttribute::HU_MOMENT_7;
			return std::nullopt; // inválido
		}

		static std::optional<StructuralAttribute> parseStructuralAttribute(const std::string& str) {
			if (str == "HEIGHT") return StructuralAttribute::HEIGHT;
			else if (str == "DEPTH") return StructuralAttribute::DEPTH;
			else if (str == "IS_LEAF") return StructuralAttribute::IS_LEAF;
			else if (str == "IS_ROOT") return StructuralAttribute::IS_ROOT;
			else if (str == "NUM_CHILDREN") return StructuralAttribute::NUM_CHILDREN;
			else if (str == "NUM_SIBLINGS") return StructuralAttribute::NUM_SIBLINGS;
			else if (str == "NUM_DESCENDANTS") return StructuralAttribute::NUM_DESCENDANTS;
			else if (str == "NUM_LEAF_DESCENDANTS") return StructuralAttribute::NUM_LEAF_DESCENDANTS;
			else if (str == "LEAF_RATIO") return StructuralAttribute::LEAF_RATIO;
			else if (str == "BALANCE") return StructuralAttribute::BALANCE;
			else if (str == "AVG_CHILD_HEIGHT") return StructuralAttribute::AVG_CHILD_HEIGHT;
			return std::nullopt; // inválido
		}
};


class AttributeComputer {
	public:
		virtual ~AttributeComputer() = default;
	
		/// Executa a computação dos atributos produzidos por essa classe
		virtual void compute(
			MorphologicalTreePtr tree,
			std::shared_ptr<float[]> buffer,
			std::shared_ptr<AttributeNames> attrNames,
			const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources = {}
		) const = 0;
	
		/// Atributos produzidos
		virtual std::vector<Attribute> attributes() const = 0;
	
		/// Atributos requeridos para o cálculo (apenas metadado)
		virtual std::vector<Attribute> requiredAttributes() const { return {}; }
	
};

using DependencyMap = std::unordered_map<Attribute, std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>>>;


class AttributeComputedIncrementally{

public:

 
    virtual void preProcessing(NodeMTPtr v);

    virtual void mergeChildren(NodeMTPtr parent, NodeMTPtr child);

    virtual void postProcessing(NodeMTPtr parent);

    void computerAttribute(NodeMTPtr root);

	static void computerAttribute(NodeMTPtr root, 
										std::function<void(NodeMTPtr)> preProcessing,
										std::function<void(NodeMTPtr, NodeMTPtr)> mergeChildren,
										std::function<void(NodeMTPtr)> postProcessing ){
		
		preProcessing(root);
			
		for(NodeMTPtr child: root->getChildren()){
			AttributeComputedIncrementally::computerAttribute(child, preProcessing, mergeChildren, postProcessing);
			mergeChildren(root, child);
		}

		postProcessing(root);
	}


	class ContoursMT{
		private:
			std::vector<std::vector<int>> contours;
			std::vector<std::vector<int>> contoursToRemove;

		public:
			ContoursMT(int numNodes): contours(numNodes), contoursToRemove(numNodes){}

			void add(NodeMTPtr node, int pixel){
				contours[node->getIndex()].push_back(pixel);
			}
			void remove(NodeMTPtr node, int pixel){
				contoursToRemove[node->getIndex()].push_back(pixel);
			}

			std::vector<std::vector<int>>& getCompactContours() {
				return contours;
			}

			std::unordered_set<int> getContour(NodeMTPtr node) {
				std::unordered_set<int> contour;
				AttributeComputedIncrementally::computerAttribute(node,
					[](NodeMTPtr node) -> void {},  // pre-processing
					[](NodeMTPtr parent, NodeMTPtr child) -> void { }, // merge-processing
					[&contour, this](NodeMTPtr node) -> void { //post-processing
						for(int p: this->contours[node->getIndex()]){
							contour.insert(p);
						}
						for(int p: this->contoursToRemove[node->getIndex()]){
							contour.erase(p);
						}
					}
				);
				return contour;
			}

			void visitContours(MorphologicalTreePtr tree, std::function<void(NodeMTPtr, const std::unordered_set<int>&)> visitor) {
				const int numNodes = tree->getNumNodes();
			
				std::vector<std::unique_ptr<std::unordered_set<int>>> contoursByNodes(numNodes);
			
				AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
					[](NodeMTPtr) -> void {},
			
					// merge: funde filhos no pai, usando o maior conjunto como base
					[&contoursByNodes](NodeMTPtr parent, NodeMTPtr child) -> void {
						auto& parentContour = contoursByNodes[parent->getIndex()];
						auto& childContour = contoursByNodes[child->getIndex()];
			
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
					[this, &contoursByNodes, &visitor](NodeMTPtr node) -> void {
						auto& contour = contoursByNodes[node->getIndex()];
						if (!contour) {
							contour = std::make_unique<std::unordered_set<int>>();
						}
			
						for (int p : this->contours[node->getIndex()]) {
							contour->insert(p);
						}
			
						for (int p : this->contoursToRemove[node->getIndex()]) {
							contour->erase(p);
						}
			
						visitor(node, *contour);
					}
				);
			}

			void visitContoursAndCCs(MorphologicalTreePtr tree, std::function<void(NodeMTPtr, const std::list<int>&, const std::unordered_set<int>&)> visitor) {
				const int numNodes = tree->getNumNodes();
			
				std::vector<std::unique_ptr<std::unordered_set<int>>> contoursByNodes(numNodes);
				std::vector<std::unique_ptr<std::list<int>>> CCsByNodes(numNodes);
			
				AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
					[](NodeMTPtr) -> void {},
			
					[&CCsByNodes, &contoursByNodes](NodeMTPtr parent, NodeMTPtr child) -> void {
						// --- Contornos ---
						auto& parentContour = contoursByNodes[parent->getIndex()];
						auto& childContour = contoursByNodes[child->getIndex()];
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
						auto& parentCC = CCsByNodes[parent->getIndex()];
						auto& childCC = CCsByNodes[child->getIndex()];
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
					[this, &contoursByNodes, &CCsByNodes, &visitor](NodeMTPtr node) -> void {
						// --- Contornos ---
						auto& contour = contoursByNodes[node->getIndex()];
						if (!contour) {
							contour = std::make_unique<std::unordered_set<int>>();
						}
						for (int p : this->contours[node->getIndex()]) {
							contour->insert(p);
						}
						for (int p : this->contoursToRemove[node->getIndex()]) {
							contour->erase(p);
						}
						
						// --- Componentes Conexos ---
						auto& cc = CCsByNodes[node->getIndex()];
						if (!cc) {
							cc = std::make_unique<std::list<int>>();
						}
						cc->insert(cc->end(), node->getCNPs().begin(), node->getCNPs().end());
			
						
						visitor(node, *cc, *contour);
					}
				);
			}
	};


	static ContoursMT extractCompactCountors(MorphologicalTreePtr tree){
		ContoursMT contoursMT(tree->getNumNodes());

		std::vector<std::vector<int>> contoursToRemoveLCA(tree->getNumNodes());
		std::vector<std::int8_t> ncount(tree->getNumRowsOfImage() * tree->getNumColsOfImage(), 0);
		AdjacencyRelationPtr adj4 = std::make_shared<AdjacencyRelation>(tree->getNumRowsOfImage(), tree->getNumColsOfImage(), 1);
		LCAEulerRMQ lca(tree);	

		AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
			[](NodeMTPtr node) -> void { // pre-processing

			},
			[](NodeMTPtr parent, NodeMTPtr child) -> void { // merge-processing
				
			},
			[&contoursMT, &contoursToRemoveLCA, &lca, &ncount, tree, adj4](NodeMTPtr nodeP) -> void { // post-processing
				std::vector<int> &NcontourToRemoveLCA = contoursToRemoveLCA[nodeP->getIndex()];

				NodeMTPtr nodeLCA = nodeP;
				for(int p: NcontourToRemoveLCA){ //pixels que sao contornos de nodes descendentes ao NodeAtual
					bool isPixelToBeRemoved = true;
					for (int r : adj4->getNeighboringPixels(p)) { //Existe um nodeQ ascendente de NodeAtual contendo p como contorno? (p, q) in A
						NodeMTPtr nodeR = tree->getSC(r); 
						if (tree->isStrictAncestor(nodeR, nodeLCA)){
							contoursToRemoveLCA[nodeR->getIndex()].push_back(p); 
							isPixelToBeRemoved = false;	
					  	}else if(!tree->isComparable(nodeLCA, nodeR)) {
							NodeMTPtr otherNodeLCA = lca.findLowestCommonAncestor(nodeLCA, nodeR);
							contoursToRemoveLCA[otherNodeLCA->getIndex()].push_back(p);
							isPixelToBeRemoved = false;
						}
					}
					if(!adj4->isBorderDomainImage(p) && isPixelToBeRemoved){
						contoursMT.remove(nodeLCA, p);
					}
				}
			
				for (int p : nodeP->getCNPs()) {
					if (adj4->isBorderDomainImage(p)){
						ncount[p]++;
					}

					for (int q : adj4->getNeighboringPixels(p)) {
						NodeMTPtr nodeQ = tree->getSC(q); 
						if(!tree->isComparable(nodeP, nodeQ)){ //se os nodeP e nodeQ não sao comparaveis, então p pode ser removido pelo LCA de nodeP e nodeQ 
							NodeMTPtr nodeLCA = lca.findLowestCommonAncestor(nodeP, nodeQ);
							contoursToRemoveLCA[nodeLCA->getIndex()].push_back(p);
							ncount[p]++;
						}
						else if(tree->isStrictDescendant(nodeP, nodeQ)){  //maxtree:  SC(p) \subset SC(q) <=> f(p) > f(q)
					  		ncount[p]++;
						}else if (tree->isStrictAncestor(nodeP, nodeQ)) { ////maxtree:  SC(q) \subset SC(p) <=> f(p) < f(q)
					  		ncount[q]--;
							if (ncount[q] == 0) {
								contoursMT.remove(nodeP, q);
							}
						}
				  	}
				  	if (ncount[p] > 0){
						contoursMT.add(nodeP, p);
					}
				}

			}
		);
				  
		return contoursMT;
	}

	static std::vector<std::unordered_set<int>> extractCountors(MorphologicalTreePtr tree){
		std::vector<std::unordered_set<int>> contours(tree->getNumNodes());
		std::vector<std::vector<int>> contoursToRemoveLCA(tree->getNumNodes());
		std::vector<std::int8_t> ncount(tree->getNumRowsOfImage() * tree->getNumColsOfImage(), 0);
		AdjacencyRelationPtr adj4 = std::make_shared<AdjacencyRelation>(tree->getNumRowsOfImage(), tree->getNumColsOfImage(), 1);
		LCAEulerRMQ lca(tree);	

		AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
			[](NodeMTPtr node) -> void { // pre-processing

			},
			[&contours, &ncount, tree, adj4](NodeMTPtr parent, NodeMTPtr child) -> void { // merge-processing
				std::unordered_set<int> &Ncontour = contours[parent->getIndex()];
				for (int p : contours[child->getIndex()]){
					Ncontour.insert(p);
				}
			},
			[&contours, &contoursToRemoveLCA, &lca, &ncount, tree, adj4](NodeMTPtr nodeP) -> void { // post-processing
				// Initialise contours of node "N"
				std::unordered_set<int> &Ncontour = contours[nodeP->getIndex()];
				std::vector<int> &NcontourToRemoveLCA = contoursToRemoveLCA[nodeP->getIndex()];
				NodeMTPtr nodeLCA = nodeP;
				for(int p: NcontourToRemoveLCA){ //pixels que sao contornos de nodes descendentes ao NodeAtual
					bool isPixelToBeRemoved = true;
					
					for (int r : adj4->getNeighboringPixels(p)) { //Existe um nodeQ ascendente de NodeAtual contendo p como contorno? (p, q) in A
						NodeMTPtr nodeR = tree->getSC(r); 
						if (tree->isStrictAncestor(nodeR, nodeLCA)){
							contoursToRemoveLCA[nodeR->getIndex()].push_back(p); 
							isPixelToBeRemoved = false;	
					  	}else if(!tree->isComparable(nodeLCA, nodeR)) {
							NodeMTPtr otherNodeLCA = lca.findLowestCommonAncestor(nodeLCA, nodeR);
							contoursToRemoveLCA[otherNodeLCA->getIndex()].push_back(p);
							isPixelToBeRemoved = false;
						}
					}
					if(!adj4->isBorderDomainImage(p) && isPixelToBeRemoved){
						Ncontour.erase(p);
					}
				}
			
				for (int p : nodeP->getCNPs()) {
					if (adj4->isBorderDomainImage(p)){
						ncount[p]++;
					}
					for (int q : adj4->getNeighboringPixels(p)) {
						NodeMTPtr nodeQ = tree->getSC(q); 
						if(!tree->isComparable(nodeP, nodeQ)){ //se os nodeP e nodeQ não sao comparaveis, então p pode ser removido pelo LCA de nodeP e nodeQ 
							NodeMTPtr nodeLCA = lca.findLowestCommonAncestor(nodeP, nodeQ);
							contoursToRemoveLCA[nodeLCA->getIndex()].push_back(p);
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



	static float* computerStructTreeAttributes(MorphologicalTreePtr tree){
		const int n = tree->getNumNodes();
		AttributeNames attrNames = AttributeNames::structural(n);
		float* attrs = new float[n * attrNames.NUM_ATTRIBUTES];

		using enum StructuralAttribute;
		AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
			[&](NodeMTPtr node) {
				int idx = node->getIndex();
				int parentDepth = node->getParent() ? attrs[attrNames.linearIndex(node->getParent()->getIndex(), DEPTH)] : 0;
	
				attrs[attrNames.linearIndex(idx, HEIGHT)] = 0.0f; // altura
				attrs[attrNames.linearIndex(idx, DEPTH)] = node->getParent() ? parentDepth + 1 : 0; // profundidade
				attrs[attrNames.linearIndex(idx, IS_LEAF)] = node->getChildren().empty() ? 1.0f : 0.0f; // é folha
				attrs[attrNames.linearIndex(idx, IS_ROOT)] = node->getParent() == nullptr ? 1.0f : 0.0f; // é raiz
				attrs[attrNames.linearIndex(idx, NUM_CHILDREN)] = node->getChildren().size();
				attrs[attrNames.linearIndex(idx, NUM_SIBLINGS)] = node->getParent() ? node->getParent()->getChildren().size() - 1 : 0;
				attrs[attrNames.linearIndex(idx, NUM_DESCENDANTS)] = 0.0f;
				attrs[attrNames.linearIndex(idx, NUM_LEAF_DESCENDANTS)] = attrs[attrNames.linearIndex(idx, IS_LEAF)];
				attrs[attrNames.linearIndex(idx, LEAF_RATIO)] = 0.0f;
				attrs[attrNames.linearIndex(idx, BALANCE)] = 0.0f;
				attrs[attrNames.linearIndex(idx, AVG_CHILD_HEIGHT)] = 0.0f;
			},
			[&](NodeMTPtr parent, NodeMTPtr child) {
				int pIdx = parent->getIndex();
				int cIdx = child->getIndex();

				// Acumulando descendentes
				attrs[attrNames.linearIndex(pIdx, NUM_DESCENDANTS)] += attrs[attrNames.linearIndex(cIdx, NUM_DESCENDANTS)] + 1;
				attrs[attrNames.linearIndex(pIdx, NUM_LEAF_DESCENDANTS)] += attrs[attrNames.linearIndex(cIdx, NUM_LEAF_DESCENDANTS)];

				// Altura
				float childHeight = attrs[attrNames.linearIndex(cIdx, HEIGHT)];
				float& parentHeight = attrs[attrNames.linearIndex(pIdx, HEIGHT)];
				parentHeight = std::max(parentHeight, childHeight + 1);

				// Balanceamento
				float& minH = attrs[attrNames.linearIndex(pIdx, BALANCE)]; // usado como mínimo temporário
				float& sumH = attrs[attrNames.linearIndex(pIdx, AVG_CHILD_HEIGHT)];
				int numChildren = parent->getChildren().size();

				if (numChildren == 1) {
					minH = childHeight;
					sumH = childHeight;
				} else {
					minH = std::min(minH, childHeight);
					sumH += childHeight;
				}

			},
			[&](NodeMTPtr node) {
				int idx = node->getIndex();
				float desc = attrs[attrNames.linearIndex(idx, NUM_DESCENDANTS)];
				float folhas = attrs[attrNames.linearIndex(idx, NUM_LEAF_DESCENDANTS)];
	
				// Razão folhas/descendentes
				attrs[attrNames.linearIndex(idx, LEAF_RATIO)] = desc > 0.0f ? folhas / (desc + 1.0f) : 1.0f;
	
				// Balanceamento e média
				if (!node->getChildren().empty()) {
					float alturaMax = attrs[attrNames.linearIndex(idx, HEIGHT)];
					float alturaMin = attrs[attrNames.linearIndex(idx, BALANCE)];
					attrs[attrNames.linearIndex(idx, BALANCE)] = alturaMax - alturaMin;
	
					attrs[attrNames.linearIndex(idx, AVG_CHILD_HEIGHT)] /= node->getChildren().size();
				}

			}
		);
	
		return attrs;
	}


	static float* computerAttributeByIndex(MorphologicalTreePtr tree, Attribute attribute){
		const int n = tree->getNumNodes();
		float* attr = new float[n];
		auto [attrNames, attrs] = AttributeComputedIncrementally::computerBasicAttributes(tree);
		for(int idx = 0; idx < n; idx++){
			attr[idx] = attrs[ attrNames.linearIndex(idx, attribute) ];
		}
		delete[] attrs;

		return attr;
	}

	static std::pair<AttributeNames, float*> computerBasicAttributes(MorphologicalTreePtr tree){
	    
		int n = tree->getNumNodes();
		AttributeNames attrNames = AttributeNames::geometric(n);
		float* attrs = new float[n * attrNames.NUM_ATTRIBUTES];

		std::unique_ptr<int[]> xmax(new int[n]);
		std::unique_ptr<int[]> ymax(new int[n]);
		std::unique_ptr<int[]> xmin(new int[n]);
		std::unique_ptr<int[]> ymin(new int[n]);
		
		//momentos geometricos para calcular o centroide
		std::unique_ptr<long int[]> sumX(new long int[n]);//sum x
		std::unique_ptr<long int[]> sumY(new long int[n]);//sum y
		

		std::unique_ptr<long[]> sumGrayLevelSquare(new long[n]);
		int numCols = tree->getNumColsOfImage();
		int numRows = tree->getNumRowsOfImage();
		
		auto indexOf = [&](int idxNode, GeometricAttribute attribute) {
			return attrNames.linearIndex(idxNode, attribute);
		};
		using enum GeometricAttribute;

		//computação dos atributos: area, volume, gray level, mean of gray level, variance of gray level, standard deviation gray level, Box width, Box height, rectangularity, ratio (Box width, Box height) e momentos geometricos 
	    AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
						[&indexOf, &attrs, n,  &xmax, &ymax, &xmin, &ymin, &sumX, &sumY, &sumGrayLevelSquare, numCols, numRows](NodeMTPtr node) -> void {
							
							attrs[ indexOf(node->getIndex(), AREA) ] = node->getCNPs().size(); //area
							attrs[ indexOf(node->getIndex(), VOLUME) ] = node->getCNPs().size() * node->getLevel(); //volume =>  \sum{ f }
							attrs[ indexOf(node->getIndex(), LEVEL) ] = node->getLevel(); //level
							if (node->isMaxtreeNode()) {
								attrs[ indexOf(node->getIndex(), DYNAMICS)] = std::numeric_limits<int>::lowest(); // Procuramos máximo para max-tree
							} else {
								attrs[ indexOf(node->getIndex(), DYNAMICS)] = std::numeric_limits<int>::max(); // Procuramos mínimo para min-tree
							}

							xmax[node->getIndex()] = 0;
							ymax[node->getIndex()] = 0;
							xmin[node->getIndex()] = numCols;
							ymin[node->getIndex()] = numRows;
							sumX[node->getIndex()] = 0;
							sumY[node->getIndex()] = 0;
							sumGrayLevelSquare[node->getIndex()] = std::pow(node->getLevel(), 2) * node->getCNPs().size(); //computando: \sum{ f^2 }
							for(int p: node->getCNPs()) {
								int x = p % numCols;
								int y = p / numCols;
								xmin[node->getIndex()] = std::min(xmin[node->getIndex()], x);
								ymin[node->getIndex()] = std::min(ymin[node->getIndex()], y);
								xmax[node->getIndex()] = std::max(xmax[node->getIndex()], x);
								ymax[node->getIndex()] = std::max(ymax[node->getIndex()], y);

								sumX[node->getIndex()] += x;
								sumY[node->getIndex()] += y;
							}
						},
						[&indexOf, &attrNames, &attrs, n, &xmax, &ymax, &xmin, &ymin, &sumX, &sumY, &sumGrayLevelSquare](NodeMTPtr parent, NodeMTPtr child) -> void {
							attrs[ indexOf(parent->getIndex(), AREA) ] += attrs[ indexOf(child->getIndex(), AREA) ]; //area
							attrs[ indexOf(parent->getIndex(), VOLUME) ] += attrs[ indexOf(child->getIndex(), VOLUME) ]; //volume
							if (parent->isMaxtreeNode()) {
								attrs[ indexOf(parent->getIndex(), DYNAMICS) ] = std::max(attrs[ indexOf(parent->getIndex(), DYNAMICS)], attrs[ indexOf(child->getIndex(), DYNAMICS)] );
							} else {
								attrs[ indexOf(parent->getIndex(), DYNAMICS) ] = std::min(attrs[ indexOf(parent->getIndex(), DYNAMICS)], attrs[ indexOf(child->getIndex(), DYNAMICS)] );
							}

							sumGrayLevelSquare[parent->getIndex()] += sumGrayLevelSquare[child->getIndex()]; //computando: \sum{ f^2 }

							ymax[parent->getIndex()] = std::max(ymax[parent->getIndex()], ymax[child->getIndex()]);
							xmax[parent->getIndex()] = std::max(xmax[parent->getIndex()], xmax[child->getIndex()]);
							ymin[parent->getIndex()] = std::min(ymin[parent->getIndex()], ymin[child->getIndex()]);
							xmin[parent->getIndex()] = std::min(xmin[parent->getIndex()], xmin[child->getIndex()]);
		
							sumX[parent->getIndex()] += sumX[child->getIndex()];
							sumY[parent->getIndex()] += sumY[child->getIndex()];
							
						},
						[&indexOf, &attrNames, &attrs, n, &xmax, &ymax, &xmin, &ymin, &sumGrayLevelSquare](NodeMTPtr node) -> void {
							
							float area = attrs[ indexOf(node->getIndex(), AREA) ];
							float volume = attrs[ indexOf(node->getIndex(), VOLUME) ];
							float width = xmax[node->getIndex()] - xmin[node->getIndex()] + 1;	
							float height = ymax[node->getIndex()] - ymin[node->getIndex()] + 1;	
							
							float meanGrayLevel = volume / area; //mean graylevel - // E(f)
							double meanGrayLevelSquare = sumGrayLevelSquare[node->getIndex()] / area; // E(f^2)
							float var = meanGrayLevelSquare - (meanGrayLevel * meanGrayLevel); //variance: E(f^2) - E(f)^2
							attrs[ indexOf(node->getIndex(), VARIANCE_LEVEL) ] = var > 0? var : 0; //variance
							
							/*
							if (attrs[indexOf(node->getIndex(), VARIANCE_LEVEL)] >= 0.0f) {
								attrs[indexOf(node->getIndex(), STANDARD_DEVIATION)] = std::sqrt(attrs[indexOf(node->getIndex(), VARIANCE_LEVEL)]);
							} else {
								attrs[indexOf(node->getIndex(), STANDARD_DEVIATION)] = 0.0f; // Se a variância for negativa, definir desvio padrão como 0
							}*/
							
							attrs[indexOf(node->getIndex(), MEAN_LEVEL)] = meanGrayLevel;
							attrs[indexOf(node->getIndex(), BOX_WIDTH)] = width;
							attrs[indexOf(node->getIndex(), BOX_HEIGHT)] = height;
							attrs[indexOf(node->getIndex(), RECTANGULARITY)] = area / (width * height);
							attrs[indexOf(node->getIndex(), RATIO_WH)] = std::max(width, height) / std::min(width, height);
							
							if (node->getChildren().empty()) {
								attrs[ indexOf(node->getIndex(), DYNAMICS)] = 0; // Folhas têm dinâmica 0
							} else {
								attrs[ indexOf(node->getIndex(), DYNAMICS)] = std::abs( node->getLevel() - attrs[ indexOf(node->getIndex(), DYNAMICS)] ) + 1;
							}


		});

		

		//Computação dos momentos centrais e momentos de Hu
		AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
			[&indexOf, &attrs, n,  &sumX, &sumY, numCols](NodeMTPtr node) -> void {				
				// Inicialização dos momentos centrais
				attrs[indexOf(node->getIndex(), CENTRAL_MOMENT_20)] = 0.0f;
				attrs[indexOf(node->getIndex(), CENTRAL_MOMENT_02)] = 0.0f;
				attrs[indexOf(node->getIndex(), CENTRAL_MOMENT_11)] = 0.0f;
				attrs[indexOf(node->getIndex(), CENTRAL_MOMENT_30)] = 0.0f;
				attrs[indexOf(node->getIndex(), CENTRAL_MOMENT_03)] = 0.0f;
				attrs[indexOf(node->getIndex(), CENTRAL_MOMENT_21)] = 0.0f;
				attrs[indexOf(node->getIndex(), CENTRAL_MOMENT_12)] = 0.0f;

				// Cálculo do centroide
				float xCentroid = sumX[node->getIndex()] / attrs[indexOf(node->getIndex(), AREA)];
				float yCentroid = sumY[node->getIndex()] / attrs[indexOf(node->getIndex(), AREA)];

				for (int p : node->getCNPs()) {
					int x = p % numCols;
					int y = p / numCols;
					float dx = x - xCentroid;
					float dy = y - yCentroid;

					// Momentos centrais de segunda ordem
					attrs[indexOf(node->getIndex(), CENTRAL_MOMENT_20)] += std::pow(dx, 2);
					attrs[indexOf(node->getIndex(), CENTRAL_MOMENT_02)] += std::pow(dy, 2);
					attrs[indexOf(node->getIndex(), CENTRAL_MOMENT_11)] += dx * dy;

					// Momentos centrais de terceira ordem
					attrs[indexOf(node->getIndex(), CENTRAL_MOMENT_30)] += std::pow(dx, 3);
					attrs[indexOf(node->getIndex(), CENTRAL_MOMENT_03)] += std::pow(dy, 3);
					attrs[indexOf(node->getIndex(), CENTRAL_MOMENT_21)] += std::pow(dx, 2) * dy;
					attrs[indexOf(node->getIndex(), CENTRAL_MOMENT_12)] += dx * std::pow(dy, 2);
				}

			},
			[&indexOf, &attrs, n](NodeMTPtr parent, NodeMTPtr child) -> void {
				attrs[indexOf(parent->getIndex(), CENTRAL_MOMENT_20)] += attrs[indexOf(child->getIndex(), CENTRAL_MOMENT_20)];
				attrs[indexOf(parent->getIndex(), CENTRAL_MOMENT_02)] += attrs[indexOf(child->getIndex(), CENTRAL_MOMENT_02)];
				attrs[indexOf(parent->getIndex(), CENTRAL_MOMENT_11)] += attrs[indexOf(child->getIndex(), CENTRAL_MOMENT_11)];
				attrs[indexOf(parent->getIndex(), CENTRAL_MOMENT_30)] += attrs[indexOf(child->getIndex(), CENTRAL_MOMENT_30)];
				attrs[indexOf(parent->getIndex(), CENTRAL_MOMENT_03)] += attrs[indexOf(child->getIndex(), CENTRAL_MOMENT_03)];
				attrs[indexOf(parent->getIndex(), CENTRAL_MOMENT_21)] += attrs[indexOf(child->getIndex(), CENTRAL_MOMENT_21)];
				attrs[indexOf(parent->getIndex(), CENTRAL_MOMENT_12)] += attrs[indexOf(child->getIndex(), CENTRAL_MOMENT_12)];			
			},
			[&indexOf, &attrs, n](NodeMTPtr node) -> void {
				int idx = node->getIndex();
				float area = attrs[indexOf(idx, AREA)];
				auto normMoment = [area](float moment, int p, int q){ 
					return moment / std::pow( area, (p + q + 2.0) / 2.0); 
				}; //função para normalizacao dos momentos				
				

				//Momentos centrais
				float mu20 = attrs[indexOf(idx, CENTRAL_MOMENT_20)];
				float mu02 = attrs[indexOf(idx, CENTRAL_MOMENT_02)];
				float mu11 = attrs[indexOf(idx, CENTRAL_MOMENT_11)];
				float mu30 = attrs[indexOf(idx, CENTRAL_MOMENT_30)];
				float mu03 = attrs[indexOf(idx, CENTRAL_MOMENT_03)];
				float mu21 = attrs[indexOf(idx, CENTRAL_MOMENT_21)];
				float mu12 = attrs[indexOf(idx, CENTRAL_MOMENT_12)];
					
				float discriminant = std::pow(mu20 - mu02, 2.0f) + 4.0f * std::pow(mu11, 2.0f);

					
				// Verificar se o denominador é zero antes de calcular atan2 para evitar divisão por zero
				if (mu20 != mu02 || mu11 != 0) {
					float radians = 0.5 * std::atan2(2 * mu11, mu20 - mu02);// orientação em radianos
					float degrees = radians * (180.0 / M_PI); // Converter para graus
					if (degrees < 0) { // Ajustar para o intervalo [0, 360] graus
						degrees += 360.0;
					}
					attrs[indexOf(idx, AXIS_ORIENTATION)] = degrees; // Armazenar a orientação em graus no intervalo [0, 360]
				} else {
					attrs[indexOf(idx, AXIS_ORIENTATION)] = 0.0; // Se não for possível calcular a orientação, definir um valor padrão
				}

				// Verificar se o discriminante é positivo para evitar raiz quadrada de números negativos
				if (discriminant < 0) {
					std::cerr << "Erro: Discriminante negativo, ajustando para zero." << std::endl;
					discriminant = 0;
				}	
				float a1 = mu20 + mu02 + std::sqrt(discriminant); // autovalores (correspondente ao eixo maior)
				float a2 = mu20 + mu02 - std::sqrt(discriminant); // autovalores (correspondente ao eixo menor)

				// Verificar se a1 e a2 são positivos antes de calcular sqrt para evitar NaN
				if (a1 > 0) {
					attrs[indexOf(idx, LENGTH_MAJOR_AXIS)] = std::sqrt((2 * a1) / area); // length major axis
				} else {
					attrs[indexOf(idx, LENGTH_MAJOR_AXIS)] = 0.0; // Definir valor padrão
				}

				if (a2 > 0) {
					attrs[indexOf(idx, LENGTH_MINOR_AXIS)] = std::sqrt((2 * a2) / area); // length minor axis
				} else {
					attrs[indexOf(idx, LENGTH_MINOR_AXIS)] = 0.0; // Definir valor padrão
				}

				// Verificar se a2 é diferente de zero antes de calcular a excentricidade
				attrs[indexOf(idx, ECCENTRICITY)] = (std::abs(a2) > std::numeric_limits<float>::epsilon()) ? a1 / a2 : a1 / 0.1; // eccentricity

				// Verificar se mu20 + mu02 é diferente de zero antes de calcular a compacidade
				if ((mu20 + mu02) > std::numeric_limits<float>::epsilon()) {
					attrs[indexOf(idx, COMPACTNESS)]  = (1.0 / (2 * PI)) * (area / (mu20 + mu02)); // compactness
				} else {
					attrs[indexOf(idx, COMPACTNESS)]  = 0.0; // Definir valor padrão
				}


				// Calcular os momentos normalizados
				float eta20 = normMoment(mu20, 2, 0);
				float eta02 = normMoment(mu02, 0, 2);
				float eta11 = normMoment(mu11, 1, 1);
				float eta30 = normMoment(mu30, 3, 0);
				float eta03 = normMoment(mu03, 0, 3);
				float eta21 = normMoment(mu21, 2, 1);
				float eta12 = normMoment(mu12, 1, 2);

				// Cálculo dos momentos de Hu
				attrs[indexOf(idx, HU_MOMENT_1)] = eta20 + eta02; // primeiro momento de Hu => inertia
				attrs[indexOf(idx, HU_MOMENT_2)]  = std::pow(eta20 - eta02, 2) + 4 * std::pow(eta11, 2);
				attrs[indexOf(idx, HU_MOMENT_3)]  = std::pow(eta30 - 3 * eta12, 2) + std::pow(3 * eta21 - eta03, 2);
				attrs[indexOf(idx, HU_MOMENT_4)]  = std::pow(eta30 + eta12, 2) + std::pow(eta21 + eta03, 2);
				
				attrs[indexOf(idx, HU_MOMENT_5)] = 
					(eta30 - 3 * eta12) * (eta30 + eta12) * (std::pow(eta30 + eta12, 2) - 3 * std::pow(eta21 + eta03, 2)) +
					(3 * eta21 - eta03) * (eta21 + eta03) * (3 * std::pow(eta30 + eta12, 2) - std::pow(eta21 + eta03, 2));
				
					attrs[indexOf(idx, HU_MOMENT_6)] = 
					(eta20 - eta02) * (std::pow(eta30 + eta12, 2) - std::pow(eta21 + eta03, 2)) + 
					4 * eta11 * (eta30 + eta12) * (eta21 + eta03);
				
					attrs[indexOf(idx, HU_MOMENT_7)] = 
					(3 * eta21 - eta03) * (eta30 + eta12) * (std::pow(eta30 + eta12, 2) - 3 * std::pow(eta21 + eta03, 2)) -
					(eta30 - 3 * eta12) * (eta21 + eta03) * (3 * std::pow(eta30 + eta12, 2) - std::pow(eta21 + eta03, 2));

				
		});
		return std::make_pair(attrNames, attrs);
    }

	struct ExtinctionValues{
		NodeMTPtr leaf;
		NodeMTPtr cutoffNode;
		float extinction;
		ExtinctionValues(NodeMTPtr leaf, NodeMTPtr cutoffNode, float extinction)
			: leaf(leaf), cutoffNode(cutoffNode), extinction(extinction) {}
		
	};

	static std::vector<AttributeComputedIncrementally::ExtinctionValues> getExtinctionValue(MorphologicalTreePtr tree, float* attr){
		std::list<NodeMTPtr> leaves = tree->getLeaves();
		std::vector<AttributeComputedIncrementally::ExtinctionValues> leavesByExtinction;
		leavesByExtinction.reserve(leaves.size());
		std::unique_ptr<bool[]> visited(new bool[tree->getNumNodes()]()); //inicializa com false
		for(NodeMTPtr leaf: leaves){
			float extinction = 0;
			NodeMTPtr cutoffNode = leaf;
			NodeMTPtr parent = cutoffNode->getParent();
			bool flag = true;
			while (flag  &&  parent != nullptr) {
				if (parent->getChildren().size() > 1) {
					for(NodeMTPtr son: parent->getChildren()){  // verifica se possui irmao com atributo maior
						if(flag){
							if (visited[son->getIndex()]  &&  son != cutoffNode  &&  attr[son->getIndex()] == attr[cutoffNode->getIndex()]) { //EMPATE Grimaud,92
								flag = false;
							}
							else if (son != cutoffNode  &&  attr[son->getIndex()] > attr[cutoffNode->getIndex()]) {
								flag = false;
							}
							visited[son->getIndex()] = true;
						}
					}
				}
				if (flag) {
					cutoffNode = parent;
					parent = cutoffNode->getParent();
				}
			}
			if(parent != nullptr)
				extinction = attr[cutoffNode->getIndex()];
			leavesByExtinction.push_back( AttributeComputedIncrementally::ExtinctionValues(leaf, cutoffNode, extinction) );
			
		}
		return leavesByExtinction;
	}

	//using DependencyMap = std::unordered_map<Attribute, std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>>>;
	//static std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>> computeSingleAttribute(MorphologicalTreePtr tree, AttributeComputer& comp, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]> >>& dependencySources={});
//	static  std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>> computeSingleAttribute(MorphologicalTreePtr tree, Attribute attr, const DependencyMap& available);

	//static std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>> computeSingleAttribute(MorphologicalTreePtr tree, Attribute attribute, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources = {});
	
	//[atributoS do nó 0] [atributoS do nó 1] [atributoS do nó 2] ...
	//static std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>> computeAttributes(MorphologicalTreePtr tree, const std::vector<Attribute>& computers, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources = {});


//static std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>> computeAttributes(MorphologicalTreePtr tree, const std::vector<Attribute>& attributes, const DependencyMap& dependencies);

	// Define isso no começo do .hpp

	// computeSingleAttribute
	static std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>> computeSingleAttribute(MorphologicalTreePtr tree, Attribute attr, const DependencyMap& available = {});

	// computeAttributes
	static std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>> computeAttributes(MorphologicalTreePtr tree, const std::vector<Attribute>& attributes, const DependencyMap& dependencies = {});
};



class AreaComputer : public AttributeComputer {
	public:
		using enum GeometricAttribute;
		std::vector<Attribute> attributes() const override {
			return {AREA};
		}
		void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources= {}) const override {
			std::cout << "\n==== AttributeComputer: Computing AREA" << std::endl;
			auto indexOf = [&](int idx) {
				return attrNames->linearIndex(idx, AREA);
			};
	
			AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
				[&](NodeMTPtr node) {
					buffer[indexOf(node->getIndex())] = node->getCNPs().size();
				},
				[&](NodeMTPtr parent, NodeMTPtr child) {
					buffer[indexOf(parent->getIndex())] += buffer[indexOf(child->getIndex())];
				},
				[](NodeMTPtr node) {}
			);
		}
};


class VolumeComputer : public AttributeComputer {
	public:
		using enum GeometricAttribute;
		std::vector<Attribute> attributes() const override {
			return {VOLUME};
		}
		void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources= {}) const override {
			std::cout << "\n==== AttributeComputer: Computing VOLUME" << std::endl;
			auto indexOf = [&](int idx) {
				return attrNames->linearIndex(idx, VOLUME);
			};
	
			AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
				[&](NodeMTPtr node) {
					buffer[indexOf(node->getIndex())] = node->getCNPs().size() * node->getLevel();
				},
				[&](NodeMTPtr parent, NodeMTPtr child) {
					buffer[indexOf(parent->getIndex())] += buffer[indexOf(child->getIndex())];
				},
				[](NodeMTPtr node) {}
			);
		}
};

class LevelComputer : public AttributeComputer {
	public:
		using enum GeometricAttribute;
		std::vector<Attribute> attributes() const override {
			return {LEVEL};
		}
		void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources= {}) const override {
			std::cout << "\n==== AttributeComputer: Computing LEVEL" << std::endl;
			auto indexOf = [&](int idx) {
				return attrNames->linearIndex(idx, LEVEL);
			};
	
			AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
				[&](NodeMTPtr node) {
					buffer[ indexOf(node->getIndex())] = node->getLevel();
				},
				[&](NodeMTPtr parent, NodeMTPtr child) { 
					
				},
				[](NodeMTPtr node) {

				}
			);
		}
};

class DynamicsComputer : public AttributeComputer {
	public:
		using enum GeometricAttribute;
		std::vector<Attribute> attributes() const override {
			return {DYNAMICS};
		}
		void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources= {}) const override {
			std::cout << "\n==== AttributeComputer: Computing DYNAMICS" << std::endl;
			auto indexOf = [&](int idx) {
				return attrNames->linearIndex(idx, DYNAMICS);
			};
	
			AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
				[&](NodeMTPtr node) {
					if (node->isMaxtreeNode()) {
						buffer[ indexOf(node->getIndex())] = std::numeric_limits<int>::lowest(); // Procuramos máximo para max-tree
					} else {
						buffer[ indexOf(node->getIndex())] = std::numeric_limits<int>::max(); // Procuramos mínimo para min-tree
					}
				},
				[&](NodeMTPtr parent, NodeMTPtr child) {
					if (parent->isMaxtreeNode()) {
						buffer[ indexOf(parent->getIndex()) ] = std::max(buffer[ indexOf(parent->getIndex())], buffer[ indexOf(child->getIndex())] );
					} else {
						buffer[ indexOf(parent->getIndex()) ] = std::min(buffer[ indexOf(parent->getIndex())], buffer[ indexOf(child->getIndex())] );
					}
				},
				[&](NodeMTPtr node) {
					if (node->getChildren().empty()) {
						buffer[ indexOf(node->getIndex()) ] = 0; // Folhas têm dinâmica 0
					} else {
						buffer[ indexOf(node->getIndex()) ] = std::abs( node->getLevel() - buffer[ indexOf(node->getIndex())] ) + 1;
					}
				}
			);
		}
};


class MeanLevelComputer : public AttributeComputer {
	public:
		using enum GeometricAttribute;
	
		std::vector<Attribute> attributes() const override {
			return { MEAN_LEVEL };
		}
	
		std::vector<Attribute> requiredAttributes() const override {
			return { AREA, VOLUME };
		}
	
		void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources = {}) const override {
			
			std::cout << "\n==== AttributeComputer: Computing MEAN_LEVEL" << std::endl;
	
			auto indexOf = [&](int idx) {
				return attrNames->linearIndex(idx, MEAN_LEVEL);
			};
			
			auto [dependencyAttrNamesArea, bufferArea] = dependencySources[0]; //area
			auto indexOfArea = [&](int idx) {
				return dependencyAttrNamesArea->linearIndex(idx, AREA);
			};

			auto [dependencyAttrNamesVol, bufferVol] = dependencySources[1]; //volume
			auto indexOfVol = [&](int idx) {
				return dependencyAttrNamesVol->linearIndex(idx, VOLUME);
			};
	
			AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
				[&](NodeMTPtr node) {
					
				},
				[&](NodeMTPtr parent, NodeMTPtr child) { 
					
				},
				[&](NodeMTPtr node) {
					buffer[ indexOf(node->getIndex()) ] = bufferVol[indexOfVol(node->getIndex())] / bufferArea[indexOfArea(node->getIndex())];
				}
			);
		}
};


class VarianceLevelComputer : public AttributeComputer {
	public:
		using enum GeometricAttribute;
	
		std::vector<Attribute> attributes() const override {
			return { VARIANCE_LEVEL };
		}
	
		std::vector<Attribute> requiredAttributes() const override {
			return { AREA, MEAN_LEVEL };
		}
	
		void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources = {}) const override {
			
			std::cout << "\n==== AttributeComputer: Computing VARIANCE_LEVEL" << std::endl;
	
			auto indexOf = [&](int idx) {
				return attrNames->linearIndex(idx, VARIANCE_LEVEL);
			};
	
			auto [namesArea, bufArea]     = dependencySources[0];
			auto indexOfArea = [&](int idx) {
				return namesArea->linearIndex(idx, AREA);
			};
	
			auto [namesMean, bufMean]     = dependencySources[1];
			auto indexOfMean = [&](int idx) {
				return namesMean->linearIndex(idx, MEAN_LEVEL);
			};
	
			const int numNodes = tree->getNumNodes();
			std::unique_ptr<long[]> sumGrayLevelSquare(new long[numNodes]);
	
			AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
				[&](NodeMTPtr node) {
					sumGrayLevelSquare[node->getIndex()] = node->getCNPs().size() * std::pow(node->getLevel(), 2);
				},
				[&](NodeMTPtr parent, NodeMTPtr child) {
					sumGrayLevelSquare[parent->getIndex()] += sumGrayLevelSquare[child->getIndex()];
				},
				[&](NodeMTPtr node) {
					int idx = node->getIndex();
					float area = bufArea[indexOfArea(idx)];
					
					float meanGrayLevel = bufMean[indexOfMean(idx)]; //mean graylevel - // E(f)
					double meanGrayLevelSquare = sumGrayLevelSquare[idx] / area; // E(f^2)
					float var = meanGrayLevelSquare - (meanGrayLevel * meanGrayLevel); //variance: E(f^2) - E(f)^2
					buffer[indexOf(idx)] = var > 0.0f ? var : 0.0f; //variance
					
					
				}
			);
		}
};

class BoundingBoxComputer : public AttributeComputer {
	public:
		using enum GeometricAttribute;
	
		std::vector<Attribute> attributes() const override {
			return { BOX_WIDTH, BOX_HEIGHT };
		}
	
		void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources = {}) const override {
			
			std::cout << "\n==== AttributeComputer: Computing BOX_WIDTH and BOX_HEIGHT" << std::endl;
	
			auto indexOfWidth  = [&](int idx) { return attrNames->linearIndex(idx, BOX_WIDTH); };
			auto indexOfHeight = [&](int idx) { return attrNames->linearIndex(idx, BOX_HEIGHT); };
	
			int n = tree->getNumNodes();
			int numCols = tree->getNumColsOfImage();
			int numRows = tree->getNumRowsOfImage();
	
			std::unique_ptr<int[]> xmin(new int[n]);
			std::unique_ptr<int[]> xmax(new int[n]);
			std::unique_ptr<int[]> ymin(new int[n]);
			std::unique_ptr<int[]> ymax(new int[n]);
	
			AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
				[&](NodeMTPtr node) {
					int idx = node->getIndex();
					xmin[idx] = numCols;
					xmax[idx] = 0;
					ymin[idx] = numRows;
					ymax[idx] = 0;
	
					for (int p : node->getCNPs()) {
						auto [y, x] = ImageUtils::to2D(p, numCols);
						xmin[idx] = std::min(xmin[idx], x);
						xmax[idx] = std::max(xmax[idx], x);
						ymin[idx] = std::min(ymin[idx], y);
						ymax[idx] = std::max(ymax[idx], y);
					}
				},
				[&](NodeMTPtr parent, NodeMTPtr child) {
					int pid = parent->getIndex();
					int cid = child->getIndex();
					xmin[pid] = std::min(xmin[pid], xmin[cid]);
					xmax[pid] = std::max(xmax[pid], xmax[cid]);
					ymin[pid] = std::min(ymin[pid], ymin[cid]);
					ymax[pid] = std::max(ymax[pid], ymax[cid]);
				},
				[&](NodeMTPtr node) {
					int idx = node->getIndex();
					buffer[indexOfWidth(idx)]  = xmax[idx] - xmin[idx] + 1;
					buffer[indexOfHeight(idx)] = ymax[idx] - ymin[idx] + 1;
				}
			);
		}
};

class RectangularityComputer : public AttributeComputer {
	public:
		using enum GeometricAttribute;
	
		std::vector<Attribute> attributes() const override {
			return { RECTANGULARITY };
		}
	
		std::vector<Attribute> requiredAttributes() const override {
			return { AREA, BOX_WIDTH, BOX_HEIGHT };
		}
	
		void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources = {}) const override {
	
			std::cout << "\n==== AttributeComputer: Computing RECTANGULARITY" << std::endl;
	
			auto indexOf = [&](int idx) {
				return attrNames->linearIndex(idx, RECTANGULARITY);
			};
	
			auto [namesArea, bufArea] = dependencySources[0];
			auto indexArea      = [&](int idx) { return namesArea->linearIndex(idx, AREA); };
			
			auto [namesBox, bufBox] = dependencySources[1];
			auto indexBoxWidth  = [&](int idx) { return namesBox->linearIndex(idx, BOX_WIDTH); };
			auto indexBoxHeight = [&](int idx) { return namesBox->linearIndex(idx, BOX_HEIGHT); };
	
			AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
				[&](NodeMTPtr) { },
				[&](NodeMTPtr, NodeMTPtr) { },
				[&](NodeMTPtr node) {
					int idx = node->getIndex();
					float area  = bufArea[indexArea(idx)];
					float width = bufBox[indexBoxWidth(idx)];
					float height = bufBox[indexBoxHeight(idx)];
					float denom = width * height;
					buffer[indexOf(idx)] = (denom > 0.0f) ? (area / denom) : 0.0f;
				}
			);
		}
};

class RatioWHComputer : public AttributeComputer {
	public:
		using enum GeometricAttribute;
	
		std::vector<Attribute> attributes() const override {
			return { RATIO_WH };
		}
	
		std::vector<Attribute> requiredAttributes() const override {
			return { BOX_WIDTH, BOX_HEIGHT };
		}
	
		void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources = {}) const override {
	
			std::cout << "\n==== AttributeComputer: Computing RATIO_WH" << std::endl;
	
			auto indexOf = [&](int idx) {
				return attrNames->linearIndex(idx, RATIO_WH);
			};
	
			auto [namesWH, bufWH] = dependencySources[0];
			auto indexBoxWidth  = [&](int idx) { return namesWH->linearIndex(idx, BOX_WIDTH); };
			auto indexBoxHeight = [&](int idx) { return namesWH->linearIndex(idx, BOX_HEIGHT); };
	
			AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
				[&](NodeMTPtr) { },
				[&](NodeMTPtr, NodeMTPtr) { },
				[&](NodeMTPtr node) {
					int idx = node->getIndex();
					float width  = bufWH[indexBoxWidth(idx)];
					float height = bufWH[indexBoxHeight(idx)];
					if (width > 0 && height > 0) {
						buffer[indexOf(idx)] = std::max(width, height) / std::min(width, height);
					} else {
						buffer[indexOf(idx)] = 0.0f;
					}
				}
			);
		}
};

class CentralMomentsComputer : public AttributeComputer {
	public:
		using enum GeometricAttribute;

		std::vector<Attribute> attributes() const override {
			return {CENTRAL_MOMENT_00,
					CENTRAL_MOMENT_20,
					CENTRAL_MOMENT_02,
					CENTRAL_MOMENT_11,
					CENTRAL_MOMENT_30,
					CENTRAL_MOMENT_03,
					CENTRAL_MOMENT_21,
					CENTRAL_MOMENT_12};
		}
		

		void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources= {}) const override {
			std::cout << "\n==== AttributeComputer: Computing CENTRAL_MOMENT" << std::endl;
			//momentos geometricos para calcular o centroide
			int numCols = tree->getNumColsOfImage();
			int n = tree->getNumColsOfImage() * tree->getNumRowsOfImage();
			std::unique_ptr<long int[]> sumX(new long int[n]);//sum x
			std::unique_ptr<long int[]> sumY(new long int[n]);//sum y
			
			auto indexOf = [&](int idx, GeometricAttribute attr) {
				return attrNames->linearIndex(idx, attr);
			};
			
			//computa sumX e sumY para calcular os centroides
			AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
				[&](NodeMTPtr node) -> void {
					buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_00)] = node->getCNPs().size(); //area
					sumX[node->getIndex()] = 0;
					sumY[node->getIndex()] = 0;
					for(int p: node->getCNPs()) {
						auto [py, px] = ImageUtils::to2D(p, numCols); 
						sumX[node->getIndex()] += px;
						sumY[node->getIndex()] += py;
					}
				},
				[&](NodeMTPtr parent, NodeMTPtr child) -> void {
					buffer[indexOf(parent->getIndex(), CENTRAL_MOMENT_00)] += buffer[indexOf(child->getIndex(), CENTRAL_MOMENT_00)];
					sumX[parent->getIndex()] += sumX[child->getIndex()];
					sumY[parent->getIndex()] += sumY[child->getIndex()];				
				},
				[](NodeMTPtr node) -> void {
					// Não é necessário fazer nada aqui
				}
			);
	
			
			//Computação dos momentos centrais
			AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
				[&](NodeMTPtr node) -> void {				
					// Inicialização dos momentos centrais
					buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_20)] = 0.0f;
					buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_02)] = 0.0f;
					buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_11)] = 0.0f;
					buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_30)] = 0.0f;
					buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_03)] = 0.0f;
					buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_21)] = 0.0f;
					buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_12)] = 0.0f;
	
					// Cálculo do centroide
					float xCentroid = sumX[node->getIndex()] / buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_00)];
					float yCentroid = sumY[node->getIndex()] / buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_00)];
	
					for (int p : node->getCNPs()) {
						auto [py, px] = ImageUtils::to2D(p, numCols); 
						float dx = px - xCentroid;
						float dy = py - yCentroid;
	
						// Momentos centrais de segunda ordem
						buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_20)] += std::pow(dx, 2);
						buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_02)] += std::pow(dy, 2);
						buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_11)] += dx * dy;
	
						// Momentos centrais de terceira ordem
						buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_30)] += std::pow(dx, 3);
						buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_03)] += std::pow(dy, 3);
						buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_21)] += std::pow(dx, 2) * dy;
						buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_12)] += dx * std::pow(dy, 2);
					}
	
				},
				[&](NodeMTPtr parent, NodeMTPtr child) -> void {
					buffer[indexOf(parent->getIndex(), CENTRAL_MOMENT_20)] += buffer[indexOf(child->getIndex(), CENTRAL_MOMENT_20)];
					buffer[indexOf(parent->getIndex(), CENTRAL_MOMENT_02)] += buffer[indexOf(child->getIndex(), CENTRAL_MOMENT_02)];
					buffer[indexOf(parent->getIndex(), CENTRAL_MOMENT_11)] += buffer[indexOf(child->getIndex(), CENTRAL_MOMENT_11)];
					buffer[indexOf(parent->getIndex(), CENTRAL_MOMENT_30)] += buffer[indexOf(child->getIndex(), CENTRAL_MOMENT_30)];
					buffer[indexOf(parent->getIndex(), CENTRAL_MOMENT_03)] += buffer[indexOf(child->getIndex(), CENTRAL_MOMENT_03)];
					buffer[indexOf(parent->getIndex(), CENTRAL_MOMENT_21)] += buffer[indexOf(child->getIndex(), CENTRAL_MOMENT_21)];
					buffer[indexOf(parent->getIndex(), CENTRAL_MOMENT_12)] += buffer[indexOf(child->getIndex(), CENTRAL_MOMENT_12)];			
				},
				[](NodeMTPtr node) -> void {
					// Não é necessário fazer nada aqui
				}
		);
	}
};

class LengthMajorAxisComputer : public AttributeComputer {
	public:
		using enum GeometricAttribute;
	
		std::vector<Attribute> attributes() const override {
			return { LENGTH_MAJOR_AXIS };
		}
	
		std::vector<Attribute> requiredAttributes() const override {
			return { CENTRAL_MOMENT_00, CENTRAL_MOMENT_20, CENTRAL_MOMENT_02, CENTRAL_MOMENT_11 };
		}
	
		void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources = {}) const override {
			
			std::cout << "\n==== AttributeComputer: Computing LENGTH_MAJOR_AXIS" << std::endl;
	
			auto indexOf = [&](int idx) {
				return attrNames->linearIndex(idx, LENGTH_MAJOR_AXIS);
			};
	
			auto [names, buf] = dependencySources[0];
			auto indexArea = [&](int idx) { return names->linearIndex(idx, CENTRAL_MOMENT_00); };
			auto indexMu20 = [&](int idx) { return names->linearIndex(idx, CENTRAL_MOMENT_20); };
			auto indexMu02 = [&](int idx) { return names->linearIndex(idx, CENTRAL_MOMENT_02); };
			auto indexMu11 = [&](int idx) { return names->linearIndex(idx, CENTRAL_MOMENT_11); };
	
			AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
				[&](NodeMTPtr) {},
				[&](NodeMTPtr, NodeMTPtr) {},
				[&](NodeMTPtr node) {
					int idx = node->getIndex();
					float mu20 = buf[indexMu20(idx)];
					float mu02 = buf[indexMu02(idx)];
					float mu11 = buf[indexMu11(idx)];
					float area = buf[indexArea(idx)];
	
					float discriminant = std::pow(mu20 - mu02, 2.0f) + 4.0f * std::pow(mu11, 2.0f);
					discriminant = std::max(discriminant, 0.0f);
	
					float lambda1 = mu20 + mu02 + std::sqrt(discriminant);
	
					if (area > 0.0f && lambda1 > 0.0f) {
						buffer[indexOf(idx)] = std::sqrt((2.0f * lambda1) / area);
					} else {
						buffer[indexOf(idx)] = 0.0f;
					}
				}
			);
		}
};

class LengthMinorAxisComputer : public AttributeComputer {
	public:
		using enum GeometricAttribute;
	
		std::vector<Attribute> attributes() const override {
			return { LENGTH_MINOR_AXIS };
		}
	
		std::vector<Attribute> requiredAttributes() const override {
			return { CENTRAL_MOMENT_00, CENTRAL_MOMENT_20, CENTRAL_MOMENT_02, CENTRAL_MOMENT_11 };
		}
	
		void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources = {}) const override {
			
			std::cout << "\n==== AttributeComputer: Computing LENGTH_MINOR_AXIS" << std::endl;
	
			auto indexOf = [&](int idx) {
				return attrNames->linearIndex(idx, LENGTH_MINOR_AXIS);
			};
	
			auto [names, buf] = dependencySources[0];
			auto indexArea = [&](int idx) { return names->linearIndex(idx, CENTRAL_MOMENT_00); };
			auto indexMu20 = [&](int idx) { return names->linearIndex(idx, CENTRAL_MOMENT_20); };
			auto indexMu02 = [&](int idx) { return names->linearIndex(idx, CENTRAL_MOMENT_02); };
			auto indexMu11 = [&](int idx) { return names->linearIndex(idx, CENTRAL_MOMENT_11); };
	
			AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
				[&](NodeMTPtr) {},
				[&](NodeMTPtr, NodeMTPtr) {},
				[&](NodeMTPtr node) {
					int idx = node->getIndex();
					float mu20 = buf[indexMu20(idx)];
					float mu02 = buf[indexMu02(idx)];
					float mu11 = buf[indexMu11(idx)];
					float area = buf[indexArea(idx)];
	
					float discriminant = std::pow(mu20 - mu02, 2.0f) + 4.0f * std::pow(mu11, 2.0f);
					discriminant = std::max(discriminant, 0.0f);
	
					float lambda2 = mu20 + mu02 - std::sqrt(discriminant);
	
					if (area > 0.0f && lambda2 > 0.0f) {
						buffer[indexOf(idx)] = std::sqrt((2.0f * lambda2) / area);
					} else {
						buffer[indexOf(idx)] = 0.0f;
					}
				}
			);
		}
};

class EccentricityComputer : public AttributeComputer {
	public:
		using enum GeometricAttribute;
	
		std::vector<Attribute> attributes() const override {
			return { ECCENTRICITY };
		}
	
		std::vector<Attribute> requiredAttributes() const override {
			return { CENTRAL_MOMENT_20, CENTRAL_MOMENT_02, CENTRAL_MOMENT_11 };
		}
	
		void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources = {}) const override {
	
			std::cout << "\n==== AttributeComputer: Computing ECCENTRICITY" << std::endl;
	
			auto indexOf = [&](int idx) {
				return attrNames->linearIndex(idx, ECCENTRICITY);
			};
	
			auto [names, buf] = dependencySources[0];
			auto indexMu20 = [&](int idx) { return names->linearIndex(idx, CENTRAL_MOMENT_20); };
			auto indexMu02 = [&](int idx) { return names->linearIndex(idx, CENTRAL_MOMENT_02); };
			auto indexMu11 = [&](int idx) { return names->linearIndex(idx, CENTRAL_MOMENT_11); };
	
			AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
				[&](NodeMTPtr) {},
				[&](NodeMTPtr, NodeMTPtr) {},
				[&](NodeMTPtr node) {
					int idx = node->getIndex();
					float mu20 = buf[indexMu20(idx)];
					float mu02 = buf[indexMu02(idx)];
					float mu11 = buf[indexMu11(idx)];
	
					float discriminant = std::pow(mu20 - mu02, 2.0f) + 4.0f * std::pow(mu11, 2.0f);
					discriminant = std::max(discriminant, 0.0f);
	
					float lambda1 = mu20 + mu02 + std::sqrt(discriminant);  // maior autovalor
					float lambda2 = mu20 + mu02 - std::sqrt(discriminant);  // menor autovalor
	
					if (std::abs(lambda2) > std::numeric_limits<float>::epsilon()) {
						buffer[indexOf(idx)] = lambda1 / lambda2;
					} else {
						buffer[indexOf(idx)] = lambda1 / 0.1f; // fallback para evitar divisão por zero
					}
				}
			);
		}
};

class CompactnessComputer : public AttributeComputer {
	public:
		using enum GeometricAttribute;
	
		std::vector<Attribute> attributes() const override {
			return { COMPACTNESS };
		}
	
		std::vector<Attribute> requiredAttributes() const override {
			return { CENTRAL_MOMENT_00, CENTRAL_MOMENT_20, CENTRAL_MOMENT_02 };
		}
	
		void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources = {}) const override {
	
			std::cout << "\n==== AttributeComputer: Computing COMPACTNESS" << std::endl;
	
			auto indexOf = [&](int idx) {
				return attrNames->linearIndex(idx, COMPACTNESS);
			};
	
			auto [names, buf] = dependencySources[0];
			auto indexArea     = [&](int idx) { return names->linearIndex(idx, CENTRAL_MOMENT_00); };
			auto indexMu20     = [&](int idx) { return names->linearIndex(idx, CENTRAL_MOMENT_20); };
			auto indexMu02     = [&](int idx) { return names->linearIndex(idx, CENTRAL_MOMENT_02); };
	
			AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
				[&](NodeMTPtr) {},
				[&](NodeMTPtr, NodeMTPtr) {},
				[&](NodeMTPtr node) {
					int idx = node->getIndex();
					float area = buf[indexArea(idx)];
					float mu20 = buf[indexMu20(idx)];
					float mu02 = buf[indexMu02(idx)];
					float denom = mu20 + mu02;
	
					if (denom > std::numeric_limits<float>::epsilon()) {
						buffer[indexOf(idx)] = (1.0f / (2.0f * static_cast<float>(M_PI))) * (area / denom);
					} else {
						buffer[indexOf(idx)] = 0.0f;
					}
				}
			);
		}
};

class OrientationComputer : public AttributeComputer {
	private:
		bool normalize0to180;
	public:
		using enum GeometricAttribute;

		OrientationComputer(): normalize0to180(false) {}
		OrientationComputer(bool normalize0to180): normalize0to180(normalize0to180) {}

		std::vector<Attribute> attributes() const override {
			return {AXIS_ORIENTATION};
		}
		std::vector<Attribute> requiredAttributes() const override{ 
			return {CENTRAL_MOMENT_00, CENTRAL_MOMENT_20, CENTRAL_MOMENT_02, CENTRAL_MOMENT_11, CENTRAL_MOMENT_30, CENTRAL_MOMENT_03, CENTRAL_MOMENT_21, CENTRAL_MOMENT_12}; 
		}

		void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources= {}) const override {
			std::cout << "\n==== AttributeComputer: Computing AXIS_ORIENTATION" << std::endl;

			int numCols = tree->getNumColsOfImage();
			auto indexOf = [&](int idx) {
				return attrNames->linearIndex(idx, AXIS_ORIENTATION);
			};
			
			auto [dependencyAttrNamesMu, bufferMu] = dependencySources[0]; //momentos centrais
			auto indexOfMu = [&](int idx, GeometricAttribute attr) {
				return dependencyAttrNamesMu->linearIndex(idx, attr);
			};
			
			//Computação da orientação em graus
			AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
				[](NodeMTPtr node) -> void { },
				[](NodeMTPtr parent, NodeMTPtr child) -> void { },
				[&](NodeMTPtr node) -> void {
					int idx = node->getIndex();
					//Momentos centrais
					float mu20 = bufferMu[indexOfMu(idx, CENTRAL_MOMENT_20)];
					float mu02 = bufferMu[indexOfMu(idx, CENTRAL_MOMENT_02)];
					float mu11 = bufferMu[indexOfMu(idx, CENTRAL_MOMENT_11)];
					
					// Verificar se o denominador é zero antes de calcular atan2 para evitar divisão por zero
					if (mu20 != mu02 || mu11 != 0) {
						float radians = 0.5 * std::atan2(2 * mu11, mu20 - mu02);// orientação em radianos
						float degrees = radians * (180.0 / M_PI); // Converter para graus
						if(normalize0to180)
							buffer[indexOf(idx)] = std::fmod(std::abs(degrees), 180.0f); // Armazenar a orientação do eixo maior [0, 180]
						else
							buffer[indexOf(idx)] = std::fmod(std::abs(degrees), 360.0f); ; // Armazenar a orientação no intervalo [0, 360]
					} else {
						buffer[indexOf(idx)] = 0.0; // Se não for possível calcular a orientação, definir um valor padrão
					}
				}
		);
	}
};


class HuMomentsComputer : public AttributeComputer {
	public:
		using enum GeometricAttribute;

		std::vector<Attribute> attributes() const override {
			return {HU_MOMENT_1,
					HU_MOMENT_2,
					HU_MOMENT_3,
					HU_MOMENT_4,
					HU_MOMENT_5,
					HU_MOMENT_6,
					HU_MOMENT_7};
		}
		std::vector<Attribute> requiredAttributes() const override{ 
			return {CENTRAL_MOMENT_00, CENTRAL_MOMENT_20, CENTRAL_MOMENT_02, CENTRAL_MOMENT_11, CENTRAL_MOMENT_30, CENTRAL_MOMENT_03, CENTRAL_MOMENT_21, CENTRAL_MOMENT_12}; 
		}

		void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources= {}) const override {
			std::cout << "\n==== AttributeComputer: Computing HU_MOMENT" << std::endl;
			int numCols = tree->getNumColsOfImage();
			auto indexOf = [&](int idx, GeometricAttribute attr) {
				return attrNames->linearIndex(idx, attr);
			};
			
			auto [dependencyAttrNamesMu, bufferMu] = dependencySources[0]; //momentos centrais
			auto indexOfMu = [&](int idx, GeometricAttribute attr) {
				return dependencyAttrNamesMu->linearIndex(idx, attr);
			};

			auto normMoment = [](int area, float moment, int p, int q){ 
				return moment / std::pow( area, (p + q + 2.0) / 2.0); 
			}; //função para normalizacao dos momentos		

			
			//Computação dos momentos de Hu
			AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
				[](NodeMTPtr node) -> void { },
				[](NodeMTPtr parent, NodeMTPtr child) -> void { },
				[&](NodeMTPtr node) -> void {
					int idx = node->getIndex();
					//Momentos centrais
					float mu20 = bufferMu[indexOfMu(idx, CENTRAL_MOMENT_20)];
					float mu02 = bufferMu[indexOfMu(idx, CENTRAL_MOMENT_02)];
					float mu11 = bufferMu[indexOfMu(idx, CENTRAL_MOMENT_11)];
					float mu30 = bufferMu[indexOfMu(idx, CENTRAL_MOMENT_30)];
					float mu03 = bufferMu[indexOfMu(idx, CENTRAL_MOMENT_03)];
					float mu21 = bufferMu[indexOfMu(idx, CENTRAL_MOMENT_21)];
					float mu12 = bufferMu[indexOfMu(idx, CENTRAL_MOMENT_12)];
					int area = bufferMu[indexOfMu(idx, CENTRAL_MOMENT_00)];

					
					// Calcular os momentos normalizados
					float eta20 = normMoment(area, mu20, 2, 0);
					float eta02 = normMoment(area, mu02, 0, 2);
					float eta11 = normMoment(area, mu11, 1, 1);
					float eta30 = normMoment(area, mu30, 3, 0);
					float eta03 = normMoment(area, mu03, 0, 3);
					float eta21 = normMoment(area, mu21, 2, 1);
					float eta12 = normMoment(area, mu12, 1, 2);

					// Cálculo dos momentos de Hu
					buffer[indexOf(idx, HU_MOMENT_1)] = eta20 + eta02; // primeiro momento de Hu => inertia
					buffer[indexOf(idx, HU_MOMENT_2)]  = std::pow(eta20 - eta02, 2) + 4 * std::pow(eta11, 2);
					buffer[indexOf(idx, HU_MOMENT_3)]  = std::pow(eta30 - 3 * eta12, 2) + std::pow(3 * eta21 - eta03, 2);
					buffer[indexOf(idx, HU_MOMENT_4)]  = std::pow(eta30 + eta12, 2) + std::pow(eta21 + eta03, 2);
					
					buffer[indexOf(idx, HU_MOMENT_5)] = (eta30 - 3 * eta12) * (eta30 + eta12) * (std::pow(eta30 + eta12, 2) - 3 * std::pow(eta21 + eta03, 2)) +
														(3 * eta21 - eta03) * (eta21 + eta03) * (3 * std::pow(eta30 + eta12, 2) - std::pow(eta21 + eta03, 2));
					
					buffer[indexOf(idx, HU_MOMENT_6)] = (eta20 - eta02) * (std::pow(eta30 + eta12, 2) - std::pow(eta21 + eta03, 2)) + 
														4 * eta11 * (eta30 + eta12) * (eta21 + eta03);
					
					buffer[indexOf(idx, HU_MOMENT_7)] = (3 * eta21 - eta03) * (eta30 + eta12) * (std::pow(eta30 + eta12, 2) - 3 * std::pow(eta21 + eta03, 2)) -
														(eta30 - 3 * eta12) * (eta21 + eta03) * (3 * std::pow(eta30 + eta12, 2) - std::pow(eta21 + eta03, 2));

				}
		);
	}
};

class InertiaComputer : public AttributeComputer {
public:
    using enum GeometricAttribute;

    std::vector<Attribute> attributes() const override {
        return { INERTIA };
    }

    std::vector<Attribute> requiredAttributes() const override {
        return { CENTRAL_MOMENT_00, CENTRAL_MOMENT_20, CENTRAL_MOMENT_02 };
    }

    void compute(MorphologicalTreePtr tree,
                 std::shared_ptr<float[]> buffer,
                 std::shared_ptr<AttributeNames> attrNames,
                 const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources = {}) const override {

        std::cout << "\n==== AttributeComputer: Computing INERTIA" << std::endl;

        auto indexOf = [&](int idx) {
            return attrNames->linearIndex(idx, INERTIA);
        };

        auto [names, buf] = dependencySources[0];
        auto indexMu00 = [&](int idx) { return names->linearIndex(idx, CENTRAL_MOMENT_00); };
        auto indexMu20 = [&](int idx) { return names->linearIndex(idx, CENTRAL_MOMENT_20); };
        auto indexMu02 = [&](int idx) { return names->linearIndex(idx, CENTRAL_MOMENT_02); };

        AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
            [&](NodeMTPtr) {},
            [&](NodeMTPtr, NodeMTPtr) {},
            [&](NodeMTPtr node) {
                int idx = node->getIndex();
                float mu00 = buf[indexMu00(idx)];
                float mu20 = buf[indexMu20(idx)];
                float mu02 = buf[indexMu02(idx)];
 
                float normMu20 = mu20 / std::pow(mu00, 2.0f);
                float normMu02 = mu02 / std::pow(mu00, 2.0f);

                buffer[indexOf(idx)] = normMu20 + normMu02;
            }
        );
    }
};


class AttributeFactory {
	public:
		static std::shared_ptr<AttributeComputer> create(Attribute attr) {
			using enum GeometricAttribute;
	
			switch (std::get<GeometricAttribute>(attr)) {
				case AREA: return std::make_shared<AreaComputer>();
				case VOLUME: return std::make_shared<VolumeComputer>();
				case LEVEL: return std::make_shared<LevelComputer>();
				case DYNAMICS: return std::make_shared<DynamicsComputer>();
				case MEAN_LEVEL: return std::make_shared<MeanLevelComputer>();
				case VARIANCE_LEVEL: return std::make_shared<VarianceLevelComputer>();
				case BOX_HEIGHT:
				case BOX_WIDTH:
				case BOUNDING_BOX: return std::make_shared<BoundingBoxComputer>();
				case RECTANGULARITY: return std::make_shared<RectangularityComputer>();
				case RATIO_WH: return std::make_shared<RatioWHComputer>();
				
				case AXIS_ORIENTATION: return std::make_shared<OrientationComputer>();
				case LENGTH_MAJOR_AXIS: return std::make_shared<LengthMajorAxisComputer>();
				case LENGTH_MINOR_AXIS: return std::make_shared<LengthMinorAxisComputer>();
				case ECCENTRICITY: return std::make_shared<EccentricityComputer>();
				case COMPACTNESS: return std::make_shared<CompactnessComputer>();

				case CENTRAL_MOMENT: 
				case CENTRAL_MOMENT_00: 
				case CENTRAL_MOMENT_20:
				case CENTRAL_MOMENT_02:
				case CENTRAL_MOMENT_11:
				case CENTRAL_MOMENT_30:
				case CENTRAL_MOMENT_03:
				case CENTRAL_MOMENT_21:
				case CENTRAL_MOMENT_12:
					return std::make_shared<CentralMomentsComputer>();

				case INERTIA:
					return std::make_shared<InertiaComputer>();	
				case HU_MOMENT:
				case HU_MOMENT_1: 
				case HU_MOMENT_2:
				case HU_MOMENT_3:
				case HU_MOMENT_4:
				case HU_MOMENT_5:
				case HU_MOMENT_6:
				case HU_MOMENT_7:
					return std::make_shared<HuMomentsComputer>();

				default:
					throw std::invalid_argument("AttributeFactory: atributo não suportado: " + AttributeNames::toString(attr));
			}
		}
};







#endif 






		