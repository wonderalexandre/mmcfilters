

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
    STANDARD_DEVIATION,
    BOX_WIDTH,
    BOX_HEIGHT,
    RECTANGULARITY,
    RATIO_WH,
    CENTRAL_MOMENT_20,
    CENTRAL_MOMENT_02,
    CENTRAL_MOMENT_11,
    CENTRAL_MOMENT_30,
    CENTRAL_MOMENT_03,
    CENTRAL_MOMENT_21,
    CENTRAL_MOMENT_12,
    ORIENTATION,
    LENGTH_MAJOR_AXIS,
    LENGTH_MINOR_AXIS,
    ECCENTRICITY,
    COMPACTNESS,
    HU_MOMENT_1_INERTIA,
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
	
		int linearIndex(const Attribute& attr, int nodeIndex) const {
			return getIndex(attr) + nodeIndex;
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
				case GeometricAttribute::STANDARD_DEVIATION: return "STANDARD_DEVIATION";
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
				case GeometricAttribute::ORIENTATION: return "ORIENTATION";
				case GeometricAttribute::LENGTH_MAJOR_AXIS: return "LENGTH_MAJOR_AXIS";
				case GeometricAttribute::LENGTH_MINOR_AXIS: return "LENGTH_MINOR_AXIS";
				case GeometricAttribute::ECCENTRICITY: return "ECCENTRICITY";
				case GeometricAttribute::COMPACTNESS: return "COMPACTNESS";
				case GeometricAttribute::HU_MOMENT_1_INERTIA: return "HU_MOMENT_1_INERTIA";
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
			else if (str == "STANDARD_DEVIATION") return GeometricAttribute::STANDARD_DEVIATION;
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
			else if (str == "ORIENTATION") return GeometricAttribute::ORIENTATION;
			else if (str == "LENGTH_MAJOR_AXIS") return GeometricAttribute::LENGTH_MAJOR_AXIS;
			else if (str == "LENGTH_MINOR_AXIS") return GeometricAttribute::LENGTH_MINOR_AXIS;
			else if (str == "ECCENTRICITY") return GeometricAttribute::ECCENTRICITY;
			else if (str == "COMPACTNESS") return GeometricAttribute::COMPACTNESS;
			else if (str == "HU_MOMENT_1_INERTIA") return GeometricAttribute::HU_MOMENT_1_INERTIA;
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
				int parentDepth = node->getParent() ? attrs[attrNames.linearIndex(DEPTH, node->getParent()->getIndex())] : 0;
	
				attrs[attrNames.linearIndex(HEIGHT, idx)] = 0.0f; // altura
				attrs[attrNames.linearIndex(DEPTH, idx)] = node->getParent() ? parentDepth + 1 : 0; // profundidade
				attrs[attrNames.linearIndex(IS_LEAF, idx)] = node->getChildren().empty() ? 1.0f : 0.0f; // é folha
				attrs[attrNames.linearIndex(IS_ROOT, idx)] = node->getParent() == nullptr ? 1.0f : 0.0f; // é raiz
				attrs[attrNames.linearIndex(NUM_CHILDREN, idx)] = node->getChildren().size();
				attrs[attrNames.linearIndex(NUM_SIBLINGS, idx)] = node->getParent() ? node->getParent()->getChildren().size() - 1 : 0;
				attrs[attrNames.linearIndex(NUM_DESCENDANTS, idx)] = 0.0f;
				attrs[attrNames.linearIndex(NUM_LEAF_DESCENDANTS, idx)] = attrs[attrNames.linearIndex(IS_LEAF, idx)];
				attrs[attrNames.linearIndex(LEAF_RATIO, idx)] = 0.0f;
				attrs[attrNames.linearIndex(BALANCE, idx)] = 0.0f;
				attrs[attrNames.linearIndex(AVG_CHILD_HEIGHT, idx)] = 0.0f;
			},
			[&](NodeMTPtr parent, NodeMTPtr child) {
				int pIdx = parent->getIndex();
				int cIdx = child->getIndex();

				// Acumulando descendentes
				attrs[attrNames.linearIndex(NUM_DESCENDANTS, pIdx)] += attrs[attrNames.linearIndex(NUM_DESCENDANTS, cIdx)] + 1;
				attrs[attrNames.linearIndex(NUM_LEAF_DESCENDANTS, pIdx)] += attrs[attrNames.linearIndex(NUM_LEAF_DESCENDANTS, cIdx)];

				// Altura
				float childHeight = attrs[attrNames.linearIndex(HEIGHT, cIdx)];
				float& parentHeight = attrs[attrNames.linearIndex(HEIGHT, pIdx)];
				parentHeight = std::max(parentHeight, childHeight + 1);

				// Balanceamento
				float& minH = attrs[attrNames.linearIndex(BALANCE, pIdx)]; // usado como mínimo temporário
				float& sumH = attrs[attrNames.linearIndex(AVG_CHILD_HEIGHT, pIdx)];
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
				float desc = attrs[attrNames.linearIndex(NUM_DESCENDANTS, idx)];
				float folhas = attrs[attrNames.linearIndex(NUM_LEAF_DESCENDANTS, idx)];
	
				// Razão folhas/descendentes
				attrs[attrNames.linearIndex(LEAF_RATIO, idx)] = desc > 0.0f ? folhas / (desc + 1.0f) : 1.0f;
	
				// Balanceamento e média
				if (!node->getChildren().empty()) {
					float alturaMax = attrs[attrNames.linearIndex(HEIGHT, idx)];
					float alturaMin = attrs[attrNames.linearIndex(BALANCE, idx)];
					attrs[attrNames.linearIndex(BALANCE, idx)] = alturaMax - alturaMin;
	
					attrs[attrNames.linearIndex(AVG_CHILD_HEIGHT, idx)] /= node->getChildren().size();
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
			attr[idx] = attrs[ attrNames.linearIndex(attribute, idx) ];
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
			return attrNames.linearIndex(attribute, idxNode);
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
							
							
							if (attrs[indexOf(node->getIndex(), VARIANCE_LEVEL)] >= 0.0f) {
								attrs[indexOf(node->getIndex(), STANDARD_DEVIATION)] = std::sqrt(attrs[indexOf(node->getIndex(), VARIANCE_LEVEL)]);
							} else {
								attrs[indexOf(node->getIndex(), STANDARD_DEVIATION)] = 0.0f; // Se a variância for negativa, definir desvio padrão como 0
							}
							
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
					attrs[indexOf(idx, ORIENTATION)] = degrees; // Armazenar a orientação em graus no intervalo [0, 360]
				} else {
					attrs[indexOf(idx, ORIENTATION)] = 0.0; // Se não for possível calcular a orientação, definir um valor padrão
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
				attrs[indexOf(idx, HU_MOMENT_1_INERTIA)] = eta20 + eta02; // primeiro momento de Hu => inertia
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
			float extinction = std::numeric_limits<float>::max();
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

};

#endif 