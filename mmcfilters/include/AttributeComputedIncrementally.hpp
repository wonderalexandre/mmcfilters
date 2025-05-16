

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
enum class Attribute {
    // Geométricos básicos
    AREA,

    // Textura / intensidade agregada
    VOLUME,
    LEVEL,
    DYNAMICS,
    MEAN_LEVEL,
    VARIANCE_LEVEL,

    // Bounding box
    BOX_WIDTH,
    BOX_HEIGHT,
	DIAGONAL_LENGTH,
    RECTANGULARITY,
    RATIO_WH,
	BOX_COL_MIN,
	BOX_COL_MAX,
	BOX_ROW_MIN,
	BOX_ROW_MAX,

    // Momentos Centrais
    CENTRAL_MOMENT_20,
    CENTRAL_MOMENT_02,
    CENTRAL_MOMENT_11,
    CENTRAL_MOMENT_30,
    CENTRAL_MOMENT_03,
    CENTRAL_MOMENT_21,
    CENTRAL_MOMENT_12,

    // Momentos de Hu
    HU_MOMENT_1,
    HU_MOMENT_2,
    HU_MOMENT_3,
    HU_MOMENT_4,
    HU_MOMENT_5,
    HU_MOMENT_6,
    HU_MOMENT_7,

    // Atributos derivados de momentos
    INERTIA,
    COMPACTNESS,
    ECCENTRICITY,
    LENGTH_MAJOR_AXIS,
    LENGTH_MINOR_AXIS,
    AXIS_ORIENTATION,

    // Estruturais (topologia da árvore)
    HEIGHT_NODE,
    DEPTH_NODE,
    IS_LEAF_NODE,
    IS_ROOT_NODE,
    NUM_CHILDREN_NODE,
    NUM_SIBLINGS_NODE,
    NUM_DESCENDANTS_NODE,
    NUM_LEAF_DESCENDANTS_NODE,
    LEAF_RATIO_NODE,
    BALANCE_NODE,
    AVG_CHILD_HEIGHT_NODE
};



enum class AttributeGroup {
    ALL,               // Todos os atributos
    GEOMETRIC,         // Forma e proporção
	MOMENT_BASED,
    BOUNDING_BOX,      // Box width/height
    CENTRAL_MOMENTS,    // Momentos centrais
    HU_MOMENTS,         // Momentos de Hu
    TEXTURE,           // Atributos baseados em níveis de cinza
    TREE_TOPOLOGY         // Topologia da árvore
};

using AttributeOrGroup = std::variant<Attribute, AttributeGroup>;
using enum Attribute;


namespace std {
    template<>
    struct hash<AttributeGroup> {
        std::size_t operator()(const AttributeGroup& group) const noexcept {
            return static_cast<std::size_t>(group);
        }
    };

    template<>
    struct hash<Attribute> {
        std::size_t operator()(const Attribute& attr) const noexcept {
            return static_cast<std::size_t>(attr);
        }
    };
    
    template <>
    struct hash<AttributeOrGroup> {
        size_t operator()(const AttributeOrGroup& attr) const {
            return std::visit([](auto&& a) -> size_t {
                return std::hash<std::decay_t<decltype(a)>>{}(a);
            }, attr);
        }
    };

}


static const std::unordered_map<AttributeGroup, std::vector<Attribute>> ATTRIBUTE_GROUPS = {
    {AttributeGroup::GEOMETRIC, {
        AREA,
		VOLUME,
        RECTANGULARITY,
		RATIO_WH,
        COMPACTNESS,
        ECCENTRICITY,
        LENGTH_MAJOR_AXIS,
        LENGTH_MINOR_AXIS,
		INERTIA,
		DIAGONAL_LENGTH,
		BOX_WIDTH,
		BOX_HEIGHT,
        AXIS_ORIENTATION
    }},
    {AttributeGroup::MOMENT_BASED, {
        COMPACTNESS,
        ECCENTRICITY,
        LENGTH_MAJOR_AXIS,
        LENGTH_MINOR_AXIS,
        AXIS_ORIENTATION,
		INERTIA
    }},
    {AttributeGroup::BOUNDING_BOX, {
        BOX_WIDTH,
        BOX_HEIGHT,
		RECTANGULARITY,
		RATIO_WH,
		BOX_COL_MIN,
		BOX_COL_MAX,
		BOX_ROW_MIN,
		BOX_ROW_MAX,
		DIAGONAL_LENGTH
    }},
    {AttributeGroup::CENTRAL_MOMENTS, {
        CENTRAL_MOMENT_20,
        CENTRAL_MOMENT_02,
        CENTRAL_MOMENT_11,
        CENTRAL_MOMENT_30,
        CENTRAL_MOMENT_03,
        CENTRAL_MOMENT_21,
        CENTRAL_MOMENT_12
    }},
    {AttributeGroup::HU_MOMENTS, {
        HU_MOMENT_1,
        HU_MOMENT_2,
        HU_MOMENT_3,
        HU_MOMENT_4,
        HU_MOMENT_5,
        HU_MOMENT_6,
        HU_MOMENT_7
    }},
    {AttributeGroup::TEXTURE, {
        VOLUME,
        LEVEL,
        DYNAMICS,
        MEAN_LEVEL,
        VARIANCE_LEVEL
    }},
    {AttributeGroup::TREE_TOPOLOGY, {
        HEIGHT_NODE,
        DEPTH_NODE,
        IS_LEAF_NODE,
        IS_ROOT_NODE,
        NUM_CHILDREN_NODE,
        NUM_SIBLINGS_NODE,
        NUM_DESCENDANTS_NODE,
        NUM_LEAF_DESCENDANTS_NODE,
        LEAF_RATIO_NODE,
        BALANCE_NODE,
        AVG_CHILD_HEIGHT_NODE
    }},
    {AttributeGroup::ALL, [] {
        std::vector<Attribute> all;
        for (int i = 0; i <= static_cast<int>(AVG_CHILD_HEIGHT_NODE); ++i)
            all.push_back(static_cast<Attribute>(i));
        return all;
    }()}
};

class AttributeNames {
public:
    std::unordered_map<Attribute, int> indexMap;
    const int NUM_ATTRIBUTES;

    AttributeNames(std::unordered_map<Attribute, int>&& map)
        : indexMap(std::move(map)), NUM_ATTRIBUTES(static_cast<int>(indexMap.size())) {}

    static AttributeNames fromList(int n, const std::vector<Attribute>& attributes) {
        std::unordered_map<Attribute, int> map;
        int i = 0;
        for (auto attr : attributes) {
            map[attr] = i++ * n;
        }
        return AttributeNames(std::move(map));
    }

    static AttributeNames fromGroup(AttributeGroup group, int n) {
        auto it = ATTRIBUTE_GROUPS.find(group);
        return fromList(n, it->second);
    }

    int getIndex(Attribute attr) const {
        return indexMap.at(attr);
    }

    int linearIndex(int nodeIndex, Attribute attr) const {
        return nodeIndex * NUM_ATTRIBUTES + getIndex(attr);
    }

    static std::string toString(Attribute attr) {
        switch (attr) {
            case AREA: return "AREA";
            case VOLUME: return "VOLUME";
            case LEVEL: return "LEVEL";
            case DYNAMICS: return "DYNAMICS";
            case MEAN_LEVEL: return "MEAN_LEVEL";
            case VARIANCE_LEVEL: return "VARIANCE_LEVEL";
            case BOX_WIDTH: return "BOX_WIDTH";
            case BOX_HEIGHT: return "BOX_HEIGHT";
            case RECTANGULARITY: return "RECTANGULARITY";
            case RATIO_WH: return "RATIO_WH";
			case DIAGONAL_LENGTH: return "DIAGONAL_LENGTH";	
			case BOX_COL_MIN: return "BOX_COL_MIN";
			case BOX_COL_MAX: return "BOX_COL_MAX";
			case BOX_ROW_MIN: return "BOX_ROW_MIN";
			case BOX_ROW_MAX: return "BOX_ROW_MAX";
            case CENTRAL_MOMENT_20: return "CENTRAL_MOMENT_20";
            case CENTRAL_MOMENT_02: return "CENTRAL_MOMENT_02";
            case CENTRAL_MOMENT_11: return "CENTRAL_MOMENT_11";
            case CENTRAL_MOMENT_30: return "CENTRAL_MOMENT_30";
            case CENTRAL_MOMENT_03: return "CENTRAL_MOMENT_03";
            case CENTRAL_MOMENT_21: return "CENTRAL_MOMENT_21";
            case CENTRAL_MOMENT_12: return "CENTRAL_MOMENT_12";
            case HU_MOMENT_1: return "HU_MOMENT_1";
            case HU_MOMENT_2: return "HU_MOMENT_2";
            case HU_MOMENT_3: return "HU_MOMENT_3";
            case HU_MOMENT_4: return "HU_MOMENT_4";
            case HU_MOMENT_5: return "HU_MOMENT_5";
            case HU_MOMENT_6: return "HU_MOMENT_6";
            case HU_MOMENT_7: return "HU_MOMENT_7";
            case INERTIA: return "INERTIA";
            case COMPACTNESS: return "COMPACTNESS";
            case ECCENTRICITY: return "ECCENTRICITY";
            case LENGTH_MAJOR_AXIS: return "LENGTH_MAJOR_AXIS";
            case LENGTH_MINOR_AXIS: return "LENGTH_MINOR_AXIS";
            case AXIS_ORIENTATION: return "AXIS_ORIENTATION";
            case HEIGHT_NODE: return "HEIGHT_NODE";
            case DEPTH_NODE: return "DEPTH_NODE";
            case IS_LEAF_NODE: return "IS_LEAF_NODE";
            case IS_ROOT_NODE: return "IS_ROOT_NODE";
            case NUM_CHILDREN_NODE: return "NUM_CHILDREN_NODE";
            case NUM_SIBLINGS_NODE: return "NUM_SIBLINGS_NODE";
            case NUM_DESCENDANTS_NODE: return "NUM_DESCENDANTS_NODE";
            case NUM_LEAF_DESCENDANTS_NODE: return "NUM_LEAF_DESCENDANTS_NODE";
            case LEAF_RATIO_NODE: return "LEAF_RATIO_NODE";
            case BALANCE_NODE: return "BALANCE_NODE";
            case AVG_CHILD_HEIGHT_NODE: return "AVG_CHILD_HEIGHT_NODE";
            default: return "UNKNOWN";
        }
    }

 	static std::string describe(Attribute attr) {
        switch (attr) {
            // Basic geometric attributes
            case Attribute::AREA: return "Area: Number of pixels in the connected component.";
            case Attribute::VOLUME: return "Volume: Sum of gray-level values of all pixels in the connected component.";
            case Attribute::LEVEL: return "Level: Threshold used to obtaned the smallest level-set containing the connected component.";
            case Attribute::DYNAMICS: return "Dynamics: Difference between the node level and the level of its darkest ancestor.";
            case Attribute::MEAN_LEVEL: return "Mean level: Average gray-level intensity of the pixels in the connected component.";
            case Attribute::VARIANCE_LEVEL: return "Variance of level: Variance of the gray-level intensities in the connected component.";

            // Bounding box attributes
            case Attribute::BOX_WIDTH: return "Bounding box width: Width of the minimum rectangle enclosing the connected component.";
            case Attribute::BOX_HEIGHT: return "Bounding box height: Height of the minimum rectangle enclosing the connected component.";
            case Attribute::RECTANGULARITY: return "Rectangularity: Ratio between the connected component area and the area of its bounding box.";
            case Attribute::RATIO_WH: return "Aspect ratio: Ratio between width and height of the bounding box.";
            case Attribute::BOX_COL_MIN: return "Bounding box column min: Minimum column index covered by the connected component.";
            case Attribute::BOX_COL_MAX: return "Bounding box column max: Maximum column index covered by the connected component.";
            case Attribute::BOX_ROW_MIN: return "Bounding box row min: Minimum row index covered by the connected component.";
            case Attribute::BOX_ROW_MAX: return "Bounding box row max: Maximum row index covered by the connected component.";
			case Attribute::DIAGONAL_LENGTH: return "Diagonal length: Length of the diagonal of the bounding box.";

            // Central moments
            case Attribute::CENTRAL_MOMENT_20: return "Central moment (2,0): Spread of pixels along the x-axis.";
            case Attribute::CENTRAL_MOMENT_02: return "Central moment (0,2): Spread of pixels along the y-axis.";
            case Attribute::CENTRAL_MOMENT_11: return "Central moment (1,1): Correlation between x and y coordinates of pixels.";
            case Attribute::CENTRAL_MOMENT_30: return "Central moment (3,0): Skewness of pixel distribution along x-axis.";
            case Attribute::CENTRAL_MOMENT_03: return "Central moment (0,3): Skewness of pixel distribution along y-axis.";
            case Attribute::CENTRAL_MOMENT_21: return "Central moment (2,1): Mixed third-order moment.";
            case Attribute::CENTRAL_MOMENT_12: return "Central moment (1,2): Mixed third-order moment.";

            // Hu moments (invariant shape descriptors)
            case Attribute::HU_MOMENT_1: return "Hu moment 1: First invariant moment, related to object size.";
            case Attribute::HU_MOMENT_2: return "Hu moment 2: Second invariant moment, related to variance.";
            case Attribute::HU_MOMENT_3: return "Hu moment 3: Third invariant moment, capturing skewness.";
            case Attribute::HU_MOMENT_4: return "Hu moment 4: Invariant moment capturing symmetry.";
            case Attribute::HU_MOMENT_5: return "Hu moment 5: Shape descriptor sensitive to orientation.";
            case Attribute::HU_MOMENT_6: return "Hu moment 6: Invariant to scale, translation and rotation.";
            case Attribute::HU_MOMENT_7: return "Hu moment 7: Captures asymmetry and fine shape variation.";

            // Derived from moments
            case Attribute::INERTIA: return "Inertia: Resistance to rotation; second-order moment around centroid.";
            case Attribute::COMPACTNESS: return "Compactness: Shape compactness; often area divided by perimeter squared.";
            case Attribute::ECCENTRICITY: return "Eccentricity: Elongation of the object; ratio of major to minor axis.";
            case Attribute::LENGTH_MAJOR_AXIS: return "Major axis length: Length of the longest diameter of the shape.";
            case Attribute::LENGTH_MINOR_AXIS: return "Minor axis length: Length of the shortest diameter of the shape.";
            case Attribute::AXIS_ORIENTATION: return "Axis orientation: Angle of the major axis with respect to the horizontal.";

            // tree topology 
            case Attribute::HEIGHT_NODE: return "Height: Longest distance to any descendant (leaf) node in the tree.";
            case Attribute::DEPTH_NODE: return "Depth: Distance from the node to the root of the tree.";
            case Attribute::IS_LEAF_NODE: return "Is leaf: Whether the node has no children.";
            case Attribute::IS_ROOT_NODE: return "Is root: Whether the node is the root of the component tree.";
            case Attribute::NUM_CHILDREN_NODE: return "Number of children: Direct child nodes of the current node.";
            case Attribute::NUM_SIBLINGS_NODE: return "Number of siblings: Other children of the node’s parent.";
            case Attribute::NUM_DESCENDANTS_NODE: return "Number of descendants: All nodes in the subtree rooted at this node.";
            case Attribute::NUM_LEAF_DESCENDANTS_NODE: return "Number of leaf descendants: Leaf nodes in the subtree rooted at this node.";
            case Attribute::LEAF_RATIO_NODE: return "Leaf ratio: Ratio of leaf descendants to all descendants.";
            case Attribute::BALANCE_NODE: return "Balance: Difference between maximum and minimum subtree height among children.";
            case Attribute::AVG_CHILD_HEIGHT_NODE: return "Average child height: Mean height of all direct children of the node.";

            default:
                return "Unknown attribute.";
    	}
	}

    static std::optional<Attribute> parse(const std::string& str) {
        static const std::unordered_map<std::string, Attribute> lookup = [] {
            std::unordered_map<std::string, Attribute> m;
            for (int i = 0; i <= static_cast<int>(AVG_CHILD_HEIGHT_NODE); ++i) {
                auto id = static_cast<Attribute>(i);
                m[toString(id)] = id;
            }
            return m;
        }();

        auto it = lookup.find(str);
        if (it != lookup.end()) return it->second;
        return std::nullopt;
    }
};




class AttributeComputer {
	public:
		virtual ~AttributeComputer() = default;
	
		/// Executa a computação dos atributos produzidos por essa classe
		virtual void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources = {}) const {
			compute(tree, buffer, attrNames, this->attributes(), dependencySources);
		}

		/// Executa a computação somente dos atributos solicitados
		virtual void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<Attribute>& requestedAttributes, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources = {}) const = 0;
	
		/// Atributos produzidos
		virtual std::vector<Attribute> attributes() const = 0;
	
		/// Atributos requeridos para o cálculo (apenas metadado)
		virtual std::vector<AttributeOrGroup> requiredAttributes() const { return {}; }
	
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



	static float* computerTreeTopologyAttributes(MorphologicalTreePtr tree){
		const int n = tree->getNumNodes();
		AttributeNames attrNames = AttributeNames::fromGroup(AttributeGroup::TREE_TOPOLOGY, n);
		float* buffer = new float[n * attrNames.NUM_ATTRIBUTES];

		using enum Attribute;
		AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
			[&](NodeMTPtr node) {
				int idx = node->getIndex();
				int parentDepth = node->getParent() ? buffer[attrNames.linearIndex(node->getParent()->getIndex(), DEPTH_NODE)] : 0;
	
				buffer[attrNames.linearIndex(idx, HEIGHT_NODE)] = 0.0f; // altura
				buffer[attrNames.linearIndex(idx, DEPTH_NODE)] = node->getParent() ? parentDepth + 1 : 0; // profundidade
				buffer[attrNames.linearIndex(idx, IS_LEAF_NODE)] = node->getChildren().empty() ? 1.0f : 0.0f; // é folha
				buffer[attrNames.linearIndex(idx, IS_ROOT_NODE)] = node->getParent() == nullptr ? 1.0f : 0.0f; // é raiz
				buffer[attrNames.linearIndex(idx, NUM_CHILDREN_NODE)] = node->getChildren().size();
				buffer[attrNames.linearIndex(idx, NUM_SIBLINGS_NODE)] = node->getParent() ? node->getParent()->getChildren().size() - 1 : 0;
				buffer[attrNames.linearIndex(idx, NUM_DESCENDANTS_NODE)] = 0.0f;
				buffer[attrNames.linearIndex(idx, NUM_LEAF_DESCENDANTS_NODE)] = buffer[attrNames.linearIndex(idx, IS_LEAF_NODE)];
				buffer[attrNames.linearIndex(idx, LEAF_RATIO_NODE)] = 0.0f;
				buffer[attrNames.linearIndex(idx, BALANCE_NODE)] = 0.0f;
				buffer[attrNames.linearIndex(idx, AVG_CHILD_HEIGHT_NODE)] = 0.0f;
			},
			[&](NodeMTPtr parent, NodeMTPtr child) {
				int pIdx = parent->getIndex();
				int cIdx = child->getIndex();

				// Acumulando descendentes
				buffer[attrNames.linearIndex(pIdx, NUM_DESCENDANTS_NODE)] += buffer[attrNames.linearIndex(cIdx, NUM_DESCENDANTS_NODE)] + 1;
				buffer[attrNames.linearIndex(pIdx, NUM_LEAF_DESCENDANTS_NODE)] += buffer[attrNames.linearIndex(cIdx, NUM_LEAF_DESCENDANTS_NODE)];

				// Altura
				float childHeight = buffer[attrNames.linearIndex(cIdx, HEIGHT_NODE)];
				float& parentHeight = buffer[attrNames.linearIndex(pIdx, HEIGHT_NODE)];
				parentHeight = std::max(parentHeight, childHeight + 1);

				// Balanceamento
				float& minH = buffer[attrNames.linearIndex(pIdx, BALANCE_NODE)]; // usado como mínimo temporário
				float& sumH = buffer[attrNames.linearIndex(pIdx, AVG_CHILD_HEIGHT_NODE)];
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
				float desc = buffer[attrNames.linearIndex(idx, NUM_DESCENDANTS_NODE)];
				float folhas = buffer[attrNames.linearIndex(idx, NUM_LEAF_DESCENDANTS_NODE)];
	
				// Razão folhas/descendentes
				buffer[attrNames.linearIndex(idx, LEAF_RATIO_NODE)] = desc > 0.0f ? folhas / (desc + 1.0f) : 1.0f;
	
				// Balanceamento e média
				if (!node->getChildren().empty()) {
					float alturaMax = buffer[attrNames.linearIndex(idx, HEIGHT_NODE)];
					float alturaMin = buffer[attrNames.linearIndex(idx, BALANCE_NODE)];
					buffer[attrNames.linearIndex(idx, BALANCE_NODE)] = alturaMax - alturaMin;
	
					buffer[attrNames.linearIndex(idx, AVG_CHILD_HEIGHT_NODE)] /= node->getChildren().size();
				}

			}
		);
	
		return buffer;
	}


	static float* computerAttributeByIndex(MorphologicalTreePtr tree, Attribute attribute){
		const int n = tree->getNumNodes();
		float* attr = new float[n];
		auto [attrNames, buffer] = AttributeComputedIncrementally::computerBasicAttributes(tree);
		for(int idx = 0; idx < n; idx++){
			attr[idx] = buffer[ attrNames.linearIndex(idx, attribute) ];
		}
		delete[] buffer;

		return attr;
	}

	static std::pair<AttributeNames, float*> computerBasicAttributes(MorphologicalTreePtr tree){
	    
		int n = tree->getNumNodes();
		AttributeNames attrNames = AttributeNames::fromGroup(AttributeGroup::ALL, n);
		float* buffer = new float[n * attrNames.NUM_ATTRIBUTES];

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
		
		auto indexOf = [&](int idxNode, Attribute attribute) {
			return attrNames.linearIndex(idxNode, attribute);
		};

		//computação dos atributos: area, volume, gray level, mean of gray level, variance of gray level, standard deviation gray level, Box width, Box height, rectangularity, ratio (Box width, Box height) e momentos geometricos 
	    AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
						[&indexOf, &buffer, n,  &xmax, &ymax, &xmin, &ymin, &sumX, &sumY, &sumGrayLevelSquare, numCols, numRows](NodeMTPtr node) -> void {
							
							buffer[ indexOf(node->getIndex(), AREA) ] = node->getCNPs().size(); //area
							buffer[ indexOf(node->getIndex(), VOLUME) ] = node->getCNPs().size() * node->getLevel(); //volume =>  \sum{ f }
							buffer[ indexOf(node->getIndex(), LEVEL) ] = node->getLevel(); //level
							if (node->isMaxtreeNode()) {
								buffer[ indexOf(node->getIndex(), DYNAMICS)] = std::numeric_limits<int>::lowest(); // Procuramos máximo para max-tree
							} else {
								buffer[ indexOf(node->getIndex(), DYNAMICS)] = std::numeric_limits<int>::max(); // Procuramos mínimo para min-tree
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
						[&indexOf, &attrNames, &buffer, n, &xmax, &ymax, &xmin, &ymin, &sumX, &sumY, &sumGrayLevelSquare](NodeMTPtr parent, NodeMTPtr child) -> void {
							buffer[ indexOf(parent->getIndex(), AREA) ] += buffer[ indexOf(child->getIndex(), AREA) ]; //area
							buffer[ indexOf(parent->getIndex(), VOLUME) ] += buffer[ indexOf(child->getIndex(), VOLUME) ]; //volume
							if (parent->isMaxtreeNode()) {
								buffer[ indexOf(parent->getIndex(), DYNAMICS) ] = std::max(buffer[ indexOf(parent->getIndex(), DYNAMICS)], buffer[ indexOf(child->getIndex(), DYNAMICS)] );
							} else {
								buffer[ indexOf(parent->getIndex(), DYNAMICS) ] = std::min(buffer[ indexOf(parent->getIndex(), DYNAMICS)], buffer[ indexOf(child->getIndex(), DYNAMICS)] );
							}

							sumGrayLevelSquare[parent->getIndex()] += sumGrayLevelSquare[child->getIndex()]; //computando: \sum{ f^2 }

							ymax[parent->getIndex()] = std::max(ymax[parent->getIndex()], ymax[child->getIndex()]);
							xmax[parent->getIndex()] = std::max(xmax[parent->getIndex()], xmax[child->getIndex()]);
							ymin[parent->getIndex()] = std::min(ymin[parent->getIndex()], ymin[child->getIndex()]);
							xmin[parent->getIndex()] = std::min(xmin[parent->getIndex()], xmin[child->getIndex()]);
		
							sumX[parent->getIndex()] += sumX[child->getIndex()];
							sumY[parent->getIndex()] += sumY[child->getIndex()];
							
						},
						[&indexOf, &attrNames, &buffer, n, &xmax, &ymax, &xmin, &ymin, &sumGrayLevelSquare](NodeMTPtr node) -> void {
							
							float area = buffer[ indexOf(node->getIndex(), AREA) ];
							float volume = buffer[ indexOf(node->getIndex(), VOLUME) ];
							float width = xmax[node->getIndex()] - xmin[node->getIndex()] + 1;	
							float height = ymax[node->getIndex()] - ymin[node->getIndex()] + 1;	
							
							float meanGrayLevel = volume / area; //mean graylevel - // E(f)
							double meanGrayLevelSquare = sumGrayLevelSquare[node->getIndex()] / area; // E(f^2)
							float var = meanGrayLevelSquare - (meanGrayLevel * meanGrayLevel); //variance: E(f^2) - E(f)^2
							buffer[ indexOf(node->getIndex(), VARIANCE_LEVEL) ] = var > 0? var : 0; //variance
							
							/*
							if (buffer[indexOf(node->getIndex(), VARIANCE_LEVEL)] >= 0.0f) {
								buffer[indexOf(node->getIndex(), STANDARD_DEVIATION)] = std::sqrt(buffer[indexOf(node->getIndex(), VARIANCE_LEVEL)]);
							} else {
								buffer[indexOf(node->getIndex(), STANDARD_DEVIATION)] = 0.0f; // Se a variância for negativa, definir desvio padrão como 0
							}*/
							
							buffer[indexOf(node->getIndex(), MEAN_LEVEL)] = meanGrayLevel;
							buffer[indexOf(node->getIndex(), BOX_WIDTH)] = width;
							buffer[indexOf(node->getIndex(), BOX_HEIGHT)] = height;
							buffer[indexOf(node->getIndex(), RECTANGULARITY)] = area / (width * height);
							buffer[indexOf(node->getIndex(), RATIO_WH)] = std::max(width, height) / std::min(width, height);
							
							if (node->getChildren().empty()) {
								buffer[ indexOf(node->getIndex(), DYNAMICS)] = 0; // Folhas têm dinâmica 0
							} else {
								buffer[ indexOf(node->getIndex(), DYNAMICS)] = std::abs( node->getLevel() - buffer[ indexOf(node->getIndex(), DYNAMICS)] ) + 1;
							}


		});

		

		//Computação dos momentos centrais e momentos de Hu
		AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
			[&indexOf, &buffer, n,  &sumX, &sumY, numCols](NodeMTPtr node) -> void {				
				// Inicialização dos momentos centrais
				buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_20)] = 0.0f;
				buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_02)] = 0.0f;
				buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_11)] = 0.0f;
				buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_30)] = 0.0f;
				buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_03)] = 0.0f;
				buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_21)] = 0.0f;
				buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_12)] = 0.0f;

				// Cálculo do centroide
				float xCentroid = sumX[node->getIndex()] / buffer[indexOf(node->getIndex(), AREA)];
				float yCentroid = sumY[node->getIndex()] / buffer[indexOf(node->getIndex(), AREA)];

				for (int p : node->getCNPs()) {
					int x = p % numCols;
					int y = p / numCols;
					float dx = x - xCentroid;
					float dy = y - yCentroid;

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
			[&indexOf, &buffer, n](NodeMTPtr parent, NodeMTPtr child) -> void {
				buffer[indexOf(parent->getIndex(), CENTRAL_MOMENT_20)] += buffer[indexOf(child->getIndex(), CENTRAL_MOMENT_20)];
				buffer[indexOf(parent->getIndex(), CENTRAL_MOMENT_02)] += buffer[indexOf(child->getIndex(), CENTRAL_MOMENT_02)];
				buffer[indexOf(parent->getIndex(), CENTRAL_MOMENT_11)] += buffer[indexOf(child->getIndex(), CENTRAL_MOMENT_11)];
				buffer[indexOf(parent->getIndex(), CENTRAL_MOMENT_30)] += buffer[indexOf(child->getIndex(), CENTRAL_MOMENT_30)];
				buffer[indexOf(parent->getIndex(), CENTRAL_MOMENT_03)] += buffer[indexOf(child->getIndex(), CENTRAL_MOMENT_03)];
				buffer[indexOf(parent->getIndex(), CENTRAL_MOMENT_21)] += buffer[indexOf(child->getIndex(), CENTRAL_MOMENT_21)];
				buffer[indexOf(parent->getIndex(), CENTRAL_MOMENT_12)] += buffer[indexOf(child->getIndex(), CENTRAL_MOMENT_12)];			
			},
			[&indexOf, &buffer, n](NodeMTPtr node) -> void {
				int idx = node->getIndex();
				float area = buffer[indexOf(idx, AREA)];
				auto normMoment = [area](float moment, int p, int q){ 
					return moment / std::pow( area, (p + q + 2.0) / 2.0); 
				}; //função para normalizacao dos momentos				
				

				//Momentos centrais
				float mu20 = buffer[indexOf(idx, CENTRAL_MOMENT_20)];
				float mu02 = buffer[indexOf(idx, CENTRAL_MOMENT_02)];
				float mu11 = buffer[indexOf(idx, CENTRAL_MOMENT_11)];
				float mu30 = buffer[indexOf(idx, CENTRAL_MOMENT_30)];
				float mu03 = buffer[indexOf(idx, CENTRAL_MOMENT_03)];
				float mu21 = buffer[indexOf(idx, CENTRAL_MOMENT_21)];
				float mu12 = buffer[indexOf(idx, CENTRAL_MOMENT_12)];
					
				float discriminant = std::pow(mu20 - mu02, 2.0f) + 4.0f * std::pow(mu11, 2.0f);

					
				// Verificar se o denominador é zero antes de calcular atan2 para evitar divisão por zero
				if (mu20 != mu02 || mu11 != 0) {
					float radians = 0.5 * std::atan2(2 * mu11, mu20 - mu02);// orientação em radianos
					float degrees = radians * (180.0 / M_PI); // Converter para graus
					if (degrees < 0) { // Ajustar para o intervalo [0, 360] graus
						degrees += 360.0;
					}
					buffer[indexOf(idx, AXIS_ORIENTATION)] = degrees; // Armazenar a orientação em graus no intervalo [0, 360]
				} else {
					buffer[indexOf(idx, AXIS_ORIENTATION)] = 0.0; // Se não for possível calcular a orientação, definir um valor padrão
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
					buffer[indexOf(idx, LENGTH_MAJOR_AXIS)] = std::sqrt((2 * a1) / area); // length major axis
				} else {
					buffer[indexOf(idx, LENGTH_MAJOR_AXIS)] = 0.0; // Definir valor padrão
				}

				if (a2 > 0) {
					buffer[indexOf(idx, LENGTH_MINOR_AXIS)] = std::sqrt((2 * a2) / area); // length minor axis
				} else {
					buffer[indexOf(idx, LENGTH_MINOR_AXIS)] = 0.0; // Definir valor padrão
				}

				// Verificar se a2 é diferente de zero antes de calcular a excentricidade
				buffer[indexOf(idx, ECCENTRICITY)] = (std::abs(a2) > std::numeric_limits<float>::epsilon()) ? a1 / a2 : a1 / 0.1; // eccentricity

				// Verificar se mu20 + mu02 é diferente de zero antes de calcular a compacidade
				if ((mu20 + mu02) > std::numeric_limits<float>::epsilon()) {
					buffer[indexOf(idx, COMPACTNESS)]  = (1.0 / (2 * PI)) * (area / (mu20 + mu02)); // compactness
				} else {
					buffer[indexOf(idx, COMPACTNESS)]  = 0.0; // Definir valor padrão
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
				buffer[indexOf(idx, HU_MOMENT_1)] = eta20 + eta02; // primeiro momento de Hu => inertia
				buffer[indexOf(idx, HU_MOMENT_2)]  = std::pow(eta20 - eta02, 2) + 4 * std::pow(eta11, 2);
				buffer[indexOf(idx, HU_MOMENT_3)]  = std::pow(eta30 - 3 * eta12, 2) + std::pow(3 * eta21 - eta03, 2);
				buffer[indexOf(idx, HU_MOMENT_4)]  = std::pow(eta30 + eta12, 2) + std::pow(eta21 + eta03, 2);
				
				buffer[indexOf(idx, HU_MOMENT_5)] = 
					(eta30 - 3 * eta12) * (eta30 + eta12) * (std::pow(eta30 + eta12, 2) - 3 * std::pow(eta21 + eta03, 2)) +
					(3 * eta21 - eta03) * (eta21 + eta03) * (3 * std::pow(eta30 + eta12, 2) - std::pow(eta21 + eta03, 2));
				
					buffer[indexOf(idx, HU_MOMENT_6)] = 
					(eta20 - eta02) * (std::pow(eta30 + eta12, 2) - std::pow(eta21 + eta03, 2)) + 
					4 * eta11 * (eta30 + eta12) * (eta21 + eta03);
				
					buffer[indexOf(idx, HU_MOMENT_7)] = 
					(3 * eta21 - eta03) * (eta30 + eta12) * (std::pow(eta30 + eta12, 2) - 3 * std::pow(eta21 + eta03, 2)) -
					(eta30 - 3 * eta12) * (eta21 + eta03) * (3 * std::pow(eta30 + eta12, 2) - std::pow(eta21 + eta03, 2));

				
		});
		return std::make_pair(attrNames, buffer);
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

	static std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>> computeAttributesByComputer(MorphologicalTreePtr tree, std::shared_ptr<AttributeComputer> comp, const DependencyMap& available = {});
	
	static std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>> computeSingleAttribute(MorphologicalTreePtr tree, AttributeOrGroup attr, const DependencyMap& available = {});

	static std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>> computeAttributes(MorphologicalTreePtr tree, const std::vector<AttributeOrGroup>& attributes,const DependencyMap& providedDependencies={});
};



class AreaComputer : public AttributeComputer {
	public:
		
		std::vector<Attribute> attributes() const override {
			return {AREA};
		}
		void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<Attribute>& requestedAttributes, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources= {}) const override {
			std::cout << "\n==== AttributeComputer: Computing AREA" << std::endl;
			auto indexOf = [&](int idx) {
				return attrNames->linearIndex(idx, AREA);
			};
	
			for(NodeMTPtr node: tree->getIndexNode()){
				buffer[indexOf(node->getIndex())] = node->getAreaCC();
			}
			/* //mesmo que: getAreaCC()
			AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
				[&](NodeMTPtr node) {
					buffer[indexOf(node->getIndex())] = node->getCNPs().size();
				},
				[&](NodeMTPtr parent, NodeMTPtr child) {
					buffer[indexOf(parent->getIndex())] += buffer[indexOf(child->getIndex())];
				},
				[](NodeMTPtr node) {}
			);*/
		}
};


class VolumeComputer : public AttributeComputer {
	public:
		std::vector<Attribute> attributes() const override {
			return {VOLUME};
		}
		void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<Attribute>& requestedAttributes, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources= {}) const override {
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

class DynamicsComputer : public AttributeComputer {
	public:
		
		std::vector<Attribute> attributes() const override {
			return {DYNAMICS};
		}
		void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<Attribute>& requestedAttributes, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources= {}) const override {
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


class GrayLevelStatsComputer : public AttributeComputer {
	public:
		
	
		std::vector<Attribute> attributes() const override {
			return { LEVEL, MEAN_LEVEL, VARIANCE_LEVEL };
		}
	
		std::vector<AttributeOrGroup> requiredAttributes() const override {
			return { VOLUME };
		}
	
		void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<Attribute>& requestedAttributes, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources = {}) const override {
			
			std::cout << "\n==== AttributeComputer: Computing GrayLevelStatsComputer " << std::endl;
	
			auto indexOfMean = [&](int idx) {
				return attrNames->linearIndex(idx, MEAN_LEVEL);
			};
			auto indexOfLevel = [&](int idx) {
				return attrNames->linearIndex(idx, LEVEL);
			};
			auto indexOfVarianceLevel = [&](int idx) {
				return attrNames->linearIndex(idx, VARIANCE_LEVEL);
			};

			bool computeMeanLevel = std::find(requestedAttributes.begin(), requestedAttributes.end(), MEAN_LEVEL) != requestedAttributes.end();
			bool computeVarianceLevel = std::find(requestedAttributes.begin(), requestedAttributes.end(), VARIANCE_LEVEL) != requestedAttributes.end();
			bool computeLevel = std::find(requestedAttributes.begin(), requestedAttributes.end(), LEVEL) != requestedAttributes.end();
			

			auto [dependencyAttrNamesVol, bufferVol] = dependencySources[0]; //volume
			auto indexOfVol = [&](int idx) {
				return dependencyAttrNamesVol->linearIndex(idx, VOLUME);
			};
	

			std::shared_ptr<long[]> sumGrayLevelSquare = nullptr;
			if(computeVarianceLevel) {
				sumGrayLevelSquare = std::shared_ptr<long[]>(new long[tree->getNumNodes()]);
			}
			AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
				[&](NodeMTPtr node) {
					if(computeVarianceLevel)
						sumGrayLevelSquare[node->getIndex()] = node->getCNPs().size() * std::pow(node->getLevel(), 2);
					if(computeLevel)
						buffer[ indexOfLevel(node->getIndex()) ] = node->getLevel();
				},
				[&](NodeMTPtr parent, NodeMTPtr child) { 
					if(computeVarianceLevel)
						sumGrayLevelSquare[parent->getIndex()] += sumGrayLevelSquare[child->getIndex()];
				},
				[&](NodeMTPtr node) {
					float area = node->getAreaCC();
					if(computeMeanLevel)
						buffer[ indexOfMean(node->getIndex()) ] = bufferVol[indexOfVol(node->getIndex())] / area;

					if(computeVarianceLevel) {
						int idx = node->getIndex();
						
						float meanGrayLevel = bufferVol[indexOfVol(node->getIndex())] / area; //mean graylevel = E(f)
						double meanGrayLevelSquare = sumGrayLevelSquare[idx] / area; // E(f^2)
						float var = meanGrayLevelSquare - (meanGrayLevel * meanGrayLevel); //variance: E(f^2) - E(f)^2
						buffer[indexOfVarianceLevel(idx)] = var > 0.0f ? var : 0.0f; //variance
					}

					
				}
			);
		}
};



class BoundingBoxComputer : public AttributeComputer {
	public:
		
	
		std::vector<Attribute> attributes() const override {
			return { BOX_WIDTH, BOX_HEIGHT, RECTANGULARITY, RATIO_WH, BOX_COL_MIN, BOX_COL_MAX, BOX_ROW_MIN, BOX_ROW_MAX, DIAGONAL_LENGTH };
		}

		
		void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<Attribute>& requestedAttributes, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources = {}) const override {
			
			std::cout << "\n==== AttributeComputer: Computing BOUNDING_BOX group" << std::endl;
	
			auto indexOfWidth  = [&](int idx) { return attrNames->linearIndex(idx, BOX_WIDTH); };
			auto indexOfHeight = [&](int idx) { return attrNames->linearIndex(idx, BOX_HEIGHT); };
			auto indexOfRectangularity = [&](int idx) { return attrNames->linearIndex(idx, RECTANGULARITY); };
			auto indexOfRatioWH = [&](int idx) { return attrNames->linearIndex(idx, RATIO_WH); };
			auto indexOfColMin = [&](int idx) { return attrNames->linearIndex(idx, BOX_COL_MIN); };
			auto indexOfColMax = [&](int idx) { return attrNames->linearIndex(idx, BOX_COL_MAX); };
			auto indexOfRowMin = [&](int idx) { return attrNames->linearIndex(idx, BOX_ROW_MIN); };
			auto indexOfRowMax = [&](int idx) { return attrNames->linearIndex(idx, BOX_ROW_MAX); };
			auto indexOfDiagonalLength = [&](int idx) { return attrNames->linearIndex(idx, DIAGONAL_LENGTH); };

			bool computeWidth  = std::find(requestedAttributes.begin(), requestedAttributes.end(), BOX_WIDTH)  != requestedAttributes.end();
			bool computeHeight = std::find(requestedAttributes.begin(), requestedAttributes.end(), BOX_HEIGHT) != requestedAttributes.end();
			bool computeRectangularity = std::find(requestedAttributes.begin(), requestedAttributes.end(), RECTANGULARITY) != requestedAttributes.end();
			bool computeRatioWH = std::find(requestedAttributes.begin(), requestedAttributes.end(), RATIO_WH) != requestedAttributes.end();
			bool computeColMin = std::find(requestedAttributes.begin(), requestedAttributes.end(), BOX_COL_MIN) != requestedAttributes.end();
			bool computeColMax = std::find(requestedAttributes.begin(), requestedAttributes.end(), BOX_COL_MAX) != requestedAttributes.end();
			bool computeRowMin = std::find(requestedAttributes.begin(), requestedAttributes.end(), BOX_ROW_MIN) != requestedAttributes.end();
			bool computeRowMax = std::find(requestedAttributes.begin(), requestedAttributes.end(), BOX_ROW_MAX) != requestedAttributes.end();
			bool computeDiagonalLength = std::find(requestedAttributes.begin(), requestedAttributes.end(), DIAGONAL_LENGTH) != requestedAttributes.end();

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
					if(computeWidth)
						buffer[indexOfWidth(idx)]  = xmax[idx] - xmin[idx] + 1;
					if(computeHeight)
						buffer[indexOfHeight(idx)] = ymax[idx] - ymin[idx] + 1;

					if(computeRectangularity) {
						float area = node->getAreaCC();
						float width = xmax[idx] - xmin[idx] + 1;
						float height = ymax[idx] - ymin[idx] + 1;
						float denom = width * height;
						buffer[indexOfRectangularity(idx)] = (denom > 0.0f) ? (area / denom) : 0.0f;
					}
					if(computeRatioWH) {
						float width  = xmax[idx] - xmin[idx] + 1;
						float height = ymax[idx] - ymin[idx] + 1;
						if (width > 0 && height > 0) {
							buffer[indexOfRatioWH(idx)] = std::max(width, height) / std::min(width, height);
						} else {
							buffer[indexOfRatioWH(idx)] = 0.0f;
						}
					}
					if(computeColMin)
						buffer[indexOfColMin(idx)]  = xmin[idx];
					if(computeColMax)
						buffer[indexOfColMax(idx)]  = xmax[idx];
					if(computeRowMin)
						buffer[indexOfRowMin(idx)]  = ymin[idx];
					if(computeRowMax)
						buffer[indexOfRowMax(idx)]  = ymax[idx];
					if(computeDiagonalLength) {
						float width  = xmax[idx] - xmin[idx] + 1;
						float height = ymax[idx] - ymin[idx] + 1;
						buffer[indexOfDiagonalLength(idx)] = std::sqrt(width*width + height*height);
					}
				}
			);
		}
};

class CentralMomentsComputer : public AttributeComputer {
	public:
		

		std::vector<Attribute> attributes() const override {
			return {CENTRAL_MOMENT_20,
					CENTRAL_MOMENT_02,
					CENTRAL_MOMENT_11,
					CENTRAL_MOMENT_30,
					CENTRAL_MOMENT_03,
					CENTRAL_MOMENT_21,
					CENTRAL_MOMENT_12};
		}
		
		void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<Attribute>& requestedAttributes, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources= {}) const override {
			std::cout << "\n==== AttributeComputer: Computing CENTRAL_MOMENT group" << std::endl;
			//momentos geometricos para calcular o centroide
			int numCols = tree->getNumColsOfImage();
			int n = tree->getNumColsOfImage() * tree->getNumRowsOfImage();
			std::unique_ptr<long int[]> sumX(new long int[n]);//sum x
			std::unique_ptr<long int[]> sumY(new long int[n]);//sum y
			
			auto indexOf = [&](int idx, Attribute attr) { return attrNames->linearIndex(idx, attr); };
			
			
			bool computeMu20 = std::find(requestedAttributes.begin(), requestedAttributes.end(), CENTRAL_MOMENT_20) != requestedAttributes.end();
			bool computeMu02 = std::find(requestedAttributes.begin(), requestedAttributes.end(), CENTRAL_MOMENT_02) != requestedAttributes.end();
			bool computeMu11 = std::find(requestedAttributes.begin(), requestedAttributes.end(), CENTRAL_MOMENT_11) != requestedAttributes.end();
			bool computeMu30 = std::find(requestedAttributes.begin(), requestedAttributes.end(), CENTRAL_MOMENT_30) != requestedAttributes.end();
			bool computeMu03 = std::find(requestedAttributes.begin(), requestedAttributes.end(), CENTRAL_MOMENT_03) != requestedAttributes.end();
			bool computeMu21 = std::find(requestedAttributes.begin(), requestedAttributes.end(), CENTRAL_MOMENT_21) != requestedAttributes.end();
			bool computeMu12 = std::find(requestedAttributes.begin(), requestedAttributes.end(), CENTRAL_MOMENT_12) != requestedAttributes.end();
			
			//computa sumX e sumY para calcular os centroides
			AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
				[&](NodeMTPtr node) -> void {
					sumX[node->getIndex()] = 0;
					sumY[node->getIndex()] = 0;
					for(int p: node->getCNPs()) {
						auto [py, px] = ImageUtils::to2D(p, numCols); 
						sumX[node->getIndex()] += px;
						sumY[node->getIndex()] += py;
					}
				},
				[&](NodeMTPtr parent, NodeMTPtr child) -> void {
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
					if(computeMu20)
						buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_20)] = 0.0f;
					if(computeMu02)
						buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_02)] = 0.0f;
					if(computeMu11)
						buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_11)] = 0.0f;
					if(computeMu30)
						buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_30)] = 0.0f;
					if(computeMu03)
						buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_03)] = 0.0f;
					if(computeMu21)
						buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_21)] = 0.0f;
					if(computeMu12)
						buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_12)] = 0.0f;
					
	
					// Cálculo do centroide
					float area = node->getAreaCC();
					float xCentroid = sumX[node->getIndex()] / area;
					float yCentroid = sumY[node->getIndex()] / area;
	
					for (int p : node->getCNPs()) {
						auto [py, px] = ImageUtils::to2D(p, numCols); 
						float dx = px - xCentroid;
						float dy = py - yCentroid;
	
						// Momentos centrais de segunda ordem
						if(computeMu20)
							buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_20)] += std::pow(dx, 2);
						if(computeMu02)
							buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_02)] += std::pow(dy, 2);
						if(computeMu11)
							buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_11)] += dx * dy;
	
						// Momentos centrais de terceira ordem
						if(computeMu30)
							buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_30)] += std::pow(dx, 3);
						if(computeMu03)
							buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_03)] += std::pow(dy, 3);
						if(computeMu21)
							buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_21)] += std::pow(dx, 2) * dy;
						if(computeMu12)
							buffer[indexOf(node->getIndex(), CENTRAL_MOMENT_12)] += dx * std::pow(dy, 2);
					}
	
				},
				[&](NodeMTPtr parent, NodeMTPtr child) -> void {
					// Acumulação dos momentos centrais
					if(computeMu20)
						buffer[indexOf(parent->getIndex(), CENTRAL_MOMENT_20)] += buffer[indexOf(child->getIndex(), CENTRAL_MOMENT_20)];
					if(computeMu02)
						buffer[indexOf(parent->getIndex(), CENTRAL_MOMENT_02)] += buffer[indexOf(child->getIndex(), CENTRAL_MOMENT_02)];
					if(computeMu11)	
						buffer[indexOf(parent->getIndex(), CENTRAL_MOMENT_11)] += buffer[indexOf(child->getIndex(), CENTRAL_MOMENT_11)];
					if(computeMu30)
						buffer[indexOf(parent->getIndex(), CENTRAL_MOMENT_30)] += buffer[indexOf(child->getIndex(), CENTRAL_MOMENT_30)];
					if(computeMu03)
						buffer[indexOf(parent->getIndex(), CENTRAL_MOMENT_03)] += buffer[indexOf(child->getIndex(), CENTRAL_MOMENT_03)];
					if(computeMu21)
						buffer[indexOf(parent->getIndex(), CENTRAL_MOMENT_21)] += buffer[indexOf(child->getIndex(), CENTRAL_MOMENT_21)];
					if(computeMu12)
						buffer[indexOf(parent->getIndex(), CENTRAL_MOMENT_12)] += buffer[indexOf(child->getIndex(), CENTRAL_MOMENT_12)];			
				},
				[](NodeMTPtr node) -> void {
					// Não é necessário fazer nada aqui
				}
		);
	}
};

class MomentBasedAttributeComputer : public AttributeComputer {
	public:
		
	
		std::vector<Attribute> attributes() const override {
			return { LENGTH_MAJOR_AXIS, LENGTH_MINOR_AXIS, ECCENTRICITY, COMPACTNESS, AXIS_ORIENTATION, INERTIA};
		}
	
		std::vector<AttributeOrGroup> requiredAttributes() const override {
			return { AttributeGroup::CENTRAL_MOMENTS };
		}
	
		void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<Attribute>& requestedAttributes, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources = {}) const override {
			
			std::cout << "\n==== AttributeComputer: Computing MOMENT_BASED group" << std::endl;
	
			auto indexOfMajorAxis = [&](int idx) { return attrNames->linearIndex(idx, LENGTH_MAJOR_AXIS); };
			auto indexOfMinorAxis = [&](int idx) { return attrNames->linearIndex(idx, LENGTH_MINOR_AXIS); };
			auto indexOfEccentricity = [&](int idx) { return attrNames->linearIndex(idx, ECCENTRICITY); };
			auto indexOfCompactness = [&](int idx) { return attrNames->linearIndex(idx, COMPACTNESS); };
			auto indexOfAxisOrientation = [&](int idx) { return attrNames->linearIndex(idx, AXIS_ORIENTATION); };
			auto indexOfInertia = [&](int idx) { return attrNames->linearIndex(idx, INERTIA); };

			bool computeMajorAxis  = std::find(requestedAttributes.begin(), requestedAttributes.end(), LENGTH_MAJOR_AXIS)  != requestedAttributes.end();
			bool computeMinorAxis = std::find(requestedAttributes.begin(), requestedAttributes.end(), LENGTH_MINOR_AXIS) != requestedAttributes.end();
			bool computeEccentricity = std::find(requestedAttributes.begin(), requestedAttributes.end(), ECCENTRICITY) != requestedAttributes.end();
			bool computeCompactness = std::find(requestedAttributes.begin(), requestedAttributes.end(), COMPACTNESS) != requestedAttributes.end();
			bool computeAxisOrientation = std::find(requestedAttributes.begin(), requestedAttributes.end(), AXIS_ORIENTATION) != requestedAttributes.end();
			bool computeInertia = std::find(requestedAttributes.begin(), requestedAttributes.end(), INERTIA) != requestedAttributes.end();
	
			auto [namesMom, bufMom] = dependencySources[0];
			auto indexMu20 = [&](int idx) { return namesMom->linearIndex(idx, CENTRAL_MOMENT_20); };
			auto indexMu02 = [&](int idx) { return namesMom->linearIndex(idx, CENTRAL_MOMENT_02); };
			auto indexMu11 = [&](int idx) { return namesMom->linearIndex(idx, CENTRAL_MOMENT_11); };
			
			
			AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
				[&](NodeMTPtr) {},
				[&](NodeMTPtr, NodeMTPtr) {},
				[&](NodeMTPtr node) {
					int idx = node->getIndex();
					float mu20 = bufMom[indexMu20(idx)];
					float mu02 = bufMom[indexMu02(idx)];
					float mu11 = bufMom[indexMu11(idx)];
					float area = node->getAreaCC();
	
					float discriminant = std::pow(mu20 - mu02, 2.0f) + 4.0f * std::pow(mu11, 2.0f);
					discriminant = std::max(discriminant, 0.0f);
					float lambda1 = mu20 + mu02 + std::sqrt(discriminant);  // maior autovalor
					float lambda2 = mu20 + mu02 - std::sqrt(discriminant);  // menor autovalor
	
					if(computeMajorAxis){
						if (area > 0.0f && lambda1 > 0.0f) {
							buffer[indexOfMajorAxis(idx)] = std::sqrt((2.0f * lambda1) / area);
						} else {
							buffer[indexOfMajorAxis(idx)] = 0.0f;
						}
					}
					if(computeMinorAxis){
						if (area > 0.0f && lambda2 > 0.0f) {
							buffer[indexOfMinorAxis(idx)] = std::sqrt((2.0f * lambda2) / area);
						} else {
							buffer[indexOfMinorAxis(idx)] = 0.0f;
						}
					}
					if(computeEccentricity){	
						if (std::abs(lambda2) > std::numeric_limits<float>::epsilon()) {
							buffer[indexOfEccentricity(idx)] = lambda1 / lambda2;
						} else {
							buffer[indexOfEccentricity(idx)] = lambda1 / 0.1f; // fallback para evitar divisão por zero
						}
					}
					if(computeCompactness){
						float denom = mu20 + mu02;
						if (denom > std::numeric_limits<float>::epsilon()) {
							buffer[indexOfCompactness(idx)] = (1.0f / (2.0f * static_cast<float>(M_PI))) * (area / denom);
						} else {
							buffer[indexOfCompactness(idx)] = 0.0f;
						}
					}
					if(computeAxisOrientation){
						// Verificar se o denominador é zero antes de calcular atan2 para evitar divisão por zero
						if (mu20 != mu02 || mu11 != 0) {
							float radians = 0.5 * std::atan2(2 * mu11, mu20 - mu02);// orientação em radianos
							float degrees = radians * (180.0 / M_PI); // Converter para graus
							buffer[indexOfAxisOrientation(idx)] = std::fmod(std::abs(degrees), 360.0f); ; // Armazenar a orientação no intervalo [0, 360]
						} else {
							buffer[indexOfAxisOrientation(idx)] = 0.0; // Se não for possível calcular a orientação, definir um valor padrão
						}
					}
					if(computeInertia){
						float normMu20 = mu20 / std::pow(area, 2.0f);
						float normMu02 = mu02 / std::pow(area, 2.0f);
						buffer[indexOfInertia(idx)] = normMu20 + normMu02;
					}

				}
			);
		}
};

class HuMomentsComputer : public AttributeComputer {
	public:
		std::vector<Attribute> attributes() const override {
			return {HU_MOMENT_1,
					HU_MOMENT_2,
					HU_MOMENT_3,
					HU_MOMENT_4,
					HU_MOMENT_5,
					HU_MOMENT_6,
					HU_MOMENT_7};
		}
		std::vector<AttributeOrGroup> requiredAttributes() const override{ 
			return { AttributeGroup::CENTRAL_MOMENTS };
		}

		void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<Attribute>& requestedAttributes, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources= {}) const override {
			std::cout << "\n==== AttributeComputer: Computing HU_MOMENT group" << std::endl;
			int numCols = tree->getNumColsOfImage();
			auto indexOf = [&](int idx, Attribute attr) {
				return attrNames->linearIndex(idx, attr);
			};
			
			auto [dependencyAttrNamesMu, bufferMu] = dependencySources[0]; //momentos centrais
			auto indexOfMu = [&](int idx, Attribute attr) {
				return dependencyAttrNamesMu->linearIndex(idx, attr);
			};
			
			auto normMoment = [](int area, float moment, int p, int q){ 
				return moment / std::pow( area, (p + q + 2.0) / 2.0); 
			}; //função para normalizacao dos momentos		

			bool computeHu1 = std::find(requestedAttributes.begin(), requestedAttributes.end(), HU_MOMENT_1) != requestedAttributes.end();
			bool computeHu2 = std::find(requestedAttributes.begin(), requestedAttributes.end(), HU_MOMENT_2) != requestedAttributes.end();
			bool computeHu3 = std::find(requestedAttributes.begin(), requestedAttributes.end(), HU_MOMENT_3) != requestedAttributes.end();
			bool computeHu4 = std::find(requestedAttributes.begin(), requestedAttributes.end(), HU_MOMENT_4) != requestedAttributes.end();
			bool computeHu5 = std::find(requestedAttributes.begin(), requestedAttributes.end(), HU_MOMENT_5) != requestedAttributes.end();
			bool computeHu6 = std::find(requestedAttributes.begin(), requestedAttributes.end(), HU_MOMENT_6) != requestedAttributes.end();
			bool computeHu7 = std::find(requestedAttributes.begin(), requestedAttributes.end(), HU_MOMENT_7) != requestedAttributes.end();
			
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
					int area = node->getAreaCC();
					
					// Calcular os momentos normalizados
					float eta20 = normMoment(area, mu20, 2, 0);
					float eta02 = normMoment(area, mu02, 0, 2);
					float eta11 = normMoment(area, mu11, 1, 1);
					float eta30 = normMoment(area, mu30, 3, 0);
					float eta03 = normMoment(area, mu03, 0, 3);
					float eta21 = normMoment(area, mu21, 2, 1);
					float eta12 = normMoment(area, mu12, 1, 2);

					// Cálculo dos momentos de Hu
					if(computeHu1)
						buffer[indexOf(idx, HU_MOMENT_1)] = eta20 + eta02; // primeiro momento de Hu => inertia
					if(computeHu2)
						buffer[indexOf(idx, HU_MOMENT_2)]  = std::pow(eta20 - eta02, 2) + 4 * std::pow(eta11, 2);
					if(computeHu3)
						buffer[indexOf(idx, HU_MOMENT_3)]  = std::pow(eta30 - 3 * eta12, 2) + std::pow(3 * eta21 - eta03, 2);
					if(computeHu4)
						buffer[indexOf(idx, HU_MOMENT_4)]  = std::pow(eta30 + eta12, 2) + std::pow(eta21 + eta03, 2);
					if(computeHu5)
						buffer[indexOf(idx, HU_MOMENT_5)] = (eta30 - 3 * eta12) * (eta30 + eta12) * (std::pow(eta30 + eta12, 2) - 3 * std::pow(eta21 + eta03, 2)) +
														(3 * eta21 - eta03) * (eta21 + eta03) * (3 * std::pow(eta30 + eta12, 2) - std::pow(eta21 + eta03, 2));
					if(computeHu6)
						buffer[indexOf(idx, HU_MOMENT_6)] = (eta20 - eta02) * (std::pow(eta30 + eta12, 2) - std::pow(eta21 + eta03, 2)) + 
														4 * eta11 * (eta30 + eta12) * (eta21 + eta03);
					if(computeHu7)
						buffer[indexOf(idx, HU_MOMENT_7)] = (3 * eta21 - eta03) * (eta30 + eta12) * (std::pow(eta30 + eta12, 2) - 3 * std::pow(eta21 + eta03, 2)) -
														(eta30 - 3 * eta12) * (eta21 + eta03) * (3 * std::pow(eta30 + eta12, 2) - std::pow(eta21 + eta03, 2));

				}
		);
	}
};

class TreeTopologyComputer : public AttributeComputer {
	public:
		std::vector<Attribute> attributes() const override {
			return { HEIGHT_NODE,
					DEPTH_NODE,
					IS_LEAF_NODE,
					IS_ROOT_NODE,
					NUM_CHILDREN_NODE,
					NUM_SIBLINGS_NODE,
					NUM_DESCENDANTS_NODE,
					NUM_LEAF_DESCENDANTS_NODE,
					LEAF_RATIO_NODE,
					BALANCE_NODE,
					AVG_CHILD_HEIGHT_NODE };
		}
		
		void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<Attribute>& requestedAttributes, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources= {}) const override {
			
			std::cout << "\n==== AttributeComputer: Computing STRUCTURE_TREE group" << std::endl;

			bool computeHeight = std::find(requestedAttributes.begin(), requestedAttributes.end(), HEIGHT_NODE) != requestedAttributes.end();
			bool computeDepth = std::find(requestedAttributes.begin(), requestedAttributes.end(), DEPTH_NODE) != requestedAttributes.end();
			bool computeIsLeaf = std::find(requestedAttributes.begin(), requestedAttributes.end(), IS_LEAF_NODE) != requestedAttributes.end();
			bool computeIsRoot = std::find(requestedAttributes.begin(), requestedAttributes.end(), IS_ROOT_NODE) != requestedAttributes.end();
			bool computeNumChildren = std::find(requestedAttributes.begin(), requestedAttributes.end(), NUM_CHILDREN_NODE) != requestedAttributes.end();
			bool computeNumSiblings = std::find(requestedAttributes.begin(), requestedAttributes.end(), NUM_SIBLINGS_NODE) != requestedAttributes.end();
			bool computeNumDescendants = std::find(requestedAttributes.begin(), requestedAttributes.end(), NUM_DESCENDANTS_NODE) != requestedAttributes.end();
			bool computeNumLeafDescendants = std::find(requestedAttributes.begin(), requestedAttributes.end(), NUM_LEAF_DESCENDANTS_NODE) != requestedAttributes.end();
			bool computeLeafRatio = std::find(requestedAttributes.begin(), requestedAttributes.end(), LEAF_RATIO_NODE) != requestedAttributes.end();
			bool computeBalance = std::find(requestedAttributes.begin(), requestedAttributes.end(), BALANCE_NODE) != requestedAttributes.end();
			bool computeAvgChildHeight = std::find(requestedAttributes.begin(), requestedAttributes.end(), AVG_CHILD_HEIGHT_NODE) != requestedAttributes.end();
			
			// constrói o buffer de altura se necessário
			std::shared_ptr<float[]> bufferHeight = nullptr;
			if(computeHeight) {
				bufferHeight = buffer;
			}else{
				bufferHeight = std::shared_ptr<float[]>(new float[tree->getNumNodes()]);
			}
			auto indexOfHeight = [&](int idx) { 
				if(computeHeight)
					return attrNames->linearIndex(idx, HEIGHT_NODE);
				else
					return idx;
			};

			// constrói o buffer de descendentes se necessário
			std::shared_ptr<float[]> bufferNumDesc = nullptr;
			if(computeNumDescendants) {
				bufferNumDesc = buffer;
			}else{
				bufferNumDesc = std::shared_ptr<float[]>(new float[tree->getNumNodes()]);
			}
			auto indexOfNumDescendants = [&](int idx) { 
				if(computeNumDescendants)
					return attrNames->linearIndex(idx, NUM_DESCENDANTS_NODE);
				else
					return idx;
			};

			// constrói o buffer de descendentes folhas se necessário
			std::shared_ptr<float[]> bufferNumLeafDesc = nullptr;
			if(computeNumLeafDescendants) {
				bufferNumLeafDesc = buffer;
			}else{
				bufferNumLeafDesc = std::shared_ptr<float[]>(new float[tree->getNumNodes()]);
			}
			auto indexOfNumLeafDescendants = [&](int idx) { 
				if(computeNumLeafDescendants)
					return attrNames->linearIndex(idx, NUM_LEAF_DESCENDANTS_NODE);
				else
					return idx;
			};
			

			AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
				[&](NodeMTPtr node) {
					int idx = node->getIndex();
					
					int parentDepth = node->getParent() ? bufferHeight[indexOfHeight(node->getParent()->getIndex())] : 0;
					bufferHeight[indexOfHeight(idx)] = node->getParent() ? parentDepth + 1 : 0; 
					
					bufferNumDesc[indexOfNumDescendants(idx)] = 0.0f;
					bufferNumLeafDesc[indexOfNumLeafDescendants(idx)] = node->getChildren().empty() ? 1.0f : 0.0f; // é folha

					if(computeHeight)
						buffer[attrNames->linearIndex(idx, HEIGHT_NODE)] = 0.0f; // altura
					if(computeIsLeaf)
						buffer[attrNames->linearIndex(idx, IS_LEAF_NODE)] = node->getChildren().empty() ? 1.0f : 0.0f; // é folha
					if(computeIsRoot)
						buffer[attrNames->linearIndex(idx, IS_ROOT_NODE)] = node->getParent() == nullptr ? 1.0f : 0.0f; // é raiz
					if(computeNumChildren)
						buffer[attrNames->linearIndex(idx, NUM_CHILDREN_NODE)] = node->getChildren().size();
					if(computeNumSiblings)
						buffer[attrNames->linearIndex(idx, NUM_SIBLINGS_NODE)] = node->getParent() ? node->getParent()->getChildren().size() - 1 : 0;
					if(computeLeafRatio)
						buffer[attrNames->linearIndex(idx, LEAF_RATIO_NODE)] = 0.0f;
					if(computeBalance)
						buffer[attrNames->linearIndex(idx, BALANCE_NODE)] = 0.0f;
					if(computeAvgChildHeight)
						buffer[attrNames->linearIndex(idx, AVG_CHILD_HEIGHT_NODE)] = 0.0f;
				},
				[&](NodeMTPtr parent, NodeMTPtr child) {
					int pIdx = parent->getIndex();
					int cIdx = child->getIndex();

					// Acumulando descendentes
					bufferNumDesc[indexOfNumDescendants(pIdx)] += bufferNumDesc[indexOfNumDescendants(cIdx)] + 1;
					bufferNumLeafDesc[indexOfNumLeafDescendants(pIdx)] += bufferNumLeafDesc[indexOfNumLeafDescendants(cIdx)];
					
					// Altura
					float childHeight = bufferHeight[indexOfHeight(cIdx)];
					float& parentHeight = bufferHeight[indexOfHeight(pIdx)];
					parentHeight = std::max(parentHeight, childHeight + 1);
					int numChildren = parent->getChildren().size();

					// Balanceamento
					if(computeBalance){
						float& minH = buffer[attrNames->linearIndex(pIdx, BALANCE_NODE)]; // usado como mínimo temporário
						if (numChildren == 1) {
							minH = childHeight;
						} else {
							minH = std::min(minH, childHeight);;
						}
					}
						
					if(computeAvgChildHeight) {
						float& sumH = buffer[attrNames->linearIndex(pIdx, AVG_CHILD_HEIGHT_NODE)];
						if (numChildren == 1) {
							sumH = childHeight;
						} else {
							sumH += childHeight;
						}
					}

				},
				[&](NodeMTPtr node) {
					int idx = node->getIndex();
					
					if(computeLeafRatio){
						float desc = bufferNumDesc[indexOfNumDescendants(idx)];
						float folhas = bufferNumLeafDesc[indexOfNumLeafDescendants(idx)]; 
						// Razão folhas/descendentes
						buffer[attrNames->linearIndex(idx, LEAF_RATIO_NODE)] = desc > 0.0f ? folhas / (desc + 1.0f) : 1.0f;
					}

					// Balanceamento e média
					if (!node->getChildren().empty()) {
						if(computeBalance){
							float alturaMax = bufferHeight[indexOfHeight(idx)];
							float alturaMin = buffer[attrNames->linearIndex(idx, BALANCE_NODE)];
							buffer[attrNames->linearIndex(idx, BALANCE_NODE)] = alturaMax - alturaMin;
						}
						
						if(computeAvgChildHeight){
							buffer[attrNames->linearIndex(idx, AVG_CHILD_HEIGHT_NODE)] = buffer[attrNames->linearIndex(idx, AVG_CHILD_HEIGHT_NODE)] / node->getChildren().size();
						}
					}

				}
			);
			
		}
};






class AttributeFactory {
	private:
		static std::shared_ptr<AttributeComputer> createImpl(Attribute attr) {
			switch (attr) {
				case AREA: return std::make_shared<AreaComputer>();
				case VOLUME: return std::make_shared<VolumeComputer>();
				
				case DYNAMICS: return std::make_shared<DynamicsComputer>();
				
				case LEVEL: 
				case MEAN_LEVEL:
				case VARIANCE_LEVEL: return std::make_shared<GrayLevelStatsComputer>();

				case BOX_COL_MIN:
				case BOX_COL_MAX:
				case BOX_ROW_MIN:
				case BOX_ROW_MAX:
				case RATIO_WH: 
				case RECTANGULARITY: 
				case DIAGONAL_LENGTH:
				case BOX_HEIGHT:
				case BOX_WIDTH: 
					return std::make_shared<BoundingBoxComputer>();
				

				case AXIS_ORIENTATION: 
				case LENGTH_MAJOR_AXIS: 
				case LENGTH_MINOR_AXIS: 
				case ECCENTRICITY: 
				case INERTIA:
				case COMPACTNESS: 
					return std::make_shared<MomentBasedAttributeComputer>();


				case CENTRAL_MOMENT_20:
				case CENTRAL_MOMENT_02:
				case CENTRAL_MOMENT_11:
				case CENTRAL_MOMENT_30:
				case CENTRAL_MOMENT_03:
				case CENTRAL_MOMENT_21:
				case CENTRAL_MOMENT_12:
					return std::make_shared<CentralMomentsComputer>();

				
				case HU_MOMENT_1: 
				case HU_MOMENT_2:
				case HU_MOMENT_3:
				case HU_MOMENT_4:
				case HU_MOMENT_5:
				case HU_MOMENT_6:
				case HU_MOMENT_7:
					return std::make_shared<HuMomentsComputer>();


				case HEIGHT_NODE:
				case DEPTH_NODE:
				case IS_LEAF_NODE:
				case IS_ROOT_NODE:
				case NUM_CHILDREN_NODE:
				case NUM_SIBLINGS_NODE:
				case NUM_DESCENDANTS_NODE:
				case NUM_LEAF_DESCENDANTS_NODE:
				case LEAF_RATIO_NODE:
				case BALANCE_NODE:
				case AVG_CHILD_HEIGHT_NODE:
					return std::make_shared<TreeTopologyComputer>();


				default:
					throw std::runtime_error("Attribute not supported.");
			}
		}

		static std::shared_ptr<AttributeComputer> createImpl(AttributeGroup group) {
			switch (group) {
				case AttributeGroup::BOUNDING_BOX:
					return std::make_shared<BoundingBoxComputer>();
				case AttributeGroup::CENTRAL_MOMENTS:
					return std::make_shared<CentralMomentsComputer>();
				case AttributeGroup::HU_MOMENTS:
					return std::make_shared<HuMomentsComputer>();
				case AttributeGroup::MOMENT_BASED:
					return std::make_shared<MomentBasedAttributeComputer>();
				case AttributeGroup::TREE_TOPOLOGY:
					return std::make_shared<TreeTopologyComputer>();
				default:
					throw std::runtime_error("Attribute group not supported.");
			}
		}

	public:
		static std::shared_ptr<AttributeComputer> create(const AttributeOrGroup& attr) {
			return std::visit([](auto&& actualAttr) -> std::shared_ptr<AttributeComputer> {
				return AttributeFactory::createImpl(actualAttr); // Correção aqui!
			}, attr);
		}



};







#endif 






		