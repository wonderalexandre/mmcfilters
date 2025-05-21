

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
#include <optional>

#define PI 3.14159265358979323846



class ContoursMT; //forward declaration
using ContoursMTPtr = std::shared_ptr<ContoursMT>;


enum class Attribute {
    // Geométricos básicos
    AREA,

    // Textura / intensidade agregada
    VOLUME,
	RELATIVE_VOLUME,
    LEVEL,
    GRAY_HEIGHT,
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

struct AttributeKey {
    Attribute attr;
    int delta = 0;

    AttributeKey(Attribute a, int d = 0) : attr(a), delta(d) {}

    bool operator==(const AttributeKey& other) const {
        return attr == other.attr && delta == other.delta;
    }
};

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

	template<>
    struct hash<AttributeKey> {
        std::size_t operator()(const AttributeKey& k) const {
            return std::hash<int>()(static_cast<int>(k.attr)) ^ (std::hash<int>()(k.delta) << 1);
        }
    };

}


static const std::unordered_map<AttributeGroup, std::vector<Attribute>> ATTRIBUTE_GROUPS = {
    {AttributeGroup::GEOMETRIC, {
        AREA,
		VOLUME,
		RELATIVE_VOLUME,
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
		RELATIVE_VOLUME,
        LEVEL,
        GRAY_HEIGHT,
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

class AttributeNamesWithDelta;  // forward declaration
using AttributeNamesWithDeltaPtr = std::shared_ptr<AttributeNamesWithDelta>;
class AttributeNamesWithDelta {
public:
    std::unordered_map<AttributeKey, int> indexMap;
    const int NUM_ATTRIBUTES;

    AttributeNamesWithDelta(std::unordered_map<AttributeKey, int>&& map)
        : indexMap(std::move(map)), NUM_ATTRIBUTES(static_cast<int>(indexMap.size())) {}


	static AttributeNamesWithDelta create(int n, int delta, const std::vector<Attribute>& attributes) {
		std::unordered_map<AttributeKey, int> map;
		int offset = 0;
		for (int d = -delta; d <= delta; ++d) {
			for (int i = 0; i < attributes.size(); ++i) {
				map[AttributeKey{attributes[i], d}] = offset++;
			}
		}
		return AttributeNamesWithDelta(std::move(map));
	}

    int getIndex(Attribute attr, int delta) const {
        return getIndex(AttributeKey{attr, delta});
    }

	int getIndex(AttributeKey attrKey) const {
		return indexMap.at(attrKey);
	}

    int linearIndex(int nodeIndex, Attribute attr, int delta) const {
        return nodeIndex * NUM_ATTRIBUTES + getIndex(attr, delta);
    }

	int linearIndex(int nodeIndex, AttributeKey attrKey) const {
		return linearIndex(nodeIndex, attrKey.attr, attrKey.delta);
	}

	static std::string toString(AttributeKey attrKey) {
		return toString(attrKey.attr, attrKey.delta);
	}

	static std::string toString(Attribute attr, int delta) {
		std::string name;
		switch (attr) {
			case AREA: name = "AREA"; break;
			case VOLUME: name = "VOLUME"; break;
			case RELATIVE_VOLUME: name = "RELATIVE_VOLUME"; break;
			case LEVEL: name = "LEVEL"; break;
			case GRAY_HEIGHT: name = "GRAY_HEIGHT"; break;
			case MEAN_LEVEL: name = "MEAN_LEVEL"; break;
			case VARIANCE_LEVEL: name = "VARIANCE_LEVEL"; break;
			case BOX_WIDTH: name = "BOX_WIDTH"; break;
			case BOX_HEIGHT: name = "BOX_HEIGHT"; break;
			case RECTANGULARITY: name = "RECTANGULARITY"; break;
			case RATIO_WH: name = "RATIO_WH"; break;
			case DIAGONAL_LENGTH: name = "DIAGONAL_LENGTH"; break;
			case BOX_COL_MIN: name = "BOX_COL_MIN"; break;
			case BOX_COL_MAX: name = "BOX_COL_MAX"; break;
			case BOX_ROW_MIN: name = "BOX_ROW_MIN"; break;
			case BOX_ROW_MAX: name = "BOX_ROW_MAX"; break;
			case CENTRAL_MOMENT_20: name = "CENTRAL_MOMENT_20"; break;
			case CENTRAL_MOMENT_02: name = "CENTRAL_MOMENT_02"; break;
			case CENTRAL_MOMENT_11: name = "CENTRAL_MOMENT_11"; break;
			case CENTRAL_MOMENT_30: name = "CENTRAL_MOMENT_30"; break;
			case CENTRAL_MOMENT_03: name = "CENTRAL_MOMENT_03"; break;
			case CENTRAL_MOMENT_21: name = "CENTRAL_MOMENT_21"; break;
			case CENTRAL_MOMENT_12: name = "CENTRAL_MOMENT_12"; break;
			case HU_MOMENT_1: name = "HU_MOMENT_1"; break;
			case HU_MOMENT_2: name = "HU_MOMENT_2"; break;
			case HU_MOMENT_3: name = "HU_MOMENT_3"; break;
			case HU_MOMENT_4: name = "HU_MOMENT_4"; break;
			case HU_MOMENT_5: name = "HU_MOMENT_5"; break;
			case HU_MOMENT_6: name = "HU_MOMENT_6"; break;
			case HU_MOMENT_7: name = "HU_MOMENT_7"; break;
			case INERTIA: name = "INERTIA"; break;
			case COMPACTNESS: name = "COMPACTNESS"; break;
			case ECCENTRICITY: name = "ECCENTRICITY"; break;
			case LENGTH_MAJOR_AXIS: name = "LENGTH_MAJOR_AXIS"; break;
			case LENGTH_MINOR_AXIS: name = "LENGTH_MINOR_AXIS"; break;
			case AXIS_ORIENTATION: name = "AXIS_ORIENTATION"; break;
			case HEIGHT_NODE: name = "HEIGHT_NODE"; break;
			case DEPTH_NODE: name = "DEPTH_NODE"; break;
			case IS_LEAF_NODE: name = "IS_LEAF_NODE"; break;
			case IS_ROOT_NODE: name = "IS_ROOT_NODE"; break;
			case NUM_CHILDREN_NODE: name = "NUM_CHILDREN_NODE"; break;
			case NUM_SIBLINGS_NODE: name = "NUM_SIBLINGS_NODE"; break;
			case NUM_DESCENDANTS_NODE: name = "NUM_DESCENDANTS_NODE"; break;
			case NUM_LEAF_DESCENDANTS_NODE: name = "NUM_LEAF_DESCENDANTS_NODE"; break;
			case LEAF_RATIO_NODE: name = "LEAF_RATIO_NODE"; break;
			case BALANCE_NODE: name = "BALANCE_NODE"; break;
			case AVG_CHILD_HEIGHT_NODE: name = "AVG_CHILD_HEIGHT_NODE"; break;
			default: name = "UNKNOWN"; break;
		}

		if (delta < 0)
			name += "_ASC_" + std::to_string(-delta);
		else if (delta > 0)
			name += "_DESC_" + std::to_string(delta);
		// Para delta == 0: permanece só o nome do atributo

		return name;
	}

};

class AttributeNames;  // forward declaration
using AttributeNamesPtr = std::shared_ptr<AttributeNames>;
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
			case RELATIVE_VOLUME: return "RELATIVE_VOLUME";
            case LEVEL: return "LEVEL";
            case GRAY_HEIGHT: return "GRAY_HEIGHT";
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
            case Attribute::VOLUME: return "Volume: Sum of the gray-level intensities of all pixels in the connected component. Interpreted as the total mass under the component, or the integral of the image function over its support.";
			case Attribute::RELATIVE_VOLUME: return "Relative volume: Sum of differences between the node level and the gray-levels of pixels in the component. Measures the amount of intensity required to fill the component to its node level.";
            case Attribute::LEVEL: return "Level: Gray-level at which the connected component appears in the threshold decomposition hierarchy; corresponds to the altitude of the node in the component tree.";
            case Attribute::GRAY_HEIGHT: return "GRAY_HEIGHT: For a node in the max-tree, the difference between its level and the maximum level among its descendants; analogously, in the min-tree, the difference to the minimum level among its descendants.";
            case Attribute::MEAN_LEVEL: return "Mean level: Average gray-level intensity of the pixels in the connected component.";
            case Attribute::VARIANCE_LEVEL: return "Variance of level: Variance of the gray-level values of the pixels in the connected component.";

            // Bounding box attributes
            case Attribute::BOX_WIDTH: return "Bounding box width: Width of the minimum rectangle enclosing the connected component.";
            case Attribute::BOX_HEIGHT: return "Bounding box height: Height of the minimum rectangle enclosing the connected component.";
            case Attribute::RECTANGULARITY: return "Rectangularity: Ratio between the connected component area and the area of its bounding box. Values closer to 1 indicate shapes that efficiently fill their bounding box.";
            case Attribute::RATIO_WH: return "Aspect ratio: Ratio of the bounding box width to its height. Describes the elongation of the component.";
            case Attribute::BOX_COL_MIN: return "Bounding box column min: Minimum column index covered by the connected component.";
            case Attribute::BOX_COL_MAX: return "Bounding box column max: Maximum column index covered by the connected component.";
            case Attribute::BOX_ROW_MIN: return "Bounding box row min: Minimum row index covered by the connected component.";
            case Attribute::BOX_ROW_MAX: return "Bounding box row max: Maximum row index covered by the connected component.";
			case Attribute::DIAGONAL_LENGTH: return "Diagonal length: Euclidean length of the diagonal of the bounding box, computed as sqrt(width^2 + height^2).";

            // Central moments
			case Attribute::CENTRAL_MOMENT_20: return "Central moment (2,0): Second-order moment about the centroid along the x-axis. Measures the horizontal spread of the component.";
			case Attribute::CENTRAL_MOMENT_02: return "Central moment (0,2): Second-order moment about the centroid along the y-axis. Measures the vertical spread of the component.";
			case Attribute::CENTRAL_MOMENT_11: return "Central moment (1,1): Mixed second-order moment about the centroid. Represents the covariance between x and y coordinates.";
			case Attribute::CENTRAL_MOMENT_30: return "Central moment (3,0): Third-order moment about the centroid along the x-axis. Describes horizontal asymmetry of the component.";
			case Attribute::CENTRAL_MOMENT_03: return "Central moment (0,3): Third-order moment about the centroid along the y-axis. Describes vertical asymmetry of the component.";
			case Attribute::CENTRAL_MOMENT_21: return "Central moment (2,1): Mixed third-order moment about the centroid. Captures joint spread and asymmetry in x and y.";
			case Attribute::CENTRAL_MOMENT_12: return "Central moment (1,2): Mixed third-order moment about the centroid. Captures joint spread and asymmetry in y and x.";

            // Hu moments (invariant shape descriptors)
            case Attribute::HU_MOMENT_1: return "Hu moment 1: Invariant to translation, scale, and rotation. Represents overall spatial variance (shape dispersion).";
			case Attribute::HU_MOMENT_2: return "Hu moment 2: Invariant capturing the difference between horizontal and vertical spread.";
			case Attribute::HU_MOMENT_3: return "Hu moment 3: Sensitive to skewness and asymmetry in the pixel distribution.";
			case Attribute::HU_MOMENT_4: return "Hu moment 4: Measures symmetry with respect to diagonal axes.";
			case Attribute::HU_MOMENT_5: return "Hu moment 5: Descriptor sensitive to orientation and reflection; captures complex asymmetries.";
			case Attribute::HU_MOMENT_6: return "Hu moment 6: Invariant capturing elliptic asymmetries, sensitive to specific shape curvature.";
			case Attribute::HU_MOMENT_7: return "Hu moment 7: Highly sensitive to irregularities and fine variations; helps discriminate mirror-symmetric shapes.";

            // Derived from moments
            case Attribute::INERTIA: return "Inertia: Sum of normalized second-order central moments (mu20 + mu02). Measures the dispersion of mass around the centroid. Higher values indicate objects with thin and elongated structures.";
            case Attribute::COMPACTNESS: return "Compactness: Area normalized by the shape's dispersion (mu20 + mu02). Higher values indicate more compact and isotropic shapes.";
            case Attribute::ECCENTRICITY: return "Eccentricity: Ratio of principal inertia axes (λ_1/λ_2). Quantifies the elongation of the shape; values near 1 indicate circularity.";
            case Attribute::LENGTH_MAJOR_AXIS: return "Major axis length: Length of the longest diameter of the shape.";
            case Attribute::LENGTH_MINOR_AXIS: return "Minor axis length: Length of the shortest diameter of the shape.";
            case Attribute::AXIS_ORIENTATION: return "Axis orientation: Angle of the principal inertia axis, computed from central moments. Indicates the dominant orientation of the shape.";

            // tree topology 
			case Attribute::HEIGHT_NODE: return "Height: Longest path from this node to any leaf in its subtree. Measures the depth of the subtree rooted at the node.";
			case Attribute::DEPTH_NODE: return "Depth: Number of steps from this node to the root of the tree. Indicates the level of embedding within the tree hierarchy.";
			case Attribute::IS_LEAF_NODE: return "Is leaf: True if the node has no children, i.e., it represents a minimal component in the hierarchy.";
			case Attribute::IS_ROOT_NODE:return "Is root: True if the node is the root of the tree, representing the entire image support.";
			case Attribute::NUM_CHILDREN_NODE:return "Number of children: Count of direct child nodes. Reflects the immediate branching factor of the node.";
			case Attribute::NUM_SIBLINGS_NODE:return "Number of siblings: Number of other nodes that share the same parent.";
			case Attribute::NUM_DESCENDANTS_NODE:return "Number of descendants: Total number of nodes in the subtree rooted at this node (excluding the node itself).";
			case Attribute::NUM_LEAF_DESCENDANTS_NODE:return "Number of leaf descendants: Number of leaf nodes in the subtree. Reflects the number of minimal patterns under this structure.";
			case Attribute::LEAF_RATIO_NODE:return "Leaf ratio: Ratio of leaf descendants to total descendants. Measures structural 'flatness' or terminal density of the subtree.";
			case Attribute::BALANCE_NODE:return "Balance: Difference between the maximum and minimum heights among the subtrees of the children. Indicates branching symmetry.";
			case Attribute::AVG_CHILD_HEIGHT_NODE:return "Average child height: Mean height of all direct child subtrees. Useful for measuring uniformity of the subtree structure.";

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



class AttributeComputer;  // forward declaration
using AttributeComputerPtr = std::shared_ptr<AttributeComputer>;
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


class PathAscendantsAndDescendants{
	private:
	MorphologicalTreePtr tree;
	std::vector<NodeMTPtr> ascendants;
	std::vector<NodeMTPtr> descendants;

	public:
	PathAscendantsAndDescendants(MorphologicalTreePtr tree): tree(tree){ }

    NodeMTPtr getNodeAscendantBasedOnHierarchical(NodeMTPtr node, int h){
		NodeMTPtr n = node;
		int i=0;
		while(i++ < h){
			n = n->getParent();
			if(n == nullptr)
				return n;
		}
		return n;
	}

	void maxAreaDescendants(NodeMTPtr nodeAsc, NodeMTPtr nodeDes){
		if(descendants[nodeAsc->getIndex()] == nullptr)
			descendants[nodeAsc->getIndex()] = nodeDes;
		
		if(descendants[nodeAsc->getIndex()]->getAreaCC() < nodeDes->getAreaCC())
			descendants[nodeAsc->getIndex()] = nodeDes;
		
	}

	void computerAscendantsAndDescendants(int delta){
		std::vector<NodeMTPtr> tmp_asc (this->tree->getNumNodes(), nullptr);
		this->ascendants = tmp_asc;

		std::vector<NodeMTPtr> tmp_des (this->tree->getNumNodes(), nullptr);
		this->descendants = tmp_des;
		
		for(NodeMTPtr node: tree->getIndexNode()){
			NodeMTPtr nodeAsc = this->getNodeAscendantBasedOnHierarchical(node, delta);
			if(nodeAsc == nullptr) continue;
			this->maxAreaDescendants(nodeAsc, node);
			if(descendants[nodeAsc->getIndex()] != nullptr){
				ascendants[node->getIndex()] = nodeAsc;
			}
		}
	}
	std::vector<NodeMTPtr>& getAscendants() {
		return this->ascendants;
	}

	std::vector<NodeMTPtr>& getDescendants() {
		return this->descendants;
	}
};


struct ExtinctionValues{
	NodeMTPtr leaf;
	NodeMTPtr cutoffNode;
	float extinction;
	ExtinctionValues(NodeMTPtr leaf, NodeMTPtr cutoffNode, float extinction)
		: leaf(leaf), cutoffNode(cutoffNode), extinction(extinction) {}
	
};
using ExtinctionValuesPtr = std::shared_ptr<ExtinctionValues>;

class AttributeComputedIncrementally;  // forward declaration
using AttributeComputedIncrementallyPtr = std::shared_ptr<AttributeComputedIncrementally>;
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



	static std::vector<ExtinctionValuesPtr> getExtinctionValue(MorphologicalTreePtr tree, std::shared_ptr<float[]> attr){
		std::list<NodeMTPtr> leaves = tree->getLeaves();
		std::vector<ExtinctionValuesPtr> leavesByExtinction;
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
			leavesByExtinction.push_back( std::make_shared<ExtinctionValues>(leaf, cutoffNode, extinction) );
			
		}
		return leavesByExtinction;
	}

	static ContoursMTPtr extractCompactCountors(MorphologicalTreePtr tree);

	static std::vector<std::unordered_set<int>> extractNonCompactCountors(MorphologicalTreePtr tree);

	static std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>> computeAttributesByComputer(MorphologicalTreePtr tree, std::shared_ptr<AttributeComputer> comp, const DependencyMap& available = {});
	
	static std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>> computeSingleAttribute(MorphologicalTreePtr tree, AttributeOrGroup attr, const DependencyMap& availableDeps = {});
	
	static std::pair<std::shared_ptr<AttributeNamesWithDelta>, std::shared_ptr<float[]>> computeSingleAttributeWithDelta(MorphologicalTreePtr tree, Attribute attribute, int delta, std::string padding="last-padding", const DependencyMap& availableDeps={});

	static std::pair<std::shared_ptr<AttributeNames>, std::shared_ptr<float[]>> computeAttributes(MorphologicalTreePtr tree, const std::vector<AttributeOrGroup>& attributes,const DependencyMap& providedDependencies={});
};



class AreaComputer : public AttributeComputer {
	public:
		
		std::vector<Attribute> attributes() const override {
			return {AREA};
		}
		void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<Attribute>& requestedAttributes, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources= {}) const override {
			if(PRINT_LOG) std::cout << "\n==== AttributeComputer: Computing AREA" << std::endl;
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
using AreaComputerPtr = std::shared_ptr<AreaComputer>;

class VolumeComputer : public AttributeComputer {
	public:
		std::vector<Attribute> attributes() const override {
			return {VOLUME, RELATIVE_VOLUME};
		}
		void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<Attribute>& requestedAttributes, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources= {}) const override {
			if(PRINT_LOG) std::cout << "\n==== AttributeComputer: Computing VOLUME" << std::endl;
			auto indexOfVol = [&](int idx) {
				return attrNames->linearIndex(idx, VOLUME);
			};
			auto indexOfVolRel = [&](int idx) {
				return attrNames->linearIndex(idx, RELATIVE_VOLUME);
			};
	
			bool computeVolume = std::find(requestedAttributes.begin(), requestedAttributes.end(), VOLUME) != requestedAttributes.end();
			bool computeRelativeVolume = std::find(requestedAttributes.begin(), requestedAttributes.end(), RELATIVE_VOLUME) != requestedAttributes.end();

			AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
				[&](NodeMTPtr node) {
					if(computeVolume)
						buffer[indexOfVol(node->getIndex())] = node->getCNPs().size() * node->getLevel();
					if(computeRelativeVolume)	
						buffer[indexOfVolRel(node->getIndex())] = 0;
				},
				[&](NodeMTPtr parent, NodeMTPtr child) {
					if(computeVolume)
						buffer[indexOfVol(parent->getIndex())] += buffer[indexOfVol(child->getIndex())];
					if(computeRelativeVolume)		
						buffer[indexOfVolRel(parent->getIndex())] += buffer[indexOfVolRel(child->getIndex())] + child->getAreaCC() * std::abs(child->getLevel() - parent->getLevel());
				},
				[&](NodeMTPtr node) {
					if(computeRelativeVolume)	
						buffer[indexOfVolRel(node->getIndex())] += node->getAreaCC();
				}
			);
		}
};
using VolumeComputerPtr = std::shared_ptr<VolumeComputer>;


class GrayLevelStatsComputer : public AttributeComputer {
	public:
		
	
		std::vector<Attribute> attributes() const override {
			return { LEVEL, MEAN_LEVEL, VARIANCE_LEVEL, GRAY_HEIGHT };
		}
	
		std::vector<AttributeOrGroup> requiredAttributes() const override {
			return { VOLUME };
		}
	
		void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<Attribute>& requestedAttributes, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources = {}) const override {
			
			if(PRINT_LOG) std::cout << "\n==== AttributeComputer: Computing GrayLevelStatsComputer " << std::endl;
	
			auto indexOfMean = [&](int idx) {
				return attrNames->linearIndex(idx, MEAN_LEVEL);
			};
			auto indexOfLevel = [&](int idx) {
				return attrNames->linearIndex(idx, LEVEL);
			};
			auto indexOfVarianceLevel = [&](int idx) {
				return attrNames->linearIndex(idx, VARIANCE_LEVEL);
			};

			auto indexOfGrayHeight = [&](int idx) {
				return attrNames->linearIndex(idx, GRAY_HEIGHT);
			};

			bool computeMeanLevel = std::find(requestedAttributes.begin(), requestedAttributes.end(), MEAN_LEVEL) != requestedAttributes.end();
			bool computeVarianceLevel = std::find(requestedAttributes.begin(), requestedAttributes.end(), VARIANCE_LEVEL) != requestedAttributes.end();
			bool computeLevel = std::find(requestedAttributes.begin(), requestedAttributes.end(), LEVEL) != requestedAttributes.end();
			bool computeGrayHeight = std::find(requestedAttributes.begin(), requestedAttributes.end(), GRAY_HEIGHT) != requestedAttributes.end();
			

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
					if(computeGrayHeight)
						buffer[ indexOfGrayHeight(node->getIndex()) ] = node->getLevel();
				},
				[&](NodeMTPtr parent, NodeMTPtr child) { 
					if(computeVarianceLevel)
						sumGrayLevelSquare[parent->getIndex()] += sumGrayLevelSquare[child->getIndex()];
					if(computeGrayHeight){
						if (tree->getTreeType()==tree->MAX_TREE || parent->isMaxtreeNode()) {
							buffer[ indexOfGrayHeight(parent->getIndex()) ] = std::max(buffer[indexOfGrayHeight(parent->getIndex())], buffer[indexOfGrayHeight(child->getIndex())] );
						} else {
							buffer[ indexOfGrayHeight(parent->getIndex()) ] = std::min(buffer[indexOfGrayHeight(parent->getIndex())], buffer[indexOfGrayHeight(child->getIndex())] );
						}
					}
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
			if(computeGrayHeight) {
				for(NodeMTPtr node: tree->getRoot()->getIteratorPostOrderTraversal()){
					if(node->isLeaf())
						buffer[indexOfGrayHeight(node->getIndex())] = 0;	
					else
						buffer[ indexOfGrayHeight(node->getIndex()) ] = std::abs( node->getLevel() - buffer[ indexOfGrayHeight(node->getIndex())] ) + 1;
				}
			}
		}
};
using GrayLevelStatsComputerPtr = std::shared_ptr<GrayLevelStatsComputer>;


class BoundingBoxComputer : public AttributeComputer {
	public:
		
	
		std::vector<Attribute> attributes() const override {
			return { BOX_WIDTH, BOX_HEIGHT, RECTANGULARITY, RATIO_WH, BOX_COL_MIN, BOX_COL_MAX, BOX_ROW_MIN, BOX_ROW_MAX, DIAGONAL_LENGTH };
		}

		
		void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<Attribute>& requestedAttributes, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources = {}) const override {
			
			if(PRINT_LOG) std::cout << "\n==== AttributeComputer: Computing BOUNDING_BOX group" << std::endl;
	
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
using BoundingBoxComputerPtr = std::shared_ptr<BoundingBoxComputer>;

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
			if(PRINT_LOG) std::cout << "\n==== AttributeComputer: Computing CENTRAL_MOMENT group" << std::endl;
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
using CentralMomentsComputerPtr = std::shared_ptr<CentralMomentsComputer>;

class MomentBasedAttributeComputer : public AttributeComputer {
	public:
		
	
		std::vector<Attribute> attributes() const override {
			return { LENGTH_MAJOR_AXIS, LENGTH_MINOR_AXIS, ECCENTRICITY, COMPACTNESS, AXIS_ORIENTATION, INERTIA};
		}
	
		std::vector<AttributeOrGroup> requiredAttributes() const override {
			return { AttributeGroup::CENTRAL_MOMENTS };
		}
	
		void compute(MorphologicalTreePtr tree, std::shared_ptr<float[]> buffer, std::shared_ptr<AttributeNames> attrNames, const std::vector<Attribute>& requestedAttributes, const std::vector<std::pair<std::shared_ptr<AttributeNames>, const std::shared_ptr<float[]>>>& dependencySources = {}) const override {
			
			if(PRINT_LOG) std::cout << "\n==== AttributeComputer: Computing MOMENT_BASED group" << std::endl;
	
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
using MomentBasedAttributeComputerPtr = std::shared_ptr<MomentBasedAttributeComputer>;

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
			if(PRINT_LOG) std::cout << "\n==== AttributeComputer: Computing HU_MOMENT group" << std::endl;
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
using HuMomentsComputerPtr = std::shared_ptr<HuMomentsComputer>;

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
			
			if(PRINT_LOG) std::cout << "\n==== AttributeComputer: Computing STRUCTURE_TREE group" << std::endl;

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
using TreeTopologyComputerPtr = std::shared_ptr<TreeTopologyComputer>;





class AttributeFactory {
	private:
		static std::shared_ptr<AttributeComputer> createImpl(Attribute attr) {
			switch (attr) {
				case AREA: return std::make_shared<AreaComputer>();
				
				case RELATIVE_VOLUME:
				case VOLUME: return std::make_shared<VolumeComputer>();
				
				case GRAY_HEIGHT: 
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






/*
Computacao incremental de countours
*/
class ContoursMT{
	private:
	std::vector<std::list<int>> contours;
	std::vector<std::list<int>> contoursToRemove;

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


	class ContourPostOrderIterator {
		private:
		using value_type = std::pair<NodeMTPtr, std::unordered_set<int>>;
		using reference = value_type&;
		using pointer = value_type*;
		using iterator_category = std::input_iterator_tag;

		ContoursMT* contoursMT;
		std::stack<NodeMTPtr> outputStack;
		std::vector<std::unique_ptr<std::unordered_set<int>>> contoursByNodes;
		value_type currentValue;

		void advance() {
			if (!outputStack.empty()) {
				NodeMTPtr node = outputStack.top(); outputStack.pop();

				// Merge dos filhos (igual antes)
				for (NodeMTPtr child : node->getChildren()) {
					auto& parentContour = contoursByNodes[node->getIndex()];
					auto& childContour = contoursByNodes[child->getIndex()];
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
				auto& contour = contoursByNodes[node->getIndex()];
				if (!contour) contour = std::make_unique<std::unordered_set<int>>();
				for (int p : contoursMT->contours[node->getIndex()])
					contour->insert(p);
				for (int p : contoursMT->contoursToRemove[node->getIndex()])
					contour->erase(p);

				currentValue = std::make_pair(node, *contour);
			}
		}

		public:
		ContourPostOrderIterator(NodeMTPtr root, ContoursMT* contoursMT) : contoursMT(contoursMT){
			// Travessia prévia para montar outputStack em pós-ordem
			if(root){
				std::stack<NodeMTPtr> tempStack;
				tempStack.push(root);
				while (!tempStack.empty()) {
					NodeMTPtr current = tempStack.top(); tempStack.pop();
					outputStack.push(current);
					for (NodeMTPtr child : current->getChildren()) {
						tempStack.push(child);
					}
				}
				int numNodes = root->getNumDescendants() + 1;
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

	class ContourPostOrderRange {
		private:
		NodeMTPtr root;
		ContoursMT* contoursMT;

		public:
		ContourPostOrderRange(NodeMTPtr root, ContoursMT* contoursMT) : root(root), contoursMT(contoursMT) {}

		ContourPostOrderIterator begin() { return ContourPostOrderIterator(root, contoursMT); }
		ContourPostOrderIterator end() { return ContourPostOrderIterator(nullptr, contoursMT); }
	};

	ContourPostOrderRange contoursLazy(NodeMTPtr root) {
		return ContourPostOrderRange(root, this);
	}

};







#endif 






		