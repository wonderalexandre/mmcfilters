#pragma once

#include "../utils/Common.hpp"

namespace mmcfilters {

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
	CIRCULARITY,

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
    AVG_CHILD_HEIGHT_NODE,

	//BitQuads
	BITQUADS_AREA,
	BITQUADS_NUMBER_EULER,
	BITQUADS_NUMBER_HOLES,
	BITQUADS_PERIMETER,
	BITQUADS_PERIMETER_CONTINUOUS,
	BITQUADS_CIRCULARITY,
	BITQUADS_PERIMETER_AVERAGE,
	BITQUADS_LENGTH_AVERAGE,
	BITQUADS_WIDTH_AVERAGE

};


enum class AttributeGroup {
    ALL,               // Todos os atributos
    GEOMETRIC,         // Forma e proporção
	MOMENT_BASED,
    BOUNDING_BOX,      // Box width/height
    CENTRAL_MOMENTS,    // Momentos centrais
    HU_MOMENTS,         // Momentos de Hu
    TEXTURE,           // Atributos baseados em níveis de cinza
    TREE_TOPOLOGY,         // Topologia da árvore
	BITQUADS          // BitQuads
	
};

using AttributeOrGroup = std::variant<Attribute, AttributeGroup>;
using enum Attribute;

/**
 * @brief Chave composta usada para indexar atributos associados a deslocamentos.
 */
struct AttributeKey {
    Attribute attr;
    int delta = 0;

    AttributeKey(Attribute a, int d = 0) : attr(a), delta(d) {}

    bool operator==(const AttributeKey& other) const {
        return attr == other.attr && delta == other.delta;
    }
};

} // namespace mmcfilters

namespace std {
    template <>
    struct hash<mmcfilters::AttributeGroup> {
        std::size_t operator()(const mmcfilters::AttributeGroup& group) const noexcept {
            return static_cast<std::size_t>(group);
        }
    };

    template <>
    struct hash<mmcfilters::Attribute> {
        std::size_t operator()(const mmcfilters::Attribute& attr) const noexcept {
            return static_cast<std::size_t>(attr);
        }
    };

    template <>
    struct hash<mmcfilters::AttributeOrGroup> {
        std::size_t operator()(const mmcfilters::AttributeOrGroup& attr) const {
            return std::visit([](auto&& a) -> std::size_t {
                return std::hash<std::decay_t<decltype(a)>>{}(a);
            }, attr);
        }
    };

    template <>
    struct hash<mmcfilters::AttributeKey> {
        std::size_t operator()(const mmcfilters::AttributeKey& k) const {
            return std::hash<int>()(static_cast<int>(k.attr)) ^ (std::hash<int>()(k.delta) << 1);
        }
    };
}

namespace mmcfilters {
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
		INERTIA,
		CIRCULARITY
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
	 {AttributeGroup::BITQUADS, {
        BITQUADS_AREA,
		BITQUADS_NUMBER_EULER,
		BITQUADS_NUMBER_HOLES,
		BITQUADS_PERIMETER,
		BITQUADS_PERIMETER_CONTINUOUS,
		BITQUADS_CIRCULARITY,
		BITQUADS_PERIMETER_AVERAGE,
		BITQUADS_LENGTH_AVERAGE,
		BITQUADS_WIDTH_AVERAGE
    }},
    {AttributeGroup::ALL, [] {
        std::vector<Attribute> all;
        for (int i = 0; i <= static_cast<int>(AVG_CHILD_HEIGHT_NODE); ++i)
            all.push_back(static_cast<Attribute>(i));
        return all;
    }()}
};


/**
 * @brief Mapeia atributos e variações de delta para índices lineares em buffers.
 */
class AttributeNamesWithDelta {
public:
    std::unordered_map<AttributeKey, int> indexMap;
    const int NUM_ATTRIBUTES;

    AttributeNamesWithDelta(std::unordered_map<AttributeKey, int>&& map)
        : indexMap(std::move(map)), NUM_ATTRIBUTES(static_cast<int>(indexMap.size())) {}


	static AttributeNamesWithDelta create(int delta, const std::vector<Attribute>& attributes) {
		std::unordered_map<AttributeKey, int> map;
		int offset = 0;
		for (int d = -delta; d <= delta; ++d) {
			for (std::size_t i = 0; i < attributes.size(); ++i) {
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
			case CIRCULARITY: name = "CIRCULARITY"; break;
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
			case BITQUADS_AREA: name = "BITQUADS_AREA"; break;
			case BITQUADS_NUMBER_EULER: name = "BITQUADS_NUMBER_EULER"; break;
			case BITQUADS_NUMBER_HOLES: name = "BITQUADS_NUMBER_HOLES"; break;
			case BITQUADS_PERIMETER: name = "BITQUADS_PERIMETER"; break;
			case BITQUADS_PERIMETER_CONTINUOUS: name = "BITQUADS_PERIMETER_CONTINUOUS"; break;
			case BITQUADS_CIRCULARITY: name = "BITQUADS_CIRCULARITY"; break;
			case BITQUADS_PERIMETER_AVERAGE: name = "BITQUADS_PERIMETER_AVERAGE"; break;
			case BITQUADS_LENGTH_AVERAGE: name = "BITQUADS_LENGTH_AVERAGE"; break;
			case BITQUADS_WIDTH_AVERAGE: name = "BITQUADS_WIDTH_AVERAGE"; break;
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


/**
 * @brief Resolve nomes de atributos em índices lineares usados pelos buffers.
 */
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
			case CIRCULARITY: return "CIRCULARITY";
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
			case BITQUADS_AREA: return "BITQUADS_AREA";
			case BITQUADS_NUMBER_EULER: return "BITQUADS_NUMBER_EULER";
			case BITQUADS_NUMBER_HOLES: return "BITQUADS_NUMBER_HOLES";
			case BITQUADS_PERIMETER: return "BITQUADS_PERIMETER";
			case BITQUADS_PERIMETER_CONTINUOUS: return "BITQUADS_PERIMETER_CONTINUOUS";
			case BITQUADS_CIRCULARITY: return "BITQUADS_CIRCULARITY";
			case BITQUADS_PERIMETER_AVERAGE: return "BITQUADS_PERIMETER_AVERAGE";
			case BITQUADS_LENGTH_AVERAGE: return "BITQUADS_LENGTH_AVERAGE";
			case BITQUADS_WIDTH_AVERAGE: return "BITQUADS_WIDTH_AVERAGE";
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
			case Attribute::ECCENTRICITY: return "Eccentricity: Ratio of principal inertia eigenvalues (λ_1/λ_2). Measures elongation; values near 1 indicate circularity, higher values indicate elongation.";
			case Attribute::CIRCULARITY: return "Circularity: Ratio of the minor to major eigenvalues of the inertia matrix (λ_2/λ_1), i.e., Inverse of eccentricity. Indicates how circular a shape is; values near 1 suggest circularity, values near 0 indicate elongation.";
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
			
			// Bitquads-based shape attributes
			case Attribute::BITQUADS_AREA: return "BitQuads area (Duda): Refined sub-pixel area estimation using fractional weights based on the geometric contribution of local 2x2 pixel patterns.";
			case Attribute::BITQUADS_NUMBER_EULER: return "BitQuads Euler number: Topological invariant computed as the number of connected components minus the number of holes, using 2x2 pattern statistics under 4- or 8-connectivity.";
			case Attribute::BITQUADS_NUMBER_HOLES: return "BitQuads number of holes: Number of interior holes in the component, derived from the Euler characteristic assuming a single connected object.";
			case Attribute::BITQUADS_PERIMETER: return "BitQuads perimeter: Discrete approximation of the shape's boundary length, calculated by summing edge-contributing patterns in the 2x2 pixel grid.";
			case Attribute::BITQUADS_PERIMETER_CONTINUOUS: return "BitQuads continuous perimeter: Smoothed estimation of the boundary length, incorporating weighted transitions across pixel edges and diagonals.";
			case Attribute::BITQUADS_CIRCULARITY: return "BitQuads circularity: Compactness measure defined as (4π x areaDuda) / perimeter². Values close to 1 indicate circular shapes; lower values suggest elongation or irregularity.";
			case Attribute::BITQUADS_PERIMETER_AVERAGE: return "BitQuads average perimeter: Mean perimeter per connected component, accounting for complex structures and holes.";
			case Attribute::BITQUADS_LENGTH_AVERAGE: return "BitQuads average length: Estimated average longitudinal extent per component, derived from the average perimeter.";
			case Attribute::BITQUADS_WIDTH_AVERAGE: return "BitQuads average width: Estimated transverse extent per component, computed as (2 x average area) / average perimeter.";

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

} // namespace mmcfilters

