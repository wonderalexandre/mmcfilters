#pragma once

#include "AttributeTypes.hpp"

#include <array>
#include <cstddef>
#include <optional>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace mmcfilters::attributes::registry {

/**
 * @brief Stable metadata attached to one scalar attribute descriptor.
 *
 * The registry is the single source used by public names, user-facing
 * descriptions, group expansion, and pipeline routing decisions. The enum value
 * stored in `attribute` must match the position occupied in
 * `ATTRIBUTE_METADATA`, because lookup is intentionally O(1) by enum ordinal.
 */
struct AttributeMetadata {
    /// Scalar attribute described by this row.
    Attribute attribute;

    /// Stable symbolic name used in Python dictionaries, string parsing, and docs.
    std::string_view name;

    /// User-facing description suitable for notebooks and `Attribute.describe`.
    std::string_view description;

    /// Whether the descriptor reads node altitude values during computation.
    bool requiresAltitude = false;

    /// Whether the descriptor can be computed from topology/support alone.
    bool topologyOnly = false;
};

/**
 * @brief Canonical metadata table indexed by `Attribute` ordinal.
 *
 * Adding an attribute requires adding one row here and keeping the enum order in
 * `AttributeTypes.hpp` consistent with the table index. `metadata(...)` rejects
 * rows whose stored enum does not match the queried ordinal, which makes enum and
 * table drift fail as an unknown attribute instead of returning incorrect
 * metadata.
 */
inline constexpr std::array<AttributeMetadata, static_cast<std::size_t>(Attribute::CONTOUR_SIDE_SOUTH) + 1> ATTRIBUTE_METADATA{{
    {AREA, "AREA", "Area: Number of pixels in the connected component.", false, false},
    {VOLUME, "VOLUME", "Volume: Sum of the gray-level intensities of all pixels in the connected component. Interpreted as the total mass under the component, or the integral of the image function over its support.", true, false},
    {RELATIVE_VOLUME, "RELATIVE_VOLUME", "Relative volume: Sum of differences between the node level and the gray-levels of pixels in the component. Measures the amount of intensity required to fill the component to its node level.", true, false},
    {LEVEL, "LEVEL", "Level: Gray-level at which the connected component appears in the threshold decomposition hierarchy; corresponds to the altitude of the node in the component tree.", true, false},
    {GRAY_HEIGHT, "GRAY_HEIGHT", "GRAY_HEIGHT: For a node in the max-tree, the absolute difference between its level and the maximum level among its descendants; analogously, in the min-tree, the difference to the minimum level among its descendants. Leaves have gray height 0.", true, false},
    {MEAN_LEVEL, "MEAN_LEVEL", "Mean level: Average gray-level intensity of the pixels in the connected component.", true, false},
    {VARIANCE_LEVEL, "VARIANCE_LEVEL", "Variance of level: Variance of the gray-level values of the pixels in the connected component.", true, false},

    {BOX_WIDTH, "BOX_WIDTH", "Bounding box width: Width of the minimum rectangle enclosing the connected component.", false, true},
    {BOX_HEIGHT, "BOX_HEIGHT", "Bounding box height: Height of the minimum rectangle enclosing the connected component.", false, true},
    {DIAGONAL_LENGTH, "DIAGONAL_LENGTH", "Diagonal length: Euclidean length of the diagonal of the bounding box, computed as sqrt(width^2 + height^2).", false, true},
    {RECTANGULARITY, "RECTANGULARITY", "Rectangularity: Ratio between the connected component area and the area of its bounding box. Values closer to 1 indicate shapes that efficiently fill their bounding box.", false, true},
    {RATIO_WH, "RATIO_WH", "Aspect ratio: Ratio of the bounding box width to its height. Describes the elongation of the component.", false, true},
    {BOX_COL_MIN, "BOX_COL_MIN", "Bounding box column min: Minimum column index covered by the connected component.", false, true},
    {BOX_COL_MAX, "BOX_COL_MAX", "Bounding box column max: Maximum column index covered by the connected component.", false, true},
    {BOX_ROW_MIN, "BOX_ROW_MIN", "Bounding box row min: Minimum row index covered by the connected component.", false, true},
    {BOX_ROW_MAX, "BOX_ROW_MAX", "Bounding box row max: Maximum row index covered by the connected component.", false, true},

    {CENTRAL_MOMENT_20, "CENTRAL_MOMENT_20", "Central moment (2,0): Second-order moment about the centroid along the x-axis. Measures the horizontal spread of the component.", false, true},
    {CENTRAL_MOMENT_02, "CENTRAL_MOMENT_02", "Central moment (0,2): Second-order moment about the centroid along the y-axis. Measures the vertical spread of the component.", false, true},
    {CENTRAL_MOMENT_11, "CENTRAL_MOMENT_11", "Central moment (1,1): Mixed second-order moment about the centroid. Represents the covariance between x and y coordinates.", false, true},
    {CENTRAL_MOMENT_30, "CENTRAL_MOMENT_30", "Central moment (3,0): Third-order moment about the centroid along the x-axis. Describes horizontal asymmetry of the component.", false, true},
    {CENTRAL_MOMENT_03, "CENTRAL_MOMENT_03", "Central moment (0,3): Third-order moment about the centroid along the y-axis. Describes vertical asymmetry of the component.", false, true},
    {CENTRAL_MOMENT_21, "CENTRAL_MOMENT_21", "Central moment (2,1): Mixed third-order moment about the centroid. Captures joint spread and asymmetry in x and y.", false, true},
    {CENTRAL_MOMENT_12, "CENTRAL_MOMENT_12", "Central moment (1,2): Mixed third-order moment about the centroid. Captures joint spread and asymmetry in y and x.", false, true},

    {HU_MOMENT_1, "HU_MOMENT_1", "Hu moment 1: Invariant to translation, scale, and rotation. Represents overall spatial variance (shape dispersion).", false, true},
    {HU_MOMENT_2, "HU_MOMENT_2", "Hu moment 2: Invariant capturing the difference between horizontal and vertical spread.", false, true},
    {HU_MOMENT_3, "HU_MOMENT_3", "Hu moment 3: Sensitive to skewness and asymmetry in the pixel distribution.", false, true},
    {HU_MOMENT_4, "HU_MOMENT_4", "Hu moment 4: Measures symmetry with respect to diagonal axes.", false, true},
    {HU_MOMENT_5, "HU_MOMENT_5", "Hu moment 5: Descriptor sensitive to orientation and reflection; captures complex asymmetries.", false, true},
    {HU_MOMENT_6, "HU_MOMENT_6", "Hu moment 6: Invariant capturing elliptic asymmetries, sensitive to specific shape curvature.", false, true},
    {HU_MOMENT_7, "HU_MOMENT_7", "Hu moment 7: Highly sensitive to irregularities and fine variations; helps discriminate mirror-symmetric shapes.", false, true},

    {INERTIA, "INERTIA", "Inertia: Sum of normalized second-order central moments (mu20 + mu02). Measures the dispersion of mass around the centroid. Higher values indicate objects with thin and elongated structures.", false, true},
    {COMPACTNESS, "COMPACTNESS", "Compactness: Area normalized by the shape's dispersion (mu20 + mu02). Higher values indicate more compact and isotropic shapes.", false, true},
    {ECCENTRICITY, "ECCENTRICITY", "Eccentricity: Ratio of principal inertia eigenvalues (λ_1/λ_2). Measures elongation; values near 1 indicate circularity, higher values indicate elongation.", false, true},
    {LENGTH_MAJOR_AXIS, "LENGTH_MAJOR_AXIS", "Major axis length: Length of the longest diameter of the shape.", false, true},
    {LENGTH_MINOR_AXIS, "LENGTH_MINOR_AXIS", "Minor axis length: Length of the shortest diameter of the shape.", false, true},
    {AXIS_ORIENTATION, "AXIS_ORIENTATION", "Axis orientation: Angle of the principal inertia axis, computed from central moments. Indicates the dominant orientation of the shape.", false, true},
    {CIRCULARITY, "CIRCULARITY", "Circularity: Ratio of the minor to major eigenvalues of the inertia matrix (λ_2/λ_1), i.e., Inverse of eccentricity. Indicates how circular a shape is; values near 1 suggest circularity, values near 0 indicate elongation.", false, true},

    {BITQUADS_AREA, "BITQUADS_AREA", "BitQuads area (Duda): Refined sub-pixel area estimation using fractional weights based on the geometric contribution of local 2x2 pixel patterns.", false, true},
    {BITQUADS_NUMBER_EULER, "BITQUADS_NUMBER_EULER", "BitQuads Euler number: Topological invariant computed as the number of connected components minus the number of holes, using 2x2 pattern statistics under 4- or 8-connectivity.", false, true},
    {BITQUADS_NUMBER_HOLES, "BITQUADS_NUMBER_HOLES", "BitQuads number of holes: Number of interior holes in the component, derived from the Euler characteristic assuming a single connected object.", false, true},
    {BITQUADS_PERIMETER, "BITQUADS_PERIMETER", "BitQuads perimeter: Discrete approximation of the shape's boundary length, calculated by summing edge-contributing patterns in the 2x2 pixel grid.", false, true},
    {BITQUADS_PERIMETER_CONTINUOUS, "BITQUADS_PERIMETER_CONTINUOUS", "BitQuads continuous perimeter: Smoothed estimation of the boundary length, incorporating weighted transitions across pixel edges and diagonals.", false, true},
    {BITQUADS_CIRCULARITY, "BITQUADS_CIRCULARITY", "BitQuads circularity: Compactness measure defined as (4π x areaDuda) / perimeter². Values close to 1 indicate circular shapes; lower values suggest elongation or irregularity.", false, true},
    {BITQUADS_PERIMETER_AVERAGE, "BITQUADS_PERIMETER_AVERAGE", "BitQuads average perimeter: Mean perimeter per connected component, accounting for complex structures and holes.", false, true},
    {BITQUADS_LENGTH_AVERAGE, "BITQUADS_LENGTH_AVERAGE", "BitQuads average length: Estimated average longitudinal extent per component, derived from the average perimeter.", false, true},
    {BITQUADS_WIDTH_AVERAGE, "BITQUADS_WIDTH_AVERAGE", "BitQuads average width: Estimated transverse extent per component, computed as (2 x average area) / average perimeter.", false, true},

    {HEIGHT_NODE, "HEIGHT_NODE", "Height: Longest path from this node to any leaf in its subtree. Measures the depth of the subtree rooted at the node.", false, true},
    {DEPTH_NODE, "DEPTH_NODE", "Depth: Number of steps from this node to the root of the tree. Indicates the level of embedding within the tree hierarchy.", false, true},
    {IS_LEAF_NODE, "IS_LEAF_NODE", "Is leaf: True if the node has no children, i.e., it represents a minimal component in the hierarchy.", false, true},
    {IS_ROOT_NODE, "IS_ROOT_NODE", "Is root: True if the node is the root of the tree, representing the entire image support.", false, true},
    {NUM_CHILDREN_NODE, "NUM_CHILDREN_NODE", "Number of children: Count of direct child nodes. Reflects the immediate branching factor of the node.", false, true},
    {NUM_SIBLINGS_NODE, "NUM_SIBLINGS_NODE", "Number of siblings: Number of other nodes that share the same parent.", false, true},
    {NUM_DESCENDANTS_NODE, "NUM_DESCENDANTS_NODE", "Number of descendants: Total number of nodes in the subtree rooted at this node (excluding the node itself).", false, true},
    {NUM_LEAF_DESCENDANTS_NODE, "NUM_LEAF_DESCENDANTS_NODE", "Number of leaf descendants: Number of leaf nodes in the subtree. Reflects the number of minimal patterns under this structure.", false, true},
    {LEAF_RATIO_NODE, "LEAF_RATIO_NODE", "Leaf ratio: Ratio of leaf descendants to total descendants. Measures structural 'flatness' or terminal density of the subtree.", false, true},
    {BALANCE_NODE, "BALANCE_NODE", "Balance: Difference between the maximum and minimum heights among the subtrees of the children. Indicates branching symmetry.", false, true},

    {MAX_DIST, "MAX_DIST", "Maximum distance transform value of the node, computed using incremental contour extraction and Differential Image Foresting Transform.", true, false},

    {AVG_CHILD_HEIGHT_NODE, "AVG_CHILD_HEIGHT_NODE", "Average child height: Mean height of all direct child subtrees. Useful for measuring uniformity of the subtree structure.", false, true},

    {CONTOUR_PIXELS, "CONTOUR_PIXELS", "Contour pixels: Number of support pixels touching the 4-neighbour complement.", false, true},
    {CONTOUR_PERIMETER, "CONTOUR_PERIMETER", "Contour perimeter: 4-neighbour exposed-side perimeter of the support.", false, true},
    {CONTOUR_SIDE_NORTH, "CONTOUR_SIDE_NORTH", "Contour north sides: Number of exposed north sides over support pixels.", false, true},
    {CONTOUR_SIDE_WEST, "CONTOUR_SIDE_WEST", "Contour west sides: Number of exposed west sides over support pixels.", false, true},
    {CONTOUR_SIDE_EAST, "CONTOUR_SIDE_EAST", "Contour east sides: Number of exposed east sides over support pixels.", false, true},
    {CONTOUR_SIDE_SOUTH, "CONTOUR_SIDE_SOUTH", "Contour south sides: Number of exposed south sides over support pixels.", false, true}
}};

/**
 * @brief Returns metadata for a scalar attribute, or `nullptr` for an invalid enum value.
 */
inline const AttributeMetadata* metadata(Attribute attribute) noexcept {
    const auto index = static_cast<std::size_t>(attribute);
    if (index >= ATTRIBUTE_METADATA.size()) {
        return nullptr;
    }
    const AttributeMetadata& item = ATTRIBUTE_METADATA[index];
    return item.attribute == attribute ? &item : nullptr;
}

/**
 * @brief Returns the stable symbolic name of `attribute`.
 *
 * Unknown values are reported as `"UNKNOWN"` rather than throwing so debugging
 * and diagnostics can still produce a string.
 */
inline std::string_view name(Attribute attribute) noexcept {
    const AttributeMetadata* item = metadata(attribute);
    return item != nullptr ? item->name : std::string_view("UNKNOWN");
}

/**
 * @brief Returns the user-facing description of `attribute`.
 */
inline std::string_view description(Attribute attribute) noexcept {
    const AttributeMetadata* item = metadata(attribute);
    return item != nullptr ? item->description : std::string_view("Unknown attribute.");
}

/**
 * @brief Tests whether `attribute` requires a weighted altitude contract.
 */
inline bool requiresAltitude(Attribute attribute) noexcept {
    const AttributeMetadata* item = metadata(attribute);
    return item != nullptr && item->requiresAltitude;
}

/**
 * @brief Tests whether `attribute` can be computed from support/topology alone.
 */
inline bool isTopologyOnly(Attribute attribute) noexcept {
    const AttributeMetadata* item = metadata(attribute);
    return item != nullptr && item->topologyOnly;
}

/**
 * @brief Tests whether `attribute` belongs to the altitude-aware pipeline path.
 *
 * `AREA` is always accepted by the pipeline because it is a common dependency
 * for altitude-aware attributes even though it does not read altitude itself.
 */
inline bool isAttributePipelineAltitudeAttribute(Attribute attribute) noexcept {
    return attribute == AREA || requiresAltitude(attribute);
}

/**
 * @brief Tests whether the ordinary attribute pipeline can materialize `attribute`.
 */
inline bool isPipelineComputed(Attribute attribute) noexcept {
    return attribute == AREA || requiresAltitude(attribute) || isTopologyOnly(attribute);
}

/**
 * @brief Parses a stable symbolic attribute name.
 *
 * Matching is exact and case-sensitive because the same names are used as public
 * dictionary keys in the Python bindings.
 */
inline std::optional<Attribute> parse(std::string_view nameToFind) noexcept {
    for (const AttributeMetadata& item : ATTRIBUTE_METADATA) {
        if (item.name == nameToFind) {
            return item.attribute;
        }
    }
    return std::nullopt;
}

/**
 * @brief Returns the canonical expansion of every public attribute group.
 *
 * The returned map is process-static and must be treated as read-only metadata.
 * Group expansion preserves the order listed here before request-level
 * deduplication is applied by higher-level helpers.
 */
inline const std::unordered_map<AttributeGroup, std::vector<Attribute>>& attributeGroups() {
    static const std::unordered_map<AttributeGroup, std::vector<Attribute>> groups = [] {
        std::vector<Attribute> all;
        all.reserve(ATTRIBUTE_METADATA.size());
        for (const AttributeMetadata& item : ATTRIBUTE_METADATA) {
            all.push_back(item.attribute);
        }

        return std::unordered_map<AttributeGroup, std::vector<Attribute>>{
            {AttributeGroup::GRAY_LEVEL, {
                VOLUME,
                RELATIVE_VOLUME,
                LEVEL,
                GRAY_HEIGHT,
                MEAN_LEVEL,
                VARIANCE_LEVEL
            }},
            {AttributeGroup::SHAPE, {
                AREA,
                BOX_WIDTH,
                BOX_HEIGHT,
                DIAGONAL_LENGTH,
                RECTANGULARITY,
                RATIO_WH,
                BOX_COL_MIN,
                BOX_COL_MAX,
                BOX_ROW_MIN,
                BOX_ROW_MAX,
                CENTRAL_MOMENT_20,
                CENTRAL_MOMENT_02,
                CENTRAL_MOMENT_11,
                CENTRAL_MOMENT_30,
                CENTRAL_MOMENT_03,
                CENTRAL_MOMENT_21,
                CENTRAL_MOMENT_12,
                HU_MOMENT_1,
                HU_MOMENT_2,
                HU_MOMENT_3,
                HU_MOMENT_4,
                HU_MOMENT_5,
                HU_MOMENT_6,
                HU_MOMENT_7,
                INERTIA,
                COMPACTNESS,
                ECCENTRICITY,
                LENGTH_MAJOR_AXIS,
                LENGTH_MINOR_AXIS,
                AXIS_ORIENTATION,
                CIRCULARITY,
                BITQUADS_AREA,
                BITQUADS_NUMBER_EULER,
                BITQUADS_NUMBER_HOLES,
                BITQUADS_PERIMETER,
                BITQUADS_PERIMETER_CONTINUOUS,
                BITQUADS_CIRCULARITY,
                BITQUADS_PERIMETER_AVERAGE,
                BITQUADS_LENGTH_AVERAGE,
                BITQUADS_WIDTH_AVERAGE,
                MAX_DIST,
                CONTOUR_PIXELS,
                CONTOUR_PERIMETER,
                CONTOUR_SIDE_NORTH,
                CONTOUR_SIDE_WEST,
                CONTOUR_SIDE_EAST,
                CONTOUR_SIDE_SOUTH
            }},
            {AttributeGroup::MOMENTS, {
                CENTRAL_MOMENT_20,
                CENTRAL_MOMENT_02,
                CENTRAL_MOMENT_11,
                CENTRAL_MOMENT_30,
                CENTRAL_MOMENT_03,
                CENTRAL_MOMENT_21,
                CENTRAL_MOMENT_12,
                HU_MOMENT_1,
                HU_MOMENT_2,
                HU_MOMENT_3,
                HU_MOMENT_4,
                HU_MOMENT_5,
                HU_MOMENT_6,
                HU_MOMENT_7,
                INERTIA,
                COMPACTNESS,
                ECCENTRICITY,
                LENGTH_MAJOR_AXIS,
                LENGTH_MINOR_AXIS,
                AXIS_ORIENTATION,
                CIRCULARITY
            }},
            {AttributeGroup::BOUNDARY, {
                BITQUADS_AREA,
                BITQUADS_NUMBER_EULER,
                BITQUADS_NUMBER_HOLES,
                BITQUADS_PERIMETER,
                BITQUADS_PERIMETER_CONTINUOUS,
                BITQUADS_CIRCULARITY,
                BITQUADS_PERIMETER_AVERAGE,
                BITQUADS_LENGTH_AVERAGE,
                BITQUADS_WIDTH_AVERAGE,
                CONTOUR_PIXELS,
                CONTOUR_PERIMETER,
                CONTOUR_SIDE_NORTH,
                CONTOUR_SIDE_WEST,
                CONTOUR_SIDE_EAST,
                CONTOUR_SIDE_SOUTH
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
            {AttributeGroup::ALL, std::move(all)}
        };
    }();
    return groups;
}

} // namespace mmcfilters::attributes::registry
