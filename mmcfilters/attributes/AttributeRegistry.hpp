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
 * @brief Adjacency capability required by one scalar attribute.
 */
enum class AttributeAdjacencyRequirement { None, UniformOrDirectional, Uniform };

/**
 * @brief Runtime capabilities required to compute one scalar attribute.
 *
 * Altitude availability is distinct from altitude monotonicity. BitQuad
 * projection, for example, is topology-only for uniform adjacency but needs
 * altitude to select between unequal directional adjacencies.
 */
struct AttributeCapabilityRequirements {
    /// Whether node altitudes must be available.
    bool altitude = false;
    /// Whether proper parts must have a regular 2D domain.
    bool gridDomain2D = false;
    /// Required form of adjacency context.
    AttributeAdjacencyRequirement adjacency = AttributeAdjacencyRequirement::None;
    /// Whether the global altitude order must be monotone.
    bool monotoneAltitudeOrder = false;
    /// Whether directional adjacency selection requires node altitudes.
    bool altitudeForDirectionalAdjacency = false;
    /// Whether only canonical 4- or 8-neighbourhoods are accepted.
    bool canonical4Or8Adjacency = false;

    /**
     * @brief Compares every capability requirement.
     *
     * @return True when the documented condition holds; otherwise false.
     */
    constexpr bool operator==(const AttributeCapabilityRequirements&) const noexcept = default;
};

inline constexpr AttributeCapabilityRequirements NO_REQUIREMENTS{};

inline constexpr AttributeCapabilityRequirements ALTITUDE_REQUIREMENTS{.altitude = true};

inline constexpr AttributeCapabilityRequirements GRID_DOMAIN_2D_REQUIREMENTS{.gridDomain2D = true};

inline constexpr AttributeCapabilityRequirements BITQUAD_REQUIREMENTS{.gridDomain2D = true,
                                                                      .adjacency = AttributeAdjacencyRequirement::UniformOrDirectional,
                                                                      .altitudeForDirectionalAdjacency = true,
                                                                      .canonical4Or8Adjacency = true};

inline constexpr AttributeCapabilityRequirements MAX_DIST_REQUIREMENTS{.gridDomain2D = true};

inline constexpr AttributeCapabilityRequirements MAX_DIST_EXACT_REQUIREMENTS{.gridDomain2D = true};

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

    /// Whether the descriptor can be computed from topology/support alone.
    bool topologyOnly = false;

    /// Complete runtime capability contract used by scheduling and validation.
    AttributeCapabilityRequirements requirements{};
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
inline constexpr std::array<AttributeMetadata, static_cast<std::size_t>(Attribute::MaxSquaredDistExact) + 1> ATTRIBUTE_METADATA{
    {{Area, "AREA", "Area: Number of pixels in the connected component.", false, NO_REQUIREMENTS},
     {Volume, "VOLUME",
      "Volume: Sum of the gray-level intensities of all pixels in the connected component. Interpreted as the total mass under the component, or the integral "
      "of the image function over its support.",
      false, ALTITUDE_REQUIREMENTS},
     {RelativeVolume, "RELATIVE_VOLUME",
      "Relative volume: Sum of differences between the node altitude and the gray levels of pixels in the component. Measures the amount of intensity required "
      "to fill the component to its node altitude.",
      false, ALTITUDE_REQUIREMENTS},
     {GrayLevelHeight, "GRAY_LEVEL_HEIGHT",
      "Gray-level height: Maximum absolute altitude difference between a node and any node in its subtree. This reduces to the traditional one-sided span on "
      "monotone max/min trees and also applies to hierarchies with unconstrained altitude order. Leaves have gray-level height 0.",
      false, ALTITUDE_REQUIREMENTS},
     {MeanGrayLevel, "MEAN_GRAY_LEVEL", "Mean gray level: Arithmetic mean of the image values over the node support.", false, ALTITUDE_REQUIREMENTS},
     {GrayLevelVariance, "GRAY_LEVEL_VARIANCE", "Gray-level variance: Population variance of the image values over the node support.", false,
      ALTITUDE_REQUIREMENTS},

     {BoxWidth, "BOX_WIDTH", "Bounding box width: Width of the minimum rectangle enclosing the connected component.", true, GRID_DOMAIN_2D_REQUIREMENTS},
     {BoundingBoxHeight, "BOUNDING_BOX_HEIGHT", "Bounding-box height: Height of the minimum rectangle enclosing the connected component.", true,
      GRID_DOMAIN_2D_REQUIREMENTS},
     {DiagonalLength, "DIAGONAL_LENGTH", "Diagonal length: Euclidean length of the diagonal of the bounding box, computed as sqrt(width^2 + height^2).", true,
      GRID_DOMAIN_2D_REQUIREMENTS},
     {Rectangularity, "RECTANGULARITY",
      "Rectangularity: Ratio between the connected component area and the area of its bounding box. Values closer to 1 indicate shapes that efficiently fill "
      "their bounding box.",
      true, GRID_DOMAIN_2D_REQUIREMENTS},
     {RatioWh, "RATIO_WH", "Aspect ratio: Ratio of the bounding box width to its height. Describes the elongation of the component.", true,
      GRID_DOMAIN_2D_REQUIREMENTS},
     {BoxColumnMin, "BOX_COLUMN_MIN", "Bounding box column min: Minimum column index covered by the connected component.", true, GRID_DOMAIN_2D_REQUIREMENTS},
     {BoxColumnMax, "BOX_COLUMN_MAX", "Bounding box column max: Maximum column index covered by the connected component.", true, GRID_DOMAIN_2D_REQUIREMENTS},
     {BoxRowMin, "BOX_ROW_MIN", "Bounding box row min: Minimum row index covered by the connected component.", true, GRID_DOMAIN_2D_REQUIREMENTS},
     {BoxRowMax, "BOX_ROW_MAX", "Bounding box row max: Maximum row index covered by the connected component.", true, GRID_DOMAIN_2D_REQUIREMENTS},

     {CentralMoment20, "CENTRAL_MOMENT_20",
      "Central moment (2,0): Second-order moment about the centroid along the x-axis. Measures the horizontal spread of the component.", true,
      GRID_DOMAIN_2D_REQUIREMENTS},
     {CentralMoment02, "CENTRAL_MOMENT_02",
      "Central moment (0,2): Second-order moment about the centroid along the y-axis. Measures the vertical spread of the component.", true,
      GRID_DOMAIN_2D_REQUIREMENTS},
     {CentralMoment11, "CENTRAL_MOMENT_11",
      "Central moment (1,1): Mixed second-order moment about the centroid. Represents the covariance between x and y coordinates.", true,
      GRID_DOMAIN_2D_REQUIREMENTS},
     {CentralMoment30, "CENTRAL_MOMENT_30",
      "Central moment (3,0): Third-order moment about the centroid along the x-axis. Describes horizontal asymmetry of the component.", true,
      GRID_DOMAIN_2D_REQUIREMENTS},
     {CentralMoment03, "CENTRAL_MOMENT_03",
      "Central moment (0,3): Third-order moment about the centroid along the y-axis. Describes vertical asymmetry of the component.", true,
      GRID_DOMAIN_2D_REQUIREMENTS},
     {CentralMoment21, "CENTRAL_MOMENT_21",
      "Central moment (2,1): Mixed third-order moment about the centroid. Captures joint spread and asymmetry in x and y.", true, GRID_DOMAIN_2D_REQUIREMENTS},
     {CentralMoment12, "CENTRAL_MOMENT_12",
      "Central moment (1,2): Mixed third-order moment about the centroid. Captures joint spread and asymmetry in y and x.", true, GRID_DOMAIN_2D_REQUIREMENTS},

     {HuMoment1, "HU_MOMENT_1", "Hu moment 1: Invariant to translation, scale, and rotation. Represents overall spatial variance (shape dispersion).", true,
      GRID_DOMAIN_2D_REQUIREMENTS},
     {HuMoment2, "HU_MOMENT_2", "Hu moment 2: Invariant capturing the difference between horizontal and vertical spread.", true, GRID_DOMAIN_2D_REQUIREMENTS},
     {HuMoment3, "HU_MOMENT_3", "Hu moment 3: Sensitive to skewness and asymmetry in the pixel distribution.", true, GRID_DOMAIN_2D_REQUIREMENTS},
     {HuMoment4, "HU_MOMENT_4", "Hu moment 4: Measures symmetry with respect to diagonal axes.", true, GRID_DOMAIN_2D_REQUIREMENTS},
     {HuMoment5, "HU_MOMENT_5", "Hu moment 5: Descriptor sensitive to orientation and reflection; captures complex asymmetries.", true,
      GRID_DOMAIN_2D_REQUIREMENTS},
     {HuMoment6, "HU_MOMENT_6", "Hu moment 6: Invariant capturing elliptic asymmetries, sensitive to specific shape curvature.", true,
      GRID_DOMAIN_2D_REQUIREMENTS},
     {HuMoment7, "HU_MOMENT_7", "Hu moment 7: Highly sensitive to irregularities and fine variations; helps discriminate mirror-symmetric shapes.", true,
      GRID_DOMAIN_2D_REQUIREMENTS},

     {Inertia, "INERTIA",
      "Inertia: Sum of normalized second-order central moments (mu20 + mu02). Measures the dispersion of mass around the centroid. Higher values indicate "
      "objects with thin and elongated structures.",
      true, GRID_DOMAIN_2D_REQUIREMENTS},
     {Compactness, "COMPACTNESS",
      "Compactness: Area normalized by the shape's dispersion (mu20 + mu02). Higher values indicate more compact and isotropic shapes.", true,
      GRID_DOMAIN_2D_REQUIREMENTS},
     {Eccentricity, "ECCENTRICITY",
      "Eccentricity: Ratio of principal inertia eigenvalues (λ_1/λ_2). Measures elongation; values near 1 indicate circularity, higher values indicate "
      "elongation. Degenerate line-like supports saturate at a finite maximum.",
      true, GRID_DOMAIN_2D_REQUIREMENTS},
     {LengthMajorAxis, "LENGTH_MAJOR_AXIS", "Major axis length: Length of the longest diameter of the shape.", true, GRID_DOMAIN_2D_REQUIREMENTS},
     {LengthMinorAxis, "LENGTH_MINOR_AXIS", "Minor axis length: Length of the shortest diameter of the shape.", true, GRID_DOMAIN_2D_REQUIREMENTS},
     {AxisOrientation, "AXIS_ORIENTATION",
      "Axis orientation: Angle of the principal inertia axis, computed from central moments. Indicates the dominant orientation of the shape.", true,
      GRID_DOMAIN_2D_REQUIREMENTS},
     {Circularity, "CIRCULARITY",
      "Circularity: Ratio of the minor to major eigenvalues of the inertia matrix (λ_2/λ_1), i.e., Inverse of eccentricity. Indicates how circular a shape is; "
      "values near 1 suggest circularity, values near 0 indicate elongation.",
      true, GRID_DOMAIN_2D_REQUIREMENTS},

     {BitquadArea, "BITQUAD_AREA",
      "Bitquad area (Duda): Refined sub-pixel area estimation using fractional weights based on the geometric contribution of local 2x2 pixel patterns.", true,
      BITQUAD_REQUIREMENTS},
     {BitquadNumberEuler, "BITQUAD_NUMBER_EULER",
      "Bitquad Euler number: Topological invariant computed as the number of connected components minus the number of holes, using 2x2 pattern statistics "
      "under 4- or 8-connectivity.",
      true, BITQUAD_REQUIREMENTS},
     {BitquadNumberHoles, "BITQUAD_NUMBER_HOLES",
      "Bitquad number of holes: Number of interior holes in the component, derived from the Euler characteristic assuming a single connected object.", true,
      BITQUAD_REQUIREMENTS},
     {BitquadPerimeter, "BITQUAD_PERIMETER",
      "Bitquad perimeter: Discrete approximation of the shape's boundary length, calculated by summing edge-contributing patterns in the 2x2 pixel grid.",
      true, BITQUAD_REQUIREMENTS},
     {BitquadPerimeterContinuous, "BITQUAD_PERIMETER_CONTINUOUS",
      "Bitquad continuous perimeter: Smoothed estimation of the boundary length, incorporating valuedTree transitions across pixel edges and diagonals.", true,
      BITQUAD_REQUIREMENTS},
     {BitquadCircularity, "BITQUAD_CIRCULARITY",
      "Bitquad circularity: Compactness measure defined as (4π x areaDuda) / perimeter². Values close to 1 indicate circular shapes; lower values suggest "
      "elongation or irregularity. Degenerate zero-perimeter supports return 0.",
      true, BITQUAD_REQUIREMENTS},
     {BitquadPerimeterAverage, "BITQUAD_PERIMETER_AVERAGE",
      "Bitquad average perimeter: Mean perimeter per connected component, accounting for complex structures and holes. Non-positive Euler component estimates "
      "return 0.",
      true, BITQUAD_REQUIREMENTS},
     {BitquadLengthAverage, "BITQUAD_LENGTH_AVERAGE",
      "Bitquad average length: Estimated average longitudinal extent per component, derived from the average perimeter with a zero fallback for non-positive "
      "Euler component estimates.",
      true, BITQUAD_REQUIREMENTS},
     {BitquadWidthAverage, "BITQUAD_WIDTH_AVERAGE",
      "Bitquad average width: Estimated transverse extent per component, computed as (2 x areaDuda) / continuous perimeter with a zero fallback for "
      "degenerate perimeter.",
      true, BITQUAD_REQUIREMENTS},

     {SubtreeHeight, "SUBTREE_HEIGHT", "Subtree height: Longest path from this node to any leaf in its subtree, measured in tree edges.", true,
      NO_REQUIREMENTS},
     {DepthNode, "DEPTH_NODE", "Depth: Number of steps from this node to the root of the tree. Indicates the level of embedding within the tree hierarchy.",
      true, NO_REQUIREMENTS},
     {IsLeafNode, "IS_LEAF_NODE", "Is leaf: True if the node has no children, i.e., it represents a minimal component in the hierarchy.", true,
      NO_REQUIREMENTS},
     {IsRootNode, "IS_ROOT_NODE", "Is root: True if the node is the root of the tree, representing the entire image support.", true, NO_REQUIREMENTS},
     {NumChildrenNode, "NUM_CHILDREN_NODE", "Number of children: Count of direct child nodes. Reflects the immediate branching factor of the node.", true,
      NO_REQUIREMENTS},
     {NumSiblingsNode, "NUM_SIBLINGS_NODE", "Number of siblings: Number of other nodes that share the same parent.", true, NO_REQUIREMENTS},
     {NumDescendantsNode, "NUM_DESCENDANTS_NODE",
      "Number of descendants: Total number of nodes in the subtree rooted at this node (excluding the node itself).", true, NO_REQUIREMENTS},
     {NumLeafDescendantsNode, "NUM_LEAF_DESCENDANTS_NODE",
      "Number of leaf descendants: Number of leaf nodes in the subtree. Reflects the number of minimal patterns under this structure.", true, NO_REQUIREMENTS},
     {LeafRatioNode, "LEAF_RATIO_NODE",
      "Leaf ratio: Ratio of leaf descendants to total descendants. Measures structural 'flatness' or terminal density of the subtree.", true, NO_REQUIREMENTS},
     {BalanceNode, "BALANCE_NODE",
      "Balance: Difference between the maximum and minimum heights among the subtrees of the children. Indicates branching symmetry.", true, NO_REQUIREMENTS},

     {MaxDist, "MAX_DIST",
      "Approximate maximum Euclidean distance from the foreground A4 contour over the node support, in pixels, computed by an adaptive A8 dynamic "
      "image foresting transform.",
      true, MAX_DIST_REQUIREMENTS},

     {AvgChildHeightNode, "AVG_CHILD_HEIGHT_NODE",
      "Average child height: Mean height of all direct child subtrees. Useful for measuring uniformity of the subtree structure.", true, NO_REQUIREMENTS},

     {ContourPixels, "CONTOUR_PIXELS", "Contour pixels: Number of support pixels touching the 4-neighbour complement.", true, GRID_DOMAIN_2D_REQUIREMENTS},
     {ContourPerimeter, "CONTOUR_PERIMETER", "Contour perimeter: 4-neighbour exposed-side perimeter of the support.", true, GRID_DOMAIN_2D_REQUIREMENTS},
     {ContourSideNorth, "CONTOUR_SIDE_NORTH", "Contour north sides: Number of exposed north sides over support pixels.", true, GRID_DOMAIN_2D_REQUIREMENTS},
     {ContourSideWest, "CONTOUR_SIDE_WEST", "Contour west sides: Number of exposed west sides over support pixels.", true, GRID_DOMAIN_2D_REQUIREMENTS},
     {ContourSideEast, "CONTOUR_SIDE_EAST", "Contour east sides: Number of exposed east sides over support pixels.", true, GRID_DOMAIN_2D_REQUIREMENTS},
     {ContourSideSouth, "CONTOUR_SIDE_SOUTH", "Contour south sides: Number of exposed south sides over support pixels.", true, GRID_DOMAIN_2D_REQUIREMENTS},

     {MaxDistExact, "MAX_DIST_EXACT",
      "Maximum exact Euclidean distance from the foreground A4 contour over the node support in the original 2D pixel domain, in pixels.", true,
      MAX_DIST_EXACT_REQUIREMENTS},

     {DistSquaredSumExact, "DIST_SQUARED_SUM_EXACT",
      "Squared-distance sum: Sum over the node support of the exact squared Euclidean distance to the foreground A4 contour.", true,
      MAX_DIST_EXACT_REQUIREMENTS},
     {DistSquaredMeanExact, "DIST_SQUARED_MEAN_EXACT",
      "Squared-distance mean: Arithmetic mean over the node support of the exact squared Euclidean distance to the foreground A4 contour.", true,
      MAX_DIST_EXACT_REQUIREMENTS},
     {DistRmsExact, "DIST_RMS_EXACT",
      "Distance RMS: Square root of the mean exact squared Euclidean distance to the foreground A4 contour, in pixel units.", true,
      MAX_DIST_EXACT_REQUIREMENTS},
     {DistSquaredVarianceExact, "DIST_SQUARED_VARIANCE_EXACT",
      "Squared-distance variance: Population variance over the node support of the exact squared Euclidean distance to the foreground A4 contour.", true,
      MAX_DIST_EXACT_REQUIREMENTS},

     {DistSquaredSum, "DIST_SQUARED_SUM",
      "Approximate squared-distance sum over the node support, computed from the adaptive A8 dynamic image foresting transform.", true,
      MAX_DIST_REQUIREMENTS},
     {DistSquaredMean, "DIST_SQUARED_MEAN",
      "Approximate squared-distance mean over the node support, computed from the adaptive A8 dynamic image foresting transform.", true,
      MAX_DIST_REQUIREMENTS},
     {DistRms, "DIST_RMS",
      "Approximate distance RMS in pixel units, computed from the adaptive A8 dynamic image foresting transform.", true,
      MAX_DIST_REQUIREMENTS},
     {DistSquaredVariance, "DIST_SQUARED_VARIANCE",
      "Approximate population variance of squared distances over the node support, computed from the adaptive A8 dynamic image foresting transform.", true,
      MAX_DIST_REQUIREMENTS},

     {MaxDistCenterRowExact, "MAX_DIST_CENTER_ROW_EXACT",
      "Maximum-distance center row: Zero-based row of the smallest row-major support pixel attaining the exact maximum squared contour distance.", true,
      MAX_DIST_EXACT_REQUIREMENTS},
     {MaxDistCenterColumnExact, "MAX_DIST_CENTER_COLUMN_EXACT",
      "Maximum-distance center column: Zero-based column of the smallest row-major support pixel attaining the exact maximum squared contour distance.",
      true, MAX_DIST_EXACT_REQUIREMENTS},
     {MaxDistCenterRow, "MAX_DIST_CENTER_ROW",
      "Approximate maximum-distance center row: Zero-based row of the smallest row-major support pixel attaining the adaptive-A8 DIFT maximum.", true,
      MAX_DIST_REQUIREMENTS},
     {MaxDistCenterColumn, "MAX_DIST_CENTER_COLUMN",
      "Approximate maximum-distance center column: Zero-based column of the smallest row-major support pixel attaining the adaptive-A8 DIFT maximum.",
      true, MAX_DIST_REQUIREMENTS},

     {MaxDistPlateauAreaExact, "MAX_DIST_PLATEAU_AREA_EXACT",
      "Maximum-distance plateau area: Number of support pixels attaining the exact maximum squared contour distance.", true,
      MAX_DIST_EXACT_REQUIREMENTS},
     {MaxDistPlateauCentroidRowExact, "MAX_DIST_PLATEAU_CENTROID_ROW_EXACT",
      "Maximum-distance plateau centroid row: Arithmetic mean of the row coordinates of all exact maximum-distance pixels.", true,
      MAX_DIST_EXACT_REQUIREMENTS},
     {MaxDistPlateauCentroidColumnExact, "MAX_DIST_PLATEAU_CENTROID_COLUMN_EXACT",
      "Maximum-distance plateau centroid column: Arithmetic mean of the column coordinates of all exact maximum-distance pixels.", true,
      MAX_DIST_EXACT_REQUIREMENTS},
     {MaxDistPlateauArea, "MAX_DIST_PLATEAU_AREA",
      "Approximate maximum-distance plateau area: Number of support pixels attaining the adaptive-A8 DIFT maximum.", true,
      MAX_DIST_REQUIREMENTS},
     {MaxDistPlateauCentroidRow, "MAX_DIST_PLATEAU_CENTROID_ROW",
      "Approximate maximum-distance plateau centroid row over all support pixels attaining the adaptive-A8 DIFT maximum.", true,
      MAX_DIST_REQUIREMENTS},
     {MaxDistPlateauCentroidColumn, "MAX_DIST_PLATEAU_CENTROID_COLUMN",
      "Approximate maximum-distance plateau centroid column over all support pixels attaining the adaptive-A8 DIFT maximum.", true,
      MAX_DIST_REQUIREMENTS},

     {DistSum, "DIST_SUM", "Sum of approximate Euclidean contour distances over the node support, in pixels.", true, MAX_DIST_REQUIREMENTS},
     {DistMean, "DIST_MEAN", "Arithmetic mean of approximate Euclidean contour distances over the node support, in pixels.", true,
      MAX_DIST_REQUIREMENTS},
     {DistVariance, "DIST_VARIANCE", "Population variance of approximate Euclidean contour distances, in squared pixels.", true,
      MAX_DIST_REQUIREMENTS},
     {DistMedian, "DIST_MEDIAN", "Lower median of the approximate Euclidean contour-distance distribution, in pixels.", true,
      MAX_DIST_REQUIREMENTS},
     {DistMode, "DIST_MODE", "Smallest approximate Euclidean contour distance among the most frequent distance levels, in pixels.", true,
      MAX_DIST_REQUIREMENTS},
     {DistQ25, "DIST_Q25", "Lower 25th percentile of the approximate Euclidean contour-distance distribution, in pixels.", true,
      MAX_DIST_REQUIREMENTS},
     {DistQ75, "DIST_Q75", "Lower 75th percentile of the approximate Euclidean contour-distance distribution, in pixels.", true,
      MAX_DIST_REQUIREMENTS},
     {DistQ90, "DIST_Q90", "Lower 90th percentile of the approximate Euclidean contour-distance distribution, in pixels.", true,
      MAX_DIST_REQUIREMENTS},
     {DistEntropy, "DIST_ENTROPY", "Shannon entropy in bits of the normalized approximate squared-distance histogram.", true,
      MAX_DIST_REQUIREMENTS},
     {DistPositiveArea, "DIST_POSITIVE_AREA", "Number of support pixels with strictly positive approximate contour distance.", true,
      MAX_DIST_REQUIREMENTS},
     {DistLevelCount, "DIST_LEVEL_COUNT", "Number of distinct approximate squared-distance levels, including zero.", true,
      MAX_DIST_REQUIREMENTS},
     {DistWeightedCentroidRow, "DIST_WEIGHTED_CENTROID_ROW",
      "Zero-based support centroid row weighted by approximate Euclidean contour distance.", true, MAX_DIST_REQUIREMENTS},
     {DistWeightedCentroidColumn, "DIST_WEIGHTED_CENTROID_COLUMN",
      "Zero-based support centroid column weighted by approximate Euclidean contour distance.", true, MAX_DIST_REQUIREMENTS},
     {DistWeightedCentralMoment20, "DIST_WEIGHTED_CENTRAL_MOMENT_20",
      "Unnormalized second-order column central moment weighted by approximate Euclidean contour distance.", true, MAX_DIST_REQUIREMENTS},
     {DistWeightedCentralMoment02, "DIST_WEIGHTED_CENTRAL_MOMENT_02",
      "Unnormalized second-order row central moment weighted by approximate Euclidean contour distance.", true, MAX_DIST_REQUIREMENTS},
     {DistWeightedCentralMoment11, "DIST_WEIGHTED_CENTRAL_MOMENT_11",
      "Unnormalized mixed row-column central moment weighted by approximate Euclidean contour distance.", true, MAX_DIST_REQUIREMENTS},
     {DistWeightedAxisOrientation, "DIST_WEIGHTED_AXIS_ORIENTATION",
      "Non-negative column-axis orientation in degrees of the approximate distance-weighted support moments.", true, MAX_DIST_REQUIREMENTS},
     {DistWeightedEccentricity, "DIST_WEIGHTED_ECCENTRICITY",
      "Major/minor eigenvalue ratio of the approximate distance-weighted second moments.", true, MAX_DIST_REQUIREMENTS},

     {DistSumExact, "DIST_SUM_EXACT", "Sum of exact Euclidean contour distances over the node support, in pixels.", true,
      MAX_DIST_EXACT_REQUIREMENTS},
     {DistMeanExact, "DIST_MEAN_EXACT", "Arithmetic mean of exact Euclidean contour distances over the node support, in pixels.", true,
      MAX_DIST_EXACT_REQUIREMENTS},
     {DistVarianceExact, "DIST_VARIANCE_EXACT", "Population variance of exact Euclidean contour distances, in squared pixels.", true,
      MAX_DIST_EXACT_REQUIREMENTS},
     {DistMedianExact, "DIST_MEDIAN_EXACT", "Lower median of the exact Euclidean contour-distance distribution, in pixels.", true,
      MAX_DIST_EXACT_REQUIREMENTS},
     {DistModeExact, "DIST_MODE_EXACT", "Smallest exact Euclidean contour distance among the most frequent distance levels, in pixels.", true,
      MAX_DIST_EXACT_REQUIREMENTS},
     {DistQ25Exact, "DIST_Q25_EXACT", "Lower 25th percentile of the exact Euclidean contour-distance distribution, in pixels.", true,
      MAX_DIST_EXACT_REQUIREMENTS},
     {DistQ75Exact, "DIST_Q75_EXACT", "Lower 75th percentile of the exact Euclidean contour-distance distribution, in pixels.", true,
      MAX_DIST_EXACT_REQUIREMENTS},
     {DistQ90Exact, "DIST_Q90_EXACT", "Lower 90th percentile of the exact Euclidean contour-distance distribution, in pixels.", true,
      MAX_DIST_EXACT_REQUIREMENTS},
     {DistEntropyExact, "DIST_ENTROPY_EXACT", "Shannon entropy in bits of the normalized exact squared-distance histogram.", true,
      MAX_DIST_EXACT_REQUIREMENTS},
     {DistPositiveAreaExact, "DIST_POSITIVE_AREA_EXACT", "Number of support pixels with strictly positive exact contour distance.", true,
      MAX_DIST_EXACT_REQUIREMENTS},
     {DistLevelCountExact, "DIST_LEVEL_COUNT_EXACT", "Number of distinct exact squared-distance levels, including zero.", true,
      MAX_DIST_EXACT_REQUIREMENTS},
     {DistWeightedCentroidRowExact, "DIST_WEIGHTED_CENTROID_ROW_EXACT",
      "Zero-based support centroid row weighted by exact Euclidean contour distance.", true, MAX_DIST_EXACT_REQUIREMENTS},
     {DistWeightedCentroidColumnExact, "DIST_WEIGHTED_CENTROID_COLUMN_EXACT",
      "Zero-based support centroid column weighted by exact Euclidean contour distance.", true, MAX_DIST_EXACT_REQUIREMENTS},
     {DistWeightedCentralMoment20Exact, "DIST_WEIGHTED_CENTRAL_MOMENT_20_EXACT",
      "Unnormalized second-order column central moment weighted by exact Euclidean contour distance.", true, MAX_DIST_EXACT_REQUIREMENTS},
     {DistWeightedCentralMoment02Exact, "DIST_WEIGHTED_CENTRAL_MOMENT_02_EXACT",
      "Unnormalized second-order row central moment weighted by exact Euclidean contour distance.", true, MAX_DIST_EXACT_REQUIREMENTS},
     {DistWeightedCentralMoment11Exact, "DIST_WEIGHTED_CENTRAL_MOMENT_11_EXACT",
      "Unnormalized mixed row-column central moment weighted by exact Euclidean contour distance.", true, MAX_DIST_EXACT_REQUIREMENTS},
     {DistWeightedAxisOrientationExact, "DIST_WEIGHTED_AXIS_ORIENTATION_EXACT",
      "Non-negative column-axis orientation in degrees of the exact distance-weighted support moments.", true, MAX_DIST_EXACT_REQUIREMENTS},
     {DistWeightedEccentricityExact, "DIST_WEIGHTED_ECCENTRICITY_EXACT",
      "Major/minor eigenvalue ratio of the exact distance-weighted second moments.", true, MAX_DIST_EXACT_REQUIREMENTS},

     {MaxSquaredDist, "MAX_SQUARED_DIST",
      "Approximate maximum squared Euclidean distance from the foreground A4 contour over the node support, in squared pixels, computed by an adaptive "
      "A8 dynamic image foresting transform.",
      true, MAX_DIST_REQUIREMENTS},
     {MaxSquaredDistExact, "MAX_SQUARED_DIST_EXACT",
      "Maximum exact squared Euclidean distance from the foreground A4 contour over the node support in the original 2D pixel domain, in squared pixels.",
      true, MAX_DIST_EXACT_REQUIREMENTS}}};

/**
 * @brief Tests whether ordinally aligned attribute metadata holds.
 *
 * @return True when ordinally aligned attribute metadata; otherwise false.
 */
inline constexpr bool hasOrdinallyAlignedAttributeMetadata() noexcept {
    for (std::size_t index = 0; index < ATTRIBUTE_METADATA.size(); ++index) {
        if (static_cast<std::size_t>(ATTRIBUTE_METADATA[index].attribute) != index) {
            return false;
        }
    }
    return true;
}

static_assert(hasOrdinallyAlignedAttributeMetadata(), "Attribute metadata rows must match Attribute enum ordinals.");

/**
 * @brief Returns metadata for a scalar attribute, or `nullptr` for an invalid enum value.
 *
 * @param attribute Attribute requested by the operation.
 * @return Metadata for a scalar attribute, or nullptr for an invalid enum value.
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
 *
 * @param attribute Attribute requested by the operation.
 * @return The stable symbolic name of attribute.
 */
inline std::string_view name(Attribute attribute) noexcept {
    const AttributeMetadata* item = metadata(attribute);
    return item != nullptr ? item->name : std::string_view("UNKNOWN");
}

/**
 * @brief Returns the user-facing description of `attribute`.
 *
 * @param attribute Attribute requested by the operation.
 * @return The user-facing description of attribute.
 */
inline std::string_view description(Attribute attribute) noexcept {
    const AttributeMetadata* item = metadata(attribute);
    return item != nullptr ? item->description : std::string_view("Unknown attribute.");
}

/**
 * @brief Tests whether `attribute` requires a valuedTree altitude contract.
 *
 * @param attribute Attribute requested by the operation.
 * @return True if attribute requires a valuedTree altitude contract; otherwise false.
 */
inline bool requiresAltitude(Attribute attribute) noexcept {
    const AttributeMetadata* item = metadata(attribute);
    return item != nullptr && item->requirements.altitude;
}

/**
 * @brief Tests whether `attribute` can be computed from support/topology alone.
 *
 * @param attribute Attribute requested by the operation.
 * @return True if attribute can be computed from support/topology alone; otherwise false.
 */
inline bool isTopologyOnly(Attribute attribute) noexcept {
    const AttributeMetadata* item = metadata(attribute);
    return item != nullptr && item->topologyOnly;
}

/**
 * @brief Returns the complete capability contract of one scalar attribute.
 *
 * This registry-level function is the canonical classification used by the
 * scheduler, diagnostics, bindings, tests, and documentation.
 *
 * @param attribute Attribute requested by the operation.
 * @return The complete capability contract of one scalar attribute.
 */
inline AttributeCapabilityRequirements capabilityRequirements(Attribute attribute) noexcept {
    const AttributeMetadata* item = metadata(attribute);
    return item != nullptr ? item->requirements : AttributeCapabilityRequirements{};
}

/**
 * @brief Tests whether `attribute` belongs to the altitude-aware pipeline path.
 *
 * `AREA` is always accepted by the pipeline because it is a common dependency
 * for altitude-aware attributes even though it does not read altitude itself.
 *
 * @param attribute Attribute requested by the operation.
 * @return True if attribute belongs to the altitude-aware pipeline path; otherwise false.
 */
inline bool isAttributePipelineAltitudeAttribute(Attribute attribute) noexcept { return attribute == Area || requiresAltitude(attribute); }

/**
 * @brief Tests whether the ordinary attribute pipeline can materialize `attribute`.
 *
 * @param attribute Attribute requested by the operation.
 * @return True if the ordinary attribute pipeline can materialize attribute; otherwise false.
 */
inline bool isPipelineComputed(Attribute attribute) noexcept { return attribute == Area || requiresAltitude(attribute) || isTopologyOnly(attribute); }

/**
 * @brief Parses a stable symbolic attribute name.
 *
 * Matching is exact and case-sensitive because the same names are used as public
 * dictionary keys in the Python bindings.
 *
 * @param nameToFind Symbolic attribute name to locate.
 * @return The parsed stable symbolic attribute name.
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
 *
 * @return The canonical expansion of every public attribute group.
 */
inline const std::unordered_map<AttributeGroup, std::vector<Attribute>>& attributeGroups() {
    static const std::unordered_map<AttributeGroup, std::vector<Attribute>> groups = [] {
        std::vector<Attribute> all;
        all.reserve(ATTRIBUTE_METADATA.size());
        for (const AttributeMetadata& item : ATTRIBUTE_METADATA) {
            all.push_back(item.attribute);
        }

        return std::unordered_map<AttributeGroup, std::vector<Attribute>>{
            {AttributeGroup::GrayLevel, {Volume, RelativeVolume, GrayLevelHeight, MeanGrayLevel, GrayLevelVariance}},
            {AttributeGroup::Shape,
             {Area,
              BoxWidth,
              BoundingBoxHeight,
              DiagonalLength,
              Rectangularity,
              RatioWh,
              BoxColumnMin,
              BoxColumnMax,
              BoxRowMin,
              BoxRowMax,
              CentralMoment20,
              CentralMoment02,
              CentralMoment11,
              CentralMoment30,
              CentralMoment03,
              CentralMoment21,
              CentralMoment12,
              HuMoment1,
              HuMoment2,
              HuMoment3,
              HuMoment4,
              HuMoment5,
              HuMoment6,
              HuMoment7,
              Inertia,
              Compactness,
              Eccentricity,
              LengthMajorAxis,
              LengthMinorAxis,
              AxisOrientation,
              Circularity,
              BitquadArea,
              BitquadNumberEuler,
              BitquadNumberHoles,
              BitquadPerimeter,
              BitquadPerimeterContinuous,
              BitquadCircularity,
              BitquadPerimeterAverage,
              BitquadLengthAverage,
              BitquadWidthAverage,
              MaxDist,
              MaxDistExact,
              ContourPixels,
              ContourPerimeter,
              ContourSideNorth,
              ContourSideWest,
              ContourSideEast,
              ContourSideSouth}},
            {AttributeGroup::Moments,
             {CentralMoment20, CentralMoment02, CentralMoment11, CentralMoment30, CentralMoment03, CentralMoment21, CentralMoment12,
              HuMoment1,       HuMoment2,       HuMoment3,       HuMoment4,       HuMoment5,       HuMoment6,       HuMoment7,
              Inertia,           Compactness,       Eccentricity,      LengthMajorAxis, LengthMinorAxis, AxisOrientation,  Circularity}},
            {AttributeGroup::Boundary,
             {BitquadArea, BitquadNumberEuler, BitquadNumberHoles, BitquadPerimeter, BitquadPerimeterContinuous, BitquadCircularity,
              BitquadPerimeterAverage, BitquadLengthAverage, BitquadWidthAverage, ContourPixels, ContourPerimeter, ContourSideNorth,
              ContourSideWest, ContourSideEast, ContourSideSouth}},
            {AttributeGroup::TreeTopology,
             {SubtreeHeight, DepthNode, IsLeafNode, IsRootNode, NumChildrenNode, NumSiblingsNode, NumDescendantsNode, NumLeafDescendantsNode,
              LeafRatioNode, BalanceNode, AvgChildHeightNode}},
            {AttributeGroup::DistTransf,
             {MaxDist,
              MaxSquaredDist,
              DistSquaredSum,
              DistSquaredMean,
              DistRms,
              DistSquaredVariance,
              MaxDistCenterRow,
              MaxDistCenterColumn,
              MaxDistPlateauArea,
              MaxDistPlateauCentroidRow,
              MaxDistPlateauCentroidColumn,
              DistSum,
              DistMean,
              DistVariance,
              DistMedian,
              DistMode,
              DistQ25,
              DistQ75,
              DistQ90,
              DistEntropy,
              DistPositiveArea,
              DistLevelCount,
              DistWeightedCentroidRow,
              DistWeightedCentroidColumn,
              DistWeightedCentralMoment20,
              DistWeightedCentralMoment02,
              DistWeightedCentralMoment11,
              DistWeightedAxisOrientation,
              DistWeightedEccentricity}},
            {AttributeGroup::DistTransfExact,
             {MaxDistExact,
              MaxSquaredDistExact,
              DistSquaredSumExact,
              DistSquaredMeanExact,
              DistRmsExact,
              DistSquaredVarianceExact,
              MaxDistCenterRowExact,
              MaxDistCenterColumnExact,
              MaxDistPlateauAreaExact,
              MaxDistPlateauCentroidRowExact,
              MaxDistPlateauCentroidColumnExact,
              DistSumExact,
              DistMeanExact,
              DistVarianceExact,
              DistMedianExact,
              DistModeExact,
              DistQ25Exact,
              DistQ75Exact,
              DistQ90Exact,
              DistEntropyExact,
              DistPositiveAreaExact,
              DistLevelCountExact,
              DistWeightedCentroidRowExact,
              DistWeightedCentroidColumnExact,
              DistWeightedCentralMoment20Exact,
              DistWeightedCentralMoment02Exact,
              DistWeightedCentralMoment11Exact,
              DistWeightedAxisOrientationExact,
              DistWeightedEccentricityExact}},
            {AttributeGroup::All, std::move(all)}};
    }();
    return groups;
}

} // namespace mmcfilters::attributes::registry
