# Attribute catalog

This reference lists the scalar attributes exposed by the public API. C++ uses
the PascalCase enumerators from `Attribute`, for example `Attribute::Area`.
Stable symbolic names use uppercase snake case in string requests, returned
layouts, and Python, for example `"AREA"` and `mmcfilters.Attribute.AREA`.

Every listed attribute produces one value per node in the returned
`AttributeNames` layout. For usage patterns and output-space rules, see
[Attributes](attributes.md). Family-specific mathematical conventions are kept
in focused guides such as
[Distance-transform attributes](distance-transform.md).

The `Contract` column classifies the input required by the attribute:

- `Topology/support`: does not read node altitudes. Some attributes use only
  finite smallest-node mapping; others additionally declare regular 2D geometry
  or adjacency requirements.
- `Altitude-aware`: reads the node altitude buffer and therefore requires a
  `ValuedMorphologicalTree<T>` or `ValuedMorphologicalTreeView<T>`.
- `Tree topology`: uses only parent/child relations in the hierarchy.

The `Groups` column lists non-`ALL` group memberships. Rows are sorted by the
first listed group; attributes that belong to more than one group appear only
once. `ALL` expands to every attribute in this table.

The complete capability matrix is:

| Attributes | Altitude | Grid 2D | Adjacency | Monotone order | Altitude for directional adjacency | Canonical 4/8 |
| --- | --- | --- | --- | --- | --- | --- |
| `AREA`, `SUBTREE_HEIGHT` through `BALANCE_NODE`, `AVG_CHILD_HEIGHT_NODE` | no | no | none | no | no | no |
| `VOLUME`, `RELATIVE_VOLUME`, `GRAY_LEVEL_HEIGHT`, `MEAN_GRAY_LEVEL`, `GRAY_LEVEL_VARIANCE` | yes | no | none | no | no | no |
| bounding boxes, central/Hu moments, moment-derived attributes, `CONTOUR_*` | no | yes | none | no | no | no |
| `BITQUAD_*` | no | yes | uniform or directional | no | yes | yes |
| `MAX_DIST_EXACT`, `MAX_DIST`, `MAX_SQUARED_DIST*`, `MAX_DIST_CENTER_*`, `DIST_*` | no | yes | none | no | no | no |

C++ callers can query
`attributes::registry::capabilityRequirements(attribute)`. Python callers use
`mmcfilters.Attribute.requirements(attribute)`. Group requests are expanded to
scalar attributes before these requirements are validated.

| C++ enumerator | Stable/Python name | Groups | Contract | Description |
| --- | --- | --- | --- | --- |
| `Attribute::Volume` | `VOLUME` | `GRAY_LEVEL` | Altitude-aware | Sum of altitude-weighted support contributions over the node subtree. It behaves like the gray-level mass or integral of the image over the connected component support. |
| `Attribute::RelativeVolume` | `RELATIVE_VOLUME` | `GRAY_LEVEL` | Altitude-aware | Recursive contrast volume `R(n) = area(n) + sum_c [R(c) + area(c) * abs(altitude(c) - altitude(n))]` over the direct children `c`. |
| `Attribute::GrayLevelHeight` | `GRAY_LEVEL_HEIGHT` | `GRAY_LEVEL` | Altitude-aware | Maximum absolute altitude difference between the node and any node in its subtree. On monotone max-tree and min-tree hierarchies this reduces to the traditional one-sided span; it also applies to hierarchies with unconstrained altitude order. Leaves have value `0`. |
| `Attribute::MeanGrayLevel` | `MEAN_GRAY_LEVEL` | `GRAY_LEVEL` | Altitude-aware | Arithmetic mean of the image values over the full node support: `sum(f(x), x in X) / card(X)`. |
| `Attribute::GrayLevelVariance` | `GRAY_LEVEL_VARIANCE` | `GRAY_LEVEL` | Altitude-aware | Population variance of the image values over the full node support, with denominator `card(X)`. |
| `Attribute::Area` | `AREA` | `SHAPE` | Topology/support | Number of pixels in the full node support. Equivalently, it is the sum of proper-part cardinalities over the node's subtree. |
| `Attribute::BoxWidth` | `BOX_WIDTH` | `SHAPE` | Topology/support | Width, in columns, of the smallest axis-aligned bounding box enclosing the node support. |
| `Attribute::BoundingBoxHeight` | `BOUNDING_BOX_HEIGHT` | `SHAPE` | Topology/support | Height, in rows, of the smallest axis-aligned bounding box enclosing the node support. |
| `Attribute::DiagonalLength` | `DIAGONAL_LENGTH` | `SHAPE` | Topology/support | Euclidean diagonal length of the bounding box, `sqrt(width^2 + height^2)`. |
| `Attribute::Rectangularity` | `RECTANGULARITY` | `SHAPE` | Topology/support | Ratio `AREA / (BOX_WIDTH * BOUNDING_BOX_HEIGHT)`. Values closer to `1` indicate that the support fills its bounding box densely. |
| `Attribute::RatioWh` | `RATIO_WH` | `SHAPE` | Topology/support | Bounding-box aspect ratio, `max(width, height) / min(width, height)` for non-degenerate boxes. Values are at least `1`. |
| `Attribute::BoxColumnMin` | `BOX_COLUMN_MIN` | `SHAPE` | Topology/support | Minimum image column index covered by the node support. |
| `Attribute::BoxColumnMax` | `BOX_COLUMN_MAX` | `SHAPE` | Topology/support | Maximum image column index covered by the node support. |
| `Attribute::BoxRowMin` | `BOX_ROW_MIN` | `SHAPE` | Topology/support | Minimum image row index covered by the node support. |
| `Attribute::BoxRowMax` | `BOX_ROW_MAX` | `SHAPE` | Topology/support | Maximum image row index covered by the node support. |
| `Attribute::MaxDistExact` | `MAX_DIST_EXACT` | `SHAPE`, `DIST_TRANSF_EXACT` | Topology/support | Maximum exact Euclidean distance, in pixels, from a support pixel to the foreground 4-connected contour. It is the square root of `MAX_SQUARED_DIST_EXACT`. |
| `Attribute::MaxDist` | `MAX_DIST` | `SHAPE`, `DIST_TRANSF` | Topology/support | Maximum of the approximate distance field, in pixels. It is the square root of `MAX_SQUARED_DIST` and may differ from `MAX_DIST_EXACT`. |
| `Attribute::MaxSquaredDistExact` | `MAX_SQUARED_DIST_EXACT` | `DIST_TRANSF_EXACT` | Topology/support | Maximum exact squared Euclidean contour distance, in squared pixels. This exposes the maximum integer cost of the exact field without applying a square root. |
| `Attribute::MaxSquaredDist` | `MAX_SQUARED_DIST` | `DIST_TRANSF` | Topology/support | Maximum approximate squared contour-distance cost, in squared pixels. |
| `Attribute::DistSquaredSumExact` | `DIST_SQUARED_SUM_EXACT` | `DIST_TRANSF_EXACT` | Topology/support | Sum of the exact squared Euclidean distance to the foreground A4 contour over every support pixel. |
| `Attribute::DistSquaredMeanExact` | `DIST_SQUARED_MEAN_EXACT` | `DIST_TRANSF_EXACT` | Topology/support | Arithmetic mean of the exact squared contour-distance field over the node support. |
| `Attribute::DistRmsExact` | `DIST_RMS_EXACT` | `DIST_TRANSF_EXACT` | Topology/support | Root mean square Euclidean contour distance, `sqrt(mean(d^2))`, in pixel units. This is not the arithmetic mean of `d`. |
| `Attribute::DistSquaredVarianceExact` | `DIST_SQUARED_VARIANCE_EXACT` | `DIST_TRANSF_EXACT` | Topology/support | Population variance of the exact squared contour-distance samples, `mean((d^2)^2) - mean(d^2)^2`, in fourth-power pixel units. |
| `Attribute::DistSquaredSum` | `DIST_SQUARED_SUM` | `DIST_TRANSF` | Topology/support | Sum of the approximate squared-distance field over the node support. |
| `Attribute::DistSquaredMean` | `DIST_SQUARED_MEAN` | `DIST_TRANSF` | Topology/support | Arithmetic mean of the approximate squared-distance field over the node support. It need not equal `DIST_SQUARED_MEAN_EXACT`. |
| `Attribute::DistRms` | `DIST_RMS` | `DIST_TRANSF` | Topology/support | Root mean square distance derived from the approximate field, `sqrt(mean(cost))`, in pixel units. |
| `Attribute::DistSquaredVariance` | `DIST_SQUARED_VARIANCE` | `DIST_TRANSF` | Topology/support | Population variance of the approximate squared-distance samples. |
| `Attribute::MaxDistCenterRowExact` | `MAX_DIST_CENTER_ROW_EXACT` | `DIST_TRANSF_EXACT` | Topology/support | Zero-based row of the exact support pixel attaining `MAX_DIST_EXACT`. If several pixels attain the maximum, the smallest row-major pixel is selected. |
| `Attribute::MaxDistCenterColumnExact` | `MAX_DIST_CENTER_COLUMN_EXACT` | `DIST_TRANSF_EXACT` | Topology/support | Zero-based column of the same deterministic exact maximum-distance center. Request row and column together to recover its position. |
| `Attribute::MaxDistCenterRow` | `MAX_DIST_CENTER_ROW` | `DIST_TRANSF` | Topology/support | Zero-based row of the support pixel attaining the approximate maximum. The same smallest-row-major tie rule is applied. |
| `Attribute::MaxDistCenterColumn` | `MAX_DIST_CENTER_COLUMN` | `DIST_TRANSF` | Topology/support | Zero-based column of the same deterministic approximate maximum-distance center. Approximate and exact centers may differ even when their maximum values agree. |
| `Attribute::MaxDistPlateauAreaExact` | `MAX_DIST_PLATEAU_AREA_EXACT` | `DIST_TRANSF_EXACT` | Topology/support | Number of support pixels attaining the exact maximum squared contour distance. |
| `Attribute::MaxDistPlateauCentroidRowExact` | `MAX_DIST_PLATEAU_CENTROID_ROW_EXACT` | `DIST_TRANSF_EXACT` | Topology/support | Arithmetic mean of the row coordinates over the exact maximum-distance plateau. It may be fractional. |
| `Attribute::MaxDistPlateauCentroidColumnExact` | `MAX_DIST_PLATEAU_CENTROID_COLUMN_EXACT` | `DIST_TRANSF_EXACT` | Topology/support | Arithmetic mean of the column coordinates over the exact maximum-distance plateau. It may be fractional. |
| `Attribute::MaxDistPlateauArea` | `MAX_DIST_PLATEAU_AREA` | `DIST_TRANSF` | Topology/support | Number of support pixels attaining the approximate maximum. It may differ from the exact plateau area. |
| `Attribute::MaxDistPlateauCentroidRow` | `MAX_DIST_PLATEAU_CENTROID_ROW` | `DIST_TRANSF` | Topology/support | Arithmetic mean of row coordinates over the approximate maximum plateau. |
| `Attribute::MaxDistPlateauCentroidColumn` | `MAX_DIST_PLATEAU_CENTROID_COLUMN` | `DIST_TRANSF` | Topology/support | Arithmetic mean of column coordinates over the approximate maximum plateau. |
| `Attribute::DistSum` | `DIST_SUM` | `DIST_TRANSF` | Topology/support | Sum of Euclidean distances `sqrt(z)` over the approximate distance field, in accumulated pixel units. |
| `Attribute::DistMean` | `DIST_MEAN` | `DIST_TRANSF` | Topology/support | Arithmetic mean of Euclidean distances over the approximate field, in pixels. |
| `Attribute::DistVariance` | `DIST_VARIANCE` | `DIST_TRANSF` | Topology/support | Population variance of Euclidean distances, `mean(z) - mean(sqrt(z))^2`, in squared pixels. |
| `Attribute::DistMedian` | `DIST_MEDIAN` | `DIST_TRANSF` | Topology/support | Lower empirical median of the approximate Euclidean-distance field. |
| `Attribute::DistMode` | `DIST_MODE` | `DIST_TRANSF` | Topology/support | Most frequent approximate Euclidean distance; ties select the smallest distance. |
| `Attribute::DistQ25` | `DIST_Q25` | `DIST_TRANSF` | Topology/support | Lower empirical 25th percentile of the approximate Euclidean-distance field. |
| `Attribute::DistQ75` | `DIST_Q75` | `DIST_TRANSF` | Topology/support | Lower empirical 75th percentile of the approximate Euclidean-distance field. |
| `Attribute::DistQ90` | `DIST_Q90` | `DIST_TRANSF` | Topology/support | Lower empirical 90th percentile of the approximate Euclidean-distance field. |
| `Attribute::DistEntropy` | `DIST_ENTROPY` | `DIST_TRANSF` | Topology/support | Shannon entropy in bits of the normalized approximate squared-distance histogram. |
| `Attribute::DistPositiveArea` | `DIST_POSITIVE_AREA` | `DIST_TRANSF` | Topology/support | Number of support pixels whose approximate squared-distance cost is positive. |
| `Attribute::DistLevelCount` | `DIST_LEVEL_COUNT` | `DIST_TRANSF` | Topology/support | Number of represented approximate squared-distance levels, including zero. |
| `Attribute::DistWeightedCentroidRow` | `DIST_WEIGHTED_CENTROID_ROW` | `DIST_TRANSF` | Topology/support | Zero-based row centroid weighted by approximate Euclidean distance. All-zero fields use the ordinary support centroid. |
| `Attribute::DistWeightedCentroidColumn` | `DIST_WEIGHTED_CENTROID_COLUMN` | `DIST_TRANSF` | Topology/support | Zero-based column centroid weighted by approximate Euclidean distance. All-zero fields use the ordinary support centroid. |
| `Attribute::DistWeightedCentralMoment20` | `DIST_WEIGHTED_CENTRAL_MOMENT_20` | `DIST_TRANSF` | Topology/support | Unnormalized approximate distance-weighted central column moment. |
| `Attribute::DistWeightedCentralMoment02` | `DIST_WEIGHTED_CENTRAL_MOMENT_02` | `DIST_TRANSF` | Topology/support | Unnormalized approximate distance-weighted central row moment. |
| `Attribute::DistWeightedCentralMoment11` | `DIST_WEIGHTED_CENTRAL_MOMENT_11` | `DIST_TRANSF` | Topology/support | Unnormalized approximate distance-weighted mixed row-column moment. |
| `Attribute::DistWeightedAxisOrientation` | `DIST_WEIGHTED_AXIS_ORIENTATION` | `DIST_TRANSF` | Topology/support | Principal axis of the approximate distance-weighted second moments, in degrees; isotropic fields return zero. |
| `Attribute::DistWeightedEccentricity` | `DIST_WEIGHTED_ECCENTRICITY` | `DIST_TRANSF` | Topology/support | Major/minor eigenvalue ratio of the approximate distance-weighted second moments, capped at `1e6`. |
| `Attribute::DistSumExact` | `DIST_SUM_EXACT` | `DIST_TRANSF_EXACT` | Topology/support | Sum of exact Euclidean contour distances over the support, in accumulated pixel units. |
| `Attribute::DistMeanExact` | `DIST_MEAN_EXACT` | `DIST_TRANSF_EXACT` | Topology/support | Arithmetic mean of exact Euclidean contour distances, in pixels. |
| `Attribute::DistVarianceExact` | `DIST_VARIANCE_EXACT` | `DIST_TRANSF_EXACT` | Topology/support | Population variance of exact Euclidean contour distances, in squared pixels. |
| `Attribute::DistMedianExact` | `DIST_MEDIAN_EXACT` | `DIST_TRANSF_EXACT` | Topology/support | Lower empirical median of the exact Euclidean-distance field. |
| `Attribute::DistModeExact` | `DIST_MODE_EXACT` | `DIST_TRANSF_EXACT` | Topology/support | Most frequent exact Euclidean distance; ties select the smallest distance. |
| `Attribute::DistQ25Exact` | `DIST_Q25_EXACT` | `DIST_TRANSF_EXACT` | Topology/support | Lower empirical 25th percentile of the exact Euclidean-distance field. |
| `Attribute::DistQ75Exact` | `DIST_Q75_EXACT` | `DIST_TRANSF_EXACT` | Topology/support | Lower empirical 75th percentile of the exact Euclidean-distance field. |
| `Attribute::DistQ90Exact` | `DIST_Q90_EXACT` | `DIST_TRANSF_EXACT` | Topology/support | Lower empirical 90th percentile of the exact Euclidean-distance field. |
| `Attribute::DistEntropyExact` | `DIST_ENTROPY_EXACT` | `DIST_TRANSF_EXACT` | Topology/support | Shannon entropy in bits of the normalized exact squared-distance histogram. |
| `Attribute::DistPositiveAreaExact` | `DIST_POSITIVE_AREA_EXACT` | `DIST_TRANSF_EXACT` | Topology/support | Number of support pixels whose exact squared contour distance is positive. |
| `Attribute::DistLevelCountExact` | `DIST_LEVEL_COUNT_EXACT` | `DIST_TRANSF_EXACT` | Topology/support | Number of represented exact squared-distance levels, including zero. |
| `Attribute::DistWeightedCentroidRowExact` | `DIST_WEIGHTED_CENTROID_ROW_EXACT` | `DIST_TRANSF_EXACT` | Topology/support | Zero-based row centroid weighted by exact Euclidean distance. All-zero fields use the ordinary support centroid. |
| `Attribute::DistWeightedCentroidColumnExact` | `DIST_WEIGHTED_CENTROID_COLUMN_EXACT` | `DIST_TRANSF_EXACT` | Topology/support | Zero-based column centroid weighted by exact Euclidean distance. All-zero fields use the ordinary support centroid. |
| `Attribute::DistWeightedCentralMoment20Exact` | `DIST_WEIGHTED_CENTRAL_MOMENT_20_EXACT` | `DIST_TRANSF_EXACT` | Topology/support | Unnormalized exact distance-weighted central column moment. |
| `Attribute::DistWeightedCentralMoment02Exact` | `DIST_WEIGHTED_CENTRAL_MOMENT_02_EXACT` | `DIST_TRANSF_EXACT` | Topology/support | Unnormalized exact distance-weighted central row moment. |
| `Attribute::DistWeightedCentralMoment11Exact` | `DIST_WEIGHTED_CENTRAL_MOMENT_11_EXACT` | `DIST_TRANSF_EXACT` | Topology/support | Unnormalized exact distance-weighted mixed row-column moment. |
| `Attribute::DistWeightedAxisOrientationExact` | `DIST_WEIGHTED_AXIS_ORIENTATION_EXACT` | `DIST_TRANSF_EXACT` | Topology/support | Principal axis of the exact distance-weighted second moments, in degrees; isotropic fields return zero. |
| `Attribute::DistWeightedEccentricityExact` | `DIST_WEIGHTED_ECCENTRICITY_EXACT` | `DIST_TRANSF_EXACT` | Topology/support | Major/minor eigenvalue ratio of the exact distance-weighted second moments, capped at `1e6`. |
| `Attribute::CentralMoment20` | `CENTRAL_MOMENT_20` | `MOMENTS`, `SHAPE` | Topology/support | Second-order central moment `mu20` around the support centroid. It measures horizontal spread using column coordinates as `x`. |
| `Attribute::CentralMoment02` | `CENTRAL_MOMENT_02` | `MOMENTS`, `SHAPE` | Topology/support | Second-order central moment `mu02` around the support centroid. It measures vertical spread using row coordinates as `y`. |
| `Attribute::CentralMoment11` | `CENTRAL_MOMENT_11` | `MOMENTS`, `SHAPE` | Topology/support | Mixed second-order central moment `mu11`. It measures covariance between column and row coordinates. |
| `Attribute::CentralMoment30` | `CENTRAL_MOMENT_30` | `MOMENTS`, `SHAPE` | Topology/support | Third-order central moment `mu30`. It captures horizontal asymmetry of the support around its centroid. |
| `Attribute::CentralMoment03` | `CENTRAL_MOMENT_03` | `MOMENTS`, `SHAPE` | Topology/support | Third-order central moment `mu03`. It captures vertical asymmetry of the support around its centroid. |
| `Attribute::CentralMoment21` | `CENTRAL_MOMENT_21` | `MOMENTS`, `SHAPE` | Topology/support | Mixed third-order central moment `mu21`. It captures combined horizontal spread and vertical asymmetry. |
| `Attribute::CentralMoment12` | `CENTRAL_MOMENT_12` | `MOMENTS`, `SHAPE` | Topology/support | Mixed third-order central moment `mu12`. It captures combined vertical spread and horizontal asymmetry. |
| `Attribute::HuMoment1` | `HU_MOMENT_1` | `MOMENTS`, `SHAPE` | Topology/support | First Hu moment invariant computed from normalized central moments. It follows the standard second-order Hu expression. |
| `Attribute::HuMoment2` | `HU_MOMENT_2` | `MOMENTS`, `SHAPE` | Topology/support | Second Hu moment invariant. It emphasizes anisotropy between horizontal and vertical spread while remaining invariant to translation, scale, and rotation. |
| `Attribute::HuMoment3` | `HU_MOMENT_3` | `MOMENTS`, `SHAPE` | Topology/support | Third Hu moment invariant. It captures third-order skewness and asymmetry patterns in the support distribution. |
| `Attribute::HuMoment4` | `HU_MOMENT_4` | `MOMENTS`, `SHAPE` | Topology/support | Fourth Hu moment invariant. It combines third-order normalized moments to describe diagonal and off-axis symmetry. |
| `Attribute::HuMoment5` | `HU_MOMENT_5` | `MOMENTS`, `SHAPE` | Topology/support | Fifth Hu moment invariant. It is sensitive to more complex orientation-dependent asymmetries and reflection-related differences. |
| `Attribute::HuMoment6` | `HU_MOMENT_6` | `MOMENTS`, `SHAPE` | Topology/support | Sixth Hu moment invariant. It combines second- and third-order normalized moments to capture elliptical asymmetry and curvature-related shape variation. |
| `Attribute::HuMoment7` | `HU_MOMENT_7` | `MOMENTS`, `SHAPE` | Topology/support | Seventh Hu moment invariant. It is highly sensitive to fine asymmetries and helps distinguish mirror-related shapes. |
| `Attribute::Inertia` | `INERTIA` | `MOMENTS`, `SHAPE` | Topology/support | Sum of normalized second-order central moments, `mu20 / area^2 + mu02 / area^2`. This is the same scalar expression as the first Hu moment invariant. |
| `Attribute::Compactness` | `COMPACTNESS` | `MOMENTS`, `SHAPE` | Topology/support | Area normalized by second-order dispersion, `(1 / (2*pi)) * area / (mu20 + mu02)` when the denominator is positive. Higher values indicate more compact supports. |
| `Attribute::Eccentricity` | `ECCENTRICITY` | `MOMENTS`, `SHAPE` | Topology/support | Ratio of the largest to smallest eigenvalue of the second-moment matrix. Values near `1` indicate isotropic shapes; larger values indicate elongation. Degenerate line-like supports saturate at the finite value `1e6`. |
| `Attribute::LengthMajorAxis` | `LENGTH_MAJOR_AXIS` | `MOMENTS`, `SHAPE` | Topology/support | Length proxy for the major axis of the equivalent second-moment ellipse, derived from the largest inertia eigenvalue and area. |
| `Attribute::LengthMinorAxis` | `LENGTH_MINOR_AXIS` | `MOMENTS`, `SHAPE` | Topology/support | Length proxy for the minor axis of the equivalent second-moment ellipse, derived from the smallest inertia eigenvalue and area. |
| `Attribute::AxisOrientation` | `AXIS_ORIENTATION` | `MOMENTS`, `SHAPE` | Topology/support | Principal-axis orientation in degrees, `0.5 * atan2(2*mu11, mu20 - mu02)`, normalized to a non-negative angle. |
| `Attribute::Circularity` | `CIRCULARITY` | `MOMENTS`, `SHAPE` | Topology/support | Ratio `lambda2 / lambda1` of second-moment eigenvalues. Values near `1` indicate circular or isotropic supports; values near `0` indicate elongation. |
| `Attribute::BitquadArea` | `BITQUAD_AREA` | `BOUNDARY`, `SHAPE` | Topology/support | Duda-style sub-pixel area estimator derived from aggregated `2x2` bitquad pattern counts. |
| `Attribute::BitquadNumberEuler` | `BITQUAD_NUMBER_EULER` | `BOUNDARY`, `SHAPE` | Topology/support | Euler characteristic estimated from bitquad family counts under an explicit connectivity projection. For a tree of shapes, non-root connectivity is selected from exact lower/upper shape polarity; the root has no polarity. |
| `Attribute::BitquadNumberHoles` | `BITQUAD_NUMBER_HOLES` | `BOUNDARY`, `SHAPE` | Topology/support | Number of holes inferred from the bitquad Euler characteristic for a single connected support. |
| `Attribute::BitquadPerimeter` | `BITQUAD_PERIMETER` | `BOUNDARY`, `SHAPE` | Topology/support | Discrete boundary-length estimate from bitquad edge-contributing patterns. |
| `Attribute::BitquadPerimeterContinuous` | `BITQUAD_PERIMETER_CONTINUOUS` | `BOUNDARY`, `SHAPE` | Topology/support | Smoothed continuous perimeter estimate from bitquad counters, using weighted transitions across local `2x2` configurations. |
| `Attribute::BitquadCircularity` | `BITQUAD_CIRCULARITY` | `BOUNDARY`, `SHAPE` | Topology/support | Bitquad compactness measure `(4*pi*BITQUAD_AREA) / BITQUAD_PERIMETER_CONTINUOUS^2`. Values closer to `1` indicate rounder supports. Degenerate zero-perimeter supports return `0`. |
| `Attribute::BitquadPerimeterAverage` | `BITQUAD_PERIMETER_AVERAGE` | `BOUNDARY`, `SHAPE` | Topology/support | Average continuous perimeter per connected component, computed from bitquad perimeter and Euler-count estimates. Euler estimates `<= 0` return `0`. |
| `Attribute::BitquadLengthAverage` | `BITQUAD_LENGTH_AVERAGE` | `BOUNDARY`, `SHAPE` | Topology/support | Average longitudinal extent proxy, derived as half of the average continuous perimeter. Euler estimates `<= 0` return `0`. |
| `Attribute::BitquadWidthAverage` | `BITQUAD_WIDTH_AVERAGE` | `BOUNDARY`, `SHAPE` | Topology/support | Average transverse extent proxy, computed as `2 * BITQUAD_AREA / BITQUAD_PERIMETER_CONTINUOUS` with a zero fallback when the continuous perimeter is degenerate. |
| `Attribute::ContourPixels` | `CONTOUR_PIXELS` | `BOUNDARY`, `SHAPE` | Topology/support | Number of support pixels that touch the 4-neighbor complement by at least one side. |
| `Attribute::ContourPerimeter` | `CONTOUR_PERIMETER` | `BOUNDARY`, `SHAPE` | Topology/support | Total number of exposed 4-neighbor sides over the support. This is the side-count perimeter, not a Euclidean perimeter estimate. |
| `Attribute::ContourSideNorth` | `CONTOUR_SIDE_NORTH` | `BOUNDARY`, `SHAPE` | Topology/support | Number of exposed north-facing sides over support pixels. |
| `Attribute::ContourSideWest` | `CONTOUR_SIDE_WEST` | `BOUNDARY`, `SHAPE` | Topology/support | Number of exposed west-facing sides over support pixels. |
| `Attribute::ContourSideEast` | `CONTOUR_SIDE_EAST` | `BOUNDARY`, `SHAPE` | Topology/support | Number of exposed east-facing sides over support pixels. |
| `Attribute::ContourSideSouth` | `CONTOUR_SIDE_SOUTH` | `BOUNDARY`, `SHAPE` | Topology/support | Number of exposed south-facing sides over support pixels. |
| `Attribute::SubtreeHeight` | `SUBTREE_HEIGHT` | `TREE_TOPOLOGY` | Tree topology | Longest child-edge path from the node to any leaf in its subtree. Leaves have height `0`. |
| `Attribute::DepthNode` | `DEPTH_NODE` | `TREE_TOPOLOGY` | Tree topology | Number of parent-edge steps from the node to the root. The root has depth `0`. |
| `Attribute::IsLeafNode` | `IS_LEAF_NODE` | `TREE_TOPOLOGY` | Tree topology | Boolean scalar encoded as `1` when the node has no children and `0` otherwise. |
| `Attribute::IsRootNode` | `IS_ROOT_NODE` | `TREE_TOPOLOGY` | Tree topology | Boolean scalar encoded as `1` for the tree root and `0` for all other nodes. |
| `Attribute::NumChildrenNode` | `NUM_CHILDREN_NODE` | `TREE_TOPOLOGY` | Tree topology | Number of direct child nodes. It measures the immediate branching factor of the hierarchy at the node. |
| `Attribute::NumSiblingsNode` | `NUM_SIBLINGS_NODE` | `TREE_TOPOLOGY` | Tree topology | Number of other nodes sharing the same parent. The root has `0` siblings. |
| `Attribute::NumDescendantsNode` | `NUM_DESCENDANTS_NODE` | `TREE_TOPOLOGY` | Tree topology | Number of internal tree nodes strictly below this node in its subtree. |
| `Attribute::NumLeafDescendantsNode` | `NUM_LEAF_DESCENDANTS_NODE` | `TREE_TOPOLOGY` | Tree topology | Number of leaf nodes in the node subtree. A leaf node contributes `1` for itself. |
| `Attribute::LeafRatioNode` | `LEAF_RATIO_NODE` | `TREE_TOPOLOGY` | Tree topology | Ratio of leaf descendants to subtree size, `leaf_descendants / (descendants + 1)`. Leaves return `1`. |
| `Attribute::BalanceNode` | `BALANCE_NODE` | `TREE_TOPOLOGY` | Tree topology | Difference between maximum and minimum child-subtree heights. It is `0` for leaves and increases when child depths are uneven. |
| `Attribute::AvgChildHeightNode` | `AVG_CHILD_HEIGHT_NODE` | `TREE_TOPOLOGY` | Tree topology | Average height of direct child subtrees. Leaves return `0`; non-leaf nodes average direct child heights. |
