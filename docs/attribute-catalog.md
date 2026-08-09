# Attribute catalog

This reference lists the scalar attributes exposed by the public API. Use the
constant directly in C++ requests, for example `AREA`, or through the Python
facade, for example `mmcfilters.Attribute.AREA`.

Every listed attribute produces one value per node in the returned
`AttributeNames` layout. For usage patterns and output-space rules, see
[Attributes](attributes.md).

The `Contract` column classifies the input required by the attribute:

- `Topology/support`: does not read node altitudes. Some attributes use only
  finite proper-part ownership; others additionally declare regular 2D geometry
  or adjacency requirements.
- `Altitude-aware`: reads the node altitude buffer and therefore requires a
  `WeightedMorphologicalTree<T>` or `WeightedTreeView<T>`.
- `Tree topology`: uses only parent/child relations in the hierarchy.

The `Groups` column lists non-`ALL` group memberships. Rows are sorted by the
first listed group; attributes that belong to more than one group appear only
once. `ALL` expands to every attribute in this table.

The complete capability matrix is:

| Attributes | Altitude | Grid 2D | Adjacency | Monotone order | Altitude for directional adjacency | Canonical 4/8 |
| --- | --- | --- | --- | --- | --- | --- |
| `AREA`, `HEIGHT_NODE` through `BALANCE_NODE`, `AVG_CHILD_HEIGHT_NODE` | no | no | none | no | no | no |
| `VOLUME`, `RELATIVE_VOLUME`, `LEVEL`, `GRAY_HEIGHT`, `MEAN_LEVEL`, `VARIANCE_LEVEL` | yes | no | none | no | no | no |
| bounding boxes, central/Hu moments, moment-derived attributes, `CONTOUR_*` | no | yes | none | no | no | no |
| `BITQUADS_*` | no | yes | uniform or directional | no | yes | yes |
| `MAX_DIST` | yes | yes | uniform | yes | no | no |

C++ callers can query
`attributes::registry::capabilityRequirements(attribute)`. Python callers use
`mmcfilters.Attribute.requirements(attribute)`. Group requests are expanded to
scalar attributes before these requirements are validated.

| Constant | Groups | Contract | Description |
| --- | --- | --- | --- |
| `VOLUME` | `GRAY_LEVEL` | Altitude-aware | Sum of altitude-weighted support contributions over the node subtree. It behaves like the gray-level mass or integral of the image over the connected component support. |
| `RELATIVE_VOLUME` | `GRAY_LEVEL` | Altitude-aware | Cumulative gray-level contrast volume: absolute parent/child altitude jumps weighted by child support area, plus the node support-area contribution. |
| `LEVEL` | `GRAY_LEVEL` | Altitude-aware | Altitude of the node in the morphological hierarchy. For component trees, this is the gray level at which the connected component appears. |
| `GRAY_HEIGHT` | `GRAY_LEVEL` | Altitude-aware | Maximum absolute altitude difference between the node and any node in its subtree. On monotone max-tree and min-tree hierarchies this reduces to the traditional one-sided span; it also applies to hierarchies with unconstrained altitude order. Leaves have value `0`. |
| `MEAN_LEVEL` | `GRAY_LEVEL` | Altitude-aware | Average gray level over the full node support. It is computed from accumulated `VOLUME / AREA`. |
| `VARIANCE_LEVEL` | `GRAY_LEVEL` | Altitude-aware | Variance of gray levels over the full node support. It uses the accumulated squared gray-level sum and the mean level. |
| `AREA` | `SHAPE` | Topology/support | Number of proper parts in the full node support, including all descendant supports. In image-domain trees this is the pixel count of the connected component represented by the node. |
| `BOX_WIDTH` | `SHAPE` | Topology/support | Width, in columns, of the smallest axis-aligned bounding box enclosing the node support. |
| `BOX_HEIGHT` | `SHAPE` | Topology/support | Height, in rows, of the smallest axis-aligned bounding box enclosing the node support. |
| `DIAGONAL_LENGTH` | `SHAPE` | Topology/support | Euclidean diagonal length of the bounding box, `sqrt(width^2 + height^2)`. |
| `RECTANGULARITY` | `SHAPE` | Topology/support | Ratio `AREA / (BOX_WIDTH * BOX_HEIGHT)`. Values closer to `1` indicate that the support fills its bounding box densely. |
| `RATIO_WH` | `SHAPE` | Topology/support | Bounding-box aspect ratio, `max(width, height) / min(width, height)` for non-degenerate boxes. Values are at least `1`. |
| `BOX_COL_MIN` | `SHAPE` | Topology/support | Minimum image column index covered by the node support. |
| `BOX_COL_MAX` | `SHAPE` | Topology/support | Maximum image column index covered by the node support. |
| `BOX_ROW_MIN` | `SHAPE` | Topology/support | Minimum image row index covered by the node support. |
| `BOX_ROW_MAX` | `SHAPE` | Topology/support | Maximum image row index covered by the node support. |
| `MAX_DIST` | `SHAPE` | Altitude-aware | Maximum squared Euclidean distance reached from the node contour during the incremental distance-transform sweep. Requires a regular 2D domain, globally monotone altitude order, and uniform adjacency; the descriptive tree kind is irrelevant. |
| `CENTRAL_MOMENT_20` | `MOMENTS`, `SHAPE` | Topology/support | Second-order central moment `mu20` around the support centroid. It measures horizontal spread using column coordinates as `x`. |
| `CENTRAL_MOMENT_02` | `MOMENTS`, `SHAPE` | Topology/support | Second-order central moment `mu02` around the support centroid. It measures vertical spread using row coordinates as `y`. |
| `CENTRAL_MOMENT_11` | `MOMENTS`, `SHAPE` | Topology/support | Mixed second-order central moment `mu11`. It measures covariance between column and row coordinates. |
| `CENTRAL_MOMENT_30` | `MOMENTS`, `SHAPE` | Topology/support | Third-order central moment `mu30`. It captures horizontal asymmetry of the support around its centroid. |
| `CENTRAL_MOMENT_03` | `MOMENTS`, `SHAPE` | Topology/support | Third-order central moment `mu03`. It captures vertical asymmetry of the support around its centroid. |
| `CENTRAL_MOMENT_21` | `MOMENTS`, `SHAPE` | Topology/support | Mixed third-order central moment `mu21`. It captures combined horizontal spread and vertical asymmetry. |
| `CENTRAL_MOMENT_12` | `MOMENTS`, `SHAPE` | Topology/support | Mixed third-order central moment `mu12`. It captures combined vertical spread and horizontal asymmetry. |
| `HU_MOMENT_1` | `MOMENTS`, `SHAPE` | Topology/support | First Hu moment invariant computed from normalized central moments. It follows the standard second-order Hu expression. |
| `HU_MOMENT_2` | `MOMENTS`, `SHAPE` | Topology/support | Second Hu moment invariant. It emphasizes anisotropy between horizontal and vertical spread while remaining invariant to translation, scale, and rotation. |
| `HU_MOMENT_3` | `MOMENTS`, `SHAPE` | Topology/support | Third Hu moment invariant. It captures third-order skewness and asymmetry patterns in the support distribution. |
| `HU_MOMENT_4` | `MOMENTS`, `SHAPE` | Topology/support | Fourth Hu moment invariant. It combines third-order normalized moments to describe diagonal and off-axis symmetry. |
| `HU_MOMENT_5` | `MOMENTS`, `SHAPE` | Topology/support | Fifth Hu moment invariant. It is sensitive to more complex orientation-dependent asymmetries and reflection-related differences. |
| `HU_MOMENT_6` | `MOMENTS`, `SHAPE` | Topology/support | Sixth Hu moment invariant. It combines second- and third-order normalized moments to capture elliptical asymmetry and curvature-related shape variation. |
| `HU_MOMENT_7` | `MOMENTS`, `SHAPE` | Topology/support | Seventh Hu moment invariant. It is highly sensitive to fine asymmetries and helps distinguish mirror-related shapes. |
| `INERTIA` | `MOMENTS`, `SHAPE` | Topology/support | Sum of normalized second-order central moments, `mu20 / area^2 + mu02 / area^2`. This is the same scalar expression as the first Hu moment invariant. |
| `COMPACTNESS` | `MOMENTS`, `SHAPE` | Topology/support | Area normalized by second-order dispersion, `(1 / (2*pi)) * area / (mu20 + mu02)` when the denominator is positive. Higher values indicate more compact supports. |
| `ECCENTRICITY` | `MOMENTS`, `SHAPE` | Topology/support | Ratio of the largest to smallest eigenvalue of the second-moment matrix. Values near `1` indicate isotropic shapes; larger values indicate elongation. Degenerate line-like supports saturate at the finite value `1e6`. |
| `LENGTH_MAJOR_AXIS` | `MOMENTS`, `SHAPE` | Topology/support | Length proxy for the major axis of the equivalent second-moment ellipse, derived from the largest inertia eigenvalue and area. |
| `LENGTH_MINOR_AXIS` | `MOMENTS`, `SHAPE` | Topology/support | Length proxy for the minor axis of the equivalent second-moment ellipse, derived from the smallest inertia eigenvalue and area. |
| `AXIS_ORIENTATION` | `MOMENTS`, `SHAPE` | Topology/support | Principal-axis orientation in degrees, `0.5 * atan2(2*mu11, mu20 - mu02)`, normalized to a non-negative angle. |
| `CIRCULARITY` | `MOMENTS`, `SHAPE` | Topology/support | Ratio `lambda2 / lambda1` of second-moment eigenvalues. Values near `1` indicate circular or isotropic supports; values near `0` indicate elongation. |
| `BITQUADS_AREA` | `BOUNDARY`, `SHAPE` | Topology/support | Duda-style sub-pixel area estimator derived from aggregated `2x2` bitquad pattern counts. |
| `BITQUADS_NUMBER_EULER` | `BOUNDARY`, `SHAPE` | Topology/support | Euler characteristic estimated from bitquad counters, representing connected components minus holes under the selected connectivity projection. |
| `BITQUADS_NUMBER_HOLES` | `BOUNDARY`, `SHAPE` | Topology/support | Number of holes inferred from the bitquad Euler characteristic for a single connected support. |
| `BITQUADS_PERIMETER` | `BOUNDARY`, `SHAPE` | Topology/support | Discrete boundary-length estimate from bitquad edge-contributing patterns. |
| `BITQUADS_PERIMETER_CONTINUOUS` | `BOUNDARY`, `SHAPE` | Topology/support | Smoothed continuous perimeter estimate from bitquad counters, using weighted transitions across local `2x2` configurations. |
| `BITQUADS_CIRCULARITY` | `BOUNDARY`, `SHAPE` | Topology/support | Bitquad compactness measure `(4*pi*BITQUADS_AREA) / BITQUADS_PERIMETER_CONTINUOUS^2`. Values closer to `1` indicate rounder supports. Degenerate zero-perimeter supports return `0`. |
| `BITQUADS_PERIMETER_AVERAGE` | `BOUNDARY`, `SHAPE` | Topology/support | Average continuous perimeter per connected component, computed from bitquad perimeter and Euler-count estimates. Euler estimates `<= 0` return `0`. |
| `BITQUADS_LENGTH_AVERAGE` | `BOUNDARY`, `SHAPE` | Topology/support | Average longitudinal extent proxy, derived as half of the average continuous perimeter. Euler estimates `<= 0` return `0`. |
| `BITQUADS_WIDTH_AVERAGE` | `BOUNDARY`, `SHAPE` | Topology/support | Average transverse extent proxy, computed as `2 * BITQUADS_AREA / BITQUADS_PERIMETER_CONTINUOUS` with a zero fallback when the continuous perimeter is degenerate. |
| `CONTOUR_PIXELS` | `BOUNDARY`, `SHAPE` | Topology/support | Number of support pixels that touch the 4-neighbor complement by at least one side. |
| `CONTOUR_PERIMETER` | `BOUNDARY`, `SHAPE` | Topology/support | Total number of exposed 4-neighbor sides over the support. This is the side-count perimeter, not a Euclidean perimeter estimate. |
| `CONTOUR_SIDE_NORTH` | `BOUNDARY`, `SHAPE` | Topology/support | Number of exposed north-facing sides over support pixels. |
| `CONTOUR_SIDE_WEST` | `BOUNDARY`, `SHAPE` | Topology/support | Number of exposed west-facing sides over support pixels. |
| `CONTOUR_SIDE_EAST` | `BOUNDARY`, `SHAPE` | Topology/support | Number of exposed east-facing sides over support pixels. |
| `CONTOUR_SIDE_SOUTH` | `BOUNDARY`, `SHAPE` | Topology/support | Number of exposed south-facing sides over support pixels. |
| `HEIGHT_NODE` | `TREE_TOPOLOGY` | Tree topology | Longest child-edge path from the node to any leaf in its subtree. Leaves have height `0`. |
| `DEPTH_NODE` | `TREE_TOPOLOGY` | Tree topology | Number of parent-edge steps from the node to the root. The root has depth `0`. |
| `IS_LEAF_NODE` | `TREE_TOPOLOGY` | Tree topology | Boolean scalar encoded as `1` when the node has no children and `0` otherwise. |
| `IS_ROOT_NODE` | `TREE_TOPOLOGY` | Tree topology | Boolean scalar encoded as `1` for the tree root and `0` for all other nodes. |
| `NUM_CHILDREN_NODE` | `TREE_TOPOLOGY` | Tree topology | Number of direct child nodes. It measures the immediate branching factor of the hierarchy at the node. |
| `NUM_SIBLINGS_NODE` | `TREE_TOPOLOGY` | Tree topology | Number of other nodes sharing the same parent. The root has `0` siblings. |
| `NUM_DESCENDANTS_NODE` | `TREE_TOPOLOGY` | Tree topology | Number of internal tree nodes strictly below this node in its subtree. |
| `NUM_LEAF_DESCENDANTS_NODE` | `TREE_TOPOLOGY` | Tree topology | Number of leaf nodes in the node subtree. A leaf node contributes `1` for itself. |
| `LEAF_RATIO_NODE` | `TREE_TOPOLOGY` | Tree topology | Ratio of leaf descendants to subtree size, `leaf_descendants / (descendants + 1)`. Leaves return `1`. |
| `BALANCE_NODE` | `TREE_TOPOLOGY` | Tree topology | Difference between maximum and minimum child-subtree heights. It is `0` for leaves and increases when child depths are uneven. |
| `AVG_CHILD_HEIGHT_NODE` | `TREE_TOPOLOGY` | Tree topology | Average height of direct child subtrees. Leaves return `0`; non-leaf nodes average direct child heights. |
