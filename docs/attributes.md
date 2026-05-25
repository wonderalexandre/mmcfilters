# Attributes

This document maps the `mmcfilters/attributes` subsystem, its execution
contracts, and the scalar attributes currently exposed by the public API.

## Purpose

The attribute layer turns public attribute requests into dense per-node buffers.
It handles:

- scalar attributes and attribute groups;
- topology-only and altitude-aware computation;
- shared intermediate quantities such as area and volume;
- optional projection from internal `NodeId` space to public export layouts.

The canonical execution layout is always the dense internal `NodeId` space of
`MorphologicalTree`. Projection is a boundary operation, not the internal
representation used by computers.

## Public API

For ordinary C++ use, include:

```cpp
#include <mmcfilters/attributes/Attributes.hpp>
```

Main entry points:

- `AttributeComputation`: public facade for computing one attribute,
  one group, or a heterogeneous request.
- `AttributeNames`: describes column layout in flat attribute buffers.
- `ComputedAttributeData` and `ComputedAttributeDataWithDelta`: owning result
  types returned by the facade.

Use `AttributeComputation` for normal application code. Concrete computers are
advanced-public implementation components, not an alternate public orchestration
path.

For the tree ownership, altitude, and `NodeId` model that attribute buffers use,
see [Morphological Trees](trees.md).
For the reconstruction operators that consume node-indexed attribute buffers,
see [Attribute Filters, Extinction Values, And UAO](filters.md).

## Tree Contracts

Topology-only requests may run on `MorphologicalTree`. Requests that read
altitude must use `WeightedMorphologicalTree<T>` or `WeightedTreeView<T>`.

The attribute layer relies on the tree contracts documented in
[Morphological Trees](trees.md):

- `WeightedMorphologicalTree<T>` owns topology plus a dense altitude buffer.
- `WeightedTreeView<T>` borrows topology plus an external altitude span.
- public computation checks that the tree is not inside an editing session;
- altitude buffers must match the number of internal node slots;
- attributes requiring adjacency must check that adjacency metadata is present.

Individual attributes may add stricter checks. For example, `MAX_DIST` requires
max-tree or min-tree topology with adjacency metadata and rejects unsupported
tree kinds explicitly.

## Common Usage Patterns

Choose the public entry point from the input contract:

- use weighted computation when any requested attribute may read altitude;
- use topology computation only when the request is explicitly support/topology
  only;
- choose non-default output spaces only at API boundaries.

Single altitude-aware attribute:

```cpp
auto [names, values] =
    AttributeComputation::computeSingleAttribute(weightedTree, LEVEL);
```

Several scalar attributes or groups in one coordinated request:

```cpp
auto [names, values] = AttributeComputation::computeAttributes(
    weightedTree,
    std::vector<AttributeOrGroup>{AREA, LEVEL, AttributeGroup::GRAY_LEVEL});
```

Topology-only attributes without requiring an altitude-bearing owner:

```cpp
auto [names, values] = AttributeComputation::computeTopologyAttributes(
    tree,
    std::vector<AttributeOrGroup>{AREA, BOX_WIDTH, BALANCE_NODE});
```

Returning values in a preserved public node-id space, for a tree created by
`createFromHigraParent(...)` and not edited since import:

```cpp
auto [names, values] = AttributeComputation::computeAttributes(
    importedTree,
    std::vector<AttributeOrGroup>{AREA, LEVEL},
    NodeIdSpace::HIGRA);
```

Delta-augmented sampling around one scalar attribute:

```cpp
auto [names, values] =
    AttributeComputation::computeSingleAttributeWithDelta(
        weightedTree,
        LEVEL,
        AltitudeDiff<std::uint8_t>{1},
        2);
```

Projecting node attributes back to another domain:

```cpp
auto mapped = AttributeComputation::computeAttributeMapping(weightedTree, AREA);

auto internal = AttributeComputation::computeAttributes(
    weightedTree,
    std::vector<AttributeOrGroup>{AREA, LEVEL});
auto exported = AttributeComputation::projectNodeValuesToExportedHigra(
    weightedTree,
    internal.attributeNames(),
    internal.values());
```

Python follows the same split between weighted and topology-only entry points:

```python
level_by_node = mmcfilters.Attribute.computeSingleAttribute(
    weightedTree,
    mmcfilters.Attribute.LEVEL,
)
names, values = mmcfilters.Attribute.computeAttributes(
    weightedTree,
    [mmcfilters.Attribute.AREA, mmcfilters.Attribute.Group.GRAY_LEVEL],
)
topology_names, topology_values = mmcfilters.Attribute.computeTopologyAttributes(
    weightedTree,
    [mmcfilters.Attribute.AREA, mmcfilters.Attribute.Group.BOUNDARY],
)
```

## Python Surface

Python keeps a smaller public surface than C++:

- image factories currently expose the canonical `np.uint8`
  `WeightedMorphologicalTree` path;
- `Attribute.computeSingleAttribute(...)` and
  `Attribute.computeAttributes(...)` are the weighted attribute entry points;
- `Attribute.computeSingleTopologyAttribute(...)` and
  `Attribute.computeTopologyAttributes(...)` are the explicit
  topology/support-only entry points;
- `NodeIdSpace` can be passed to the Python attribute methods when a preserved
  output node-id space is needed;
- `AttributePipeline`, concrete C++ computers, local-event deltas, and bitquad
  delta buffers are not Python API.

## Result Layout And Output Spaces

Attribute results are dense flat buffers interpreted by `AttributeNames`.
The canonical internal layout is:

```text
values[node_id * num_attributes + attribute_column]
```

where `node_id` is an internal dense `MorphologicalTree` node slot. Dead
internal slots keep the default buffer value; consumers that reason about tree
nodes should iterate `tree.getAliveNodeIds()`.

`ComputedAttributeData` and `ComputedAttributeDataWithDelta` also store the
`NodeIdSpace` of the returned buffer. Public computation methods can request a
different public node-id space, but projection always happens after the internal
pipeline has computed the result in `NodeIdSpace::MORPHOLOGICAL_TREE`.

`NodeIdSpace::HIGRA` means the preserved imported Higra node-id domain. It is
available only for trees imported from Higra whose original node-id space has
not been invalidated by edits. Direct projection to this space copies live
internal-node rows and fills proper-part rows with unit-component values for
the requested attributes.

For a compact Higra layout exported from the current tree, use
`AttributeComputation::projectNodeValuesToExportedHigra(...)`. That
helper emits the `[proper parts | live internal nodes]` layout produced by
hierarchy export and fills unit proper-part rows through the responsible
attribute computers. `computeAttributeMapping(...)` is the image-domain helper:
each proper part receives the value stored at its owner node.

For the distinction between preserved imported Higra ids and exported compact
Higra snapshots, see [Higra Interoperability](higra-interoperability.md).

## Execution Model

At a high level, a request follows this path:

```text
request -> expand groups -> validate support -> materialize dependencies
        -> compute buffers -> assemble requested result -> project if needed
```

The important rules are:

- group requests are expanded to scalar attributes before result assembly;
- hidden dependencies may be computed without being returned to the caller;
- internal buffers are row-major by internal node id and attribute column;
- unsupported requests fail explicitly instead of returning placeholder values.

The intended internal orchestration path is `detail::AttributePipeline`. It
computes typed altitude attributes through typed kernels and delegates
altitude-independent geometry/topology attributes to the topology backend.
New code should route ordinary attribute requests through this path instead of
adding another top-level execution pipeline.

## Bottom-Up Accumulation And Delta Computation

Several attribute computers use bottom-up accumulation over the tree, but the
public facade should not be read as a live incremental-update engine for every
attribute. The shared traversal skeleton exposes three phases:

- preprocess the current node;
- merge each child result into the parent;
- finalize the current node after all children are processed.

`computeSingleAttributeWithDelta` builds on the normal attribute computation:
it computes the base attribute, then materializes ancestor/descendant offsets
using a typed altitude step, a radius, and a padding policy.

After structural edits, call the public computation facade again unless the
specific operator owns a separate edit-aware update path.

## Dependencies

Dependencies are ordinary attribute results consumed by another attribute
computer. They are passed as `DependencySource`, a non-owning pair of:

- `AttributeNames`, describing the dependency layout;
- `const float*`, pointing to the dependency buffer.

The current implementation has one orchestration path:

- `AttributePipeline` materializes common hidden dependencies directly, such as
  area or volume when needed by another attribute;
- the topology backend computes topology-only families explicitly, including
  bounding-box, tree-topology, moments, bitquads, and contour attributes;
- `AttributeDependencyCache` stores reusable internal-space buffers only for
  dependencies that have already been materialized by the pipeline/backend.

Dependency buffers are reusable only when they are valid, contain the requested
attributes, and are expressed in `NodeIdSpace::MORPHOLOGICAL_TREE`.

## Computers

Concrete computers live under `mmcfilters/attributes/computers/`. A computer
owns a coherent attribute family, not necessarily a single scalar descriptor.

Each computer is expected to:

- declare the attributes it can produce;
- validate any dependency sources required by the requested descriptors;
- respect the requested subset passed by the pipeline;
- write into the caller-owned flat result buffer;
- define unit-component values for compact Higra-style export projection.

Computers are treated as stateless function objects. Implementations must not
keep request-specific mutable state inside the computer object.

Public bucket types that may be useful to direct computer tests live outside
`detail`, such as `BitquadFamilyCounts` and `ContourSideCounts`. Local deltas,
state histograms, traversal policies, and other implementation-only helpers
remain under `computers/detail` and are not public contracts.

## Header-Only Boundary

The C++ library is header-only. Installed packages therefore include internal
`detail` headers whenever public headers need them as transitive implementation
dependencies. Those files are shipped so public headers compile for downstream
consumers; they are not part of the compatibility contract.

Ordinary consumer code should include documented public facade headers such as
`mmcfilters/attributes/Attributes.hpp`. Advanced extension or test code may
include concrete computers under `mmcfilters/attributes/computers/` directly
when it is intentionally targeting the computer protocol. Public bucket types
that are intentionally shareable live outside `detail`, for example
`attributes::computers::BitquadFamilyCounts` and
`attributes::computers::ContourSideCounts`.

## Registry

`AttributeRegistry.hpp` is the source of attribute metadata:

- public name;
- description;
- group membership;
- whether an attribute requires altitude;
- whether an attribute is topology-only.

Group membership is registry metadata. Public requests may mix scalar
attributes and groups; the pipeline expands groups, deduplicates scalars, and
assembles a result containing the requested scalar attributes. Current public
groups are `GRAY_LEVEL`, `SHAPE`, `MOMENTS`, `BOUNDARY`, `TREE_TOPOLOGY`, and
`ALL`.

## Adding An Attribute

Use this checklist for new attributes:

1. Add the scalar enum and registry metadata.
2. Decide whether the attribute is topology-only, altitude-aware, adjacency
   dependent, or tree-kind specific.
3. Add the attribute to any coherent group only if the group semantics still
   hold.
4. Choose whether it belongs in a typed pipeline kernel or a concrete computer.
5. Declare and validate dependencies explicitly.
6. Implement unit-component behavior for exported layouts.
7. Register the attribute in the pipeline/topology backend route.
8. Add small C++ tests with hand-checkable trees/images.
9. Update Python bindings, notebooks, or examples only if the public surface
   changes.

Prefer the smallest coherent implementation path. If an attribute shares a
traversal or intermediate state with an existing family, extend that family
instead of adding a separate computer.

## Validation

The attribute boundary is covered by focused unit tests plus the installed
consumer test. Useful checks while changing this subsystem are:

```bash
cmake --build build --target \
  unit_public_attribute_api \
  unit_attribute_plumbing \
  unit_attribute_unit_values \
  unit_attributes_on_morphological_tree \
  unit_local_event_computations \
  unit_maxdist_support

ctest --test-dir build --output-on-failure -R \
  "unit_(public_attribute_api|attribute_plumbing|attribute_unit_values|attributes_on_morphological_tree|local_event_computations|maxdist_support|installed_consumer)"
```

Run the Python tests when the Python facade or bindings change.

## Attribute Catalog

Use the constant directly in C++ requests, for example `AREA`, or through the
Python facade, for example `mmcfilters.Attribute.AREA`. Every row is materialized
per node in the returned `AttributeNames` layout.

The `Contract` column classifies the input contract required by the attribute:

- `Topology/support`: the attribute uses the tree support in the image domain,
  such as proper-part ownership, pixel coordinates, contours, bounding boxes, or
  shape moments. It does not read node altitudes.
- `Altitude-aware`: the attribute reads the node altitude buffer and therefore
  requires a `WeightedMorphologicalTree<T>` or `WeightedTreeView<T>`. Some
  altitude-aware attributes also depend on support attributes such as `AREA`.
- `Tree topology`: the attribute uses only parent/child relations in the
  hierarchy. It is a stricter topology-only case: no image geometry and no
  altitude values are required.

The `Groups` column lists non-`ALL` group memberships. Rows are sorted by the
first listed group; attributes that belong to more than one group appear only
once. `ALL` expands to every attribute in this table.

| Constant | Groups | Contract | Description |
| --- | --- | --- | --- |
| `VOLUME` | `GRAY_LEVEL` | Altitude-aware | Sum of altitude-weighted support contributions over the node subtree. It behaves like the grey-level mass or integral of the image over the connected component support. |
| `RELATIVE_VOLUME` | `GRAY_LEVEL` | Altitude-aware | Cumulative grey-level contrast volume. It accumulates absolute parent/child altitude jumps weighted by child support area, then adds the node area term used by the current implementation. |
| `LEVEL` | `GRAY_LEVEL` | Altitude-aware | Altitude of the node in the morphological hierarchy. For component trees, this is the grey level at which the connected component appears. |
| `GRAY_HEIGHT` | `GRAY_LEVEL` | Altitude-aware | Grey-level span from the node altitude to the most extreme descendant altitude: maximum descendant level in a max-tree, minimum descendant level in a min-tree. Leaves have value `0`. |
| `MEAN_LEVEL` | `GRAY_LEVEL` | Altitude-aware | Average grey level over the full node support. It is computed from accumulated `VOLUME / AREA`. |
| `VARIANCE_LEVEL` | `GRAY_LEVEL` | Altitude-aware | Variance of grey levels over the full node support. It uses the accumulated squared grey-level sum and the mean level. |
| `AREA` | `SHAPE` | Topology/support | Number of proper parts in the full node support, including all descendant supports. In image-domain trees this is the pixel count of the connected component represented by the node. |
| `BOX_WIDTH` | `SHAPE` | Topology/support | Width, in columns, of the smallest axis-aligned bounding box enclosing the node support. |
| `BOX_HEIGHT` | `SHAPE` | Topology/support | Height, in rows, of the smallest axis-aligned bounding box enclosing the node support. |
| `DIAGONAL_LENGTH` | `SHAPE` | Topology/support | Euclidean diagonal length of the bounding box, `sqrt(width^2 + height^2)`. |
| `RECTANGULARITY` | `SHAPE` | Topology/support | Ratio `AREA / (BOX_WIDTH * BOX_HEIGHT)`. Values closer to `1` indicate that the support fills its bounding box densely. |
| `RATIO_WH` | `SHAPE` | Topology/support | Bounding-box aspect ratio. The implementation returns `max(width, height) / min(width, height)` for non-degenerate boxes, so values are at least `1`. |
| `BOX_COL_MIN` | `SHAPE` | Topology/support | Minimum image column index covered by the node support. |
| `BOX_COL_MAX` | `SHAPE` | Topology/support | Maximum image column index covered by the node support. |
| `BOX_ROW_MIN` | `SHAPE` | Topology/support | Minimum image row index covered by the node support. |
| `BOX_ROW_MAX` | `SHAPE` | Topology/support | Maximum image row index covered by the node support. |
| `MAX_DIST` | `SHAPE` | Altitude-aware | Maximum squared Euclidean distance reached from the node contour during the incremental distance-transform sweep. It is currently defined for max-trees and min-trees with valid adjacency metadata. |
| `CENTRAL_MOMENT_20` | `MOMENTS`, `SHAPE` | Topology/support | Second-order central moment `mu20` around the support centroid. It measures horizontal spread using column coordinates as `x`. |
| `CENTRAL_MOMENT_02` | `MOMENTS`, `SHAPE` | Topology/support | Second-order central moment `mu02` around the support centroid. It measures vertical spread using row coordinates as `y`. |
| `CENTRAL_MOMENT_11` | `MOMENTS`, `SHAPE` | Topology/support | Mixed second-order central moment `mu11`. It measures covariance between column and row coordinates. |
| `CENTRAL_MOMENT_30` | `MOMENTS`, `SHAPE` | Topology/support | Third-order central moment `mu30`. It captures horizontal asymmetry of the support around its centroid. |
| `CENTRAL_MOMENT_03` | `MOMENTS`, `SHAPE` | Topology/support | Third-order central moment `mu03`. It captures vertical asymmetry of the support around its centroid. |
| `CENTRAL_MOMENT_21` | `MOMENTS`, `SHAPE` | Topology/support | Mixed third-order central moment `mu21`. It captures combined horizontal spread and vertical asymmetry. |
| `CENTRAL_MOMENT_12` | `MOMENTS`, `SHAPE` | Topology/support | Mixed third-order central moment `mu12`. It captures combined vertical spread and horizontal asymmetry. |
| `HU_MOMENT_1` | `MOMENTS`, `SHAPE` | Topology/support | First Hu invariant computed from normalized central moments. In the current implementation it is numerically the same scalar as `INERTIA`, but it is materialized by `HuMomentsComputer`; `INERTIA` is materialized by `MomentBasedAttributeComputer`. |
| `HU_MOMENT_2` | `MOMENTS`, `SHAPE` | Topology/support | Second Hu invariant. It emphasizes anisotropy between horizontal and vertical spread while remaining invariant to translation, scale, and rotation. |
| `HU_MOMENT_3` | `MOMENTS`, `SHAPE` | Topology/support | Third Hu invariant. It captures third-order skewness and asymmetry patterns in the support distribution. |
| `HU_MOMENT_4` | `MOMENTS`, `SHAPE` | Topology/support | Fourth Hu invariant. It combines third-order normalized moments to describe diagonal and off-axis symmetry. |
| `HU_MOMENT_5` | `MOMENTS`, `SHAPE` | Topology/support | Fifth Hu invariant. It is sensitive to more complex orientation-dependent asymmetries and reflection-related differences. |
| `HU_MOMENT_6` | `MOMENTS`, `SHAPE` | Topology/support | Sixth Hu invariant. It combines second- and third-order normalized moments to capture elliptical asymmetry and curvature-related shape variation. |
| `HU_MOMENT_7` | `MOMENTS`, `SHAPE` | Topology/support | Seventh Hu invariant. It is highly sensitive to fine asymmetries and helps distinguish mirror-related shapes. |
| `INERTIA` | `MOMENTS`, `SHAPE` | Topology/support | Sum of normalized second-order central moments, `mu20 / area^2 + mu02 / area^2`. In the current implementation this equals `HU_MOMENT_1`, but it is materialized by `MomentBasedAttributeComputer` and belongs to the moment-derived descriptor family. |
| `COMPACTNESS` | `MOMENTS`, `SHAPE` | Topology/support | Area normalized by second-order dispersion, currently `(1 / (2*pi)) * area / (mu20 + mu02)` when the denominator is positive. Higher values indicate more compact supports. |
| `ECCENTRICITY` | `MOMENTS`, `SHAPE` | Topology/support | Ratio of the largest to smallest eigenvalue of the second-moment matrix. Values near `1` indicate isotropic shapes; larger values indicate elongation. |
| `LENGTH_MAJOR_AXIS` | `MOMENTS`, `SHAPE` | Topology/support | Length proxy for the major axis of the equivalent second-moment ellipse, derived from the largest inertia eigenvalue and area. |
| `LENGTH_MINOR_AXIS` | `MOMENTS`, `SHAPE` | Topology/support | Length proxy for the minor axis of the equivalent second-moment ellipse, derived from the smallest inertia eigenvalue and area. |
| `AXIS_ORIENTATION` | `MOMENTS`, `SHAPE` | Topology/support | Principal-axis orientation in degrees, computed as `0.5 * atan2(2*mu11, mu20 - mu02)` and normalized by the implementation to a non-negative angle. |
| `CIRCULARITY` | `MOMENTS`, `SHAPE` | Topology/support | Ratio `lambda2 / lambda1` of second-moment eigenvalues. Values near `1` indicate circular or isotropic supports; values near `0` indicate elongation. |
| `BITQUADS_AREA` | `BOUNDARY`, `SHAPE` | Topology/support | Duda-style sub-pixel area estimator derived from aggregated `2x2` bitquad pattern counts. |
| `BITQUADS_NUMBER_EULER` | `BOUNDARY`, `SHAPE` | Topology/support | Euler characteristic estimated from bitquad counters, representing connected components minus holes under the selected connectivity projection. |
| `BITQUADS_NUMBER_HOLES` | `BOUNDARY`, `SHAPE` | Topology/support | Number of holes inferred from the bitquad Euler characteristic for a single connected support. |
| `BITQUADS_PERIMETER` | `BOUNDARY`, `SHAPE` | Topology/support | Discrete boundary-length estimate from bitquad edge-contributing patterns. |
| `BITQUADS_PERIMETER_CONTINUOUS` | `BOUNDARY`, `SHAPE` | Topology/support | Smoothed continuous perimeter estimate from bitquad counters, using weighted transitions across local `2x2` configurations. |
| `BITQUADS_CIRCULARITY` | `BOUNDARY`, `SHAPE` | Topology/support | Bitquad compactness measure `(4*pi*BITQUADS_AREA) / BITQUADS_PERIMETER_CONTINUOUS^2`. Values closer to `1` indicate rounder supports. |
| `BITQUADS_PERIMETER_AVERAGE` | `BOUNDARY`, `SHAPE` | Topology/support | Average continuous perimeter per connected component, computed from bitquad perimeter and Euler-count estimates. |
| `BITQUADS_LENGTH_AVERAGE` | `BOUNDARY`, `SHAPE` | Topology/support | Average longitudinal extent proxy, derived as half of the average continuous perimeter. |
| `BITQUADS_WIDTH_AVERAGE` | `BOUNDARY`, `SHAPE` | Topology/support | Average transverse extent proxy, derived from average bitquad area and average perimeter. |
| `CONTOUR_PIXELS` | `BOUNDARY`, `SHAPE` | Topology/support | Number of support pixels that touch the 4-neighbour complement by at least one side. |
| `CONTOUR_PERIMETER` | `BOUNDARY`, `SHAPE` | Topology/support | Total number of exposed 4-neighbour sides over the support. This is the side-count perimeter, not an Euclidean perimeter estimate. |
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
| `LEAF_RATIO_NODE` | `TREE_TOPOLOGY` | Tree topology | Ratio of leaf descendants to the subtree size used by the implementation, `leaf_descendants / (descendants + 1)`. Leaves return `1`. |
| `BALANCE_NODE` | `TREE_TOPOLOGY` | Tree topology | Difference between maximum and minimum child-subtree heights. It is `0` for leaves and increases when child depths are uneven. |
| `AVG_CHILD_HEIGHT_NODE` | `TREE_TOPOLOGY` | Tree topology | Average height of direct child subtrees. Leaves return `0`; non-leaf nodes average the child height values computed in the topology pass. |
