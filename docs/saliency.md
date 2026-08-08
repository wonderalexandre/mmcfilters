# Hierarchy Saliency Maps

This is the canonical scientific guide to saliency in `mmcfilters`. It defines
the mathematical objects, distinguishes the available operators, and records
the validity, determinism, and complexity contracts shared by the C++ and
Python APIs. The tree, filter, and Python guides contain API-specific examples
and link back here instead of redefining these concepts.

## Mathematical model

Let \f$E\f$ be a finite proper-part domain, usually the row-major pixels of an
image, and let \f$G=(E,\mathcal{A})\f$ be the undirected graph induced by a
`RegularGridAdjacency2D`. A rooted morphological tree is

\f[
\mathcal{T}=(\mathcal{V},r,{\rm par}),
\f]

where \f$r\f$ is the root. Every proper part \f$p\in E\f$ has one direct owner
\f$o(p)\in\mathcal{V}\f$. A hierarchy valuation
\f$h:\mathcal{V}\rightarrow\mathbb{R}_{\ge 0}\f$ is compatible when it is
non-decreasing toward the root:

\f[
h({\rm par}(v))\ge h(v).
\f]

For an adjacency edge \f$\{p,q\}\f$, let
\f$a={\rm LCA}(o(p),o(q))\f$. Under
`HierarchyLevelConvention::EdgeSaliencyValue`, the full edge map is

\f[
\Phi_h(\{p,q\})=
\left\{\begin{array}{ll}
0, & o(p)=o(q),\\
h(a), & o(p)\ne o(q).
\end{array}\right.
\f]

Under `PartitionAppearanceLevel`, positive partition-appearance levels are
interpreted literally as in Cousty's construction and a transition edge
receives \f$h(a)-1\f$. This is the `Phi(H)` representation of a connected
hierarchy. It is distinct from the full `Psi(w) = Phi(QFZ(G,w))` pipeline, which
starts from arbitrary graph-edge weights.

## Choose the intended operator

| Goal | API | Mathematical operation | Result domain |
| --- | --- | --- | --- |
| Project a known monotone hierarchy valuation | `HierarchySaliencyMap::computeSaliencyEdgeMap` | LCA projection `Phi(H)` | every graph edge |
| Obtain dense ranks only for effective transitions | `computeCanonicalRankedSaliencyEdgeMap` | LCA projection plus effective-edge ranking | every graph edge |
| Use a structural hierarchy scale | `computeTopologicalLevelEdgeMap` | canonically ranked topological levels | every graph edge |
| Use max/min-tree gray-level altitude | `computeNormalizedAltitudeEdgeMap` | polarity-aware normalization followed by LCA | every graph edge |
| Build a Cousty hierarchical watershed from extinction values | `ExtinctionValues::computeFormalSaliencyEdgeMap` / `HierarchicalWatershedSaliency` | altitude-ordered MST, persistence, QFZ dendrogram, full-graph projection | every graph edge |
| Reproduce the former extinction experiment | `computeMonotoneExtinctionProjection` | max-descendant extinction valuation followed by LCA | every graph edge |
| Draw selected extinction cutoff contours | `ExtinctionValues::contourMap` | raster visualization of selected node contours | pixels, not a formal saliency map |
| Apply Xu shaping to a non-monotone node attribute | `ShapeSpaceSaliency` | second component hierarchy on the tree-node graph, extinction, maximum on original contours | every graph edge |

The formal extinction method and the monotone projection are not equivalent.
They can return different edge values for the same extinction records.

## Primary references and implementation correspondence

The saliency subsystem combines three algorithms from distinct primary
publications. The references below are normative for the scientific operations;
the implementation-specific contracts in this guide remain authoritative for
accepted inputs, deterministic ties, storage, and exceptions.

1. **Component-tree extinction records.** Alexandre Gonçalves Silva and Roberto
   de Alencar Lotufo, “Efficient computation of new extinction values from
   extended component tree,” *Pattern Recognition Letters*, 32(1), 79–90,
   2011. [DOI 10.1016/j.patrec.2010.07.019](https://doi.org/10.1016/j.patrec.2010.07.019).
   `ExtinctionValues::initialize` follows the leaf-to-root branch traversal of
   Algorithm 1 for an increasing attribute. In `mmcfilters`, the tree and its
   dense attribute buffer already exist; they are not built incrementally with
   the extinction computation. The finite dominant-extremum sentinel and the
   exact equal-strength tie order are library conventions.
2. **Connected hierarchies, QFZs, MSTs, and saliency maps.** Jean Cousty,
   Laurent Najman, Yukiko Kenmochi, and Silvio Guimarães, “Hierarchical
   segmentations with graphs: quasi-flat zones, minimum spanning trees, and
   saliency maps,” *Journal of Mathematical Imaging and Vision*, 60(4),
   479–502, 2018.
   [DOI 10.1007/s10851-017-0768-7](https://doi.org/10.1007/s10851-017-0768-7);
   [author manuscript, HAL hal-01344727v2](https://hal.science/hal-01344727v2).
   `HierarchySaliencyMap` implements the edge-indexed `Phi(H)` correspondence
   of Section 4, Equations (5)–(6), through the LCA algorithm in Section 7,
   Algorithm 1. `HierarchicalWatershedSaliency` concretizes the BPTAO,
   extinction extension, and persistence construction sketched in Section 8.1,
   then uses the paper's QFZ/saliency correspondence to return a full-graph map.
   Cousty et al. start from an arbitrary edge-weighted graph `(G,w)`. This API
   instead starts from a committed max-tree or min-tree: same-owner graph edges
   are ordered first and transition edges are ordered by the altitude of the LCA
   of their endpoint owners. The resulting component-tree-derived edge order is
   then used by the two Kruskal constructions. This is an implementation adapter,
   not an additional construction stated in the paper.
3. **Extinction-based shape-space hierarchy transformation.** Yongchao Xu,
   Edwin Carlinet, Thierry Géraud, and Laurent Najman, “Hierarchical
   Segmentation Using Tree-Based Shape Spaces,” *IEEE Transactions on Pattern
   Analysis and Machine Intelligence*, 39(3), 457–469, 2017.
   [DOI 10.1109/TPAMI.2016.2554550](https://doi.org/10.1109/TPAMI.2016.2554550);
   [author manuscript, HAL hal-01301966v1](https://hal.science/hal-01301966v1).
   `ShapeSpaceSaliency` implements the construction in Section 4.3: extrema in
   the node-weighted tree graph receive extinction values, and every original
   region boundary receives the maximum score of the extrema whose regions
   contain that boundary. The paper describes local minima on a Khalimsky grid;
   the library also exposes the order-dual maxima case and stores the result on
   regular-grid graph edges. Those are explicit generalizations, not claims made
   by the publication.

For the computational background of the contour projection, see Yongchao Xu,
Edwin Carlinet, Thierry Géraud, and Laurent Najman, “Efficient Computation of
Attributes and Saliency Maps on Tree-Based Image Representations,” in
*Mathematical Morphology and Its Applications to Signal and Image Processing*
(ISMM 2015), Lecture Notes in Computer Science 9082, 693–704, Springer, 2015,
[DOI 10.1007/978-3-319-18720-4_58](https://doi.org/10.1007/978-3-319-18720-4_58).
The current implementation uses owner-to-LCA paths and binary lifting rather
than reproducing that paper's incremental storage layout.

## Direct hierarchy projection

`HierarchySaliencyMap` enumerates each undirected edge once through the forward
half of the selected adjacency, obtains the endpoint owners, computes their LCA,
and writes the selected hierarchy level. `EdgeSaliencyMap<T>` stores parallel
`sources`, `targets`, and `values` arrays plus the image-domain metadata.

Same-owner edges are explicitly present with value zero. Consequently, a formal
valuation must be finite and non-negative. Equal parent/child levels are accepted
with `AllowLevelCollapse`; `RequireStrictHierarchy` rejects them. Use
`rankHierarchyValuation` when dense node ranks are needed, and use
`computeCanonicalRankedSaliencyEdgeMap` when only values realized on graph edges
should define the rank scale.

For a component tree, `ComponentTreePartitionHierarchyAdapter` makes the
partial-to-complete partition interpretation explicit: pixel singletons form
partition zero, connected same-owner edges form direct-proper-part atoms, and
component nodes merge those atoms with child supports.

## Extinction hierarchical watershed

`ExtinctionValues<T, Real>` first computes one record for every component-tree
leaf. It climbs the leaf branch until the first merge with a stronger branch or
with an already visited equal-strength branch. The record contains the leaf,
its cutoff node, and the extinction value. The dominant extremum reaches the
root and receives the finite ordering sentinel
`numeric_limits<Real>::max()`.

`HierarchicalWatershedSaliency` then performs the constructive persistence path:

1. Enumerate graph edges and order them by finest-region status, component-tree
   altitude, and row-major endpoint ids.
2. Select a deterministic altitude-ordered MST with Kruskal's algorithm.
3. Traverse its binary merges. If \f$M_L\f$ and \f$M_R\f$ are the maximum descendant
   extinctions on the two sides, assign

   \f[
   {\rm pers}(L,R)=\min(M_L,M_R),
   \qquad M_{L\cup R}=\max(M_L,M_R).
   \f]

4. Build the QFZ dendrogram of the persistence-weighted MST.
5. Project that dendrogram back to every edge of \f$G\f$ by LCA.

The dominant-extremum sentinel is therefore used for ordering records but does
not become an ordinary persistence merger level. The ranked overload ranks only
values present on the final full-graph map.

`computeMonotoneExtinctionProjection` intentionally exposes the earlier
max-descendant valuation followed by direct LCA projection. It is a valid
monotone hierarchy, but it is not the persistence construction above.

## Xu shape-space saliency

`ShapeSpaceSaliency` accepts an arbitrary finite floating-point attribute on the
live nodes of the original tree. It treats those nodes as graph vertices and the
original parent/child relations as graph edges. Equal-valued connected nodes
form plateaus; regional minima or maxima are selected in this second graph and
receive finite extinction values.

The final projection is maximum-on-contours, not LCA valuation. If
\f$s(v)\ge0\f$ is the sparse score of an original-tree region, then

\f[
w(\{p,q\})=
\max\{s(v)\mid \{p,q\}\subseteq\partial C(v)\}.
\f]

Equivalently, the maximum is taken on the two owner-to-LCA paths, excluding the
LCA. Use this API when the input attribute is not a hierarchy valuation.

## Validity contracts

All public saliency paths require a committed rooted topology. Cached views and
`ExtinctionValues` objects reject use after a topology mutation.

- Node-valued buffers contain one slot per dense internal `NodeId`; live-node
  values used by an operation must be finite.
- The image-domain graph must match the tree's row/column and proper-part domain.
  A stored adjacency is used only when it is unambiguous; otherwise the caller
  must pass an explicit relation.
- Formal LCA projection requires non-negative values that are non-decreasing
  toward the root. Strict mode also rejects equal parent/child levels.
- Formal hierarchy connectivity means every node support is connected in the
  selected graph, not merely that the parent array has one root.
- `HierarchicalWatershedSaliency` additionally requires a globally monotone
  max-tree/min-tree altitude, finite non-negative leaf extinctions, a connected
  projection graph, and at least one directly owned proper part for every
  component-tree leaf.
- `ExtinctionValues` models regional extrema as leaves and therefore accepts
  hierarchies with a globally monotone altitude order. Standard trees of shapes
  and self-dual residual trees declare unconstrained altitude order and are
  rejected by this operator.
- `ShapeSpaceSaliency` requires one proper part per image pixel for contour
  projection. The input attribute may be negative, but projected node scores
  must be finite and non-negative.

Disconnected adjacencies, mismatched domains, stale topology, invalid LCAs,
negative formal values, and non-finite values are errors; they are not silently
repaired.

## Determinism and ties

The implementations make all algorithmic ties reproducible:

- Extinction records are sorted by decreasing extinction, then increasing
  cutoff `NodeId`, then increasing leaf `NodeId`. `byTopK(k)` keeps the first
  `k`; `byThreshold(t)` keeps every record with `extinction >= t`.
- Equal-strength extinction branches are resolved by the deterministic leaf and
  visitation order of the dense topology.
- The hierarchical-watershed MST orders equal edges by row-major endpoint ids.
  Equal persistence values retain that stable first-pass order.
- Shape-space plateaus use their outermost node as representative, with smaller
  `NodeId` as the final tie break. Equal-strength extrema also prefer the smaller
  representative id.

These rules are implementation contracts for reproducibility; they do not claim
uniqueness of an MST or of a hierarchy under a different tie convention.

## Cuts and visualization

A threshold cut at \f$\lambda\f$ selects contour edges with
\f$w(e)\ge\lambda\f$. The complementary edges with \f$w(e)<\lambda\f$ define the
quasi-flat-zone connected components. `edgeMapToPixelImage` only aggregates
incident edge values for display and is not the formal representation.
Likewise, `ExtinctionValues::contourMap` draws selected cutoff-node contours into
a raster and may combine several selections on a pixel. Do not compare either
rasterization to a formal edge map without stating the aggregation rule.

## Complexity and memory

Let \f$m=|\mathcal{V}|\f$ be the number of live tree nodes, \f$p=|E|\f$ the number
of proper parts, and \f$e=|\mathcal{A}|\f$ the number of undirected graph edges.

| Operation | Time | Dominant auxiliary memory |
| --- | --- | --- |
| hierarchy-connectivity validation | `O(m + p + e)` | `O(m + p + e)` |
| direct LCA saliency projection | linear preprocessing plus `O(e)` projection | `O(m + e)` including output |
| `HierarchicalWatershedSaliency` | `O(e log e + p log p)` plus linear validation/projection terms | `O(m + p + e)` |
| shape-space extinction | `O(m log m)` for the ordered sweep and tie handling | `O(m)` excluding results |
| shape-space contour projection | `O(m log m + e log m)` | `O(m log m + e)` including binary-lifting tables and output |

For `ExtinctionValues` record construction, the current branch-climbing
implementation can revisit ancestors from different leaves; its conservative
worst-case bound is `O(l m)` for (l) leaves, with `O(m + l)` storage. In
typical component trees, tree construction and adjacency-edge materialization
may dominate the later projections.

## C++ and Python entry points

| Concept | C++ | Python |
| --- | --- | --- |
| formal LCA projection | `HierarchySaliencyMap::computeSaliencyEdgeMap` | `mmcfilters.HierarchySaliencyMap.computeSaliencyEdgeMap` |
| canonical edge ranks | `computeCanonicalRankedSaliencyEdgeMap` | same method name |
| Cousty persistence | `ExtinctionValues::computeFormalSaliencyEdgeMap` | `extinction.computeFormalSaliencyEdgeMap` |
| explicit former behavior | `computeMonotoneExtinctionProjection` | `extinction.computeMonotoneExtinctionProjection` |
| raster cutoff contours | `ExtinctionValues::contourMap` | `extinction.contourMap` |
| Xu shaping | `ShapeSpaceSaliency` | `mmcfilters.ShapeSpaceSaliency` |
| cuts and display projections | `HierarchySaliencyMapProjection` | `mmcfilters.HierarchySaliencyMapProjection` |

The maintained English tutorial is
[`Saliency_Maps_Tutorial.ipynb`](../notebooks/Saliency_Maps_Tutorial.ipynb).
Release-to-release API changes are listed in [`CHANGELOG.md`](../CHANGELOG.md).
