# Saliency maps

This guide distinguishes the edge-indexed hierarchy and shape-space saliency
operators exposed by `mmcfilters`. It defines their inputs, outputs, validity,
determinism, and complexity contracts.

## Choose an operator

| Need | API | Result |
| --- | --- | --- |
| project a monotone node valuation | `HierarchySaliencyMap::computeSaliencyEdgeMap` | value on every graph edge |
| rank only values realized on graph transitions | `computeCanonicalRankedSaliencyEdgeMap` | ranked value on every graph edge |
| use structural hierarchy levels | `computeTopologicalLevelEdgeMap` | ranked value on every graph edge |
| normalize max-tree or min-tree altitude | `computeNormalizedAltitudeEdgeMap` | value in `[0, 1]` on every graph edge |
| build a persistence hierarchy from extinction values | `ExtinctionValues::computeFormalSaliencyEdgeMap` | value on every graph edge |
| project the max-propagated extinction valuation directly | `computeMonotoneExtinctionProjection` | value on every graph edge |
| draw selected extinction cutoff contours | `ExtinctionValues::contourMap` | pixel visualization |
| compute extinction from an arbitrary tree-node attribute | `ShapeSpaceSaliency` | value on every graph edge |

`contourMap` is a raster visualization; it is not an edge-indexed hierarchy
saliency map.

## Hierarchy model

Let \f$E\f$ be a finite pixel domain and
\f$G=(E,\mathcal{A})\f$ an undirected graph. A rooted hierarchy is

\f[
\mathcal{T}=(\mathcal{V},r,{\rm par}),
\f]

where every pixel \f$p\in E\f$ has one smallest node
\f$P(p)\in\mathcal{V}\f$. A compatible valuation
\f$h:\mathcal{V}\rightarrow\mathbb{R}_{\ge 0}\f$ is non-decreasing toward the
root:

\f[
h({\rm par}(v))\ge h(v).
\f]

For an edge \f$\{p,q\}\f$, let
\f$a={\rm LCA}(P(p),P(q))\f$ be the lowest common ancestor (LCA). Under
`HierarchyLevelConvention::EdgeSaliencyValue`,

\f[
\Phi_h(\{p,q\})=
\left\{\begin{array}{ll}
0, & P(p)=P(q),\\
h(a), & P(p)\ne P(q).
\end{array}\right.
\f]

Under `PartitionAppearanceLevel`, a transition edge receives \f$h(a)-1\f$.
These conventions represent a supplied hierarchy on the graph; they do not
construct a hierarchy from arbitrary graph-edge weights.

## Direct hierarchy projection

`HierarchySaliencyMap` enumerates each undirected adjacency edge once, obtains
its endpoints' smallest nodes, computes their LCA, and writes the selected hierarchy
level. `EdgeSaliencyMap<T>` stores parallel `sources`, `targets`, and `values`
arrays with image-domain metadata.

Edges whose endpoints have the same smallest node are present with value zero.
Formal valuations must be finite and non-negative. Equal parent/child levels
are accepted by `AllowLevelCollapse`; `RequireStrictHierarchy` rejects them.

Use:

- `rankHierarchyValuation` for dense ranks over all live nodes;
- `computeCanonicalRankedSaliencyEdgeMap` for dense ranks over values realized
  on graph transitions;
- `computeTopologicalLevelEdgeMap` for a structural scale;
- `computeNormalizedAltitudeEdgeMap` for polarity-aware max-tree or min-tree
  altitude in `[0, 1]`.

The normalized-altitude operation requires one global max-tree or min-tree
polarity. A self-dual residual tree or a tree of shapes may instead use a
supplied non-decreasing structural or attribute valuation.

For component trees, `ComponentTreePartitionHierarchyAdapter` completes the
partial-partition interpretation: pixel singletons form partition zero,
zero-valued edges connect pixels within proper parts, and component nodes merge
those components with child supports.

## Extinction persistence

`ExtinctionValues<T, Real>` computes one record per component-tree leaf. Each
record contains the leaf, its cutoff node, and its extinction. The dominant
extremum receives the finite ordering sentinel
`std::numeric_limits<Real>::max()`.

`computeFormalSaliencyEdgeMap` constructs an edge hierarchy in five steps:

1. order graph edges by finest-region status, component-tree altitude, and
   endpoint IDs;
2. select a deterministic altitude-ordered minimum spanning tree;
3. assign each binary merge the persistence

   \f[
   {\rm pers}(L,R)=\min(M_L,M_R),
   \qquad M_{L\cup R}=\max(M_L,M_R),
   \f]

   where \f$M_L\f$ and \f$M_R\f$ are maximum descendant extinctions;
4. build the quasi-flat-zone dendrogram of the persistence-weighted tree;
5. project the dendrogram to every graph edge by LCA.

`computeMonotoneExtinctionProjection` performs a different operation: it
max-propagates extinction values over the supplied component tree and projects
that monotone valuation directly by LCA.

## Shape-space saliency

`ShapeSpaceSaliency` accepts an arbitrary finite floating-point attribute on the
live nodes of an input tree. It treats those nodes as graph vertices and the
parent/child relations as graph edges. Equal-valued connected nodes form
plateaus; regional minima or maxima receive finite extinction values in this
second graph.

Let \f$C(v)\subseteq E\f$ be the support of node \f$v\f$. Its graph boundary is
the cut

\f[
\partial_{\mathcal A}C(v)=
\{\{p,q\}\in\mathcal A : |\{p,q\}\cap C(v)|=1\}.
\f]

An edge belongs to this boundary exactly when one endpoint is in the support.
For non-negative node scores \f$s(v)\f$, projection takes the maximum over all
boundaries containing the edge:

\f[
w(\{p,q\})=
\max\bigl(\{0\}\cup
\{s(v) : v\in\mathcal V,\ \{p,q\}\in\partial_{\mathcal A}C(v)\}\bigr).
\f]

Equivalently, the maximum is taken over the two smallest-node-to-LCA paths, excluding
the LCA. An edge whose endpoints have the same smallest node receives zero.
Use this operation when the input attribute is not a monotone valuation of the
original hierarchy.

## Validity contracts

All operations require a committed rooted topology. Cached views and helper
objects reject use after topology mutation.

- Node buffers contain one slot per internal `NodeId`; live values consumed by
  an operation must be finite.
- The graph domain must match the tree's rows, columns, and pixels.
- A stored adjacency is used only when it defines one unambiguous graph;
  otherwise the caller supplies an explicit relation.
- Formal LCA projection requires a finite, non-negative valuation that is
  non-decreasing toward the root.
- Formal hierarchy connectivity requires every node support to be connected in
  the selected graph.
- Extinction persistence additionally requires globally monotone component-tree
  altitude and a non-empty proper part for each leaf.
- Shape-space contour projection requires the tree pixel domain to match the image grid. Its
  input attribute may be negative, but projected scores must be finite and
  non-negative.

Invalid domains, disconnected supports, stale topology, invalid LCAs, negative
formal values, and non-finite values are rejected.

## Determinism

- Extinction records are ordered by decreasing extinction, increasing cutoff
  `NodeId`, then increasing leaf `NodeId`.
- Equal-strength extinction branches follow deterministic leaf and traversal
  order.
- Equal hierarchy edges use row-major endpoint IDs.
- Shape-space plateaus use their outermost node as representative. Ties between
  extrema and result ordering use the minimum pixel index in the node support,
  then support cardinality and depth, independently of internal `NodeId` assignment.

These rules make results reproducible without claiming mathematical uniqueness
under other valid tie policies.

## Cuts and display

A threshold cut at \f$\lambda\f$ selects edges with
\f$w(e)\ge\lambda\f$. Edges with \f$w(e)<\lambda\f$ define the associated
quasi-flat-zone components.

`HierarchySaliencyMapProjection::edgeMapToPixelImage` aggregates incident edge
values for display. `ExtinctionValues::contourMap` draws selected cutoff-node
contours. Neither raster should be compared with an edge map without stating
the aggregation rule.

## Complexity

Let `N` count internal node slots, \f$m=|\mathcal{V}|\f$ live nodes,
\f$p=|E|\f$ pixels, \f$e=|\mathcal{A}|\f$ graph edges, and \f$l\f$
component-tree leaves. Edits can leave dead slots, so `m <= N`; dense buffers
are sized by `N`. Let `alpha` denote the inverse Ackermann function.

The following bounds exclude source-tree construction and its LCA caches.
Projection rows also exclude hierarchy-connectivity validation, listed
separately. Direct projection validates connectivity by default; extinction
persistence always validates it.

| Operation | Time | Additional storage, including output |
| --- | --- | --- |
| hierarchy-connectivity validation | `O(N + (m + p + e) alpha(p))` | `O(N + p + e)` |
| direct LCA projection | `O(N + p + e)` | `O(N + e)` |
| extinction persistence map | `O(N + e log(e + 1) + p log(p + 1))` | `O(N + e + p log(p + 1))` |
| extinction record construction | conservative `O(N + lm)` | `O(N + l)` |
| shape-space extinction | `O(N + p + m log(m + 1))` | `O(N)` |
| shape-space contour projection | `O(N log(N + 1) + p + e log(N + 1))` | `O(N log(N + 1) + e)` |

These projection paths and connectivity validation use scalar
`MorphologicalTree::lowestCommonAncestor` queries. DFS intervals take `O(N)`
time and storage to prepare and resolve comparable pairs in constant time.
The first incomparable pair may additionally build an Euler/RMQ cache with
`O(N log(N + 1))` time and storage. Later queries reuse it in constant time.
Add this preparation cost when the source cache is absent; retained cache
storage belongs to the tree. See [Morphological trees](trees.md).

Each extinction-persistence call constructs a new dendrogram with `O(p)` nodes;
the table includes its possible `O(p log(p + 1))` LCA cache. Extinction-record
construction is a separate cost paid when initializing `ExtinctionValues`.
Shape-space extinction includes preparation of the support cardinality and
minimum pixel index in the node support for each node. Retaining multiple maps
adds their `O(e)` outputs.

## C++ and Python entry points

| Operation | C++ | Python |
| --- | --- | --- |
| direct LCA projection | `HierarchySaliencyMap::computeSaliencyEdgeMap` | `HierarchySaliencyMap.compute_saliency_edge_map` |
| canonical transition ranks | `computeCanonicalRankedSaliencyEdgeMap` | `compute_canonical_ranked_saliency_edge_map` |
| extinction persistence | `ExtinctionValues::computeFormalSaliencyEdgeMap` | `extinction.compute_formal_saliency_edge_map` |
| direct extinction projection | `computeMonotoneExtinctionProjection` | `extinction.compute_monotone_extinction_projection` |
| cutoff contour visualization | `ExtinctionValues::contourMap` | `extinction.contour_map` |
| shape-space extinction | `ShapeSpaceSaliency` | `mmcfilters.ShapeSpaceSaliency` |
| cuts and display | `HierarchySaliencyMapProjection` | `mmcfilters.HierarchySaliencyMapProjection` |

## Scientific references

- A. G. Silva and R. A. Lotufo, “Efficient computation of new extinction values
  from extended component tree,” *Pattern Recognition Letters* 32(1), 79–90,
  2011. [DOI](https://doi.org/10.1016/j.patrec.2010.07.019).
- J. Cousty, L. Najman, Y. Kenmochi, and S. Guimarães, “Hierarchical
  segmentations with graphs: quasi-flat zones, minimum spanning trees, and
  saliency maps,” *JMIV* 60(4), 479–502, 2018.
  [DOI](https://doi.org/10.1007/s10851-017-0768-7).
- Y. Xu, E. Carlinet, T. Géraud, and L. Najman, “Hierarchical Segmentation Using
  Tree-Based Shape Spaces,” *IEEE TPAMI* 39(3), 457–469, 2017.
  [DOI](https://doi.org/10.1109/TPAMI.2016.2554550).
