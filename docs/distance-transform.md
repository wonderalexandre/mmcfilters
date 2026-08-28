# Distance-transform attributes

The distance-transform subsystem is topology-only and requires a non-empty
regular 2D tree domain. It measures every node support against its foreground
A4 contour in the original image grid. Tree-construction adjacency and node
altitudes do not change this metric.

Two coordinated families are public:

- unsuffixed names, such as `MAX_DIST`, select the adaptive-A8 DIFT
  approximation;
- names ending in `_EXACT`, such as `MAX_DIST_EXACT`, select the exact
  Euclidean distance transform, whose internal cost field is squared.

There is deliberately no `_APPROX` suffix. `AttributeGroup::DistTransf` in C++
and `Attribute.Group.DIST_TRANSF` in Python contain the 29 approximate scalars.
`AttributeGroup::DistTransfExact` and `Attribute.Group.DIST_TRANSF_EXACT` contain
the 29 matching exact scalars. The exact group enumerator was appended, so the
numeric ordinals of pre-existing groups remain unchanged. Request both groups
explicitly for a paired 58-column study. As with every multi-attribute request,
a computed result layout uses the API's canonical scalar-ordinal order rather
than the registry's presentation order.

## Mathematical contract

For a node support `X`, write `z(x) = d(x)^2`, where `d(x)` is the Euclidean
distance in pixels from `x` to the foreground A4 contour of `X`. In the
approximate family, `z` denotes the corresponding adaptive-A8 DIFT cost. Let
`n = |X|`.

The maximum and squared-distance summaries are:

- `MAX_DIST = sqrt(max z) = max d`;
- `MAX_SQUARED_DIST = max z`;
- `DIST_SQUARED_SUM = sum z`;
- `DIST_SQUARED_MEAN = (sum z) / n`;
- `DIST_RMS = sqrt(DIST_SQUARED_MEAN)`;
- `DIST_SQUARED_VARIANCE = mean(z^2) - mean(z)^2`.

The real-distance summaries use `d = sqrt(z)`:

- `DIST_SUM = sum d`;
- `DIST_MEAN = (sum d) / n`;
- `DIST_VARIANCE = mean(d^2) - mean(d)^2`.

The maximum center is the smallest row-major support pixel among tied
maximizers. Maximum-plateau area counts every maximizer, and plateau centroid
coordinates are arithmetic means, so they may be fractional.

The sparse histogram stores `(z, count)` bins in increasing squared-distance
order. `DIST_MEDIAN`, `DIST_Q25`, `DIST_Q75`, and `DIST_Q90` are lower empirical
quantiles of `d`: the smallest distance whose cumulative count reaches
`ceil(p*n)`. `DIST_MODE` chooses the smallest distance when frequencies tie.
`DIST_ENTROPY` is Shannon entropy in bits. `DIST_POSITIVE_AREA` counts samples
with `z > 0`, and `DIST_LEVEL_COUNT` counts represented squared-distance levels,
including zero.

## Distance-weighted geometry

Spatial descriptors use weight `w(x) = d(x)`. The row and column centroids are
`sum(w*coordinate) / sum(w)`. If every distance is zero, the centroid falls
back to the ordinary support centroid and all weighted central moments are
zero.

`DIST_WEIGHTED_CENTRAL_MOMENT_20` is the unnormalized weighted column spread,
`DIST_WEIGHTED_CENTRAL_MOMENT_02` is the unnormalized weighted row spread,
and `DIST_WEIGHTED_CENTRAL_MOMENT_11` is the unnormalized row-column mixed
moment. Axis orientation is reported in degrees. Eccentricity is the
major/minor second-moment eigenvalue ratio: isotropic or point-degenerate
fields return `1`, and line-degenerate fields saturate at `1e6`.

Every name above has an exact counterpart obtained by appending `_EXACT`.

## Units

| Quantity | Unit |
| --- | --- |
| maximum distance | pixels |
| maximum squared distance, squared mean | pixel² |
| squared sum | pixel² accumulated over support samples |
| squared-distance variance | pixel⁴ |
| real-distance sum, mean, quantiles, mode, RMS | pixels |
| real-distance variance | pixel² |
| entropy | bits |
| areas and level counts | counts |
| centroids | zero-based pixel coordinates |
| weighted central moments | pixel³ |
| orientation | degrees |
| eccentricity | dimensionless |

## Performance behavior

Maximum-only requests remain specialized. `MAX_SQUARED_DIST` exposes the
internal maximum cost directly; `MAX_DIST` adds only one square root per live
node. Both use the lightweight approximate maximum observer, and squared
moments use the smaller dynamic-moment observer. Histogram and spatial state
are independently selected at compile time: distributional requests enable the
sparse histogram only when needed, and spatial-only requests do not allocate or
update a histogram. Mixed scalar requests retain exactly the union of the
required state. Exact bundles share one EDT sample stream, while approximate
bundles share one DIFT traversal.

The approximate DIFT has one production queue: the PQueue32-style `GFT_FAST`
intrusive-node layout. It is not selected through CMake or the public API.
After the DIFT establishes valid dense element and squared-cost domains, the
queue does not repeat public boundary checks in `insert`, `update`, `erase`,
`contains`, or `pop`; its constructor still rejects invalid allocation
domains. FIFO ordering is preserved within each squared-distance bucket.

Production first checks a topology-only condition: whether the
inclusion-smallest owners of every A8 domain edge are comparable. When it
holds, no sibling can propagate into another before their join, so all domain
edges remain active and the hot kernel is the same unrestricted form used by
Opt3. Otherwise, each edge is activated at the LCA of its endpoint owners.
This dispatch depends on the hierarchy itself rather than max-tree, min-tree,
tree-of-shapes, residual, or generic type labels. The production policy also
compiles internal-operation validation out of the hot path; audited policies
remain available to tests.

See [Scientific benchmark builds](scientific-benchmark-builds.md) for the
benchmark entry points that isolate maximum-only, distributional, and spatial
scalar costs.
