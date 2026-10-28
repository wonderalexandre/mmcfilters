"""Python bindings for mmcfilters.

Preferred usage is the NodeId-based tree query API:
- `root`
- `alive_node_ids`
- `children`
- `parents`
- `proper_part`
- `node_support`
- `smallest_node`

Factory constructors return `ValuedMorphologicalTree`:
- `MorphologicalTreeFactory.create_max_tree(image)`
- `MorphologicalTreeFactory.create_min_tree(image)`
- `MorphologicalTreeFactory.create_extremum_residual_tree(image)`
- `MorphologicalTreeFactory.create_tree_of_shapes(image)`
- `MorphologicalTreeFactory.create_from_higra_parent(parent, node_altitudes, rows, columns, MorphologicalTreeKind.MAX_TREE)`

Images passed to Python factories must be 2D C-contiguous `np.uint8` arrays.
Altitudes passed from Python must be integer sequences in `[0, 255]` or 1D
C-contiguous `np.uint8` arrays.

The C++ API has typed altitude owners intended for future Python expansion, but
the current Python API is deliberately limited to the canonical 8-bit image and
altitude contract. Python factories reject `np.int32`, `np.int64`,
`np.float32`, `np.float64`, `bool`, `object`, and non-contiguous image arrays
instead of converting them implicitly.

Python construction is centralized in `MorphologicalTreeFactory`. Factory and
import paths return `ValuedMorphologicalTree`, which also exposes the topology
queries listed above. A separate `MorphologicalTree` topology class is not
exported by the public Python API.

Valued-tree helpers live on `ValuedMorphologicalTree`:
- `valued_tree.node_altitude(node_id)`
- `valued_tree.set_node_altitude(node_id, value)`
- `valued_tree.node_altitudes = values`
- `valued_tree.reconstruct_from_node_altitudes()`
- `valued_tree.reconstruct_from_node_contributions(node_contributions)`

Filtering rules are available on `AttributeFilters(tree)`:
- `AttributeFilters(tree).filtering_by_direct_rule(mask)`
- `AttributeFilters(tree).filtering_by_subtractive_rule(mask)`
- `AttributeFilters(tree).filtering_by_soft_subtractive_rule(node_preservation_scores)`

Hierarchy edge saliency maps live on `HierarchySaliencyMap`:
- `HierarchySaliencyMap.compute_saliency_edge_map(valued_tree, valuation, radius=None, strict=False, level_convention=..., validate_connectivity=True)`
- `HierarchySaliencyMap.compute_canonical_ranked_saliency_edge_map(valued_tree, valuation)`
- `HierarchySaliencyMap.compute_topological_level_edge_map(valued_tree)`
- `HierarchySaliencyMap.compute_normalized_altitude_edge_map(valued_tree)`

Component-tree completion is explicit on `ComponentTreePartitionHierarchyAdapter`:
- `ComponentTreePartitionHierarchyAdapter.validate(valued_tree, radius=None)`
- `ComponentTreePartitionHierarchyAdapter.compute_partition_appearance_levels(valued_tree)`
- `ComponentTreePartitionHierarchyAdapter.compute_saliency_edge_map(valued_tree, levels, radius=None, strict=False)`

Hierarchy saliency valuation helpers live on `HierarchySaliencyMapValidation`:
- `HierarchySaliencyMapValidation.validate_hierarchy_connectivity(valued_tree, radius=None)`
- `HierarchySaliencyMapValidation.validate_hierarchy_valuation(valued_tree, valuation, strict=False, nonnegative=False)`
- `HierarchySaliencyMapValidation.rank_hierarchy_valuation(valued_tree, valuation, strict=False)`
- `HierarchySaliencyMapValidation.compute_normalized_scores(valued_tree, valuation, strict=False, nonnegative=False)`

Edge-map operations and hierarchy saliency projections:
- `EdgeMap.to_pixel_image(reducer=EdgeToPixelReducer.MAX)`
- `EdgeMap.threshold(threshold)` returns `EdgeContours`
- `HierarchySaliencyMapProjection.node_contour_edges(valued_tree)`
- `HierarchySaliencyMapProjection.compute_incremental_node_contours(valued_tree)`
- `IncrementalNodeContours.project_node_valuation(node_valuation)`
- `IncrementalNodeContours.threshold_by_node_valuation(node_valuation, threshold)`

Xu shape-space extinction saliency lives on the separate `ShapeSpaceSaliency`
API:
- `ShapeSpaceSaliency.compute_extinction_values(valued_tree, attribute, polarity)`
- `ShapeSpaceSaliency.project_contour_scores(valued_tree, node_scores, radius=None)`
- `ShapeSpaceSaliency.compute(valued_tree, attribute, polarity, radius=None)`

Python exposes only the safe altitude setters. They validate max-tree/min-tree
altitude order before publishing changes; the C++ `Unchecked` setters and
editor sessions are intentionally not part of the Python API.

Attribute computation is centered on `Attribute`:
- `Attribute.compute(tree, attributes, output_space=NodeIdSpace.MORPHOLOGICAL_TREE,
  dtype=np.float32)` selects the topology or altitude-aware path from the
  registry requirements. A scalar value or name returns a 1D array; a group or
  iterable returns `(layout, values)` with a 2D array;
- `Attribute.compute_sampled_node_attribute(...)` uses
  `NodeAttributeSamplingPolicy` and `MissingNodeAttributeSamplePolicy` for
  altitude-based ancestor/current/representative-descendant samples.

Distance-transform fields and profiles are explicit topology-only results:
- `DistanceTransform.compute_node_field(valued_tree, node_id, accuracy="exact")`
- `DistanceTransform.compute_profiles(valued_tree, accuracy="exact")`

Depth-stability helpers expose topological depth windows through
`DepthStableRegionComputer.compute_by_depth(depth_window_radius)` and report variation
scores with `get_variation`/`get_variations`. Result getters raise until the first
successful `compute_by_depth` call.

Extinction-value helpers expose explicit selection policies, contour
visualization, and the formal edge-indexed saliency path:
- `ExtinctionSelectionPolicy.by_top_k(extrema_to_keep)`
- `ExtinctionSelectionPolicy.by_threshold(threshold)`
- `ExtinctionContourScorePolicy.RANK_SCORE`
- `ExtinctionContourScorePolicy.EXTINCTION_VALUE`
- `ExtinctionValues.get_regional_extrema()`
- `ExtinctionValues.get_extinction_value_attribute()`
- `ExtinctionValues.compute_ranked_extinction_value_attribute()`
- `ExtinctionValues.filtering(selection)`
- `ExtinctionValues.contour_map(selection, score_policy)`
- `ExtinctionValues.compute_formal_saliency_edge_map(radius=None, ranked=False)`
- `ExtinctionValues.compute_monotone_extinction_projection(radius=None, ranked=False)`
"""

# Version (generated by setuptools-scm at build time when packaging)
try:
    from ._version import version as __version__
except ImportError:  # pragma: no cover - local CMake builds may not generate it
    __version__ = "0+unknown"

# Import native pybind11 module packaged under this package
from . import mmcfilters as _native

# Re-export native symbols at package level
MorphologicalTreeFactory = _native.MorphologicalTreeFactory
ValuedMorphologicalTree = _native.ValuedMorphologicalTree
HierarchySaliencyMap = _native.HierarchySaliencyMap
HierarchySaliencyMapValidation = _native.HierarchySaliencyMapValidation
HierarchySaliencyMapProjection = _native.HierarchySaliencyMapProjection
ComponentTreePartitionHierarchyAdapter = _native.ComponentTreePartitionHierarchyAdapter
ShapeSpaceExtremaPolarity = _native.ShapeSpaceExtremaPolarity
ShapeSpaceSaliency = _native.ShapeSpaceSaliency
EdgeMap = _native.EdgeMap
EdgeContours = _native.EdgeContours
NodeContourEdges = _native.NodeContourEdges
IncrementalNodeContours = _native.IncrementalNodeContours
ShapeSpaceExtremum = _native.ShapeSpaceExtremum
ShapeSpaceExtinctionResult = _native.ShapeSpaceExtinctionResult
ShapeSpaceSaliencyResult = _native.ShapeSpaceSaliencyResult
EdgeToPixelReducer = _native.EdgeToPixelReducer
HierarchyLevelConvention = _native.HierarchyLevelConvention
NodeIdSpace = _native.NodeIdSpace
MorphologicalTreeKind = _native.MorphologicalTreeKind
NodeAltitudeOrder = _native.NodeAltitudeOrder
GridDomain2D = _native.GridDomain2D
MorphologicalTreeSemantics = _native.MorphologicalTreeSemantics
NoConstructionContext = _native.NoConstructionContext
SharedAdjacencyContext = _native.SharedAdjacencyContext
ComplementaryAdjacencies = _native.ComplementaryAdjacencies
SelfDualSpanImmersion = _native.SelfDualSpanImmersion
ComplementaryPairing = _native.ComplementaryPairing
CanonicalComplementaryGridImmersion = _native.CanonicalComplementaryGridImmersion
ComplementaryGridImmersion = _native.ComplementaryGridImmersion
TopographicDomainExtension = _native.TopographicDomainExtension
TopographicAltitudeEncoding = _native.TopographicAltitudeEncoding
TopographicConvention = _native.TopographicConvention
self_dual_span_convention = _native.self_dual_span_convention
RegularGridAdjacency2D = _native.RegularGridAdjacency2D
RegularGridAdjacencyShape = _native.RegularGridAdjacencyShape
AttributeFilters = _native.AttributeFilters
IncompleteStabilityWindowPolicy = _native.IncompleteStabilityWindowPolicy
adjust_node_preservation_mask_by_altitude_stability = _native.adjust_node_preservation_mask_by_altitude_stability
adjust_node_preservation_mask_by_depth_stability = _native.adjust_node_preservation_mask_by_depth_stability
DirectReconstruction = _native.DirectReconstruction
SubtractiveResidueModulation = _native.SubtractiveResidueModulation
DepthStableRegionComputer = _native.DepthStableRegionComputer
UltimateAttributeOpening = _native.UltimateAttributeOpening
Attribute = _native.Attribute
DistanceTransformAccuracy = _native.DistanceTransformAccuracy
NodeDistanceField = _native.NodeDistanceField
DistanceTransformProfiles = _native.DistanceTransformProfiles
DistanceTransform = _native.DistanceTransform
NodeAttributeSamplingPolicy = _native.NodeAttributeSamplingPolicy
MissingNodeAttributeSamplePolicy = _native.MissingNodeAttributeSamplePolicy
ContourComputation = _native.ContourComputation
ContourTraceComputation = _native.ContourTraceComputation
ContourTrace = _native.ContourTrace
ContourEdge = _native.ContourEdge
ContourBoundary = _native.ContourBoundary
ContourSide = _native.ContourSide
ContourBoundaryKind = _native.ContourBoundaryKind
ExtinctionValues = _native.ExtinctionValues
ExtinctionSelectionPolicy = _native.ExtinctionSelectionPolicy
ExtinctionContourScorePolicy = _native.ExtinctionContourScorePolicy
DualMinMaxTreeIncrementalFilter = _native.DualMinMaxTreeIncrementalFilter
CasfComponentTrees = _native.CasfComponentTrees
CasfComponentTreesAttribute = _native.CasfComponentTreesAttribute

__all__ = [
    "__version__",
    "MorphologicalTreeFactory",
    "ValuedMorphologicalTree",
    "HierarchySaliencyMap",
    "HierarchySaliencyMapValidation",
    "HierarchySaliencyMapProjection",
    "ComponentTreePartitionHierarchyAdapter",
    "ShapeSpaceExtremaPolarity",
    "ShapeSpaceSaliency",
    "EdgeMap",
    "EdgeContours",
    "NodeContourEdges",
    "IncrementalNodeContours",
    "ShapeSpaceExtremum",
    "ShapeSpaceExtinctionResult",
    "ShapeSpaceSaliencyResult",
    "EdgeToPixelReducer",
    "HierarchyLevelConvention",
    "NodeIdSpace",
    "MorphologicalTreeKind",
    "NodeAltitudeOrder",
    "GridDomain2D",
    "MorphologicalTreeSemantics",
    "NoConstructionContext",
    "SharedAdjacencyContext",
    "ComplementaryAdjacencies",
    "SelfDualSpanImmersion",
    "ComplementaryPairing",
    "CanonicalComplementaryGridImmersion",
    "ComplementaryGridImmersion",
    "TopographicDomainExtension",
    "TopographicAltitudeEncoding",
    "TopographicConvention",
    "self_dual_span_convention",
    "RegularGridAdjacency2D",
    "RegularGridAdjacencyShape",
    "AttributeFilters",
    "IncompleteStabilityWindowPolicy",
    "adjust_node_preservation_mask_by_altitude_stability",
    "adjust_node_preservation_mask_by_depth_stability",
    "DirectReconstruction",
    "SubtractiveResidueModulation",
    "DepthStableRegionComputer",
    "Attribute",
    "DistanceTransformAccuracy",
    "NodeDistanceField",
    "DistanceTransformProfiles",
    "DistanceTransform",
    "NodeAttributeSamplingPolicy",
    "MissingNodeAttributeSamplePolicy",
    "UltimateAttributeOpening",
    "ContourComputation",
    "ContourTraceComputation",
    "ContourTrace",
    "ContourEdge",
    "ContourBoundary",
    "ContourSide",
    "ContourBoundaryKind",
    "ExtinctionValues",
    "ExtinctionSelectionPolicy",
    "ExtinctionContourScorePolicy",
    "DualMinMaxTreeIncrementalFilter",
    "CasfComponentTrees",
    "CasfComponentTreesAttribute",
]
