"""Python bindings for mmcfilters.

Preferred usage is the NodeId-based tree query API:
- `root`
- `alive_node_ids`
- `children`
- `proper_part`
- `node_support`
- `smallest_node`

Factory constructors return `ValuedMorphologicalTree`:
- `MorphologicalTreeFactory.create_max_tree(image)`
- `MorphologicalTreeFactory.create_min_tree(image)`
- `MorphologicalTreeFactory.create_unrestricted_residual_tree(image)`
- `MorphologicalTreeFactory.create_saturated_residual_tree(image, infinity_pixel=0)`
- `MorphologicalTreeFactory.create_tree_of_shapes(image)`
- `MorphologicalTreeFactory.create_from_native_topology(parent, smallest_node_map, node_altitudes, root, semantics)`
- `MorphologicalTreeFactory.create_from_native_topology(parent, smallest_node_map, node_altitudes, root, rows, columns, semantics)`
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

Direct and subtractive reconstruction are separate filter families:
- `DirectAttributeFilter(tree).apply_direct_attribute_filter(NodePreservationMask(...))`
- `SubtractiveAttributeFilter(tree).apply_subtractive_attribute_filter(NodePreservationMask(...))`
- `SoftSubtractiveAttributeFilter(tree).apply_soft_subtractive_attribute_filter(node_preservation_scores)`

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

Hierarchy saliency projections live on `HierarchySaliencyMapProjection`:
- `HierarchySaliencyMapProjection.edge_map_to_pixel_image(edge_map, reducer=EdgeToPixelReducer.MAX)`
- `HierarchySaliencyMapProjection.threshold_cut(edge_map, threshold)`
- `HierarchySaliencyMapProjection.node_contour_edges(valued_tree)`
- `HierarchySaliencyMapProjection.compute_incremental_node_contours(valued_tree)`
- `HierarchySaliencyMapProjection.project_node_valuation(contours, node_valuation)`
- `HierarchySaliencyMapProjection.threshold_by_node_valuation(contours, node_valuation, threshold)`

Xu shape-space extinction saliency lives on the separate `ShapeSpaceSaliency`
API:
- `ShapeSpaceSaliency.compute_extinction_values(valued_tree, attribute, polarity)`
- `ShapeSpaceSaliency.project_contour_scores(valued_tree, node_scores, radius=None)`
- `ShapeSpaceSaliency.compute(valued_tree, attribute, polarity, radius=None)`

Python exposes only the safe altitude setters. They validate max-tree/min-tree
altitude order before publishing changes; the C++ `Unchecked` setters and
editor sessions are intentionally not part of the Python API.

Attribute computation is centered on `Attribute`:
- `Attribute.compute_single_attribute(valued_tree, attr)` and
  `Attribute.compute_attributes(valued_tree, attrs)` require a
  `ValuedMorphologicalTree`;
- `Attribute.compute_single_topology_attribute(tree_or_valued_tree, attr)` and
  `Attribute.compute_topology_attributes(tree_or_valued_tree, attrs)` are the
  explicit topology/support-only entry points;
- `Attribute.compute_sampled_node_attribute(...)` uses
  `NodeAttributeSamplingPolicy` and `MissingNodeAttributeSamplePolicy` for
  altitude-based ancestor/current/representative-descendant samples.

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
EdgeToPixelReducer = _native.EdgeToPixelReducer
HierarchyLevelConvention = _native.HierarchyLevelConvention
NodeIdSpace = _native.NodeIdSpace
MorphologicalTreeKind = _native.MorphologicalTreeKind
Polarity = _native.Polarity
SpatialOrder = _native.SpatialOrder
RowMajorSpatialOrder = _native.RowMajorSpatialOrder
SelfDualResidualKey = _native.SelfDualResidualKey
SelfDualResidualOrder = _native.SelfDualResidualOrder
SelfDualResidualSchedule = _native.SelfDualResidualSchedule
NodeAltitudeOrder = _native.NodeAltitudeOrder
GridDomain2D = _native.GridDomain2D
MorphologicalTreeSemantics = _native.MorphologicalTreeSemantics
NoConstructionContext = _native.NoConstructionContext
SharedAdjacencyContext = _native.SharedAdjacencyContext
SaturatedResidualContext = _native.SaturatedResidualContext
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
NodePreservationMask = _native.NodePreservationMask
NodePruningMask = _native.NodePruningMask
to_node_pruning_mask = _native.to_node_pruning_mask
to_node_preservation_mask = _native.to_node_preservation_mask
IncompleteStabilityWindowPolicy = _native.IncompleteStabilityWindowPolicy
compute_node_preservation_mask = _native.compute_node_preservation_mask
adjust_node_preservation_mask_by_altitude_stability = _native.adjust_node_preservation_mask_by_altitude_stability
adjust_node_preservation_mask_by_depth_stability = _native.adjust_node_preservation_mask_by_depth_stability
DirectReconstruction = _native.DirectReconstruction
SubtractiveResidueModulation = _native.SubtractiveResidueModulation
DirectAttributeFilter = _native.DirectAttributeFilter
SubtractiveAttributeFilter = _native.SubtractiveAttributeFilter
SoftSubtractiveAttributeFilter = _native.SoftSubtractiveAttributeFilter
DepthStableRegionComputer = _native.DepthStableRegionComputer
UltimateAttributeOpening = _native.UltimateAttributeOpening
Attribute = _native.Attribute
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
    "EdgeToPixelReducer",
    "HierarchyLevelConvention",
    "NodeIdSpace",
    "MorphologicalTreeKind",
    "Polarity",
    "SpatialOrder",
    "RowMajorSpatialOrder",
    "SelfDualResidualKey",
    "SelfDualResidualOrder",
    "SelfDualResidualSchedule",
    "NodeAltitudeOrder",
    "GridDomain2D",
    "MorphologicalTreeSemantics",
    "NoConstructionContext",
    "SharedAdjacencyContext",
    "SaturatedResidualContext",
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
    "NodePreservationMask",
    "NodePruningMask",
    "to_node_pruning_mask",
    "to_node_preservation_mask",
    "IncompleteStabilityWindowPolicy",
    "compute_node_preservation_mask",
    "adjust_node_preservation_mask_by_altitude_stability",
    "adjust_node_preservation_mask_by_depth_stability",
    "DirectReconstruction",
    "SubtractiveResidueModulation",
    "DirectAttributeFilter",
    "SubtractiveAttributeFilter",
    "SoftSubtractiveAttributeFilter",
    "DepthStableRegionComputer",
    "Attribute",
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
