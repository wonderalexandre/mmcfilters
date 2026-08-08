"""Python bindings for MorphologicalAttributeFilters.

Preferred usage is the NodeId-based tree query API:
- `getRoot()`
- `getAliveNodeIds`
- `getChildren`
- `getProperParts`
- `getProperPartOwner`

Factory constructors return `WeightedMorphologicalTree`:
- `MorphologicalTreeFactory.createMaxTree(image)`
- `MorphologicalTreeFactory.createMinTree(image)`
- `MorphologicalTreeFactory.createSelfDualResidualTree(image)`
- `MorphologicalTreeFactory.createSaturatedSelfDualResidualTree(image, infinityPixel=0)`
- `MorphologicalTreeFactory.createTreeOfShapes(image)`
- `MorphologicalTreeFactory.createFromNativeTopology(node_parent, proper_part_owner, altitude, root, semantics)`
- `MorphologicalTreeFactory.createFromNativeTopology(node_parent, proper_part_owner, altitude, root, rows, cols, semantics)`
- `MorphologicalTreeFactory.createFromHigraParent(parent, altitude, rows, cols, MorphologicalTreeKind.MAX_TREE)`

Images passed to Python factories must be 2D C-contiguous `np.uint8` arrays.
Altitudes passed from Python must be integer sequences in `[0, 255]` or 1D
C-contiguous `np.uint8` arrays.

The C++ API has typed altitude owners intended for future Python expansion, but
the current Python API is deliberately limited to the canonical 8-bit image and
altitude contract. Python factories reject `np.int32`, `np.int64`,
`np.float32`, `np.float64`, `bool`, `object`, and non-contiguous image arrays
instead of converting them implicitly.

Python construction is centralized in `MorphologicalTreeFactory`. Factory and
import paths return `WeightedMorphologicalTree`, which also exposes the topology
queries listed above. A separate `MorphologicalTree` topology class is not
exported by the public Python API.

Weighted helpers live on `WeightedMorphologicalTree`:
- `weighted.getAltitude(node_id)`
- `weighted.setAltitude(node_id, value)`
- `weighted.setAltitudeBuffer(values)` / `weighted.altitude = values`
- `weighted.reconstructionImage()`

Hierarchy edge saliency maps live on `HierarchySaliencyMap`:
- `HierarchySaliencyMap.computeSaliencyEdgeMap(weighted, valuation, radius=None, strict=False, levelConvention=..., validateConnectivity=True)`
- `HierarchySaliencyMap.computeCanonicalRankedSaliencyEdgeMap(weighted, valuation)`
- `HierarchySaliencyMap.computeTopologicalLevelEdgeMap(weighted)`
- `HierarchySaliencyMap.computeNormalizedAltitudeEdgeMap(weighted)`

Component-tree completion is explicit on `ComponentTreePartitionHierarchyAdapter`:
- `ComponentTreePartitionHierarchyAdapter.validate(weighted, radius=None)`
- `ComponentTreePartitionHierarchyAdapter.computePartitionAppearanceLevels(weighted)`
- `ComponentTreePartitionHierarchyAdapter.computeSaliencyEdgeMap(weighted, levels, radius=None, strict=False)`

Hierarchy saliency valuation helpers live on `HierarchySaliencyMapValidation`:
- `HierarchySaliencyMapValidation.validateHierarchyConnectivity(weighted, radius=None)`
- `HierarchySaliencyMapValidation.validateHierarchyValuation(weighted, valuation, strict=False, nonnegative=False)`
- `HierarchySaliencyMapValidation.rankHierarchyValuation(weighted, valuation, strict=False)`
- `HierarchySaliencyMapValidation.computeNormalizedScores(weighted, valuation, strict=False, nonnegative=False)`

Hierarchy saliency projections live on `HierarchySaliencyMapProjection`:
- `HierarchySaliencyMapProjection.edgeMapToPixelImage(edge_map, reducer=EdgeToPixelReducer.Max)`
- `HierarchySaliencyMapProjection.thresholdCut(edge_map, threshold)`
- `HierarchySaliencyMapProjection.nodeContourEdges(weighted)`
- `HierarchySaliencyMapProjection.computeIncrementalNodeContours(weighted)`
- `HierarchySaliencyMapProjection.projectNodeValuation(contours, node_valuation)`
- `HierarchySaliencyMapProjection.thresholdByNodeValuation(contours, node_valuation, threshold)`

Xu shape-space extinction saliency lives on the separate `ShapeSpaceSaliency`
API:
- `ShapeSpaceSaliency.computeExtinctionValues(weighted, attribute, polarity)`
- `ShapeSpaceSaliency.projectContourScores(weighted, node_scores, radius=None)`
- `ShapeSpaceSaliency.compute(weighted, attribute, polarity, radius=None)`

Python exposes only the safe altitude setters. They validate max-tree/min-tree
altitude order before publishing changes; the C++ `Unchecked` setters and
editor sessions are intentionally not part of the Python API.

Attribute computation is centered on `Attribute`:
- `Attribute.computeSingleAttribute(weighted, attr)` and
  `Attribute.computeAttributes(weighted, attrs)` require a
  `WeightedMorphologicalTree`;
- `Attribute.computeSingleTopologyAttribute(tree_or_weighted, attr)` and
  `Attribute.computeTopologyAttributes(tree_or_weighted, attrs)` are the
  explicit topology/support-only entry points.

Depth-stability helpers expose topological depth windows through
`DepthStableRegionComputer.computeByDepth(depthDelta)` and report variation
scores with `getVariation`/`getVariations`. Result getters raise until the first
successful `computeByDepth` call.

Extinction-value helpers expose explicit selection policies, contour
visualization, and the formal edge-indexed saliency path:
- `ExtinctionSelectionPolicy.byTopK(extremaToKeep)`
- `ExtinctionSelectionPolicy.byThreshold(threshold)`
- `ExtinctionContourScorePolicy.RankScore`
- `ExtinctionContourScorePolicy.ExtinctionValue`
- `ExtinctionValues.getRegionalExtrema()`
- `ExtinctionValues.getExtinctionValueAttribute()`
- `ExtinctionValues.computeRankedExtinctionValueAttribute()`
- `ExtinctionValues.filtering(selection)`
- `ExtinctionValues.contourMap(selection, scorePolicy)`
- `ExtinctionValues.computeFormalSaliencyEdgeMap(radius=None, ranked=False)`
- `ExtinctionValues.computeMonotoneExtinctionProjection(radius=None, ranked=False)`
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
WeightedMorphologicalTree = _native.WeightedMorphologicalTree
HierarchySaliencyMap = _native.HierarchySaliencyMap
HierarchySaliencyMapValidation = _native.HierarchySaliencyMapValidation
HierarchySaliencyMapProjection = _native.HierarchySaliencyMapProjection
ComponentTreePartitionHierarchyAdapter = _native.ComponentTreePartitionHierarchyAdapter
ShapeSpaceExtremaPolarity = _native.ShapeSpaceExtremaPolarity
ShapeSpaceSaliency = _native.ShapeSpaceSaliency
EdgeToPixelReducer = _native.EdgeToPixelReducer
HierarchyLevelConvention = _native.HierarchyLevelConvention
ToSInterpolation = _native.ToSInterpolation
ToSPaddingPolicy = _native.ToSPaddingPolicy
TreeOfShapesProducerOptions = _native.TreeOfShapesProducerOptions
NodeIdSpace = _native.NodeIdSpace
MorphologicalTreeKind = _native.MorphologicalTreeKind
SdrtTiePolicy = _native.SdrtTiePolicy
AltitudeOrder = _native.AltitudeOrder
AdjacencyMode = _native.AdjacencyMode
GridDomain2D = _native.GridDomain2D
HierarchySemantics = _native.HierarchySemantics
RegularGridAdjacency2D = _native.RegularGridAdjacency2D
RegularGridAdjacencyShape = _native.RegularGridAdjacencyShape
DirectionalGridAdjacency2D = _native.DirectionalGridAdjacency2D
AttributeFilters = _native.AttributeFilters
DepthStableRegionComputer = _native.DepthStableRegionComputer
UltimateAttributeOpening = _native.UltimateAttributeOpening
Attribute = _native.Attribute
ContourComputation = _native.ContourComputation
Contours = _native.Contours
ContourRange = _native.ContourRange
ContoursIterator = _native.ContoursIterator
ContourTraceComputation = _native.ContourTraceComputation
ContourTraces = _native.ContourTraces
ContourTraceEdge = _native.ContourTraceEdge
ContourTraceLoop = _native.ContourTraceLoop
ContourTraceSide = _native.ContourTraceSide
ContourLoopKind = _native.ContourLoopKind
ExtinctionValues = _native.ExtinctionValues
ExtinctionSelectionPolicy = _native.ExtinctionSelectionPolicy
ExtinctionContourScorePolicy = _native.ExtinctionContourScorePolicy
DualMinMaxTreeIncrementalFilter = _native.DualMinMaxTreeIncrementalFilter
CasfComponentTrees = _native.CasfComponentTrees
CasfComponentTreesAttribute = _native.CasfComponentTreesAttribute

__all__ = [
    "__version__",
    "MorphologicalTreeFactory",
    "WeightedMorphologicalTree",
    "HierarchySaliencyMap",
    "HierarchySaliencyMapValidation",
    "HierarchySaliencyMapProjection",
    "ComponentTreePartitionHierarchyAdapter",
    "ShapeSpaceExtremaPolarity",
    "ShapeSpaceSaliency",
    "EdgeToPixelReducer",
    "HierarchyLevelConvention",
    "ToSInterpolation",
    "ToSPaddingPolicy",
    "TreeOfShapesProducerOptions",
    "NodeIdSpace",
    "MorphologicalTreeKind",
    "SdrtTiePolicy",
    "AltitudeOrder",
    "AdjacencyMode",
    "GridDomain2D",
    "HierarchySemantics",
    "RegularGridAdjacency2D",
    "RegularGridAdjacencyShape",
    "DirectionalGridAdjacency2D",
    "AttributeFilters",
    "DepthStableRegionComputer",
    "Attribute",
    "UltimateAttributeOpening",
    "ContourComputation",
    "Contours",
    "ContourRange",
    "ContoursIterator",
    "ContourTraceComputation",
    "ContourTraces",
    "ContourTraceEdge",
    "ContourTraceLoop",
    "ContourTraceSide",
    "ContourLoopKind",
    "ExtinctionValues",
    "ExtinctionSelectionPolicy",
    "ExtinctionContourScorePolicy",
    "DualMinMaxTreeIncrementalFilter",
    "CasfComponentTrees",
    "CasfComponentTreesAttribute",
]
