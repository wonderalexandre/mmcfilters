# Stable headers that users may include directly. Headers outside this list can
# still be installed below when a public header needs their template
# definitions, but those support headers are not independent API entry points.
set(MMCFILTERS_PUBLIC_HEADER_ENTRYPOINTS
    ${MMCFILTERS_ROOT}/attributes/AttributeComputation.hpp
    ${MMCFILTERS_ROOT}/attributes/AttributeNames.hpp
    ${MMCFILTERS_ROOT}/attributes/AttributeRegistry.hpp
    ${MMCFILTERS_ROOT}/attributes/AttributeResultTypes.hpp
    ${MMCFILTERS_ROOT}/attributes/AttributeTypes.hpp
    ${MMCFILTERS_ROOT}/attributes/Attributes.hpp
    ${MMCFILTERS_ROOT}/contours/ContourTraceComputation.hpp
    ${MMCFILTERS_ROOT}/contours/ContoursComputedIncrementally.hpp
    ${MMCFILTERS_ROOT}/filters/AttributeFilters.hpp
    ${MMCFILTERS_ROOT}/filters/AttributeReconstructionFilters.hpp
    ${MMCFILTERS_ROOT}/filters/NodeDecisionMasks.hpp
    ${MMCFILTERS_ROOT}/filters/NodePreservationStability.hpp
    ${MMCFILTERS_ROOT}/filters/DepthStableRegionComputer.hpp
    ${MMCFILTERS_ROOT}/filters/ExtinctionValues.hpp
    ${MMCFILTERS_ROOT}/filters/MSERComputer.hpp
    ${MMCFILTERS_ROOT}/filters/UltimateAttributeOpening.hpp
    ${MMCFILTERS_ROOT}/localAttributes/FiniteWindowLocalAttributeComputer.hpp
    ${MMCFILTERS_ROOT}/trees/MorphologicalTreeSemantics.hpp
    ${MMCFILTERS_ROOT}/trees/MorphologicalTree.hpp
    ${MMCFILTERS_ROOT}/trees/MorphologicalTreeFactory.hpp
    ${MMCFILTERS_ROOT}/trees/NativeHierarchy.hpp
    ${MMCFILTERS_ROOT}/trees/ProperPartDomain.hpp
    ${MMCFILTERS_ROOT}/trees/TreeAltitudeAlgorithms.hpp
    ${MMCFILTERS_ROOT}/trees/TreeEditor.hpp
    ${MMCFILTERS_ROOT}/trees/TreeOfShapesProducer.hpp
    ${MMCFILTERS_ROOT}/trees/ValuedMorphologicalTree.hpp
    ${MMCFILTERS_ROOT}/trees/ValuedMorphologicalTreeView.hpp
    ${MMCFILTERS_ROOT}/trees/adjust/CasfComponentTrees.hpp
    ${MMCFILTERS_ROOT}/trees/saliency/HierarchySaliencyMap.hpp
    ${MMCFILTERS_ROOT}/trees/saliency/HierarchySaliencyMapProjection.hpp
    ${MMCFILTERS_ROOT}/trees/saliency/HierarchySaliencyMapValidation.hpp
    ${MMCFILTERS_ROOT}/trees/saliency/HierarchicalWatershedSaliency.hpp
    ${MMCFILTERS_ROOT}/trees/saliency/ShapeSpaceSaliency.hpp
    ${MMCFILTERS_ROOT}/trees/sdrt/SaturatedResidualTreeBuilder.hpp
    ${MMCFILTERS_ROOT}/trees/sdrt/ResidualTreeBuildStatistics.hpp
    ${MMCFILTERS_ROOT}/trees/sdrt/ResidualEvolution.hpp
    ${MMCFILTERS_ROOT}/trees/sdrt/ResidualTreePolicies.hpp
    ${MMCFILTERS_ROOT}/trees/sdrt/UnrestrictedResidualTreeBuilder.hpp
    ${MMCFILTERS_ROOT}/utils/Altitude.hpp
    ${MMCFILTERS_ROOT}/utils/Common.hpp
    ${MMCFILTERS_ROOT}/utils/Contract.hpp
    ${MMCFILTERS_ROOT}/utils/Image.hpp
    ${MMCFILTERS_ROOT}/utils/RegularGridAdjacency2D.hpp
)

# Explicit installation manifest. Besides the entry points above, this contains
# only the implementation support reached by their header-only definitions.
# Unreferenced internal helpers are deliberately absent.
set(MMCFILTERS_INSTALL_HEADERS
    ${MMCFILTERS_ROOT}/attributes/AttributeComputation.hpp
    ${MMCFILTERS_ROOT}/attributes/AttributeNames.hpp
    ${MMCFILTERS_ROOT}/attributes/AttributeRegistry.hpp
    ${MMCFILTERS_ROOT}/attributes/AttributeResultTypes.hpp
    ${MMCFILTERS_ROOT}/attributes/AttributeTypes.hpp
    ${MMCFILTERS_ROOT}/attributes/Attributes.hpp
    ${MMCFILTERS_ROOT}/attributes/computers/AreaComputer.hpp
    ${MMCFILTERS_ROOT}/attributes/computers/AttributeComputerDomain.hpp
    ${MMCFILTERS_ROOT}/attributes/computers/AttributeComputerFamily.hpp
    ${MMCFILTERS_ROOT}/attributes/computers/AttributeComputerRegistry.hpp
    ${MMCFILTERS_ROOT}/attributes/computers/BitquadAttributeComputer.hpp
    ${MMCFILTERS_ROOT}/attributes/computers/BoundingBoxComputer.hpp
    ${MMCFILTERS_ROOT}/attributes/computers/ContourSideAttributeComputer.hpp
    ${MMCFILTERS_ROOT}/attributes/computers/GrayLevelStatsComputer.hpp
    ${MMCFILTERS_ROOT}/attributes/computers/MaxDistComputer.hpp
    ${MMCFILTERS_ROOT}/attributes/computers/MomentBasedAttributeComputer.hpp
    ${MMCFILTERS_ROOT}/attributes/computers/TreeTopologyComputer.hpp
    ${MMCFILTERS_ROOT}/attributes/computers/VolumeComputer.hpp
    ${MMCFILTERS_ROOT}/attributes/computers/detail/BitquadAttributeData.hpp
    ${MMCFILTERS_ROOT}/attributes/computers/detail/BitquadConnectivityPolicy.hpp
    ${MMCFILTERS_ROOT}/attributes/computers/detail/BitquadAttributeProjection.hpp
    ${MMCFILTERS_ROOT}/attributes/computers/detail/BitquadFiniteWindowComputation.hpp
    ${MMCFILTERS_ROOT}/attributes/computers/detail/ContourSideAttributeData.hpp
    ${MMCFILTERS_ROOT}/attributes/computers/detail/ContourSideAttributeMaterialization.hpp
    ${MMCFILTERS_ROOT}/attributes/computers/detail/ContourSideFiniteWindowComputation.hpp
    ${MMCFILTERS_ROOT}/attributes/computers/detail/maxdist/EdtDIFT.hpp
    ${MMCFILTERS_ROOT}/attributes/computers/detail/maxdist/Geometry.hpp
    ${MMCFILTERS_ROOT}/attributes/computers/detail/maxdist/PQueue.hpp
    ${MMCFILTERS_ROOT}/attributes/detail/AttributeCapabilityValidation.hpp
    ${MMCFILTERS_ROOT}/attributes/detail/NodeAttributeSampleMaterialization.hpp
    ${MMCFILTERS_ROOT}/attributes/detail/AttributeDependencyCache.hpp
    ${MMCFILTERS_ROOT}/attributes/detail/AttributeFamilyScheduler.hpp
    ${MMCFILTERS_ROOT}/attributes/detail/AttributeKernelSupport.hpp
    ${MMCFILTERS_ROOT}/attributes/detail/AttributeNumericPolicy.hpp
    ${MMCFILTERS_ROOT}/attributes/detail/AttributePipeline.hpp
    ${MMCFILTERS_ROOT}/attributes/detail/AttributeProjection.hpp
    ${MMCFILTERS_ROOT}/attributes/detail/AttributeRequestUtils.hpp
    ${MMCFILTERS_ROOT}/attributes/detail/TopologyAttributeBackend.hpp
    ${MMCFILTERS_ROOT}/contours/ContourTraceComputation.hpp
    ${MMCFILTERS_ROOT}/contours/ContoursComputedIncrementally.hpp
    ${MMCFILTERS_ROOT}/contours/detail/ContourDeltaStore.hpp
    ${MMCFILTERS_ROOT}/contours/detail/ContourTraceDeltaStore.hpp
    ${MMCFILTERS_ROOT}/contours/detail/PendingPixelLists.hpp
    ${MMCFILTERS_ROOT}/dataStructure/FastQueue.hpp
    ${MMCFILTERS_ROOT}/filters/AttributeFilters.hpp
    ${MMCFILTERS_ROOT}/filters/AttributeReconstructionFilters.hpp
    ${MMCFILTERS_ROOT}/filters/NodeDecisionMasks.hpp
    ${MMCFILTERS_ROOT}/filters/NodePreservationStability.hpp
    ${MMCFILTERS_ROOT}/filters/DepthStableRegionComputer.hpp
    ${MMCFILTERS_ROOT}/filters/ExtinctionValues.hpp
    ${MMCFILTERS_ROOT}/filters/MSERComputer.hpp
    ${MMCFILTERS_ROOT}/filters/UltimateAttributeOpening.hpp
    ${MMCFILTERS_ROOT}/filters/detail/VariationMeasure.hpp
    ${MMCFILTERS_ROOT}/filters/detail/ViterbiDecision.hpp
    ${MMCFILTERS_ROOT}/localAttributes/FiniteWindowLocalAttributeComputer.hpp
    ${MMCFILTERS_ROOT}/trees/MorphologicalTreeSemantics.hpp
    ${MMCFILTERS_ROOT}/trees/MorphologicalTree.hpp
    ${MMCFILTERS_ROOT}/trees/MorphologicalTreeFactory.hpp
    ${MMCFILTERS_ROOT}/trees/NativeHierarchy.hpp
    ${MMCFILTERS_ROOT}/trees/ProperPartDomain.hpp
    ${MMCFILTERS_ROOT}/trees/TreeAltitudeAlgorithms.hpp
    ${MMCFILTERS_ROOT}/trees/TreeEditor.hpp
    ${MMCFILTERS_ROOT}/trees/TreeOfShapesProducer.hpp
    ${MMCFILTERS_ROOT}/trees/ValuedMorphologicalTree.hpp
    ${MMCFILTERS_ROOT}/trees/ValuedMorphologicalTreeView.hpp
    ${MMCFILTERS_ROOT}/trees/adjust/CasfComponentTrees.hpp
    ${MMCFILTERS_ROOT}/trees/adjust/DualMinMaxTreeIncrementalFilterLeaf.hpp
    ${MMCFILTERS_ROOT}/trees/adjust/DualMinMaxTreeIncrementalFilter.hpp
    ${MMCFILTERS_ROOT}/trees/adjust/DynamicTreeAttributeComputer.hpp
    ${MMCFILTERS_ROOT}/trees/detail/ComponentTreeProducerDetail.hpp
    ${MMCFILTERS_ROOT}/trees/detail/CommittedTreeAccess.hpp
    ${MMCFILTERS_ROOT}/trees/detail/MorphologicalTreeConstructionContextQueries.hpp
    ${MMCFILTERS_ROOT}/trees/detail/ComponentTreeUnionFind.hpp
    ${MMCFILTERS_ROOT}/trees/detail/HierarchyCapabilityValidation.hpp
    ${MMCFILTERS_ROOT}/trees/detail/HigraExportLayoutDetail.hpp
    ${MMCFILTERS_ROOT}/trees/detail/HigraHierarchyAdapterDetail.hpp
    ${MMCFILTERS_ROOT}/trees/detail/HigraImportLayoutDetail.hpp
    ${MMCFILTERS_ROOT}/trees/detail/MorphologicalTreeConstructionTag.hpp
    ${MMCFILTERS_ROOT}/trees/detail/NativeHierarchyValidationDetail.hpp
    ${MMCFILTERS_ROOT}/trees/detail/TreeAttributeSamplingNeighborhood.hpp
    ${MMCFILTERS_ROOT}/trees/detail/TreeStabilityNeighborhood.hpp
    ${MMCFILTERS_ROOT}/trees/detail/TreeTraversalDetail.hpp
    ${MMCFILTERS_ROOT}/trees/saliency/HierarchySaliencyMap.hpp
    ${MMCFILTERS_ROOT}/trees/saliency/HierarchySaliencyMapProjection.hpp
    ${MMCFILTERS_ROOT}/trees/saliency/HierarchySaliencyMapValidation.hpp
    ${MMCFILTERS_ROOT}/trees/saliency/HierarchicalWatershedSaliency.hpp
    ${MMCFILTERS_ROOT}/trees/saliency/ShapeSpaceSaliency.hpp
    ${MMCFILTERS_ROOT}/trees/sdrt/ResidualTreeBuildStatistics.hpp
    ${MMCFILTERS_ROOT}/trees/sdrt/ResidualEvolution.hpp
    ${MMCFILTERS_ROOT}/trees/sdrt/ResidualTreePolicies.hpp
    ${MMCFILTERS_ROOT}/trees/sdrt/SaturatedResidualTreeBuilder.hpp
    ${MMCFILTERS_ROOT}/trees/sdrt/UnrestrictedResidualTreeBuilder.hpp
    ${MMCFILTERS_ROOT}/trees/sdrt/detail/SynchronizedResidualTreeEvolution.hpp
    ${MMCFILTERS_ROOT}/trees/sdrt/detail/FlatZonePartition.hpp
    ${MMCFILTERS_ROOT}/trees/sdrt/detail/ResidualTreeCandidateAgenda.hpp
    ${MMCFILTERS_ROOT}/trees/sdrt/detail/ResidualTreeCandidateContext.hpp
    ${MMCFILTERS_ROOT}/trees/sdrt/detail/ResidualTreeCandidatePreparation.hpp
    ${MMCFILTERS_ROOT}/trees/sdrt/detail/ResidualTreeCandidateTypes.hpp
    ${MMCFILTERS_ROOT}/trees/sdrt/detail/ResidualTreeEventAssembler.hpp
    ${MMCFILTERS_ROOT}/trees/sdrt/detail/ResidualTreeMaterialization.hpp
    ${MMCFILTERS_ROOT}/trees/sdrt/detail/ResidualTreeRegionTypes.hpp
    ${MMCFILTERS_ROOT}/trees/sdrt/detail/SaturatedDynamicLca.hpp
    ${MMCFILTERS_ROOT}/trees/sdrt/detail/SaturatedResidualEligibility.hpp
    ${MMCFILTERS_ROOT}/utils/Altitude.hpp
    ${MMCFILTERS_ROOT}/utils/Assert.hpp
    ${MMCFILTERS_ROOT}/utils/Common.hpp
    ${MMCFILTERS_ROOT}/utils/CommittedGridAccess.hpp
    ${MMCFILTERS_ROOT}/utils/CommittedImageAccess.hpp
    ${MMCFILTERS_ROOT}/utils/Contract.hpp
    ${MMCFILTERS_ROOT}/utils/GenerationStampSet.hpp
    ${MMCFILTERS_ROOT}/utils/Image.hpp
    ${MMCFILTERS_ROOT}/utils/RegularGridAdjacency2D.hpp
)

list(REMOVE_DUPLICATES MMCFILTERS_PUBLIC_HEADER_ENTRYPOINTS)
list(REMOVE_DUPLICATES MMCFILTERS_INSTALL_HEADERS)

set(MMCFILTERS_PUBLIC_HEADER_SUPPORT ${MMCFILTERS_INSTALL_HEADERS})
list(REMOVE_ITEM
    MMCFILTERS_PUBLIC_HEADER_SUPPORT
    ${MMCFILTERS_PUBLIC_HEADER_ENTRYPOINTS})

foreach(public_header IN LISTS MMCFILTERS_PUBLIC_HEADER_ENTRYPOINTS)
    if(NOT public_header IN_LIST MMCFILTERS_INSTALL_HEADERS)
        message(FATAL_ERROR
            "Public header is absent from the installation manifest: "
            "${public_header}")
    endif()
endforeach()

foreach(installed_header IN LISTS MMCFILTERS_INSTALL_HEADERS)
    if(NOT EXISTS "${installed_header}")
        message(FATAL_ERROR
            "Installed-header manifest contains a missing file: "
            "${installed_header}")
    endif()
endforeach()
