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
    ${MMCFILTERS_ROOT}/filters/DepthStableRegionComputer.hpp
    ${MMCFILTERS_ROOT}/filters/ExtinctionValues.hpp
    ${MMCFILTERS_ROOT}/filters/MSERComputer.hpp
    ${MMCFILTERS_ROOT}/filters/UltimateAttributeOpening.hpp
    ${MMCFILTERS_ROOT}/localEvents/EventEngine.hpp
    ${MMCFILTERS_ROOT}/trees/HierarchySemantics.hpp
    ${MMCFILTERS_ROOT}/trees/MorphologicalTree.hpp
    ${MMCFILTERS_ROOT}/trees/MorphologicalTreeFactory.hpp
    ${MMCFILTERS_ROOT}/trees/NativeHierarchy.hpp
    ${MMCFILTERS_ROOT}/trees/ProperPartDomain.hpp
    ${MMCFILTERS_ROOT}/trees/TreeAltitudeAlgorithms.hpp
    ${MMCFILTERS_ROOT}/trees/TreeEditor.hpp
    ${MMCFILTERS_ROOT}/trees/TreeOfShapesProducer.hpp
    ${MMCFILTERS_ROOT}/trees/WeightedMorphologicalTree.hpp
    ${MMCFILTERS_ROOT}/trees/WeightedTreeView.hpp
    ${MMCFILTERS_ROOT}/trees/adjust/CasfComponentTrees.hpp
    ${MMCFILTERS_ROOT}/trees/sdrt/MinMaxResidualTreeBuilder.hpp
    ${MMCFILTERS_ROOT}/trees/sdrt/SdrtTiePolicy.hpp
    ${MMCFILTERS_ROOT}/utils/Altitude.hpp
    ${MMCFILTERS_ROOT}/utils/Common.hpp
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
    ${MMCFILTERS_ROOT}/attributes/computers/detail/BitquadAttributeMaterialization.hpp
    ${MMCFILTERS_ROOT}/attributes/computers/detail/BitquadLocalEventComputation.hpp
    ${MMCFILTERS_ROOT}/attributes/computers/detail/ContourSideAttributeData.hpp
    ${MMCFILTERS_ROOT}/attributes/computers/detail/ContourSideAttributeMaterialization.hpp
    ${MMCFILTERS_ROOT}/attributes/computers/detail/ContourSideLocalEventComputation.hpp
    ${MMCFILTERS_ROOT}/attributes/computers/detail/maxdist/EdtDIFT.hpp
    ${MMCFILTERS_ROOT}/attributes/computers/detail/maxdist/Geometry.hpp
    ${MMCFILTERS_ROOT}/attributes/computers/detail/maxdist/PQueue.hpp
    ${MMCFILTERS_ROOT}/attributes/detail/AttributeCapabilityValidation.hpp
    ${MMCFILTERS_ROOT}/attributes/detail/AttributeDeltaMaterialization.hpp
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
    ${MMCFILTERS_ROOT}/filters/DepthStableRegionComputer.hpp
    ${MMCFILTERS_ROOT}/filters/ExtinctionValues.hpp
    ${MMCFILTERS_ROOT}/filters/MSERComputer.hpp
    ${MMCFILTERS_ROOT}/filters/UltimateAttributeOpening.hpp
    ${MMCFILTERS_ROOT}/filters/detail/VariationMeasure.hpp
    ${MMCFILTERS_ROOT}/filters/detail/ViterbiDecision.hpp
    ${MMCFILTERS_ROOT}/localEvents/EventEngine.hpp
    ${MMCFILTERS_ROOT}/trees/HierarchySemantics.hpp
    ${MMCFILTERS_ROOT}/trees/MorphologicalTree.hpp
    ${MMCFILTERS_ROOT}/trees/MorphologicalTreeFactory.hpp
    ${MMCFILTERS_ROOT}/trees/NativeHierarchy.hpp
    ${MMCFILTERS_ROOT}/trees/ProperPartDomain.hpp
    ${MMCFILTERS_ROOT}/trees/TreeAltitudeAlgorithms.hpp
    ${MMCFILTERS_ROOT}/trees/TreeEditor.hpp
    ${MMCFILTERS_ROOT}/trees/TreeOfShapesProducer.hpp
    ${MMCFILTERS_ROOT}/trees/WeightedMorphologicalTree.hpp
    ${MMCFILTERS_ROOT}/trees/WeightedTreeView.hpp
    ${MMCFILTERS_ROOT}/trees/adjust/CasfComponentTrees.hpp
    ${MMCFILTERS_ROOT}/trees/adjust/DualMinMaxTreeIncrementalFilterLeaf.hpp
    ${MMCFILTERS_ROOT}/trees/adjust/DualMinMaxTreeIncrementalFilter.hpp
    ${MMCFILTERS_ROOT}/trees/adjust/DynamicTreeAttributeComputer.hpp
    ${MMCFILTERS_ROOT}/trees/detail/ComponentTreeProducerDetail.hpp
    ${MMCFILTERS_ROOT}/trees/detail/ComponentTreeUnionFind.hpp
    ${MMCFILTERS_ROOT}/trees/detail/HierarchyCapabilityValidation.hpp
    ${MMCFILTERS_ROOT}/trees/detail/HigraExportLayoutDetail.hpp
    ${MMCFILTERS_ROOT}/trees/detail/HigraHierarchyAdapterDetail.hpp
    ${MMCFILTERS_ROOT}/trees/detail/HigraImportLayoutDetail.hpp
    ${MMCFILTERS_ROOT}/trees/detail/MorphologicalTreeConstructionTag.hpp
    ${MMCFILTERS_ROOT}/trees/detail/NativeHierarchyValidationDetail.hpp
    ${MMCFILTERS_ROOT}/trees/detail/ProperPartEntryNode.hpp
    ${MMCFILTERS_ROOT}/trees/detail/TreeAltitudeDeltaNeighborhood.hpp
    ${MMCFILTERS_ROOT}/trees/detail/TreeStabilityNeighborhood.hpp
    ${MMCFILTERS_ROOT}/trees/detail/TreeTraversalDetail.hpp
    ${MMCFILTERS_ROOT}/trees/sdrt/MinMaxResidualTreeBuilder.hpp
    ${MMCFILTERS_ROOT}/trees/sdrt/SdrtTiePolicy.hpp
    ${MMCFILTERS_ROOT}/trees/sdrt/detail/UnionFindRegionTypes.hpp
    ${MMCFILTERS_ROOT}/trees/sdrt/detail/UnionFindResidualTreeAssembler.hpp
    ${MMCFILTERS_ROOT}/utils/Altitude.hpp
    ${MMCFILTERS_ROOT}/utils/Assert.hpp
    ${MMCFILTERS_ROOT}/utils/Common.hpp
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
