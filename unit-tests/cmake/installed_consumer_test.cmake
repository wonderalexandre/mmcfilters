cmake_minimum_required(VERSION 3.20)
include(GNUInstallDirs)

if(NOT DEFINED MMCFILTERS_BUILD_DIR)
    message(FATAL_ERROR "MMCFILTERS_BUILD_DIR is required")
endif()
if(NOT DEFINED MMCFILTERS_PROJECT_ROOT)
    message(FATAL_ERROR "MMCFILTERS_PROJECT_ROOT is required")
endif()
if(NOT DEFINED MMCFILTERS_PYTHON_EXECUTABLE)
    message(FATAL_ERROR "MMCFILTERS_PYTHON_EXECUTABLE is required")
endif()

set(work_dir "${MMCFILTERS_BUILD_DIR}/installed-consumer-test")
set(prefix "${work_dir}/prefix")
set(consumer_source_dir "${work_dir}/consumer")
set(consumer_build_dir "${work_dir}/consumer-build")
set(consumer_config "Release")

file(REMOVE_RECURSE "${work_dir}")
file(MAKE_DIRECTORY "${consumer_source_dir}")

set(install_command
    "${CMAKE_COMMAND}" --install "${MMCFILTERS_BUILD_DIR}" --prefix "${prefix}")

if(DEFINED MMCFILTERS_CTEST_CONFIG AND NOT "${MMCFILTERS_CTEST_CONFIG}" STREQUAL "")
    list(APPEND install_command --config "${MMCFILTERS_CTEST_CONFIG}")
endif()

execute_process(
    COMMAND ${install_command}
    RESULT_VARIABLE install_result
    OUTPUT_VARIABLE install_output
    ERROR_VARIABLE install_error)

if(NOT install_result EQUAL 0)
    message(FATAL_ERROR
        "Failed to install mmcfilters for consumer test.\n"
        "stdout:\n${install_output}\n"
        "stderr:\n${install_error}")
endif()

if(EXISTS "${prefix}/${CMAKE_INSTALL_LIBDIR}/mmcfilters")
    message(FATAL_ERROR
        "Plain CMake install must not publish a partial Python package. "
        "Install Python with pip/scikit-build-core instead.")
endif()

set(vocabulary_audit_command
    "${MMCFILTERS_PYTHON_EXECUTABLE}"
    "${MMCFILTERS_PROJECT_ROOT}/unit-tests/cmake/audit_api_refactoring.py"
    --root "${MMCFILTERS_PROJECT_ROOT}"
)
set(refactoring_status
    "${MMCFILTERS_PROJECT_ROOT}/docs/refactoring/api-migration-status.json")
set(refactoring_vocabulary
    "${MMCFILTERS_PROJECT_ROOT}/docs/refactoring/api-vocabulary-audit.json")
if(EXISTS "${refactoring_status}" AND EXISTS "${refactoring_vocabulary}")
    list(APPEND vocabulary_audit_command
        --status "${refactoring_status}"
        --vocabulary "${refactoring_vocabulary}"
        --artifact-root "${prefix}/${CMAKE_INSTALL_INCLUDEDIR}"
    )
endif()

execute_process(
    COMMAND ${vocabulary_audit_command}
    RESULT_VARIABLE vocabulary_audit_result
    OUTPUT_VARIABLE vocabulary_audit_output
    ERROR_VARIABLE vocabulary_audit_error)

if(NOT vocabulary_audit_result EQUAL 0)
    message(FATAL_ERROR
        "Installed-header vocabulary audit failed.\n"
        "stdout:\n${vocabulary_audit_output}\n"
        "stderr:\n${vocabulary_audit_error}")
endif()

if(EXISTS "${prefix}/${CMAKE_INSTALL_INCLUDEDIR}/mmcfilters/filters/AttributeOpeningPrimitivesFamily.hpp")
    message(FATAL_ERROR "Removed AttributeOpeningPrimitivesFamily header must not be installed.")
endif()

if(EXISTS "${prefix}/${CMAKE_INSTALL_INCLUDEDIR}/mmcfilters/attributes/computers/AttributeComputerTraits.hpp")
    message(FATAL_ERROR "Removed AttributeComputerTraits header must not be installed.")
endif()

if(EXISTS "${prefix}/${CMAKE_INSTALL_INCLUDEDIR}/mmcfilters/attributes/computers/AttributeComputerProtocol.hpp")
    message(FATAL_ERROR "Removed AttributeComputerProtocol header must not be installed.")
endif()

foreach(internal_header IN ITEMS
    "attributes/computers/detail/maxdist/EdtDIFT.hpp"
    "attributes/computers/detail/maxdist/Geometry.hpp"
    "attributes/computers/detail/maxdist/PQueue.hpp"
    "dataStructure/FastStack.hpp"
    "trees/BuilderMorphologicalTreeByUnionFind.hpp"
    "trees/ComponentTreeKind.hpp"
    "trees/detail/TreeKindValidation.hpp"
    "utils/AdjacencyRelation.hpp")
    if(EXISTS "${prefix}/${CMAKE_INSTALL_INCLUDEDIR}/mmcfilters/${internal_header}")
        message(FATAL_ERROR
            "Internal implementation header must not be installed: ${internal_header}")
    endif()
endforeach()

file(WRITE "${consumer_source_dir}/CMakeLists.txt"
"cmake_minimum_required(VERSION 3.20)\n"
"project(mmcfilters_installed_consumer LANGUAGES CXX)\n"
"find_package(mmcfilters CONFIG REQUIRED)\n"
"if(NOT TARGET mmcfilters::core)\n"
"    message(FATAL_ERROR \"mmcfilters::core target was not exported\")\n"
"endif()\n"
"add_executable(consumer main.cpp)\n"
"set_property(TARGET consumer PROPERTY NO_SYSTEM_FROM_IMPORTED TRUE)\n"
"target_link_libraries(consumer PRIVATE mmcfilters::core)\n")

file(WRITE "${consumer_source_dir}/main.cpp"
"#include <array>\n"
"#include <cassert>\n"
"#include <cstdint>\n"
"#include <optional>\n"
"#include <span>\n"
"#include <type_traits>\n"
"#include <vector>\n"
"#include <mmcfilters/attributes/AttributeRegistry.hpp>\n"
"#include <mmcfilters/attributes/Attributes.hpp>\n"
"#include <mmcfilters/attributes/computers/AttributeComputerRegistry.hpp>\n"
"#include <mmcfilters/attributes/computers/MaxDistExactComputer.hpp>\n"
"#include <mmcfilters/attributes/computers/MaxDistComputer.hpp>\n"
"#include <mmcfilters/attributes/computers/BitquadAttributeComputer.hpp>\n"
"#include <mmcfilters/attributes/computers/ContourSideAttributeComputer.hpp>\n"
"#include <mmcfilters/filters/MSERComputer.hpp>\n"
"#include <mmcfilters/filters/AttributeFilters.hpp>\n"
"#include <mmcfilters/filters/NodePreservationStability.hpp>\n"
"#include <mmcfilters/filters/ExtinctionValues.hpp>\n"
"#include <mmcfilters/filters/UltimateAttributeOpening.hpp>\n"
"#include <mmcfilters/localAttributes/FiniteWindowLocalAttributeComputer.hpp>\n"
"#include <mmcfilters/trees/MorphologicalTreeFactory.hpp>\n"
"#include <mmcfilters/trees/NativeHierarchy.hpp>\n"
"#include <mmcfilters/trees/TreeOfShapesProducer.hpp>\n"
"#include <mmcfilters/trees/saliency/HierarchySaliencyMap.hpp>\n"
"#include <mmcfilters/trees/saliency/HierarchySaliencyMapProjection.hpp>\n"
"#include <mmcfilters/trees/saliency/HierarchySaliencyMapValidation.hpp>\n"
"#include <mmcfilters/trees/saliency/ShapeSpaceSaliency.hpp>\n"
"#include <mmcfilters/trees/sdrt/SaturatedResidualTreeBuilder.hpp>\n"
"#include <mmcfilters/trees/sdrt/ResidualTreePolicies.hpp>\n"
"#include <mmcfilters/trees/sdrt/UnrestrictedResidualTreeBuilder.hpp>\n"
"#include <mmcfilters/trees/TreeAltitudeAlgorithms.hpp>\n"
"#include <mmcfilters/trees/ValuedMorphologicalTreeView.hpp>\n"
"#include <mmcfilters/trees/adjust/CasfComponentTrees.hpp>\n"
"#include <mmcfilters/utils/Image.hpp>\n"
"#include <mmcfilters/utils/RegularGridAdjacency2D.hpp>\n"
"\n"
"int main()\n"
"{\n"
"    std::array<std::uint8_t, 4> pixels{1, 2, 3, 4};\n"
"    auto image = mmcfilters::ImageUInt8::fromExternal(pixels.data(), 2, 2);\n"
"    auto tree = mmcfilters::MorphologicalTreeFactory::createMaxTree(image);\n"
"    const mmcfilters::RegularGridAdjacency2D residualAdjacency(2, 2, 1.0);\n"
"    mmcfilters::sdrt::UnrestrictedResidualTreeBuilder<std::uint8_t> unrestrictedBuilder(residualAdjacency);\n"
"    mmcfilters::sdrt::SaturatedResidualTreeBuilder<std::uint8_t> saturatedBuilder(residualAdjacency, mmcfilters::PixelId{0});\n"
"    mmcfilters::sdrt::UnrestrictedResidualTreeOptions unrestrictedOptions;\n"
"    mmcfilters::sdrt::SaturatedResidualTreeOptions saturatedOptions;\n"
"    assert(unrestrictedBuilder.spatialOrder().isRowMajor());\n"
"    assert(saturatedBuilder.infinityPixel() == 0);\n"
"    auto residualTree = mmcfilters::MorphologicalTreeFactory::createUnrestrictedResidualTree(image, 1.0, unrestrictedOptions);\n"
"    auto saturatedResidualTree = mmcfilters::MorphologicalTreeFactory::createSaturatedResidualTree(image, mmcfilters::PixelId{0}, 1.0, saturatedOptions);\n"
"    assert(tree.topology().numNodes() > 0);\n"
"    assert(residualTree.topology().kind() == mmcfilters::MorphologicalTreeKind::UnrestrictedResidualTree);\n"
"    assert(saturatedResidualTree.topology().kind() == mmcfilters::MorphologicalTreeKind::SaturatedResidualTree);\n"
"    assert(residualTree.reconstructFromNodeAltitudes()->getSize() == image->getSize());\n"
"    assert(saturatedResidualTree.reconstructFromNodeAltitudes()->getSize() == image->getSize());\n"
"    assert(tree.topology().numRows() == 2);\n"
"    assert(tree.topology().numColumns() == 2);\n"
"    assert(tree.topology().sharedAdjacencyContext() != nullptr);\n"
"    const auto neighbors = tree.topology().sharedAdjacencyContext()->adjacency.getNeighborIndices(0);\n"
"    assert(neighbors.begin() != neighbors.end());\n"
"    const auto maxDistExactRequirements =\n"
"        mmcfilters::attributes::registry::capabilityRequirements(\n"
"            mmcfilters::Attribute::MaxDistExact);\n"
"    assert(!maxDistExactRequirements.altitude);\n"
"    assert(maxDistExactRequirements.gridDomain2D);\n"
"    assert(!maxDistExactRequirements.monotoneAltitudeOrder);\n"
"    assert(maxDistExactRequirements.adjacency ==\n"
"        mmcfilters::attributes::registry::AttributeAdjacencyRequirement::None);\n"
"    const auto maxDistRequirements =\n"
"        mmcfilters::attributes::registry::capabilityRequirements(\n"
"            mmcfilters::Attribute::MaxDist);\n"
"    assert(!maxDistRequirements.altitude);\n"
"    assert(maxDistRequirements.gridDomain2D);\n"
"    assert(maxDistRequirements.adjacency ==\n"
"        mmcfilters::attributes::registry::AttributeAdjacencyRequirement::None);\n"
"    const std::array<mmcfilters::NodeId, 2> nativeParent{0, 0};\n"
"    const std::array<mmcfilters::NodeId, 1> nativeSmallestNodeMap{1};\n"
"    const std::array<std::uint8_t, 2> nativeAltitude{10, 3};\n"
"    auto abstractTree =\n"
"        mmcfilters::MorphologicalTreeFactory::createFromNativeHierarchy<std::uint8_t>(\n"
"            mmcfilters::NativeHierarchyView<std::uint8_t>{\n"
"                std::span<const mmcfilters::NodeId>(nativeParent),\n"
"                std::span<const mmcfilters::NodeId>(nativeSmallestNodeMap),\n"
"                std::span<const std::uint8_t>(nativeAltitude),\n"
"                0,\n"
"                std::nullopt,\n"
"                mmcfilters::MorphologicalTreeSemantics{}});\n"
"    assert(!abstractTree.topology().hasGridDomain2D());\n"
"    auto exactTreeOfShapes =\n"
"        mmcfilters::MorphologicalTreeFactory::createTreeOfShapes<mmcfilters::ToSGrayLevel>(\n"
"            image,\n"
"            mmcfilters::selfDualSpanConvention(\n"
"                mmcfilters::TopographicDomainExtension::None,\n"
"                0));\n"
"    static_assert(std::is_same_v<\n"
"        decltype(exactTreeOfShapes),\n"
"        mmcfilters::ValuedMorphologicalTree<mmcfilters::ToSGrayLevel>>);\n"
"    assert(exactTreeOfShapes.topology().numRows() == 2);\n"
"    assert(exactTreeOfShapes.topology().numColumns() == 2);\n"
"    assert(exactTreeOfShapes.topology().topographicConvention() != nullptr);\n"
"    auto defaultTreeOfShapes =\n"
"        mmcfilters::MorphologicalTreeFactory::createTreeOfShapes(image);\n"
"    static_assert(std::is_same_v<\n"
"        decltype(defaultTreeOfShapes),\n"
"        mmcfilters::ValuedMorphologicalTree<std::uint8_t>>);\n"
"    assert(defaultTreeOfShapes.topology().topographicConvention()->altitudeEncoding ==\n"
"           mmcfilters::TopographicAltitudeEncoding::UInt8);\n"
"    auto reconstructed = tree.reconstructFromNodeAltitudes();\n"
"    assert(reconstructed->getSize() == 4);\n"
"    auto hierarchyLevels = mmcfilters::HierarchySaliencyMap::computeTopologicalLevels(tree.topology());\n"
"    mmcfilters::HierarchySaliencyMapValidation::validateHierarchyValuation(\n"
"        tree.topology(),\n"
"        std::span<const int>(hierarchyLevels),\n"
"        mmcfilters::HierarchyValuationPolicy::RequireStrictHierarchy,\n"
"        mmcfilters::HierarchyValuationRangePolicy::RequireNonNegative);\n"
"    auto hierarchyEdgeMap = mmcfilters::HierarchySaliencyMap::computeTopologicalLevelEdgeMap(tree.topology());\n"
"    assert(hierarchyEdgeMap.sources.size() == hierarchyEdgeMap.targets.size());\n"
"    assert(hierarchyEdgeMap.sources.size() == hierarchyEdgeMap.values.size());\n"
"    auto hierarchyCut = mmcfilters::HierarchySaliencyMapProjection::thresholdCut(hierarchyEdgeMap, 1);\n"
"    assert(hierarchyCut.sources.size() == hierarchyCut.targets.size());\n"
"    std::vector<double> shapeAttribute(\n"
"        static_cast<std::size_t>(tree.topology().numInternalNodeSlots()), 0.0);\n"
"    auto shapeSpaceMap = mmcfilters::ShapeSpaceSaliency::compute(\n"
"        tree.topology(),\n"
"        std::span<const double>(shapeAttribute),\n"
"        mmcfilters::ShapeSpaceExtremaPolarity::Maxima);\n"
"    assert(shapeSpaceMap.nodeScores.size() == shapeAttribute.size());\n"
"    assert(shapeSpaceMap.edgeMap.sources.size() == shapeSpaceMap.edgeMap.values.size());\n"
"    auto publicAttributes = mmcfilters::AttributeComputation::computeAttributes(\n"
"        tree,\n"
"        std::vector<mmcfilters::AttributeOrGroup>{mmcfilters::Attribute::Area, mmcfilters::Attribute::MeanGrayLevel});\n"
"    assert(publicAttributes.first.NUM_ATTRIBUTES == 2);\n"
"    assert(publicAttributes.second.size() == static_cast<std::size_t>(\n"
"        tree.topology().numInternalNodeSlots() * publicAttributes.first.NUM_ATTRIBUTES));\n"
"    auto publicTopologyAttributes = mmcfilters::AttributeComputation::computeTopologyAttributes(\n"
"        tree.topology(),\n"
"        std::vector<mmcfilters::AttributeOrGroup>{mmcfilters::Attribute::Area, mmcfilters::Attribute::BoxWidth, mmcfilters::Attribute::MaxDistExact});\n"
"    assert(publicTopologyAttributes.first.NUM_ATTRIBUTES == 3);\n"
"    assert(publicTopologyAttributes.second.size() == static_cast<std::size_t>(\n"
"        tree.topology().numInternalNodeSlots() * publicTopologyAttributes.first.NUM_ATTRIBUTES));\n"
"    assert(!mmcfilters::attributes::computers::runtimeProducedAttributes<mmcfilters::attributes::computers::BitquadAttributeComputer>().empty());\n"
"    assert(!mmcfilters::attributes::computers::runtimeProducedAttributes<mmcfilters::attributes::computers::ContourSideAttributeComputer>().empty());\n"
"    assert(!mmcfilters::attributes::computers::runtimeProducedAttributes<mmcfilters::attributes::computers::MaxDistExactComputer>().empty());\n"
"    assert(!mmcfilters::attributes::computers::runtimeProducedAttributes<mmcfilters::attributes::computers::MaxDistComputer>().empty());\n"
"    auto publicMaxDistExact = mmcfilters::AttributeComputation::computeSingleAttribute(tree, mmcfilters::Attribute::MaxDistExact);\n"
"    assert(publicMaxDistExact.first.NUM_ATTRIBUTES == 1);\n"
"    assert(publicMaxDistExact.second.size() == static_cast<std::size_t>(tree.topology().numInternalNodeSlots()));\n"
"    auto publicMaxDist =\n"
"        mmcfilters::AttributeComputation::computeSingleTopologyAttribute(\n"
"            tree.topology(), mmcfilters::Attribute::MaxDist);\n"
"    assert(publicMaxDist.first.NUM_ATTRIBUTES == 1);\n"
"    assert(publicMaxDist.second.size() == static_cast<std::size_t>(tree.topology().numInternalNodeSlots()));\n"
"    mmcfilters::local_attributes::WindowOffset offset{0, 0};\n"
"    assert(offset.rowOffset == 0 && offset.columnOffset == 0);\n"
"    mmcfilters::local_attributes::ObservationWindow observationWindow{{0, 0}, {0, 1}};\n"
"    assert(observationWindow.size() == 2);\n"
"    struct InstalledLocalDecision {\n"
"        using Value = int;\n"
"        Value evaluateLocalDecision(mmcfilters::local_attributes::BinaryVisibilityState state) const { return static_cast<int>(state.bits()); }\n"
"    };\n"
"    struct InstalledEventAlgebra {\n"
"        int additiveIdentity() const { return 0; }\n"
"        void addAssign(int& target, const int& source) const { target += source; }\n"
"        void subtractAssign(int& target, const int& source) const { target -= source; }\n"
"    };\n"
"    static_assert(mmcfilters::local_attributes::LocalDecision<InstalledLocalDecision>);\n"
"    static_assert(mmcfilters::local_attributes::EventAlgebra<InstalledEventAlgebra, int>);\n"
"    const InstalledLocalDecision installedLocalDecision;\n"
"    const InstalledEventAlgebra installedEventAlgebra;\n"
"    auto eventDeltas = mmcfilters::local_attributes::FiniteWindowLocalAttributeComputer::computeEventDeltas(\n"
"        tree.topology(), 0, observationWindow, installedLocalDecision, installedEventAlgebra);\n"
"    auto localAttributeIncrements = mmcfilters::local_attributes::FiniteWindowLocalAttributeComputer::computeLocalAttributeIncrements(\n"
"        tree.topology(), observationWindow, installedLocalDecision, installedEventAlgebra);\n"
"    auto nodeAttributes = mmcfilters::local_attributes::FiniteWindowLocalAttributeComputer::aggregateEventIncrements(\n"
"        tree.topology(), std::span<const mmcfilters::local_attributes::LocalAttributeIncrement<int>>(localAttributeIncrements), installedEventAlgebra);\n"
"    auto computedNodeAttributes = mmcfilters::local_attributes::FiniteWindowLocalAttributeComputer::compute(\n"
"        tree.topology(), observationWindow, installedLocalDecision, installedEventAlgebra);\n"
"    static_assert(std::is_same_v<decltype(eventDeltas), std::vector<mmcfilters::local_attributes::EventDelta<int>>>);\n"
"    static_assert(std::is_same_v<decltype(nodeAttributes), std::vector<mmcfilters::local_attributes::NodeAttribute<int>>>);\n"
"    assert(!eventDeltas.empty());\n"
"    assert(nodeAttributes.size() == static_cast<std::size_t>(tree.topology().numInternalNodeSlots()));\n"
"    assert(computedNodeAttributes == nodeAttributes);\n"
"    auto publicSamples = mmcfilters::AttributeComputation::computeSampledNodeAttribute(\n"
"        tree,\n"
"        mmcfilters::Attribute::MeanGrayLevel,\n"
"        mmcfilters::AltitudeDifference<std::uint8_t>{1},\n"
"        1);\n"
"    assert(publicSamples.first.NUM_ATTRIBUTES == 3);\n"
"    assert(publicSamples.second.size() == static_cast<std::size_t>(\n"
"        tree.topology().numInternalNodeSlots() * publicSamples.first.NUM_ATTRIBUTES));\n"
"\n"
"    std::array<std::int32_t, 4> intPixels{1, 2, 3, 4};\n"
"    auto intImage = mmcfilters::ImageInt32::fromExternal(intPixels.data(), 2, 2);\n"
"    auto typedTree = mmcfilters::MorphologicalTreeFactory::createMaxTree(intImage);\n"
"    static_assert(std::is_same_v<decltype(typedTree), mmcfilters::ValuedMorphologicalTree<std::int32_t>>);\n"
"    auto typedReconstruction = typedTree.reconstructFromNodeAltitudes();\n"
"    static_assert(std::is_same_v<decltype(typedReconstruction), mmcfilters::ImageInt32Ptr>);\n"
"    assert(typedReconstruction->getSize() == intImage->getSize());\n"
"    for (int i = 0; i < intImage->getSize(); ++i) {\n"
"        assert((*typedReconstruction)[i] == (*intImage)[i]);\n"
"    }\n"
"    auto typedView = typedTree.asView();\n"
"    static_assert(std::is_same_v<decltype(typedView), mmcfilters::ValuedMorphologicalTreeView<std::int32_t>>);\n"
"    mmcfilters::TreeAltitudeAlgorithms::validateNodeAltitudeBufferShape(typedView.topology(), typedView.nodeAltitudes());\n"
"    auto [typedParent, typedAltitude] = typedTree.exportHigraHierarchy();\n"
"    static_assert(std::is_same_v<decltype(typedAltitude), std::vector<std::int32_t>>);\n"
"    assert(typedParent.size() == typedAltitude.size());\n"
"    auto typedAttributes = mmcfilters::AttributeComputation::computeAttributes(\n"
"        typedTree,\n"
"        std::vector<mmcfilters::AttributeOrGroup>{mmcfilters::Attribute::MeanGrayLevel, mmcfilters::Attribute::GrayLevelHeight});\n"
"    assert(typedAttributes.first.NUM_ATTRIBUTES == 2);\n"
"    assert(typedAttributes.second.size() == static_cast<std::size_t>(\n"
"        typedView.topology().numInternalNodeSlots() * typedAttributes.first.NUM_ATTRIBUTES));\n"
"    auto mappedMeanGrayLevel = mmcfilters::AttributeComputation::computeAttributeMapping(\n"
"        typedTree,\n"
"        mmcfilters::Attribute::MeanGrayLevel);\n"
"    assert(mappedMeanGrayLevel->getSize() == intImage->getSize());\n"
"    mmcfilters::NodePreservationMask keepAll(typedTree.topology().numInternalNodeSlots(), true);\n"
"    std::vector<float> thresholdAttributes(static_cast<std::size_t>(typedTree.topology().numInternalNodeSlots()), 1.0f);\n"
"    auto thresholdMask = mmcfilters::computeNodePreservationMask<float>(thresholdAttributes, 0.5f);\n"
"    auto altitudeAdjusted = mmcfilters::adjustNodePreservationMaskByAltitudeStability(\n"
"        typedTree, thresholdMask, mmcfilters::AltitudeDifference<std::int32_t>{1});\n"
"    auto depthAdjusted = mmcfilters::adjustNodePreservationMaskByDepthStability(typedTree, thresholdMask, 1);\n"
"    static_assert(std::is_same_v<decltype(altitudeAdjusted), mmcfilters::NodePreservationMask>);\n"
"    static_assert(std::is_same_v<decltype(depthAdjusted), mmcfilters::NodePreservationMask>);\n"
"    assert(altitudeAdjusted.decisions() == thresholdMask.decisions());\n"
"    assert(depthAdjusted.decisions() == thresholdMask.decisions());\n"
"    mmcfilters::DirectAttributeFilter<std::int32_t> typedFilter(typedTree);\n"
"    auto typedFiltered = typedFilter.applyDirectAttributeFilter(keepAll);\n"
"    static_assert(std::is_same_v<decltype(typedFiltered), mmcfilters::ImageInt32Ptr>);\n"
"    assert(typedFiltered->getSize() == intImage->getSize());\n"
"    for (int i = 0; i < intImage->getSize(); ++i) {\n"
"        assert((*typedFiltered)[i] == (*intImage)[i]);\n"
"    }\n"
"    auto typedMeanGrayLevel = mmcfilters::AttributeComputation::computeSingleAttribute(\n"
"        typedTree,\n"
"        mmcfilters::Attribute::MeanGrayLevel);\n"
"    auto typedMeanGrayLevelSamples = mmcfilters::AttributeComputation::computeSampledNodeAttribute(\n"
"        typedTree,\n"
"        mmcfilters::Attribute::MeanGrayLevel,\n"
"        mmcfilters::AltitudeDifference<std::int32_t>{1},\n"
"        1);\n"
"    assert(typedMeanGrayLevelSamples.first.NUM_ATTRIBUTES == 3);\n"
"    assert(typedMeanGrayLevelSamples.second.size() == static_cast<std::size_t>(\n"
"        typedTree.topology().numInternalNodeSlots() * typedMeanGrayLevelSamples.first.NUM_ATTRIBUTES));\n"
"    mmcfilters::ExtinctionValues<std::int32_t> typedExtinction(typedTree, typedMeanGrayLevel.second);\n"
"    auto typedExtinctionFiltered = typedExtinction.filtering(mmcfilters::ExtinctionSelectionPolicy<float>::byTopK(1024));\n"
"    static_assert(std::is_same_v<decltype(typedExtinctionFiltered), mmcfilters::ImageInt32Ptr>);\n"
"    assert(typedExtinctionFiltered->getSize() == intImage->getSize());\n"
"    mmcfilters::UltimateAttributeOpening<std::int32_t> typedUao(typedTree, typedMeanGrayLevel.second);\n"
"    typedUao.execute(intImage->getSize());\n"
"    auto typedMaxContrast = typedUao.getMaxContrastImage();\n"
"    static_assert(std::is_same_v<decltype(typedMaxContrast), mmcfilters::ImageInt32Ptr>);\n"
"    assert(typedMaxContrast->getSize() == intImage->getSize());\n"
"    mmcfilters::MSERComputer<std::int32_t> typedMser(typedTree);\n"
"    auto typedMserFlags = typedMser.computeMSER(mmcfilters::AltitudeDifference<std::int32_t>{1});\n"
"    assert(typedMserFlags.size() == static_cast<std::size_t>(typedTree.topology().numInternalNodeSlots()));\n"
"    typedUao.executeWithMSER(intImage->getSize(), mmcfilters::AltitudeDifference<std::int32_t>{1});\n"
"    auto typedAssociated = typedUao.getAssociatedImage();\n"
"    static_assert(std::is_same_v<decltype(typedAssociated), mmcfilters::ImageInt32Ptr>);\n"
"    assert(typedAssociated->getSize() == intImage->getSize());\n"
"    mmcfilters::adjust::CasfComponentTrees<std::int32_t> typedCasf(\n"
"        intImage,\n"
"        mmcfilters::adjust::CasfComponentTreesAttribute::Area);\n"
	"    auto typedCasfFiltered = typedCasf.filter({1.0});\n"
	"    static_assert(std::is_same_v<decltype(typedCasfFiltered), mmcfilters::ImageInt32Ptr>);\n"
	"    assert(typedCasfFiltered->getSize() == intImage->getSize());\n"
	"    auto typedProjectionAttrs = mmcfilters::AttributeComputation::computeAttributes(\n"
	"        typedTree,\n"
	"        std::vector<mmcfilters::AttributeOrGroup>{mmcfilters::Attribute::MeanGrayLevel, mmcfilters::Attribute::Volume});\n"
	"    auto typedProjectedHigra = mmcfilters::AttributeComputation::projectNodeValuesToExportedHigra(\n"
	"        typedTree,\n"
	"        typedProjectionAttrs.first,\n"
	"        typedProjectionAttrs.second);\n"
	"    assert(typedProjectedHigra.size() == static_cast<std::size_t>(\n"
	"        (typedTree.topology().numPixels() + typedTree.topology().numNodes()) * typedProjectionAttrs.first.NUM_ATTRIBUTES));\n"
"\n"
"    std::array<float, 4> floatPixels{1.0f, 2.5f, 3.0f, 4.5f};\n"
"    auto floatImage = mmcfilters::ImageFloat::fromExternal(floatPixels.data(), 2, 2);\n"
"    auto floatTree = mmcfilters::MorphologicalTreeFactory::createMaxTree(floatImage);\n"
"    static_assert(std::is_same_v<decltype(floatTree), mmcfilters::ValuedMorphologicalTree<float>>);\n"
"    auto floatReconstruction = floatTree.reconstructFromNodeAltitudes();\n"
"    static_assert(std::is_same_v<decltype(floatReconstruction), mmcfilters::ImageFloatPtr>);\n"
"    assert(floatReconstruction->getSize() == floatImage->getSize());\n"
"    for (int i = 0; i < floatImage->getSize(); ++i) {\n"
"        assert((*floatReconstruction)[i] == (*floatImage)[i]);\n"
"    }\n"
"    auto floatView = floatTree.asView();\n"
"    static_assert(std::is_same_v<decltype(floatView), mmcfilters::ValuedMorphologicalTreeView<float>>);\n"
"    mmcfilters::TreeAltitudeAlgorithms::validateNodeAltitudeBufferShape(floatView.topology(), floatView.nodeAltitudes());\n"
"    auto [floatParent, floatAltitude] = floatTree.exportHigraHierarchy();\n"
"    static_assert(std::is_same_v<decltype(floatAltitude), std::vector<float>>);\n"
"    assert(floatParent.size() == floatAltitude.size());\n"
"    auto floatMeanGrayLevel = mmcfilters::AttributeComputation::computeSingleAttribute(\n"
"        floatTree,\n"
"        mmcfilters::Attribute::MeanGrayLevel);\n"
"    auto floatMeanGrayLevelSamples = mmcfilters::AttributeComputation::computeSampledNodeAttribute(\n"
"        floatTree,\n"
"        mmcfilters::Attribute::MeanGrayLevel,\n"
"        mmcfilters::AltitudeDifference<float>{0.5f},\n"
"        1);\n"
"    assert(floatMeanGrayLevelSamples.first.NUM_ATTRIBUTES == 3);\n"
"    assert(floatMeanGrayLevelSamples.second.size() == static_cast<std::size_t>(\n"
"        floatTree.topology().numInternalNodeSlots() * floatMeanGrayLevelSamples.first.NUM_ATTRIBUTES));\n"
"    mmcfilters::NodePreservationMask keepAllFloat(floatTree.topology().numInternalNodeSlots(), true);\n"
"    mmcfilters::DirectAttributeFilter<float> floatFilter(floatTree);\n"
"    auto floatFiltered = floatFilter.applyDirectAttributeFilter(keepAllFloat);\n"
"    static_assert(std::is_same_v<decltype(floatFiltered), mmcfilters::ImageFloatPtr>);\n"
"    assert(floatFiltered->getSize() == floatImage->getSize());\n"
"    mmcfilters::ExtinctionValues<float> floatExtinction(floatTree, floatMeanGrayLevel.second);\n"
"    auto floatExtinctionFiltered = floatExtinction.filtering(mmcfilters::ExtinctionSelectionPolicy<float>::byTopK(1024));\n"
"    static_assert(std::is_same_v<decltype(floatExtinctionFiltered), mmcfilters::ImageFloatPtr>);\n"
"    assert(floatExtinctionFiltered->getSize() == floatImage->getSize());\n"
"    mmcfilters::UltimateAttributeOpening<float> floatUao(floatTree, floatMeanGrayLevel.second);\n"
"    floatUao.execute(floatImage->getSize());\n"
"    auto floatMaxContrast = floatUao.getMaxContrastImage();\n"
"    static_assert(std::is_same_v<decltype(floatMaxContrast), mmcfilters::ImageFloatPtr>);\n"
"    assert(floatMaxContrast->getSize() == floatImage->getSize());\n"
"    mmcfilters::MSERComputer<float> floatMser(floatTree);\n"
"    auto floatMserFlags = floatMser.computeMSER(0.5f);\n"
"    assert(floatMserFlags.size() == static_cast<std::size_t>(floatTree.topology().numInternalNodeSlots()));\n"
"    floatUao.executeWithMSER(floatImage->getSize(), 0.5f);\n"
"    mmcfilters::adjust::CasfComponentTrees<float> floatCasf(\n"
"        floatImage,\n"
"        mmcfilters::adjust::CasfComponentTreesAttribute::Area);\n"
	"    auto floatCasfFiltered = floatCasf.filter({1.0});\n"
	"    static_assert(std::is_same_v<decltype(floatCasfFiltered), mmcfilters::ImageFloatPtr>);\n"
	"    assert(floatCasfFiltered->getSize() == floatImage->getSize());\n"
	"    auto floatProjectionAttrs = mmcfilters::AttributeComputation::computeAttributes(\n"
	"        floatTree,\n"
	"        std::vector<mmcfilters::AttributeOrGroup>{mmcfilters::Attribute::MeanGrayLevel, mmcfilters::Attribute::Volume});\n"
	"    auto floatProjectedHigra = mmcfilters::AttributeComputation::projectNodeValuesToExportedHigra(\n"
	"        floatTree,\n"
	"        floatProjectionAttrs.first,\n"
	"        floatProjectionAttrs.second);\n"
	"    assert(floatProjectedHigra.size() == static_cast<std::size_t>(\n"
	"        (floatTree.topology().numPixels() + floatTree.topology().numNodes()) * floatProjectionAttrs.first.NUM_ATTRIBUTES));\n"
	"    return 0;\n"
"}\n")

execute_process(
    COMMAND
        "${CMAKE_COMMAND}"
        -S "${consumer_source_dir}"
        -B "${consumer_build_dir}"
        "-DCMAKE_PREFIX_PATH=${prefix}"
        "-DCMAKE_BUILD_TYPE=${consumer_config}"
    RESULT_VARIABLE configure_result
    OUTPUT_VARIABLE configure_output
    ERROR_VARIABLE configure_error)

if(NOT configure_result EQUAL 0)
    message(FATAL_ERROR
        "Failed to configure installed mmcfilters consumer.\n"
        "stdout:\n${configure_output}\n"
        "stderr:\n${configure_error}")
endif()

set(consumer_build_command
    "${CMAKE_COMMAND}" --build "${consumer_build_dir}" --parallel
    --config "${consumer_config}")

execute_process(
    COMMAND ${consumer_build_command}
    RESULT_VARIABLE build_result
    OUTPUT_VARIABLE build_output
    ERROR_VARIABLE build_error)

if(NOT build_result EQUAL 0)
    message(FATAL_ERROR
        "Failed to build installed mmcfilters consumer.\n"
        "stdout:\n${build_output}\n"
        "stderr:\n${build_error}")
endif()

if(WIN32)
    set(consumer_executable_candidates
        "${consumer_build_dir}/${consumer_config}/consumer.exe"
        "${consumer_build_dir}/consumer.exe")
else()
    set(consumer_executable_candidates
        "${consumer_build_dir}/consumer"
        "${consumer_build_dir}/${consumer_config}/consumer")
endif()

set(consumer_executable "")
foreach(candidate IN LISTS consumer_executable_candidates)
    if(EXISTS "${candidate}")
        set(consumer_executable "${candidate}")
        break()
    endif()
endforeach()

if("${consumer_executable}" STREQUAL "")
    message(FATAL_ERROR
        "Installed mmcfilters consumer executable was not found.\n"
        "Checked:\n${consumer_executable_candidates}")
endif()

execute_process(
    COMMAND "${consumer_executable}"
    RESULT_VARIABLE run_result
    OUTPUT_VARIABLE run_output
    ERROR_VARIABLE run_error)

if(NOT run_result EQUAL 0)
    message(FATAL_ERROR
        "Installed mmcfilters consumer executable failed.\n"
        "stdout:\n${run_output}\n"
        "stderr:\n${run_error}")
endif()
