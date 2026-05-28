cmake_minimum_required(VERSION 3.20)
include(GNUInstallDirs)

if(NOT DEFINED MMCFILTERS_BUILD_DIR)
    message(FATAL_ERROR "MMCFILTERS_BUILD_DIR is required")
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

if(EXISTS "${prefix}/${CMAKE_INSTALL_INCLUDEDIR}/mmcfilters/filters/AttributeOpeningPrimitivesFamily.hpp")
    message(FATAL_ERROR "Removed AttributeOpeningPrimitivesFamily header must not be installed.")
endif()

file(WRITE "${consumer_source_dir}/CMakeLists.txt"
"cmake_minimum_required(VERSION 3.20)\n"
"project(mmcfilters_installed_consumer LANGUAGES CXX)\n"
"find_package(mmcfilters CONFIG REQUIRED)\n"
"if(NOT TARGET mmcfilters::core)\n"
"    message(FATAL_ERROR \"mmcfilters::core target was not exported\")\n"
"endif()\n"
"add_executable(consumer main.cpp)\n"
"target_link_libraries(consumer PRIVATE mmcfilters::core)\n")

file(WRITE "${consumer_source_dir}/main.cpp"
"#include <array>\n"
"#include <cassert>\n"
"#include <cstdint>\n"
"#include <span>\n"
"#include <type_traits>\n"
"#include <vector>\n"
"#include <mmcfilters/attributes/Attributes.hpp>\n"
"#include <mmcfilters/attributes/computers/AttributeComputerTraits.hpp>\n"
"#include <mmcfilters/attributes/computers/MaxDistComputer.hpp>\n"
"#include <mmcfilters/attributes/computers/BitquadAttributeComputer.hpp>\n"
"#include <mmcfilters/attributes/computers/ContourSideAttributeComputer.hpp>\n"
"#include <mmcfilters/filters/MSERComputer.hpp>\n"
"#include <mmcfilters/filters/AttributeFilters.hpp>\n"
"#include <mmcfilters/filters/ExtinctionValues.hpp>\n"
"#include <mmcfilters/filters/UltimateAttributeOpening.hpp>\n"
"#include <mmcfilters/localEvents/LocalEvents.hpp>\n"
"#include <mmcfilters/trees/MorphologicalTreeFactory.hpp>\n"
"#include <mmcfilters/trees/TreeAltitudeAlgorithms.hpp>\n"
"#include <mmcfilters/trees/WeightedTreeView.hpp>\n"
"#include <mmcfilters/trees/adjust/CasfComponentTrees.hpp>\n"
"#include <mmcfilters/utils/Image.hpp>\n"
"\n"
"int main()\n"
"{\n"
"    std::array<std::uint8_t, 4> pixels{1, 2, 3, 4};\n"
"    auto image = mmcfilters::ImageUInt8::fromExternal(pixels.data(), 2, 2);\n"
"    auto tree = mmcfilters::MorphologicalTreeFactory::createMaxTree(image);\n"
"    assert(tree.topology().getNumNodes() > 0);\n"
"    auto reconstructed = tree.reconstructionImage();\n"
"    assert(reconstructed->getSize() == 4);\n"
"    auto publicAttributes = mmcfilters::AttributeComputation::computeAttributes(\n"
"        tree,\n"
"        std::vector<mmcfilters::AttributeOrGroup>{mmcfilters::Attribute::AREA, mmcfilters::Attribute::LEVEL});\n"
"    assert(publicAttributes.first.NUM_ATTRIBUTES == 2);\n"
"    assert(publicAttributes.second.size() == static_cast<std::size_t>(\n"
"        tree.topology().getNumInternalNodeSlots() * publicAttributes.first.NUM_ATTRIBUTES));\n"
"    auto publicTopologyAttributes = mmcfilters::AttributeComputation::computeTopologyAttributes(\n"
"        tree.topology(),\n"
"        std::vector<mmcfilters::AttributeOrGroup>{mmcfilters::Attribute::AREA, mmcfilters::Attribute::BOX_WIDTH});\n"
"    assert(publicTopologyAttributes.first.NUM_ATTRIBUTES == 2);\n"
"    assert(publicTopologyAttributes.second.size() == static_cast<std::size_t>(\n"
"        tree.topology().getNumInternalNodeSlots() * publicTopologyAttributes.first.NUM_ATTRIBUTES));\n"
"    assert(!mmcfilters::attributes::computers::runtimeProducedAttributes<mmcfilters::attributes::computers::BitquadAttributeComputer>().empty());\n"
"    assert(!mmcfilters::attributes::computers::runtimeProducedAttributes<mmcfilters::attributes::computers::ContourSideAttributeComputer>().empty());\n"
"    assert(!mmcfilters::attributes::computers::runtimeProducedAttributes<mmcfilters::attributes::computers::MaxDistComputer>().empty());\n"
"    auto publicMaxDist = mmcfilters::AttributeComputation::computeSingleAttribute(tree, mmcfilters::Attribute::MAX_DIST);\n"
"    assert(publicMaxDist.first.NUM_ATTRIBUTES == 1);\n"
"    assert(publicMaxDist.second.size() == static_cast<std::size_t>(tree.topology().getNumInternalNodeSlots()));\n"
"    mmcfilters::local_events::WindowOffset offset{0, 0};\n"
"    assert(offset.rowOffset == 0 && offset.colOffset == 0);\n"
"    auto publicDelta = mmcfilters::AttributeComputation::computeSingleAttributeWithDelta(\n"
"        tree,\n"
"        mmcfilters::Attribute::LEVEL,\n"
"        mmcfilters::AltitudeDiff<std::uint8_t>{1},\n"
"        1);\n"
"    assert(publicDelta.first.NUM_ATTRIBUTES == 3);\n"
"    assert(publicDelta.second.size() == static_cast<std::size_t>(\n"
"        tree.topology().getNumInternalNodeSlots() * publicDelta.first.NUM_ATTRIBUTES));\n"
"\n"
"    std::array<std::int32_t, 4> intPixels{1, 2, 3, 4};\n"
"    auto intImage = mmcfilters::ImageInt32::fromExternal(intPixels.data(), 2, 2);\n"
"    auto typedTree = mmcfilters::MorphologicalTreeFactory::createMaxTree(intImage);\n"
"    static_assert(std::is_same_v<decltype(typedTree), mmcfilters::WeightedMorphologicalTree<std::int32_t>>);\n"
"    auto typedReconstruction = typedTree.reconstructionImage();\n"
"    static_assert(std::is_same_v<decltype(typedReconstruction), mmcfilters::ImageInt32Ptr>);\n"
"    assert(typedReconstruction->getSize() == intImage->getSize());\n"
"    for (int i = 0; i < intImage->getSize(); ++i) {\n"
"        assert((*typedReconstruction)[i] == (*intImage)[i]);\n"
"    }\n"
"    auto typedView = typedTree.asView();\n"
"    static_assert(std::is_same_v<decltype(typedView), mmcfilters::WeightedTreeView<std::int32_t>>);\n"
"    mmcfilters::TreeAltitudeAlgorithms::validateAltitudeBufferShape(typedView.topology(), typedView.altitude());\n"
"    auto [typedParent, typedAltitude] = typedTree.exportHigraHierarchy();\n"
"    static_assert(std::is_same_v<decltype(typedAltitude), std::vector<std::int32_t>>);\n"
"    assert(typedParent.size() == typedAltitude.size());\n"
"    auto typedAttributes = mmcfilters::AttributeComputation::computeAttributes(\n"
"        typedTree,\n"
"        std::vector<mmcfilters::AttributeOrGroup>{mmcfilters::Attribute::LEVEL, mmcfilters::Attribute::GRAY_HEIGHT});\n"
"    assert(typedAttributes.first.NUM_ATTRIBUTES == 2);\n"
"    assert(typedAttributes.second.size() == static_cast<std::size_t>(\n"
"        typedView.topology().getNumInternalNodeSlots() * typedAttributes.first.NUM_ATTRIBUTES));\n"
"    auto mappedLevel = mmcfilters::AttributeComputation::computeAttributeMapping(\n"
"        typedTree,\n"
"        mmcfilters::Attribute::LEVEL);\n"
"    assert(mappedLevel->getSize() == intImage->getSize());\n"
"    std::vector<bool> keepAll(typedTree.topology().getNumInternalNodeSlots(), true);\n"
"    mmcfilters::AttributeFilters<std::int32_t> typedFilters(typedTree);\n"
"    auto typedFiltered = typedFilters.filteringByDirectRule(keepAll);\n"
"    static_assert(std::is_same_v<decltype(typedFiltered), mmcfilters::ImageInt32Ptr>);\n"
"    assert(typedFiltered->getSize() == intImage->getSize());\n"
"    for (int i = 0; i < intImage->getSize(); ++i) {\n"
"        assert((*typedFiltered)[i] == (*intImage)[i]);\n"
"    }\n"
"    auto typedLevel = mmcfilters::AttributeComputation::computeSingleAttribute(\n"
"        typedTree,\n"
"        mmcfilters::Attribute::LEVEL);\n"
"    auto typedLevelDelta = mmcfilters::AttributeComputation::computeSingleAttributeWithDelta(\n"
"        typedTree,\n"
"        mmcfilters::Attribute::LEVEL,\n"
"        mmcfilters::AltitudeDiff<std::int32_t>{1},\n"
"        1);\n"
"    assert(typedLevelDelta.first.NUM_ATTRIBUTES == 3);\n"
"    assert(typedLevelDelta.second.size() == static_cast<std::size_t>(\n"
"        typedTree.topology().getNumInternalNodeSlots() * typedLevelDelta.first.NUM_ATTRIBUTES));\n"
"    mmcfilters::ExtinctionValues<std::int32_t> typedExtinction(typedTree, typedLevel.second);\n"
"    auto typedExtinctionFiltered = typedExtinction.filtering(1024);\n"
"    static_assert(std::is_same_v<decltype(typedExtinctionFiltered), mmcfilters::ImageInt32Ptr>);\n"
"    assert(typedExtinctionFiltered->getSize() == intImage->getSize());\n"
"    mmcfilters::UltimateAttributeOpening<std::int32_t> typedUao(typedTree, typedLevel.second);\n"
"    typedUao.execute(intImage->getSize());\n"
"    auto typedMaxContrast = typedUao.getMaxContrastImage();\n"
"    static_assert(std::is_same_v<decltype(typedMaxContrast), mmcfilters::ImageInt32Ptr>);\n"
"    assert(typedMaxContrast->getSize() == intImage->getSize());\n"
"    mmcfilters::MSERComputer<std::int32_t> typedMser(typedTree);\n"
"    auto typedMserFlags = typedMser.computeMSER(mmcfilters::AltitudeDiff<std::int32_t>{1});\n"
"    assert(typedMserFlags.size() == static_cast<std::size_t>(typedTree.topology().getNumInternalNodeSlots()));\n"
"    typedUao.executeWithMSER(intImage->getSize(), mmcfilters::AltitudeDiff<std::int32_t>{1});\n"
"    auto typedAssociated = typedUao.getAssociatedImage();\n"
"    static_assert(std::is_same_v<decltype(typedAssociated), mmcfilters::ImageInt32Ptr>);\n"
"    assert(typedAssociated->getSize() == intImage->getSize());\n"
"    mmcfilters::adjust::CasfComponentTrees<std::int32_t> typedCasf(\n"
"        intImage,\n"
"        mmcfilters::adjust::CasfComponentTreesAttribute::AREA);\n"
	"    auto typedCasfFiltered = typedCasf.filter({1.0});\n"
	"    static_assert(std::is_same_v<decltype(typedCasfFiltered), mmcfilters::ImageInt32Ptr>);\n"
	"    assert(typedCasfFiltered->getSize() == intImage->getSize());\n"
	"    auto typedSdrt = mmcfilters::MorphologicalTreeFactory::createSelfDualResidualTree(intImage);\n"
	"    static_assert(std::is_same_v<decltype(typedSdrt), mmcfilters::WeightedMorphologicalTree<std::int32_t>>);\n"
	"    auto typedSdrtReconstruction = typedSdrt.reconstructionImage();\n"
	"    static_assert(std::is_same_v<decltype(typedSdrtReconstruction), mmcfilters::ImageInt32Ptr>);\n"
	"    assert(typedSdrtReconstruction->getSize() == intImage->getSize());\n"
	"    auto typedProjectionAttrs = mmcfilters::AttributeComputation::computeAttributes(\n"
	"        typedTree,\n"
	"        std::vector<mmcfilters::AttributeOrGroup>{mmcfilters::Attribute::LEVEL, mmcfilters::Attribute::VOLUME});\n"
	"    auto typedProjectedHigra = mmcfilters::AttributeComputation::projectNodeValuesToExportedHigra(\n"
	"        typedTree,\n"
	"        typedProjectionAttrs.first,\n"
	"        typedProjectionAttrs.second);\n"
	"    assert(typedProjectedHigra.size() == static_cast<std::size_t>(\n"
	"        (typedTree.topology().getNumTotalProperParts() + typedTree.topology().getNumNodes()) * typedProjectionAttrs.first.NUM_ATTRIBUTES));\n"
"\n"
"    std::array<float, 4> floatPixels{1.0f, 2.5f, 3.0f, 4.5f};\n"
"    auto floatImage = mmcfilters::ImageFloat::fromExternal(floatPixels.data(), 2, 2);\n"
"    auto floatTree = mmcfilters::MorphologicalTreeFactory::createMaxTree(floatImage);\n"
"    static_assert(std::is_same_v<decltype(floatTree), mmcfilters::WeightedMorphologicalTree<float>>);\n"
"    auto floatReconstruction = floatTree.reconstructionImage();\n"
"    static_assert(std::is_same_v<decltype(floatReconstruction), mmcfilters::ImageFloatPtr>);\n"
"    assert(floatReconstruction->getSize() == floatImage->getSize());\n"
"    for (int i = 0; i < floatImage->getSize(); ++i) {\n"
"        assert((*floatReconstruction)[i] == (*floatImage)[i]);\n"
"    }\n"
"    auto floatView = floatTree.asView();\n"
"    static_assert(std::is_same_v<decltype(floatView), mmcfilters::WeightedTreeView<float>>);\n"
"    mmcfilters::TreeAltitudeAlgorithms::validateAltitudeBufferShape(floatView.topology(), floatView.altitude());\n"
"    auto [floatParent, floatAltitude] = floatTree.exportHigraHierarchy();\n"
"    static_assert(std::is_same_v<decltype(floatAltitude), std::vector<float>>);\n"
"    assert(floatParent.size() == floatAltitude.size());\n"
"    auto floatLevel = mmcfilters::AttributeComputation::computeSingleAttribute(\n"
"        floatTree,\n"
"        mmcfilters::Attribute::LEVEL);\n"
"    auto floatLevelDelta = mmcfilters::AttributeComputation::computeSingleAttributeWithDelta(\n"
"        floatTree,\n"
"        mmcfilters::Attribute::LEVEL,\n"
"        mmcfilters::AltitudeDiff<float>{0.5f},\n"
"        1);\n"
"    assert(floatLevelDelta.first.NUM_ATTRIBUTES == 3);\n"
"    assert(floatLevelDelta.second.size() == static_cast<std::size_t>(\n"
"        floatTree.topology().getNumInternalNodeSlots() * floatLevelDelta.first.NUM_ATTRIBUTES));\n"
"    std::vector<bool> keepAllFloat(floatTree.topology().getNumInternalNodeSlots(), true);\n"
"    mmcfilters::AttributeFilters<float> floatFilters(floatTree);\n"
"    auto floatFiltered = floatFilters.filteringByDirectRule(keepAllFloat);\n"
"    static_assert(std::is_same_v<decltype(floatFiltered), mmcfilters::ImageFloatPtr>);\n"
"    assert(floatFiltered->getSize() == floatImage->getSize());\n"
"    mmcfilters::ExtinctionValues<float> floatExtinction(floatTree, floatLevel.second);\n"
"    auto floatExtinctionFiltered = floatExtinction.filtering(1024);\n"
"    static_assert(std::is_same_v<decltype(floatExtinctionFiltered), mmcfilters::ImageFloatPtr>);\n"
"    assert(floatExtinctionFiltered->getSize() == floatImage->getSize());\n"
"    mmcfilters::UltimateAttributeOpening<float> floatUao(floatTree, floatLevel.second);\n"
"    floatUao.execute(floatImage->getSize());\n"
"    auto floatMaxContrast = floatUao.getMaxContrastImage();\n"
"    static_assert(std::is_same_v<decltype(floatMaxContrast), mmcfilters::ImageFloatPtr>);\n"
"    assert(floatMaxContrast->getSize() == floatImage->getSize());\n"
"    mmcfilters::MSERComputer<float> floatMser(floatTree);\n"
"    auto floatMserFlags = floatMser.computeMSER(0.5f);\n"
"    assert(floatMserFlags.size() == static_cast<std::size_t>(floatTree.topology().getNumInternalNodeSlots()));\n"
"    floatUao.executeWithMSER(floatImage->getSize(), 0.5f);\n"
"    mmcfilters::adjust::CasfComponentTrees<float> floatCasf(\n"
"        floatImage,\n"
"        mmcfilters::adjust::CasfComponentTreesAttribute::AREA);\n"
	"    auto floatCasfFiltered = floatCasf.filter({1.0});\n"
	"    static_assert(std::is_same_v<decltype(floatCasfFiltered), mmcfilters::ImageFloatPtr>);\n"
	"    assert(floatCasfFiltered->getSize() == floatImage->getSize());\n"
	"    auto floatSdrt = mmcfilters::MorphologicalTreeFactory::createSelfDualResidualTree(floatImage);\n"
	"    static_assert(std::is_same_v<decltype(floatSdrt), mmcfilters::WeightedMorphologicalTree<float>>);\n"
	"    auto floatSdrtReconstruction = floatSdrt.reconstructionImage();\n"
	"    static_assert(std::is_same_v<decltype(floatSdrtReconstruction), mmcfilters::ImageFloatPtr>);\n"
	"    assert(floatSdrtReconstruction->getSize() == floatImage->getSize());\n"
	"    auto floatProjectionAttrs = mmcfilters::AttributeComputation::computeAttributes(\n"
	"        floatTree,\n"
	"        std::vector<mmcfilters::AttributeOrGroup>{mmcfilters::Attribute::LEVEL, mmcfilters::Attribute::VOLUME});\n"
	"    auto floatProjectedHigra = mmcfilters::AttributeComputation::projectNodeValuesToExportedHigra(\n"
	"        floatTree,\n"
	"        floatProjectionAttrs.first,\n"
	"        floatProjectionAttrs.second);\n"
	"    assert(floatProjectedHigra.size() == static_cast<std::size_t>(\n"
	"        (floatTree.topology().getNumTotalProperParts() + floatTree.topology().getNumNodes()) * floatProjectionAttrs.first.NUM_ATTRIBUTES));\n"
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
