cmake_minimum_required(VERSION 3.20)

if(NOT DEFINED MMCFILTERS_BUILD_DIR)
    message(FATAL_ERROR "MMCFILTERS_BUILD_DIR is required")
endif()

set(work_dir "${MMCFILTERS_BUILD_DIR}/installed-consumer-test")
set(prefix "${work_dir}/prefix")
set(consumer_source_dir "${work_dir}/consumer")
set(consumer_build_dir "${work_dir}/consumer-build")

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
"#include <mmcfilters/trees/WeightedMorphologicalTree.hpp>\n"
"#include <mmcfilters/utils/Image.hpp>\n"
"\n"
"int main()\n"
"{\n"
"    std::array<std::uint8_t, 4> pixels{1, 2, 3, 4};\n"
"    auto image = mmcfilters::ImageUInt8::fromExternal(pixels.data(), 2, 2);\n"
"    auto tree = mmcfilters::WeightedMorphologicalTree::createComponentTree(image, true);\n"
"    assert(tree.topology().getNumNodes() > 0);\n"
"    auto reconstructed = tree.reconstructionImage();\n"
"    assert(reconstructed->getSize() == 4);\n"
"    return 0;\n"
"}\n")

execute_process(
    COMMAND
        "${CMAKE_COMMAND}"
        -S "${consumer_source_dir}"
        -B "${consumer_build_dir}"
        "-DCMAKE_PREFIX_PATH=${prefix}"
    RESULT_VARIABLE configure_result
    OUTPUT_VARIABLE configure_output
    ERROR_VARIABLE configure_error)

if(NOT configure_result EQUAL 0)
    message(FATAL_ERROR
        "Failed to configure installed mmcfilters consumer.\n"
        "stdout:\n${configure_output}\n"
        "stderr:\n${configure_error}")
endif()

set(consumer_build_command "${CMAKE_COMMAND}" --build "${consumer_build_dir}" --parallel)

if(DEFINED MMCFILTERS_CTEST_CONFIG AND NOT "${MMCFILTERS_CTEST_CONFIG}" STREQUAL "")
    list(APPEND consumer_build_command --config "${MMCFILTERS_CTEST_CONFIG}")
endif()

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

set(consumer_executable "${consumer_build_dir}/consumer")

if(DEFINED MMCFILTERS_CTEST_CONFIG AND NOT "${MMCFILTERS_CTEST_CONFIG}" STREQUAL "")
    set(consumer_executable "${consumer_build_dir}/${MMCFILTERS_CTEST_CONFIG}/consumer")
endif()

if(WIN32)
    set(consumer_executable "${consumer_executable}.exe")
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
