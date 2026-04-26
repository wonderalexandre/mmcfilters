# -----------------------------------------------------------------------------
# unit-tests/CMakeLists.txt
# -----------------------------------------------------------------------------
# Testes unitários focados na API atual da MorphologicalTree.
# Esta pasta concentra a suíte principal do projeto.

cmake_minimum_required(VERSION 3.14)
project(mmcfilters_unit_tests LANGUAGES CXX)

set(MMCFILTERS_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../mmcfilters")

if(NOT TARGET mmcfilters::core)
    add_subdirectory(${MMCFILTERS_SOURCE_DIR} ${CMAKE_BINARY_DIR}/mmcfilters_build)
endif()

if(NOT TARGET mmcfilters::core)
    message(FATAL_ERROR "O alvo mmcfilters::core não foi encontrado.")
endif()

function(mmcfilters_add_unit_test target)
    add_executable(${target} ${ARGN})
    target_link_libraries(${target} PRIVATE mmcfilters::core)
    target_compile_features(${target} PRIVATE cxx_std_20)
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
    target_include_directories(${target}
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${CMAKE_CURRENT_SOURCE_DIR}/support
            ${CMAKE_CURRENT_SOURCE_DIR}/..
    )
    add_test(NAME ${target} COMMAND ${target})
endfunction()

mmcfilters_add_unit_test(unit_morphological_tree_build
    trees/test_morphological_tree_build.cpp
)

mmcfilters_add_unit_test(unit_morphological_tree_traversals
    trees/test_morphological_tree_traversals.cpp
)

mmcfilters_add_unit_test(unit_morphological_tree_tos
    trees/test_morphological_tree_tos.cpp
)

mmcfilters_add_unit_test(unit_morphological_tree_topology
    trees/test_morphological_tree_topology.cpp
)

mmcfilters_add_unit_test(unit_morphological_tree_proper_parts
    trees/test_morphological_tree_proper_parts.cpp
)

mmcfilters_add_unit_test(unit_morphological_tree_mutations
    trees/test_morphological_tree_mutations.cpp
)

mmcfilters_add_unit_test(unit_tree_editor
    trees/test_tree_editor.cpp
)

mmcfilters_add_unit_test(unit_edit_api_contracts
    trees/test_edit_api_contracts.cpp
)

mmcfilters_add_unit_test(unit_weighted_morphological_tree
    trees/test_weighted_morphological_tree.cpp
)

mmcfilters_add_unit_test(unit_component_tree_adjust
    trees/test_component_tree_adjust.cpp
)

mmcfilters_add_unit_test(unit_morphological_tree_invariants
    trees/test_morphological_tree_invariants.cpp
)

mmcfilters_add_unit_test(unit_morphological_tree_internals
    trees/test_morphological_tree_internals.cpp
)

mmcfilters_add_unit_test(unit_attributes_on_morphological_tree
    attributes/test_attributes_on_morphological_tree.cpp
)

mmcfilters_add_unit_test(unit_attribute_plumbing
    attributes/test_attribute_plumbing.cpp
)

mmcfilters_add_unit_test(unit_secondary_attribute_consumers
    attributes/test_secondary_attribute_consumers.cpp
)

mmcfilters_add_unit_test(unit_contours_on_morphological_tree
    contours/test_contours_on_morphological_tree.cpp
)

mmcfilters_add_unit_test(unit_extinction_values_on_morphological_tree
    filters/test_extinction_values_on_morphological_tree.cpp
)

mmcfilters_add_unit_test(unit_attribute_filters_and_uao
    filters/test_attribute_filters_and_uao.cpp
)

if(TARGET mmcfilters)
    find_package(Python3 COMPONENTS Interpreter REQUIRED)
    set(MMCFILTERS_PYTHON_TEST_EXECUTABLE "${Python3_EXECUTABLE}")

    execute_process(
        COMMAND "${MMCFILTERS_PYTHON_TEST_EXECUTABLE}" -c "import numpy"
        RESULT_VARIABLE MMCFILTERS_PYTHON_HAS_NUMPY
        OUTPUT_QUIET
        ERROR_QUIET
    )

    if(NOT MMCFILTERS_PYTHON_HAS_NUMPY EQUAL 0 AND EXISTS "/opt/anaconda3/bin/python3.12")
        execute_process(
            COMMAND /opt/anaconda3/bin/python3.12 -c "import numpy"
            RESULT_VARIABLE MMCFILTERS_CONDA_PYTHON_HAS_NUMPY
            OUTPUT_QUIET
            ERROR_QUIET
        )
        if(MMCFILTERS_CONDA_PYTHON_HAS_NUMPY EQUAL 0)
            set(MMCFILTERS_PYTHON_TEST_EXECUTABLE "/opt/anaconda3/bin/python3.12")
        endif()
    endif()

    add_test(
        NAME unit_python_nodeid_api
        COMMAND ${MMCFILTERS_PYTHON_TEST_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/python/test_python_nodeid_api.py
            ${CMAKE_BINARY_DIR}
    )

    add_test(
        NAME unit_python_filter_wrappers
        COMMAND ${MMCFILTERS_PYTHON_TEST_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/python/test_python_filter_wrappers.py
            ${CMAKE_BINARY_DIR}
    )
endif()
