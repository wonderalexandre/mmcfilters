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
    target_compile_options(${target}
        PRIVATE
            $<$<CXX_COMPILER_ID:MSVC>:/W4>
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall>
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wextra>
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wpedantic>
    )
    target_include_directories(${target}
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${CMAKE_CURRENT_SOURCE_DIR}/support
            ${CMAKE_CURRENT_SOURCE_DIR}/..
    )
    add_test(NAME ${target} COMMAND ${target})
endfunction()

mmcfilters_add_unit_test(unit_altitude_traits
    utils/test_altitude_traits.cpp
)

mmcfilters_add_unit_test(unit_regular_grid_adjacency_2d
    utils/test_regular_grid_adjacency_2d.cpp
)

mmcfilters_add_unit_test(unit_validation_boundaries
    utils/test_validation_boundaries.cpp
)

find_package(Python3 COMPONENTS Interpreter REQUIRED)
add_test(
    NAME unit_validation_boundary_audit
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/audit_validation_boundaries.py
        ${CMAKE_CURRENT_SOURCE_DIR}/..
)

set(MMCFILTERS_API_REFACTORING_AUDIT_ARGS
    --root ${CMAKE_CURRENT_SOURCE_DIR}/..
)
if(MMCFILTERS_REFACTORING_METADATA_AVAILABLE)
    list(APPEND MMCFILTERS_API_REFACTORING_AUDIT_ARGS
        --status ${MMCFILTERS_REFACTORING_METADATA_DIR}/api-migration-status.json
        --vocabulary ${MMCFILTERS_REFACTORING_METADATA_DIR}/api-vocabulary-audit.json
        --contextual ${MMCFILTERS_REFACTORING_METADATA_DIR}/api-contextual-vocabulary-review.json
        --context-report ${CMAKE_BINARY_DIR}/refactoring-contextual-vocabulary-report.json
    )
endif()
add_test(
    NAME unit_api_refactoring_audit
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/audit_api_refactoring.py
        ${MMCFILTERS_API_REFACTORING_AUDIT_ARGS}
)
add_test(
    NAME unit_api_refactoring_audit_self_test
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/audit_api_refactoring.py
        --self-test
)

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

mmcfilters_add_unit_test(unit_valued_morphological_tree
    trees/test_valued_morphological_tree.cpp
)

mmcfilters_add_unit_test(unit_valued_altitude_templates
    trees/test_valued_altitude_templates.cpp
)

mmcfilters_add_unit_test(unit_hierarchy_saliency_map
    trees/test_hierarchy_saliency_map.cpp
)

mmcfilters_add_unit_test(unit_min_max_residual_tree_factory
    trees/sdrt/test_min_max_residual_tree_factory.cpp
)

mmcfilters_add_unit_test(unit_shape_space_saliency
    trees/test_shape_space_saliency.cpp
)

mmcfilters_add_unit_test(unit_owner_state_contracts
    trees/test_owner_state_contracts.cpp
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

mmcfilters_add_unit_test(unit_finite_window_local_attribute_computations
    attributes/test_finite_window_local_attribute_computations.cpp
)

mmcfilters_add_unit_test(unit_attribute_plumbing
    attributes/test_attribute_plumbing.cpp
)

mmcfilters_add_unit_test(unit_public_attribute_api
    attributes/test_public_attribute_api.cpp
)

mmcfilters_add_unit_test(unit_attribute_unit_values
    attributes/test_attribute_unit_values.cpp
)

mmcfilters_add_unit_test(unit_all_attributes_numeric_validation
    attributes/test_all_attributes_numeric_validation.cpp
)
target_compile_definitions(unit_all_attributes_numeric_validation
    PRIVATE
        MMCFILTERS_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../dat"
)

mmcfilters_add_unit_test(unit_maxdist_support
    attributes/test_maxdist_support.cpp
)

mmcfilters_add_unit_test(unit_secondary_attribute_consumers
    attributes/test_secondary_attribute_consumers.cpp
)

mmcfilters_add_unit_test(unit_contours_on_morphological_tree
    contours/test_contours_on_morphological_tree.cpp
)

mmcfilters_add_unit_test(unit_contour_traces_on_morphological_tree
    contours/test_contour_traces_on_morphological_tree.cpp
)
target_compile_definitions(unit_contour_traces_on_morphological_tree
    PRIVATE
        MMCFILTERS_ENABLE_CONTOUR_TRACE_PROFILE
)

mmcfilters_add_unit_test(unit_extinction_values_on_morphological_tree
    filters/test_extinction_values_on_morphological_tree.cpp
)

mmcfilters_add_unit_test(unit_attribute_filters_and_uao
    filters/test_attribute_filters_and_uao.cpp
)

mmcfilters_add_unit_test(unit_stability_computers
    filters/test_stability_computers.cpp
)

set(MMCFILTERS_INSTALLED_CONSUMER_TEST_ARGS
    -DMMCFILTERS_BUILD_DIR=${CMAKE_BINARY_DIR}
    -DMMCFILTERS_PROJECT_ROOT=${CMAKE_CURRENT_SOURCE_DIR}/..
    -DMMCFILTERS_PYTHON_EXECUTABLE=${Python3_EXECUTABLE}
)
if(CMAKE_CONFIGURATION_TYPES)
    list(APPEND MMCFILTERS_INSTALLED_CONSUMER_TEST_ARGS
        -DMMCFILTERS_CTEST_CONFIG=$<CONFIG>
    )
endif()

add_test(
    NAME unit_installed_consumer
    COMMAND
        ${CMAKE_COMMAND}
            ${MMCFILTERS_INSTALLED_CONSUMER_TEST_ARGS}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/installed_consumer_test.cmake)

if(TARGET mmcfilters)
    if(DEFINED Python_EXECUTABLE)
        set(MMCFILTERS_PYTHON_TEST_EXECUTABLE "${Python_EXECUTABLE}")
    else()
        find_package(Python3 COMPONENTS Interpreter REQUIRED)
        set(MMCFILTERS_PYTHON_TEST_EXECUTABLE "${Python3_EXECUTABLE}")
    endif()

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
        NAME unit_python_import_smoke
        COMMAND ${CMAKE_COMMAND} -E env
            "MMCFILTERS_BUILD_PYTHON=${CMAKE_BINARY_DIR}/python"
            "PYTHONPATH=${CMAKE_BINARY_DIR}/python"
            ${MMCFILTERS_PYTHON_TEST_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/python/test_python_import_smoke.py
    )

    set(MMCFILTERS_PYTHON_API_STYLE_AUDIT_ARGS
        --build-dir ${CMAKE_BINARY_DIR}
    )
    if(EXISTS "${MMCFILTERS_REFACTORING_METADATA_DIR}/python-public-api-style-baseline.json")
        list(APPEND MMCFILTERS_PYTHON_API_STYLE_AUDIT_ARGS
            --baseline ${MMCFILTERS_REFACTORING_METADATA_DIR}/python-public-api-style-baseline.json
        )
    endif()
    add_test(
        NAME unit_python_public_api_style_audit
        COMMAND ${MMCFILTERS_PYTHON_TEST_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/python/audit_python_public_api_style.py
            ${MMCFILTERS_PYTHON_API_STYLE_AUDIT_ARGS}
    )

    add_test(
        NAME unit_python_nodeid_api
        COMMAND ${MMCFILTERS_PYTHON_TEST_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/python/test_python_nodeid_api.py
            ${CMAKE_BINARY_DIR}
    )

    add_test(
        NAME unit_python_attribute_numeric_validation
        COMMAND ${CMAKE_COMMAND} -E env
            "MMCFILTERS_TEST_DATA_DIR=${CMAKE_CURRENT_SOURCE_DIR}/../dat"
            ${MMCFILTERS_PYTHON_TEST_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/python/test_python_attribute_numeric_validation.py
            ${CMAKE_BINARY_DIR}
    )

    add_test(
        NAME unit_python_filter_wrappers
        COMMAND ${MMCFILTERS_PYTHON_TEST_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/python/test_python_filter_wrappers.py
            ${CMAKE_BINARY_DIR}
    )

    add_test(
        NAME unit_python_hierarchy_saliency_map
        COMMAND ${MMCFILTERS_PYTHON_TEST_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/python/test_python_hierarchy_saliency_map.py
            ${CMAKE_BINARY_DIR}
    )

    add_test(
        NAME unit_python_shape_space_saliency
        COMMAND ${MMCFILTERS_PYTHON_TEST_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/python/test_python_shape_space_saliency.py
            ${CMAKE_BINARY_DIR}
    )

    add_test(
        NAME unit_python_contour_traces
        COMMAND ${MMCFILTERS_PYTHON_TEST_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/python/test_python_contour_traces.py
            ${CMAKE_BINARY_DIR}
    )
endif()
