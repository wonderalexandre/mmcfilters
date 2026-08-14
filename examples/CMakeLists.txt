set(MMCFILTERS_EXAMPLES_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
set(MMCFILTERS_INTERNAL_EXAMPLES_ROOT ${MMCFILTERS_EXAMPLES_ROOT}/internal)

function(mmcfilters_configure_example target)
    target_link_libraries(${target} PRIVATE mmcfilters::core)
    target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
    target_compile_features(${target} PRIVATE cxx_std_20)
endfunction()

function(mmcfilters_configure_stb_example target)
    mmcfilters_configure_example(${target})
    target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../external/stb)
endfunction()

# Internal examples may include detail headers to exercise/debug implementation
# kernels directly; public examples should stay on stable headers.
add_executable(mmcfilters_internal_edt_dift
    ${MMCFILTERS_INTERNAL_EXAMPLES_ROOT}/edt_dift_example.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../external/stb/stb_image.cpp
)

add_executable(mmcfilters_example_editing_api
    ${MMCFILTERS_EXAMPLES_ROOT}/editing_api_example.cpp
)

add_executable(mmcfilters_contour_benchmark
    ${MMCFILTERS_EXAMPLES_ROOT}/contour_benchmark.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../external/stb/stb_image.cpp
)

add_executable(mmcfilters_contour_trace_benchmark
    ${MMCFILTERS_EXAMPLES_ROOT}/contour_trace_benchmark.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../external/stb/stb_image.cpp
)

add_executable(mmcfilters_contour_trace_profile
    ${MMCFILTERS_INTERNAL_EXAMPLES_ROOT}/contour_trace_profile.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../external/stb/stb_image.cpp
)

add_executable(mmcfilters_finite_window_bitquad_benchmark
    ${MMCFILTERS_INTERNAL_EXAMPLES_ROOT}/finite_window_bitquad_benchmark.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../external/stb/stb_image.cpp
)

add_executable(mmcfilters_tos_bitquad_projection_export
    ${MMCFILTERS_INTERNAL_EXAMPLES_ROOT}/tos_bitquad_projection_export.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../external/stb/stb_image.cpp
)

add_executable(mmcfilters_tos_bitquad_naive_validation
    ${MMCFILTERS_INTERNAL_EXAMPLES_ROOT}/tos_bitquad_naive_validation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../external/stb/stb_image.cpp
)

add_executable(mmcfilters_maxdist_benchmark
    ${MMCFILTERS_EXAMPLES_ROOT}/maxdist_benchmark.cpp
)

add_executable(mmcfilters_altitude_span_attribute_benchmark
    ${MMCFILTERS_EXAMPLES_ROOT}/altitude_span_attribute_benchmark.cpp
)

add_executable(mmcfilters_typed_altitude_benchmark
    ${MMCFILTERS_EXAMPLES_ROOT}/typed_altitude_benchmark.cpp
)

mmcfilters_configure_stb_example(mmcfilters_internal_edt_dift)
mmcfilters_configure_example(mmcfilters_example_editing_api)
mmcfilters_configure_stb_example(mmcfilters_contour_benchmark)
mmcfilters_configure_stb_example(mmcfilters_contour_trace_benchmark)
mmcfilters_configure_stb_example(mmcfilters_contour_trace_profile)
mmcfilters_configure_stb_example(mmcfilters_finite_window_bitquad_benchmark)
mmcfilters_configure_stb_example(mmcfilters_tos_bitquad_projection_export)
mmcfilters_configure_stb_example(mmcfilters_tos_bitquad_naive_validation)
mmcfilters_configure_example(mmcfilters_maxdist_benchmark)
mmcfilters_configure_example(mmcfilters_altitude_span_attribute_benchmark)
mmcfilters_configure_example(mmcfilters_typed_altitude_benchmark)

source_group(TREE ${MMCFILTERS_EXAMPLES_ROOT} FILES
    ${MMCFILTERS_EXAMPLES_ROOT}/contour_benchmark.cpp
    ${MMCFILTERS_EXAMPLES_ROOT}/contour_trace_benchmark.cpp
    ${MMCFILTERS_INTERNAL_EXAMPLES_ROOT}/contour_trace_profile.cpp
    ${MMCFILTERS_INTERNAL_EXAMPLES_ROOT}/edt_dift_example.cpp
    ${MMCFILTERS_INTERNAL_EXAMPLES_ROOT}/finite_window_bitquad_benchmark.cpp
    ${MMCFILTERS_INTERNAL_EXAMPLES_ROOT}/tos_bitquad_projection_export.cpp
    ${MMCFILTERS_INTERNAL_EXAMPLES_ROOT}/tos_bitquad_naive_validation.cpp
    ${MMCFILTERS_EXAMPLES_ROOT}/editing_api_example.cpp
    ${MMCFILTERS_EXAMPLES_ROOT}/altitude_span_attribute_benchmark.cpp
    ${MMCFILTERS_EXAMPLES_ROOT}/maxdist_benchmark.cpp
    ${MMCFILTERS_EXAMPLES_ROOT}/typed_altitude_benchmark.cpp
)
