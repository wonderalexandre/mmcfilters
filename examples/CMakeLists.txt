set(MMCFILTERS_EXAMPLES_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
set(MMCFILTERS_INTERNAL_EXAMPLES_ROOT ${MMCFILTERS_EXAMPLES_ROOT}/internal)

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

add_executable(mmcfilters_local_events_bitquads_benchmark
    ${MMCFILTERS_INTERNAL_EXAMPLES_ROOT}/local_events_bitquads_benchmark.cpp
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

target_link_libraries(mmcfilters_internal_edt_dift PRIVATE mmcfilters::core)
target_link_libraries(mmcfilters_example_editing_api PRIVATE mmcfilters::core)
target_link_libraries(mmcfilters_contour_benchmark PRIVATE mmcfilters::core)
target_link_libraries(mmcfilters_local_events_bitquads_benchmark PRIVATE mmcfilters::core)
target_link_libraries(mmcfilters_tos_bitquad_projection_export PRIVATE mmcfilters::core)
target_link_libraries(mmcfilters_tos_bitquad_naive_validation PRIVATE mmcfilters::core)
target_link_libraries(mmcfilters_maxdist_benchmark PRIVATE mmcfilters::core)
target_link_libraries(mmcfilters_altitude_span_attribute_benchmark PRIVATE mmcfilters::core)
target_link_libraries(mmcfilters_typed_altitude_benchmark PRIVATE mmcfilters::core)
target_include_directories(mmcfilters_internal_edt_dift PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${CMAKE_CURRENT_SOURCE_DIR}/../external/stb
)
target_include_directories(mmcfilters_example_editing_api PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)
target_include_directories(mmcfilters_contour_benchmark PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${CMAKE_CURRENT_SOURCE_DIR}/../external/stb
)
target_include_directories(mmcfilters_local_events_bitquads_benchmark PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${CMAKE_CURRENT_SOURCE_DIR}/../external/stb
)
target_include_directories(mmcfilters_tos_bitquad_projection_export PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${CMAKE_CURRENT_SOURCE_DIR}/../external/stb
)
target_include_directories(mmcfilters_tos_bitquad_naive_validation PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${CMAKE_CURRENT_SOURCE_DIR}/../external/stb
)
target_include_directories(mmcfilters_maxdist_benchmark PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)
target_include_directories(mmcfilters_altitude_span_attribute_benchmark PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)
target_include_directories(mmcfilters_typed_altitude_benchmark PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)
target_compile_features(mmcfilters_internal_edt_dift PRIVATE cxx_std_20)
target_compile_features(mmcfilters_example_editing_api PRIVATE cxx_std_20)
target_compile_features(mmcfilters_contour_benchmark PRIVATE cxx_std_20)
target_compile_features(mmcfilters_local_events_bitquads_benchmark PRIVATE cxx_std_20)
target_compile_features(mmcfilters_tos_bitquad_projection_export PRIVATE cxx_std_20)
target_compile_features(mmcfilters_tos_bitquad_naive_validation PRIVATE cxx_std_20)
target_compile_features(mmcfilters_maxdist_benchmark PRIVATE cxx_std_20)
target_compile_features(mmcfilters_altitude_span_attribute_benchmark PRIVATE cxx_std_20)
target_compile_features(mmcfilters_typed_altitude_benchmark PRIVATE cxx_std_20)

source_group(TREE ${MMCFILTERS_EXAMPLES_ROOT} FILES
    ${MMCFILTERS_EXAMPLES_ROOT}/contour_benchmark.cpp
    ${MMCFILTERS_INTERNAL_EXAMPLES_ROOT}/edt_dift_example.cpp
    ${MMCFILTERS_INTERNAL_EXAMPLES_ROOT}/local_events_bitquads_benchmark.cpp
    ${MMCFILTERS_INTERNAL_EXAMPLES_ROOT}/tos_bitquad_projection_export.cpp
    ${MMCFILTERS_INTERNAL_EXAMPLES_ROOT}/tos_bitquad_naive_validation.cpp
    ${MMCFILTERS_EXAMPLES_ROOT}/editing_api_example.cpp
    ${MMCFILTERS_EXAMPLES_ROOT}/altitude_span_attribute_benchmark.cpp
    ${MMCFILTERS_EXAMPLES_ROOT}/maxdist_benchmark.cpp
    ${MMCFILTERS_EXAMPLES_ROOT}/typed_altitude_benchmark.cpp
)
