set(MMCFILTERS_EXAMPLES_ROOT ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(mmcfilters_example_edt_dift
    ${MMCFILTERS_EXAMPLES_ROOT}/edt_dift_example.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../external/stb/stb_image.cpp
)

add_executable(mmcfilters_example_editing_api
    ${MMCFILTERS_EXAMPLES_ROOT}/editing_api_example.cpp
)

add_executable(mmcfilters_contour_benchmark
    ${MMCFILTERS_EXAMPLES_ROOT}/contour_benchmark.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../external/stb/stb_image.cpp
)

target_link_libraries(mmcfilters_example_edt_dift PRIVATE mmcfilters::core)
target_link_libraries(mmcfilters_example_editing_api PRIVATE mmcfilters::core)
target_link_libraries(mmcfilters_contour_benchmark PRIVATE mmcfilters::core)
target_include_directories(mmcfilters_example_edt_dift PRIVATE
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
target_compile_features(mmcfilters_example_edt_dift PRIVATE cxx_std_20)
target_compile_features(mmcfilters_example_editing_api PRIVATE cxx_std_20)
target_compile_features(mmcfilters_contour_benchmark PRIVATE cxx_std_20)

source_group(TREE ${MMCFILTERS_EXAMPLES_ROOT} FILES
    ${MMCFILTERS_EXAMPLES_ROOT}/contour_benchmark.cpp
    ${MMCFILTERS_EXAMPLES_ROOT}/edt_dift_example.cpp
    ${MMCFILTERS_EXAMPLES_ROOT}/editing_api_example.cpp
)
