add_executable(mmcfilters_tree_validation_benchmark
    ${CMAKE_CURRENT_SOURCE_DIR}/tree_validation_benchmark.cpp
)

add_executable(mmcfilters_adjacency_iteration_benchmark
    ${CMAKE_CURRENT_SOURCE_DIR}/adjacency_iteration_benchmark.cpp
)

add_executable(mmcfilters_min_max_residual_tree_benchmark
    ${CMAKE_CURRENT_SOURCE_DIR}/min_max_residual_tree_benchmark.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../external/stb/stb_image.cpp
)

add_executable(mmcfilters_tree_construction_comparison_benchmark
    ${CMAKE_CURRENT_SOURCE_DIR}/tree_construction_comparison_benchmark.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../external/stb/stb_image.cpp
)

foreach(target IN ITEMS
        mmcfilters_tree_validation_benchmark
        mmcfilters_adjacency_iteration_benchmark
        mmcfilters_min_max_residual_tree_benchmark
        mmcfilters_tree_construction_comparison_benchmark)
    target_link_libraries(${target} PRIVATE mmcfilters::core)
    target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
    target_compile_features(${target} PRIVATE cxx_std_20)
endforeach()

target_include_directories(mmcfilters_min_max_residual_tree_benchmark
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../external/stb)

target_include_directories(mmcfilters_tree_construction_comparison_benchmark
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../external/stb)

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/tree_validation_benchmark.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/adjacency_iteration_benchmark.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/min_max_residual_tree_benchmark.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tree_construction_comparison_benchmark.cpp
)
