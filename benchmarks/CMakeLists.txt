set(MMCFILTERS_STB_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/../external/stb/stb_image.cpp)

# stb is vendored third-party code. Keep platform deprecation diagnostics out
# of project builds without weakening warnings for mmcfilters sources.
if(MSVC)
    set_source_files_properties(${MMCFILTERS_STB_SOURCE} PROPERTIES COMPILE_OPTIONS "/wd4996")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    set_source_files_properties(${MMCFILTERS_STB_SOURCE} PROPERTIES COMPILE_OPTIONS "-Wno-deprecated-declarations")
endif()

add_executable(mmcfilters_tree_validation_benchmark
    ${CMAKE_CURRENT_SOURCE_DIR}/tree_validation_benchmark.cpp
)

add_executable(mmcfilters_adjacency_iteration_benchmark
    ${CMAKE_CURRENT_SOURCE_DIR}/adjacency_iteration_benchmark.cpp
)

add_executable(mmcfilters_min_max_residual_tree_benchmark
    ${CMAKE_CURRENT_SOURCE_DIR}/min_max_residual_tree_benchmark.cpp
    ${MMCFILTERS_STB_SOURCE}
)

add_executable(mmcfilters_tree_construction_comparison_benchmark
    ${CMAKE_CURRENT_SOURCE_DIR}/tree_construction_comparison_benchmark.cpp
    ${MMCFILTERS_STB_SOURCE}
)

add_executable(mmcfilters_scientific_pipeline_benchmark
    ${CMAKE_CURRENT_SOURCE_DIR}/scientific_pipeline_benchmark.cpp
)

add_executable(mmcfilters_validation_sensitive_algorithms_benchmark
    ${CMAKE_CURRENT_SOURCE_DIR}/validation_sensitive_algorithms_benchmark.cpp
)

set(MMCFILTERS_API_BENCHMARK_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/api/api_benchmark.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/api/ConstructionScenarios.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/api/AttributeScenarios.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/api/FilterScenarios.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/api/InteroperabilityScenarios.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/api/EditingScenarios.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/api/CasfScenarios.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/api/PipelineScenarios.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/api/WorkloadManifest.cpp
)

add_executable(mmcfilters_api_benchmark
    ${MMCFILTERS_API_BENCHMARK_SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/support/stb_image_read.cpp
)

target_include_directories(mmcfilters_api_benchmark
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../external/stb)

if(MMCFILTERS_BUILD_TESTS)
    add_test(
        NAME benchmark_api_manifest_smoke
        COMMAND mmcfilters_api_benchmark
            --manifest ${CMAKE_CURRENT_SOURCE_DIR}/workloads/scientific-api.ini
            --workload smoke_structured
            --repetitions 1
            --suite attributes,casf
            --emit-samples
    )
    set_tests_properties(benchmark_api_manifest_smoke PROPERTIES LABELS "benchmark;scientific")
endif()

foreach(target IN ITEMS
        mmcfilters_tree_validation_benchmark
        mmcfilters_adjacency_iteration_benchmark
        mmcfilters_min_max_residual_tree_benchmark
        mmcfilters_tree_construction_comparison_benchmark
        mmcfilters_scientific_pipeline_benchmark
        mmcfilters_validation_sensitive_algorithms_benchmark
        mmcfilters_api_benchmark)
    target_link_libraries(${target} PRIVATE mmcfilters::core)
    target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
    target_compile_features(${target} PRIVATE cxx_std_20)
endforeach()

target_include_directories(mmcfilters_min_max_residual_tree_benchmark
    SYSTEM PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../external/stb)

target_include_directories(mmcfilters_tree_construction_comparison_benchmark
    SYSTEM PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../external/stb)

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/tree_validation_benchmark.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/adjacency_iteration_benchmark.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/min_max_residual_tree_benchmark.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tree_construction_comparison_benchmark.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/scientific_pipeline_benchmark.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/validation_sensitive_algorithms_benchmark.cpp
    ${MMCFILTERS_API_BENCHMARK_SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/api/ApiBenchmark.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/support/ScientificBenchmark.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/support/stb_image_read.cpp
)
