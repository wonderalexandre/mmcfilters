# Definir flags de compilação condicionalmente (Release ou Debug)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(MMCFILTERS_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
set(MMCFILTERS_COMPONENT_DIRS
    ${MMCFILTERS_ROOT}/trees
    ${MMCFILTERS_ROOT}/filters
    ${MMCFILTERS_ROOT}/attributes
    ${MMCFILTERS_ROOT}/contours
    ${MMCFILTERS_ROOT}/utils
)

file(GLOB_RECURSE HEADER_FILES
    "${MMCFILTERS_ROOT}/trees/*.hpp"
    "${MMCFILTERS_ROOT}/trees/*.tpp"
    "${MMCFILTERS_ROOT}/filters/*.hpp"
    "${MMCFILTERS_ROOT}/attributes/*.hpp"
    "${MMCFILTERS_ROOT}/contours/*.hpp"
    "${MMCFILTERS_ROOT}/utils/*.hpp"
)

file(GLOB_RECURSE SOURCE_FILES
    "${MMCFILTERS_ROOT}/trees/*.cpp"
    "${MMCFILTERS_ROOT}/filters/*.cpp"
    "${MMCFILTERS_ROOT}/attributes/*.cpp"
    "${MMCFILTERS_ROOT}/contours/*.cpp"
    "${MMCFILTERS_ROOT}/utils/*.cpp"
)

add_library(mmcfilters_lib SHARED ${HEADER_FILES} ${SOURCE_FILES})

foreach(dir IN LISTS MMCFILTERS_COMPONENT_DIRS)
    target_include_directories(mmcfilters_lib PRIVATE ${dir})
endforeach()

if(BUILD_PYBIND)
    if(DEFINED MMC_PYBINDS_DIR)
        get_filename_component(PROJECT_PYBIND_DIR "${CMAKE_SOURCE_DIR}/${MMC_PYBINDS_DIR}" ABSOLUTE)
    else()
        get_filename_component(PROJECT_PYBIND_DIR "${MMCFILTERS_ROOT}/../pybinds" ABSOLUTE)
    endif()
    file(GLOB_RECURSE PYBIND_HEADER_FILES "${PROJECT_PYBIND_DIR}/*.hpp")
    file(GLOB_RECURSE PYBIND_SOURCE_FILES "${PROJECT_PYBIND_DIR}/*.cpp")

    target_sources(mmcfilters_lib PRIVATE ${PYBIND_HEADER_FILES} ${PYBIND_SOURCE_FILES})
    target_include_directories(mmcfilters_lib PRIVATE ${PROJECT_PYBIND_DIR})
endif()

source_group(TREE ${MMCFILTERS_ROOT} FILES ${HEADER_FILES} ${SOURCE_FILES})

if(BUILD_PYBIND)
    source_group(TREE ${PROJECT_PYBIND_DIR} PREFIX "Pybind" FILES ${PYBIND_HEADER_FILES} ${PYBIND_SOURCE_FILES})
endif()

install(TARGETS mmcfilters_lib DESTINATION lib)

foreach(dir IN LISTS MMCFILTERS_COMPONENT_DIRS)
    get_filename_component(dir_name ${dir} NAME)
    install(DIRECTORY ${dir}/ DESTINATION include/mmcfilters/${dir_name}
            FILES_MATCHING PATTERN "*.hpp" PATTERN "*.tpp")
endforeach()

if(BUILD_PYBIND)
    install(DIRECTORY ${PROJECT_PYBIND_DIR}/ DESTINATION include/mmcfilters/pybinds
            FILES_MATCHING PATTERN "*.hpp")
endif()
