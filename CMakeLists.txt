# -----------------------------------------------------------------------------
# CMakeLists.txt (raiz)
# -----------------------------------------------------------------------------
# Como usar:
#   1. Configure normalmente: `cmake -S . -B build`.
#   2. Para compilar apenas o núcleo C++ use `-DMMCFILTERS_BUILD_PYTHON=OFF`.
#   3. Para habilitar os testes use `-DMMCFILTERS_BUILD_TESTS=ON` e então
#      `cmake --build build --target testeQualquer` ou `ctest`.
#   4. O CMake produz um pacote Python importável na árvore de build para testes
#      locais. A instalação Python oficial usa `pip`/`scikit-build-core`; uma
#      instalação CMake publica apenas o pacote C++.

cmake_minimum_required(VERSION 3.20)
project(mmcfilters LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(MMCFILTERS_BUILD_PYTHON "Build the pybind11 extension module" ON)
option(MMCFILTERS_BUILD_TESTS "Build C++ test executables" OFF)
option(MMCFILTERS_BUILD_EXAMPLES "Build example executables" OFF)
option(MMCFILTERS_BUILD_BENCHMARKS "Build internal benchmark executables" OFF)
option(MMCFILTERS_ENABLE_ASSERTS "Keep runtime asserts enabled even in Release" OFF)
option(MMCFILTERS_ENABLE_DOC_TARGETS "Enable Doxygen documentation targets" ON)

# Núcleo C++ compartilhado entre as saídas (biblioteca + bindings Python).
add_subdirectory(mmcfilters)

if(MMCFILTERS_ENABLE_ASSERTS)
    target_compile_definitions(mmcfilters_lib INTERFACE MMCFILTERS_ENABLE_ASSERTS)
    target_compile_options(mmcfilters_lib
        INTERFACE
            $<$<CXX_COMPILER_ID:MSVC>:/UNDEBUG>
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-UNDEBUG>
    )
endif()

if(MMCFILTERS_BUILD_PYTHON)
    include(pybind11.cmake)

    if(NOT DEFINED MMC_PYBINDS_DIR)
        set(MMC_PYBINDS_DIR "pybinds")
    endif()

    if(NOT IS_ABSOLUTE "${MMC_PYBINDS_DIR}")
        set(MMC_PYBINDS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${MMC_PYBINDS_DIR}")
    endif()

    set(MMCFILTERS_PYBIND_MAIN "${MMC_PYBINDS_DIR}/mmcfilters.cpp")
    set(MMCFILTERS_PYTHON_PACKAGE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/python/mmcfilters")
    set(MMCFILTERS_PYTHON_BUILD_PACKAGE_DIR "${CMAKE_BINARY_DIR}/python/mmcfilters")

    if(NOT EXISTS "${MMCFILTERS_PYBIND_MAIN}")
        message(FATAL_ERROR
            "Arquivo de binding não encontrado: ${MMCFILTERS_PYBIND_MAIN}.\n"
            "Verifique se o diretório 'pybinds' está presente.")
    endif()

    file(GLOB MMCFILTERS_PYBIND_CPP CONFIGURE_DEPENDS "${MMC_PYBINDS_DIR}/*.cpp")
    file(GLOB_RECURSE MMCFILTERS_PYBIND_HEADERS "${MMC_PYBINDS_DIR}/*.hpp")
    file(GLOB MMCFILTERS_PYTHON_PACKAGE_FILES CONFIGURE_DEPENDS
        "${MMCFILTERS_PYTHON_PACKAGE_DIR}/*.py")

    file(MAKE_DIRECTORY "${MMCFILTERS_PYTHON_BUILD_PACKAGE_DIR}")
    foreach(_package_file IN LISTS MMCFILTERS_PYTHON_PACKAGE_FILES)
        get_filename_component(_package_name "${_package_file}" NAME)
        configure_file(
            "${_package_file}"
            "${MMCFILTERS_PYTHON_BUILD_PACKAGE_DIR}/${_package_name}"
            COPYONLY
        )
    endforeach()

    pybind11_add_module(mmcfilters MODULE
        ${MMCFILTERS_PYBIND_CPP}
        ${MMCFILTERS_PYBIND_HEADERS}
    )

    target_link_libraries(mmcfilters PRIVATE mmcfilters::core)
    target_include_directories(mmcfilters PRIVATE ${MMC_PYBINDS_DIR})
    target_compile_features(mmcfilters PRIVATE cxx_std_20)
    set_target_properties(mmcfilters PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY "${MMCFILTERS_PYTHON_BUILD_PACKAGE_DIR}"
        RUNTIME_OUTPUT_DIRECTORY "${MMCFILTERS_PYTHON_BUILD_PACKAGE_DIR}"
        ARCHIVE_OUTPUT_DIRECTORY "${MMCFILTERS_PYTHON_BUILD_PACKAGE_DIR}")

    source_group(TREE ${MMC_PYBINDS_DIR} PREFIX "Pybind" FILES
        ${MMCFILTERS_PYBIND_CPP}
        ${MMCFILTERS_PYBIND_HEADERS}
    )

    # scikit-build-core stages the native extension beside the Python sources
    # selected by wheel.packages. A plain CMake install intentionally remains a
    # complete C++ installation instead of publishing an incomplete Python
    # package containing only the extension module.
    if(DEFINED SKBUILD_PLATLIB_DIR AND NOT SKBUILD_PLATLIB_DIR STREQUAL "")
        install(TARGETS mmcfilters
            LIBRARY DESTINATION "${SKBUILD_PLATLIB_DIR}/mmcfilters"
            RUNTIME DESTINATION "${SKBUILD_PLATLIB_DIR}/mmcfilters"
            ARCHIVE DESTINATION "${SKBUILD_PLATLIB_DIR}/mmcfilters")
    endif()
endif()

if(MMCFILTERS_BUILD_TESTS)
    enable_testing()
    add_subdirectory(unit-tests)
endif()

if(MMCFILTERS_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(MMCFILTERS_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

if(MMCFILTERS_ENABLE_DOC_TARGETS)
    find_package(Doxygen QUIET)
    if(DOXYGEN_FOUND)
        # Doxygen canonicalizes Markdown link targets. Canonicalize INPUT paths as
        # well so documentation builds remain warning-free when the checkout is
        # reached through a symlink (for example /Users/... -> /Volumes/...).
        file(REAL_PATH "${CMAKE_CURRENT_SOURCE_DIR}" MMCFILTERS_DOCS_SOURCE_DIR)
        set(MMCFILTERS_DOCS_PUBLIC_DIR "${CMAKE_CURRENT_BINARY_DIR}/docs/public")
        set(MMCFILTERS_DOCS_INTERNAL_DIR "${CMAKE_CURRENT_BINARY_DIR}/docs/internal")
        set(MMCFILTERS_DOCS_PUBLIC_WARN_LOG "${CMAKE_CURRENT_BINARY_DIR}/docs/doxygen-public.log")
        set(MMCFILTERS_DOCS_INTERNAL_WARN_LOG "${CMAKE_CURRENT_BINARY_DIR}/docs/doxygen-internal.log")
        file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/docs")

        configure_file(
            "${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.public.in"
            "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile.public"
            @ONLY
        )
        configure_file(
            "${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.internal.in"
            "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile.internal"
            @ONLY
        )

        add_custom_target(docs-public
            COMMAND "${DOXYGEN_EXECUTABLE}" "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile.public"
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
            COMMENT "Generating public API documentation"
            VERBATIM
        )
        add_custom_target(docs-internal
            COMMAND "${DOXYGEN_EXECUTABLE}" "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile.internal"
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
            COMMENT "Generating internal contributor documentation"
            VERBATIM
        )
        add_custom_target(docs
            DEPENDS docs-public docs-internal
        )
    else()
        message(STATUS "Doxygen not found; docs-public and docs-internal targets are unavailable.")
    endif()
endif()

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/mmcfiltersConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/mmcfiltersConfig.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/mmcfilters")

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/mmcfiltersConfig.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/mmcfilters")
