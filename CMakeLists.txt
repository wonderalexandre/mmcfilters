# -----------------------------------------------------------------------------
# CMakeLists.txt (raiz)
# -----------------------------------------------------------------------------
# Como usar:
#   1. Configure normalmente: `cmake -S . -B build`.
#   2. Para compilar apenas o núcleo C++ use `-DMMCFILTERS_BUILD_PYTHON=OFF`.
#   3. Para habilitar os testes use `-DMMCFILTERS_BUILD_TESTS=ON` e então
#      `cmake --build build --target testeQualquer` ou `ctest`.
#   4. O módulo Python reutiliza a biblioteca C++ (alvo `mmcfilters::core`) e
#      pode ser instalado via `cmake --install build` ou ferramentas como
#      `scikit-build-core`.

cmake_minimum_required(VERSION 3.20)
project(mmcfilters LANGUAGES CXX)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(MMCFILTERS_BUILD_PYTHON "Build the pybind11 extension module" ON)
option(MMCFILTERS_BUILD_TESTS "Build C++ test executables" OFF)
option(MMCFILTERS_ENABLE_ASSERTS "Keep runtime asserts enabled even in Release" OFF)

# Núcleo C++ compartilhado entre as saídas (biblioteca + bindings Python).
add_subdirectory(mmcfilters)

if(MMCFILTERS_ENABLE_ASSERTS)
    target_compile_definitions(mmcfilters_lib PUBLIC MMCFILTERS_ENABLE_ASSERTS)
    target_compile_options(mmcfilters_lib
        PUBLIC
            $<$<CXX_COMPILER_ID:MSVC>:/UNDEBUG>
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-UNDEBUG>
    )
else()
    target_compile_definitions(mmcfilters_lib PUBLIC NDEBUG)
endif()

if(MMCFILTERS_BUILD_PYTHON)
    include(pybind11.cmake)

    if(NOT DEFINED MMC_PYBINDS_DIR)
        set(MMC_PYBINDS_DIR "pybinds")
    endif()

    if(NOT IS_ABSOLUTE "${MMC_PYBINDS_DIR}")
        set(MMC_PYBINDS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${MMC_PYBINDS_DIR}")
    endif()

    set(MMCFILTERS_PYBIND_CPP "${MMC_PYBINDS_DIR}/mmcfilters.cpp")

    if(NOT EXISTS "${MMCFILTERS_PYBIND_CPP}")
        message(FATAL_ERROR
            "Arquivo de binding não encontrado: ${MMCFILTERS_PYBIND_CPP}.\n"
            "Verifique se o diretório 'pybinds' está presente.")
    endif()

    file(GLOB_RECURSE MMCFILTERS_PYBIND_HEADERS "${MMC_PYBINDS_DIR}/*.hpp")

    pybind11_add_module(mmcfilters MODULE
        ${MMCFILTERS_PYBIND_CPP}
        ${MMCFILTERS_PYBIND_HEADERS}
    )

    target_link_libraries(mmcfilters PRIVATE mmcfilters::core)
    target_include_directories(mmcfilters PRIVATE ${MMC_PYBINDS_DIR})
    target_compile_features(mmcfilters PRIVATE cxx_std_20)

    source_group(TREE ${MMC_PYBINDS_DIR} PREFIX "Pybind" FILES
        ${MMCFILTERS_PYBIND_CPP}
        ${MMCFILTERS_PYBIND_HEADERS}
    )

    set(_PY_INSTALL_DIR "${SKBUILD_PLATLIB_DIR}")
    if(NOT _PY_INSTALL_DIR)
        if(DEFINED PYTHON_LIBRARY_DIR)
            set(_PY_INSTALL_DIR "${PYTHON_LIBRARY_DIR}")
        else()
            include(GNUInstallDirs)
            set(_PY_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}")
        endif()
    endif()

    install(TARGETS mmcfilters
        LIBRARY DESTINATION "${_PY_INSTALL_DIR}"
        RUNTIME DESTINATION "${_PY_INSTALL_DIR}"
        ARCHIVE DESTINATION "${_PY_INSTALL_DIR}")
endif()

if(MMCFILTERS_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
