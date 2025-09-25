# Instruções de compilação
#
# Compilar somente os arquivos de teste:
#    cmake -S tests -B tests/build
#    cmake --build tests/build
#
# Compilar o modulo Python:
#    cmake -DPYTHON_LIBRARY_DIR=$(python -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))") -DPYTHON_EXECUTABLE=$(python -c "import sys; print(sys.executable)") ../
#
# Depois:
#   make
#   make install

cmake_minimum_required(VERSION 3.20)

project(mmcfilters)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Opções configuráveis
option(BUILD_PYBIND "Enable building Python bindings with pybind11" ON)
option(BUILD_DEBUG "Enable Debug mode" OFF)
option(BUILD_TESTS "Enable building of tests" OFF)

# Configurações de compilação (Release por padrão)
if(BUILD_DEBUG)
    set(CMAKE_BUILD_TYPE Debug)
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
else()
    set(CMAKE_BUILD_TYPE Release)
    set(CMAKE_CXX_FLAGS "-O3")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()

# C++20 sem extensões específicas de compilador
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Mensagens úteis
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ Standard: C++${CMAKE_CXX_STANDARD}")

# Incluir pybind11 se BUILD_PYBIND estiver ativado
if(BUILD_PYBIND)
    
    include(pybind11.cmake)

    # Se MMC_PYBINDS_DIR não vier, use 'pybinds'
    if(NOT DEFINED MMC_PYBINDS_DIR)
    set(MMC_PYBINDS_DIR "pybinds")
    endif()

    # Se for relativo, torne-o absoluto baseado no CMakeLists atual
    if(NOT IS_ABSOLUTE "${MMC_PYBINDS_DIR}")
    set(MMC_PYBINDS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${MMC_PYBINDS_DIR}")
    endif()

    set(PYBIND_CPP_FILE "${MMC_PYBINDS_DIR}/mmcfilters.cpp")

    if(NOT EXISTS "${PYBIND_CPP_FILE}")
    message(FATAL_ERROR
        "pybind source não encontrado: ${PYBIND_CPP_FILE}\n"
        "Verifique se 'pybinds/mmcfilters.cpp' está no sdist (MANIFEST.in).")
    endif()

    set(MMCFILTERS_COMPONENT_DIRS
        mmcfilters/trees
        mmcfilters/filters
        mmcfilters/attributes
        mmcfilters/contours
        mmcfilters/utils
    )

    file(GLOB_RECURSE SOURCE_FILES
        "mmcfilters/trees/*.cpp"
        "mmcfilters/filters/*.cpp"
        "mmcfilters/attributes/*.cpp"
        "mmcfilters/contours/*.cpp"
        "mmcfilters/utils/*.cpp"
    )

    file(GLOB_RECURSE HEADER_FILES
        "mmcfilters/trees/*.hpp"
        "mmcfilters/trees/*.tpp"
        "mmcfilters/filters/*.hpp"
        "mmcfilters/attributes/*.hpp"
        "mmcfilters/contours/*.hpp"
        "mmcfilters/utils/*.hpp"
        "${MMC_PYBINDS_DIR}/*.hpp"
    )

    # Set up such that XCode organizes the files
    source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCE_FILES} ${HEADER_FILES} )

    # Arquivo fonte principal do módulo pybind
    set(PYBIND_CPP_FILE "${MMC_PYBINDS_DIR}/mmcfilters.cpp")
    
    pybind11_add_module(mmcfilters 
      ${PYBIND_CPP_FILE}
      ${SOURCE_FILES}
      ${HEADER_FILES}
    )

    # Install into Python's platlib when building with scikit-build-core.
    # Fallbacks: PYTHON_LIBRARY_DIR (legacy) -> standard CMake libdir.
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

# Condicionalmente adicionar a compilação de testes
# (alternativamente, use "cmake -S tests -B tests/build" para um build dedicado)
if(BUILD_TESTS)
    add_subdirectory(tests)
endif()
