cmake_minimum_required(VERSION 3.14)
project(mmcfilters_tests LANGUAGES CXX)


# C++20 sem extensões específicas de compilador
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Define build padrão, se não for especificado
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

# Flags específicas por tipo de build
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -Wextra -Wpedantic")

# Mensagens úteis
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ Standard: C++${CMAKE_CXX_STANDARD}")

set(MMCFILTERS_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../mmcfilters")
add_subdirectory(${MMCFILTERS_SOURCE_DIR} ${CMAKE_BINARY_DIR}/mmcfilters_build)

function(mmcfilters_add_test target)
    add_executable(${target} ${ARGN})
    target_link_libraries(${target} PRIVATE mmcfilters_lib)
    target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
endfunction()

mmcfilters_add_test(testeResidualTree testeResidualTree.cpp)
mmcfilters_add_test(TestMorphologicalTree TestMorphologicalTree.cpp)
mmcfilters_add_test(testeTreeOfShapes testeTreeOfShapes.cpp)

mmcfilters_add_test(TesteContours TesteContours.cpp ../external/stb/stb_image.cpp)
target_compile_options(TesteContours PRIVATE -Wno-deprecated-declarations)

mmcfilters_add_test(TesteContours2 TesteContours2.cpp ../external/stb/stb_image.cpp)
target_compile_options(TesteContours2 PRIVATE -Wno-deprecated-declarations)

mmcfilters_add_test(TesteContours3 TesteContours3.cpp ../external/stb/stb_image.cpp)
target_compile_options(TesteContours3 PRIVATE -Wno-deprecated-declarations)

mmcfilters_add_test(TesteAdjacencyUC TesteAdjacencyUC.cpp ../external/stb/stb_image.cpp)
target_compile_options(TesteAdjacencyUC PRIVATE -Wno-deprecated-declarations)

mmcfilters_add_test(testeComputerAttribute testeComputerAttribute.cpp ../external/stb/stb_image.cpp)
target_compile_options(testeComputerAttribute PRIVATE -Wno-deprecated-declarations)

mmcfilters_add_test(TesteBitQuads TesteBitQuads.cpp ../external/stb/stb_image.cpp)
target_compile_options(TesteBitQuads PRIVATE -Wno-deprecated-declarations)
