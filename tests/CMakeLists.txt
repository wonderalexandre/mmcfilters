# -----------------------------------------------------------------------------
# tests/CMakeLists.txt
# -----------------------------------------------------------------------------
# Como usar:
#   • Quando chamado a partir do projeto raiz, o alvo `mmcfilters::core` já terá
#     sido criado por `add_subdirectory(mmcfilters)` e os executáveis apenas o usam.
#   • Para rodar isoladamente execute `cmake -S tests -B tests/build` – este
#     script garante que o núcleo seja adicionado se ainda não existir.
#   • Novos executáveis podem ser adicionados via `mmcfilters_add_program(nome arquivo.cpp)`,
#     recebendo automaticamente C++20 e warnings padrão.

cmake_minimum_required(VERSION 3.14)
project(mmcfilters_tests LANGUAGES CXX)

set(MMCFILTERS_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../mmcfilters")

# Se o alvo principal ainda não existir (execução isolada), inclui o núcleo.
if(NOT TARGET mmcfilters::core)
    add_subdirectory(${MMCFILTERS_SOURCE_DIR} ${CMAKE_BINARY_DIR}/mmcfilters_build)
endif()

if(NOT TARGET mmcfilters::core)
    message(FATAL_ERROR "O alvo mmcfilters::core não foi encontrado. Verifique se o diretório 'mmcfilters' existe." )
endif()


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



function(mmcfilters_add_program target)
    add_executable(${target} ${ARGN})
    target_link_libraries(${target} PRIVATE mmcfilters::core)
    target_compile_features(${target} PRIVATE cxx_std_20)
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
    target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
endfunction()

mmcfilters_add_program(testeResidualTree testeResidualTree.cpp)
mmcfilters_add_program(TestMorphologicalTree TestMorphologicalTree.cpp)
mmcfilters_add_program(testeTreeOfShapes testeTreeOfShapes.cpp)

mmcfilters_add_program(TesteContours TesteContours.cpp ../external/stb/stb_image.cpp)
target_compile_options(TesteContours PRIVATE -Wno-deprecated-declarations)

mmcfilters_add_program(TesteContours2 TesteContours2.cpp ../external/stb/stb_image.cpp)
target_compile_options(TesteContours2 PRIVATE -Wno-deprecated-declarations)

mmcfilters_add_program(TesteContours3 TesteContours3.cpp ../external/stb/stb_image.cpp)
target_compile_options(TesteContours3 PRIVATE -Wno-deprecated-declarations)

mmcfilters_add_program(TesteAdjacencyUC TesteAdjacencyUC.cpp ../external/stb/stb_image.cpp)
target_compile_options(TesteAdjacencyUC PRIVATE -Wno-deprecated-declarations)

mmcfilters_add_program(testeComputerAttribute testeComputerAttribute.cpp ../external/stb/stb_image.cpp)
target_compile_options(testeComputerAttribute PRIVATE -Wno-deprecated-declarations)

mmcfilters_add_program(TesteBitQuads TesteBitQuads.cpp ../external/stb/stb_image.cpp)
target_compile_options(TesteBitQuads PRIVATE -Wno-deprecated-declarations)
