# Instruções de compilação
#
# Compilar somente os arquivos de teste: 
#    cmake -DBUILD_PYBIND=OFF -DBUILD_DEBUG=ON -DBUILD_TESTS=ON ../
#    cmake -DBUILD_PYBIND=OFF -DBUILD_DEBUG=OFF -DBUILD_TESTS=ON  ../
#
# Compilar o modulo Python:
#    cmake -DPYTHON_LIBRARY_DIR=$(python -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))") -DPYTHON_EXECUTABLE=$(python -c "import sys; print(sys.executable)") ../
#
# Depois:
#   make
#   make install

cmake_minimum_required(VERSION 3.14)
project(mmcfilters)

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

# Padrões de C++
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Incluir pybind11 se BUILD_PYBIND estiver ativado
if(BUILD_PYBIND)
    
    include(pybind11.cmake)
    file (GLOB SOURCE_FILES "mmcfilters/src/*.cpp")
    file (GLOB HEADER_FILES "mmcfilters/include/*.hpp" "mmcfilters/pybind/*.hpp")

    # Set up such that XCode organizes the files
    source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCE_FILES} ${HEADER_FILES} )

    # Arquivo fonte principal do módulo pybind
    set(PYBIND_CPP_FILE "mmcfilters/mmcfilters.cpp")
    
    pybind11_add_module(mmcfilters ${PYBIND_CPP_FILE}
      ${SOURCE_FILES}
      ${HEADER_FILES}
    )

    install(TARGETS mmcfilters
      COMPONENT python
      LIBRARY DESTINATION "${PYTHON_LIBRARY_DIR}"
      )

endif()

# Condicionalmente adicionar a compilação de testes
if(BUILD_TESTS)
    #crei uma lib mmcfilters_lib
    add_subdirectory(mmcfilters)

    # Adicionar o executável de testes: testeResidualTree.cpp
    add_executable(testeResidualTree 
        mmcfilters/test/testeResidualTree.cpp
    )
    target_link_libraries(testeResidualTree PRIVATE mmcfilters_lib)


    # Adicionar o executável de testes: 
    add_executable(testeTreeOfShapes 
        mmcfilters/test/testeTreeOfShapes.cpp
    )
    target_link_libraries(testeTreeOfShapes PRIVATE mmcfilters_lib)

    # Adicionar o executável de testes: 
    add_executable(TesteCountors 
        mmcfilters/test/TesteCountors.cpp
        external/stb/stb_image.cpp
    )
    target_link_libraries(TesteCountors PRIVATE mmcfilters_lib)
    target_compile_options(TesteCountors PRIVATE -Wno-deprecated-declarations)

        # Adicionar o executável de testes: 
    add_executable(TesteCountors2 
        mmcfilters/test/TesteCountors2.cpp
        external/stb/stb_image.cpp
    )
    target_link_libraries(TesteCountors2 PRIVATE mmcfilters_lib)
    target_compile_options(TesteCountors2 PRIVATE -Wno-deprecated-declarations)
        
    add_executable(TesteAdjacencyUC 
        mmcfilters/test/TesteAdjacencyUC.cpp
        external/stb/stb_image.cpp
    )
    target_link_libraries(TesteAdjacencyUC PRIVATE mmcfilters_lib)
    target_compile_options(TesteAdjacencyUC PRIVATE -Wno-deprecated-declarations)

    add_executable(testeComputerAttribute 
        mmcfilters/test/testeComputerAttribute.cpp
        external/stb/stb_image.cpp
    )
    target_link_libraries(testeComputerAttribute PRIVATE mmcfilters_lib)
    target_compile_options(testeComputerAttribute PRIVATE -Wno-deprecated-declarations)



endif()
