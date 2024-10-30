# Definir flags de compilação condicionalmente (Release ou Debug)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Diretórios do projeto
set(PROJECT_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(PROJECT_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")

# Usar GLOB para obter automaticamente todos os arquivos .hpp e .cpp
file(GLOB_RECURSE HEADER_FILES "${PROJECT_INCLUDE_DIR}/*.hpp")
file(GLOB_RECURSE SOURCE_FILES "${PROJECT_SOURCE_DIR}/*.cpp")

# Criação da biblioteca compartilhada com o nome mmcfilters_lib
add_library(mmcfilters_lib SHARED ${HEADER_FILES} ${SOURCE_FILES})

# Diretórios de inclusão
target_include_directories(mmcfilters_lib PRIVATE ${PROJECT_INCLUDE_DIR})

# Se pybind estiver ativado, incluir os arquivos de binding Python
if(BUILD_PYBIND)
    set(PROJECT_PYBIND_DIR "${CMAKE_CURRENT_SOURCE_DIR}/pybind")
    file(GLOB_RECURSE PYBIND_HEADER_FILES "${PROJECT_PYBIND_DIR}/*.hpp")
    file(GLOB_RECURSE PYBIND_SOURCE_FILES "${PROJECT_PYBIND_DIR}/*.cpp")

    # Adicionar os arquivos de pybind aos arquivos de compilação
    target_sources(mmcfilters_lib PRIVATE ${PYBIND_HEADER_FILES} ${PYBIND_SOURCE_FILES})

    # Incluir diretórios do pybind
    target_include_directories(mmcfilters_lib PRIVATE ${PROJECT_PYBIND_DIR})
endif()

# Agrupamento para IDEs
source_group(TREE ${PROJECT_INCLUDE_DIR} PREFIX "Headers" FILES ${HEADER_FILES})
source_group(TREE ${PROJECT_SOURCE_DIR} PREFIX "Sources" FILES ${SOURCE_FILES})

if(BUILD_PYBIND)
    source_group(TREE ${PROJECT_PYBIND_DIR} PREFIX "Pybind Headers" FILES ${PYBIND_HEADER_FILES})
    source_group(TREE ${PROJECT_PYBIND_DIR} PREFIX "Pybind Sources" FILES ${PYBIND_SOURCE_FILES})
endif()

# Instalação da biblioteca
install(TARGETS mmcfilters_lib DESTINATION lib)

# Instalação de headers
install(DIRECTORY ${PROJECT_INCLUDE_DIR}/ DESTINATION include/mmcfilters)

if(BUILD_PYBIND)
    install(DIRECTORY ${PROJECT_PYBIND_DIR}/ DESTINATION include/mmcfilters/pybind)
endif()
