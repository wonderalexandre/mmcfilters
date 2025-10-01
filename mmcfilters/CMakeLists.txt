# -----------------------------------------------------------------------------
# mmcfilters/CMakeLists.txt
# -----------------------------------------------------------------------------
# Como usar:
#   1. Adicione este diretório no seu projeto com `add_subdirectory(mmcfilters)`.
#   2. Consuma a biblioteca através do alvo `mmcfilters::core`:
#         target_link_libraries(meu_alvo PRIVATE mmcfilters::core)
#   3. Inclua os headers normalmente (ex.: `#include "mmcfilters/trees/..."`).
#   4. Instale/exporte usando os comandos abaixo ou complemente conforme
#      necessário no seu projeto principal.
#
# Este CMakeLists gera apenas o núcleo C++ independente de python/pybind. A
# criação de módulos Python deve ocorrer em nível superior reutilizando este
# alvo para evitar recompilações e dependências desnecessárias.

set(MMCFILTERS_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
set(MMCFILTERS_COMPONENT_DIRS
    ${MMCFILTERS_ROOT}/trees
    ${MMCFILTERS_ROOT}/filters
    ${MMCFILTERS_ROOT}/attributes
    ${MMCFILTERS_ROOT}/contours
    ${MMCFILTERS_ROOT}/utils
)

# Agrupa headers e fontes do núcleo. Usamos GLOB para manter listas atualizadas
# automaticamente ao adicionar novos arquivos.
file(GLOB_RECURSE MMCFILTERS_HEADER_FILES
    "${MMCFILTERS_ROOT}/trees/*.hpp"
    "${MMCFILTERS_ROOT}/trees/*.tpp"
    "${MMCFILTERS_ROOT}/filters/*.hpp"
    "${MMCFILTERS_ROOT}/attributes/*.hpp"
    "${MMCFILTERS_ROOT}/contours/*.hpp"
    "${MMCFILTERS_ROOT}/utils/*.hpp"
)

file(GLOB_RECURSE MMCFILTERS_SOURCE_FILES
    "${MMCFILTERS_ROOT}/trees/*.cpp"
    "${MMCFILTERS_ROOT}/filters/*.cpp"
    "${MMCFILTERS_ROOT}/attributes/*.cpp"
    "${MMCFILTERS_ROOT}/contours/*.cpp"
    "${MMCFILTERS_ROOT}/utils/*.cpp"
)

if(NOT TARGET mmcfilters_lib)
    add_library(mmcfilters_lib)
    target_sources(mmcfilters_lib
        PRIVATE
            ${MMCFILTERS_HEADER_FILES}
            ${MMCFILTERS_SOURCE_FILES}
    )
endif()

# Exponha os diretórios de include para consumidores externos.
target_include_directories(mmcfilters_lib
    PUBLIC
        $<BUILD_INTERFACE:${MMCFILTERS_ROOT}>
        $<INSTALL_INTERFACE:include/mmcfilters>
)

# Garante C++20 para tudo que linkar nesta biblioteca.
target_compile_features(mmcfilters_lib PUBLIC cxx_std_20)

# Ajusta flags de otimização dependendo do compilador: evita passar -O3 ao MSVC.
target_compile_options(mmcfilters_lib
    PUBLIC
        $<$<CXX_COMPILER_ID:MSVC>:/O2>
        $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-O3>
)

# Alias canônico para consumo moderno via target namespace.
add_library(mmcfilters::core ALIAS mmcfilters_lib)

# Organiza os arquivos por árvore no Xcode/Visual Studio.
source_group(TREE ${MMCFILTERS_ROOT} FILES
    ${MMCFILTERS_HEADER_FILES}
    ${MMCFILTERS_SOURCE_FILES}
)

# Regras de instalação padrão: lib + headers.
install(TARGETS mmcfilters_lib
        EXPORT mmcfiltersTargets
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

foreach(dir IN LISTS MMCFILTERS_COMPONENT_DIRS)
    get_filename_component(dir_name ${dir} NAME)
    install(DIRECTORY ${dir}/
            DESTINATION include/mmcfilters/${dir_name}
            FILES_MATCHING PATTERN "*.hpp" PATTERN "*.tpp")
endforeach()

install(EXPORT mmcfiltersTargets
        NAMESPACE mmcfilters::
        DESTINATION lib/cmake/mmcfilters)
