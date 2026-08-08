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

include(GNUInstallDirs)

set(MMCFILTERS_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
get_filename_component(MMCFILTERS_INCLUDE_ROOT "${MMCFILTERS_ROOT}" DIRECTORY)
include("${MMCFILTERS_INCLUDE_ROOT}/cmake/mmcfiltersPublicHeaders.cmake")

# Headers are intentionally explicit: adding a file does not make it part of
# the installed API by accident.
set(MMCFILTERS_HEADER_FILES ${MMCFILTERS_INSTALL_HEADERS})

if(NOT TARGET mmcfilters_lib)
    add_library(mmcfilters_lib INTERFACE)
endif()
set_target_properties(mmcfilters_lib PROPERTIES EXPORT_NAME core)

# Exponha os diretórios de include para consumidores externos.
target_include_directories(mmcfilters_lib
    INTERFACE
        $<BUILD_INTERFACE:${MMCFILTERS_INCLUDE_ROOT}>
        $<BUILD_INTERFACE:${MMCFILTERS_ROOT}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/mmcfilters>
)

# Garante C++20 para tudo que linkar nesta biblioteca.
target_compile_features(mmcfilters_lib INTERFACE cxx_std_20)

# Alias canônico para consumo moderno via target namespace.
add_library(mmcfilters::core ALIAS mmcfilters_lib)

# Organiza os headers por árvore no Xcode/Visual Studio.
source_group(TREE ${MMCFILTERS_ROOT} FILES ${MMCFILTERS_HEADER_FILES})

# Regras de instalação padrão: lib + headers.
install(TARGETS mmcfilters_lib
        EXPORT mmcfiltersTargets
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

foreach(header IN LISTS MMCFILTERS_INSTALL_HEADERS)
    file(RELATIVE_PATH
        header_relative
        "${MMCFILTERS_INCLUDE_ROOT}"
        "${header}")
    get_filename_component(header_relative_dir "${header_relative}" DIRECTORY)
    install(
        FILES "${header}"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${header_relative_dir}")
endforeach()

install(EXPORT mmcfiltersTargets
        NAMESPACE mmcfilters::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/mmcfilters)
