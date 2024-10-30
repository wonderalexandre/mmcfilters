include(FetchContent)

set(FETCHCONTENT_QUIET FALSE)

FetchContent_Declare(
    pybind11
    URL https://github.com/pybind/pybind11/archive/v2.3.0.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE  # Adicione essa linha
)

FetchContent_MakeAvailable(pybind11)  # Substitui FetchContent_Populate e add_subdirectory

#FetchContent_GetProperties(pybind11)
#if(NOT pybind11_POPULATED)
#    FetchContent_Populate(pybind11)
#    add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
#endif()
