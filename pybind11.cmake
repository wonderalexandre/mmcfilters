set(PYBIND11_FINDPYTHON ON)
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

set(_MMCFILTERS_PYBIND11_FALLBACK_VERSION "3.0.1")
set(_MMCFILTERS_PYBIND11_VERSION_RANGE "${_MMCFILTERS_PYBIND11_FALLBACK_VERSION}...<4")

# A Python build frontend already installs pybind11 from pyproject.toml. Prefer
# that exact environment so wheel and CMake dependency resolution cannot drift.
execute_process(
    COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
    RESULT_VARIABLE _mmcfilters_pybind11_cmakedir_result
    OUTPUT_VARIABLE _mmcfilters_pybind11_cmakedir
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)
if(_mmcfilters_pybind11_cmakedir_result EQUAL 0 AND
   IS_DIRECTORY "${_mmcfilters_pybind11_cmakedir}")
    find_package(pybind11 ${_MMCFILTERS_PYBIND11_VERSION_RANGE} CONFIG QUIET
        HINTS "${_mmcfilters_pybind11_cmakedir}")
else()
    find_package(pybind11 ${_MMCFILTERS_PYBIND11_VERSION_RANGE} CONFIG QUIET)
endif()

if(NOT TARGET pybind11::module)
    if(NOT MMCFILTERS_FETCH_PYBIND11)
        message(FATAL_ERROR
            "A compatible pybind11 (${_MMCFILTERS_PYBIND11_VERSION_RANGE}) was not found. "
            "Install it in ${Python_EXECUTABLE}'s environment or enable "
            "MMCFILTERS_FETCH_PYBIND11.")
    endif()

    include(FetchContent)
    set(FETCHCONTENT_QUIET FALSE)
    FetchContent_Declare(
        pybind11
        URL https://github.com/pybind/pybind11/archive/refs/tags/v${_MMCFILTERS_PYBIND11_FALLBACK_VERSION}.tar.gz
        URL_HASH SHA256=741633da746b7c738bb71f1854f957b9da660bcd2dce68d71949037f0969d0ca
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    FetchContent_MakeAvailable(pybind11)
    message(STATUS "Using bundled pybind11 ${_MMCFILTERS_PYBIND11_FALLBACK_VERSION}")
else()
    message(STATUS "Using local pybind11 ${pybind11_VERSION}")
endif()

unset(_MMCFILTERS_PYBIND11_FALLBACK_VERSION)
unset(_MMCFILTERS_PYBIND11_VERSION_RANGE)
unset(_mmcfilters_pybind11_cmakedir)
unset(_mmcfilters_pybind11_cmakedir_result)
