cmake_minimum_required(VERSION 3.14) 
include(FetchContent)

set(FETCHCONTENT_QUIET FALSE)
FetchContent_Declare(
    pybind11
    URL https://github.com/pybind/pybind11/archive/refs/tags/v2.11.1.tar.gz
)

FetchContent_MakeAvailable(pybind11)