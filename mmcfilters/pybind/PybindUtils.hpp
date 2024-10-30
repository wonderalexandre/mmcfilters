
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

#ifndef PYBIND_UTILS_H
#define PYBIND_UTILS_H

class PybindUtils{
    public:
        static py::array_t<int> toNumpy(int *data, int size){
            return py::array(py::buffer_info(
                data,            /* Pointer to data (nullptr -> ask NumPy to allocate!) */
                sizeof(int),     /* Size of one item */
                py::format_descriptor<int>::value, /* Buffer format */
                1,          /* How many dimensions? */
                { size },  /* Number of elements for each dimension */
                { sizeof(int) }  /* Strides for each dimension */
            ));
        }

        static py::array_t<float> toNumpyFloat(float *data, int size){
            return py::array(py::buffer_info(
                data,            /* Pointer to data (nullptr -> ask NumPy to allocate!) */
                sizeof(float),     /* Size of one item */
                py::format_descriptor<float>::value, /* Buffer format */
                1,          /* How many dimensions? */
                { size },  /* Number of elements for each dimension */
                { sizeof(float) }  /* Strides for each dimension */
            ));
        }





};

#endif
