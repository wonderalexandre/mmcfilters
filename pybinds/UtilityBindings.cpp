#include "ModuleBindings.hpp"

#include "../mmcfilters/utils/AdjacencyRelation.hpp"

#include <pybind11/pybind11.h>

namespace mmcfilters::pybindings {

namespace py = pybind11;
using namespace pybind11::literals;

void initAdjacencyRelation(py::module_& m) {
    py::class_<AdjacencyRelation>(m, "AdjacencyRelation", py::module_local(false),
        "Reusable 2D image adjacency relation with row-major pixel iteration.")
        .def(py::init<int, int, double>(),
            "rows"_a,
            "cols"_a,
            "radius"_a,
            "Create an adjacency relation for a 2D image domain.")
        .def_property_readonly("size", &AdjacencyRelation::getSize,
            "Number of offsets in the adjacency stencil.")
        .def_property_readonly("radius", &AdjacencyRelation::getRadius,
            "Neighbourhood radius used to build the stencil.")
        .def("getAdjPixels", py::overload_cast<int, int>(&AdjacencyRelation::getAdjPixels),
            "row"_a,
            "col"_a,
            "Prepare iteration over adjacent row-major pixel indices.");
}

} // namespace mmcfilters::pybindings
