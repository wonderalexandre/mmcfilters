#include "ModuleBindings.hpp"
#include "PythonValuedMorphologicalTree.hpp"

#include "../mmcfilters/contours/ContourComputation.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <algorithm>
#include <memory>
#include <stdexcept>

namespace mmcfilters::pybindings {

namespace py = pybind11;
using namespace pybind11::literals;

namespace {

py::array_t<PixelId> copyContourPixels(std::span<const PixelId> pixels) {
    py::array_t<PixelId> pixelsCopy(static_cast<py::ssize_t>(pixels.size()));
    std::copy(pixels.begin(), pixels.end(), pixelsCopy.mutable_data());
    return pixelsCopy;
}

struct ContourIterator {
    ContourComputation::iterator iterator;
    bool hasStarted = false;

    explicit ContourIterator(const ContourComputation& contours) : iterator(contours.begin()) {}

    py::tuple next() {
        if (hasStarted && iterator != std::default_sentinel) {
            ++iterator;
        }
        hasStarted = true;
        if (iterator == std::default_sentinel) {
            throw py::stop_iteration();
        }
        const auto [node, pixels] = *iterator;
        return py::make_tuple(node, copyContourPixels(pixels));
    }
};

} // namespace

/** @brief Registers sequential pixel-contour access and queries for one node. */
void initContourComputation(py::module_& m) {
    py::class_<ContourIterator>(m, "_ContourIterator", "Single-pass iterator with independently owned NumPy contour arrays.")
        .def("__iter__", [](ContourIterator& self) -> ContourIterator& { return self; }, py::return_value_policy::reference_internal)
        .def("__next__", &ContourIterator::next);

    py::class_<ContourComputation>(m, "ContourComputation", R"doc(Incremental foreground A4 pixel contours for a stable valued tree.

Iterate (node_id, pixels) in post-order, or query one node with contour.
Each result is an independent NumPy array. No contour history is cached.)doc")
        .def(py::init([](const std::shared_ptr<PythonValuedMorphologicalTree>& tree) {
                 if (!tree) {
                     throw std::invalid_argument("ContourComputation requires a non-null valued tree.");
                 }
                 return tree->visit([](const auto& concreteTree) { return ContourComputation(concreteTree->asView()); });
             }),
             py::keep_alive<1, 2>(), "tree"_a.none(false))
        .def("__iter__", [](const ContourComputation& self) { return ContourIterator(self); }, py::keep_alive<0, 1>())
        .def("contour", [](const ContourComputation& self, NodeId node) { return copyContourPixels(self.contour(node)); }, "node_id"_a,
             "Return an owned contour array for one live internal node by scanning its support, without caching.")
        .def("for_each_contour", [](const ContourComputation& self, const py::function& consumer) {
            self.forEachContour([&](NodeId node, std::span<const PixelId> pixels) { consumer(node, copyContourPixels(pixels)); });
        }, "consumer"_a, "Call consumer(node_id, pixels) in post-order; each array owns its data.");
}

} // namespace mmcfilters::pybindings
