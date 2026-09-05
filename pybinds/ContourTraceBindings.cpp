#include "ModuleBindings.hpp"
#include "PythonValuedMorphologicalTree.hpp"

#include "../mmcfilters/contours/ContourTraceComputation.hpp"
#include "../mmcfilters/trees/ValuedMorphologicalTree.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>
#include <stdexcept>
#include <vector>

namespace mmcfilters::pybindings {

namespace py = pybind11;
using namespace pybind11::literals;

namespace {

std::vector<ContourEdge> copyEdges(ContourEdgeRange range) {
    return std::vector<ContourEdge>(range.begin(), range.end());
}

std::vector<PixelId> copyPixels(ContourPixelRange range) {
    return std::vector<PixelId>(range.begin(), range.end());
}

std::vector<ContourBoundary> copyBoundaries(std::span<const ContourBoundary> boundaries) {
    return std::vector<ContourBoundary>(boundaries.begin(), boundaries.end());
}

struct ContourTraceIterator {
    ContourTraceComputation::iterator iterator;
    bool hasStarted = false;

    explicit ContourTraceIterator(const ContourTraceComputation& traces) : iterator(traces.begin()) {}

    py::tuple next() {
        if (hasStarted && iterator != std::default_sentinel) {
            ++iterator;
        }
        hasStarted = true;
        if (iterator == std::default_sentinel) {
            throw py::stop_iteration();
        }
        const auto [node, traceView] = *iterator;
        return py::make_tuple(node, ContourTrace(traceView));
    }
};

} // namespace

/** @brief Registers sequential and node-local contour trace access. */
void initContourTraceComputation(py::module_& m) {
    py::enum_<ContourSide>(m, "ContourSide", py::module_local(false))
        .value("NORTH", ContourSide::North)
        .value("WEST", ContourSide::West)
        .value("EAST", ContourSide::East)
        .value("SOUTH", ContourSide::South)
        .export_values();

    py::enum_<ContourBoundaryKind>(m, "ContourBoundaryKind", py::module_local(false))
        .value("EXTERNAL", ContourBoundaryKind::External)
        .value("INTERNAL", ContourBoundaryKind::Internal)
        .export_values();

    py::class_<ContourEdge>(m, "ContourEdge", py::module_local(false), "One contour edge identified by its support pixel and side.")
        .def(py::init<>())
        .def_readwrite("pixel", &ContourEdge::pixel)
        .def_readwrite("side", &ContourEdge::side)
        .def("__repr__", [](const ContourEdge& edge) { return "<ContourEdge pixel=" + std::to_string(edge.pixel) + ">"; });

    py::class_<ContourBoundary>(m, "ContourBoundary", py::module_local(false), "Metadata for one ordered contour boundary.")
        .def(py::init<>())
        .def_readwrite("kind", &ContourBoundary::kind)
        .def_readwrite("edge_offset", &ContourBoundary::edgeOffset)
        .def_readwrite("edge_count", &ContourBoundary::edgeCount)
        .def_readwrite("doubled_signed_area", &ContourBoundary::doubledSignedArea)
        .def("__repr__", [](const ContourBoundary& boundary) { return "<ContourBoundary edge_count=" + std::to_string(boundary.edgeCount) + ">"; });

    py::class_<ContourTrace>(m, "ContourTrace", py::module_local(false), "Owned ordered contour trace for one tree node.")
        .def("boundaries", [](const ContourTrace& self) { return copyBoundaries(self.boundaries()); },
             "Return independent boundary descriptors.")
        .def("external_boundary", &ContourTrace::externalBoundary, "Return the unique external boundary.")
        .def("edges", [](const ContourTrace& self) { return copyEdges(self.edges()); }, "Return every contour edge in boundary order.")
        .def("boundary_edges", [](const ContourTrace& self, const ContourBoundary& boundary) { return copyEdges(self.boundaryEdges(boundary)); },
             "boundary"_a, "Return the ordered edges of one boundary.")
        .def("boundary_pixels", [](const ContourTrace& self, const ContourBoundary& boundary) { return copyPixels(self.boundaryPixels(boundary)); },
             "boundary"_a, "Return one support pixel for each ordered boundary edge, including repetitions.");

    py::class_<ContourTraceIterator>(m, "_ContourTraceIterator", "Single-pass iterator returning independently owned node traces.")
        .def("__iter__", [](ContourTraceIterator& self) -> ContourTraceIterator& { return self; }, py::return_value_policy::reference_internal)
        .def("__next__", &ContourTraceIterator::next);

    using Computation = ContourTraceComputation;
    py::class_<Computation, std::shared_ptr<Computation>>(m, "ContourTraceComputation", py::module_local(false),
               R"doc(Incremental ordered contour traces for a stable valued tree.

Iterate (node_id, trace) in post-order, use for_each_trace for callback access,
or query one node with trace. Python trace results own their data.)doc")
        .def(py::init([](const std::shared_ptr<PythonValuedMorphologicalTree>& tree) {
                 if (!tree) {
                     throw std::invalid_argument("ContourTraceComputation requires a non-null valued tree.");
                 }
                 return tree->visit([](const auto& concreteTree) { return ContourTraceComputation(concreteTree->asView()); });
             }), py::keep_alive<1, 2>(), "tree"_a.none(false))
        .def("__iter__", [](const Computation& self) { return ContourTraceIterator(self); }, py::keep_alive<0, 1>())
        .def("trace", &Computation::trace, "node_id"_a, "Return an owned ordered trace for one live internal node.")
        .def("for_each_trace", [](const Computation& self, const py::function& consumer) {
            self.forEachTrace([&](NodeId node, ContourTraceView traceView) { consumer(node, ContourTrace(traceView)); });
        }, "consumer"_a, "Call consumer(node_id, trace) in post-order; each trace owns its data.");
}

} // namespace mmcfilters::pybindings
