#include "ModuleBindings.hpp"
#include "PythonValuedMorphologicalTree.hpp"

#include "../mmcfilters/contours/ContourTraceComputation.hpp"
#include "../mmcfilters/trees/ValuedMorphologicalTree.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <vector>

namespace mmcfilters::pybindings {

namespace py = pybind11;
using namespace pybind11::literals;

namespace {

std::vector<ContourEdge> copyEdges(ContourTraceComputation::EdgeRange range) {
    std::vector<ContourEdge> edges;
    edges.reserve(range.size());
    for (ContourEdge edge : range) {
        edges.push_back(edge);
    }
    return edges;
}

std::vector<PixelId> copyPixels(ContourTraceComputation::PixelRange range) {
    std::vector<PixelId> pixels;
    pixels.reserve(range.size());
    for (PixelId pixel : range) {
        pixels.push_back(pixel);
    }
    return pixels;
}

} // namespace

/**
 * @brief Registers contour trace computation bindings in the Python module.
 *
 * @param m Python module receiving the bindings.
 */
void initContourTraceComputation(py::module_& m) {
    using Computation = ContourTraceComputation;

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

    py::class_<Computation, std::shared_ptr<Computation>>(m, "ContourTraceComputation", py::module_local(false),
               R"doc(Lazy contour-edge construction and ordered boundary tracing.

Edges are exposed as support-pixel sides. Boundaries are traced from oriented edges,
with external and internal boundaries distinguished by orientation. The computation
keeps its source tree alive.)doc")
        .def(py::init([](const std::shared_ptr<PythonValuedMorphologicalTree>& tree) {
                 if (!tree) {
                     throw std::invalid_argument("Contour trace computation requires a non-null ValuedMorphologicalTree.");
                 }
                 return tree->visit([](const auto& concreteTree) { return ContourTraceComputation(concreteTree->asView()); });
             }), py::keep_alive<1, 2>(), "tree"_a.none(false))
        .def(
            "edges", [](const Computation& self, NodeId nodeId) { return copyEdges(self.edges(nodeId)); }, "node_id"_a,
            "Return the unordered contour edges of one live internal node.")
        .def("boundaries", &Computation::boundaries, "node_id"_a, "Return an independent list of contour boundary metadata for one live internal node.")
        .def("external_boundary", &Computation::externalBoundary, "node_id"_a, "Return the unique external boundary of one live internal node.")
        .def(
            "boundary_edges", [](const Computation& self, const ContourBoundary& boundary) { return copyEdges(self.boundaryEdges(boundary)); }, "boundary"_a,
            "Return the ordered contour edges of one boundary.")
        .def(
            "boundary_pixels", [](const Computation& self, const ContourBoundary& boundary) { return copyPixels(self.boundaryPixels(boundary)); }, "boundary"_a,
            "Return one support pixel per ordered boundary edge, including repetitions.")
        .def("trace_all", &Computation::traceAll, "Trace the contours of every live node.")
        .def_property_readonly("has_traced_all_boundaries", &Computation::hasTracedAllBoundaries,
                               "True when the ordered boundaries of every live node are cached.")
        .def("has_cached_edges", &Computation::hasCachedEdges, "node_id"_a, "Return whether the contour edges of one live node are cached.")
        .def("has_traced_boundaries", &Computation::hasTracedBoundaries, "node_id"_a,
             "Return whether the ordered boundaries of one live node are cached.");
}

} // namespace mmcfilters::pybindings
