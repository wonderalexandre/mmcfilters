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

struct ContourTraceComputationBinding {};

std::vector<ContourTraceEdge> collectEdges(ContourTraceComputation::IncrementalContourTraces::EdgeRange range) {
    std::vector<ContourTraceEdge> values;
    values.reserve(range.size());
    for (ContourTraceEdge edge : range) {
        values.push_back(edge);
    }
    return values;
}

} // namespace

/**
 * @brief Registers contour trace computation bindings in the Python module.
 *
 * @param m Python module receiving the bindings.
 */
void initContourTraceComputation(py::module_& m) {
    using Traces = ContourTraceComputation::IncrementalContourTraces;

    py::enum_<ContourTraceSide>(m, "ContourTraceSide", py::module_local(false))
        .value("NORTH", ContourTraceSide::North)
        .value("WEST", ContourTraceSide::West)
        .value("EAST", ContourTraceSide::East)
        .value("SOUTH", ContourTraceSide::South)
        .export_values();

    py::enum_<ContourLoopKind>(m, "ContourLoopKind", py::module_local(false))
        .value("EXTERNAL", ContourLoopKind::External)
        .value("INTERNAL", ContourLoopKind::Internal)
        .export_values();

    py::class_<ContourTraceEdge>(m, "ContourTraceEdge", py::module_local(false), "One side-level contour edge attached to a support pixel.")
        .def(py::init<>())
        .def_readwrite("pixel", &ContourTraceEdge::pixel)
        .def_readwrite("side", &ContourTraceEdge::side)
        .def("__repr__", [](const ContourTraceEdge& edge) { return "<ContourTraceEdge pixel=" + std::to_string(edge.pixel) + ">"; });

    py::class_<ContourTraceLoop>(m, "ContourTraceLoop", py::module_local(false), "Metadata for one ordered contour loop.")
        .def(py::init<>())
        .def_readwrite("kind", &ContourTraceLoop::kind)
        .def_readwrite("edge_offset", &ContourTraceLoop::edgeOffset)
        .def_readwrite("edge_count", &ContourTraceLoop::edgeCount)
        .def_readwrite("signed_area2", &ContourTraceLoop::signedArea2)
        .def("__repr__", [](const ContourTraceLoop& loop) { return "<ContourTraceLoop edge_count=" + std::to_string(loop.edgeCount) + ">"; });

    py::class_<Traces, std::shared_ptr<Traces>>(m, "ContourTraces", py::module_local(false),
                                                R"doc(Incremental contour trace result with lazy edge and loop materialization.

Edges are exposed as support-pixel sides. Loops are traced from oriented edges,
with external and internal boundaries separated by orientation. The source tree
is retained for the lifetime of this result.)doc")
        .def(
            "get_edges", [](const Traces& self, NodeId nodeId) { return collectEdges(self.getEdges(nodeId)); }, "node_id"_a,
            "Return unordered side-level contour edges for one live internal node.")
        .def("get_loops", &Traces::getLoops, "node_id"_a, "Return an independent list of contour loop metadata for one live internal node.")
        .def(
            "get_loop_edges", [](const Traces& self, const ContourTraceLoop& loop) { return collectEdges(self.getLoopEdges(loop)); }, "loop"_a,
            "Return ordered side-level contour edges for one loop.")
        .def("materialize_all", &Traces::materializeAll, "Materialize and trace every live-node contour.")
        .def_property_readonly("is_materialized", &Traces::isMaterialized, "Whether every live-node trace has been materialized.")
        .def("is_edge_materialized", &Traces::isEdgeMaterialized, "node_id"_a, "Return whether one live node's boundary edges have been materialized.")
        .def("is_node_traced", &Traces::isNodeTraced, "node_id"_a, "Return whether one live node's loops have been traced.");

    py::class_<ContourTraceComputationBinding>(m, "ContourTraceComputation", py::module_local(false),
                                               "Factory for definitive geometric contour trace extraction on morphological trees.")
        .def_static(
            "extraction",
            [](const std::shared_ptr<PythonValuedMorphologicalTree>& valuedTree) {
                if (!valuedTree) {
                    throw std::invalid_argument("Contour trace extraction requires a non-null ValuedMorphologicalTree.");
                }
                return valuedTree->visit([](const auto& concreteTree) { return ContourTraceComputation::extract(concreteTree->asView()); });
            },
            py::keep_alive<0, 1>(), "tree"_a.none(false), "Extract lazy contour traces from a valued tree.");
}

} // namespace mmcfilters::pybindings
