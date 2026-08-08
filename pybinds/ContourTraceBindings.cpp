#include "ModuleBindings.hpp"

#include "../mmcfilters/contours/ContourTraceComputation.hpp"
#include "../mmcfilters/trees/WeightedMorphologicalTree.hpp"

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
        .value("North", ContourTraceSide::North)
        .value("West", ContourTraceSide::West)
        .value("East", ContourTraceSide::East)
        .value("South", ContourTraceSide::South)
        .export_values();

    py::enum_<ContourLoopKind>(m, "ContourLoopKind", py::module_local(false))
        .value("External", ContourLoopKind::External)
        .value("Internal", ContourLoopKind::Internal)
        .export_values();

    py::class_<ContourTraceEdge>(m, "ContourTraceEdge", py::module_local(false), "One side-level contour edge attached to a support pixel.")
        .def(py::init<>())
        .def_readwrite("pixel", &ContourTraceEdge::pixel)
        .def_readwrite("side", &ContourTraceEdge::side)
        .def("__repr__", [](const ContourTraceEdge& edge) { return "<ContourTraceEdge pixel=" + std::to_string(edge.pixel) + ">"; });

    py::class_<ContourTraceLoop>(m, "ContourTraceLoop", py::module_local(false), "Metadata for one ordered contour loop.")
        .def(py::init<>())
        .def_readwrite("kind", &ContourTraceLoop::kind)
        .def_readwrite("edgeOffset", &ContourTraceLoop::edgeOffset)
        .def_readwrite("edgeCount", &ContourTraceLoop::edgeCount)
        .def_readwrite("signedArea2", &ContourTraceLoop::signedArea2)
        .def_property_readonly("edge_offset", [](const ContourTraceLoop& loop) { return loop.edgeOffset; })
        .def_property_readonly("edge_count", [](const ContourTraceLoop& loop) { return loop.edgeCount; })
        .def_property_readonly("signed_area2", [](const ContourTraceLoop& loop) { return loop.signedArea2; })
        .def("__repr__", [](const ContourTraceLoop& loop) { return "<ContourTraceLoop edgeCount=" + std::to_string(loop.edgeCount) + ">"; });

    py::class_<Traces, std::shared_ptr<Traces>>(m, "ContourTraces", py::module_local(false),
                                                R"doc(Incremental contour trace result with lazy edge and loop materialization.

Edges are exposed as support-pixel sides. Loops are traced from oriented edges,
with external and internal boundaries separated by orientation. The source tree
is retained for the lifetime of this result.)doc")
        .def(
            "getEdges", [](const Traces& self, NodeId nodeId) { return collectEdges(self.getEdges(nodeId)); }, "nodeId"_a,
            "Return unordered side-level contour edges for one live internal node.")
        .def("getLoops", &Traces::getLoops, "nodeId"_a, "Return an independent list of contour loop metadata for one live internal node.")
        .def(
            "getLoopEdges", [](const Traces& self, const ContourTraceLoop& loop) { return collectEdges(self.getLoopEdges(loop)); }, "loop"_a,
            "Return ordered side-level contour edges for one loop.")
        .def("materializeAll", &Traces::materializeAll, "Materialize and trace every live-node contour.")
        .def_property_readonly("isMaterialized", &Traces::isMaterialized, "Whether every live-node trace has been materialized.")
        .def("isEdgeMaterialized", &Traces::isEdgeMaterialized, "nodeId"_a, "Return whether one live node's boundary edges have been materialized.")
        .def("isNodeTraced", &Traces::isNodeTraced, "nodeId"_a, "Return whether one live node's loops have been traced.");

    py::class_<ContourTraceComputationBinding>(m, "ContourTraceComputation", py::module_local(false),
                                               "Factory for definitive geometric contour trace extraction on morphological trees.")
        .def_static(
            "extraction",
            [](const std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>>& weighted) {
                if (!weighted) {
                    throw std::invalid_argument("Contour trace extraction requires a non-null WeightedMorphologicalTree.");
                }
                return ContourTraceComputation::extract(weighted->asView());
            },
            py::keep_alive<0, 1>(), "tree"_a.none(false), "Extract lazy contour traces from a weighted tree.");
}

} // namespace mmcfilters::pybindings
