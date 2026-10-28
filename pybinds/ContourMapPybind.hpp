#pragma once

#include "EdgeMapPybind.hpp"

namespace mmcfilters::pybindings {

/** @brief Owns transition edges and their boundary-owning hierarchy nodes. */
class NodeContourEdgesPybind {
    NodeContourEdgeMap contours_; ///< Native edge endpoints, nodes, and domain metadata.

  public:
    /**
     * @brief Takes ownership of projected contour edges.
     * @param contours Native contour projection transferred without copying.
     */
    explicit NodeContourEdgesPybind(NodeContourEdgeMap&& contours) : contours_(std::move(contours)) {}

    /**
     * @brief Reads the immutable contour projection.
     * @return Edge endpoints, owning nodes, and domain metadata.
     */
    const NodeContourEdgeMap& data() const { return contours_; }
};

/**
 * @brief Dispatches a dense node valuation without copying contiguous inputs.
 * @param values Input NumPy valuation retained during the native operation.
 * @param numNodeSlots Required size of the dense node-slot domain.
 * @param visitor Operation applied to a typed read-only valuation span.
 * @return Result of the operation.
 */
template <class Visitor> auto visitNodeValuation(py::array values, int numNodeSlots, Visitor&& visitor) {
    py::array contiguous = py::array::ensure(values, py::array::c_style);
    if (!contiguous || (!contiguous.dtype().is(py::dtype::of<float>()) && !contiguous.dtype().is(py::dtype::of<double>()))) {
        throw std::invalid_argument("node_valuation must be a 1D C-contiguous float32 or float64 array");
    }
    pybind_utils::require1DArray(contiguous.request(), numNodeSlots, "node_valuation");
    if (contiguous.dtype().is(py::dtype::of<float>())) {
        return visitor(std::span<const float>(static_cast<const float*>(contiguous.data()), static_cast<std::size_t>(numNodeSlots)));
    }
    return visitor(std::span<const double>(static_cast<const double*>(contiguous.data()), static_cast<std::size_t>(numNodeSlots)));
}

/** @brief Owns immutable per-node transition contour slices. */
class IncrementalNodeContoursPybind {
    IncrementalNodeContourMap contours_; ///< Native contour slices and image-domain metadata.

  public:
    /**
     * @brief Takes ownership of incremental contour buffers.
     * @param contours Native contour data transferred without copying.
     */
    explicit IncrementalNodeContoursPybind(IncrementalNodeContourMap&& contours) : contours_(std::move(contours)) {}

    /**
     * @brief Reads the immutable contour slices.
     * @return Node offsets, edge endpoints, and domain metadata.
     */
    const IncrementalNodeContourMap& data() const { return contours_; }

    /**
     * @brief Projects a dense node valuation onto the stored transition edges.
     * @param values Float32 or float64 valuation with one value per node slot.
     * @return Owned sparse edge map retaining the input value dtype.
     */
    EdgeMapPybind projectNodeValuation(py::array values) const {
        return visitNodeValuation(std::move(values), contours_.numNodeSlots, [&](auto valuation) {
            return EdgeMapPybind(HierarchySaliencyMapProjection::projectNodeValuation(contours_, valuation));
        });
    }

    /**
     * @brief Selects contour slices whose node valuations reach a threshold.
     * @param values Float32 or float64 valuation with one value per node slot.
     * @param threshold Minimum node valuation retained in the contour set.
     * @return Owned contour edge set selected with greater-than-or-equal comparison.
     */
    EdgeContoursPybind thresholdByNodeValuation(py::array values, double threshold) const {
        return visitNodeValuation(std::move(values), contours_.numNodeSlots, [&](auto valuation) {
            return EdgeContoursPybind(HierarchySaliencyMapProjection::thresholdByNodeValuation(contours_, valuation, threshold));
        });
    }
};

/**
 * @brief Registers immutable projected and incremental contour results.
 * @param m Python module receiving the result types.
 */
inline void initContourMapBindings(py::module_& m) {
    using namespace pybind11::literals;
    py::class_<NodeContourEdgesPybind>(m, "NodeContourEdges", py::module_local(false),
                                     "Owned transition edge endpoints and their boundary-owning hierarchy nodes.")
        .def_property_readonly("num_rows", [](const NodeContourEdgesPybind& self) { return self.data().numRows; })
        .def_property_readonly("num_columns", [](const NodeContourEdgesPybind& self) { return self.data().numColumns; })
        .def_property_readonly("adjacency_radius", [](const NodeContourEdgesPybind& self) { return self.data().adjacencyRadius; })
        .def_property_readonly("sources", [](py::object self) {
            return readOnlyResultArray(self.cast<const NodeContourEdgesPybind&>().data().sources, self);
        })
        .def_property_readonly("targets", [](py::object self) {
            return readOnlyResultArray(self.cast<const NodeContourEdgesPybind&>().data().targets, self);
        })
        .def_property_readonly("nodes", [](py::object self) {
            return readOnlyResultArray(self.cast<const NodeContourEdgesPybind&>().data().nodes, self);
        })
        .def("__len__", [](const NodeContourEdgesPybind& self) { return self.data().size(); });

    py::class_<IncrementalNodeContoursPybind>(m, "IncrementalNodeContours", py::module_local(false),
                                            "Owned per-node transition edge slices with read-only NumPy views.")
        .def_property_readonly("num_rows", [](const IncrementalNodeContoursPybind& self) { return self.data().numRows; })
        .def_property_readonly("num_columns", [](const IncrementalNodeContoursPybind& self) { return self.data().numColumns; })
        .def_property_readonly("num_node_slots", [](const IncrementalNodeContoursPybind& self) { return self.data().numNodeSlots; })
        .def_property_readonly("adjacency_radius", [](const IncrementalNodeContoursPybind& self) { return self.data().adjacencyRadius; })
        .def_property_readonly("offsets", [](py::object self) {
            return readOnlyResultArray(self.cast<const IncrementalNodeContoursPybind&>().data().offsets, self);
        })
        .def_property_readonly("sources", [](py::object self) {
            return readOnlyResultArray(self.cast<const IncrementalNodeContoursPybind&>().data().sources, self);
        })
        .def_property_readonly("targets", [](py::object self) {
            return readOnlyResultArray(self.cast<const IncrementalNodeContoursPybind&>().data().targets, self);
        })
        .def("__len__", [](const IncrementalNodeContoursPybind& self) { return self.data().size(); })
        .def("project_node_valuation", &IncrementalNodeContoursPybind::projectNodeValuation, "node_valuation"_a,
             "Return an EdgeMap over transition edges; omitted finest-region edges have implicit value zero.")
        .def("threshold_by_node_valuation", &IncrementalNodeContoursPybind::thresholdByNodeValuation,
             "node_valuation"_a, "threshold"_a,
             "Return EdgeContours selecting node slices with valuation >= threshold; validate formal hierarchy valuations separately.");
}

} // namespace mmcfilters::pybindings
