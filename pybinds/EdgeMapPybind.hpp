#pragma once

#include "PybindConversions.hpp"
#include "../mmcfilters/trees/saliency/HierarchySaliencyMapProjection.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <utility>
#include <variant>
#include <vector>

namespace mmcfilters::pybindings {

namespace py = pybind11;

/** @brief Owns an immutable contour edge set returned to Python. */
class EdgeContoursPybind {
    EdgeContourMap contours_; ///< Contour endpoints and image-domain metadata.

  public:
    /**
     * @brief Takes ownership of contour buffers without copying them.
     * @param contours Contour edge set produced by the native algorithm.
     */
    explicit EdgeContoursPybind(EdgeContourMap&& contours) : contours_(std::move(contours)) {}

    /**
     * @brief Returns the owned contour data for read-only queries.
     * @return Immutable contour endpoints and domain metadata.
     */
    const EdgeContourMap& data() const { return contours_; }
};

/** @brief Owns an immutable edge map while retaining its native value type. */
class EdgeMapPybind {
    /// Edge buffers and metadata in the native value type.
    std::variant<EdgeSaliencyMap<int>, EdgeSaliencyMap<float>, EdgeSaliencyMap<double>> map_;

  public:
    /**
     * @brief Takes ownership of an edge map without copying its buffers.
     * @param map Edge saliency map produced by the native algorithm.
     */
    template <class Value> explicit EdgeMapPybind(EdgeSaliencyMap<Value>&& map) : map_(std::move(map)) {}

    /**
     * @brief Reads the map using its native value type.
     * @param visitor Function applied to the immutable native map.
     * @return Result of the visitor.
     */
    template <class Visitor> decltype(auto) visit(Visitor&& visitor) const {
        return std::visit(std::forward<Visitor>(visitor), map_);
    }

    /**
     * @brief Projects values to pixels directly from the owned buffers.
     * @param reducer Reduction applied to incident edge values at each pixel.
     * @return Independently owned pixel image used for display.
     */
    py::array toPixelImage(EdgeToPixelReducer reducer) const {
        return visit([&](const auto& map) -> py::array {
            return pybind_utils::toNumpy(HierarchySaliencyMapProjection::edgeMapToPixelImage(map, reducer));
        });
    }

    /**
     * @brief Selects contour edges whose values reach the threshold.
     * @param threshold Minimum edge value retained in the contour set.
     * @return Owned contour set containing the selected endpoints.
     */
    EdgeContoursPybind threshold(double threshold) const {
        return visit([&](const auto& map) {
            return EdgeContoursPybind(HierarchySaliencyMapProjection::thresholdCut(map, threshold));
        });
    }
};

/**
 * @brief Returns a read-only array that retains the Python buffer owner.
 * @param values Native buffer exposed without a copy.
 * @param owner Python object that owns the buffer.
 * @return Read-only NumPy view retaining the owner.
 */
template <class Value> py::array readOnlyResultArray(const std::vector<Value>& values, py::handle owner) {
    // A non-null pointer prevents NumPy from allocating storage for empty arrays.
    static const Value emptyValue{};
    py::array_t<Value> array({static_cast<py::ssize_t>(values.size())},
                            {static_cast<py::ssize_t>(sizeof(Value))},
                            values.empty() ? &emptyValue : values.data(), owner);
    array.attr("setflags")(false);
    return array;
}

/**
 * @brief Registers immutable edge-map results and their operations.
 * @param m Python module receiving the result types.
 */
inline void initEdgeMapBindings(py::module_& m) {
    using namespace pybind11::literals;

    py::class_<EdgeContoursPybind>(m, "EdgeContours", py::module_local(false),
                                   "Owned contour edge set with read-only NumPy endpoint views.")
        .def_property_readonly("num_rows", [](const EdgeContoursPybind& self) { return self.data().numRows; })
        .def_property_readonly("num_columns", [](const EdgeContoursPybind& self) { return self.data().numColumns; })
        .def_property_readonly("adjacency_radius", [](const EdgeContoursPybind& self) { return self.data().adjacencyRadius; })
        .def_property_readonly("sources", [](py::object self) {
            return readOnlyResultArray(self.cast<const EdgeContoursPybind&>().data().sources, self);
        })
        .def_property_readonly("targets", [](py::object self) {
            return readOnlyResultArray(self.cast<const EdgeContoursPybind&>().data().targets, self);
        })
        .def("__len__", [](const EdgeContoursPybind& self) { return self.data().size(); });

    py::class_<EdgeMapPybind>(m, "EdgeMap", py::module_local(false),
                              "Owned edge saliency map with read-only NumPy views independent of the source tree.")
        .def_property_readonly("num_rows", [](const EdgeMapPybind& self) {
            return self.visit([](const auto& map) { return map.numRows; });
        })
        .def_property_readonly("num_columns", [](const EdgeMapPybind& self) {
            return self.visit([](const auto& map) { return map.numColumns; });
        })
        .def_property_readonly("adjacency_radius", [](const EdgeMapPybind& self) {
            return self.visit([](const auto& map) { return map.adjacencyRadius; });
        })
        .def_property_readonly("sources", [](py::object self) {
            return self.cast<const EdgeMapPybind&>().visit([&](const auto& map) {
                return readOnlyResultArray(map.sources, self);
            });
        })
        .def_property_readonly("targets", [](py::object self) {
            return self.cast<const EdgeMapPybind&>().visit([&](const auto& map) {
                return readOnlyResultArray(map.targets, self);
            });
        })
        .def_property_readonly("values", [](py::object self) {
            return self.cast<const EdgeMapPybind&>().visit([&](const auto& map) {
                return readOnlyResultArray(map.values, self);
            });
        })
        .def("__len__", [](const EdgeMapPybind& self) {
            return self.visit([](const auto& map) { return map.size(); });
        })
        .def("to_pixel_image", &EdgeMapPybind::toPixelImage, "reducer"_a = EdgeToPixelReducer::Max,
             "Project incident edge values to pixels using MAX or MEAN; this is a display projection.")
        .def("threshold", &EdgeMapPybind::threshold, "threshold"_a,
             "Return EdgeContours selecting edges whose values are greater than or equal to threshold.");
}

} // namespace mmcfilters::pybindings
