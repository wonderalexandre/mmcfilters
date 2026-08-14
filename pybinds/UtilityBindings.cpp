#include "ModuleBindings.hpp"

#include "../mmcfilters/utils/RegularGridAdjacency2D.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <span>
#include <utility>
#include <vector>

namespace mmcfilters::pybindings {

namespace py = pybind11;
using namespace pybind11::literals;

/**
 * @brief Registers regular grid adjacency2 d bindings in the Python module.
 *
 * @param m Python module receiving the bindings.
 */
void initRegularGridAdjacency2D(py::module_& m) {
    py::enum_<RegularGridAdjacencyShape>(m, "RegularGridAdjacencyShape", py::module_local(false))
        .value("EUCLIDEAN_DISK", RegularGridAdjacencyShape::EuclideanDisk)
        .value("STRUCTURING_ELEMENT", RegularGridAdjacencyShape::StructuringElement);

    py::class_<RegularGridAdjacency2D>(m, "RegularGridAdjacency2D", py::module_local(false),
                                       "Reusable regular-grid 2D adjacency with row-major index iteration.")
        .def(py::init<int, int, double>(), "rows"_a, "columns"_a, "radius"_a, "Create an adjacency relation for a regular 2D grid domain.")
        .def_static(
            "from_structuring_element",
            [](int rows, int columns, const std::vector<std::pair<int, int>>& offsets) {
                std::vector<GridOffset2D> converted;
                converted.reserve(offsets.size());
                for (const auto [row, column] : offsets) {
                    converted.push_back({row, column});
                }
                return RegularGridAdjacency2D::fromStructuringElement(rows, columns, converted);
            },
            "rows"_a, "columns"_a, "offsets"_a, "Create undirected adjacency from symmetric (row, column) offsets.")
        .def_static("rectangular", &RegularGridAdjacency2D::rectangular, "rows"_a, "columns"_a, "row_radius"_a, "column_radius"_a,
                    "Create a centered rectangular adjacency.")
        .def_static("line", &RegularGridAdjacency2D::line, "rows"_a, "columns"_a, "row_extent"_a, "column_extent"_a, "Create a centered digital-line adjacency.")
        .def_static("horizontal_line", &RegularGridAdjacency2D::horizontalLine, "rows"_a, "columns"_a, "half_length"_a,
                    "Create a centered horizontal-line adjacency.")
        .def_static("vertical_line", &RegularGridAdjacency2D::verticalLine, "rows"_a, "columns"_a, "half_length"_a, "Create a centered vertical-line adjacency.")
        .def_property_readonly("size", &RegularGridAdjacency2D::getSize, "Number of offsets in the adjacency stencil.")
        .def_property_readonly("radius", &RegularGridAdjacency2D::getRadius, "Configured or bounding Euclidean stencil radius.")
        .def_property_readonly("shape", &RegularGridAdjacency2D::getShape, "Stencil construction family.")
        .def_property_readonly(
            "offsets",
            [](const RegularGridAdjacency2D& self) {
                std::vector<std::pair<int, int>> offsets;
                offsets.reserve(static_cast<std::size_t>(self.getSize()));
                for (int index = 0; index < self.getSize(); ++index) {
                    offsets.emplace_back(self.getOffsetRow(index), self.getOffsetColumn(index));
                }
                return offsets;
            },
            "Canonical clockwise (row, column) stencil offsets.")
        .def(
            "adjacent_indices",
            [](const RegularGridAdjacency2D& self, int row, int column) {
                const auto range = self.getAdjacentIndices(row, column);
                return std::vector<int>(range.begin(), range.end());
            },
            "row"_a, "column"_a, "Return adjacent row-major grid indices, including the origin.")
        .def(
            "neighbor_indices",
            [](const RegularGridAdjacency2D& self, int row, int column) {
                const auto range = self.getNeighborIndices(row, column);
                return std::vector<int>(range.begin(), range.end());
            },
            "row"_a, "column"_a, "Return neighboring row-major grid indices.")
        .def(
            "forward_neighbor_indices",
            [](const RegularGridAdjacency2D& self, int row, int column) {
                const auto range = self.getForwardNeighborIndices(row, column);
                return std::vector<int>(range.begin(), range.end());
            },
            "row"_a, "column"_a, "Return one directed half of the neighboring grid indices.");
}

} // namespace mmcfilters::pybindings
