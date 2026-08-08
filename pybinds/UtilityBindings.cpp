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
        .value("EuclideanDisk", RegularGridAdjacencyShape::EuclideanDisk)
        .value("StructuringElement", RegularGridAdjacencyShape::StructuringElement);

    py::class_<RegularGridAdjacency2D>(m, "RegularGridAdjacency2D", py::module_local(false),
                                       "Reusable regular-grid 2D adjacency with row-major index iteration.")
        .def(py::init<int, int, double>(), "rows"_a, "cols"_a, "radius"_a, "Create an adjacency relation for a regular 2D grid domain.")
        .def_static(
            "fromStructuringElement",
            [](int rows, int cols, const std::vector<std::pair<int, int>>& offsets) {
                std::vector<GridOffset2D> converted;
                converted.reserve(offsets.size());
                for (const auto [row, col] : offsets) {
                    converted.push_back({row, col});
                }
                return RegularGridAdjacency2D::fromStructuringElement(rows, cols, converted);
            },
            "rows"_a, "cols"_a, "offsets"_a, "Create undirected adjacency from symmetric (row, col) offsets.")
        .def_static("rectangular", &RegularGridAdjacency2D::rectangular, "rows"_a, "cols"_a, "rowRadius"_a, "colRadius"_a,
                    "Create a centered rectangular adjacency.")
        .def_static("line", &RegularGridAdjacency2D::line, "rows"_a, "cols"_a, "rowExtent"_a, "colExtent"_a, "Create a centered digital-line adjacency.")
        .def_static("horizontalLine", &RegularGridAdjacency2D::horizontalLine, "rows"_a, "cols"_a, "halfLength"_a,
                    "Create a centered horizontal-line adjacency.")
        .def_static("verticalLine", &RegularGridAdjacency2D::verticalLine, "rows"_a, "cols"_a, "halfLength"_a, "Create a centered vertical-line adjacency.")
        .def_property_readonly("size", &RegularGridAdjacency2D::getSize, "Number of offsets in the adjacency stencil.")
        .def_property_readonly("radius", &RegularGridAdjacency2D::getRadius, "Configured or bounding Euclidean stencil radius.")
        .def_property_readonly("shape", &RegularGridAdjacency2D::getShape, "Stencil construction family.")
        .def_property_readonly(
            "offsets",
            [](const RegularGridAdjacency2D& self) {
                std::vector<std::pair<int, int>> offsets;
                offsets.reserve(static_cast<std::size_t>(self.getSize()));
                for (int index = 0; index < self.getSize(); ++index) {
                    offsets.emplace_back(self.getOffsetRow(index), self.getOffsetCol(index));
                }
                return offsets;
            },
            "Canonical clockwise (row, col) stencil offsets.")
        .def(
            "adjacentIndices",
            [](const RegularGridAdjacency2D& self, int row, int col) {
                const auto range = self.getAdjacentIndices(row, col);
                return std::vector<int>(range.begin(), range.end());
            },
            "row"_a, "col"_a, "Return adjacent row-major grid indices, including the origin.")
        .def(
            "neighborIndices",
            [](const RegularGridAdjacency2D& self, int row, int col) {
                const auto range = self.getNeighborIndices(row, col);
                return std::vector<int>(range.begin(), range.end());
            },
            "row"_a, "col"_a, "Return neighboring row-major grid indices.")
        .def(
            "forwardNeighborIndices",
            [](const RegularGridAdjacency2D& self, int row, int col) {
                const auto range = self.getForwardNeighborIndices(row, col);
                return std::vector<int>(range.begin(), range.end());
            },
            "row"_a, "col"_a, "Return one directed half of the neighboring grid indices.");
}

} // namespace mmcfilters::pybindings
