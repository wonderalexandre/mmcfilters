#include "ModuleBindings.hpp"

#include "../mmcfilters/contours/ContoursComputedIncrementally.hpp"
#include "../mmcfilters/trees/WeightedMorphologicalTree.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace mmcfilters::pybindings {

namespace py = pybind11;
using namespace pybind11::literals;

namespace {

struct ContourComputationBinding {};

} // namespace

/**
 * @brief Registers contours computed incrementally bindings in the Python module.
 *
 * @param m Python module receiving the bindings.
 */
void initContoursComputedIncrementally(py::module_& m) {
    using Contours = ContoursComputedIncrementally::IncrementalContours;
    using ContourRange = Contours::ContourRange;
    using Range = decltype(std::declval<Contours&>().contoursByNode());
    using Iter = decltype(std::declval<Range&>().begin());

    py::class_<ContourRange>(m, "ContourRange", py::module_local(false), "Lazy iterable range of row-major pixel indices for one internal-node contour.")
        .def(
            "__iter__", [](const ContourRange& p) { return py::make_iterator(p.begin(), p.end()); }, py::keep_alive<0, 1>(),
            "Iterate contour pixels as row-major integer indices.")
        .def("empty", &ContourRange::empty, "Return true when the contour has no pixels.");

    struct ContoursIterator {
        Range range;
        Iter it;
        Iter itEnd;

        explicit ContoursIterator(Contours& self) : range(self.contoursByNode()), it(range.begin()), itEnd(range.end()) {}
    };

    py::class_<ContoursIterator>(m, "ContoursIterator", py::module_local(false), "Iterator over `(node_id, ContourRange)` pairs for all live nodes.")
        .def(py::init<Contours&>(), "contours"_a)
        .def(
            "__iter__", [](ContoursIterator& self) -> ContoursIterator& { return self; }, py::return_value_policy::reference_internal, "Return this iterator.")
        .def(
            "__next__",
            [](ContoursIterator& self) -> py::object {
                if (self.it == self.itEnd) {
                    throw py::stop_iteration();
                }
                auto entry = *self.it++;
                auto nodeId = std::get<0>(entry);
                auto proxy = std::get<1>(entry);
                return py::make_tuple(nodeId, proxy);
            },
            "Return the next `(node_id, ContourRange)` pair.");

    py::class_<Contours, std::shared_ptr<Contours>>(m, "Contours", py::module_local(false),
                                                    R"doc(Incremental contour result with lazy materialization.

Contours are exposed as row-major pixel-index ranges. The first read of a node
may materialize and cache its subtree contour; later reads use cached storage.)doc")
        .def(
            "contoursByNode", [](Contours& self) { return ContoursIterator(self); }, py::keep_alive<0, 1>(),
            "Iterate `(node_id, contour_range)` for every live internal node.")
        .def("getContour", &Contours::getContour, "nodeId"_a, "Return a lazy contour range for one live internal node.")
        .def("materializeAll", &Contours::materializeAll, "Materialize and cache contours for every live node.")
        .def_property_readonly("isMaterialized", &Contours::isMaterialized, "Whether all live-node contours have been materialized.")
        .def("isContourMaterialized", &Contours::isContourMaterialized, "nodeId"_a, "Return whether one live node contour has been materialized.");

    py::class_<ContourComputationBinding>(m, "ContourComputation", py::module_local(false),
                                          "Factory for incremental contour extraction on morphological trees.")
        .def_static(
            "extraction",
            [](const std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>>& weighted) {
                if (!weighted) {
                    throw std::invalid_argument("Contour extraction requires a non-null WeightedMorphologicalTree.");
                }
                return ContoursComputedIncrementally::extractCompactContours(weighted->asView());
            },
            "tree"_a.none(false), "Extract compact contours from a weighted tree.");
}

} // namespace mmcfilters::pybindings
