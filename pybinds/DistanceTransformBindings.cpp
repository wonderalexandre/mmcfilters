#include "ModuleBindings.hpp"
#include "PythonValuedMorphologicalTree.hpp"

#include "../mmcfilters/distanceTransform/DistanceTransformComputation.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace mmcfilters::pybindings {

namespace py = pybind11;
using namespace pybind11::literals;

namespace {

struct DistanceTransformComputationBinding {};

DistanceTransformAccuracy parseDistanceTransformAccuracy(const py::handle& value) {
    if (py::isinstance<py::str>(value)) {
        std::string name = py::cast<std::string>(value);
        std::transform(name.begin(), name.end(), name.begin(), [](unsigned char character) { return static_cast<char>(std::tolower(character)); });
        if (name == "exact") {
            return DistanceTransformAccuracy::Exact;
        }
        if (name == "approximate") {
            return DistanceTransformAccuracy::Approximate;
        }
        throw std::invalid_argument("Distance-transform accuracy must be 'exact' or 'approximate'.");
    }

    try {
        return py::cast<DistanceTransformAccuracy>(value);
    } catch (const py::cast_error&) {
        throw std::invalid_argument("Distance-transform accuracy must be DistanceTransformAccuracy.EXACT, "
                                    "DistanceTransformAccuracy.APPROXIMATE, 'exact', or 'approximate'.");
    }
}

} // namespace

/** @brief Registers public distance-field and profile bindings. */
void initDistanceTransformComputation(py::module_& m) {
    py::enum_<DistanceTransformAccuracy>(m, "DistanceTransformAccuracy", py::module_local(false))
        .value("APPROXIMATE", DistanceTransformAccuracy::Approximate)
        .value("EXACT", DistanceTransformAccuracy::Exact);

    py::class_<NodeDistanceField>(m, "NodeDistanceField", py::module_local(false), "Self-contained compact distance field for one live internal node.")
        .def_property_readonly("node_id", &NodeDistanceField::node)
        .def_property_readonly("rows", &NodeDistanceField::rows)
        .def_property_readonly("columns", &NodeDistanceField::columns)
        .def_property_readonly("accuracy", &NodeDistanceField::accuracy)
        .def_property_readonly("support_pixels",
                               [](const NodeDistanceField& field) { return std::vector<PixelId>(field.supportPixels().begin(), field.supportPixels().end()); })
        .def_property_readonly(
            "squared_distances",
            [](const NodeDistanceField& field) { return std::vector<std::int64_t>(field.squaredDistances().begin(), field.squaredDistances().end()); })
        .def(
            "to_dense",
            [](const NodeDistanceField& field, std::int64_t outsideValue) {
                const std::vector<std::int64_t> dense = field.denseSquaredDistances(outsideValue);
                py::array_t<std::int64_t> output({field.rows(), field.columns()});
                std::copy(dense.begin(), dense.end(), output.mutable_data());
                return output;
            },
            "outside_value"_a = -1, "Return a dense squared-distance image using the supplied value outside the node support.");

    py::class_<DistanceTransformProfiles>(m, "DistanceTransformProfiles", py::module_local(false),
                                          "Compact sparse distance histograms for all live internal nodes.")
        .def_property_readonly("accuracy", &DistanceTransformProfiles::accuracy)
        .def_property_readonly("num_node_slots", &DistanceTransformProfiles::numNodeSlots)
        .def_property_readonly("num_bins", &DistanceTransformProfiles::numBins)
        .def(
            "histogram",
            [](const DistanceTransformProfiles& profiles, NodeId node) {
                py::list result;
                for (const DistanceTransformHistogramBin& bin : profiles.histogram(node)) {
                    result.append(py::make_tuple(bin.squaredDistance, bin.count));
                }
                return result;
            },
            "node_id"_a, "Return sorted `(squared_distance, count)` bins for one live node.")
        .def(
            "survival_profile",
            [](const DistanceTransformProfiles& profiles, NodeId node) {
                py::list result;
                for (const DistanceTransformSurvivalBin& bin : profiles.survivalProfile(node)) {
                    result.append(py::make_tuple(bin.squaredDistance, bin.countAtOrAbove));
                }
                return result;
            },
            "node_id"_a, "Return sorted `(squared_distance, count_at_or_above)` levels for one live node.");

    py::class_<DistanceTransformComputationBinding>(m, "DistanceTransform", py::module_local(false),
                                                    "Topology-only node distance-field and profile computation.")
        .def_static(
            "compute_node_field",
            [](const std::shared_ptr<PythonValuedMorphologicalTree>& valuedTree, NodeId node, const py::object& accuracyObject) {
                if (!valuedTree) {
                    throw std::invalid_argument("Distance-transform field computation requires a non-null valued tree.");
                }
                const DistanceTransformAccuracy accuracy = parseDistanceTransformAccuracy(accuracyObject);
                return valuedTree->visit([node, accuracy](const auto& concreteTree) {
                    return DistanceTransformComputation::computeNodeField(concreteTree->topology(), node, accuracy);
                });
            },
            "tree"_a.none(false), "node_id"_a, "accuracy"_a = py::str("exact"), "Compute the compact distance field of one live internal node.")
        .def_static(
            "compute_profiles",
            [](const std::shared_ptr<PythonValuedMorphologicalTree>& valuedTree, const py::object& accuracyObject) {
                if (!valuedTree) {
                    throw std::invalid_argument("Distance-transform profile computation requires a non-null valued tree.");
                }
                const DistanceTransformAccuracy accuracy = parseDistanceTransformAccuracy(accuracyObject);
                return valuedTree->visit(
                    [accuracy](const auto& concreteTree) { return DistanceTransformComputation::computeProfiles(concreteTree->topology(), accuracy); });
            },
            "tree"_a.none(false), "accuracy"_a = py::str("exact"), "Compute compact sparse distance histograms for every live internal node.");
}

} // namespace mmcfilters::pybindings
