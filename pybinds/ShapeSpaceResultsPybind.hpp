#pragma once

#include "EdgeMapPybind.hpp"
#include "../mmcfilters/trees/saliency/ShapeSpaceSaliency.hpp"

namespace mmcfilters::pybindings {

/** @brief Immutable scalar record for one shape-space regional extremum. */
struct ShapeSpaceExtremumPybind {
    NodeId representative; ///< Stable representative hierarchy node.
    double birthLevel; ///< Attribute level at which the extremum appears.
    double deathLevel; ///< Attribute level at which the extremum is absorbed.
    double extinction; ///< Absolute difference between death and birth levels.

    /**
     * @brief Copies scalar extremum data without changing its numeric value.
     * @param value Native extremum record in its computed precision.
     */
    template <std::floating_point Real> explicit ShapeSpaceExtremumPybind(const ShapeSpaceExtremum<Real>& value)
        : representative(value.representative), birthLevel(value.birthLevel), deathLevel(value.deathLevel), extinction(value.extinction) {}
};

/**
 * @brief Materializes immutable extremum records once for a computed result.
 * @param extrema Native extrema in computation order.
 * @return Immutable tuple of independent scalar records.
 */
template <std::floating_point Real> py::tuple shapeSpaceExtremaTuple(const std::vector<ShapeSpaceExtremum<Real>>& extrema) {
    py::tuple result(extrema.size());
    for (std::size_t index = 0; index < extrema.size(); ++index) {
        result[index] = py::cast(ShapeSpaceExtremumPybind(extrema[index]));
    }
    return result;
}

/** @brief Owns immutable extrema and dense-domain shape-space scores. */
class ShapeSpaceExtinctionResultPybind {
    std::variant<std::vector<float>, std::vector<double>> nodeScores_; ///< Scores in their computed precision.
    py::tuple extrema_; ///< Cached immutable scalar records in computation order.

  public:
    /**
     * @brief Takes ownership of computed shape-space scores.
     * @param result Native extrema and scores transferred to the result object.
     */
    template <std::floating_point Real> explicit ShapeSpaceExtinctionResultPybind(ShapeSpaceExtinctionResult<Real>&& result)
        : nodeScores_(std::move(result.nodeScores)), extrema_(shapeSpaceExtremaTuple(result.extrema)) {}

    /**
     * @brief Reads the cached immutable extrema.
     * @return Tuple of independent scalar records.
     */
    py::tuple extrema() const { return extrema_; }

    /**
     * @brief Reads scores using their computed numeric type.
     * @param visitor Function applied to the immutable score buffer.
     * @return Result of the visitor.
     */
    template <class Visitor> decltype(auto) visitScores(Visitor&& visitor) const {
        return std::visit(std::forward<Visitor>(visitor), nodeScores_);
    }
};

/** @brief Owns Xu extrema, dense scores, and their contour saliency map. */
class ShapeSpaceSaliencyResultPybind {
    ShapeSpaceExtinctionResultPybind extinction_; ///< Owned extrema and node scores.
    EdgeMapPybind edgeMap_; ///< Owned maximum-on-contours saliency map.

  public:
    /**
     * @brief Takes ownership of a complete Xu saliency result.
     * @param result Native extrema, scores, and edge map transferred without buffer copies.
     */
    template <std::floating_point Real> explicit ShapeSpaceSaliencyResultPybind(ShapeSpaceSaliencyResult<Real>&& result)
        : extinction_(ShapeSpaceExtinctionResult<Real>{std::move(result.extrema), std::move(result.nodeScores)}),
          edgeMap_(std::move(result.edgeMap)) {}

    /**
     * @brief Reads the immutable extrema and score component.
     * @return Owned extinction result component.
     */
    const ShapeSpaceExtinctionResultPybind& extinction() const { return extinction_; }

    /**
     * @brief Reads the owned saliency map without copying its buffers.
     * @return Immutable edge-map component.
     */
    const EdgeMapPybind& edgeMap() const { return edgeMap_; }
};

/**
 * @brief Registers immutable Xu scalar records and composite results.
 * @param m Python module receiving the result types.
 */
inline void initShapeSpaceResultBindings(py::module_& m) {
    py::class_<ShapeSpaceExtremumPybind>(m, "ShapeSpaceExtremum", py::module_local(false),
                                       "Immutable representative node, birth level, death level, and extinction.")
        .def_readonly("representative", &ShapeSpaceExtremumPybind::representative)
        .def_readonly("birth_level", &ShapeSpaceExtremumPybind::birthLevel)
        .def_readonly("death_level", &ShapeSpaceExtremumPybind::deathLevel)
        .def_readonly("extinction", &ShapeSpaceExtremumPybind::extinction);

    py::class_<ShapeSpaceExtinctionResultPybind>(m, "ShapeSpaceExtinctionResult", py::module_local(false),
                                               "Owned shape-space extrema and read-only NumPy node scores.")
        .def_property_readonly("extrema", &ShapeSpaceExtinctionResultPybind::extrema)
        .def_property_readonly("node_scores", [](py::object self) {
            return self.cast<const ShapeSpaceExtinctionResultPybind&>().visitScores([&](const auto& scores) {
                return readOnlyResultArray(scores, self);
            });
        });

    py::class_<ShapeSpaceSaliencyResultPybind>(m, "ShapeSpaceSaliencyResult", py::module_local(false),
                                             "Owned Xu extrema, read-only NumPy node scores, and an EdgeMap.")
        .def_property_readonly("extrema", [](const ShapeSpaceSaliencyResultPybind& self) { return self.extinction().extrema(); })
        .def_property_readonly("node_scores", [](py::object self) {
            return self.cast<const ShapeSpaceSaliencyResultPybind&>().extinction().visitScores([&](const auto& scores) {
                return readOnlyResultArray(scores, self);
            });
        })
        .def_property_readonly("edge_map", &ShapeSpaceSaliencyResultPybind::edgeMap, py::return_value_policy::reference_internal);
}

} // namespace mmcfilters::pybindings
