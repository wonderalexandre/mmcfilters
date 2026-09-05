#pragma once

#include "../mmcfilters/filters/ExtinctionValues.hpp"
#include "PybindConversions.hpp"
#include "PythonValuedMorphologicalTree.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <concepts>
#include <optional>
#include <type_traits>
#include <variant>
namespace mmcfilters {

namespace py = pybind11;

/**
 * @brief Python-facing extinction selection policy converted to the native Real type on dispatch.
 */
struct ExtinctionSelectionPolicyPybind {
    /** @brief Enumerates the supported kind values. */
    enum class Kind { TopK, MinimumExtinction };

    /** @brief Selection mode represented by this policy. */
    Kind kind = Kind::TopK;
    /** @brief Number of top-ranked extrema retained by `TopK`. */
    int extremaToKeep = 0;
    /** @brief Minimum extinction accepted by `MinimumExtinction`. */
    double threshold = 0.0;

    /**
     * @brief Creates a policy that retains the top-k ranked extrema.
     *
     * @param extremaToKeep Number of highest-ranked extrema to retain.
     * @return Configured top-k selection policy.
     */
    static ExtinctionSelectionPolicyPybind byTopK(int extremaToKeep) {
        ExtinctionSelectionPolicyPybind policy;
        policy.kind = Kind::TopK;
        policy.extremaToKeep = extremaToKeep;
        return policy;
    }

    /**
     * @brief Creates a policy that retains extrema above a minimum extinction value.
     *
     * @param threshold Threshold applied by the operation.
     * @return Configured minimum-extinction selection policy.
     */
    static ExtinctionSelectionPolicyPybind byThreshold(double threshold) {
        ExtinctionSelectionPolicyPybind policy;
        policy.kind = Kind::MinimumExtinction;
        policy.threshold = threshold;
        return policy;
    }

    /**
     * @brief Converts to native.
     *
     * @return Converted to native.
     */
    template <std::floating_point Real> [[nodiscard]] ExtinctionSelectionPolicy<Real> toNative() const {
        switch (kind) {
        case Kind::TopK:
            return ExtinctionSelectionPolicy<Real>::byTopK(extremaToKeep);
        case Kind::MinimumExtinction:
            return ExtinctionSelectionPolicy<Real>::byThreshold(static_cast<Real>(threshold));
        }
        ExtinctionSelectionPolicy<Real> invalid;
        return invalid;
    }
};

/**
 * @brief Pybind11 wrapper exposing extinction-value computation to Python.
 */
class ExtinctionValuesPybind {
    /** @brief References the valued-tree owner used by the component. */
    std::shared_ptr<pybindings::PythonValuedMorphologicalTree> valuedTreeOwner_;
    /** @brief Defines all native extinction evaluators supported at runtime. */
    using ExtinctionStorage =
        std::variant<ExtinctionValues<std::uint8_t, float>, ExtinctionValues<std::uint8_t, double>, ExtinctionValues<ToSGrayLevel, float>,
                     ExtinctionValues<ToSGrayLevel, double>>;
    ExtinctionStorage extinction_; ///< Concrete extinction-value implementation.

    /**
     * @brief Creates an extinction-value computer for the runtime attribute type.
     *
     * @param valuedTree Valued tree.
     * @param attribute Attribute requested by the operation.
     * @return Extinction-value computer variant for the array element type.
     */
    template <std::floating_point Real, AltitudeValue Altitude>
    static ExtinctionValues<Altitude, Real> makeExtinction(ValuedMorphologicalTree<Altitude>& valuedTree, py::array attribute) {
        auto typed = pybind_utils::requireNodeAttributeArray<Real>(std::move(attribute), valuedTree.topology(), "attribute");
        return ExtinctionValues<Altitude, Real>(valuedTree, static_cast<const Real*>(typed.request().ptr));
    }

    /**
     * @brief Converts an edge-indexed saliency map to a Python dictionary.
     *
     * @param edgeMap Saliency values indexed by graph edge.
     * @return Dictionary containing the edge map values and grid metadata.
     */
    template <class Value> static py::dict edgeSaliencyMapToDict(EdgeSaliencyMap<Value>&& edgeMap) {
        const int numEdges = static_cast<int>(edgeMap.values.size());
        py::dict result;
        result["num_rows"] = edgeMap.numRows;
        result["num_columns"] = edgeMap.numColumns;
        result["adjacency_radius"] = edgeMap.adjacencyRadius;
        result["sources"] = pybind_utils::toNumpyOwned(std::move(edgeMap.sources), numEdges);
        result["targets"] = pybind_utils::toNumpyOwned(std::move(edgeMap.targets), numEdges);
        result["values"] = pybind_utils::toNumpyOwned(std::move(edgeMap.values), numEdges);
        return result;
    }

    /**
     * @brief Creates an extinction-value computer for the runtime attribute type.
     *
     * @param valuedTree Valued tree.
     * @param attribute Attribute requested by the operation.
     * @return Extinction-value computer variant for the array element type.
     */
    static ExtinctionStorage makeExtinction(pybindings::PythonValuedMorphologicalTree& valuedTree, py::array attribute) {
        if (pybind_utils::parseFloatingArrayDType(attribute, "attribute") == pybind_utils::FloatingDType::Float64) {
            return valuedTree.visit([&](auto& concreteTree) -> ExtinctionStorage {
                return makeExtinction<double>(*concreteTree, std::move(attribute));
            });
        }
        return valuedTree.visit([&](auto& concreteTree) -> ExtinctionStorage {
            return makeExtinction<float>(*concreteTree, std::move(attribute));
        });
    }

  public:
    /**
     * @brief Constructs `ExtinctionValuesPybind` from the supplied inputs.
     *
     * @param valuedTree Valued tree.
     * @param attribute Attribute requested by the operation.
     */
    ExtinctionValuesPybind(std::shared_ptr<pybindings::PythonValuedMorphologicalTree> valuedTree, py::array attribute)
        : valuedTreeOwner_(std::move(valuedTree)), extinction_(makeExtinction(*valuedTreeOwner_, std::move(attribute))) {}

    /**
     * @brief Returns map.
     *
     * @param selection Policy used to select extrema from the ranked candidates.
     * @param scorePolicy Policy used to score and rank candidate extrema.
     * @return Map.
     */
    py::array contourMap(const ExtinctionSelectionPolicyPybind& selection, ExtinctionContourScorePolicy scorePolicy) {
        return std::visit(
            [&selection, scorePolicy](auto& extinction) -> py::array {
                using Extinction = std::decay_t<decltype(extinction)>;
                using Real = typename Extinction::value_type;
                return pybind_utils::toNumpy(extinction.contourMap(selection.toNative<Real>(), scorePolicy));
            },
            extinction_);
    }

    /**
     * @brief Returns regional extrema py.
     *
     * @return Regional extrema py.
     */
    std::vector<py::tuple> getRegionalExtremaPy() {
        return std::visit(
            [](auto& extinction) {
                const auto& vec = extinction.getRegionalExtrema();
                std::vector<py::tuple> out;
                out.reserve(vec.size());
                for (const auto& item : vec) {
                    out.push_back(py::make_tuple(item.leaf, item.cutoffNode, item.extinction));
                }
                return out;
            },
            extinction_);
    }

    /**
     * @brief Returns extinction value attribute.
     *
     * @return Extinction value attribute.
     */
    py::array getExtinctionValueAttribute() {
        return std::visit(
            [](auto& extinction) -> py::array {
                const auto& cached = extinction.getExtinctionValueAttribute();
                std::vector<typename std::decay_t<decltype(extinction)>::value_type> valuation(cached.begin(), cached.end());
                const int numValues = static_cast<int>(valuation.size());
                return pybind_utils::toNumpyOwned(std::move(valuation), numValues);
            },
            extinction_);
    }

    /**
     * @brief Computes ranked extinction value attribute.
     *
     * @return Computed ranked extinction value attribute.
     */
    py::array computeRankedExtinctionValueAttribute() {
        return std::visit(
            [](auto& extinction) -> py::array {
                std::vector<int> valuation = extinction.computeRankedExtinctionValueAttribute();
                const int numValues = static_cast<int>(valuation.size());
                return pybind_utils::toNumpyOwned(std::move(valuation), numValues);
            },
            extinction_);
    }

    /**
     * @brief Computes formal saliency edge map.
     *
     * @param radius Neighbourhood radius.
     * @param ranked Ranked extrema from which the selection is produced.
     * @return Computed formal saliency edge map.
     */
    py::dict computeFormalSaliencyEdgeMap(std::optional<double> radius = std::nullopt, bool ranked = false) {
        const MorphologicalTree& tree = valuedTreeOwner_->topology();
        std::optional<RegularGridAdjacency2D> adjacency;
        if (radius.has_value()) {
            adjacency.emplace(pybind_utils::makeRegularGridAdjacency2D(tree.numRows(), tree.numColumns(), *radius,
                                                                       "ExtinctionValues.compute_formal_saliency_edge_map"));
        }

        return std::visit(
            [&](auto& extinction) -> py::dict {
                if (ranked) {
                    if (adjacency.has_value()) {
                        return edgeSaliencyMapToDict(extinction.computeRankedFormalSaliencyEdgeMap(*adjacency));
                    }
                    return edgeSaliencyMapToDict(extinction.computeRankedFormalSaliencyEdgeMap());
                }

                if (adjacency.has_value()) {
                    return edgeSaliencyMapToDict(extinction.computeFormalSaliencyEdgeMap(*adjacency));
                }
                return edgeSaliencyMapToDict(extinction.computeFormalSaliencyEdgeMap());
            },
            extinction_);
    }

    /**
     * @brief Computes the max-propagated extinction LCA projection.
     *
     * @param radius Optional projection-neighbourhood radius. When omitted, the
     * tree must carry one unambiguous stored adjacency.
     * @param ranked Whether to return dense effective-edge ranks instead of raw
     * extinction values.
     * @return Python dictionary containing the edge domain and projected values.
     */
    py::dict computeMonotoneExtinctionProjection(std::optional<double> radius = std::nullopt, bool ranked = false) {
        const MorphologicalTree& tree = valuedTreeOwner_->topology();
        std::optional<RegularGridAdjacency2D> adjacency;
        if (radius.has_value()) {
            adjacency.emplace(pybind_utils::makeRegularGridAdjacency2D(tree.numRows(), tree.numColumns(), *radius,
                                                                       "ExtinctionValues.compute_monotone_extinction_projection"));
        }

        return std::visit(
            [&](auto& extinction) -> py::dict {
                if (ranked) {
                    if (adjacency.has_value()) {
                        return edgeSaliencyMapToDict(extinction.computeRankedMonotoneExtinctionProjection(*adjacency));
                    }
                    return edgeSaliencyMapToDict(extinction.computeRankedMonotoneExtinctionProjection());
                }
                if (adjacency.has_value()) {
                    return edgeSaliencyMapToDict(extinction.computeMonotoneExtinctionProjection(*adjacency));
                }
                return edgeSaliencyMapToDict(extinction.computeMonotoneExtinctionProjection());
            },
            extinction_);
    }

    /**
     * @brief Applies filtering.
     *
     * @param selection Policy used to select extrema from the ranked candidates.
     * @return Image or array produced by the operation.
     */
    py::array filtering(const ExtinctionSelectionPolicyPybind& selection) {
        return std::visit(
            [&selection](auto& extinction) -> py::array {
                using Extinction = std::decay_t<decltype(extinction)>;
                using Real = typename Extinction::value_type;
                return pybind_utils::toNumpy(extinction.filtering(selection.toNative<Real>()));
            },
            extinction_);
    }
};

} // namespace mmcfilters
