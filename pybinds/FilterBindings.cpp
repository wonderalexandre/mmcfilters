#include "ModuleBindings.hpp"

#include "AttributeFiltersPybind.hpp"
#include "DepthStableRegionComputerPybind.hpp"
#include "ExtinctionValuesPybind.hpp"
#include "UltimateAttributeOpeningPybind.hpp"
#include "../mmcfilters/filters/NodePreservationStability.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <concepts>
#include <cstdint>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace mmcfilters::pybindings {

namespace py = pybind11;
using namespace pybind11::literals;

namespace {

std::shared_ptr<PythonValuedMorphologicalTree> requireStabilityTreeOwner(
    std::shared_ptr<PythonValuedMorphologicalTree> tree, const char* context) {
    if (!tree) {
        throw std::invalid_argument(std::string(context) + " requires a non-null ValuedMorphologicalTree.");
    }
    return tree;
}

py::array_t<bool> toNumpyPreservationMask(const NodePreservationMask& mask) {
    py::array_t<bool> result(mask.size());
    bool* decisions = result.mutable_data();
    for (std::size_t index = 0; index < mask.size(); ++index) {
        decisions[index] = mask[index];
    }
    return result;
}

py::array_t<bool> adjustNodePreservationMaskByAltitudeStabilityPy(
    std::shared_ptr<PythonValuedMorphologicalTree> valuedTree, py::array nodePreservationMask,
    std::int64_t altitudeWindowRadius, IncompleteStabilityWindowPolicy incompleteWindowPolicy) {
    auto owner = requireStabilityTreeOwner(std::move(valuedTree), "adjust_node_preservation_mask_by_altitude_stability");
    const auto decisions = pybind_utils::requireNodePreservationArray(nodePreservationMask, owner->topology());
    return owner->visit([&](const auto& native) {
        using Tree = typename std::remove_cvref_t<decltype(native)>::element_type;
        using Altitude = typename Tree::AltitudeType;
        return toNumpyPreservationMask(adjustNodePreservationMaskByAltitudeStability(
            *native, decisions, static_cast<AltitudeDifference<Altitude>>(altitudeWindowRadius), incompleteWindowPolicy));
    });
}

py::array_t<bool> adjustNodePreservationMaskByDepthStabilityPy(
    std::shared_ptr<PythonValuedMorphologicalTree> valuedTree, py::array nodePreservationMask,
    int depthWindowRadius, IncompleteStabilityWindowPolicy incompleteWindowPolicy) {
    auto owner = requireStabilityTreeOwner(std::move(valuedTree), "adjust_node_preservation_mask_by_depth_stability");
    const auto decisions = pybind_utils::requireNodePreservationArray(nodePreservationMask, owner->topology());
    return toNumpyPreservationMask(adjustNodePreservationMaskByDepthStability(owner->topology(), decisions, depthWindowRadius, incompleteWindowPolicy));
}

void initExtinctionPolicyTypes(py::module_& m) {
    py::class_<ExtinctionSelectionPolicyPybind>(m, "ExtinctionSelectionPolicy", py::module_local(false),
                                                R"doc(Selection policy shared by extinction filtering and contour maps.)doc")
        .def_static("by_top_k", &ExtinctionSelectionPolicyPybind::byTopK, "extrema_to_keep"_a, "Select the strongest extrema by decreasing extinction ranking.")
        .def_static("by_threshold", &ExtinctionSelectionPolicyPybind::byThreshold, "threshold"_a,
                    "Select every extremum whose extinction value is greater than or equal to threshold.")
        .def_readonly("extrema_to_keep", &ExtinctionSelectionPolicyPybind::extremaToKeep)
        .def_readonly("threshold", &ExtinctionSelectionPolicyPybind::threshold);

    py::enum_<ExtinctionContourScorePolicy>(m, "ExtinctionContourScorePolicy", py::module_local(false))
        .value("RANK_SCORE", ExtinctionContourScorePolicy::RankScore)
        .value("EXTINCTION_VALUE", ExtinctionContourScorePolicy::ExtinctionValue)
        .export_values();
}

} // namespace

/**
 * @brief Registers attribute filters bindings in the Python module.
 *
 * @param m Python module receiving the bindings.
 */
void initAttributeFilters(py::module_& m) {
    initExtinctionPolicyTypes(m);

    py::enum_<IncompleteStabilityWindowPolicy>(m, "IncompleteStabilityWindowPolicy", py::module_local(false))
        .value("PRESERVE_INPUT_DECISION", IncompleteStabilityWindowPolicy::PreserveInputDecision)
        .export_values();

    m.def("adjust_node_preservation_mask_by_altitude_stability", &adjustNodePreservationMaskByAltitudeStabilityPy,
          "tree"_a, "node_preservation_mask"_a.noconvert(), "altitude_window_radius"_a,
          "incomplete_window_policy"_a = IncompleteStabilityWindowPolicy::PreserveInputDecision,
          "Relocate input rejections using a monotone altitude-distance stability window.");
    m.def("adjust_node_preservation_mask_by_depth_stability", &adjustNodePreservationMaskByDepthStabilityPy,
          "tree"_a, "node_preservation_mask"_a.noconvert(), "depth_window_radius"_a,
          "incomplete_window_policy"_a = IncompleteStabilityWindowPolicy::PreserveInputDecision,
          "Relocate input rejections using a topology-only edge-count stability window.");

    py::class_<DirectReconstruction>(m, "DirectReconstruction", py::module_local(false)).def(py::init<>());
    py::class_<SubtractiveResidueModulation>(m, "SubtractiveResidueModulation", py::module_local(false)).def(py::init<>());

    py::class_<AttributeFiltersPybind>(m, "AttributeFilters", py::module_local(false),
                                       R"doc(Attribute-based filtering operators over a valued morphological tree.

Preservation masks are C-contiguous 1D Boolean NumPy arrays with true meaning
preserve. Attribute arrays are 1D `np.float32` or `np.float64` buffers indexed by dense
internal `NodeId`, C-contiguous, and with length `tree.num_internal_node_slots`.
Returned images are 2D NumPy arrays on the original image domain.)doc")
        .def(py::init<std::shared_ptr<PythonValuedMorphologicalTree>>(), "tree"_a, "Create filtering operators over a `ValuedMorphologicalTree`.")
        .def("filtering_by_direct_rule", &AttributeFiltersPybind::filteringByDirectRule, "node_preservation_mask"_a.noconvert(),
             "Apply direct reconstruction; the root altitude is always retained.")
        .def("filtering_by_subtractive_rule", &AttributeFiltersPybind::filteringBySubtractiveRule, "node_preservation_mask"_a.noconvert(),
             "Gate every zero-baseline node residue with a Boolean preservation mask.")
        .def("filtering_by_soft_subtractive_rule", &AttributeFiltersPybind::filteringBySoftSubtractiveRule, "node_preservation_scores"_a.noconvert(),
             "Gate every zero-baseline node residue with finite scores in [0, 1].")
        .def("filtering_by_pruning_min_rule", py::overload_cast<py::array>(&AttributeFiltersPybind::filteringByPruningMin),
             "node_preservation_mask"_a.noconvert(), "Apply pruning-min filtering from explicit node-preservation decisions.")
        .def(
            "filtering_by_pruning_min_rule",
            [](AttributeFiltersPybind& self, py::array attribute, double threshold) { return self.filteringByPruningMin(std::move(attribute), threshold); }, "attribute"_a.noconvert(),
            "threshold"_a, "Apply pruning-min filtering from a node-indexed floating-point attribute buffer.")
        .def("filtering_by_pruning_max_rule", py::overload_cast<py::array>(&AttributeFiltersPybind::filteringByPruningMax),
             "node_preservation_mask"_a.noconvert(), "Apply pruning-max filtering from explicit node-preservation decisions.")
        .def(
            "filtering_by_pruning_max_rule",
            [](AttributeFiltersPybind& self, py::array attribute, double threshold) { return self.filteringByPruningMax(std::move(attribute), threshold); }, "attribute"_a.noconvert(),
            "threshold"_a, "Apply pruning-max filtering from a node-indexed floating-point attribute buffer.")
        .def(
            "filtering_by_viterbi_rule",
            [](AttributeFiltersPybind& self, py::array attribute, double threshold) { return self.filteringByViterbiRule(std::move(attribute), threshold); }, "attribute"_a.noconvert(),
            "threshold"_a, "Apply connected Viterbi filtering from a node-indexed floating-point attribute buffer.")
        .def(
            "filtering_by_extinction_rule",
            [](AttributeFiltersPybind& self, py::array attribute, const ExtinctionSelectionPolicyPybind& selection) {
                return self.filteringByExtinctionValue(std::move(attribute), selection);
            },
            "attribute"_a.noconvert(), "selection"_a, "Filter by extinction using an explicit selection policy.")
        .def(
            "contour_map_by_extinction",
            [](AttributeFiltersPybind& self, py::array attribute, const ExtinctionSelectionPolicyPybind& selection, ExtinctionContourScorePolicy scorePolicy) {
                return self.contourMapByExtinctionValue(std::move(attribute), selection, scorePolicy);
            },
            "attribute"_a.noconvert(), "selection"_a, "score_policy"_a, "Build a contour-valued image from extinction values.");
}

/**
 * @brief Registers depth-stable-region bindings in the Python module.
 *
 * @param m Python module receiving the bindings.
 */
void initDepthStableRegionComputer(py::module_& m) {
    py::class_<DepthStableRegionComputerPybind>(m, "DepthStableRegionComputer", py::module_local(false),
                                                R"doc(Topological depth-stability helper over a valued morphological tree.

The optional attribute array must be 1D `np.float32` or `np.float64`,
C-contiguous, and indexed by dense internal `NodeId`. Without an explicit
attribute, the helper computes topology-only AREA internally. The reported
numeric score is variation; lower finite values are more stable. Result getters
require a successful `compute_by_depth` call.)doc")
        .def(py::init<std::shared_ptr<PythonValuedMorphologicalTree>>(), "tree"_a, "Create a depth-stability helper using topology-only AREA.")
        .def(py::init<std::shared_ptr<PythonValuedMorphologicalTree>, py::array>(), "tree"_a, "attribute"_a.noconvert(),
             "Create a depth-stability helper from a node-indexed increasing attribute buffer.")
        .def("compute_by_depth", &DepthStableRegionComputerPybind::computeByDepth, "depth_window_radius"_a, "Return a dense uint8 mask of strict local variation minima.")
        .def("get_variation", &DepthStableRegionComputerPybind::getVariation, "node_id"_a, "Return the variation score for one node after `compute_by_depth`.")
        .def("get_variations", &DepthStableRegionComputerPybind::getVariations, "Return the dense variation array after `compute_by_depth`.")
        .def("node_with_minimum_variation_in_window", &DepthStableRegionComputerPybind::nodeWithMinimumVariationInWindow, "node_id"_a,
             "Return the node with minimum finite variation in the current depth window.")
        .def("ancestor_in_stability_window", &DepthStableRegionComputerPybind::ancestorInStabilityWindow, "node_id"_a,
             "Return the ancestor used in the current depth window.")
        .def("descendant_in_stability_window", &DepthStableRegionComputerPybind::descendantInStabilityWindow, "node_id"_a,
             "Return the descendant used in the current depth window.")
        .def("num_nodes", &DepthStableRegionComputerPybind::numNodes, "Return the number of nodes selected by the last `compute_by_depth` call.")
        .def("set_max_variation", &DepthStableRegionComputerPybind::setMaxVariation, "value"_a, "Set the maximum accepted variation value.")
        .def("set_min_attribute", &DepthStableRegionComputerPybind::setMinAttribute, "value"_a, "Set the lower accepted attribute bound.")
        .def("set_max_attribute", &DepthStableRegionComputerPybind::setMaxAttribute, "value"_a, "Set the upper accepted attribute bound.");
}

/**
 * @brief Registers extinction-value bindings in the Python module.
 *
 * @param m Python module receiving the bindings.
 */
void initExtinctionValues(py::module_& m) {
    py::class_<ExtinctionValuesPybind>(m, "ExtinctionValues", py::module_local(false),
                                       R"doc(Extinction value utilities over a valued morphological tree.

This class requires a globally monotone altitude-order capability, which the
standard max-tree and min-tree producers provide. The constructor attribute
array must be 1D `np.float32` or `np.float64`, C-contiguous, and indexed by
dense internal `NodeId`. The dominant extremum is reported with
`numpy.finfo(dtype).max`/`std::numeric_limits<Real>::max()`. The `contour_map`
method returns an image-domain contour visualization; use
`compute_formal_saliency_edge_map` for the persistence-based hierarchical-watershed
QFZ saliency projection.

Primary extinction reference: Alexandre Gonçalves Silva and Roberto de Alencar
Lotufo, "Efficient computation of new extinction values from extended component
tree," Pattern Recognition Letters 32(1):79-90, 2011,
https://doi.org/10.1016/j.patrec.2010.07.019. The formal hierarchical-watershed
projection follows Section 8.1 of Jean Cousty, Laurent Najman, Yukiko Kenmochi,
and Silvio Guimarães, "Hierarchical segmentations with graphs: quasi-flat zones,
minimum spanning trees, and saliency maps," Journal of Mathematical Imaging and
Vision 60(4):479-502, 2018, https://doi.org/10.1007/s10851-017-0768-7.)doc")
        .def(py::init<std::shared_ptr<PythonValuedMorphologicalTree>, py::array>(), "tree"_a, "attribute"_a.noconvert(),
             "Compute extinction values from a node-indexed attribute buffer.")
        .def("filtering", &ExtinctionValuesPybind::filtering, "selection"_a, "Reconstruct an image by applying an extinction selection policy.")
        .def("contour_map", &ExtinctionValuesPybind::contourMap, "selection"_a, "score_policy"_a,
             "Return a 2D contour visualization from an extinction selection policy.")
        .def("get_extinction_value_attribute", &ExtinctionValuesPybind::getExtinctionValueAttribute,
             R"doc(Return the dense node extinction attribute induced by extinction records.

Each regional-extremum leaf receives its raw extinction value and every non-leaf
node receives the maximum extinction value among the extrema contained in its
subtree. The result is compatible with
`HierarchySaliencyMap.compute_saliency_edge_map`.)doc")
        .def("compute_ranked_extinction_value_attribute", &ExtinctionValuesPybind::computeRankedExtinctionValueAttribute,
             "Return dense non-negative integer ranks for the extinction-value attribute.")
        .def("compute_formal_saliency_edge_map", &ExtinctionValuesPybind::computeFormalSaliencyEdgeMap, "radius"_a = py::none(), "ranked"_a = false,
             R"doc(Return the edge-indexed hierarchical-watershed saliency map induced by extinction.

This method follows the Cousty persistence construction: it builds an
altitude-ordered MST/BPTAO, assigns each binary merge the minimum of its two
max-descendant extinction values, and returns the full-graph QFZ saliency map.
Set `ranked=True` for the canonical dense edge scale.)doc")
        .def("compute_monotone_extinction_projection", &ExtinctionValuesPybind::computeMonotoneExtinctionProjection, "radius"_a = py::none(), "ranked"_a = false,
             R"doc(Project the max-propagated extinction node attribute directly by LCA.

Each transition edge receives the maximum descendant extinction of its
endpoints' lowest common ancestor. Edges within one smallest node receive zero.
Set `ranked=True` to return dense ranks of the realized transition values.)doc")
        .def("get_regional_extrema", &ExtinctionValuesPybind::getRegionalExtremaPy, "Return extinction tuples as (leaf_node_id, cutoff_node_id, value).");
}

/**
 * @brief Registers ultimate-attribute-opening bindings in the Python module.
 *
 * @param m Python module receiving the bindings.
 */
void initUltimateAttributeOpening(py::module_& m) {
    py::class_<UltimateAttributeOpeningPybind>(m, "UltimateAttributeOpening", py::module_local(false),
                                               R"doc(Ultimate Attribute Opening over a valued morphological tree.

The constructor attribute array must be 1D `np.float32` or `np.float64`,
C-contiguous, and indexed by dense internal `NodeId`. Call `execute` before
reading output images.)doc")
        .def(py::init<std::shared_ptr<PythonValuedMorphologicalTree>, py::array>(), "tree"_a, "attribute"_a.noconvert(),
             "Create a UAO computation from a valued tree and increasing attribute buffer.")
        .def("execute", &UltimateAttributeOpeningPybind::execute, "maximum_attribute_threshold"_a, "Run UAO using all nodes as selectable candidates.")
        .def("execute_with_mser", &UltimateAttributeOpeningPybind::executeWithMSER, "maximum_attribute_threshold"_a, "altitude_window_radius"_a,
             "Run UAO using an MSER-derived node-selection mask.")
        .def("execute_with_depth_stability", &UltimateAttributeOpeningPybind::executeWithDepthStability, "maximum_attribute_threshold"_a, "depth_window_radius"_a,
             "Run UAO using a depth-stability node-selection mask.")
        .def("get_max_contrast_image", &UltimateAttributeOpeningPybind::getMaxContrastImage, "Return the 2D uint8 maximum-contrast image.")
        .def("get_associated_image", &UltimateAttributeOpeningPybind::getAssociatedImage, "Return the 2D int32 associated attribute-index image.")
        .def("get_associated_colored_image", &UltimateAttributeOpeningPybind::getAssociatedColorImage,
             "Return a uint8 color rendering of the associated-index image.");
}

} // namespace mmcfilters::pybindings
