#include "ModuleBindings.hpp"

#include "AttributeComputationBindings.hpp"
#include "AttributeFiltersPybind.hpp"
#include "AttributeReconstructionFiltersPybind.hpp"
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
#include <type_traits>
#include <utility>
#include <vector>

namespace mmcfilters::pybindings {

namespace py = pybind11;
using namespace pybind11::literals;

namespace {

template <std::floating_point Real>
NodePreservationMask computeNodePreservationMaskPyTyped(py::array nodeAttributes, double threshold) {
    const py::buffer_info info = nodeAttributes.request();
    if (info.ndim != 1) {
        throw std::invalid_argument("node_attributes must be a 1D array.");
    }
    if (info.strides[0] != static_cast<py::ssize_t>(sizeof(Real))) {
        throw std::invalid_argument("node_attributes must be C-contiguous.");
    }
    const auto* values = static_cast<const Real*>(info.ptr);
    return computeNodePreservationMask(std::span<const Real>(values, static_cast<std::size_t>(info.shape[0])), static_cast<Real>(threshold));
}

NodePreservationMask computeNodePreservationMaskPy(py::array nodeAttributes, double threshold) {
    if (pybind_utils::parseFloatingArrayDType(nodeAttributes, "node_attributes") == pybind_utils::FloatingDType::Float64) {
        return computeNodePreservationMaskPyTyped<double>(std::move(nodeAttributes), threshold);
    }
    return computeNodePreservationMaskPyTyped<float>(std::move(nodeAttributes), threshold);
}

NodePreservationMask adjustNodePreservationMaskByAltitudeStabilityPy(
    std::shared_ptr<PythonValuedMorphologicalTree> valuedTree, const NodePreservationMask& nodePreservationMask,
    std::int64_t altitudeWindowRadius, IncompleteStabilityWindowPolicy incompleteWindowPolicy) {
    auto owner = requireReconstructionFilterOwner(std::move(valuedTree), "adjust_node_preservation_mask_by_altitude_stability");
    return owner->visit([&](const auto& native) {
        using Tree = typename std::remove_cvref_t<decltype(native)>::element_type;
        using Altitude = typename Tree::AltitudeType;
        return adjustNodePreservationMaskByAltitudeStability(
            *native, nodePreservationMask, static_cast<AltitudeDifference<Altitude>>(altitudeWindowRadius), incompleteWindowPolicy);
    });
}

NodePreservationMask adjustNodePreservationMaskByDepthStabilityPy(
    std::shared_ptr<PythonValuedMorphologicalTree> valuedTree, const NodePreservationMask& nodePreservationMask,
    int depthWindowRadius, IncompleteStabilityWindowPolicy incompleteWindowPolicy) {
    auto owner = requireReconstructionFilterOwner(std::move(valuedTree), "adjust_node_preservation_mask_by_depth_stability");
    return adjustNodePreservationMaskByDepthStability(owner->topology(), nodePreservationMask, depthWindowRadius, incompleteWindowPolicy);
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

    py::class_<NodePreservationMask>(m, "NodePreservationMask", py::module_local(false),
                                     "Dense Boolean node decisions where true means preserve the node.")
        .def(py::init<std::vector<bool>>(), "decisions"_a)
        .def("__len__", &NodePreservationMask::size)
        .def("__getitem__", [](const NodePreservationMask& mask, std::size_t index) {
            if (index >= mask.size()) {
                throw py::index_error();
            }
            return mask[index];
        })
        .def("to_list", [](const NodePreservationMask& mask) { return mask.decisions(); });

    py::class_<NodePruningMask>(m, "NodePruningMask", py::module_local(false),
                                "Dense Boolean node decisions where true means prune the node.")
        .def(py::init<std::vector<bool>>(), "decisions"_a)
        .def("__len__", &NodePruningMask::size)
        .def("__getitem__", [](const NodePruningMask& mask, std::size_t index) {
            if (index >= mask.size()) {
                throw py::index_error();
            }
            return mask[index];
        })
        .def("to_list", [](const NodePruningMask& mask) { return mask.decisions(); });

    m.def("to_node_pruning_mask", &toNodePruningMask, "node_preservation_mask"_a,
          "Explicitly complement preservation decisions into pruning decisions.");
    m.def("to_node_preservation_mask", &toNodePreservationMask, "node_pruning_mask"_a,
          "Explicitly complement pruning decisions into preservation decisions.");

    py::enum_<IncompleteStabilityWindowPolicy>(m, "IncompleteStabilityWindowPolicy", py::module_local(false))
        .value("PRESERVE_INPUT_DECISION", IncompleteStabilityWindowPolicy::PreserveInputDecision)
        .export_values();

    m.def("compute_node_preservation_mask", &computeNodePreservationMaskPy, "node_attributes"_a, "threshold"_a,
          "Threshold scalar node attributes with the inclusive >= preservation rule.");
    m.def("adjust_node_preservation_mask_by_altitude_stability", &adjustNodePreservationMaskByAltitudeStabilityPy,
          "tree"_a, "node_preservation_mask"_a, "altitude_window_radius"_a,
          "incomplete_window_policy"_a = IncompleteStabilityWindowPolicy::PreserveInputDecision,
          "Relocate input rejections using a monotone altitude-distance stability window.");
    m.def("adjust_node_preservation_mask_by_depth_stability", &adjustNodePreservationMaskByDepthStabilityPy,
          "tree"_a, "node_preservation_mask"_a, "depth_window_radius"_a,
          "incomplete_window_policy"_a = IncompleteStabilityWindowPolicy::PreserveInputDecision,
          "Relocate input rejections using a topology-only edge-count stability window.");

    py::class_<DirectReconstruction>(m, "DirectReconstruction", py::module_local(false)).def(py::init<>());
    py::class_<SubtractiveResidueModulation>(m, "SubtractiveResidueModulation", py::module_local(false)).def(py::init<>());

    py::class_<DirectAttributeFilterPybind>(m, "DirectAttributeFilter", py::module_local(false))
        .def(py::init<std::shared_ptr<PythonValuedMorphologicalTree>>(), "tree"_a)
        .def("apply_direct_attribute_filter", &DirectAttributeFilterPybind::applyDirectAttributeFilter, "node_preservation_mask"_a,
             "Apply direct reconstruction; the root must be preserved.")
        .def("apply", &DirectAttributeFilterPybind::applyDirectAttributeFilter, "node_preservation_mask"_a,
             "Short alias of `apply_direct_attribute_filter`.");

    py::class_<SubtractiveAttributeFilterPybind>(m, "SubtractiveAttributeFilter", py::module_local(false))
        .def(py::init<std::shared_ptr<PythonValuedMorphologicalTree>>(), "tree"_a)
        .def("apply_subtractive_attribute_filter", &SubtractiveAttributeFilterPybind::applySubtractiveAttributeFilter,
             "node_preservation_mask"_a, "Gate every zero-baseline node residue with a Boolean preservation mask.")
        .def("apply", &SubtractiveAttributeFilterPybind::applySubtractiveAttributeFilter, "node_preservation_mask"_a,
             "Short alias of `apply_subtractive_attribute_filter`.");

    py::class_<SoftSubtractiveAttributeFilterPybind>(m, "SoftSubtractiveAttributeFilter", py::module_local(false))
        .def(py::init<std::shared_ptr<PythonValuedMorphologicalTree>>(), "tree"_a)
        .def("apply_soft_subtractive_attribute_filter", &SoftSubtractiveAttributeFilterPybind::applySoftSubtractiveAttributeFilter,
             "node_preservation_scores"_a, "Gate every zero-baseline node residue with finite scores in [0, 1].")
        .def("apply", &SoftSubtractiveAttributeFilterPybind::applySoftSubtractiveAttributeFilter, "node_preservation_scores"_a,
             "Short alias of `apply_soft_subtractive_attribute_filter`.");

    py::class_<AttributeFiltersPybind>(m, "AttributeFilters", py::module_local(false),
                                       R"doc(Attribute-based filtering operators over a valued morphological tree.

Attribute arrays are 1D `np.float32` or `np.float64` buffers indexed by dense
internal `NodeId`, C-contiguous, and with length `tree.num_internal_node_slots`.
Returned images are 2D NumPy arrays on the original image domain.)doc")
        .def(py::init<std::shared_ptr<PythonValuedMorphologicalTree>>(), "tree"_a, "Create filtering operators over a `ValuedMorphologicalTree`.")
        .def("filtering_by_pruning_min", py::overload_cast<const NodePreservationMask&>(&AttributeFiltersPybind::filteringByPruningMin),
             "node_preservation_mask"_a, "Apply pruning-min filtering from explicit node-preservation decisions.")
        .def(
            "filtering_by_pruning_min",
            [](AttributeFiltersPybind& self, py::array attr, double threshold) { return self.filteringByPruningMin(std::move(attr), threshold); }, "attr"_a,
            "threshold"_a, "Apply pruning-min filtering from a node-indexed floating-point attribute buffer.")
        .def(
            "filtering_by_pruning_min",
            [](AttributeFiltersPybind& self, const py::object& attr, double threshold) {
                return self.filteringByPruningMin(attribute_computation::attributeBufferFor(self.treeOwner(), attr), threshold);
            },
            "attr"_a, "threshold"_a, "Accepts an `Attribute` value or its symbolic name and computes the buffer internally.")
        .def("filtering_by_pruning_max", py::overload_cast<const NodePreservationMask&>(&AttributeFiltersPybind::filteringByPruningMax),
             "node_preservation_mask"_a, "Apply pruning-max filtering from explicit node-preservation decisions.")
        .def(
            "filtering_by_pruning_max",
            [](AttributeFiltersPybind& self, py::array attr, double threshold) { return self.filteringByPruningMax(std::move(attr), threshold); }, "attr"_a,
            "threshold"_a, "Apply pruning-max filtering from a node-indexed floating-point attribute buffer.")
        .def(
            "filtering_by_pruning_max",
            [](AttributeFiltersPybind& self, const py::object& attr, double threshold) {
                return self.filteringByPruningMax(attribute_computation::attributeBufferFor(self.treeOwner(), attr), threshold);
            },
            "attr"_a, "threshold"_a, "Accepts an `Attribute` value or its symbolic name and computes the buffer internally.")
        .def(
            "filtering_by_viterbi_rule",
            [](AttributeFiltersPybind& self, py::array attr, double threshold) { return self.filteringByViterbiRule(std::move(attr), threshold); }, "attr"_a,
            "threshold"_a, "Apply connected Viterbi filtering from a node-indexed floating-point attribute buffer.")
        .def(
            "filtering_by_viterbi_rule",
            [](AttributeFiltersPybind& self, const py::object& attr, double threshold) {
                return self.filteringByViterbiRule(attribute_computation::attributeBufferFor(self.treeOwner(), attr), threshold);
            },
            "attr"_a, "threshold"_a, "Accepts an `Attribute` value or its symbolic name and computes the buffer internally.")
        .def(
            "filtering_by_extinction",
            [](AttributeFiltersPybind& self, py::array attr, const ExtinctionSelectionPolicyPybind& selection) {
                return self.filteringByExtinctionValue(std::move(attr), selection);
            },
            "attr"_a, "selection"_a, "Filter by extinction using an explicit selection policy.")
        .def(
            "filtering_by_extinction",
            [](AttributeFiltersPybind& self, const py::object& attr, const ExtinctionSelectionPolicyPybind& selection) {
                return self.filteringByExtinctionValue(attribute_computation::attributeBufferFor(self.treeOwner(), attr), selection);
            },
            "attr"_a, "selection"_a, "Accepts an `Attribute` value or its symbolic name and computes the buffer internally.")
        .def(
            "contour_map_by_extinction",
            [](AttributeFiltersPybind& self, py::array attr, const ExtinctionSelectionPolicyPybind& selection, ExtinctionContourScorePolicy scorePolicy) {
                return self.contourMapByExtinctionValue(std::move(attr), selection, scorePolicy);
            },
            "attr"_a, "selection"_a, "score_policy"_a, "Build a contour-valued image from extinction values.");
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
        .def(
            "contour_map_by_extinction",
            [](AttributeFiltersPybind& self, const py::object& attr, const ExtinctionSelectionPolicyPybind& selection,
               ExtinctionContourScorePolicy scorePolicy) {
                return self.contourMapByExtinctionValue(attribute_computation::attributeBufferFor(self.treeOwner(), attr), selection, scorePolicy);
            },
            "attr"_a, "selection"_a, "score_policy"_a, "Accepts an `Attribute` value or its symbolic name and computes the buffer internally.")
        .def(py::init<std::shared_ptr<PythonValuedMorphologicalTree>>(), "tree"_a, "Create a depth-stability helper using topology-only AREA.")
        .def(py::init<std::shared_ptr<PythonValuedMorphologicalTree>, py::array>(), "tree"_a, "attribute"_a,
             "Create a depth-stability helper from a node-indexed increasing attribute buffer.")
        .def(py::init([](std::shared_ptr<PythonValuedMorphologicalTree> tree, const py::object& attribute) {
                 py::array buffer = attribute_computation::attributeBufferFor(tree, attribute);
                 return std::make_unique<DepthStableRegionComputerPybind>(std::move(tree), std::move(buffer));
             }),
             "tree"_a, "attribute"_a, "Accepts an `Attribute` value or its symbolic name and computes the buffer internally.")
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
        .def(py::init<std::shared_ptr<PythonValuedMorphologicalTree>, py::array>(), "tree"_a, "attribute"_a,
             "Compute extinction values from a node-indexed attribute buffer.")
        .def(py::init([](std::shared_ptr<PythonValuedMorphologicalTree> tree, const py::object& attribute) {
                 py::array buffer = attribute_computation::attributeBufferFor(tree, attribute);
                 return std::make_unique<ExtinctionValuesPybind>(std::move(tree), std::move(buffer));
             }),
             "tree"_a, "attribute"_a, "Accepts an `Attribute` value or its symbolic name and computes the buffer internally.")
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

This method preserves the former `compute_formal_saliency_edge_map` behavior for
experiments that intentionally use the monotone node valuation. It is not the
Cousty hierarchical-watershed persistence construction.)doc")
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
        .def(py::init<std::shared_ptr<PythonValuedMorphologicalTree>, py::array>(), "tree"_a, "attr"_a,
             "Create a UAO computation from a valued tree and increasing attribute buffer.")
        .def(py::init([](std::shared_ptr<PythonValuedMorphologicalTree> tree, const py::object& attr) {
                 py::array buffer = attribute_computation::attributeBufferFor(tree, attr);
                 return std::make_unique<UltimateAttributeOpeningPybind>(std::move(tree), std::move(buffer));
             }),
             "tree"_a, "attr"_a, "Accepts an `Attribute` value or its symbolic name and computes the buffer internally.")
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
