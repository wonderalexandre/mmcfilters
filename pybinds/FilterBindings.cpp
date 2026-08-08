#include "ModuleBindings.hpp"

#include "AttributeFiltersPybind.hpp"
#include "DepthStableRegionComputerPybind.hpp"
#include "ExtinctionValuesPybind.hpp"
#include "UltimateAttributeOpeningPybind.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>
#include <utility>
#include <vector>

namespace mmcfilters::pybindings {

namespace py = pybind11;
using namespace pybind11::literals;

namespace {

void initExtinctionPolicyTypes(py::module_& m) {
    py::class_<ExtinctionSelectionPolicyPybind>(m, "ExtinctionSelectionPolicy", py::module_local(false),
                                                R"doc(Selection policy shared by extinction filtering and contour maps.)doc")
        .def_static("byTopK", &ExtinctionSelectionPolicyPybind::byTopK, "extremaToKeep"_a, "Select the strongest extrema by decreasing extinction ranking.")
        .def_static("byThreshold", &ExtinctionSelectionPolicyPybind::byThreshold, "threshold"_a,
                    "Select every extremum whose extinction value is greater than or equal to threshold.")
        .def_readonly("extremaToKeep", &ExtinctionSelectionPolicyPybind::extremaToKeep)
        .def_readonly("threshold", &ExtinctionSelectionPolicyPybind::threshold);

    py::enum_<ExtinctionContourScorePolicy>(m, "ExtinctionContourScorePolicy", py::module_local(false))
        .value("RankScore", ExtinctionContourScorePolicy::RankScore)
        .value("ExtinctionValue", ExtinctionContourScorePolicy::ExtinctionValue)
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

    py::class_<AttributeFiltersPybind>(m, "AttributeFilters", py::module_local(false),
                                       R"doc(Attribute-based filtering operators over a weighted morphological tree.

Attribute arrays are 1D `np.float32` or `np.float64` buffers indexed by dense
internal `NodeId`, C-contiguous, and with length `tree.numInternalNodeSlots`.
Returned images are 2D NumPy arrays on the original image domain.)doc")
        .def(py::init<std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>>>(), "tree"_a, "Create filtering operators over a `WeightedMorphologicalTree`.")
        .def("filteringByPruningMin", py::overload_cast<std::vector<bool>&>(&AttributeFiltersPybind::filteringByPruningMin), "criterion"_a,
             "Apply pruning-min filtering from a dense boolean keep/remove criterion.")
        .def(
            "filteringByPruningMin",
            [](AttributeFiltersPybind& self, py::array attr, double threshold) { return self.filteringByPruningMin(std::move(attr), threshold); }, "attr"_a,
            "threshold"_a, "Apply pruning-min filtering from a node-indexed floating-point attribute buffer.")
        .def("filteringByPruningMax", py::overload_cast<std::vector<bool>&>(&AttributeFiltersPybind::filteringByPruningMax), "criterion"_a,
             "Apply pruning-max filtering from a dense boolean keep/remove criterion.")
        .def("filteringDirectRule", py::overload_cast<std::vector<bool>&>(&AttributeFiltersPybind::filteringByDirectRule), "criterion"_a,
             "Apply the direct filtering rule from a dense boolean criterion.")
        .def("filteringSubtractiveRule", py::overload_cast<std::vector<bool>&>(&AttributeFiltersPybind::filteringBySubtractiveRule), "criterion"_a,
             "Apply the subtractive filtering rule from a dense boolean criterion.")
        .def("filteringSubtractiveScoreRule", py::overload_cast<std::vector<float>&>(&AttributeFiltersPybind::filteringBySubtractiveScoreRule), "scores"_a,
             "Apply subtractive-score filtering from dense per-node float scores.")
        .def(
            "filteringByPruningMax",
            [](AttributeFiltersPybind& self, py::array attr, double threshold) { return self.filteringByPruningMax(std::move(attr), threshold); }, "attr"_a,
            "threshold"_a, "Apply pruning-max filtering from a node-indexed floating-point attribute buffer.")
        .def(
            "filteringByViterbiRule",
            [](AttributeFiltersPybind& self, py::array attr, double threshold) { return self.filteringByViterbiRule(std::move(attr), threshold); }, "attr"_a,
            "threshold"_a, "Apply connected Viterbi filtering from a node-indexed floating-point attribute buffer.")
        .def(
            "filteringByExtinction",
            [](AttributeFiltersPybind& self, py::array attr, const ExtinctionSelectionPolicyPybind& selection) {
                return self.filteringByExtinctionValue(std::move(attr), selection);
            },
            "attr"_a, "selection"_a, "Filter by extinction using an explicit selection policy.")
        .def(
            "contourMapByExtinction",
            [](AttributeFiltersPybind& self, py::array attr, const ExtinctionSelectionPolicyPybind& selection, ExtinctionContourScorePolicy scorePolicy) {
                return self.contourMapByExtinctionValue(std::move(attr), selection, scorePolicy);
            },
            "attr"_a, "selection"_a, "scorePolicy"_a, "Build a contour-valued image from extinction values.")
        .def("getAdaptiveCriterion", &AttributeFiltersPybind::getAdaptiveCriterion, "criterion"_a, "delta"_a,
             "Expand a dense boolean criterion by an altitude ancestor/descendant delta.")
        .def("getAdaptiveCriterionByDepth", &AttributeFiltersPybind::getAdaptiveCriterionByDepth, "criterion"_a, "depthDelta"_a,
             "Expand a dense boolean criterion by a topological depth stability window.");
}

/**
 * @brief Registers depth-stable-region bindings in the Python module.
 *
 * @param m Python module receiving the bindings.
 */
void initDepthStableRegionComputer(py::module_& m) {
    py::class_<DepthStableRegionComputerPybind>(m, "DepthStableRegionComputer", py::module_local(false),
                                                R"doc(Topological depth-stability helper over a weighted morphological tree.

The optional attribute array must be 1D `np.float32` or `np.float64`,
C-contiguous, and indexed by dense internal `NodeId`. Without an explicit
attribute, the helper computes topology-only AREA internally. The reported
numeric score is variation; lower finite values are more stable. Result getters
require a successful computeByDepth call.)doc")
        .def(py::init<std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>>>(), "tree"_a, "Create a depth-stability helper using topology-only AREA.")
        .def(py::init<std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>>, py::array>(), "tree"_a, "attribute"_a,
             "Create a depth-stability helper from a node-indexed increasing attribute buffer.")
        .def("computeByDepth", &DepthStableRegionComputerPybind::computeByDepth, "depthDelta"_a, "Return a dense uint8 mask of strict local variation minima.")
        .def("getVariation", &DepthStableRegionComputerPybind::getVariation, "nodeId"_a, "Return the variation score for one node after computeByDepth.")
        .def("getVariations", &DepthStableRegionComputerPybind::getVariations, "Return the dense variation array after computeByDepth.")
        .def("nodeWithMinimumVariationInWindow", &DepthStableRegionComputerPybind::nodeWithMinimumVariationInWindow, "nodeId"_a,
             "Return the node with minimum finite variation in the current depth window.")
        .def("ascendantInStabilityWindow", &DepthStableRegionComputerPybind::ascendantInStabilityWindow, "nodeId"_a,
             "Return the ascendant used in the current depth window.")
        .def("descendantInStabilityWindow", &DepthStableRegionComputerPybind::descendantInStabilityWindow, "nodeId"_a,
             "Return the descendant used in the current depth window.")
        .def("getNumNodes", &DepthStableRegionComputerPybind::getNumNodes, "Return the number of nodes selected by the last computeByDepth call.")
        .def("setMaxVariation", &DepthStableRegionComputerPybind::setMaxVariation, "value"_a, "Set the maximum accepted variation value.")
        .def("setMinAttribute", &DepthStableRegionComputerPybind::setMinAttribute, "value"_a, "Set the lower accepted attribute bound.")
        .def("setMaxAttribute", &DepthStableRegionComputerPybind::setMaxAttribute, "value"_a, "Set the upper accepted attribute bound.");
}

/**
 * @brief Registers extinction-value bindings in the Python module.
 *
 * @param m Python module receiving the bindings.
 */
void initExtinctionValues(py::module_& m) {
    py::class_<ExtinctionValuesPybind>(m, "ExtinctionValues", py::module_local(false),
                                       R"doc(Extinction value utilities over a weighted morphological tree.

This class requires a globally monotone altitude-order capability, which the
standard max-tree and min-tree producers provide. The constructor attribute
array must be 1D `np.float32` or `np.float64`, C-contiguous, and indexed by
dense internal `NodeId`. The dominant extremum is reported with
`numpy.finfo(dtype).max`/`std::numeric_limits<Real>::max()`. The `contourMap`
method returns an image-domain contour visualization; use
`computeFormalSaliencyEdgeMap` for the persistence-based hierarchical-watershed
QFZ saliency projection.

Primary extinction reference: Alexandre Gonçalves Silva and Roberto de Alencar
Lotufo, "Efficient computation of new extinction values from extended component
tree," Pattern Recognition Letters 32(1):79-90, 2011,
https://doi.org/10.1016/j.patrec.2010.07.019. The formal hierarchical-watershed
projection follows Section 8.1 of Jean Cousty, Laurent Najman, Yukiko Kenmochi,
and Silvio Guimarães, "Hierarchical segmentations with graphs: quasi-flat zones,
minimum spanning trees, and saliency maps," Journal of Mathematical Imaging and
Vision 60(4):479-502, 2018, https://doi.org/10.1007/s10851-017-0768-7.)doc")
        .def(py::init<std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>>, py::array>(), "tree"_a, "attribute"_a,
             "Compute extinction values from a node-indexed attribute buffer.")
        .def("filtering", &ExtinctionValuesPybind::filtering, "selection"_a, "Reconstruct an image by applying an extinction selection policy.")
        .def("contourMap", &ExtinctionValuesPybind::contourMap, "selection"_a, "scorePolicy"_a,
             "Return a 2D contour visualization from an extinction selection policy.")
        .def("getExtinctionValueAttribute", &ExtinctionValuesPybind::getExtinctionValueAttribute,
             R"doc(Return the dense node extinction attribute induced by extinction records.

Each regional-extremum leaf receives its raw extinction value and every non-leaf
node receives the maximum extinction value among the extrema contained in its
subtree. The result is compatible with
`HierarchySaliencyMap.computeSaliencyEdgeMap`.)doc")
        .def("computeRankedExtinctionValueAttribute", &ExtinctionValuesPybind::computeRankedExtinctionValueAttribute,
             "Return dense non-negative integer ranks for the extinction-value attribute.")
        .def("computeFormalSaliencyEdgeMap", &ExtinctionValuesPybind::computeFormalSaliencyEdgeMap, "radius"_a = py::none(), "ranked"_a = false,
             R"doc(Return the edge-indexed hierarchical-watershed saliency map induced by extinction.

This method follows the Cousty persistence construction: it builds an
altitude-ordered MST/BPTAO, assigns each binary merge the minimum of its two
max-descendant extinction values, and returns the full-graph QFZ saliency map.
Set `ranked=True` for the canonical dense edge scale.)doc")
        .def("computeMonotoneExtinctionProjection", &ExtinctionValuesPybind::computeMonotoneExtinctionProjection, "radius"_a = py::none(), "ranked"_a = false,
             R"doc(Project the max-propagated extinction node attribute directly by LCA.

This method preserves the former `computeFormalSaliencyEdgeMap` behavior for
experiments that intentionally use the monotone node valuation. It is not the
Cousty hierarchical-watershed persistence construction.)doc")
        .def("getRegionalExtrema", &ExtinctionValuesPybind::getRegionalExtremaPy, "Return extinction tuples as (leafNodeId, cutoffNodeId, value).");
}

/**
 * @brief Registers ultimate-attribute-opening bindings in the Python module.
 *
 * @param m Python module receiving the bindings.
 */
void initUltimateAttributeOpening(py::module_& m) {
    py::class_<UltimateAttributeOpeningPybind>(m, "UltimateAttributeOpening", py::module_local(false),
                                               R"doc(Ultimate Attribute Opening over a weighted morphological tree.

The constructor attribute array must be 1D `np.float32` or `np.float64`,
C-contiguous, and indexed by dense internal `NodeId`. Call `execute` before
reading output images.)doc")
        .def(py::init<std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>>, py::array>(), "tree"_a, "attr"_a,
             "Create a UAO computation from a weighted tree and increasing attribute buffer.")
        .def("execute", &UltimateAttributeOpeningPybind::execute, "maxCriterion"_a, "Run UAO using all nodes as selectable candidates.")
        .def("executeWithMSER", &UltimateAttributeOpeningPybind::executeWithMSER, "maxCriterion"_a, "deltaMSER"_a,
             "Run UAO using an MSER-derived node-selection mask.")
        .def("executeWithDepthStability", &UltimateAttributeOpeningPybind::executeWithDepthStability, "maxCriterion"_a, "depthDelta"_a,
             "Run UAO using a depth-stability node-selection mask.")
        .def("getMaxContrastImage", &UltimateAttributeOpeningPybind::getMaxContrastImage, "Return the 2D uint8 maximum-contrast image.")
        .def("getAssociatedImage", &UltimateAttributeOpeningPybind::getAssociatedImage, "Return the 2D int32 associated attribute-index image.")
        .def("getAssociatedColoredImage", &UltimateAttributeOpeningPybind::getAssociatedColorImage,
             "Return a uint8 color rendering of the associated-index image.");
}

} // namespace mmcfilters::pybindings
