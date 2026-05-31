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

void initAttributeFilters(py::module_& m) {
    py::class_<AttributeFiltersPybind>(m, "AttributeFilters", py::module_local(false),
        R"doc(Attribute-based filtering operators over a weighted morphological tree.

Attribute arrays are 1D `np.float32` or `np.float64` buffers indexed by dense
internal `NodeId`, C-contiguous, and with length `tree.numInternalNodeSlots`.
Returned images are 2D NumPy arrays on the original image domain.)doc")
        .def(py::init<std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>>>(),
            "tree"_a,
            "Create filtering operators over a `WeightedMorphologicalTree`.")
        .def("filteringByPruningMin", py::overload_cast<std::vector<bool>&>(&AttributeFiltersPybind::filteringByPruningMin),
            "criterion"_a,
            "Apply pruning-min filtering from a dense boolean keep/remove criterion.")
        .def("filteringByPruningMin", [](AttributeFiltersPybind& self, py::array attr, double threshold) {
            return self.filteringByPruningMin(std::move(attr), threshold);
        }, "attr"_a, "threshold"_a,
            "Apply pruning-min filtering from a node-indexed floating-point attribute buffer.")
        .def("filteringByPruningMax", py::overload_cast<std::vector<bool>&>(&AttributeFiltersPybind::filteringByPruningMax),
            "criterion"_a,
            "Apply pruning-max filtering from a dense boolean keep/remove criterion.")
        .def("filteringDirectRule", py::overload_cast<std::vector<bool>&>(&AttributeFiltersPybind::filteringByDirectRule),
            "criterion"_a,
            "Apply the direct filtering rule from a dense boolean criterion.")
        .def("filteringSubtractiveRule", py::overload_cast<std::vector<bool>&>(&AttributeFiltersPybind::filteringBySubtractiveRule),
            "criterion"_a,
            "Apply the subtractive filtering rule from a dense boolean criterion.")
        .def("filteringSubtractiveScoreRule", py::overload_cast<std::vector<float>&>(&AttributeFiltersPybind::filteringBySubtractiveScoreRule),
            "scores"_a,
            "Apply subtractive-score filtering from dense per-node float scores.")
        .def("filteringByPruningMax", [](AttributeFiltersPybind& self, py::array attr, double threshold) {
            return self.filteringByPruningMax(std::move(attr), threshold);
        }, "attr"_a, "threshold"_a,
            "Apply pruning-max filtering from a node-indexed floating-point attribute buffer.")
        .def("filteringByViterbiRule", [](AttributeFiltersPybind& self, py::array attr, double threshold) {
            return self.filteringByViterbiRule(std::move(attr), threshold);
        }, "attr"_a, "threshold"_a,
            "Apply connected Viterbi filtering from a node-indexed floating-point attribute buffer.")
        .def("filteringByExtinction", [](AttributeFiltersPybind& self, py::array attr, int leafToKeep) {
            return self.filteringByExtinctionValue(std::move(attr), leafToKeep);
        }, "attr"_a, "leafToKeep"_a,
            "Filter by keeping the strongest extinction extrema.")
        .def("saliencyMapByExtinction", [](AttributeFiltersPybind& self, py::array attr, int leafToKeep, bool unweighted) {
            return self.saliencyMapByExtinctionValue(std::move(attr), leafToKeep, unweighted);
        }, "attr"_a, "leafToKeep"_a, "unweighted"_a = false,
            "Build a contour saliency map from extinction values.")
        .def("getAdaptiveCriterion", &AttributeFiltersPybind::getAdaptiveCriterion,
            "criterion"_a,
            "delta"_a,
            "Expand a dense boolean criterion by an altitude ancestor/descendant delta.")
        .def("getAdaptiveCriterionByDepth", &AttributeFiltersPybind::getAdaptiveCriterionByDepth,
            "criterion"_a,
            "depthDelta"_a,
            "Expand a dense boolean criterion by a topological depth stability window.");
	}

	void initDepthStableRegionComputer(py::module_& m) {
	    py::class_<DepthStableRegionComputerPybind>(m, "DepthStableRegionComputer", py::module_local(false),
	        R"doc(Topological depth-stability helper over a weighted morphological tree.

The optional attribute array must be 1D `np.float32` or `np.float64`,
C-contiguous, and indexed by dense internal `NodeId`. Without an explicit
attribute, the helper computes topology-only AREA internally. The reported
numeric score is variation; lower finite values are more stable.)doc")
	        .def(py::init<std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>>>(),
	            "tree"_a,
	            "Create a depth-stability helper using topology-only AREA.")
	        .def(py::init<std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>>, py::array>(),
	            "tree"_a,
	            "attribute"_a,
	            "Create a depth-stability helper from a node-indexed increasing attribute buffer.")
	        .def("computeByDepth", &DepthStableRegionComputerPybind::computeByDepth,
	            "depthDelta"_a,
	            "Return a dense uint8 mask of strict local variation minima.")
	        .def("getVariation", &DepthStableRegionComputerPybind::getVariation,
	            "nodeId"_a,
	            "Return the variation score for one node after computeByDepth.")
	        .def("getVariations", &DepthStableRegionComputerPybind::getVariations,
	            "Return the dense variation array after computeByDepth.")
	        .def("nodeWithMinimumVariationInWindow", &DepthStableRegionComputerPybind::nodeWithMinimumVariationInWindow,
	            "nodeId"_a,
	            "Return the node with minimum finite variation in the current depth window.")
	        .def("ascendantInStabilityWindow", &DepthStableRegionComputerPybind::ascendantInStabilityWindow,
	            "nodeId"_a,
	            "Return the ascendant used in the current depth window.")
	        .def("descendantInStabilityWindow", &DepthStableRegionComputerPybind::descendantInStabilityWindow,
	            "nodeId"_a,
	            "Return the descendant used in the current depth window.")
	        .def("getNumNodes", &DepthStableRegionComputerPybind::getNumNodes,
	            "Return the number of nodes selected by the last computeByDepth call.")
	        .def("setMaxVariation", &DepthStableRegionComputerPybind::setMaxVariation,
	            "value"_a,
	            "Set the maximum accepted variation value.")
	        .def("setMinAttribute", &DepthStableRegionComputerPybind::setMinAttribute,
	            "value"_a,
	            "Set the lower accepted attribute bound.")
	        .def("setMaxAttribute", &DepthStableRegionComputerPybind::setMaxAttribute,
	            "value"_a,
	            "Set the upper accepted attribute bound.");
	}

	void initExtinctionValues(py::module_& m) {
	    py::class_<ExtinctionValuesPybind>(m, "ExtinctionValues", py::module_local(false),
	        R"doc(Extinction value utilities over a weighted morphological tree.

	This class is currently defined for max-trees and min-trees. The constructor
	attribute array must be 1D `np.float32` or `np.float64`, C-contiguous, and
	indexed by dense internal `NodeId`. The dominant extremum is reported with
	`numpy.finfo(dtype).max`/`std::numeric_limits<Real>::max()`.)doc")
        .def(py::init<std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>>, py::array>(),
            "tree"_a,
            "attribute"_a,
            "Compute extinction values from a node-indexed attribute buffer.")
        .def("filtering", &ExtinctionValuesPybind::filtering,
            "leafToKeep"_a,
            "Reconstruct an image by keeping the strongest extrema.")
        .def("saliencyMap", &ExtinctionValuesPybind::saliencyMap, "leafToKeep"_a, "unweighted"_a = true,
            "Return a 2D float saliency map from the strongest extrema.")
        .def("getExtinctionValues", &ExtinctionValuesPybind::getExtinctionValuesPy,
            "Return extinction tuples as (leafNodeId, cutoffNodeId, value).");
}

void initUltimateAttributeOpening(py::module_& m) {
    py::class_<UltimateAttributeOpeningPybind>(m, "UltimateAttributeOpening", py::module_local(false),
        R"doc(Ultimate Attribute Opening over a weighted morphological tree.

The constructor attribute array must be 1D `np.float32` or `np.float64`,
C-contiguous, and indexed by dense internal `NodeId`. Call `execute` before
reading output images.)doc")
        .def(py::init<std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>>, py::array>(),
            "tree"_a,
            "attr"_a,
            "Create a UAO computation from a weighted tree and increasing attribute buffer.")
        .def("execute", &UltimateAttributeOpeningPybind::execute,
            "maxCriterion"_a,
            "Run UAO using all nodes as selectable candidates.")
        .def("executeWithMSER", &UltimateAttributeOpeningPybind::executeWithMSER,
            "maxCriterion"_a,
            "deltaMSER"_a,
            "Run UAO using an MSER-derived node-selection mask.")
        .def("executeWithDepthStability", &UltimateAttributeOpeningPybind::executeWithDepthStability,
            "maxCriterion"_a,
            "depthDelta"_a,
            "Run UAO using a depth-stability node-selection mask.")
        .def("getMaxContrastImage", &UltimateAttributeOpeningPybind::getMaxContrastImage,
            "Return the 2D uint8 maximum-contrast image.")
        .def("getAssociatedImage", &UltimateAttributeOpeningPybind::getAssociatedImage,
            "Return the 2D int32 associated attribute-index image.")
        .def("getAssociatedColoredImage", &UltimateAttributeOpeningPybind::getAssociatedColorImage,
            "Return a uint8 color rendering of the associated-index image.");
}

} // namespace mmcfilters::pybindings
