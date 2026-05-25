#include "ModuleBindings.hpp"

#include "AttributeFiltersPybind.hpp"
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
    using FloatArray = py::array_t<float, py::array::c_style | py::array::forcecast>;

    py::class_<AttributeFiltersPybind>(m, "AttributeFilters", py::module_local(false),
        R"doc(Attribute-based filtering operators over a weighted morphological tree.

All attribute arrays are 1D `np.float32` buffers indexed by dense internal
`NodeId`, C-contiguous, and with length `tree.numInternalNodeSlots`. Returned
images are 2D NumPy arrays on the original image domain.)doc")
        .def(py::init<std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>>>(),
            "tree"_a,
            "Create filtering operators over a `WeightedMorphologicalTree`.")
        .def("filteringMin", [](AttributeFiltersPybind& self, FloatArray attr, float threshold) {
            return self.filteringByPruningMin(std::move(attr), threshold);
        }, "attr"_a, "threshold"_a,
            "Apply pruning-min filtering from a node-indexed float attribute buffer.")
        .def("filteringMin", py::overload_cast<std::vector<bool>&>(&AttributeFiltersPybind::filteringByPruningMin),
            "criterion"_a,
            "Apply pruning-min filtering from a dense boolean keep/remove criterion.")
        .def("filteringMax", py::overload_cast<std::vector<bool>&>(&AttributeFiltersPybind::filteringByPruningMax),
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
        .def("filteringMax", [](AttributeFiltersPybind& self, FloatArray attr, float threshold) {
            return self.filteringByPruningMax(std::move(attr), threshold);
        }, "attr"_a, "threshold"_a,
            "Apply pruning-max filtering from a node-indexed float attribute buffer.")
        .def("filteringByExtinction", [](AttributeFiltersPybind& self, FloatArray attr, int leafToKeep) {
            return self.filteringByExtinctionValue(std::move(attr), leafToKeep);
        }, "attr"_a, "leafToKeep"_a,
            "Filter by keeping the strongest extinction extrema.")
        .def("saliencyMapByExtinction", [](AttributeFiltersPybind& self, FloatArray attr, int leafToKeep, bool unweighted) {
            return self.saliencyMapByExtinctionValue(std::move(attr), leafToKeep, unweighted);
        }, "attr"_a, "leafToKeep"_a, "unweighted"_a = false,
            "Build a contour saliency map from extinction values.")
        .def("getAdaptiveCriterion", &AttributeFiltersPybind::getAdaptiveCriterion,
            "criterion"_a,
            "delta"_a,
            "Expand a dense boolean criterion by an ancestor/descendant delta.");
}

void initExtinctionValues(py::module_& m) {
    using FloatArray = py::array_t<float, py::array::c_style | py::array::forcecast>;

    py::class_<ExtinctionValuesPybind>(m, "ExtinctionValues", py::module_local(false),
        R"doc(Extinction value utilities over a weighted morphological tree.

The constructor attribute array must be 1D `np.float32`, C-contiguous, and
indexed by dense internal `NodeId`.)doc")
        .def(py::init<std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>>, FloatArray>(),
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
    using FloatArray = py::array_t<float, py::array::c_style | py::array::forcecast>;

    py::class_<UltimateAttributeOpeningPybind>(m, "UltimateAttributeOpening", py::module_local(false),
        R"doc(Ultimate Attribute Opening over a weighted morphological tree.

The constructor attribute array must be 1D `np.float32`, C-contiguous, and
indexed by dense internal `NodeId`. Call `execute` before reading output images.)doc")
        .def(py::init<std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>>, FloatArray>(),
            "tree"_a,
            "attr"_a,
            "Create a UAO computation from a weighted tree and increasing attribute buffer.")
        .def("execute", py::overload_cast<int>(&UltimateAttributeOpeningPybind::execute),
            "maxCriterion"_a,
            "Run UAO using all nodes as selectable candidates.")
        .def("executeWithMSER", &UltimateAttributeOpeningPybind::executeWithMSER,
            "maxCriterion"_a,
            "deltaMSER"_a,
            "Run UAO using an MSER-derived node-selection mask.")
        .def("getMaxContrastImage", &UltimateAttributeOpeningPybind::getMaxContrastImage,
            "Return the 2D uint8 maximum-contrast image.")
        .def("getAssociatedImage", &UltimateAttributeOpeningPybind::getAssociatedImage,
            "Return the 2D int32 associated attribute-index image.")
        .def("getAssociatedColoredImage", &UltimateAttributeOpeningPybind::getAssociatedColorImage,
            "Return a uint8 color rendering of the associated-index image.");
}

} // namespace mmcfilters::pybindings
