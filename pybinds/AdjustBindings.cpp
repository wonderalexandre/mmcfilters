#include "ModuleBindings.hpp"

#include "PybindConversions.hpp"

#include "../mmcfilters/trees/WeightedMorphologicalTree.hpp"
#include "../mmcfilters/trees/adjust/CasfComponentTrees.hpp"
#include "../mmcfilters/trees/adjust/DualMinMaxTreeIncrementalFilter.hpp"
#include "../mmcfilters/utils/RegularGridAdjacency2D.hpp"
#include "../mmcfilters/utils/Common.hpp"
#include "../mmcfilters/utils/Image.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcfilters::pybindings {

namespace py = pybind11;
using namespace pybind11::literals;

using UInt8InputArray = py::array;

namespace {

ImageUInt8Ptr imageFromArray(const UInt8InputArray& input) {
    if (!input.dtype().is(py::dtype::of<uint8_t>())) {
        throw std::invalid_argument("input must be a 2D uint8 array");
    }
    auto buf = input.request();
    if (buf.ndim != 2) {
        throw std::invalid_argument("input must be a 2D uint8 array");
    }
    const int rows = static_cast<int>(buf.shape[0]);
    const int cols = static_cast<int>(buf.shape[1]);
    if (buf.strides[1] != static_cast<py::ssize_t>(sizeof(uint8_t)) || buf.strides[0] != static_cast<py::ssize_t>(cols * sizeof(uint8_t))) {
        throw std::invalid_argument("input must be a C-contiguous 2D uint8 array");
    }
    return ImageUInt8::fromExternal(static_cast<uint8_t*>(buf.ptr), rows, cols);
}

class DualMinMaxTreeIncrementalFilterPybind {
  private:
    std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> minTree_;
    std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> maxTree_;
    RegularGridAdjacency2D adjacency_;
    mmcfilters::adjust::DualMinMaxTreeIncrementalFilter<std::uint8_t> adjust_;

    static RegularGridAdjacency2D makeAdjacency(const WeightedMorphologicalTree<std::uint8_t>& minTree,
                                                const WeightedMorphologicalTree<std::uint8_t>& maxTree) {
        const RegularGridAdjacency2D* minAdjacency = minTree.topology().getUniformGridAdjacency2D();
        const RegularGridAdjacency2D* maxAdjacency = maxTree.topology().getUniformGridAdjacency2D();
        if (minAdjacency == nullptr || maxAdjacency == nullptr) {
            throw std::invalid_argument("DualMinMaxTreeIncrementalFilter requires weighted component trees with adjacency information.");
        }
        if (minTree.topology().getNumRowsOfGridDomain2D() != maxTree.topology().getNumRowsOfGridDomain2D() ||
            minTree.topology().getNumColsOfGridDomain2D() != maxTree.topology().getNumColsOfGridDomain2D() ||
            minTree.topology().getNumTotalProperParts() != maxTree.topology().getNumTotalProperParts()) {
            throw std::invalid_argument("DualMinMaxTreeIncrementalFilter requires trees built on the same image domain.");
        }
        return *minAdjacency;
    }

  public:
    DualMinMaxTreeIncrementalFilterPybind(std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> minTree,
                                          std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> maxTree)
        : minTree_(std::move(minTree)), maxTree_(std::move(maxTree)), adjacency_(makeAdjacency(*minTree_, *maxTree_)),
          adjust_(minTree_.get(), maxTree_.get(), adjacency_) {}

    void pruneMaxTreeAndUpdateMinTree(const std::vector<NodeId>& nodesToPrune) { adjust_.pruneMaxTreeAndUpdateMinTree(nodesToPrune); }

    void pruneMinTreeAndUpdateMaxTree(const std::vector<NodeId>& nodesToPrune) { adjust_.pruneMinTreeAndUpdateMaxTree(nodesToPrune); }

    std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> minTree() const { return minTree_; }

    std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>> maxTree() const { return maxTree_; }
};

} // namespace

/**
 * @brief Registers component tree adjust bindings in the Python module.
 *
 * @param m Python module receiving the bindings.
 */
void initComponentTreeAdjust(py::module_& m) {
    namespace adjust = mmcfilters::adjust;

    py::enum_<adjust::CasfComponentTreesAttribute>(m, "CasfComponentTreesAttribute", py::module_local(false))
        .value("AREA", adjust::CasfComponentTreesAttribute::AREA)
        .value("BOUNDING_BOX_WIDTH", adjust::CasfComponentTreesAttribute::BOUNDING_BOX_WIDTH)
        .value("BOUNDING_BOX_HEIGHT", adjust::CasfComponentTreesAttribute::BOUNDING_BOX_HEIGHT)
        .value("BOUNDING_BOX_DIAGONAL", adjust::CasfComponentTreesAttribute::BOUNDING_BOX_DIAGONAL)
        .export_values();

    py::class_<DualMinMaxTreeIncrementalFilterPybind, std::shared_ptr<DualMinMaxTreeIncrementalFilterPybind>>(
        m, "DualMinMaxTreeIncrementalFilter", py::module_local(false),
        R"doc(Incremental updater for paired weighted min-tree/max-tree component trees.

Both trees must be weighted component trees built on the same 2D uint8 image
domain and must carry adjacency metadata.)doc")
        .def(py::init<std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>>, std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>>>(), "minTree"_a,
             "maxTree"_a, "Create the paired incremental updater.")
        .def("pruneMaxTreeAndUpdateMinTree", &DualMinMaxTreeIncrementalFilterPybind::pruneMaxTreeAndUpdateMinTree, "nodesToPrune"_a,
             "Update the min-tree incrementally and then prune the requested max-tree subtrees.")
        .def("pruneMinTreeAndUpdateMaxTree", &DualMinMaxTreeIncrementalFilterPybind::pruneMinTreeAndUpdateMaxTree, "nodesToPrune"_a,
             "Update the max-tree incrementally and then prune the requested min-tree subtrees.")
        .def_property_readonly("minTree", &DualMinMaxTreeIncrementalFilterPybind::minTree, "Borrowed min-tree owner updated by paired pruning operations.")
        .def_property_readonly("maxTree", &DualMinMaxTreeIncrementalFilterPybind::maxTree, "Borrowed max-tree owner updated by paired pruning operations.");

    py::class_<adjust::CasfComponentTrees<std::uint8_t>, std::shared_ptr<adjust::CasfComponentTrees<std::uint8_t>>>(
        m, "CasfComponentTrees", py::module_local(false),
        R"doc(Connected alternating sequential filter using paired component trees.

The input image must be a 2D C-contiguous `np.uint8` array. Thresholds passed to
`filter` are interpreted in increasing order for the selected CASF attribute.)doc")
        .def(py::init([](UInt8InputArray input, adjust::CasfComponentTreesAttribute attribute, double radius) {
                 return std::make_shared<adjust::CasfComponentTrees<std::uint8_t>>(imageFromArray(input), attribute, radius);
             }),
             "input"_a, "attribute"_a = adjust::CasfComponentTreesAttribute::AREA, "radius"_a = 1.5,
             "Create a CASF object and its paired min/max component trees.")
        .def(
            "filter",
            [](adjust::CasfComponentTrees<std::uint8_t>& self, const std::vector<double>& thresholds) {
                return pybind_utils::toNumpy(self.filter(thresholds));
            },
            "thresholds"_a, "Run the connected alternating sequential filter and return a 2D uint8 image.")
        .def_property_readonly(
            "minTree", [](adjust::CasfComponentTrees<std::uint8_t>& self) -> const WeightedMorphologicalTree<std::uint8_t>& { return self.minTree(); },
            py::return_value_policy::reference_internal, "Internal min-tree maintained by the CASF.")
        .def_property_readonly(
            "maxTree", [](adjust::CasfComponentTrees<std::uint8_t>& self) -> const WeightedMorphologicalTree<std::uint8_t>& { return self.maxTree(); },
            py::return_value_policy::reference_internal, "Internal max-tree maintained by the CASF.")
        .def_property_readonly("attribute", &adjust::CasfComponentTrees<std::uint8_t>::attribute, "CASF attribute used to evaluate thresholds.")
        .def("exportMinTree", &adjust::CasfComponentTrees<std::uint8_t>::exportMinTree,
             "Export the current min-tree as a compact Higra `(parent, altitude)` pair.")
        .def("exportMaxTree", &adjust::CasfComponentTrees<std::uint8_t>::exportMaxTree,
             "Export the current max-tree as a compact Higra `(parent, altitude)` pair.");
}

} // namespace mmcfilters::pybindings
