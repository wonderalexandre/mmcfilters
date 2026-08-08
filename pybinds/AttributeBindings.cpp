#include "ModuleBindings.hpp"

#include "AttributeComputationBindings.hpp"

#include "../mmcfilters/attributes/AttributeRegistry.hpp"
#include "../mmcfilters/trees/WeightedMorphologicalTree.hpp"
#include "../mmcfilters/utils/Common.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace mmcfilters::pybindings {

namespace py = pybind11;

namespace {

struct AttributeBinding {};

} // namespace

/**
 * @brief Registers attribute computation bindings in the Python module.
 *
 * @param m Python module receiving the bindings.
 */
void initAttributeComputation(py::module_& m) {
    auto cls =
        py::class_<AttributeBinding>(m, "Attribute", py::module_local(false),
                                     R"doc(Attribute computation utilities based on dense node-id buffers.

Methods returning a single attribute return a 1D floating-point array indexed by
the selected `NodeIdSpace`. Methods returning several attributes return
`(layout, values)`, where `layout` maps stable attribute names to columns and
`values` is a 2D floating-point array with one row per output node. The optional
`dtype` keyword accepts `np.float32` or `np.float64` and defaults to
`np.float32`; it controls only the returned NumPy buffer dtype.)doc")
            .def_static(
                "computeAttributes",
                py::overload_cast<std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>>, const std::vector<AttributeOrGroup>&, NodeIdSpace, py::object>(
                    &attribute_computation::computeAttributesFromList),
                py::arg("tree"), py::arg("attributes"), py::arg("outputSpace") = NodeIdSpace::MORPHOLOGICAL_TREE, py::kw_only(), py::arg("dtype") = py::none(),
                R"doc(Compute one or more altitude-dependent attributes.

Parameters:
    tree: `WeightedMorphologicalTree` with uint8 altitudes.
    attributes: Sequence of `Attribute` or `Attribute.Group` values.
    outputSpace: Node-id domain of the returned rows.
    dtype: Returned attribute value dtype, either `np.float32` or `np.float64`.

Returns:
    `(layout, values)` where `values.shape == (num_output_nodes, len(layout))`
    and `values.dtype` matches `dtype`.)doc")
            .def_static(
                "computeTopologyAttributes",
                py::overload_cast<std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>>, const std::vector<AttributeOrGroup>&, NodeIdSpace, py::object>(
                    &attribute_computation::computeTopologyAttributesFromList),
                py::arg("tree"), py::arg("attributes"), py::arg("outputSpace") = NodeIdSpace::MORPHOLOGICAL_TREE, py::kw_only(), py::arg("dtype") = py::none(),
                R"doc(Compute topology/support-only attributes from a weighted tree.

The weighted tree is accepted for convenience, but altitude-dependent
attributes must use `computeAttributes` or `computeSingleAttribute`.)doc")
            .def_static("computeSingleAttribute",
                        py::overload_cast<std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>>, Attribute, NodeIdSpace, py::object>(
                            &attribute_computation::computeSingleAttribute),
                        py::arg("tree"), py::arg("attribute"), py::arg("outputSpace") = NodeIdSpace::MORPHOLOGICAL_TREE, py::kw_only(),
                        py::arg("dtype") = py::none(),
                        R"doc(Compute one altitude-dependent attribute.

Returns a 1D array indexed by `outputSpace`. `dtype` controls the returned
buffer dtype and accepts `np.float32` or `np.float64`; the default is
`np.float32`.)doc")
            .def_static("computeSingleTopologyAttribute",
                        py::overload_cast<std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>>, Attribute, NodeIdSpace, py::object>(
                            &attribute_computation::computeSingleTopologyAttribute),
                        py::arg("tree"), py::arg("attribute"), py::arg("outputSpace") = NodeIdSpace::MORPHOLOGICAL_TREE, py::kw_only(),
                        py::arg("dtype") = py::none(),
                        R"doc(Compute one topology/support-only attribute from a weighted tree.

Returns a 1D floating-point array indexed by `outputSpace`.)doc")
            .def_static("computeSingleAttributeWithDelta",
                        py::overload_cast<std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>>, Attribute, int, std::string, NodeIdSpace, py::object>(
                            &attribute_computation::computeSingleAttributeWithDelta),
                        py::arg("tree"), py::arg("attribute"), py::arg("delta"), py::arg("padding") = "last-padding",
                        py::arg("outputSpace") = NodeIdSpace::MORPHOLOGICAL_TREE, py::kw_only(), py::arg("dtype") = py::none(),
                        R"doc(Compute one attribute at ancestor/descendant delta offsets.

Parameters:
    delta: Non-negative neighbourhood radius. Columns cover `[-delta, delta]`.
    padding: Boundary strategy. Supported values are `"last-padding"`,
        `"nan-padding"`, `"null-padding"`, and `"zero-padding"`.
    dtype: Returned attribute value dtype, either `np.float32` or `np.float64`.

Returns:
    `(layout, values)` where layout keys include suffixes such as `_ASC_1`
    and `_DESC_1`.)doc")
            .def_static("describe", &attribute_computation::describeAttribute, py::arg("attribute"), "Return the human-readable description of an attribute.")
            .def_static(
                "requirements",
                [](Attribute attribute) {
                    const auto requirements = attributes::registry::capabilityRequirements(attribute);
                    const char* adjacency = "none";
                    switch (requirements.adjacency) {
                    case attributes::registry::AttributeAdjacencyRequirement::NONE:
                        adjacency = "none";
                        break;
                    case attributes::registry::AttributeAdjacencyRequirement::UNIFORM_OR_DIRECTIONAL:
                        adjacency = "uniform-or-directional";
                        break;
                    case attributes::registry::AttributeAdjacencyRequirement::UNIFORM:
                        adjacency = "uniform";
                        break;
                    }

                    py::dict result;
                    result["altitude"] = requirements.altitude;
                    result["gridDomain2D"] = requirements.gridDomain2D;
                    result["adjacency"] = adjacency;
                    result["monotoneAltitudeOrder"] = requirements.monotoneAltitudeOrder;
                    result["altitudeForDirectionalAdjacency"] = requirements.altitudeForDirectionalAdjacency;
                    result["canonical4Or8Adjacency"] = requirements.canonical4Or8Adjacency;
                    return result;
                },
                py::arg("attribute"), "Return the declared capability requirements of an attribute.")
            .def_static("computeAttributeMapping",
                        py::overload_cast<std::shared_ptr<WeightedMorphologicalTree<std::uint8_t>>, Attribute, py::object>(
                            &attribute_computation::computeAttributeMapping),
                        py::arg("tree"), py::arg("attribute"), py::kw_only(), py::arg("dtype") = py::none(),
                        R"doc(Compute an attribute and project it back to the image domain.

Returns a 2D floating-point array with the same shape as the tree image domain.)doc");

    py::enum_<AttributeGroup>(cls, "Group", py::module_local(false))
        .value("ALL", AttributeGroup::ALL)
        .value("GRAY_LEVEL", AttributeGroup::GRAY_LEVEL)
        .value("SHAPE", AttributeGroup::SHAPE)
        .value("MOMENTS", AttributeGroup::MOMENTS)
        .value("BOUNDARY", AttributeGroup::BOUNDARY)
        .value("TREE_TOPOLOGY", AttributeGroup::TREE_TOPOLOGY)
        .export_values();

    py::enum_<Attribute>(cls, "Type", py::module_local(false))
        .value("AREA", Attribute::AREA)
        .value("VOLUME", Attribute::VOLUME)
        .value("RELATIVE_VOLUME", Attribute::RELATIVE_VOLUME)
        .value("LEVEL", Attribute::LEVEL)
        .value("GRAY_HEIGHT", Attribute::GRAY_HEIGHT)
        .value("MEAN_LEVEL", Attribute::MEAN_LEVEL)
        .value("VARIANCE_LEVEL", Attribute::VARIANCE_LEVEL)
        .value("BOX_WIDTH", Attribute::BOX_WIDTH)
        .value("BOX_HEIGHT", Attribute::BOX_HEIGHT)
        .value("RECTANGULARITY", Attribute::RECTANGULARITY)
        .value("DIAGONAL_LENGTH", Attribute::DIAGONAL_LENGTH)
        .value("BOX_COL_MIN", Attribute::BOX_COL_MIN)
        .value("BOX_COL_MAX", Attribute::BOX_COL_MAX)
        .value("BOX_ROW_MIN", Attribute::BOX_ROW_MIN)
        .value("BOX_ROW_MAX", Attribute::BOX_ROW_MAX)
        .value("RATIO_WH", Attribute::RATIO_WH)
        .value("CENTRAL_MOMENT_20", Attribute::CENTRAL_MOMENT_20)
        .value("CENTRAL_MOMENT_02", Attribute::CENTRAL_MOMENT_02)
        .value("CENTRAL_MOMENT_11", Attribute::CENTRAL_MOMENT_11)
        .value("CENTRAL_MOMENT_30", Attribute::CENTRAL_MOMENT_30)
        .value("CENTRAL_MOMENT_03", Attribute::CENTRAL_MOMENT_03)
        .value("CENTRAL_MOMENT_21", Attribute::CENTRAL_MOMENT_21)
        .value("CENTRAL_MOMENT_12", Attribute::CENTRAL_MOMENT_12)
        .value("AXIS_ORIENTATION", Attribute::AXIS_ORIENTATION)
        .value("LENGTH_MAJOR_AXIS", Attribute::LENGTH_MAJOR_AXIS)
        .value("LENGTH_MINOR_AXIS", Attribute::LENGTH_MINOR_AXIS)
        .value("ECCENTRICITY", Attribute::ECCENTRICITY)
        .value("CIRCULARITY", Attribute::CIRCULARITY)
        .value("COMPACTNESS", Attribute::COMPACTNESS)
        .value("INERTIA", Attribute::INERTIA)
        .value("HU_MOMENT_1", Attribute::HU_MOMENT_1)
        .value("HU_MOMENT_2", Attribute::HU_MOMENT_2)
        .value("HU_MOMENT_3", Attribute::HU_MOMENT_3)
        .value("HU_MOMENT_4", Attribute::HU_MOMENT_4)
        .value("HU_MOMENT_5", Attribute::HU_MOMENT_5)
        .value("HU_MOMENT_6", Attribute::HU_MOMENT_6)
        .value("HU_MOMENT_7", Attribute::HU_MOMENT_7)
        .value("HEIGHT_NODE", Attribute::HEIGHT_NODE)
        .value("DEPTH_NODE", Attribute::DEPTH_NODE)
        .value("IS_LEAF_NODE", Attribute::IS_LEAF_NODE)
        .value("IS_ROOT_NODE", Attribute::IS_ROOT_NODE)
        .value("NUM_CHILDREN_NODE", Attribute::NUM_CHILDREN_NODE)
        .value("NUM_SIBLINGS_NODE", Attribute::NUM_SIBLINGS_NODE)
        .value("NUM_DESCENDANTS_NODE", Attribute::NUM_DESCENDANTS_NODE)
        .value("NUM_LEAF_DESCENDANTS_NODE", Attribute::NUM_LEAF_DESCENDANTS_NODE)
        .value("LEAF_RATIO_NODE", Attribute::LEAF_RATIO_NODE)
        .value("BALANCE_NODE", Attribute::BALANCE_NODE)
        .value("AVG_CHILD_HEIGHT_NODE", Attribute::AVG_CHILD_HEIGHT_NODE)
        .value("BITQUADS_AREA", Attribute::BITQUADS_AREA)
        .value("BITQUADS_NUMBER_EULER", Attribute::BITQUADS_NUMBER_EULER)
        .value("BITQUADS_NUMBER_HOLES", Attribute::BITQUADS_NUMBER_HOLES)
        .value("BITQUADS_PERIMETER", Attribute::BITQUADS_PERIMETER)
        .value("BITQUADS_PERIMETER_CONTINUOUS", Attribute::BITQUADS_PERIMETER_CONTINUOUS)
        .value("BITQUADS_CIRCULARITY", Attribute::BITQUADS_CIRCULARITY)
        .value("BITQUADS_PERIMETER_AVERAGE", Attribute::BITQUADS_PERIMETER_AVERAGE)
        .value("BITQUADS_LENGTH_AVERAGE", Attribute::BITQUADS_LENGTH_AVERAGE)
        .value("BITQUADS_WIDTH_AVERAGE", Attribute::BITQUADS_WIDTH_AVERAGE)
        .value("MAX_DIST", Attribute::MAX_DIST)
        .value("CONTOUR_PIXELS", Attribute::CONTOUR_PIXELS)
        .value("CONTOUR_PERIMETER", Attribute::CONTOUR_PERIMETER)
        .value("CONTOUR_SIDE_NORTH", Attribute::CONTOUR_SIDE_NORTH)
        .value("CONTOUR_SIDE_WEST", Attribute::CONTOUR_SIDE_WEST)
        .value("CONTOUR_SIDE_EAST", Attribute::CONTOUR_SIDE_EAST)
        .value("CONTOUR_SIDE_SOUTH", Attribute::CONTOUR_SIDE_SOUTH)
        .export_values();
}

} // namespace mmcfilters::pybindings
