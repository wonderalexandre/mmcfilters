#include "ModuleBindings.hpp"

#include "AttributeComputationBindings.hpp"
#include "PythonValuedMorphologicalTree.hpp"

#include "../mmcfilters/attributes/AttributeRegistry.hpp"
#include "../mmcfilters/trees/ValuedMorphologicalTree.hpp"
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
    py::enum_<NodeAttributeSamplingPolicy>(m, "NodeAttributeSamplingPolicy", py::module_local(false))
        .value("LARGEST_SUPPORT_DESCENDANT", NodeAttributeSamplingPolicy::LargestSupportDescendant)
        .export_values();

    py::enum_<MissingNodeAttributeSamplePolicy>(m, "MissingNodeAttributeSamplePolicy", py::module_local(false))
        .value("REPEAT_NEAREST", MissingNodeAttributeSamplePolicy::RepeatNearest)
        .value("NOT_A_NUMBER", MissingNodeAttributeSamplePolicy::NotANumber)
        .value("ZERO", MissingNodeAttributeSamplePolicy::Zero)
        .export_values();

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
                "compute_attributes",
                [](std::shared_ptr<PythonValuedMorphologicalTree> tree, const py::object& attributes, NodeIdSpace outputSpace, py::object dtype) {
                    return attribute_computation::computeAttributesFromList(std::move(tree), attribute_computation::resolveAttributeOrGroupList(attributes),
                                                                            outputSpace, std::move(dtype));
                },
                py::arg("tree"), py::arg("attributes"), py::arg("output_space") = NodeIdSpace::MorphologicalTree, py::kw_only(), py::arg("dtype") = py::none(),
                R"doc(Compute one or more altitude-dependent attributes.

Parameters:
    tree: `ValuedMorphologicalTree`.
    attributes: Sequence of `Attribute` or `Attribute.Group` values.
    output_space: Node-id domain of the returned rows.
    dtype: Returned attribute value dtype, either `np.float32` or `np.float64`.

Returns:
    `(layout, values)` where `values.shape == (num_output_nodes, len(layout))`
    and `values.dtype` matches `dtype`.)doc")
            .def_static(
                "compute_topology_attributes",
                [](std::shared_ptr<PythonValuedMorphologicalTree> tree, const py::object& attributes, NodeIdSpace outputSpace, py::object dtype) {
                    return attribute_computation::computeTopologyAttributesFromList(
                        std::move(tree), attribute_computation::resolveAttributeOrGroupList(attributes), outputSpace, std::move(dtype));
                },
                py::arg("tree"), py::arg("attributes"), py::arg("output_space") = NodeIdSpace::MorphologicalTree, py::kw_only(), py::arg("dtype") = py::none(),
                R"doc(Compute topology/support-only attributes from a valued tree.

The valued tree is accepted for convenience, but altitude-dependent
attributes must use `compute_attributes` or `compute_single_attribute`.)doc")
            .def_static("compute_single_attribute",
                        [](std::shared_ptr<PythonValuedMorphologicalTree> tree, const py::object& attribute, NodeIdSpace outputSpace, py::object dtype) {
                            return attribute_computation::computeSingleAttribute(std::move(tree), attribute_computation::resolveAttribute(attribute),
                                                                                 outputSpace, std::move(dtype));
                        },
                        py::arg("tree"), py::arg("attribute"), py::arg("output_space") = NodeIdSpace::MorphologicalTree, py::kw_only(),
                        py::arg("dtype") = py::none(),
                        R"doc(Compute one altitude-dependent attribute.

Returns a 1D array indexed by `output_space`. `dtype` controls the returned
buffer dtype and accepts `np.float32` or `np.float64`; the default is
`np.float32`.)doc")
            .def_static("compute_single_topology_attribute",
                        [](std::shared_ptr<PythonValuedMorphologicalTree> tree, const py::object& attribute, NodeIdSpace outputSpace, py::object dtype) {
                            return attribute_computation::computeSingleTopologyAttribute(std::move(tree), attribute_computation::resolveAttribute(attribute),
                                                                                        outputSpace, std::move(dtype));
                        },
                        py::arg("tree"), py::arg("attribute"), py::arg("output_space") = NodeIdSpace::MorphologicalTree, py::kw_only(),
                        py::arg("dtype") = py::none(),
                        R"doc(Compute one topology/support-only attribute from a valued tree.

Returns a 1D floating-point array indexed by `output_space`.)doc")
            .def_static("compute_sampled_node_attribute",
                        [](std::shared_ptr<PythonValuedMorphologicalTree> tree, const py::object& attribute, std::int64_t altitudeStep, int samplingRadius,
                           NodeAttributeSamplingPolicy samplingPolicy, MissingNodeAttributeSamplePolicy missingSamplePolicy, NodeIdSpace outputSpace,
                           py::object dtype) {
                            return attribute_computation::computeSampledNodeAttribute(std::move(tree), attribute_computation::resolveAttribute(attribute),
                                                                                      altitudeStep, samplingRadius, samplingPolicy, missingSamplePolicy,
                                                                                      outputSpace, std::move(dtype));
                        },
                        py::arg("tree"), py::arg("attribute"), py::arg("altitude_step"), py::arg("sampling_radius"),
                        py::arg("sampling_policy") = NodeAttributeSamplingPolicy::LargestSupportDescendant,
                        py::arg("missing_sample_policy") = MissingNodeAttributeSamplePolicy::RepeatNearest,
                        py::arg("output_space") = NodeIdSpace::MorphologicalTree, py::kw_only(), py::arg("dtype") = py::none(),
                        R"doc(Sample one node attribute at altitude-based ancestor/current/descendant positions.

Parameters:
    altitude_step: Positive altitude distance between adjacent sample positions.
    sampling_radius: Non-negative radius. Columns cover every integer sample
        offset in `[-sampling_radius, sampling_radius]`.
    sampling_policy: Representative-descendant selection policy.
    missing_sample_policy: Policy for unavailable ancestor or descendant samples.
    dtype: Returned attribute value dtype, either `np.float32` or `np.float64`.

Returns:
    `(layout, values)` where layout keys include suffixes such as
    `_ANCESTOR_1` and `_DESCENDANT_1`.)doc")
            .def_static(
                "describe", [](const py::object& attribute) { return attribute_computation::describeAttribute(attribute_computation::resolveAttribute(attribute)); },
                py::arg("attribute"), "Return the human-readable description of an attribute.")
            .def_static(
                "requirements",
                [](const py::object& attributeRequest) {
                    const Attribute attribute = attribute_computation::resolveAttribute(attributeRequest);
                    const auto requirements = attributes::registry::capabilityRequirements(attribute);
                    const char* adjacency = "none";
                    switch (requirements.adjacency) {
                    case attributes::registry::AttributeAdjacencyRequirement::None:
                        adjacency = "none";
                        break;
                    case attributes::registry::AttributeAdjacencyRequirement::UniformOrDirectional:
                        adjacency = "uniform-or-directional";
                        break;
                    case attributes::registry::AttributeAdjacencyRequirement::Uniform:
                        adjacency = "uniform";
                        break;
                    }

                    py::dict result;
                    result["altitude"] = requirements.altitude;
                    result["grid_domain_2d"] = requirements.gridDomain2D;
                    result["adjacency"] = adjacency;
                    result["monotone_altitude_order"] = requirements.monotoneAltitudeOrder;
                    result["altitude_for_directional_adjacency"] = requirements.altitudeForDirectionalAdjacency;
                    result["canonical_4_or_8_adjacency"] = requirements.canonical4Or8Adjacency;
                    return result;
                },
                py::arg("attribute"), "Return the declared capability requirements of an attribute.")
            .def_static("compute_attribute_mapping",
                        [](std::shared_ptr<PythonValuedMorphologicalTree> tree, const py::object& attribute, py::object dtype) {
                            return attribute_computation::computeAttributeMapping(std::move(tree), attribute_computation::resolveAttribute(attribute),
                                                                                  std::move(dtype));
                        },
                        py::arg("tree"), py::arg("attribute"), py::kw_only(), py::arg("dtype") = py::none(),
                        R"doc(Compute an attribute and project it back to the image domain.

Returns a 2D floating-point array with the same shape as the tree image domain.)doc");

    py::enum_<AttributeGroup>(cls, "Group", py::module_local(false))
        .value("ALL", AttributeGroup::All)
        .value("GRAY_LEVEL", AttributeGroup::GrayLevel)
        .value("SHAPE", AttributeGroup::Shape)
        .value("MOMENTS", AttributeGroup::Moments)
        .value("BOUNDARY", AttributeGroup::Boundary)
        .value("TREE_TOPOLOGY", AttributeGroup::TreeTopology)
        .value("DIST_TRANSF", AttributeGroup::DistTransf)
        .value("DIST_TRANSF_EXACT", AttributeGroup::DistTransfExact)
        .export_values();

    py::enum_<Attribute>(cls, "Type", py::module_local(false))
        .value("AREA", Attribute::Area)
        .value("VOLUME", Attribute::Volume)
        .value("RELATIVE_VOLUME", Attribute::RelativeVolume)
        .value("GRAY_LEVEL_HEIGHT", Attribute::GrayLevelHeight)
        .value("MEAN_GRAY_LEVEL", Attribute::MeanGrayLevel)
        .value("GRAY_LEVEL_VARIANCE", Attribute::GrayLevelVariance)
        .value("BOX_WIDTH", Attribute::BoxWidth)
        .value("BOUNDING_BOX_HEIGHT", Attribute::BoundingBoxHeight)
        .value("RECTANGULARITY", Attribute::Rectangularity)
        .value("DIAGONAL_LENGTH", Attribute::DiagonalLength)
        .value("BOX_COLUMN_MIN", Attribute::BoxColumnMin)
        .value("BOX_COLUMN_MAX", Attribute::BoxColumnMax)
        .value("BOX_ROW_MIN", Attribute::BoxRowMin)
        .value("BOX_ROW_MAX", Attribute::BoxRowMax)
        .value("RATIO_WH", Attribute::RatioWh)
        .value("CENTRAL_MOMENT_20", Attribute::CentralMoment20)
        .value("CENTRAL_MOMENT_02", Attribute::CentralMoment02)
        .value("CENTRAL_MOMENT_11", Attribute::CentralMoment11)
        .value("CENTRAL_MOMENT_30", Attribute::CentralMoment30)
        .value("CENTRAL_MOMENT_03", Attribute::CentralMoment03)
        .value("CENTRAL_MOMENT_21", Attribute::CentralMoment21)
        .value("CENTRAL_MOMENT_12", Attribute::CentralMoment12)
        .value("AXIS_ORIENTATION", Attribute::AxisOrientation)
        .value("LENGTH_MAJOR_AXIS", Attribute::LengthMajorAxis)
        .value("LENGTH_MINOR_AXIS", Attribute::LengthMinorAxis)
        .value("ECCENTRICITY", Attribute::Eccentricity)
        .value("CIRCULARITY", Attribute::Circularity)
        .value("COMPACTNESS", Attribute::Compactness)
        .value("INERTIA", Attribute::Inertia)
        .value("HU_MOMENT_1", Attribute::HuMoment1)
        .value("HU_MOMENT_2", Attribute::HuMoment2)
        .value("HU_MOMENT_3", Attribute::HuMoment3)
        .value("HU_MOMENT_4", Attribute::HuMoment4)
        .value("HU_MOMENT_5", Attribute::HuMoment5)
        .value("HU_MOMENT_6", Attribute::HuMoment6)
        .value("HU_MOMENT_7", Attribute::HuMoment7)
        .value("SUBTREE_HEIGHT", Attribute::SubtreeHeight)
        .value("DEPTH_NODE", Attribute::DepthNode)
        .value("IS_LEAF_NODE", Attribute::IsLeafNode)
        .value("IS_ROOT_NODE", Attribute::IsRootNode)
        .value("NUM_CHILDREN_NODE", Attribute::NumChildrenNode)
        .value("NUM_SIBLINGS_NODE", Attribute::NumSiblingsNode)
        .value("NUM_DESCENDANTS_NODE", Attribute::NumDescendantsNode)
        .value("NUM_LEAF_DESCENDANTS_NODE", Attribute::NumLeafDescendantsNode)
        .value("LEAF_RATIO_NODE", Attribute::LeafRatioNode)
        .value("BALANCE_NODE", Attribute::BalanceNode)
        .value("AVG_CHILD_HEIGHT_NODE", Attribute::AvgChildHeightNode)
        .value("BITQUAD_AREA", Attribute::BitquadArea)
        .value("BITQUAD_NUMBER_EULER", Attribute::BitquadNumberEuler)
        .value("BITQUAD_NUMBER_HOLES", Attribute::BitquadNumberHoles)
        .value("BITQUAD_PERIMETER", Attribute::BitquadPerimeter)
        .value("BITQUAD_PERIMETER_CONTINUOUS", Attribute::BitquadPerimeterContinuous)
        .value("BITQUAD_CIRCULARITY", Attribute::BitquadCircularity)
        .value("BITQUAD_PERIMETER_AVERAGE", Attribute::BitquadPerimeterAverage)
        .value("BITQUAD_LENGTH_AVERAGE", Attribute::BitquadLengthAverage)
        .value("BITQUAD_WIDTH_AVERAGE", Attribute::BitquadWidthAverage)
        .value("MAX_DIST", Attribute::MaxDist)
        .value("CONTOUR_PIXELS", Attribute::ContourPixels)
        .value("CONTOUR_PERIMETER", Attribute::ContourPerimeter)
        .value("CONTOUR_SIDE_NORTH", Attribute::ContourSideNorth)
        .value("CONTOUR_SIDE_WEST", Attribute::ContourSideWest)
        .value("CONTOUR_SIDE_EAST", Attribute::ContourSideEast)
        .value("CONTOUR_SIDE_SOUTH", Attribute::ContourSideSouth)
        .value("MAX_DIST_EXACT", Attribute::MaxDistExact)
        .value("DIST_SQUARED_SUM_EXACT", Attribute::DistSquaredSumExact)
        .value("DIST_SQUARED_MEAN_EXACT", Attribute::DistSquaredMeanExact)
        .value("DIST_RMS_EXACT", Attribute::DistRmsExact)
        .value("DIST_SQUARED_VARIANCE_EXACT", Attribute::DistSquaredVarianceExact)
        .value("DIST_SQUARED_SUM", Attribute::DistSquaredSum)
        .value("DIST_SQUARED_MEAN", Attribute::DistSquaredMean)
        .value("DIST_RMS", Attribute::DistRms)
        .value("DIST_SQUARED_VARIANCE", Attribute::DistSquaredVariance)
        .value("MAX_DIST_CENTER_ROW_EXACT", Attribute::MaxDistCenterRowExact)
        .value("MAX_DIST_CENTER_COLUMN_EXACT", Attribute::MaxDistCenterColumnExact)
        .value("MAX_DIST_CENTER_ROW", Attribute::MaxDistCenterRow)
        .value("MAX_DIST_CENTER_COLUMN", Attribute::MaxDistCenterColumn)
        .value("MAX_DIST_PLATEAU_AREA_EXACT", Attribute::MaxDistPlateauAreaExact)
        .value("MAX_DIST_PLATEAU_CENTROID_ROW_EXACT", Attribute::MaxDistPlateauCentroidRowExact)
        .value("MAX_DIST_PLATEAU_CENTROID_COLUMN_EXACT", Attribute::MaxDistPlateauCentroidColumnExact)
        .value("MAX_DIST_PLATEAU_AREA", Attribute::MaxDistPlateauArea)
        .value("MAX_DIST_PLATEAU_CENTROID_ROW", Attribute::MaxDistPlateauCentroidRow)
        .value("MAX_DIST_PLATEAU_CENTROID_COLUMN", Attribute::MaxDistPlateauCentroidColumn)
        .value("DIST_SUM", Attribute::DistSum)
        .value("DIST_MEAN", Attribute::DistMean)
        .value("DIST_VARIANCE", Attribute::DistVariance)
        .value("DIST_MEDIAN", Attribute::DistMedian)
        .value("DIST_MODE", Attribute::DistMode)
        .value("DIST_Q25", Attribute::DistQ25)
        .value("DIST_Q75", Attribute::DistQ75)
        .value("DIST_Q90", Attribute::DistQ90)
        .value("DIST_ENTROPY", Attribute::DistEntropy)
        .value("DIST_POSITIVE_AREA", Attribute::DistPositiveArea)
        .value("DIST_LEVEL_COUNT", Attribute::DistLevelCount)
        .value("DIST_WEIGHTED_CENTROID_ROW", Attribute::DistWeightedCentroidRow)
        .value("DIST_WEIGHTED_CENTROID_COLUMN", Attribute::DistWeightedCentroidColumn)
        .value("DIST_WEIGHTED_CENTRAL_MOMENT_20", Attribute::DistWeightedCentralMoment20)
        .value("DIST_WEIGHTED_CENTRAL_MOMENT_02", Attribute::DistWeightedCentralMoment02)
        .value("DIST_WEIGHTED_CENTRAL_MOMENT_11", Attribute::DistWeightedCentralMoment11)
        .value("DIST_WEIGHTED_AXIS_ORIENTATION", Attribute::DistWeightedAxisOrientation)
        .value("DIST_WEIGHTED_ECCENTRICITY", Attribute::DistWeightedEccentricity)
        .value("DIST_SUM_EXACT", Attribute::DistSumExact)
        .value("DIST_MEAN_EXACT", Attribute::DistMeanExact)
        .value("DIST_VARIANCE_EXACT", Attribute::DistVarianceExact)
        .value("DIST_MEDIAN_EXACT", Attribute::DistMedianExact)
        .value("DIST_MODE_EXACT", Attribute::DistModeExact)
        .value("DIST_Q25_EXACT", Attribute::DistQ25Exact)
        .value("DIST_Q75_EXACT", Attribute::DistQ75Exact)
        .value("DIST_Q90_EXACT", Attribute::DistQ90Exact)
        .value("DIST_ENTROPY_EXACT", Attribute::DistEntropyExact)
        .value("DIST_POSITIVE_AREA_EXACT", Attribute::DistPositiveAreaExact)
        .value("DIST_LEVEL_COUNT_EXACT", Attribute::DistLevelCountExact)
        .value("DIST_WEIGHTED_CENTROID_ROW_EXACT", Attribute::DistWeightedCentroidRowExact)
        .value("DIST_WEIGHTED_CENTROID_COLUMN_EXACT", Attribute::DistWeightedCentroidColumnExact)
        .value("DIST_WEIGHTED_CENTRAL_MOMENT_20_EXACT", Attribute::DistWeightedCentralMoment20Exact)
        .value("DIST_WEIGHTED_CENTRAL_MOMENT_02_EXACT", Attribute::DistWeightedCentralMoment02Exact)
        .value("DIST_WEIGHTED_CENTRAL_MOMENT_11_EXACT", Attribute::DistWeightedCentralMoment11Exact)
        .value("DIST_WEIGHTED_AXIS_ORIENTATION_EXACT", Attribute::DistWeightedAxisOrientationExact)
        .value("DIST_WEIGHTED_ECCENTRICITY_EXACT", Attribute::DistWeightedEccentricityExact)
        .value("MAX_SQUARED_DIST", Attribute::MaxSquaredDist)
        .value("MAX_SQUARED_DIST_EXACT", Attribute::MaxSquaredDistExact)
        .export_values();
}

} // namespace mmcfilters::pybindings
