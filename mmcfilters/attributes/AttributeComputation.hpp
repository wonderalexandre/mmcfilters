#pragma once

#include "../utils/Image.hpp"
#include "../attributes/AttributeResultTypes.hpp"
#include "../trees/MorphologicalTree.hpp"
#include "../trees/WeightedMorphologicalTree.hpp"
#include "../trees/WeightedTreeView.hpp"

#include <concepts>
#include <span>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace mmcfilters {

/**
 * @brief Public facade for attribute computation.
 *
 * @details
 * This class is the public facade of the attribute framework. It provides:
 * - routing to `AttributePipeline` for ordinary weighted/topology requests;
 * - cache-aware materialisation of one attribute, one group, or several
 *   heterogeneous requests at once inside the implementation;
 * - projection from the internal `MorphologicalTree` node-id space to other
 *   public node-id spaces such as the preserved imported Higra convention;
 * - convenience helpers that project node attributes back to the image domain.
 *
 * The surrounding attribute subsystem is organised around:
 * - `AttributeNames` / `AttributeNamesWithDelta` for flat-buffer layout;
 * - typed attribute kernels for concrete attribute families;
 * - `AttributeResultTypes.hpp` for owning public results;
 * - `AttributePipeline` / `TopologyAttributeBackend` for ordinary requests.
 *
 * Ordinary public attribute computation is expressed through
 * `WeightedMorphologicalTree<T>` or `WeightedTreeView<T>`, so altitude-dependent
 * requests always carry an explicit altitude contract. Plain `MorphologicalTree`
 * remains public only for topology-only requests.
 *
 * The canonical execution space is always the tree's dense internal node-id
 * space. Projection to the preserved imported Higra space only happens at the
 * boundary of the public API. Attribute value buffers are materialized as
 * `float` by default; public methods can select any `std::floating_point`
 * output `Real` type, with `double` supported by the packaged Python bindings.
 * The public `Real` controls result storage. Ordinary facade computations run
 * their internal attribute pipeline in `double` and cast only at the API
 * boundary.
 */
class AttributeComputation {
public:

    /**
     * @brief Computes one topology/support-only scalar attribute or group.
     *
     * @param tree Topology whose dense internal `NodeId` domain indexes the result.
     * @param attr Scalar attribute or group that must not require altitude.
     * @param outputSpace Node-id space requested for the returned buffer.
     * @return Owning attribute data and layout projected to `outputSpace`.
     *
     * @throws std::invalid_argument If `attr` requires altitude or `outputSpace`
     * is unavailable for `tree`.
     */
    template <std::floating_point Real = float>
    [[nodiscard]] static ComputedAttributeData<Real> computeSingleTopologyAttribute(const MorphologicalTree& tree, AttributeOrGroup attr, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE);

    /**
     * @brief Computes several topology/support-only attributes or groups.
     *
     * @param tree Topology whose dense internal `NodeId` domain indexes the result.
     * @param attributes Scalar attributes and groups that must not require altitude.
     * @param outputSpace Node-id space requested for the returned buffer.
     * @return Owning node-major attribute data projected to `outputSpace`.
     *
     * @throws std::invalid_argument If any request requires altitude or the output
     * space is unavailable for `tree`.
     */
    template <std::floating_point Real = float>
    [[nodiscard]] static ComputedAttributeData<Real> computeTopologyAttributes(const MorphologicalTree& tree, const std::vector<AttributeOrGroup>& attributes, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE);

    /**
     * @brief Computes one topology/support-only scalar attribute or group on a weighted owner.
     *
     * @tparam T Altitude type owned by `tree`.
     * @param tree Weighted tree whose topology defines the node-id domain.
     * @param attr Scalar attribute or group that must not require altitude.
     * @param outputSpace Node-id space requested for the returned buffer.
     * @return Owning attribute data and layout projected to `outputSpace`.
     */
    template<std::floating_point Real = float, AltitudeValue T>
    [[nodiscard]] static ComputedAttributeData<Real> computeSingleTopologyAttribute(const WeightedMorphologicalTree<T>& tree, AttributeOrGroup attr, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE);

    /**
     * @brief Computes one topology/support-only scalar attribute or group on a weighted view.
     *
     * @tparam T Altitude type borrowed by `tree`.
     * @param tree Weighted view whose topology defines the node-id domain.
     * @param attr Scalar attribute or group that must not require altitude.
     * @param outputSpace Node-id space requested for the returned buffer.
     * @return Owning attribute data and layout projected to `outputSpace`.
     *
     * @throws std::runtime_error If the borrowed topology changed since view construction.
     */
    template<std::floating_point Real = float, AltitudeValue T>
    [[nodiscard]] static ComputedAttributeData<Real> computeSingleTopologyAttribute(const WeightedTreeView<T>& tree, AttributeOrGroup attr, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE);

    /**
     * @brief Computes several topology/support-only attributes or groups on a weighted owner.
     *
     * @tparam T Altitude type owned by `tree`.
     * @param tree Weighted tree whose topology defines the node-id domain.
     * @param attributes Scalar attributes and groups that must not require altitude.
     * @param outputSpace Node-id space requested for the returned buffer.
     * @return Owning node-major attribute data projected to `outputSpace`.
     */
    template<std::floating_point Real = float, AltitudeValue T>
    [[nodiscard]] static ComputedAttributeData<Real> computeTopologyAttributes(const WeightedMorphologicalTree<T>& tree, const std::vector<AttributeOrGroup>& attributes, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE);

    /**
     * @brief Computes several topology/support-only attributes or groups on a weighted view.
     *
     * @tparam T Altitude type borrowed by `tree`.
     * @param tree Weighted view whose topology defines the node-id domain.
     * @param attributes Scalar attributes and groups that must not require altitude.
     * @param outputSpace Node-id space requested for the returned buffer.
     * @return Owning node-major attribute data projected to `outputSpace`.
     *
     * @throws std::runtime_error If the borrowed topology changed since view construction.
     */
    template<std::floating_point Real = float, AltitudeValue T>
    [[nodiscard]] static ComputedAttributeData<Real> computeTopologyAttributes(const WeightedTreeView<T>& tree, const std::vector<AttributeOrGroup>& attributes, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE);

    /**
     * @brief Computes a single scalar attribute or a full attribute group.
     *
     * @details
     * Requests are routed through the typed weighted-tree attribute pipeline.
     *
     * @tparam T Altitude type owned by `tree`.
     * @param tree Weighted tree supplying topology and dense altitude data.
     * @param attr Scalar attribute or group to compute.
     * @param outputSpace Node-id space requested for the returned buffer.
     * @return Owning attribute data and layout projected to `outputSpace`.
     */
    template<std::floating_point Real = float, AltitudeValue T>
    [[nodiscard]] static ComputedAttributeData<Real> computeSingleAttribute(const WeightedMorphologicalTree<T>& tree, AttributeOrGroup attr, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE);

    /**
     * @brief Computes a single scalar attribute or group on a non-owning weighted view.
     *
     * @tparam T Altitude type borrowed by `tree`.
     * @param tree Weighted view supplying topology and dense altitude data.
     * @param attr Scalar attribute or group to compute.
     * @param outputSpace Node-id space requested for the returned buffer.
     * @return Owning attribute data and layout projected to `outputSpace`.
     *
     * @throws std::runtime_error If the borrowed topology changed since view construction.
     */
    template<std::floating_point Real = float, AltitudeValue T>
    [[nodiscard]] static ComputedAttributeData<Real> computeSingleAttribute(const WeightedTreeView<T>& tree, AttributeOrGroup attr, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE);

    /**
     * @brief Computes a delta-augmented version of one scalar attribute.
     *
     * @tparam T Altitude type owned by `tree`.
     * @param tree Weighted tree supplying topology and dense altitude data.
     * @param attribute Scalar attribute to compute before delta sampling.
     * @param deltaStep Typed altitude distance between consecutive materialised
     * offsets.
     * @param radius Maximum signed ancestor/descendant offset to materialise.
     * @param padding Strategy used when the requested offset leaves the valid
     * ancestor/descendant chain.
     * @param outputSpace Node-id space requested for the returned buffer.
     * @return Owning delta-augmented data and layout projected to `outputSpace`.
     */
    template<std::floating_point Real = float, AltitudeValue T>
    [[nodiscard]] static ComputedAttributeDataWithDelta<Real> computeSingleAttributeWithDelta(const WeightedMorphologicalTree<T>& tree, Attribute attribute, AltitudeDiff<T> deltaStep, int radius, std::string padding="last-padding", NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE);

    /**
     * @brief Computes a delta-augmented scalar attribute on a non-owning weighted view.
     *
     * @tparam T Altitude type borrowed by `tree`.
     * @param tree Weighted view supplying topology and dense altitude data.
     * @param attribute Scalar attribute to compute before delta sampling.
     * @param deltaStep Typed altitude distance between consecutive materialised offsets.
     * @param radius Maximum signed ancestor/descendant offset to materialise.
     * @param padding Strategy used when the requested offset leaves the valid chain.
     * @param outputSpace Node-id space requested for the returned buffer.
     * @return Owning delta-augmented data and layout projected to `outputSpace`.
     *
     * @throws std::runtime_error If the borrowed topology changed since view construction.
     */
    template<std::floating_point Real = float, AltitudeValue T>
    [[nodiscard]] static ComputedAttributeDataWithDelta<Real> computeSingleAttributeWithDelta(const WeightedTreeView<T>& tree, Attribute attribute, AltitudeDiff<T> deltaStep, int radius, std::string padding="last-padding", NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE);

    /**
     * @brief Computes a heterogeneous set of scalar attributes and attribute
     * groups in one coordinated run.
     *
     * @details
     * Requests are routed to `AttributePipeline`. That path combines typed altitude kernels with the topology
     * backend and returns a single layout containing exactly the requested
     * scalar attributes.
     *
     * @tparam T Altitude type owned by `tree`.
     * @param tree Weighted tree supplying topology and dense altitude data.
     * @param attributes Scalar attributes and groups to compute.
     * @param outputSpace Node-id space requested for the returned buffer.
     * @return Owning node-major attribute data projected to `outputSpace`.
     */
    template<std::floating_point Real = float, AltitudeValue T>
    [[nodiscard]] static ComputedAttributeData<Real> computeAttributes(const WeightedMorphologicalTree<T>& tree, const std::vector<AttributeOrGroup>& attributes, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE);

    /**
     * @brief Computes several scalar attributes and groups on a non-owning weighted view.
     *
     * @tparam T Altitude type borrowed by `tree`.
     * @param tree Weighted view supplying topology and dense altitude data.
     * @param attributes Scalar attributes and groups to compute.
     * @param outputSpace Node-id space requested for the returned buffer.
     * @return Owning node-major attribute data projected to `outputSpace`.
     *
     * @throws std::runtime_error If the borrowed topology changed since view construction.
     */
    template<std::floating_point Real = float, AltitudeValue T>
    [[nodiscard]] static ComputedAttributeData<Real> computeAttributes(const WeightedTreeView<T>& tree, const std::vector<AttributeOrGroup>& attributes, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE);

    /**
     * @brief Computes attributes over an external altitude span.
     *
     * @details
     * Attributes that read altitude use templated kernels. Attributes that are
     * independent from altitude are delegated to the topology backend, avoiding
     * an altitude copy or conversion. `WeightedTreeView<T>` is the named
     * non-owning form of the same topology plus external-altitude contract.
     *
     * @tparam T Altitude type borrowed by `weighted`.
     * @param weighted Weighted view supplying topology and dense altitude data.
     * @param attributes Scalar attributes and groups to compute.
     * @param outputSpace Node-id space requested for the returned buffer.
     * @return Owning node-major attribute data projected to `outputSpace`.
     */
    template<std::floating_point Real = float, AltitudeValue T>
    [[nodiscard]] static ComputedAttributeData<Real> computeAttributesFromAltitudeSpan(const WeightedTreeView<T>& weighted, const std::vector<AttributeOrGroup>& attributes, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE);

    /**
     * @brief Projects an already computed internal-node attribute buffer to the
     * compact Higra layout produced by `TreeAltitudeAlgorithms::exportHigraHierarchy()`.
     *
     * @details
     * `nodeValues` must be row-major by dense internal `NodeId`, with one column
     * per attribute described by `attrNames`. The exported buffer follows the
     * compact Higra convention `[proper parts | live internal nodes]`. Unit
     * proper-part rows are computed by the export projection path for each
     * requested attribute, and internal-node rows are copied from `nodeValues`
     * using the same layout helper as hierarchy export.
     *
     * @tparam T Altitude type owned by `tree`.
     * @param tree Weighted tree defining the current exported Higra hierarchy.
     * @param attrNames Layout of the input `nodeValues` columns.
     * @param nodeValues Dense internal-node buffer in row-major layout.
     * @return Row-major values in compact exported Higra layout.
     */
    template<std::floating_point Real = float, AltitudeValue T>
    [[nodiscard]] static std::vector<Real> projectNodeValuesToExportedHigra(const WeightedMorphologicalTree<T>& tree, const AttributeNames& attrNames, std::span<const Real> nodeValues);

    /**
     * @brief Projects a vector-backed internal-node attribute buffer to exported Higra layout.
     *
     * @tparam T Altitude type owned by `tree`.
     * @param tree Weighted tree defining the current exported Higra hierarchy.
     * @param attrNames Layout of the input `nodeValues` columns.
     * @param nodeValues Dense internal-node buffer in row-major layout.
     * @return Row-major values in compact exported Higra layout.
     */
    template<std::floating_point Real = float, AltitudeValue T>
    [[nodiscard]] static std::vector<Real> projectNodeValuesToExportedHigra(const WeightedMorphologicalTree<T>& tree, const AttributeNames& attrNames, const std::vector<Real>& nodeValues);

    /**
     * @brief Projects an internal-node attribute buffer from a weighted view to exported Higra layout.
     *
     * @tparam T Altitude type borrowed by `tree`.
     * @param tree Weighted view defining topology and unit-component altitude data.
     * @param attrNames Layout of the input `nodeValues` columns.
     * @param nodeValues Dense internal-node buffer in row-major layout.
     * @return Row-major values in compact exported Higra layout.
     */
    template<std::floating_point Real = float, AltitudeValue T>
    [[nodiscard]] static std::vector<Real> projectNodeValuesToExportedHigra(const WeightedTreeView<T>& tree, const AttributeNames& attrNames, std::span<const Real> nodeValues);

    /**
     * @brief Projects a vector-backed internal-node attribute buffer from a weighted view to exported Higra layout.
     *
     * @tparam T Altitude type borrowed by `tree`.
     * @param tree Weighted view defining topology and unit-component altitude data.
     * @param attrNames Layout of the input `nodeValues` columns.
     * @param nodeValues Dense internal-node buffer in row-major layout.
     * @return Row-major values in compact exported Higra layout.
     */
    template<std::floating_point Real = float, AltitudeValue T>
    [[nodiscard]] static std::vector<Real> projectNodeValuesToExportedHigra(const WeightedTreeView<T>& tree, const AttributeNames& attrNames, const std::vector<Real>& nodeValues);

    /**
     * @brief Projects a node attribute to a proper-part image in the original domain.
     *
     * @tparam T Altitude type owned by `tree`.
     * @param tree Weighted tree supplying topology and dense altitude data.
     * @param attribute Scalar attribute to compute and map.
     * @return Image where each proper part receives its owner-node attribute value.
     */
    template<std::floating_point Real = float, AltitudeValue T>
    [[nodiscard]] static ImagePtr<Real> computeAttributeMapping(const WeightedMorphologicalTree<T>& tree, Attribute attribute);

    /**
     * @brief Projects a node attribute from a weighted view to a proper-part image.
     *
     * @tparam T Altitude type borrowed by `tree`.
     * @param tree Weighted view supplying topology and dense altitude data.
     * @param attribute Scalar attribute to compute and map.
     * @return Image where each proper part receives its owner-node attribute value.
     */
    template<std::floating_point Real = float, AltitudeValue T>
    [[nodiscard]] static ImagePtr<Real> computeAttributeMapping(const WeightedTreeView<T>& tree, Attribute attribute);

};

} // namespace mmcfilters



#include "../attributes/detail/AttributeDeltaMaterialization.hpp"
#include "../attributes/detail/AttributePipeline.hpp"
#include "../attributes/detail/AttributeProjection.hpp"



namespace mmcfilters {

namespace detail {

template<std::floating_point OutputReal, std::floating_point InternalReal>
[[nodiscard]] inline ComputedAttributeData<OutputReal> castComputedAttributeData(ComputedAttributeData<InternalReal> computed) {
    if constexpr (std::is_same_v<OutputReal, InternalReal>) {
        return computed;
    } else {
        std::vector<OutputReal> output;
        output.reserve(computed.second.size());
        for (const InternalReal value : computed.second) {
            output.push_back(static_cast<OutputReal>(value));
        }
        return ComputedAttributeData<OutputReal>{
            std::move(computed.first),
            std::move(output),
            computed.nodeIdSpace};
    }
}

template<std::floating_point OutputReal, std::floating_point InternalReal>
[[nodiscard]] inline ComputedAttributeDataWithDelta<OutputReal> castComputedAttributeDataWithDelta(ComputedAttributeDataWithDelta<InternalReal> computed) {
    if constexpr (std::is_same_v<OutputReal, InternalReal>) {
        return computed;
    } else {
        std::vector<OutputReal> output;
        output.reserve(computed.second.size());
        for (const InternalReal value : computed.second) {
            output.push_back(static_cast<OutputReal>(value));
        }
        return ComputedAttributeDataWithDelta<OutputReal>{
            std::move(computed.first),
            std::move(output),
            computed.nodeIdSpace};
    }
}

template<std::floating_point OutputReal, std::floating_point InternalReal>
[[nodiscard]] inline std::vector<OutputReal> castAttributeValues(std::vector<InternalReal> values) {
    if constexpr (std::is_same_v<OutputReal, InternalReal>) {
        return values;
    } else {
        std::vector<OutputReal> output;
        output.reserve(values.size());
        for (const InternalReal value : values) {
            output.push_back(static_cast<OutputReal>(value));
        }
        return output;
    }
}

template<std::floating_point OutputReal, std::floating_point InternalReal>
[[nodiscard]] inline ImagePtr<OutputReal> mapNodeAttributeToImageCast(
    const MorphologicalTree& tree,
    const AttributeNames& attrNames,
    std::span<const InternalReal> nodeValues,
    Attribute attribute)
{
    ImagePtr<OutputReal> imgPtr = Image<OutputReal>::create(tree.getNumRowsOfImage(), tree.getNumColsOfImage());
    OutputReal* img = imgPtr->rawData();
    for (int p = 0; p < imgPtr->getSize(); ++p) {
        const NodeId nodeId = tree.getProperPartOwner(p);
        img[p] = static_cast<OutputReal>(nodeValues[attrNames.linearIndex(nodeId, attribute)]);
    }
    return imgPtr;
}

} // namespace detail

template <std::floating_point Real>
inline ComputedAttributeData<Real> AttributeComputation::computeSingleTopologyAttribute(const MorphologicalTree& tree, AttributeOrGroup attrOrGroup, NodeIdSpace outputSpace) {
    return detail::castComputedAttributeData<Real>(
        detail::materializeAttributesWithoutAltitude<double>(tree, {attrOrGroup}, outputSpace));
}

template <std::floating_point Real>
inline ComputedAttributeData<Real> AttributeComputation::computeTopologyAttributes(const MorphologicalTree& tree, const std::vector<AttributeOrGroup>& attributes, NodeIdSpace outputSpace) {
    return detail::castComputedAttributeData<Real>(
        detail::materializeAttributesWithoutAltitude<double>(tree, attributes, outputSpace));
}

template<std::floating_point Real, AltitudeValue T>
inline ComputedAttributeData<Real> AttributeComputation::computeSingleTopologyAttribute(const WeightedMorphologicalTree<T>& tree, AttributeOrGroup attrOrGroup, NodeIdSpace outputSpace) {
    return computeSingleTopologyAttribute<Real>(tree.asView(), attrOrGroup, outputSpace);
}

template<std::floating_point Real, AltitudeValue T>
inline ComputedAttributeData<Real> AttributeComputation::computeSingleTopologyAttribute(const WeightedTreeView<T>& tree, AttributeOrGroup attrOrGroup, NodeIdSpace outputSpace) {
    tree.requireTopologyUnchanged("AttributeComputation::computeSingleTopologyAttribute");
    return detail::castComputedAttributeData<Real>(
        detail::materializeTopologyAttributeRequest<double>(
            tree.topology(),
            tree.altitude(),
            {attrOrGroup},
            detail::DependencyMapT<double>{},
            outputSpace));
}

template<std::floating_point Real, AltitudeValue T>
inline ComputedAttributeData<Real> AttributeComputation::computeTopologyAttributes(const WeightedMorphologicalTree<T>& tree, const std::vector<AttributeOrGroup>& attributes, NodeIdSpace outputSpace) {
    return computeTopologyAttributes<Real>(tree.asView(), attributes, outputSpace);
}

template<std::floating_point Real, AltitudeValue T>
inline ComputedAttributeData<Real> AttributeComputation::computeTopologyAttributes(const WeightedTreeView<T>& tree, const std::vector<AttributeOrGroup>& attributes, NodeIdSpace outputSpace) {
    tree.requireTopologyUnchanged("AttributeComputation::computeTopologyAttributes");
    return detail::castComputedAttributeData<Real>(
        detail::materializeTopologyAttributeRequest<double>(
            tree.topology(),
            tree.altitude(),
            attributes,
            detail::DependencyMapT<double>{},
            outputSpace));
}

template<std::floating_point Real, AltitudeValue T>
inline ComputedAttributeData<Real> AttributeComputation::computeSingleAttribute(const WeightedMorphologicalTree<T>& tree, AttributeOrGroup attrOrGroup, NodeIdSpace outputSpace) {
    return computeSingleAttribute<Real>(tree.asView(), attrOrGroup, outputSpace);
}

template<std::floating_point Real, AltitudeValue T>
inline ComputedAttributeData<Real> AttributeComputation::computeSingleAttribute(const WeightedTreeView<T>& tree, AttributeOrGroup attrOrGroup, NodeIdSpace outputSpace) {
    return computeAttributesFromAltitudeSpan<Real>(tree, {attrOrGroup}, outputSpace);
}

template<std::floating_point Real, AltitudeValue T>
inline ComputedAttributeDataWithDelta<Real> AttributeComputation::computeSingleAttributeWithDelta(const WeightedMorphologicalTree<T>& tree, Attribute attribute, AltitudeDiff<T> deltaStep, int radius, std::string padding, NodeIdSpace outputSpace) {
    return computeSingleAttributeWithDelta<Real>(tree.asView(), attribute, deltaStep, radius, std::move(padding), outputSpace);
}

template<std::floating_point Real, AltitudeValue T>
inline ComputedAttributeDataWithDelta<Real> AttributeComputation::computeSingleAttributeWithDelta(const WeightedTreeView<T>& tree, Attribute attribute, AltitudeDiff<T> deltaStep, int radius, std::string padding, NodeIdSpace outputSpace) {
    tree.requireTopologyUnchanged("AttributeComputation::computeSingleAttributeWithDelta");
    auto base = computeSingleAttribute<double>(tree, attribute, NodeIdSpace::MORPHOLOGICAL_TREE);

    return detail::castComputedAttributeDataWithDelta<Real>(
        detail::materializeSingleAttributeWithTypedDelta<double>(
            tree.topology(),
            tree.altitude(),
            std::move(base),
            attribute,
            deltaStep,
            radius,
            std::move(padding),
            outputSpace));
}

template<std::floating_point Real, AltitudeValue T>
inline ComputedAttributeData<Real> AttributeComputation::computeAttributes(const WeightedMorphologicalTree<T>& tree, const std::vector<AttributeOrGroup>& attributes, NodeIdSpace outputSpace) {
    return computeAttributes<Real>(tree.asView(), attributes, outputSpace);
}

template<std::floating_point Real, AltitudeValue T>
inline ComputedAttributeData<Real> AttributeComputation::computeAttributes(const WeightedTreeView<T>& tree, const std::vector<AttributeOrGroup>& attributes, NodeIdSpace outputSpace) {
    tree.requireTopologyUnchanged("AttributeComputation::computeAttributes");
    return computeAttributesFromAltitudeSpan<Real>(tree, attributes, outputSpace);
}

template<std::floating_point Real, AltitudeValue T>
inline ComputedAttributeData<Real> AttributeComputation::computeAttributesFromAltitudeSpan(
    const WeightedTreeView<T>& weighted,
    const std::vector<AttributeOrGroup>& attributes,
    NodeIdSpace outputSpace) {
    weighted.requireTopologyUnchanged("AttributeComputation::computeAttributesFromAltitudeSpan");
    return detail::castComputedAttributeData<Real>(
        detail::materializeAttributes<double>(weighted.topology(), weighted.altitude(), attributes, outputSpace));
}

template<std::floating_point Real, AltitudeValue T>
inline std::vector<Real> AttributeComputation::projectNodeValuesToExportedHigra(
    const WeightedMorphologicalTree<T>& tree,
    const AttributeNames& attrNames,
    std::span<const Real> nodeValues) {
    return projectNodeValuesToExportedHigra<Real>(tree.asView(), attrNames, nodeValues);
}

template<std::floating_point Real, AltitudeValue T>
inline std::vector<Real> AttributeComputation::projectNodeValuesToExportedHigra(
    const WeightedMorphologicalTree<T>& tree,
    const AttributeNames& attrNames,
    const std::vector<Real>& nodeValues) {
    return projectNodeValuesToExportedHigra<Real>(tree, attrNames, std::span<const Real>(nodeValues));
}

template<std::floating_point Real, AltitudeValue T>
inline std::vector<Real> AttributeComputation::projectNodeValuesToExportedHigra(
    const WeightedTreeView<T>& tree,
    const AttributeNames& attrNames,
    std::span<const Real> nodeValues) {
    tree.topology().requireNotEditing("AttributeComputation::projectNodeValuesToExportedHigra");
    tree.requireTopologyUnchanged("AttributeComputation::projectNodeValuesToExportedHigra");
    if constexpr (std::is_same_v<Real, double>) {
        return detail::projectNodeValuesToExportedHigraTyped<double>(
            tree.topology(),
            tree.altitude(),
            attrNames,
            nodeValues);
    } else {
        std::vector<double> internalValues;
        internalValues.reserve(nodeValues.size());
        for (const Real value : nodeValues) {
            internalValues.push_back(static_cast<double>(value));
        }
        return detail::castAttributeValues<Real>(
            detail::projectNodeValuesToExportedHigraTyped<double>(
                tree.topology(),
                tree.altitude(),
                attrNames,
                std::span<const double>(internalValues)));
    }
}

template<std::floating_point Real, AltitudeValue T>
inline std::vector<Real> AttributeComputation::projectNodeValuesToExportedHigra(
    const WeightedTreeView<T>& tree,
    const AttributeNames& attrNames,
    const std::vector<Real>& nodeValues) {
    return projectNodeValuesToExportedHigra<Real>(tree, attrNames, std::span<const Real>(nodeValues));
}

template<std::floating_point Real, AltitudeValue T>
inline ImagePtr<Real> AttributeComputation::computeAttributeMapping(const WeightedMorphologicalTree<T>& tree, Attribute attribute) {
    auto [attrNames, buffer] = AttributeComputation::computeSingleAttribute<double>(tree, attribute, NodeIdSpace::MORPHOLOGICAL_TREE);
    return detail::mapNodeAttributeToImageCast<Real>(tree.topology(), attrNames, std::span<const double>(buffer), attribute);
}

template<std::floating_point Real, AltitudeValue T>
inline ImagePtr<Real> AttributeComputation::computeAttributeMapping(const WeightedTreeView<T>& tree, Attribute attribute) {
    tree.requireTopologyUnchanged("AttributeComputation::computeAttributeMapping");
    auto [attrNames, buffer] = AttributeComputation::computeSingleAttribute<double>(tree, attribute, NodeIdSpace::MORPHOLOGICAL_TREE);
    return detail::mapNodeAttributeToImageCast<Real>(tree.topology(), attrNames, std::span<const double>(buffer), attribute);
}

} // namespace mmcfilters
