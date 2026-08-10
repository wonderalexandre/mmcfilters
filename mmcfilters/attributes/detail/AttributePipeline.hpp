#pragma once

#include "AttributeFamilyScheduler.hpp"
#include "AttributeCapabilityValidation.hpp"
#include "AttributeProjection.hpp"
#include "AttributeRequestUtils.hpp"
#include "TopologyAttributeBackend.hpp"
#include "../AttributeNames.hpp"
#include "../AttributeRegistry.hpp"
#include "../AttributeResultTypes.hpp"
#include "../computers/AreaComputer.hpp"
#include "../computers/GrayLevelStatsComputer.hpp"
#include "../computers/VolumeComputer.hpp"
#include "../computers/MaxDistComputer.hpp"
#include "../../trees/MorphologicalTree.hpp"
#include "../../trees/TreeAltitudeAlgorithms.hpp"
#include "../../utils/Altitude.hpp"
#include "../../utils/Common.hpp"
#include "../../utils/Contract.hpp"

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mmcfilters::detail {

/**
 * @brief Preferred internal orchestration path for public attribute requests.
 *
 * @details
 * This pipeline handles ordinary requests without injected dependencies. It
 * computes altitude-dependent attributes with typed kernels and delegates
 * altitude-independent geometry/topology attributes to `TopologyAttributeBackend`.
 * New public-facing attribute orchestration should be integrated here instead
 * of adding a parallel execution pipeline.
 *
 * @param attribute Attribute requested by the operation.
 * @return True when the documented condition holds; otherwise false.
 */

inline bool isAttributePipelineAltitudeAttribute(Attribute attribute) noexcept { return attributes::registry::isAttributePipelineAltitudeAttribute(attribute); }

/**
 * @brief Tests whether attribute pipeline attribute holds.
 *
 * @param attribute Attribute requested by the operation.
 * @return True when attribute pipeline attribute; otherwise false.
 */
inline bool isAttributePipelineAttribute(Attribute attribute) noexcept { return attributes::registry::isPipelineComputed(attribute); }

/**
 * @brief Checks whether an attribute-pipeline request requires hierarchy altitudes.
 *
 * @param attribute Attribute requested by the operation.
 * @return True when the documented condition holds; otherwise false.
 */
inline bool requiresAltitudeForAttributePipeline(Attribute attribute) noexcept { return attributes::registry::requiresAltitude(attribute); }

/**
 * @brief Tests whether attribute pipeline request holds.
 *
 * @param attributes Attributes requested by the operation.
 * @return True when attribute pipeline request; otherwise false.
 */
inline bool isAttributePipelineRequest(const std::vector<AttributeOrGroup>& attributes) {
    return requestContainsOnlyAttributes(attributes, [](Attribute attribute) { return isAttributePipelineAttribute(attribute); });
}

/**
 * @brief Expands attribute pipeline attributes.
 *
 * @param attributes Attributes requested by the operation.
 * @return Values produced by the operation.
 */
inline std::vector<Attribute> expandAttributePipelineAttributes(const std::vector<AttributeOrGroup>& attributes) {
    return expandUniqueAttributeRequest(
        attributes, [](Attribute attribute) { return isAttributePipelineAttribute(attribute); },
        "computeAttributesFromAltitudeSpan received an unsupported attribute.",
        "computeAttributesFromAltitudeSpan received an attribute group containing unsupported attributes.");
}

/**
 * @brief Tests whether altitude dependent attribute holds.
 *
 * @param attributes Attributes requested by the operation.
 * @return True when altitude dependent attribute; otherwise false.
 */
inline bool containsAltitudeDependentAttribute(std::span<const Attribute> attributes) {
    return std::any_of(attributes.begin(), attributes.end(), [](Attribute attribute) { return requiresAltitudeForAttributePipeline(attribute); });
}

/** @brief Builds a scalar layout after request expansion established uniqueness. */
inline AttributeNames makeAttributeNamesForExpandedRequest(std::span<const Attribute> attributes) {
    std::unordered_map<Attribute, int> offsets;
    offsets.reserve(attributes.size());
    for (std::size_t index = 0; index < attributes.size(); ++index) {
        offsets.emplace(attributes[index], static_cast<int>(index));
    }
    return AttributeNames(std::move(offsets));
}

/**
 * @brief Executes an altitude-aware attribute plan into the caller result buffer.
 *
 * @details
 * The function assumes `resultNames` describes only the public requested
 * attributes, while `plan` may also contain hidden dependencies. Hidden buffers
 * are materialized in internal dense `NodeId` space and are used only as
 * dependency sources for downstream families. Public attributes are written into
 * `resultBuffer` according to `resultNames`.
 *
 * Execution order follows the scheduler closure: area first when needed,
 * volume before gray-level statistics, `MAX_DIST` with altitude, then
 * topology/support families through `TopologyAttributeBackend`.
 *
 * @param tree Tree topology used by the operation.
 * @param altitude Altitude data indexed by node identifier.
 * @param plan Attribute computation plan.
 * @param resultNames Destination represented by `resultNames`.
 * @param resultBuffer Destination represented by `resultBuffer`.
 */
template <std::floating_point Real, AltitudeValue T>
inline void executeAttributeComputationPlan(const MorphologicalTree& tree, std::span<const T> altitude, const AttributeComputationPlan& plan,
                                            const AttributeNames& resultNames, std::span<Real> resultBuffer) {
    const int numNodes = tree.getNumInternalNodeSlots();

    const AttributeNames areaNames = AttributeNames::fromList({AREA});
    std::vector<Real> areaBuffer;
    DependencySourceT<Real> areaSource{};
    if (plan.materializes(AREA)) {
        if (plan.requests(AREA)) {
            attributes::computers::detail::kernel::computeArea(
                AttributeComputeContext<Real>{tree, resultBuffer, resultNames, std::span<const Attribute>{}});
            areaSource = DependencySourceT<Real>{&resultNames, resultBuffer.data()};
        } else {
            areaBuffer.assign(static_cast<std::size_t>(numNodes), Real{0});
            attributes::computers::detail::kernel::computeArea(
                AttributeComputeContext<Real>{tree, std::span<Real>(areaBuffer), areaNames, std::span<const Attribute>{}});
            areaSource = DependencySourceT<Real>{&areaNames, areaBuffer.data()};
        }
    }

    std::vector<Real> volumeDependencyBuffer;
    std::unique_ptr<AttributeNames> volumeDependencyNames;
    DependencySourceT<Real> volumeSource{};
    const std::vector<Attribute> volumeAttributes = plan.materializedForFamily(attributes::computers::AttributeComputerFamily::Volume);
    if (!volumeAttributes.empty()) {
        const bool hasHiddenVolumeDependency =
            std::any_of(volumeAttributes.begin(), volumeAttributes.end(), [&](Attribute attribute) { return plan.hides(attribute); });

        const bool canWriteVolumeDirectly = !hasHiddenVolumeDependency;
        const AttributeNames* volumeNames = &resultNames;
        std::span<Real> volumeBuffer(resultBuffer);
        if (!canWriteVolumeDirectly) {
            volumeDependencyNames = std::make_unique<AttributeNames>(AttributeNames::fromList(volumeAttributes));
            volumeDependencyBuffer.assign(static_cast<std::size_t>(numNodes) * static_cast<std::size_t>(volumeDependencyNames->NUM_ATTRIBUTES), Real{0});
            volumeNames = volumeDependencyNames.get();
            volumeBuffer = std::span<Real>(volumeDependencyBuffer);
        }

        std::array<DependencySourceT<Real>, 1> volumeDependencies{{areaSource}};
        const bool needsAreaForVolume = anyAttributeRequiresDependency(std::span<const Attribute>(volumeAttributes), AREA);
        const std::span<const DependencySourceT<Real>> volumeDependencySpan =
            needsAreaForVolume ? std::span<const DependencySourceT<Real>>(volumeDependencies) : std::span<const DependencySourceT<Real>>();

        const auto volumeContext = AltitudeAttributeComputeContext<Real, T>{
            tree, altitude, volumeBuffer, *volumeNames, std::span<const Attribute>(volumeAttributes), volumeDependencySpan};
        const auto volumeRequest = ::mmcfilters::attributes::computers::detail::VolumeRequest::from(volumeContext.requestedAttributes);
        const DependencySourceT<Real>* areaDependency = needsAreaForVolume ? &volumeDependencies[0] : nullptr;
        ::mmcfilters::attributes::computers::detail::kernel::computeVolume(volumeContext, volumeRequest, areaDependency);

        if (!canWriteVolumeDirectly) {
            for (Attribute attribute : volumeAttributes) {
                if (plan.requests(attribute)) {
                    copyAttributeValuesBetweenLayouts<Real>(tree, *volumeNames, volumeDependencyBuffer, resultNames, resultBuffer, attribute);
                }
            }
        }

        if (plan.materializes(VOLUME)) {
            volumeSource = canWriteVolumeDirectly ? DependencySourceT<Real>{&resultNames, resultBuffer.data()}
                                                  : DependencySourceT<Real>{volumeNames, volumeDependencyBuffer.data()};
        }
    }

    const std::vector<Attribute> grayAttributes = plan.requestedForFamily(attributes::computers::AttributeComputerFamily::GrayLevelStats);
    if (!grayAttributes.empty()) {
        std::array<DependencySourceT<Real>, 2> grayDependencies{{volumeSource, areaSource}};
        const bool needsGrayAggregateDependencies = anyAttributeRequiresDependency(std::span<const Attribute>(grayAttributes), VOLUME) ||
                                                    anyAttributeRequiresDependency(std::span<const Attribute>(grayAttributes), AREA);
        const std::span<const DependencySourceT<Real>> grayDependencySpan =
            needsGrayAggregateDependencies ? std::span<const DependencySourceT<Real>>(grayDependencies) : std::span<const DependencySourceT<Real>>();

        const auto grayContext = AltitudeAttributeComputeContext<Real, T>{
            tree, altitude, std::span<Real>(resultBuffer), resultNames, std::span<const Attribute>(grayAttributes), grayDependencySpan};
        const auto grayRequest = ::mmcfilters::attributes::computers::detail::GrayLevelStatsRequest::from(grayContext.requestedAttributes);
        const DependencySourceT<Real>* volumeDependency = needsGrayAggregateDependencies ? &grayDependencies[0] : nullptr;
        const DependencySourceT<Real>* grayAreaDependency = needsGrayAggregateDependencies ? &grayDependencies[1] : nullptr;
        ::mmcfilters::attributes::computers::detail::kernel::computeGrayLevelStats(grayContext, grayRequest, volumeDependency, grayAreaDependency);
    }

    const std::vector<Attribute> maxDistAttributes = plan.requestedForFamily(attributes::computers::AttributeComputerFamily::MaxDist);
    if (!maxDistAttributes.empty()) {
        ::mmcfilters::attributes::computers::detail::kernel::computeMaxDist(
            AltitudeAttributeComputeContext<Real, T>{tree, altitude, resultBuffer, resultNames, std::span<const Attribute>(maxDistAttributes)});
    }

    DependencyMapT<Real> topologyOnlyDependencies;
    if (areaSource.attrNames != nullptr && areaSource.buffer != nullptr) {
        topologyOnlyDependencies[AREA] = ComputedAttributeViewT<Real>{areaSource.attrNames, areaSource.buffer, NodeIdSpace::MORPHOLOGICAL_TREE};
    }
    computeTopologyOnlyAttributesIntoResult<Real>(tree, std::span<const Attribute>(plan.requestedAttributes), std::move(topologyOnlyDependencies), resultNames,
                                                  resultBuffer, altitude);
}

/**
 * @brief Materializes a topology-only request from an unweighted tree.
 *
 * @param tree Tree topology used by the operation.
 * @param attributes Attributes requested by the operation.
 * @param outputSpace Node-id domain used to index the output.
 * @return The materialized topology-only request from an unweighted tree.
 *
 * @throws std::invalid_argument if the expanded request contains any
 * altitude-dependent attribute.
 *
 */
template <std::floating_point Real = float>
inline ComputedAttributeData<Real> materializeAttributesWithoutAltitude(const MorphologicalTree& tree, const std::vector<AttributeOrGroup>& attributes,
                                                                        NodeIdSpace outputSpace) {
    const std::vector<Attribute> requestedAttributes = expandAttributePipelineAttributes(attributes);
    const std::span<const Attribute> requestedSpan(requestedAttributes);
    MMCFILTERS_CONTRACT_CHECKED_ONLY(tree.requireNotEditing("AttributeComputation::computeAttributes");
                                     validateAttributeCapabilities(tree, requestedSpan, false, "AttributeComputation");
                                     if (containsAltitudeDependentAttribute(requestedSpan)) {
                                         throw std::invalid_argument(
                                             "Requested attributes require altitude data. Use a WeightedMorphologicalTree<T>, WeightedTreeView, or the overload "
                                             "that receives an altitude buffer.");
                                     });

    const AttributeNames resultNames = makeAttributeNamesForExpandedRequest(requestedAttributes);
    std::vector<Real> resultBuffer = makeAttributeValueBuffer<Real>(tree, resultNames);

    DependencyMapT<Real> topologyOnlyDependencies;
    if (containsAttribute(requestedSpan, AREA)) {
        attributes::computers::detail::kernel::computeArea(
            AttributeComputeContext<Real>{tree, std::span<Real>(resultBuffer), resultNames, std::span<const Attribute>{}});
        topologyOnlyDependencies[AREA] = ComputedAttributeViewT<Real>{&resultNames, resultBuffer.data(), NodeIdSpace::MORPHOLOGICAL_TREE};
    }

    computeTopologyOnlyAttributesIntoResult<Real>(tree, requestedSpan, std::move(topologyOnlyDependencies), resultNames, resultBuffer);

    return projectComputedDataToNodeIdSpace(tree, ComputedAttributeData<Real>{std::move(resultNames), std::move(resultBuffer), NodeIdSpace::MORPHOLOGICAL_TREE},
                                            outputSpace);
}

/**
 * @brief Materializes a mixed topology/altitude request from a typed altitude span.
 *
 * @details
 * This is the internal implementation shared by weighted-tree and
 * altitude-span public facades. It validates altitude shape, builds the
 * dependency plan, executes the plan in internal node-id space, and projects the
 * result only at the API boundary.
 *
 * @param tree Tree topology used by the operation.
 * @param altitude Altitude data indexed by node identifier.
 * @param attributes Attributes requested by the operation.
 * @param outputSpace Node-id domain used to index the output.
 * @return The materialized mixed topology/altitude request from a typed altitude span.
 */
template <std::floating_point Real = float, AltitudeValue T>
inline ComputedAttributeData<Real> materializeAttributesForExpandedRequest(const MorphologicalTree& tree, std::span<const T> altitude,
                                                                           const std::vector<Attribute>& requestedAttributes, NodeIdSpace outputSpace) {
    const AttributeNames resultNames = makeAttributeNamesForExpandedRequest(requestedAttributes);
    std::vector<Real> resultBuffer = makeAttributeValueBuffer<Real>(tree, resultNames);
    const AttributeComputationPlan plan = makeAttributeComputationPlan(std::span<const Attribute>(requestedAttributes));

    executeAttributeComputationPlan<Real, T>(tree, altitude, plan, resultNames, std::span<Real>(resultBuffer));

    return projectComputedDataToNodeIdSpace(
        tree, altitude, ComputedAttributeData<Real>{std::move(resultNames), std::move(resultBuffer), NodeIdSpace::MORPHOLOGICAL_TREE}, outputSpace);
}

/**
 * @brief Validates a public weighted request once and enters the established pipeline.
 */
template <std::floating_point Real = float, AltitudeValue T>
inline ComputedAttributeData<Real> materializeAttributes(const MorphologicalTree& tree, std::span<const T> altitude,
                                                         const std::vector<AttributeOrGroup>& attributes, NodeIdSpace outputSpace) {
    const std::vector<Attribute> requestedAttributes = expandAttributePipelineAttributes(attributes);
    MMCFILTERS_CONTRACT_CHECKED_ONLY(tree.requireNotEditing("AttributeComputation::computeAttributesFromAltitudeSpan");
                                     TreeAltitudeAlgorithms::validateAltitudeBufferShape(tree, altitude);
                                     validateAttributeCapabilities(tree, std::span<const Attribute>(requestedAttributes), true, "AttributeComputation");
                                     if (containsAttribute(std::span<const Attribute>(requestedAttributes), MAX_DIST)) {
                                         ::mmcfilters::attributes::computers::detail::validateFiniteMaxDistAltitude(altitude);
                                     });
    return materializeAttributesForExpandedRequest<Real, T>(tree, altitude, requestedAttributes, outputSpace);
}

} // namespace mmcfilters::detail
