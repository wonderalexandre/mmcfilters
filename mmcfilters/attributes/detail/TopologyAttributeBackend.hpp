#pragma once

#include "AttributeDependencyCache.hpp"
#include "AttributeFamilyScheduler.hpp"
#include "AttributeProjection.hpp"
#include "AttributeRequestUtils.hpp"
#include "../AttributeNames.hpp"
#include "../AttributeRegistry.hpp"
#include "../AttributeResultTypes.hpp"
#include "../computers/AttributeComputerRegistry.hpp"
#include "../computers/AreaComputer.hpp"
#include "../computers/BoundingBoxComputer.hpp"
#include "../computers/BitquadAttributeComputer.hpp"
#include "../computers/ContourSideAttributeComputer.hpp"
#include "../computers/MomentBasedAttributeComputer.hpp"
#include "../computers/TreeTopologyComputer.hpp"
#include "../../trees/MorphologicalTree.hpp"
#include "../../utils/Altitude.hpp"
#include "../../utils/Common.hpp"

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcfilters::detail {

/**
 * @brief Backend for altitude-independent attributes used by the attribute pipeline.
 *
 * @details
 * The backend owns topology-only dependency materialisation for attributes used
 * by `AttributePipeline`. Topology-only families are handled explicitly here so
 * public requests do not rely on a second dispatch/orchestration path.
 *
 * The backend always computes in dense internal `NodeId` space. Output-space
 * projection is handled after materialization by the public facade helpers.
 */

template <std::floating_point Real = float>
inline ComputedAttributeData<Real> materializeTopologyAttributeRequest(const MorphologicalTree& tree, const std::vector<AttributeOrGroup>& attributes,
                                                                       const DependencyMapT<Real>& availableDeps, NodeIdSpace outputSpace);

template <std::floating_point Real = float, AltitudeValue T>
inline ComputedAttributeData<Real> materializeTopologyAttributeRequest(const MorphologicalTree& tree, std::span<const T> altitude,
                                                                       const std::vector<AttributeOrGroup>& attributes,
                                                                       const DependencyMapT<Real>& availableDeps, NodeIdSpace outputSpace);

/**
 * @brief Tests whether topology only attribute holds.
 *
 * @param attribute Attribute requested by the operation.
 * @return True when topology only attribute; otherwise false.
 */
inline bool isTopologyOnlyAttribute(Attribute attribute) noexcept { return attributes::registry::isTopologyOnly(attribute); }

/**
 * @brief Tests whether topology backend attribute holds.
 *
 * @param attribute Attribute requested by the operation.
 * @return True when topology backend attribute; otherwise false.
 */
inline bool isTopologyBackendAttribute(Attribute attribute) noexcept { return attribute == AREA || isTopologyOnlyAttribute(attribute); }

/**
 * @brief Expands topology backend attributes.
 *
 * @param attributes Attributes requested by the operation.
 * @return Values produced by the operation.
 */
inline std::vector<Attribute> expandTopologyBackendAttributes(const std::vector<AttributeOrGroup>& attributes) {
    return expandUniqueAttributeRequest(
        attributes, [](Attribute attribute) { return isTopologyBackendAttribute(attribute); },
        "TopologyAttributeBackend received an attribute that requires altitude.",
        "TopologyAttributeBackend received an attribute group containing altitude-dependent attributes.");
}

/**
 * @brief Sentinel used when bitquad projection does not have altitude data.
 */
struct NoTopologyBackendAltitude {};

/**
 * @brief Computes bitquad scalar attributes for trees with regular adjacency.
 *
 * @param tree Tree topology used by the operation.
 * @param bitquadAttributes Attribute information represented by `bitquadAttributes`.
 * @param resultNames Destination represented by `resultNames`.
 * @param resultBuffer Destination represented by `resultBuffer`.
 */
template <std::floating_point Real>
inline void computeBitquadBackendAttributesIntoResult(const MorphologicalTree& tree, std::span<const Attribute> bitquadAttributes,
                                                      const AttributeNames& resultNames, std::span<Real> resultBuffer, NoTopologyBackendAltitude) {
    ::mmcfilters::attributes::computers::BitquadAttributeComputer::compute(AttributeComputeContext<Real>{tree, resultBuffer, resultNames, bitquadAttributes});
}

/**
 * @brief Computes bitquad scalar attributes using altitude when ToS polarity is needed.
 *
 * @param tree Tree topology used by the operation.
 * @param bitquadAttributes Attribute information represented by `bitquadAttributes`.
 * @param resultNames Destination represented by `resultNames`.
 * @param resultBuffer Destination represented by `resultBuffer`.
 * @param altitude Altitude data indexed by node identifier.
 */
template <std::floating_point Real, AltitudeValue T>
inline void computeBitquadBackendAttributesIntoResult(const MorphologicalTree& tree, std::span<const Attribute> bitquadAttributes,
                                                      const AttributeNames& resultNames, std::span<Real> resultBuffer, std::span<const T> altitude) {
    ::mmcfilters::attributes::computers::BitquadAttributeComputer::compute(
        AltitudeAttributeComputeContext<Real, T>{tree, altitude, resultBuffer, resultNames, bitquadAttributes});
}

/**
 * @brief Executes topology/support families into an already allocated result buffer.
 *
 * @details
 * `plan.requestedAttributes` determines the public columns to fill. Hidden
 * dependencies are materialized in owned scratch buffers and registered in the
 * local dependency cache, but they are not copied to the public result unless
 * requested directly. `BitquadAltitude` is either `NoTopologyBackendAltitude`
 * or a typed altitude span for Tree of Shapes bitquad projection.
 *
 * @param tree Tree topology used by the operation.
 * @param plan Attribute computation plan.
 * @param available Available computed-attribute views.
 * @param resultNames Destination represented by `resultNames`.
 * @param resultBuffer Destination represented by `resultBuffer`.
 * @param bitquadAltitude Altitude or level represented by `bitquadAltitude`.
 */
template <std::floating_point Real, class BitquadAltitude>
inline void executeTopologyAttributeComputationPlan(const MorphologicalTree& tree, const AttributeComputationPlan& plan, DependencyMapT<Real> available,
                                                    const AttributeNames& resultNames, std::span<Real> resultBuffer, BitquadAltitude bitquadAltitude) {
    OwnedComputedResultsT<Real> ownedResults;

    auto ensureInternalAreaDependency = [&]() -> DependencySourceT<Real> {
        const auto it = available.find(AREA);
        if (it != available.end() && isReusableDependencyData(it->second, {AREA})) {
            return it->second.dependencySource();
        }

        AttributeNames areaNames = AttributeNames::fromList({AREA});
        std::vector<Real> areaBuffer = makeAttributeValueBuffer<Real>(tree, areaNames);
        attributes::computers::AreaComputer::compute(AttributeComputeContext<Real>{tree, std::span<Real>(areaBuffer), areaNames, std::span<const Attribute>{}});
        stashComputedAttributes(ownedResults, available, ComputedAttributeData<Real>(std::move(areaNames), std::move(areaBuffer)));
        return available.at(AREA).dependencySource();
    };

    auto ensureInternalCentralMomentDependency = [&](const std::vector<Attribute>& requiredAttributes) -> DependencySourceT<Real> {
        for (const Attribute attribute : requiredAttributes) {
            const auto it = available.find(attribute);
            if (it != available.end() && isReusableDependencyData(it->second, requiredAttributes)) {
                return it->second.dependencySource();
            }
        }

        AttributeNames momentNames = AttributeNames::fromList(requiredAttributes);
        std::vector<Real> momentBuffer = makeAttributeValueBuffer<Real>(tree, momentNames);
        attributes::computers::CentralMomentsComputer::compute(
            AttributeComputeContext<Real>{tree, std::span<Real>(momentBuffer), momentNames, std::span<const Attribute>(requiredAttributes)});
        stashComputedAttributes(ownedResults, available, ComputedAttributeData<Real>(std::move(momentNames), std::move(momentBuffer)));
        return available.at(requiredAttributes.front()).dependencySource();
    };

    const std::vector<Attribute> boundingBoxAttributes = plan.requestedForFamily(attributes::computers::AttributeComputerFamily::BoundingBox);
    const std::vector<Attribute> treeTopologyAttributes = plan.requestedForFamily(attributes::computers::AttributeComputerFamily::TreeTopology);
    const std::vector<Attribute> centralMomentAttributes = plan.requestedForFamily(attributes::computers::AttributeComputerFamily::CentralMoments);
    const std::vector<Attribute> huMomentAttributes = plan.requestedForFamily(attributes::computers::AttributeComputerFamily::HuMoments);
    const std::vector<Attribute> momentBasedAttributes = plan.requestedForFamily(attributes::computers::AttributeComputerFamily::MomentDerived);
    const std::vector<Attribute> bitquadAttributes = plan.requestedForFamily(attributes::computers::AttributeComputerFamily::Bitquad);
    const std::vector<Attribute> contourAttributes = plan.requestedForFamily(attributes::computers::AttributeComputerFamily::ContourSide);
    std::vector<AttributeOrGroup> topologyOnlyRequests;

    for (const Attribute attribute : plan.requestedAttributes) {
        if (isTopologyOnlyAttribute(attribute) && familyForAttribute(attribute) == attributes::computers::AttributeComputerFamily::Unsupported) {
            topologyOnlyRequests.emplace_back(attribute);
        }
    }

    if (!boundingBoxAttributes.empty()) {
        const bool needsAreaDependency = anyAttributeRequiresDependency(std::span<const Attribute>(boundingBoxAttributes), AREA);
        std::array<DependencySourceT<Real>, 1> areaDependency{{}};
        std::span<const DependencySourceT<Real>> dependencySources;
        if (needsAreaDependency) {
            areaDependency[0] = ensureInternalAreaDependency();
            dependencySources = std::span<const DependencySourceT<Real>>(areaDependency);
        }

        attributes::computers::BoundingBoxComputer::compute(
            AttributeComputeContext<Real>{tree, resultBuffer, resultNames, std::span<const Attribute>(boundingBoxAttributes), dependencySources});

        for (const Attribute attribute : boundingBoxAttributes) {
            available[attribute] = ComputedAttributeViewT<Real>{&resultNames, resultBuffer.data(), NodeIdSpace::MORPHOLOGICAL_TREE};
        }
    }

    if (!treeTopologyAttributes.empty()) {
        attributes::computers::TreeTopologyComputer::compute(
            AttributeComputeContext<Real>{tree, resultBuffer, resultNames, std::span<const Attribute>(treeTopologyAttributes)});

        for (const Attribute attribute : treeTopologyAttributes) {
            available[attribute] = ComputedAttributeViewT<Real>{&resultNames, resultBuffer.data(), NodeIdSpace::MORPHOLOGICAL_TREE};
        }
    }

    if (!centralMomentAttributes.empty()) {
        attributes::computers::CentralMomentsComputer::compute(
            AttributeComputeContext<Real>{tree, resultBuffer, resultNames, std::span<const Attribute>(centralMomentAttributes)});

        for (const Attribute attribute : centralMomentAttributes) {
            available[attribute] = ComputedAttributeViewT<Real>{&resultNames, resultBuffer.data(), NodeIdSpace::MORPHOLOGICAL_TREE};
        }
    }

    if (!huMomentAttributes.empty()) {
        static const std::vector<Attribute> requiredCentralMoments{CENTRAL_MOMENT_20, CENTRAL_MOMENT_02, CENTRAL_MOMENT_11, CENTRAL_MOMENT_30,
                                                                   CENTRAL_MOMENT_03, CENTRAL_MOMENT_21, CENTRAL_MOMENT_12};
        std::array<DependencySourceT<Real>, 2> dependencies{{ensureInternalCentralMomentDependency(requiredCentralMoments), ensureInternalAreaDependency()}};

        attributes::computers::HuMomentsComputer::compute(AttributeComputeContext<Real>{
            tree, resultBuffer, resultNames, std::span<const Attribute>(huMomentAttributes), std::span<const DependencySourceT<Real>>(dependencies)});

        for (const Attribute attribute : huMomentAttributes) {
            available[attribute] = ComputedAttributeViewT<Real>{&resultNames, resultBuffer.data(), NodeIdSpace::MORPHOLOGICAL_TREE};
        }
    }

    if (!momentBasedAttributes.empty()) {
        static const std::vector<Attribute> requiredCentralMoments{CENTRAL_MOMENT_20, CENTRAL_MOMENT_02, CENTRAL_MOMENT_11};
        std::array<DependencySourceT<Real>, 2> dependencies{{ensureInternalCentralMomentDependency(requiredCentralMoments), ensureInternalAreaDependency()}};

        attributes::computers::MomentBasedAttributeComputer::compute(AttributeComputeContext<Real>{
            tree, resultBuffer, resultNames, std::span<const Attribute>(momentBasedAttributes), std::span<const DependencySourceT<Real>>(dependencies)});

        for (const Attribute attribute : momentBasedAttributes) {
            available[attribute] = ComputedAttributeViewT<Real>{&resultNames, resultBuffer.data(), NodeIdSpace::MORPHOLOGICAL_TREE};
        }
    }

    if (!bitquadAttributes.empty()) {
        computeBitquadBackendAttributesIntoResult(tree, std::span<const Attribute>(bitquadAttributes), resultNames, resultBuffer, bitquadAltitude);

        for (const Attribute attribute : bitquadAttributes) {
            available[attribute] = ComputedAttributeViewT<Real>{&resultNames, resultBuffer.data(), NodeIdSpace::MORPHOLOGICAL_TREE};
        }
    }

    if (!contourAttributes.empty()) {
        ::mmcfilters::attributes::computers::ContourSideAttributeComputer::compute(
            AttributeComputeContext<Real>{tree, resultBuffer, resultNames, std::span<const Attribute>(contourAttributes)});

        for (const Attribute attribute : contourAttributes) {
            available[attribute] = ComputedAttributeViewT<Real>{&resultNames, resultBuffer.data(), NodeIdSpace::MORPHOLOGICAL_TREE};
        }
    }

    if (!topologyOnlyRequests.empty()) {
        throw std::runtime_error("TopologyAttributeBackend received a topology-only attribute without an explicit pipeline family.");
    }

    for (const Attribute attribute : plan.requestedAttributes) {
        if (!isTopologyOnlyAttribute(attribute)) {
            continue;
        }
        const auto it = available.find(attribute);
        if (it == available.end()) {
            throw std::runtime_error("Requested topology-only attribute was not materialised.");
        }
        copyAttributesIntoBuffer(tree, it->second, {attribute}, resultNames, resultBuffer.data());
    }
}

/**
 * @brief Builds a scheduler plan and runs topology/support families in-place.
 *
 * @param tree Tree topology used by the operation.
 * @param requestedAttributes Attributes requested for materialization.
 * @param available Available computed-attribute views.
 * @param resultNames Destination represented by `resultNames`.
 * @param resultBuffer Destination represented by `resultBuffer`.
 * @param bitquadAltitude Altitude or level represented by `bitquadAltitude`.
 */
template <std::floating_point Real, class BitquadAltitude>
inline void computeTopologyOnlyAttributesIntoResultImpl(const MorphologicalTree& tree, std::span<const Attribute> requestedAttributes,
                                                        DependencyMapT<Real> available, const AttributeNames& resultNames, std::span<Real> resultBuffer,
                                                        BitquadAltitude bitquadAltitude) {
    executeTopologyAttributeComputationPlan(tree, makeAttributeComputationPlan(requestedAttributes), std::move(available), resultNames, resultBuffer,
                                            bitquadAltitude);
}

/**
 * @brief Computes topology/support attributes into an existing result buffer.
 *
 * @param tree Tree topology used by the operation.
 * @param requestedAttributes Attributes requested for materialization.
 * @param available Available computed-attribute views.
 * @param resultNames Destination represented by `resultNames`.
 * @param resultBuffer Destination represented by `resultBuffer`.
 */
template <std::floating_point Real>
inline void computeTopologyOnlyAttributesIntoResult(const MorphologicalTree& tree, std::span<const Attribute> requestedAttributes,
                                                    DependencyMapT<Real> available, const AttributeNames& resultNames, std::span<Real> resultBuffer) {
    computeTopologyOnlyAttributesIntoResultImpl(tree, requestedAttributes, std::move(available), resultNames, resultBuffer, NoTopologyBackendAltitude{});
}

/**
 * @brief Computes topology/support attributes with altitude available for ToS bitquads.
 *
 * @param tree Tree topology used by the operation.
 * @param requestedAttributes Attributes requested for materialization.
 * @param available Available computed-attribute views.
 * @param resultNames Destination represented by `resultNames`.
 * @param resultBuffer Destination represented by `resultBuffer`.
 * @param altitude Altitude data indexed by node identifier.
 */
template <std::floating_point Real, AltitudeValue T>
inline void computeTopologyOnlyAttributesIntoResult(const MorphologicalTree& tree, std::span<const Attribute> requestedAttributes,
                                                    DependencyMapT<Real> available, const AttributeNames& resultNames, std::span<Real> resultBuffer,
                                                    std::span<const T> altitude) {
    computeTopologyOnlyAttributesIntoResultImpl(tree, requestedAttributes, std::move(available), resultNames, resultBuffer, altitude);
}

/**
 * @brief Materializes a public topology/support-only request.
 *
 * @details
 * `availableDeps` can provide reusable internal-space dependencies, commonly
 * `AREA`, from an enclosing attribute pipeline execution. Missing dependencies
 * are computed locally and kept internal unless part of the public request.
 *
 * @param tree Tree topology used by the operation.
 * @param attributes Attributes requested by the operation.
 * @param availableDeps Available dependency views.
 * @param outputSpace Node-id domain used to index the output.
 * @return Computed attribute data produced by the backend.
 */
template <std::floating_point Real>
inline ComputedAttributeData<Real> materializeTopologyAttributeRequest(const MorphologicalTree& tree, const std::vector<AttributeOrGroup>& attributes,
                                                                       const DependencyMapT<Real>& availableDeps, NodeIdSpace outputSpace) {
    const std::vector<Attribute> requestedAttributes = expandTopologyBackendAttributes(attributes);
    const std::span<const Attribute> requestedSpan(requestedAttributes);
    const AttributeNames resultNames = AttributeNames::fromList(requestedAttributes);
    std::vector<Real> resultBuffer = makeAttributeValueBuffer<Real>(tree, resultNames);

    DependencyMapT<Real> available = availableDeps;
    if (containsAttribute(requestedSpan, AREA)) {
        const auto it = available.find(AREA);
        if (it != available.end() && isReusableDependencyData(it->second, {AREA})) {
            copyAttributesIntoBuffer(tree, it->second, {AREA}, resultNames, resultBuffer.data());
        } else {
            attributes::computers::AreaComputer::compute(
                AttributeComputeContext<Real>{tree, std::span<Real>(resultBuffer), resultNames, std::span<const Attribute>{}});
        }
        available[AREA] = ComputedAttributeViewT<Real>{&resultNames, resultBuffer.data(), NodeIdSpace::MORPHOLOGICAL_TREE};
    }

    computeTopologyOnlyAttributesIntoResult(tree, requestedSpan, std::move(available), resultNames, resultBuffer);

    return projectComputedDataToNodeIdSpace(tree, ComputedAttributeData<Real>{std::move(resultNames), std::move(resultBuffer), NodeIdSpace::MORPHOLOGICAL_TREE},
                                            outputSpace);
}

/**
 * @brief Materializes topology/support attributes with altitude available for projection.
 *
 * @details
 * The altitude span is not used by ordinary topology families, but is forwarded
 * to bitquad projection for Tree of Shapes inputs.
 *
 * @param tree Tree topology used by the operation.
 * @param altitude Altitude data indexed by node identifier.
 * @param attributes Attributes requested by the operation.
 * @param availableDeps Available dependency views.
 * @param outputSpace Node-id domain used to index the output.
 * @return The materialized topology/support attributes with altitude available for projection.
 */
template <std::floating_point Real, AltitudeValue T>
inline ComputedAttributeData<Real> materializeTopologyAttributeRequest(const MorphologicalTree& tree, std::span<const T> altitude,
                                                                       const std::vector<AttributeOrGroup>& attributes,
                                                                       const DependencyMapT<Real>& availableDeps, NodeIdSpace outputSpace) {
    const std::vector<Attribute> requestedAttributes = expandTopologyBackendAttributes(attributes);
    const std::span<const Attribute> requestedSpan(requestedAttributes);
    const AttributeNames resultNames = AttributeNames::fromList(requestedAttributes);
    std::vector<Real> resultBuffer = makeAttributeValueBuffer<Real>(tree, resultNames);

    DependencyMapT<Real> available = availableDeps;
    if (containsAttribute(requestedSpan, AREA)) {
        const auto it = available.find(AREA);
        if (it != available.end() && isReusableDependencyData(it->second, {AREA})) {
            copyAttributesIntoBuffer(tree, it->second, {AREA}, resultNames, resultBuffer.data());
        } else {
            attributes::computers::AreaComputer::compute(
                AttributeComputeContext<Real>{tree, std::span<Real>(resultBuffer), resultNames, std::span<const Attribute>{}});
        }
        available[AREA] = ComputedAttributeViewT<Real>{&resultNames, resultBuffer.data(), NodeIdSpace::MORPHOLOGICAL_TREE};
    }

    computeTopologyOnlyAttributesIntoResult<Real>(tree, requestedSpan, std::move(available), resultNames, std::span<Real>(resultBuffer), altitude);

    return projectComputedDataToNodeIdSpace(
        tree, altitude, ComputedAttributeData<Real>{std::move(resultNames), std::move(resultBuffer), NodeIdSpace::MORPHOLOGICAL_TREE}, outputSpace);
}

} // namespace mmcfilters::detail
