#pragma once

#include "AttributeDependencyCache.hpp"
#include "AttributeProjection.hpp"
#include "AttributeRequestUtils.hpp"
#include "../AttributeNames.hpp"
#include "../AttributeRegistry.hpp"
#include "../AttributeResultTypes.hpp"
#include "../computers/AreaComputer.hpp"
#include "../computers/BoundingBoxComputer.hpp"
#include "../computers/BitquadAttributeComputer.hpp"
#include "../computers/ContourSideAttributeComputer.hpp"
#include "../computers/detail/BitquadLocalEventComputation.hpp"
#include "../computers/MomentBasedAttributeComputer.hpp"
#include "../computers/TreeTopologyComputer.hpp"
#include "../../trees/MorphologicalTree.hpp"
#include "../../utils/Altitude.hpp"
#include "../../utils/Common.hpp"

#include <algorithm>
#include <array>
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
 */

inline ComputedAttributeData materializeTopologyAttributeRequest(const MorphologicalTree& tree, const std::vector<AttributeOrGroup>& attributes, const DependencyMap& availableDeps, NodeIdSpace outputSpace);

template<AltitudeValue T>
inline ComputedAttributeData materializeTopologyAttributeRequest(const MorphologicalTree& tree, std::span<const T> altitude, const std::vector<AttributeOrGroup>& attributes, const DependencyMap& availableDeps, NodeIdSpace outputSpace);

inline bool isTopologyOnlyAttribute(Attribute attribute) noexcept {
    return attributes::registry::isTopologyOnly(attribute);
}

inline bool isTopologyBackendAttribute(Attribute attribute) noexcept {
    return attribute == AREA || isTopologyOnlyAttribute(attribute);
}

inline bool isBoundingBoxBackendAttribute(Attribute attribute) noexcept {
    switch (attribute) {
        case BOX_WIDTH:
        case BOX_HEIGHT:
        case DIAGONAL_LENGTH:
        case RECTANGULARITY:
        case RATIO_WH:
        case BOX_COL_MIN:
        case BOX_COL_MAX:
        case BOX_ROW_MIN:
        case BOX_ROW_MAX:
            return true;
        default:
            return false;
    }
}

inline bool isTreeTopologyBackendAttribute(Attribute attribute) noexcept {
    switch (attribute) {
        case HEIGHT_NODE:
        case DEPTH_NODE:
        case IS_LEAF_NODE:
        case IS_ROOT_NODE:
        case NUM_CHILDREN_NODE:
        case NUM_SIBLINGS_NODE:
        case NUM_DESCENDANTS_NODE:
        case NUM_LEAF_DESCENDANTS_NODE:
        case LEAF_RATIO_NODE:
        case BALANCE_NODE:
        case AVG_CHILD_HEIGHT_NODE:
            return true;
        default:
            return false;
    }
}

inline bool isCentralMomentBackendAttribute(Attribute attribute) noexcept {
    switch (attribute) {
        case CENTRAL_MOMENT_20:
        case CENTRAL_MOMENT_02:
        case CENTRAL_MOMENT_11:
        case CENTRAL_MOMENT_30:
        case CENTRAL_MOMENT_03:
        case CENTRAL_MOMENT_21:
        case CENTRAL_MOMENT_12:
            return true;
        default:
            return false;
    }
}

inline bool isHuMomentBackendAttribute(Attribute attribute) noexcept {
    switch (attribute) {
        case HU_MOMENT_1:
        case HU_MOMENT_2:
        case HU_MOMENT_3:
        case HU_MOMENT_4:
        case HU_MOMENT_5:
        case HU_MOMENT_6:
        case HU_MOMENT_7:
            return true;
        default:
            return false;
    }
}

inline bool isMomentBasedBackendAttribute(Attribute attribute) noexcept {
    switch (attribute) {
        case COMPACTNESS:
        case ECCENTRICITY:
        case LENGTH_MAJOR_AXIS:
        case LENGTH_MINOR_AXIS:
        case AXIS_ORIENTATION:
        case INERTIA:
        case CIRCULARITY:
            return true;
        default:
            return false;
    }
}

inline bool isBitquadBackendAttribute(Attribute attribute) noexcept {
    switch (attribute) {
        case BITQUADS_AREA:
        case BITQUADS_NUMBER_EULER:
        case BITQUADS_NUMBER_HOLES:
        case BITQUADS_PERIMETER:
        case BITQUADS_PERIMETER_CONTINUOUS:
        case BITQUADS_CIRCULARITY:
        case BITQUADS_PERIMETER_AVERAGE:
        case BITQUADS_LENGTH_AVERAGE:
        case BITQUADS_WIDTH_AVERAGE:
            return true;
        default:
            return false;
    }
}

inline bool isContourBackendAttribute(Attribute attribute) noexcept {
    switch (attribute) {
        case CONTOUR_PIXELS:
        case CONTOUR_PERIMETER:
        case CONTOUR_SIDE_NORTH:
        case CONTOUR_SIDE_WEST:
        case CONTOUR_SIDE_EAST:
        case CONTOUR_SIDE_SOUTH:
            return true;
        default:
            return false;
    }
}

inline std::vector<Attribute> expandTopologyBackendAttributes(const std::vector<AttributeOrGroup>& attributes){
    return expandUniqueAttributeRequest(
        attributes,
        [](Attribute attribute) {
            return isTopologyBackendAttribute(attribute);
        },
        "TopologyAttributeBackend received an attribute that requires altitude.",
        "TopologyAttributeBackend received an attribute group containing altitude-dependent attributes.");
}

struct NoTopologyBackendAltitude {};

inline void computeBitquadBackendAttributesIntoResult(const MorphologicalTree& tree, std::span<const Attribute> bitquadAttributes, const AttributeNames& resultNames, std::span<float> resultBuffer, NoTopologyBackendAltitude) {
    ::mmcfilters::attributes::computers::BitquadAttributeComputer computer;
    computer.compute(tree, AttributeAltitudeView{}, resultBuffer, resultNames, bitquadAttributes, {});
}

template<AltitudeValue T>
inline void computeBitquadBackendAttributesIntoResult(const MorphologicalTree& tree, std::span<const Attribute> bitquadAttributes, const AttributeNames& resultNames, std::span<float> resultBuffer, std::span<const T> altitude) {
    using BitquadLocalEventComputation = ::mmcfilters::attributes::computers::detail::BitquadLocalEventComputation;
    const auto familyDeltas = BitquadLocalEventComputation::computeBitquadFamilyDeltas(tree);
    const auto familyCounts = BitquadLocalEventComputation::aggregateBitquadFamilyDeltas(tree, familyDeltas);
    ::mmcfilters::attributes::computers::BitquadAttributeComputer::materializeAttributesFromBitquadFamilyCounts(tree, altitude, familyCounts, resultBuffer, resultNames, bitquadAttributes);
}

template<class BitquadAltitude>
inline void computeTopologyOnlyAttributesIntoResultImpl(const MorphologicalTree& tree, std::span<const Attribute> requestedAttributes, DependencyMap available, const AttributeNames& resultNames, std::span<float> resultBuffer, BitquadAltitude bitquadAltitude) {
    OwnedComputedResults ownedResults;

    auto ensureInternalAreaDependency = [&]() -> DependencySource {
        const auto it = available.find(AREA);
        if (it != available.end() && isReusableDependencyData(it->second, {AREA})) {
            return it->second.dependencySource();
        }

        AttributeNames areaNames = AttributeNames::fromList({AREA});
        std::vector<float> areaBuffer = makeAttributeValueBuffer(tree, areaNames);
        attributes::computers::AreaComputer::computeAreaAttribute(tree, areaBuffer, areaNames);
        stashComputedAttributes(ownedResults, available, ComputedAttributeData(std::move(areaNames), std::move(areaBuffer)));
        return available.at(AREA).dependencySource();
    };

    auto ensureInternalCentralMomentDependency =
        [&](const std::vector<Attribute>& requiredAttributes) -> DependencySource {
            for (const Attribute attribute : requiredAttributes) {
                const auto it = available.find(attribute);
                if (it != available.end() && isReusableDependencyData(it->second, requiredAttributes)) {
                    return it->second.dependencySource();
                }
            }

            AttributeNames momentNames = AttributeNames::fromList(requiredAttributes);
            std::vector<float> momentBuffer = makeAttributeValueBuffer(tree, momentNames);
            attributes::computers::CentralMomentsComputer computer;
            computer.compute(tree, AttributeAltitudeView{}, momentBuffer, momentNames, std::span<const Attribute>(requiredAttributes), {});
            stashComputedAttributes(ownedResults, available, ComputedAttributeData(std::move(momentNames), std::move(momentBuffer)));
            return available.at(requiredAttributes.front()).dependencySource();
        };

    std::vector<Attribute> boundingBoxAttributes;
    std::vector<Attribute> treeTopologyAttributes;
    std::vector<Attribute> centralMomentAttributes;
    std::vector<Attribute> huMomentAttributes;
    std::vector<Attribute> momentBasedAttributes;
    std::vector<Attribute> bitquadAttributes;
    std::vector<Attribute> contourAttributes;
    std::vector<AttributeOrGroup> topologyOnlyRequests;

    for (const Attribute attribute : requestedAttributes) {
        if (!isTopologyOnlyAttribute(attribute)) {
            continue;
        }
        if (isBoundingBoxBackendAttribute(attribute)) {
            boundingBoxAttributes.push_back(attribute);
            continue;
        }
        if (isTreeTopologyBackendAttribute(attribute)) {
            treeTopologyAttributes.push_back(attribute);
            continue;
        }
        if (isCentralMomentBackendAttribute(attribute)) {
            centralMomentAttributes.push_back(attribute);
            continue;
        }
        if (isHuMomentBackendAttribute(attribute)) {
            huMomentAttributes.push_back(attribute);
            continue;
        }
        if (isMomentBasedBackendAttribute(attribute)) {
            momentBasedAttributes.push_back(attribute);
            continue;
        }
        if (isBitquadBackendAttribute(attribute)) {
            bitquadAttributes.push_back(attribute);
            continue;
        }
        if (isContourBackendAttribute(attribute)) {
            contourAttributes.push_back(attribute);
            continue;
        }
        topologyOnlyRequests.emplace_back(attribute);
    }

    if (!boundingBoxAttributes.empty()) {
        const bool needsAreaDependency = containsAttribute(boundingBoxAttributes, RECTANGULARITY);
        std::array<DependencySource, 1> areaDependency{{}};
        std::span<const DependencySource> dependencySources;
        if (needsAreaDependency) {
            areaDependency[0] = ensureInternalAreaDependency();
            dependencySources = std::span<const DependencySource>(areaDependency);
        }

        attributes::computers::BoundingBoxComputer computer;
        computer.compute(tree, AttributeAltitudeView{}, resultBuffer, resultNames, std::span<const Attribute>(boundingBoxAttributes), dependencySources);

        for (const Attribute attribute : boundingBoxAttributes) {
            available[attribute] = ComputedAttributeView{
                &resultNames,
                resultBuffer.data(),
                NodeIdSpace::MORPHOLOGICAL_TREE};
        }
    }

    if (!treeTopologyAttributes.empty()) {
        attributes::computers::TreeTopologyComputer computer;
        computer.compute(tree, AttributeAltitudeView{}, resultBuffer, resultNames, std::span<const Attribute>(treeTopologyAttributes), {});

        for (const Attribute attribute : treeTopologyAttributes) {
            available[attribute] = ComputedAttributeView{
                &resultNames,
                resultBuffer.data(),
                NodeIdSpace::MORPHOLOGICAL_TREE};
        }
    }

    if (!centralMomentAttributes.empty()) {
        attributes::computers::CentralMomentsComputer computer;
        computer.compute(tree, AttributeAltitudeView{}, resultBuffer, resultNames, std::span<const Attribute>(centralMomentAttributes), {});

        for (const Attribute attribute : centralMomentAttributes) {
            available[attribute] = ComputedAttributeView{
                &resultNames,
                resultBuffer.data(),
                NodeIdSpace::MORPHOLOGICAL_TREE};
        }
    }

    if (!huMomentAttributes.empty()) {
        static const std::vector<Attribute> requiredCentralMoments{
            CENTRAL_MOMENT_20,
            CENTRAL_MOMENT_02,
            CENTRAL_MOMENT_11,
            CENTRAL_MOMENT_30,
            CENTRAL_MOMENT_03,
            CENTRAL_MOMENT_21,
            CENTRAL_MOMENT_12};
        std::array<DependencySource, 2> dependencies{{
            ensureInternalCentralMomentDependency(requiredCentralMoments),
            ensureInternalAreaDependency()}};

        attributes::computers::HuMomentsComputer computer;
        computer.compute(tree, AttributeAltitudeView{}, resultBuffer, resultNames, std::span<const Attribute>(huMomentAttributes), std::span<const DependencySource>(dependencies));

        for (const Attribute attribute : huMomentAttributes) {
            available[attribute] = ComputedAttributeView{
                &resultNames,
                resultBuffer.data(),
                NodeIdSpace::MORPHOLOGICAL_TREE};
        }
    }

    if (!momentBasedAttributes.empty()) {
        static const std::vector<Attribute> requiredCentralMoments{
            CENTRAL_MOMENT_20,
            CENTRAL_MOMENT_02,
            CENTRAL_MOMENT_11};
        std::array<DependencySource, 2> dependencies{{
            ensureInternalCentralMomentDependency(requiredCentralMoments),
            ensureInternalAreaDependency()}};

        attributes::computers::MomentBasedAttributeComputer computer;
        computer.compute(tree, AttributeAltitudeView{}, resultBuffer, resultNames, std::span<const Attribute>(momentBasedAttributes), std::span<const DependencySource>(dependencies));

        for (const Attribute attribute : momentBasedAttributes) {
            available[attribute] = ComputedAttributeView{
                &resultNames,
                resultBuffer.data(),
                NodeIdSpace::MORPHOLOGICAL_TREE};
        }
    }

    if (!bitquadAttributes.empty()) {
        computeBitquadBackendAttributesIntoResult(tree, std::span<const Attribute>(bitquadAttributes), resultNames, resultBuffer, bitquadAltitude);

        for (const Attribute attribute : bitquadAttributes) {
            available[attribute] = ComputedAttributeView{
                &resultNames,
                resultBuffer.data(),
                NodeIdSpace::MORPHOLOGICAL_TREE};
        }
    }

    if (!contourAttributes.empty()) {
        ::mmcfilters::attributes::computers::ContourSideAttributeComputer computer;
        computer.compute(tree, AttributeAltitudeView{}, resultBuffer, resultNames, std::span<const Attribute>(contourAttributes), {});

        for (const Attribute attribute : contourAttributes) {
            available[attribute] = ComputedAttributeView{
                &resultNames,
                resultBuffer.data(),
                NodeIdSpace::MORPHOLOGICAL_TREE};
        }
    }

    if (!topologyOnlyRequests.empty()) {
        throw std::runtime_error("TopologyAttributeBackend received a topology-only attribute without an explicit pipeline family.");
    }

    for (const Attribute attribute : requestedAttributes) {
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

inline void computeTopologyOnlyAttributesIntoResult(const MorphologicalTree& tree, std::span<const Attribute> requestedAttributes, DependencyMap available, const AttributeNames& resultNames, std::span<float> resultBuffer) {
    computeTopologyOnlyAttributesIntoResultImpl(tree, requestedAttributes, std::move(available), resultNames, resultBuffer, NoTopologyBackendAltitude{});
}

template<AltitudeValue T>
inline void computeTopologyOnlyAttributesIntoResult( const MorphologicalTree& tree, std::span<const Attribute> requestedAttributes, DependencyMap available, const AttributeNames& resultNames, std::span<float> resultBuffer, std::span<const T> altitude) {
    computeTopologyOnlyAttributesIntoResultImpl(tree, requestedAttributes, std::move(available), resultNames, resultBuffer, altitude);
}

inline ComputedAttributeData materializeTopologyAttributeRequest(const MorphologicalTree& tree, const std::vector<AttributeOrGroup>& attributes, const DependencyMap& availableDeps, NodeIdSpace outputSpace) {
    const std::vector<Attribute> requestedAttributes = expandTopologyBackendAttributes(attributes);
    const std::span<const Attribute> requestedSpan(requestedAttributes);
    const AttributeNames resultNames = AttributeNames::fromList(requestedAttributes);
    std::vector<float> resultBuffer = makeAttributeValueBuffer(tree, resultNames);

    DependencyMap available = availableDeps;
    if (containsAttribute(requestedSpan, AREA)) {
        const auto it = available.find(AREA);
        if (it != available.end() && isReusableDependencyData(it->second, {AREA})) {
            copyAttributesIntoBuffer(tree, it->second, {AREA}, resultNames, resultBuffer.data());
        } else {
            attributes::computers::AreaComputer::computeAreaAttribute(tree, resultBuffer, resultNames);
        }
        available[AREA] = ComputedAttributeView{
            &resultNames,
            resultBuffer.data(),
            NodeIdSpace::MORPHOLOGICAL_TREE};
    }

    computeTopologyOnlyAttributesIntoResult(tree, requestedSpan, std::move(available), resultNames, resultBuffer);

    return projectComputedDataToNodeIdSpace(tree, {std::move(resultNames), std::move(resultBuffer), NodeIdSpace::MORPHOLOGICAL_TREE}, outputSpace);
}

template<AltitudeValue T>
inline ComputedAttributeData materializeTopologyAttributeRequest(const MorphologicalTree& tree, std::span<const T> altitude, const std::vector<AttributeOrGroup>& attributes, const DependencyMap& availableDeps, NodeIdSpace outputSpace) {
    const std::vector<Attribute> requestedAttributes = expandTopologyBackendAttributes(attributes);
    const std::span<const Attribute> requestedSpan(requestedAttributes);
    const AttributeNames resultNames = AttributeNames::fromList(requestedAttributes);
    std::vector<float> resultBuffer = makeAttributeValueBuffer(tree, resultNames);

    DependencyMap available = availableDeps;
    if (containsAttribute(requestedSpan, AREA)) {
        const auto it = available.find(AREA);
        if (it != available.end() && isReusableDependencyData(it->second, {AREA})) {
            copyAttributesIntoBuffer(tree, it->second, {AREA}, resultNames, resultBuffer.data());
        } else {
            attributes::computers::AreaComputer::computeAreaAttribute(tree, resultBuffer, resultNames);
        }
        available[AREA] = ComputedAttributeView{
            &resultNames,
            resultBuffer.data(),
            NodeIdSpace::MORPHOLOGICAL_TREE};
    }

    computeTopologyOnlyAttributesIntoResult(tree, requestedSpan, std::move(available), resultNames, resultBuffer, altitude);

    return projectComputedDataToNodeIdSpace(tree, altitude, {std::move(resultNames), std::move(resultBuffer), NodeIdSpace::MORPHOLOGICAL_TREE}, outputSpace);
}

} // namespace mmcfilters::detail
