#pragma once

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

#include <algorithm>
#include <array>
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
 */

inline bool isAttributePipelineAltitudeAttribute(Attribute attribute) noexcept {
    return attributes::registry::isAttributePipelineAltitudeAttribute(attribute);
}

inline bool isAttributePipelineAttribute(Attribute attribute) noexcept {
    return attributes::registry::isPipelineComputed(attribute);
}

inline bool requiresAltitudeForAttributePipeline(Attribute attribute) noexcept {
    return attributes::registry::requiresAltitude(attribute);
}

inline bool isAttributePipelineRequest(const std::vector<AttributeOrGroup>& attributes) {
    return requestContainsOnlyAttributes(
        attributes,
        [](Attribute attribute) {
            return isAttributePipelineAttribute(attribute);
        });
}

inline std::vector<Attribute> expandAttributePipelineAttributes(const std::vector<AttributeOrGroup>& attributes) {
    return expandUniqueAttributeRequest(
        attributes,
        [](Attribute attribute) {
            return isAttributePipelineAttribute(attribute);
        },
        "computeAttributesFromAltitudeSpan received an unsupported attribute.",
        "computeAttributesFromAltitudeSpan received an attribute group containing unsupported attributes.");
}

inline bool containsAltitudeDependentAttribute(std::span<const Attribute> attributes) {
    return std::any_of(
        attributes.begin(),
        attributes.end(),
        [](Attribute attribute) {
            return requiresAltitudeForAttributePipeline(attribute);
        });
}

template<AltitudeValue T>
inline void computeMaxDistAttributes(const MorphologicalTree& tree, std::span<const T> altitude, std::span<float> buffer, const AttributeNames& attrNames) {
    ::mmcfilters::attributes::computers::MaxDistComputer::requireSupportedTreeKind(tree);
    if (!tree.hasAdjacencyRelation()) {
        throw std::invalid_argument("MAX_DIST requires an adjacency relation.");
    }

    ::mmcfilters::attributes::computers::MaxDistComputer maxDistComputer(tree);
    const std::vector<float> maxDist = maxDistComputer.getAttributes(altitude);
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        buffer[attrNames.linearIndex(nodeId, MAX_DIST)] = maxDist[static_cast<std::size_t>(nodeId)];
    }
}

inline ComputedAttributeData materializeAttributesWithoutAltitude(const MorphologicalTree& tree, const std::vector<AttributeOrGroup>& attributes, NodeIdSpace outputSpace){
    tree.requireNotEditing("AttributeComputation::computeAttributes");

    const std::vector<Attribute> requestedAttributes = expandAttributePipelineAttributes(attributes);
    const std::span<const Attribute> requestedSpan(requestedAttributes);
    if (containsAltitudeDependentAttribute(requestedSpan)) {
        throw std::invalid_argument(
            "Requested attributes require altitude data. Use a WeightedMorphologicalTree<T>, "
            "WeightedTreeView, or the overload that receives an altitude buffer.");
    }

    const AttributeNames resultNames = AttributeNames::fromList(requestedAttributes);
    std::vector<float> resultBuffer = makeAttributeValueBuffer(tree, resultNames);

    DependencyMap topologyOnlyDependencies;
    if (containsAttribute(requestedSpan, AREA)) {
        attributes::computers::AreaComputer::computeAreaAttribute(tree, resultBuffer, resultNames);
        topologyOnlyDependencies[AREA] = ComputedAttributeView{
            &resultNames,
            resultBuffer.data(),
            NodeIdSpace::MORPHOLOGICAL_TREE};
    }

    computeTopologyOnlyAttributesIntoResult(tree, requestedSpan, std::move(topologyOnlyDependencies), resultNames, resultBuffer);

    return projectComputedDataToNodeIdSpace(tree, {std::move(resultNames), std::move(resultBuffer), NodeIdSpace::MORPHOLOGICAL_TREE}, outputSpace);
}

template<AltitudeValue T>
inline ComputedAttributeData materializeAttributes(const MorphologicalTree& tree, std::span<const T> altitude, const std::vector<AttributeOrGroup>& attributes, NodeIdSpace outputSpace){
    tree.requireNotEditing("AttributeComputation::computeAttributesFromAltitudeSpan");
    TreeAltitudeAlgorithms::validateAltitudeBufferShape(tree, altitude);

    const std::vector<Attribute> requestedAttributes = expandAttributePipelineAttributes(attributes);
    const AttributeNames resultNames = AttributeNames::fromList(requestedAttributes);
    const int numNodes = tree.getNumInternalNodeSlots();
    std::vector<float> resultBuffer = makeAttributeValueBuffer(tree, resultNames);

    const std::span<const Attribute> requestedSpan(requestedAttributes);
    const bool requestsArea = containsAttribute(requestedSpan, AREA);
    const bool requestsVolume = containsAttribute(requestedSpan, VOLUME);
    const bool requestsRelativeVolume = containsAttribute(requestedSpan, RELATIVE_VOLUME);
    const bool requestsLevel = containsAttribute(requestedSpan, LEVEL);
    const bool requestsMeanLevel = containsAttribute(requestedSpan, MEAN_LEVEL);
    const bool requestsVarianceLevel = containsAttribute(requestedSpan, VARIANCE_LEVEL);
    const bool requestsGrayHeight = containsAttribute(requestedSpan, GRAY_HEIGHT);
    const bool requestsMaxDist = containsAttribute(requestedSpan, MAX_DIST);

    const bool needsAreaDependency = requestsRelativeVolume || requestsMeanLevel || requestsVarianceLevel;
    const bool needsVolumeDependency = requestsMeanLevel || requestsVarianceLevel;

    const AttributeNames areaNames = AttributeNames::fromList({AREA});
    std::vector<float> areaBuffer;
    DependencySource areaSource{};
    if (requestsArea || needsAreaDependency) {
        if (requestsArea) {
            attributes::computers::AreaComputer::computeAreaAttribute(tree, resultBuffer, resultNames);
            areaSource = DependencySource{&resultNames, resultBuffer.data()};
        } else {
            areaBuffer.assign(static_cast<std::size_t>(numNodes), 0.0f);
            attributes::computers::AreaComputer::computeAreaAttribute(tree, areaBuffer, areaNames);
            areaSource = DependencySource{&areaNames, areaBuffer.data()};
        }
    }

    std::vector<float> volumeDependencyBuffer;
    std::unique_ptr<AttributeNames> volumeDependencyNames;
    DependencySource volumeSource{};
    const bool hasHiddenVolumeDependency = needsVolumeDependency && !requestsVolume;
    const bool hasRequestedVolumeFamily = requestsVolume || requestsRelativeVolume;
    if (hasHiddenVolumeDependency || hasRequestedVolumeFamily) {
        std::vector<Attribute> volumeAttributes;
        if (hasHiddenVolumeDependency) {
            volumeAttributes.push_back(VOLUME);
        }
        if (requestsVolume) {
            volumeAttributes.push_back(VOLUME);
        }
        if (requestsRelativeVolume) {
            volumeAttributes.push_back(RELATIVE_VOLUME);
        }
        std::sort(volumeAttributes.begin(), volumeAttributes.end());
        volumeAttributes.erase(std::unique(volumeAttributes.begin(), volumeAttributes.end()), volumeAttributes.end());

        const bool canWriteVolumeDirectly = !hasHiddenVolumeDependency;
        const AttributeNames* volumeNames = &resultNames;
        std::span<float> volumeBuffer(resultBuffer);
        if (!canWriteVolumeDirectly) {
            volumeDependencyNames = std::make_unique<AttributeNames>(AttributeNames::fromList(volumeAttributes));
            volumeDependencyBuffer.assign(static_cast<std::size_t>(numNodes) * static_cast<std::size_t>(volumeDependencyNames->NUM_ATTRIBUTES), 0.0f);
            volumeNames = volumeDependencyNames.get();
            volumeBuffer = std::span<float>(volumeDependencyBuffer);
        }

        std::array<DependencySource, 1> volumeDependencies{{areaSource}};
        const std::span<const DependencySource> volumeDependencySpan = requestsRelativeVolume ? std::span<const DependencySource>(volumeDependencies) : std::span<const DependencySource>();

        ::mmcfilters::attributes::computers::detail::computeVolumeAttributes(tree, altitude, volumeBuffer, *volumeNames, std::span<const Attribute>(volumeAttributes), volumeDependencySpan);

        if (!canWriteVolumeDirectly) {
            if (requestsRelativeVolume) {
                copyAttributeValuesBetweenLayouts(tree, *volumeNames, volumeDependencyBuffer, resultNames, resultBuffer, RELATIVE_VOLUME);
            }
            volumeSource = DependencySource{volumeNames, volumeDependencyBuffer.data()};
        } else if (needsVolumeDependency) {
            volumeSource = DependencySource{&resultNames, resultBuffer.data()};
        }
    }

    std::vector<Attribute> grayAttributes;
    if (requestsLevel) {
        grayAttributes.push_back(LEVEL);
    }
    if (requestsGrayHeight) {
        grayAttributes.push_back(GRAY_HEIGHT);
    }
    if (requestsMeanLevel) {
        grayAttributes.push_back(MEAN_LEVEL);
    }
    if (requestsVarianceLevel) {
        grayAttributes.push_back(VARIANCE_LEVEL);
    }
    if (!grayAttributes.empty()) {
        std::sort(grayAttributes.begin(), grayAttributes.end());
        std::array<DependencySource, 2> grayDependencies{{volumeSource, areaSource}};
        const std::span<const DependencySource> grayDependencySpan = needsVolumeDependency ? std::span<const DependencySource>(grayDependencies) : std::span<const DependencySource>();

        ::mmcfilters::attributes::computers::detail::computeGrayLevelStatsAttributes(tree, altitude, resultBuffer, resultNames, std::span<const Attribute>(grayAttributes), grayDependencySpan);
    }

    if (requestsMaxDist) {
        computeMaxDistAttributes(tree, altitude, resultBuffer, resultNames);
    }

    DependencyMap topologyOnlyDependencies;
    if (areaSource.attrNames != nullptr && areaSource.buffer != nullptr) {
        topologyOnlyDependencies[AREA] = ComputedAttributeView{
            areaSource.attrNames,
            areaSource.buffer,
            NodeIdSpace::MORPHOLOGICAL_TREE};
    }
    computeTopologyOnlyAttributesIntoResult(tree, requestedSpan, std::move(topologyOnlyDependencies), resultNames, resultBuffer, altitude);

    return projectComputedDataToNodeIdSpace(tree, altitude, {std::move(resultNames), std::move(resultBuffer), NodeIdSpace::MORPHOLOGICAL_TREE}, outputSpace);
}

} // namespace mmcfilters::detail
