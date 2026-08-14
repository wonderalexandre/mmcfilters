#pragma once

/**
 * @file NodeAttributeSampleMaterialization.hpp
 * @internal
 * @brief Internal helpers for altitude-based node-attribute sampling.
 *
 * Public callers should use `AttributeComputation::computeSampledNodeAttribute(...)`.
 */

#include "../../trees/detail/TreeAttributeSamplingNeighborhood.hpp"
#include "AttributeProjection.hpp"
#include "../AttributeResultTypes.hpp"
#include "../../trees/MorphologicalTree.hpp"
#include "../../trees/TreeAltitudeAlgorithms.hpp"
#include "../../utils/Altitude.hpp"
#include "../../utils/Common.hpp"
#include "../../utils/Contract.hpp"

#include <cmath>
#include <concepts>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mmcfilters::detail {

/**
 * @brief Rejects a negative ancestor/descendant sampling radius.
 * @param samplingRadius Requested radius.
 * @param context Diagnostic operation name.
 */
inline void validateSamplingRadius(int samplingRadius, const char* context) {
    MMCFILTERS_CONTRACT_REQUIRE(samplingRadius >= 0,
                                throw std::invalid_argument(std::string(context) + " requires a non-negative sampling radius."));
}

/**
 * @brief Materializes one scalar node attribute at altitude-based sample positions.
 *
 * The output contains the current-node value at offset zero, ancestor samples
 * at negative offsets, and policy-selected representative-descendant samples
 * at positive offsets. Missing positions are finalized independently of the
 * values stored by valid samples, so a valid `NaN` attribute value is not
 * mistaken for an unavailable sample.
 * @tparam Real Public attribute storage type.
 * @tparam T Node-altitude type.
 * @param tree Tree topology.
 * @param altitude Dense node-altitude buffer.
 * @param base Base scalar node attribute.
 * @param attribute Attribute sampled from `base`.
 * @param altitudeStep Positive distance between samples.
 * @param samplingRadius Number of ancestor/descendant positions.
 * @param samplingPolicy Representative-descendant policy.
 * @param missingSamplePolicy Missing-position materialization policy.
 * @param outputSpace Requested result node-id space.
 * @return Sample layout and materialized node values.
 */
template <std::floating_point Real, AltitudeValue T>
inline SampledNodeAttributeData<Real>
materializeNodeAttributeSamples(const MorphologicalTree& tree, std::span<const T> altitude, ComputedAttributeData<Real> base,
                                Attribute attribute, AltitudeDifference<T> altitudeStep, int samplingRadius,
                                NodeAttributeSamplingPolicy samplingPolicy, MissingNodeAttributeSamplePolicy missingSamplePolicy,
                                NodeIdSpace outputSpace) {
    validateSamplingRadius(samplingRadius, "computeSampledNodeAttribute");
    validatePositiveAltitudeStep<T>(altitudeStep, "computeSampledNodeAttribute");
    TreeAltitudeAlgorithms::validateNodeAltitudeBufferShape(tree, altitude);
    validateNodeAttributeSamplingPolicy(samplingPolicy);
    switch (missingSamplePolicy) {
    case MissingNodeAttributeSamplePolicy::RepeatNearest:
    case MissingNodeAttributeSamplePolicy::NotANumber:
    case MissingNodeAttributeSamplePolicy::Zero:
        break;
    default:
        throw std::invalid_argument("Unknown missing node-attribute sample policy.");
    }
    if (base.nodeIdSpace != NodeIdSpace::MorphologicalTree) {
        throw std::logic_error("Node-attribute sampling requires base attributes in internal node-id space.");
    }

    const AttributeNames& baseLayout = base.first;
    const std::vector<Real>& baseValues = base.second;
    const int numSlots = tree.numInternalNodeSlots();
    NodeAttributeSampleLayout sampleLayout = NodeAttributeSampleLayout::create(samplingRadius, {attribute});
    const std::size_t outputSize = static_cast<std::size_t>(numSlots) * static_cast<std::size_t>(sampleLayout.NUM_ATTRIBUTES);
    std::vector<Real> sampledValues(outputSize, std::numeric_limits<Real>::quiet_NaN());
    std::vector<std::uint8_t> samplePresent(outputSize, std::uint8_t{0});
    const auto aliveNodes = tree.aliveNodeIds();

    NodeSupportSamplingMetadata supportMetadata;
    NodeAttributeSamplingNeighborhood neighborhood;
    if (samplingRadius > 0) {
        validateGlobalMonotoneAltitudeOrder(tree, "computeSampledNodeAttribute");
        supportMetadata = computeNodeSupportSamplingMetadata(tree);
        neighborhood = {
            std::vector<NodeId>(static_cast<std::size_t>(numSlots), InvalidNode),
            std::vector<NodeId>(static_cast<std::size_t>(numSlots), InvalidNode),
        };
    }

    for (NodeId nodeId : aliveNodes) {
        const int outputIndex = sampleLayout.linearIndex(nodeId, attribute, 0);
        const int baseIndex = baseLayout.linearIndex(nodeId, attribute);
        sampledValues[static_cast<std::size_t>(outputIndex)] = baseValues[static_cast<std::size_t>(baseIndex)];
        samplePresent[static_cast<std::size_t>(outputIndex)] = std::uint8_t{1};
    }

    for (int sampleDistance = 1; sampleDistance <= samplingRadius; ++sampleDistance) {
        const AltitudeDifference<T> altitudeDistance = altitudeStep * static_cast<AltitudeDifference<T>>(sampleDistance);
        validatePositiveAltitudeStep<T>(altitudeDistance, "computeSampledNodeAttribute");
        kernel::fillNodeAttributeSamplingNeighborhood(tree, aliveNodes, altitude, altitudeDistance, supportMetadata, neighborhood);

        for (NodeId nodeId : aliveNodes) {
            const NodeId ancestor = neighborhood.ancestors[static_cast<std::size_t>(nodeId)];
            if (ancestor != InvalidNode) {
                const int outputIndex = sampleLayout.linearIndex(nodeId, attribute, -sampleDistance);
                const int baseIndex = baseLayout.linearIndex(ancestor, attribute);
                sampledValues[static_cast<std::size_t>(outputIndex)] = baseValues[static_cast<std::size_t>(baseIndex)];
                samplePresent[static_cast<std::size_t>(outputIndex)] = std::uint8_t{1};
            }

            const NodeId descendant = neighborhood.representativeDescendants[static_cast<std::size_t>(nodeId)];
            if (descendant != InvalidNode) {
                const int outputIndex = sampleLayout.linearIndex(nodeId, attribute, sampleDistance);
                const int baseIndex = baseLayout.linearIndex(descendant, attribute);
                sampledValues[static_cast<std::size_t>(outputIndex)] = baseValues[static_cast<std::size_t>(baseIndex)];
                samplePresent[static_cast<std::size_t>(outputIndex)] = std::uint8_t{1};
            }
        }
    }

    for (NodeId nodeId : aliveNodes) {
        for (const int direction : {-1, 1}) {
            for (int sampleDistance = 1; sampleDistance <= samplingRadius; ++sampleDistance) {
                const int sampleOffset = direction * sampleDistance;
                const int outputIndex = sampleLayout.linearIndex(nodeId, attribute, sampleOffset);
                const std::size_t outputPosition = static_cast<std::size_t>(outputIndex);
                if (samplePresent[outputPosition] != 0) {
                    continue;
                }

                switch (missingSamplePolicy) {
                case MissingNodeAttributeSamplePolicy::RepeatNearest: {
                    const int nearestOffset = direction * (sampleDistance - 1);
                    const int nearestIndex = sampleLayout.linearIndex(nodeId, attribute, nearestOffset);
                    sampledValues[outputPosition] = sampledValues[static_cast<std::size_t>(nearestIndex)];
                    break;
                }
                case MissingNodeAttributeSamplePolicy::NotANumber:
                    sampledValues[outputPosition] = std::numeric_limits<Real>::quiet_NaN();
                    break;
                case MissingNodeAttributeSamplePolicy::Zero:
                    sampledValues[outputPosition] = Real{0};
                    break;
                default:
                    throw std::invalid_argument("Unknown missing node-attribute sample policy.");
                }
            }
        }
    }

    return projectComputedDataToNodeIdSpace(
        tree, altitude,
        SampledNodeAttributeData<Real>{std::move(sampleLayout), std::move(sampledValues), NodeIdSpace::MorphologicalTree}, outputSpace);
}

} // namespace mmcfilters::detail
