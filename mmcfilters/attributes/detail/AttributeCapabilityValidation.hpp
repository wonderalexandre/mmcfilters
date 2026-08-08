#pragma once

#include "../AttributeRegistry.hpp"
#include "../../trees/MorphologicalTree.hpp"

#include <span>
#include <stdexcept>
#include <string>

namespace mmcfilters::detail {

/**
 * @brief Tests whether unequal directional adjacency needs altitude selection.
 *
 * @param tree Tree topology used by the operation.
 * @return True if unequal directional adjacency needs altitude selection; otherwise false.
 */
inline bool directionalAdjacencyNeedsAltitude(const MorphologicalTree& tree) noexcept {
    const RegularGridAdjacency2D* decreasing = tree.getDecreasingGridAdjacency2D();
    const RegularGridAdjacency2D* increasing = tree.getIncreasingGridAdjacency2D();
    return decreasing != nullptr && increasing != nullptr && decreasing->is4connectivity() != increasing->is4connectivity();
}

/**
 * @brief Validates expanded scalar-attribute requirements before scheduling.
 *
 * @param tree Tree topology used by the operation.
 * @param attributes Attributes requested by the operation.
 * @param altitudeAvailable Altitude or level represented by `altitudeAvailable`.
 * @param context Operation context or diagnostic label.
 */
inline void validateAttributeCapabilities(const MorphologicalTree& tree, std::span<const Attribute> attributes, bool altitudeAvailable, const char* context) {
    for (Attribute attribute : attributes) {
        const auto requirements = attributes::registry::capabilityRequirements(attribute);
        const std::string prefix = std::string(context) + " cannot compute " + std::string(attributes::registry::name(attribute)) + ": ";

        if (requirements.altitude && !altitudeAvailable) {
            throw std::invalid_argument(prefix + "a dense node-altitude buffer is required.");
        }
        if (requirements.gridDomain2D && !tree.hasGridDomain2D()) {
            throw std::invalid_argument(prefix + "a regular 2D proper-part domain is required.");
        }
        if (requirements.monotoneAltitudeOrder && tree.getAltitudeOrder() == AltitudeOrder::UNCONSTRAINED) {
            throw std::invalid_argument(prefix + "a globally monotone altitude order is required.");
        }
        if (requirements.adjacency == attributes::registry::AttributeAdjacencyRequirement::UNIFORM && !tree.hasUniformGridAdjacency2D()) {
            throw std::invalid_argument(prefix + "uniform adjacency is required.");
        }
        if (requirements.adjacency == attributes::registry::AttributeAdjacencyRequirement::UNIFORM_OR_DIRECTIONAL && !tree.hasUniformGridAdjacency2D() &&
            !tree.hasDirectionalGridAdjacency2D()) {
            throw std::invalid_argument(prefix + "uniform or directional adjacency is required.");
        }
        if (requirements.altitudeForDirectionalAdjacency && tree.hasDirectionalGridAdjacency2D() && directionalAdjacencyNeedsAltitude(tree) &&
            !altitudeAvailable) {
            throw std::invalid_argument(prefix + "altitude is required to select between unequal directional adjacencies.");
        }
        if (requirements.canonical4Or8Adjacency) {
            const RegularGridAdjacency2D* uniform = tree.getUniformGridAdjacency2D();
            const RegularGridAdjacency2D* decreasing = tree.getDecreasingGridAdjacency2D();
            const RegularGridAdjacency2D* increasing = tree.getIncreasingGridAdjacency2D();
            if ((uniform != nullptr && !uniform->isCanonical4Or8Connectivity()) || (decreasing != nullptr && !decreasing->isCanonical4Or8Connectivity()) ||
                (increasing != nullptr && !increasing->isCanonical4Or8Connectivity())) {
                throw std::invalid_argument(
                    prefix + "canonical 4- or 8-connectivity is required; the stored structuring-element adjacency is unsupported by BitQuad formulas.");
            }
        }
    }
}

} // namespace mmcfilters::detail
