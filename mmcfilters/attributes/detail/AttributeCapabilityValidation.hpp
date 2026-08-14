#pragma once

#include "../AttributeRegistry.hpp"
#include "../../trees/MorphologicalTree.hpp"
#include "../../trees/detail/MorphologicalTreeConstructionContextQueries.hpp"

#include <span>
#include <stdexcept>
#include <string>

namespace mmcfilters::detail {

/**
 * @brief Tests whether bitquad projection needs non-root shape polarity.
 *
 * @param tree Tree topology.
 * @return True if lower- and upper-shape connectivity differ; otherwise false.
 */
inline bool bitquadProjectionNeedsShapePolarity(const MorphologicalTree& tree) {
    const auto adjacencies = currentBitquadProjectionAdjacencies(tree);
    return adjacencies && adjacencies->minAdjacency.is4connectivity() != adjacencies->maxAdjacency.is4connectivity();
}

/**
 * @brief Validates expanded scalar-attribute requirements before scheduling.
 *
 * @param tree Tree topology.
 * @param attributes Attributes requested by the operation.
 * @param altitudeAvailable Altitude or level.
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
            throw std::invalid_argument(prefix + "a regular 2D pixel domain is required.");
        }
        if (requirements.monotoneAltitudeOrder && tree.nodeAltitudeOrder() == NodeAltitudeOrder::Unconstrained) {
            throw std::invalid_argument(prefix + "a globally monotone altitude order is required.");
        }
        const RegularGridAdjacency2D* constructionAdjacency = ::mmcfilters::detail::constructionAdjacency(tree);
        const auto projectionAdjacencies = currentBitquadProjectionAdjacencies(tree);
        if (requirements.adjacency == attributes::registry::AttributeAdjacencyRequirement::Uniform && constructionAdjacency == nullptr) {
            throw std::invalid_argument(prefix + "a shared or saturated construction adjacency is required.");
        }
        if (requirements.adjacency == attributes::registry::AttributeAdjacencyRequirement::UniformOrDirectional && constructionAdjacency == nullptr &&
            !projectionAdjacencies) {
            throw std::invalid_argument(prefix + "a construction adjacency or topographic convention is required.");
        }
        if (requirements.altitudeForDirectionalAdjacency && bitquadProjectionNeedsShapePolarity(tree) && !altitudeAvailable) {
            throw std::invalid_argument(prefix + "node altitudes are required to derive lower/upper shape polarity for bitquad connectivity selection.");
        }
        if (requirements.canonical4Or8Adjacency) {
            if ((constructionAdjacency != nullptr && !constructionAdjacency->isCanonical4Or8Connectivity()) ||
                (projectionAdjacencies && (!projectionAdjacencies->minAdjacency.isCanonical4Or8Connectivity() ||
                                           !projectionAdjacencies->maxAdjacency.isCanonical4Or8Connectivity()))) {
                throw std::invalid_argument(
                    prefix + "canonical 4- or 8-connectivity is required; the retained adjacency is unsupported by bitquad formulas.");
            }
        }
    }
}

} // namespace mmcfilters::detail
