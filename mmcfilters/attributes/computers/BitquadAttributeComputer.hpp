#pragma once

#include "AttributeComputerDomain.hpp"
#include "../detail/AttributeKernelSupport.hpp"
#include "detail/BitquadAttributeMaterialization.hpp"
#include "detail/BitquadLocalEventComputation.hpp"

#include <array>
#include <numbers>
#include <span>
#include <stdexcept>
#include <string_view>

namespace mmcfilters::attributes::computers {

/**
 * @brief Bitquad scalar computer backed by local events.
 *
 * @details
 * This class obtains bitquad family deltas from its local-event computation,
 * aggregates them into per-node family counts, and materializes the public
 * `BITQUADS_*` descriptors. The public attribute pipeline/topology backend
 * invokes this computer for bitquad requests.
 *
 * The computer exposes the standard scalar bitquad descriptors. For component
 * trees, connectivity comes from the tree's adjacency relation. For Tree of
 * Shapes inputs, scalar projection additionally requires node altitudes so the
 * node can be classified as
 * min-tree or max-tree relative to its parent and the appropriate auxiliary
 * adjacency can be selected.
 */
class BitquadAttributeComputer {
public:
    /// Family name used in dependency-plan diagnostics.
    static constexpr std::string_view familyName = "bitquad";

    /// Execution domain required by the computer.
    static constexpr AttributeComputerDomain domain = AttributeComputerDomain::Topology;

    /**
     * @brief Canonical list of scalar bitquad attributes materialized by this computer.
     */
    inline static constexpr std::array<Attribute, 9> producedAttributes{
        BITQUADS_AREA,
        BITQUADS_NUMBER_EULER,
        BITQUADS_NUMBER_HOLES,
        BITQUADS_PERIMETER,
        BITQUADS_PERIMETER_CONTINUOUS,
        BITQUADS_CIRCULARITY,
        BITQUADS_PERIMETER_AVERAGE,
        BITQUADS_LENGTH_AVERAGE,
        BITQUADS_WIDTH_AVERAGE};

    /**
     * @brief Computes requested scalar bitquad attributes for every live node.
     *
     * @details
     * This entry point first computes family deltas, aggregates them into final
     * per-node family counts, and then projects only the requested scalar
     * descriptors into the caller-provided flat buffer. Rows are dense internal
     * node ids interpreted by `context.attrNames`. This topology-only overload
     * is valid for component trees with a regular adjacency relation.
     */
    template <std::floating_point Real>
    static void compute(const AttributeComputeContext<Real>& context) {
        requireAttributeBufferShape(context.tree, context.buffer, context.attrNames);
        const auto familyDeltas = detail::BitquadLocalEventComputation::computeBitquadFamilyDeltas(context.tree);
        const auto familyCounts = detail::BitquadLocalEventComputation::aggregateBitquadFamilyDeltas(context.tree, familyDeltas);
        detail::BitquadAttributeMaterialization::materializeAttributesFromBitquadFamilyCounts(
            context.tree,
            familyCounts,
            context.buffer,
            context.attrNames,
            context.requestedAttributes);
    }

    /**
     * @brief Computes requested bitquad attributes using a generic altitude span.
     *
     * @details
     * This overload exists for weighted Tree of Shapes callers whose altitude
     * buffer is not stored as `uint8_t`. Component-tree inputs still derive
     * scalar connectivity from the tree adjacency relation.
     */
    template<std::floating_point Real, AltitudeValue T>
    static void compute(const AltitudeAttributeComputeContext<Real, T>& context) {
        requireAttributeBufferShape(context.tree, context.buffer, context.attrNames);
        const auto familyDeltas = detail::BitquadLocalEventComputation::computeBitquadFamilyDeltas(context.tree);
        const auto familyCounts = detail::BitquadLocalEventComputation::aggregateBitquadFamilyDeltas(context.tree, familyDeltas);
        detail::BitquadAttributeMaterialization::materializeAttributesFromBitquadFamilyCounts(
            context.tree,
            context.altitude,
            familyCounts,
            context.buffer,
            context.attrNames,
            context.requestedAttributes);
    }

    /**
     * @brief Materializes bitquad attributes for unit proper-part supports.
     *
     * @details
     * The unit path defines the scalar bitquad constants for one-pixel
     * supports. It is used by the attribute pipeline when values must be
     * supplied for unit leaves outside the internal tree-node accumulation path.
     *
     * Unit constants depend on the tree adjacency relation. The method rejects
     * trees without regular adjacency metadata.
     */
    template <std::floating_point Real>
    static void computeUnitRows(const UnitAttributeComputeContext<Real>& context) {
        requireUnitAttributeBufferShape(context.tree, context.unitProperParts, context.buffer, context.attrNames);
        const AdjacencyRelation* adjacency = context.tree.getAdjacencyRelation();
        if (adjacency == nullptr) {
            throw std::invalid_argument("Local-event BitQuads attributes require an adjacency relation.");
        }

        const bool is4Connectivity = adjacency->is4connectivity();
        auto unitValue = [&](Attribute attribute) -> Real {
            if (is4Connectivity) {
                switch (attribute) {
                    case BITQUADS_AREA: return Real{0};
                    case BITQUADS_NUMBER_EULER: return Real{1};
                    case BITQUADS_NUMBER_HOLES: return Real{0};
                    case BITQUADS_PERIMETER: return Real{0};
                    case BITQUADS_PERIMETER_CONTINUOUS: return Real{0};
                    case BITQUADS_CIRCULARITY: return Real{0};
                    case BITQUADS_PERIMETER_AVERAGE: return Real{0};
                    case BITQUADS_LENGTH_AVERAGE: return Real{0};
                    case BITQUADS_WIDTH_AVERAGE: return Real{0};
                    default: break;
                }
            }

            switch (attribute) {
                case BITQUADS_AREA: return Real{1};
                case BITQUADS_NUMBER_EULER: return Real{0};
                case BITQUADS_NUMBER_HOLES: return Real{1};
                case BITQUADS_PERIMETER: return Real{4};
                case BITQUADS_PERIMETER_CONTINUOUS: return Real{8} / Real{3};
                case BITQUADS_CIRCULARITY: return Real{9} * std::numbers::pi_v<Real> / Real{16};
                case BITQUADS_PERIMETER_AVERAGE: return Real{0};
                case BITQUADS_LENGTH_AVERAGE: return Real{0};
                case BITQUADS_WIDTH_AVERAGE: return static_cast<Real>(0.75);
                default: break;
            }
            throw std::runtime_error("Unsupported local-event BitQuads unit attribute.");
        };

        constexpr std::array<Attribute, 9> unitAttributes{
            BITQUADS_AREA,
            BITQUADS_NUMBER_EULER,
            BITQUADS_NUMBER_HOLES,
            BITQUADS_PERIMETER,
            BITQUADS_PERIMETER_CONTINUOUS,
            BITQUADS_CIRCULARITY,
            BITQUADS_PERIMETER_AVERAGE,
            BITQUADS_LENGTH_AVERAGE,
            BITQUADS_WIDTH_AVERAGE};

        for (Attribute attribute : unitAttributes) {
            if (!requestsAttribute(context.requestedAttributes, attribute)) {
                continue;
            }
            const Real value = unitValue(attribute);
            for (NodeId leafIndex = 0; leafIndex < static_cast<NodeId>(context.unitProperParts.size()); ++leafIndex) {
                context.buffer[context.attrNames.linearIndex(leafIndex, attribute)] = value;
            }
        }
    }
};

} // namespace mmcfilters::attributes::computers
