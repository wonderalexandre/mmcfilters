#pragma once

#include "AttributeComputerDomain.hpp"
#include "AttributeComputerFamily.hpp"
#include "../detail/AttributeKernelSupport.hpp"
#include "detail/BitquadAttributeProjection.hpp"
#include "detail/BitquadFiniteWindowComputation.hpp"
#include "../../trees/TreeAltitudeAlgorithms.hpp"
#include "../../trees/detail/MorphologicalTreeConstructionContextQueries.hpp"
#include "../../utils/Contract.hpp"

#include <array>
#include <numbers>
#include <span>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace mmcfilters::attributes::computers {

namespace detail {

inline void validateBitquadAdjacency(const MorphologicalTree& tree, bool altitudeAvailable) {
    static_cast<void>(makeBitquadConnectivityPolicy(tree, altitudeAvailable));
}

template <std::floating_point Real> inline void validateBitquadContext(const AttributeComputeContext<Real>& context) {
    requireAttributeBufferShape(context.tree, context.buffer, context.attrNames);
    requireRequestedAttributeColumns(context);
    local_attributes::detail::validateFiniteWindowLocalAttributeInput(context.tree);
    validateBitquadAdjacency(context.tree, false);
}

template <std::floating_point Real, AltitudeValue T> inline void validateBitquadContext(const AltitudeAttributeComputeContext<Real, T>& context) {
    requireAttributeBufferShape(context.tree, context.buffer, context.attrNames);
    requireRequestedAttributeColumns(context);
    local_attributes::detail::validateFiniteWindowLocalAttributeInput(context.tree);
    validateBitquadAdjacency(context.tree, true);
    TreeAltitudeAlgorithms::validateNodeAltitudeBufferShape(context.tree, context.altitude);

    const BitquadConnectivityPolicy connectivityPolicy = makeBitquadConnectivityPolicy(context.tree, true);
    if (!connectivityPolicy.requiresShapePolarity()) {
        return;
    }
    for (NodeId node : context.tree.aliveNodeIds()) {
        static_cast<void>(shapePolarity(context.tree, context.altitude, node));
    }
}

} // namespace detail

/**
 * @brief Bitquad scalar computer backed by finite-window computation.
 *
 * @details
 * This class obtains bitquad family increments from its finite-window computation,
 * aggregates them into per-node family counts, and materializes the public
 * `BITQUAD_*` descriptors. The public attribute pipeline/topology backend
 * invokes this computer for bitquad requests.
 *
 * The computer exposes the standard scalar bitquad descriptors. For component
 * trees, connectivity comes from the tree's adjacency relation. For a tree of
 * shapes whose lower- and upper-shape connectivity differ, scalar projection
 * derives each non-root shape's polarity from its exact node altitude and its
 * parent's altitude. The root has no shape polarity and is handled explicitly
 * by the projection policy.
 */
class BitquadAttributeComputer {
  public:
    /// Family name used in dependency-plan diagnostics.
    static constexpr std::string_view familyName = "bitquad";

    /// Stable family id used by the scheduler.
    static constexpr AttributeComputerFamily family = AttributeComputerFamily::Bitquad;

    /// Execution domain required by the computer.
    static constexpr AttributeComputerDomain domain = AttributeComputerDomain::Topology;

    /**
     * @brief Canonical list of scalar bitquad attributes materialized by this computer.
     */
    inline static constexpr std::array<Attribute, 9> producedAttributes{
        BitquadArea,        BitquadNumberEuler,      BitquadNumberHoles,   BitquadPerimeter,    BitquadPerimeterContinuous,
        BitquadCircularity, BitquadPerimeterAverage, BitquadLengthAverage, BitquadWidthAverage};

    /**
     * @brief Computes requested scalar bitquad attributes for every live node.
     *
     * @details
     * This entry point first computes family increments, aggregates them into final
     * per-node family counts, and then projects only the requested scalar
     * descriptors into the caller-provided flat buffer. Rows are dense internal
     * node ids interpreted by `context.attrNames`. This topology-only overload
     * is valid for component trees with a regular adjacency relation.
     *
     * @param context Operation context or diagnostic label.
     */
    template <std::floating_point Real> static void compute(const AttributeComputeContext<Real>& context) {
        const detail::BitquadRequest request = detail::BitquadRequest::from(context.requestedAttributes);
        MMCFILTERS_CONTRACT_CHECKED_ONLY(detail::validateBitquadContext(context));
        if (!request.any()) {
            return;
        }
        const std::vector<detail::BitquadFamilyCounts> familyCounts = detail::kernel::computeBitquadFamilyCounts(context.tree);
        const detail::BitquadConnectivityPolicy connectivityPolicy = detail::makeBitquadConnectivityPolicy(context.tree, false);
        detail::BitquadAttributeProjection::materializeBitquadAttributes(context.tree, std::span<const detail::BitquadFamilyCounts>(familyCounts),
                                                                         connectivityPolicy, context.buffer, context.attrNames, context.requestedAttributes);
    }

    /**
     * @brief Computes requested bitquad attributes using a generic altitude span.
     *
     * @details
     * This overload exists for valued-tree callers with a generic exact
     * altitude type. Component-tree inputs still derive scalar connectivity
     * from the tree adjacency relation.
     *
     * @param context Operation context or diagnostic label.
     */
    template <std::floating_point Real, AltitudeValue T> static void compute(const AltitudeAttributeComputeContext<Real, T>& context) {
        const detail::BitquadRequest request = detail::BitquadRequest::from(context.requestedAttributes);
        MMCFILTERS_CONTRACT_CHECKED_ONLY(detail::validateBitquadContext(context));
        if (!request.any()) {
            return;
        }
        const std::vector<detail::BitquadFamilyCounts> familyCounts = detail::kernel::computeBitquadFamilyCounts(context.tree);
        const detail::BitquadConnectivityPolicy connectivityPolicy = detail::makeBitquadConnectivityPolicy(context.tree, true);
        detail::BitquadAttributeProjection::materializeBitquadAttributes(context.tree, context.altitude,
                                                                         std::span<const detail::BitquadFamilyCounts>(familyCounts), connectivityPolicy,
                                                                         context.buffer, context.attrNames, context.requestedAttributes);
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
     *
     * @param context Operation context or diagnostic label.
     */
    template <std::floating_point Real> static void computeUnitRows(const UnitAttributeComputeContext<Real>& context) {
        requireUnitAttributeBufferShape(context.tree, context.unitPixels, context.buffer, context.attrNames);
        const RegularGridAdjacency2D* adjacency = ::mmcfilters::detail::constructionAdjacency(context.tree);
        if (adjacency == nullptr) {
            throw std::invalid_argument("Bitquad attributes require an adjacency relation.");
        }
        if (!adjacency->isCanonical4Or8Connectivity()) {
            throw std::invalid_argument("Bitquad attributes require canonical 4- or 8-connectivity.");
        }

        const bool is4Connectivity = adjacency->is4connectivity();
        auto unitValue = [&](Attribute attribute) -> Real {
            if (is4Connectivity) {
                switch (attribute) {
                case BitquadArea:
                    return Real{0};
                case BitquadNumberEuler:
                    return Real{1};
                case BitquadNumberHoles:
                    return Real{0};
                case BitquadPerimeter:
                    return Real{0};
                case BitquadPerimeterContinuous:
                    return Real{0};
                case BitquadCircularity:
                    return Real{0};
                case BitquadPerimeterAverage:
                    return Real{0};
                case BitquadLengthAverage:
                    return Real{0};
                case BitquadWidthAverage:
                    return Real{0};
                default:
                    break;
                }
            }

            switch (attribute) {
            case BitquadArea:
                return Real{1};
            case BitquadNumberEuler:
                return Real{0};
            case BitquadNumberHoles:
                return Real{1};
            case BitquadPerimeter:
                return Real{4};
            case BitquadPerimeterContinuous:
                return Real{8} / Real{3};
            case BitquadCircularity:
                return Real{9} * std::numbers::pi_v<Real> / Real{16};
            case BitquadPerimeterAverage:
                return Real{0};
            case BitquadLengthAverage:
                return Real{0};
            case BitquadWidthAverage:
                return static_cast<Real>(0.75);
            default:
                break;
            }
            throw std::runtime_error("Unsupported finite-window Bitquad unit attribute.");
        };

        constexpr std::array<Attribute, 9> unitAttributes{
            BitquadArea,        BitquadNumberEuler,      BitquadNumberHoles,   BitquadPerimeter,    BitquadPerimeterContinuous,
            BitquadCircularity, BitquadPerimeterAverage, BitquadLengthAverage, BitquadWidthAverage};

        for (Attribute attribute : unitAttributes) {
            if (!requestsAttribute(context.requestedAttributes, attribute)) {
                continue;
            }
            const Real value = unitValue(attribute);
            for (NodeId leafIndex = 0; leafIndex < static_cast<NodeId>(context.unitPixels.size()); ++leafIndex) {
                context.buffer[context.attrNames.linearIndex(leafIndex, attribute)] = value;
            }
        }
    }
};

} // namespace mmcfilters::attributes::computers
