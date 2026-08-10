#pragma once

#include "AttributeComputerDomain.hpp"
#include "AttributeComputerFamily.hpp"
#include "../detail/AttributeKernelSupport.hpp"
#include "detail/BitquadAttributeMaterialization.hpp"
#include "detail/BitquadLocalEventComputation.hpp"
#include "../../trees/TreeAltitudeAlgorithms.hpp"
#include "../../trees/detail/CommittedTreeAccess.hpp"
#include "../../utils/Contract.hpp"

#include <array>
#include <numbers>
#include <span>
#include <stdexcept>
#include <string_view>

namespace mmcfilters::attributes::computers {

namespace detail {

namespace kernel {

/**
 * @brief Computes requested bitquad attributes with established uniform connectivity.
 * @param context Established tree, output layout, and output buffer.
 * @param request Bitquad columns to materialize.
 */
template <std::floating_point Real>
inline void computeBitquadAttributes(const AttributeComputeContext<Real>& context, const BitquadRequest& request) {
    if (!request.any()) {
        return;
    }
    const std::vector<BitquadFamilyCounts> familyCounts = computeBitquadFamilyCounts(context.tree);
    const RegularGridAdjacency2D* uniformAdjacency = context.tree.getUniformGridAdjacency2D();
    const bool is4Connectivity = uniformAdjacency != nullptr ? uniformAdjacency->is4connectivity()
                                                               : context.tree.getDecreasingGridAdjacency2D()->is4connectivity();
    materializeBitquadAttributes(context, request, familyCounts, [is4Connectivity](NodeId) { return is4Connectivity; });
}

/**
 * @brief Computes requested bitquad attributes with established typed altitudes.
 * @param context Established tree, altitude span, output layout, and output buffer.
 * @param request Bitquad columns to materialize.
 */
template <std::floating_point Real, AltitudeValue T>
inline void computeBitquadAttributes(const AltitudeAttributeComputeContext<Real, T>& context, const BitquadRequest& request) {
    if (!request.any()) {
        return;
    }
    const std::vector<BitquadFamilyCounts> familyCounts = computeBitquadFamilyCounts(context.tree);
    const AttributeComputeContext<Real> outputContext{context.tree, context.buffer, context.attrNames, context.requestedAttributes,
                                                      context.dependencySources};
    const RegularGridAdjacency2D* uniformAdjacency = context.tree.getUniformGridAdjacency2D();
    if (uniformAdjacency != nullptr) {
        const bool is4Connectivity = uniformAdjacency->is4connectivity();
        materializeBitquadAttributes(outputContext, request, familyCounts, [is4Connectivity](NodeId) { return is4Connectivity; });
        return;
    }

    const bool decreasingIs4Connectivity = context.tree.getDecreasingGridAdjacency2D()->is4connectivity();
    const bool increasingIs4Connectivity = context.tree.getIncreasingGridAdjacency2D()->is4connectivity();
    const NodeId root = context.tree.getRoot();
    materializeBitquadAttributes(outputContext, request, familyCounts, [&](NodeId node) {
        if (node == root) {
            return decreasingIs4Connectivity == increasingIs4Connectivity ? decreasingIs4Connectivity : false;
        }
        const NodeId parent = ::mmcfilters::detail::CommittedTreeAccess::nodeParent(context.tree, node);
        if (context.altitude[static_cast<std::size_t>(node)] > context.altitude[static_cast<std::size_t>(parent)]) {
            return increasingIs4Connectivity;
        }
        if (context.altitude[static_cast<std::size_t>(node)] < context.altitude[static_cast<std::size_t>(parent)]) {
            return decreasingIs4Connectivity;
        }
        return decreasingIs4Connectivity;
    });
}

} // namespace kernel

inline void validateBitquadAdjacency(const MorphologicalTree& tree, bool altitudeAvailable) {
    const RegularGridAdjacency2D* uniformAdjacency = tree.getUniformGridAdjacency2D();
    if (uniformAdjacency != nullptr) {
        if (!uniformAdjacency->isCanonical4Or8Connectivity()) {
            throw std::invalid_argument("Local-event BitQuads scalar projection requires canonical 4- or 8-connectivity.");
        }
        return;
    }

    const RegularGridAdjacency2D* decreasingAdjacency = tree.getDecreasingGridAdjacency2D();
    const RegularGridAdjacency2D* increasingAdjacency = tree.getIncreasingGridAdjacency2D();
    if (decreasingAdjacency == nullptr || increasingAdjacency == nullptr) {
        throw std::invalid_argument("Local-event BitQuads scalar attributes require a regular or dual adjacency context.");
    }
    if (!decreasingAdjacency->isCanonical4Or8Connectivity() || !increasingAdjacency->isCanonical4Or8Connectivity()) {
        throw std::invalid_argument("Local-event BitQuads scalar projection requires canonical 4- or 8-connectivity.");
    }
    if (!altitudeAvailable && decreasingAdjacency->is4connectivity() != increasingAdjacency->is4connectivity()) {
        throw std::invalid_argument("Local-event BitQuads scalar projection requires an altitude buffer when decreasing and increasing connectivity differ.");
    }
}

template <std::floating_point Real> inline void validateBitquadContext(const AttributeComputeContext<Real>& context) {
    requireAttributeBufferShape(context.tree, context.buffer, context.attrNames);
    requireRequestedAttributeColumns(context);
    local_events::detail::validateEventEngineInput(context.tree, kernel::bitquadWindows[0]);
    validateBitquadAdjacency(context.tree, false);
}

template <std::floating_point Real, AltitudeValue T>
inline void validateBitquadContext(const AltitudeAttributeComputeContext<Real, T>& context) {
    requireAttributeBufferShape(context.tree, context.buffer, context.attrNames);
    requireRequestedAttributeColumns(context);
    local_events::detail::validateEventEngineInput(context.tree, kernel::bitquadWindows[0]);
    validateBitquadAdjacency(context.tree, true);
    TreeAltitudeAlgorithms::validateAltitudeBufferShape(context.tree, context.altitude);

    const RegularGridAdjacency2D* decreasingAdjacency = context.tree.getDecreasingGridAdjacency2D();
    const RegularGridAdjacency2D* increasingAdjacency = context.tree.getIncreasingGridAdjacency2D();
    if (context.tree.getUniformGridAdjacency2D() != nullptr || decreasingAdjacency->is4connectivity() == increasingAdjacency->is4connectivity()) {
        return;
    }
    for (NodeId node : context.tree.getAliveNodeIds()) {
        if (context.tree.isRoot(node)) {
            continue;
        }
        const NodeId parent = context.tree.getNodeParent(node);
        if (context.altitude[static_cast<std::size_t>(node)] == context.altitude[static_cast<std::size_t>(parent)]) {
            throw std::runtime_error("Local-event BitQuads scalar projection cannot select decreasing or increasing connectivity from equal node and parent altitudes.");
        }
    }
}

} // namespace detail

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

    /// Stable family id used by the scheduler.
    static constexpr AttributeComputerFamily family = AttributeComputerFamily::Bitquad;

    /// Execution domain required by the computer.
    static constexpr AttributeComputerDomain domain = AttributeComputerDomain::Topology;

    /**
     * @brief Canonical list of scalar bitquad attributes materialized by this computer.
     */
    inline static constexpr std::array<Attribute, 9> producedAttributes{
        BITQUADS_AREA,        BITQUADS_NUMBER_EULER,      BITQUADS_NUMBER_HOLES,   BITQUADS_PERIMETER,    BITQUADS_PERIMETER_CONTINUOUS,
        BITQUADS_CIRCULARITY, BITQUADS_PERIMETER_AVERAGE, BITQUADS_LENGTH_AVERAGE, BITQUADS_WIDTH_AVERAGE};

    /**
     * @brief Computes requested scalar bitquad attributes for every live node.
     *
     * @details
     * This entry point first computes family deltas, aggregates them into final
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
        detail::kernel::computeBitquadAttributes(context, request);
    }

    /**
     * @brief Computes requested bitquad attributes using a generic altitude span.
     *
     * @details
     * This overload exists for weighted Tree of Shapes callers whose altitude
     * buffer is not stored as `uint8_t`. Component-tree inputs still derive
     * scalar connectivity from the tree adjacency relation.
     *
     * @param context Operation context or diagnostic label.
     */
    template <std::floating_point Real, AltitudeValue T> static void compute(const AltitudeAttributeComputeContext<Real, T>& context) {
        const detail::BitquadRequest request = detail::BitquadRequest::from(context.requestedAttributes);
        MMCFILTERS_CONTRACT_CHECKED_ONLY(detail::validateBitquadContext(context));
        detail::kernel::computeBitquadAttributes(context, request);
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
        requireUnitAttributeBufferShape(context.tree, context.unitProperParts, context.buffer, context.attrNames);
        const RegularGridAdjacency2D* adjacency = context.tree.getUniformGridAdjacency2D();
        if (adjacency == nullptr) {
            throw std::invalid_argument("Local-event BitQuads attributes require an adjacency relation.");
        }
        if (!adjacency->isCanonical4Or8Connectivity()) {
            throw std::invalid_argument("Local-event BitQuads attributes require canonical 4- or 8-connectivity.");
        }

        const bool is4Connectivity = adjacency->is4connectivity();
        auto unitValue = [&](Attribute attribute) -> Real {
            if (is4Connectivity) {
                switch (attribute) {
                case BITQUADS_AREA:
                    return Real{0};
                case BITQUADS_NUMBER_EULER:
                    return Real{1};
                case BITQUADS_NUMBER_HOLES:
                    return Real{0};
                case BITQUADS_PERIMETER:
                    return Real{0};
                case BITQUADS_PERIMETER_CONTINUOUS:
                    return Real{0};
                case BITQUADS_CIRCULARITY:
                    return Real{0};
                case BITQUADS_PERIMETER_AVERAGE:
                    return Real{0};
                case BITQUADS_LENGTH_AVERAGE:
                    return Real{0};
                case BITQUADS_WIDTH_AVERAGE:
                    return Real{0};
                default:
                    break;
                }
            }

            switch (attribute) {
            case BITQUADS_AREA:
                return Real{1};
            case BITQUADS_NUMBER_EULER:
                return Real{0};
            case BITQUADS_NUMBER_HOLES:
                return Real{1};
            case BITQUADS_PERIMETER:
                return Real{4};
            case BITQUADS_PERIMETER_CONTINUOUS:
                return Real{8} / Real{3};
            case BITQUADS_CIRCULARITY:
                return Real{9} * std::numbers::pi_v<Real> / Real{16};
            case BITQUADS_PERIMETER_AVERAGE:
                return Real{0};
            case BITQUADS_LENGTH_AVERAGE:
                return Real{0};
            case BITQUADS_WIDTH_AVERAGE:
                return static_cast<Real>(0.75);
            default:
                break;
            }
            throw std::runtime_error("Unsupported local-event BitQuads unit attribute.");
        };

        constexpr std::array<Attribute, 9> unitAttributes{
            BITQUADS_AREA,        BITQUADS_NUMBER_EULER,      BITQUADS_NUMBER_HOLES,   BITQUADS_PERIMETER,    BITQUADS_PERIMETER_CONTINUOUS,
            BITQUADS_CIRCULARITY, BITQUADS_PERIMETER_AVERAGE, BITQUADS_LENGTH_AVERAGE, BITQUADS_WIDTH_AVERAGE};

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
