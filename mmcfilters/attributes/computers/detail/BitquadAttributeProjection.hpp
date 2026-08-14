#pragma once

#include "BitquadAttributeData.hpp"
#include "BitquadConnectivityPolicy.hpp"
#include "../../detail/AttributeKernelSupport.hpp"
#include "../../../trees/MorphologicalTree.hpp"
#include "../../../trees/ValuedMorphologicalTreeView.hpp"
#include "../../../trees/detail/TreeTraversalDetail.hpp"
#include "../../../utils/Contract.hpp"

#include <concepts>
#include <cstddef>
#include <numbers>
#include <optional>
#include <span>
#include <stdexcept>

namespace mmcfilters::attributes::computers::detail {

/** @brief Selection mask for scalar attributes derived from bitquad-family counts. */
struct BitquadRequest {
    bool area = false;                ///< Whether bitquad area is requested.
    bool numberEuler = false;         ///< Whether Euler number is requested.
    bool numberHoles = false;         ///< Whether the number of holes is requested.
    bool perimeter = false;           ///< Whether discrete perimeter is requested.
    bool perimeterContinuous = false; ///< Whether continuous perimeter is requested.
    bool circularity = false;         ///< Whether bitquad circularity is requested.
    bool perimeterAverage = false;    ///< Whether average perimeter is requested.
    bool lengthAverage = false;       ///< Whether average length is requested.
    bool widthAverage = false;        ///< Whether average width is requested.

    /** @brief Reports whether at least one scalar projection is requested. @return True when any request flag is set. */
    [[nodiscard]] bool any() const noexcept {
        return area || numberEuler || numberHoles || perimeter || perimeterContinuous || circularity || perimeterAverage || lengthAverage || widthAverage;
    }

    /**
     * @brief Builds the selection mask from requested scalar attributes.
     * @param requestedAttributes Requested scalar attributes.
     * @return Bitquad scalar selection mask.
     */
    [[nodiscard]] static BitquadRequest from(std::span<const Attribute> requestedAttributes) {
        return {.area = requestsAttribute(requestedAttributes, BitquadArea),
                .numberEuler = requestsAttribute(requestedAttributes, BitquadNumberEuler),
                .numberHoles = requestsAttribute(requestedAttributes, BitquadNumberHoles),
                .perimeter = requestsAttribute(requestedAttributes, BitquadPerimeter),
                .perimeterContinuous = requestsAttribute(requestedAttributes, BitquadPerimeterContinuous),
                .circularity = requestsAttribute(requestedAttributes, BitquadCircularity),
                .perimeterAverage = requestsAttribute(requestedAttributes, BitquadPerimeterAverage),
                .lengthAverage = requestsAttribute(requestedAttributes, BitquadLengthAverage),
                .widthAverage = requestsAttribute(requestedAttributes, BitquadWidthAverage)};
    }
};

namespace kernel {

/** @brief Connectivity-adjusted scalar counters used by bitquad formulas. */
struct BitquadScalarCounters {
    int c1c4 = 0; ///< Q1/QD combination used by the 4-connected Euler formula.
    int c1 = 0;   ///< Q1 count used by perimeter and 8-connected formulas.
    int c2 = 0;   ///< Q2 family count.
    int cd = 0;   ///< Diagonal QD family count.
    int c3 = 0;   ///< Q3 family count.
    int c4 = 0;   ///< Q4 family count.
};

/**
 * @brief Converts family counts into connectivity-specific scalar counters.
 * @param counts Per-node bitquad-family counts.
 * @param is4Connectivity Whether foreground connectivity is four-neighbour.
 * @return Counters used by scalar bitquad formulas.
 */
inline BitquadScalarCounters bitquadScalarCounters(const BitquadFamilyCounts& counts, bool is4Connectivity) noexcept {
    BitquadScalarCounters counters;
    if (is4Connectivity) {
        counters.c1c4 = counts.q1 + 2 * counts.qd;
    } else {
        counters.c1 = counts.q1;
        counters.cd = counts.qd;
    }
    counters.c2 = counts.q2;
    counters.c3 = counts.q3;
    counters.c4 = counts.q4;
    return counters;
}

/** @brief Computes Euler number from scalar counters. @param counters Scalar counters. @return Euler number. */
inline int bitquadEulerNumber(const BitquadScalarCounters& counters) noexcept { return (counters.c1c4 - counters.c3) / 4; }
/** @brief Computes discrete perimeter from scalar counters. @param counters Scalar counters. @return Discrete perimeter. */
inline int bitquadPerimeter(const BitquadScalarCounters& counters) noexcept { return counters.c1 + counters.c2 + counters.c3 + 2 * counters.cd; }
/** @brief Computes bitquad area from scalar counters. @param counters Scalar counters. @return Estimated area. */
inline double bitquadArea(const BitquadScalarCounters& counters) noexcept {
    return 0.25 * counters.c1 + 0.5 * counters.c2 + 0.875 * counters.c3 + counters.c4 + 0.75 * counters.cd;
}
/**
 * @brief Computes continuous perimeter from scalar counters.
 * @param counters Scalar counters.
 * @return Continuous perimeter estimate.
 */
inline double bitquadContinuousPerimeter(const BitquadScalarCounters& counters) noexcept { return counters.c2 + (counters.c1 + counters.c3) / 1.5; }

/**
 * @brief Materializes requested bitquad scalars from established family counts.
 * @param context Established tree, output layout, and output buffer.
 * @param request Scalar columns to materialize.
 * @param familyCounts Per-node bitquad-family counts.
 * @param connectivityPolicy Explicit root/lower/upper connectivity policy.
 * @param polarityOf Callable returning no polarity for the root and the derived polarity for non-root shapes.
 */
template <std::floating_point Real, class ShapePolarityResolver>
inline void materializeBitquadAttributes(const AttributeComputeContext<Real>& context, const BitquadRequest& request,
                                         std::span<const BitquadFamilyCounts> familyCounts, const BitquadConnectivityPolicy& connectivityPolicy,
                                         ShapePolarityResolver polarityOf) {
    if (!request.any()) {
        return;
    }

    const int stride = context.attrNames.NUM_ATTRIBUTES;
    const auto offsetOf = [&](Attribute attribute) { return context.attrNames.indexMap.find(attribute)->second; };
    const int areaOffset = request.area ? offsetOf(BitquadArea) : 0;
    const int eulerOffset = request.numberEuler ? offsetOf(BitquadNumberEuler) : 0;
    const int holesOffset = request.numberHoles ? offsetOf(BitquadNumberHoles) : 0;
    const int perimeterOffset = request.perimeter ? offsetOf(BitquadPerimeter) : 0;
    const int continuousPerimeterOffset = request.perimeterContinuous ? offsetOf(BitquadPerimeterContinuous) : 0;
    const int circularityOffset = request.circularity ? offsetOf(BitquadCircularity) : 0;
    const int perimeterAverageOffset = request.perimeterAverage ? offsetOf(BitquadPerimeterAverage) : 0;
    const int lengthAverageOffset = request.lengthAverage ? offsetOf(BitquadLengthAverage) : 0;
    const int widthAverageOffset = request.widthAverage ? offsetOf(BitquadWidthAverage) : 0;
    const auto outputIndex = [&](NodeId node, int offset) { return static_cast<std::size_t>(node * stride + offset); };

    ::mmcfilters::detail::kernel::traversePostOrder(
        context.tree, context.tree.root(), [](NodeId) {}, [](NodeId, NodeId) {},
        [&](NodeId node) {
            const BitquadScalarCounters counters =
                bitquadScalarCounters(familyCounts[static_cast<std::size_t>(node)], connectivityPolicy.uses4Connectivity(polarityOf(node)));
            const int euler = bitquadEulerNumber(counters);
            const double area = bitquadArea(counters);
            const double continuousPerimeter = bitquadContinuousPerimeter(counters);
            if (request.area)
                context.buffer[outputIndex(node, areaOffset)] = static_cast<Real>(area);
            if (request.numberEuler)
                context.buffer[outputIndex(node, eulerOffset)] = static_cast<Real>(euler);
            if (request.numberHoles)
                context.buffer[outputIndex(node, holesOffset)] = static_cast<Real>(1 - euler);
            if (request.perimeter)
                context.buffer[outputIndex(node, perimeterOffset)] = static_cast<Real>(bitquadPerimeter(counters));
            if (request.perimeterContinuous)
                context.buffer[outputIndex(node, continuousPerimeterOffset)] = static_cast<Real>(continuousPerimeter);
            if (request.circularity)
                context.buffer[outputIndex(node, circularityOffset)] =
                    static_cast<Real>(::mmcfilters::attributes::numeric::safeDivide(4.0 * std::numbers::pi * area, continuousPerimeter * continuousPerimeter));
            if (request.perimeterAverage)
                context.buffer[outputIndex(node, perimeterAverageOffset)] =
                    static_cast<Real>(::mmcfilters::attributes::numeric::safeDivide(continuousPerimeter, static_cast<double>(euler)));
            if (request.lengthAverage)
                context.buffer[outputIndex(node, lengthAverageOffset)] =
                    static_cast<Real>(::mmcfilters::attributes::numeric::safeDivide(continuousPerimeter, 2.0 * static_cast<double>(euler)));
            if (request.widthAverage)
                context.buffer[outputIndex(node, widthAverageOffset)] =
                    static_cast<Real>(::mmcfilters::attributes::numeric::safeDivide(2.0 * area, continuousPerimeter));
        });
}

} // namespace kernel

/**
 * @brief Connectivity-dependent scalar projection for precomputed bitquad-family counts.
 *
 * @details
 * Finite-window counting stores hierarchy-independent support information as
 * dense `BitquadFamilyCounts` buffers. This helper applies an explicit
 * connectivity policy and converts those counts into public scalar attributes.
 * Counting and connectivity-dependent projection remain independently testable.
 */
class BitquadAttributeProjection {
  public:
    /**
     * @brief Projects family counters for component-tree inputs.
     *
     * A regular adjacency applies to every node. A dual adjacency policy can
     * also be used without altitudes when both branch directions have the same
     * connectivity.
     *
     * @param tree Tree topology.
     * @param familyCounts Per-node counters for the attribute family.
     * @param connectivityPolicy Explicit scalar-projection connectivity policy.
     * @param buffer Buffer read or written by the operation.
     * @param attrNames Layout that maps attributes to buffer columns.
     * @param requestedAttributes Attributes requested for materialization.
     */
    template <std::floating_point Real>
    static void materializeBitquadAttributes(const MorphologicalTree& tree, std::span<const BitquadFamilyCounts> familyCounts,
                                             const BitquadConnectivityPolicy& connectivityPolicy, std::span<Real> buffer, const AttributeNames& attrNames,
                                             std::span<const Attribute> requestedAttributes) {
        const std::size_t numNodeSlots = static_cast<std::size_t>(tree.numInternalNodeSlots());
        MMCFILTERS_CONTRACT_REQUIRE(familyCounts.size() >= numNodeSlots,
                                    throw std::invalid_argument("Bitquad family counts do not cover all tree node slots."));
        MMCFILTERS_CONTRACT_REQUIRE(
            !connectivityPolicy.requiresShapePolarity(),
            throw std::invalid_argument(
                "Bitquad scalar materialization requires node altitudes for a connectivity policy that distinguishes lower and upper shapes."));

        const AttributeComputeContext<Real> context{tree, buffer, attrNames, requestedAttributes};
        const BitquadRequest request = BitquadRequest::from(requestedAttributes);
        kernel::materializeBitquadAttributes(context, request, familyCounts, connectivityPolicy, [](NodeId) { return std::optional<ShapePolarity>{}; });
    }

    /**
     * @brief Projects family counters using a generic altitude span.
     *
     * @details
     * Trees with a regular adjacency ignore `altitude`. Trees whose
     * topographic convention has unequal lower- and upper-shape connectivity
     * use the span to derive the polarity of every non-root shape. The root has
     * no shape polarity and is handled explicitly by the connectivity policy.
     *
     * @param tree Tree topology.
     * @param altitude Altitude data indexed by node identifier.
     * @param familyCounts Per-node counters for the attribute family.
     * @param connectivityPolicy Explicit scalar-projection connectivity policy.
     * @param buffer Buffer read or written by the operation.
     * @param attrNames Layout that maps attributes to buffer columns.
     * @param requestedAttributes Attributes requested for materialization.
     *
     * @throws std::invalid_argument If family counts or altitude shape are
     * invalid or if the topology lacks regular and dual adjacency context.
     *
     */
    template <AltitudeValue T, std::floating_point Real>
    static void materializeBitquadAttributes(const MorphologicalTree& tree, std::span<const T> altitude, std::span<const BitquadFamilyCounts> familyCounts,
                                             const BitquadConnectivityPolicy& connectivityPolicy, std::span<Real> buffer, const AttributeNames& attrNames,
                                             std::span<const Attribute> requestedAttributes) {
        const std::size_t numNodeSlots = static_cast<std::size_t>(tree.numInternalNodeSlots());
        MMCFILTERS_CONTRACT_REQUIRE(familyCounts.size() >= numNodeSlots,
                                    throw std::invalid_argument("Bitquad family counts do not cover all tree node slots."));

        const AttributeComputeContext<Real> context{tree, buffer, attrNames, requestedAttributes};
        const BitquadRequest request = BitquadRequest::from(requestedAttributes);
        if (connectivityPolicy.requiresShapePolarity()) {
            TreeAltitudeAlgorithms::validateNodeAltitudeBufferShape(tree, altitude);
        }
        kernel::materializeBitquadAttributes(context, request, familyCounts, connectivityPolicy, [&](NodeId node) {
            if (!connectivityPolicy.requiresShapePolarity()) {
                return std::optional<ShapePolarity>{};
            }
            return shapePolarity(tree, altitude, node);
        });
    }

    /**
     * @brief Projects family counters using a non-owning valued-tree view.
     *
     * The view must remain valid for the duration of the call and its altitude
     * span must cover the topology's dense internal node-slot domain.
     *
     * @param tree Tree topology.
     * @param familyCounts Per-node counters for the attribute family.
     * @param connectivityPolicy Explicit scalar-projection connectivity policy.
     * @param buffer Buffer read or written by the operation.
     * @param attrNames Layout that maps attributes to buffer columns.
     * @param requestedAttributes Attributes requested for materialization.
     */
    template <AltitudeValue T, std::floating_point Real>
    static void materializeBitquadAttributes(const ValuedMorphologicalTreeView<T>& tree, std::span<const BitquadFamilyCounts> familyCounts,
                                             const BitquadConnectivityPolicy& connectivityPolicy, std::span<Real> buffer, const AttributeNames& attrNames,
                                             std::span<const Attribute> requestedAttributes) {
        materializeBitquadAttributes(tree.topology(), tree.nodeAltitudes(), familyCounts, connectivityPolicy, buffer, attrNames, requestedAttributes);
    }
};

} // namespace mmcfilters::attributes::computers::detail
