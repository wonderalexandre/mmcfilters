#pragma once

#include "BitquadAttributeData.hpp"
#include "../../detail/AttributeKernelSupport.hpp"
#include "../../../trees/MorphologicalTree.hpp"
#include "../../../trees/TreeAltitudeAlgorithms.hpp"
#include "../../../trees/WeightedTreeView.hpp"
#include "../../../trees/detail/CommittedTreeAccess.hpp"
#include "../../../trees/detail/TreeTraversalDetail.hpp"
#include "../../../utils/Contract.hpp"

#include <concepts>
#include <cstddef>
#include <numbers>
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
        return {.area = requestsAttribute(requestedAttributes, BITQUADS_AREA),
                .numberEuler = requestsAttribute(requestedAttributes, BITQUADS_NUMBER_EULER),
                .numberHoles = requestsAttribute(requestedAttributes, BITQUADS_NUMBER_HOLES),
                .perimeter = requestsAttribute(requestedAttributes, BITQUADS_PERIMETER),
                .perimeterContinuous = requestsAttribute(requestedAttributes, BITQUADS_PERIMETER_CONTINUOUS),
                .circularity = requestsAttribute(requestedAttributes, BITQUADS_CIRCULARITY),
                .perimeterAverage = requestsAttribute(requestedAttributes, BITQUADS_PERIMETER_AVERAGE),
                .lengthAverage = requestsAttribute(requestedAttributes, BITQUADS_LENGTH_AVERAGE),
                .widthAverage = requestsAttribute(requestedAttributes, BITQUADS_WIDTH_AVERAGE)};
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
inline double bitquadContinuousPerimeter(const BitquadScalarCounters& counters) noexcept {
    return counters.c2 + (counters.c1 + counters.c3) / 1.5;
}

/**
 * @brief Materializes requested bitquad scalars from established family counts.
 * @param context Established tree, output layout, and output buffer.
 * @param request Scalar columns to materialize.
 * @param familyCounts Per-node bitquad-family counts.
 * @param uses4Connectivity Callable selecting connectivity for each node.
 */
template <std::floating_point Real, class ConnectivitySelector>
inline void materializeBitquadAttributes(const AttributeComputeContext<Real>& context, const BitquadRequest& request,
                                         std::span<const BitquadFamilyCounts> familyCounts, ConnectivitySelector uses4Connectivity) {
    if (!request.any()) {
        return;
    }

    const int stride = context.attrNames.NUM_ATTRIBUTES;
    const auto offsetOf = [&](Attribute attribute) { return context.attrNames.indexMap.find(attribute)->second; };
    const int areaOffset = request.area ? offsetOf(BITQUADS_AREA) : 0;
    const int eulerOffset = request.numberEuler ? offsetOf(BITQUADS_NUMBER_EULER) : 0;
    const int holesOffset = request.numberHoles ? offsetOf(BITQUADS_NUMBER_HOLES) : 0;
    const int perimeterOffset = request.perimeter ? offsetOf(BITQUADS_PERIMETER) : 0;
    const int continuousPerimeterOffset = request.perimeterContinuous ? offsetOf(BITQUADS_PERIMETER_CONTINUOUS) : 0;
    const int circularityOffset = request.circularity ? offsetOf(BITQUADS_CIRCULARITY) : 0;
    const int perimeterAverageOffset = request.perimeterAverage ? offsetOf(BITQUADS_PERIMETER_AVERAGE) : 0;
    const int lengthAverageOffset = request.lengthAverage ? offsetOf(BITQUADS_LENGTH_AVERAGE) : 0;
    const int widthAverageOffset = request.widthAverage ? offsetOf(BITQUADS_WIDTH_AVERAGE) : 0;
    const auto outputIndex = [&](NodeId node, int offset) { return static_cast<std::size_t>(node * stride + offset); };

    ::mmcfilters::detail::kernel::traversePostOrder(
        context.tree, context.tree.getRoot(), [](NodeId) {}, [](NodeId, NodeId) {},
        [&](NodeId node) {
            const BitquadScalarCounters counters =
                bitquadScalarCounters(familyCounts[static_cast<std::size_t>(node)], uses4Connectivity(node));
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
                context.buffer[outputIndex(node, circularityOffset)] = static_cast<Real>(
                    ::mmcfilters::attributes::numeric::safeDivide(4.0 * std::numbers::pi * area, continuousPerimeter * continuousPerimeter));
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
 * @brief Internal scalar projection for precomputed bitquad-family buckets.
 *
 * @details
 * Local-event code stores topology support information as dense
 * `BitquadFamilyCounts` buffers. This helper owns the conversion from those
 * implementation buckets to public `BITQUADS_*` scalar attributes. Keeping the
 * bucket-aware projection in `detail` prevents the public computer API from
 * exposing local-event storage types.
 */
class BitquadAttributeMaterialization {
  public:
    /**
     * @brief Projects family counters for component-tree inputs.
     *
     * A regular adjacency applies to every node. A dual adjacency policy can
     * also be used without altitudes when both branch directions have the same
     * connectivity.
     *
     * @param tree Tree topology used by the operation.
     * @param familyCounts Per-node counters for the attribute family.
     * @param buffer Buffer read or written by the operation.
     * @param attrNames Layout that maps attributes to buffer columns.
     * @param requestedAttributes Attributes requested for materialization.
     */
    template <std::floating_point Real>
    static void materializeAttributesFromBitquadFamilyCounts(const MorphologicalTree& tree, std::span<const BitquadFamilyCounts> familyCounts,
                                                             std::span<Real> buffer, const AttributeNames& attrNames,
                                                             std::span<const Attribute> requestedAttributes) {
        const std::size_t numNodeSlots = static_cast<std::size_t>(tree.getNumInternalNodeSlots());
        MMCFILTERS_CONTRACT_REQUIRE(familyCounts.size() >= numNodeSlots,
                                    throw std::invalid_argument("Local-event bitquad family counts do not cover all tree node slots."));

        const RegularGridAdjacency2D* adjacency = tree.getUniformGridAdjacency2D();
        if (adjacency != nullptr) {
            MMCFILTERS_CONTRACT_REQUIRE(
                adjacency->isCanonical4Or8Connectivity(),
                throw std::invalid_argument("Local-event BitQuads scalar projection requires canonical 4- or 8-connectivity."));
            const AttributeComputeContext<Real> context{tree, buffer, attrNames, requestedAttributes};
            const BitquadRequest request = BitquadRequest::from(requestedAttributes);
            kernel::materializeBitquadAttributes(context, request, familyCounts,
                                                 [is4Connectivity = adjacency->is4connectivity()](NodeId) { return is4Connectivity; });
            return;
        }

        const RegularGridAdjacency2D* decreasingAdjacency = tree.getDecreasingGridAdjacency2D();
        const RegularGridAdjacency2D* increasingAdjacency = tree.getIncreasingGridAdjacency2D();
        MMCFILTERS_CONTRACT_REQUIRE(decreasingAdjacency != nullptr && increasingAdjacency != nullptr,
                                    throw std::invalid_argument("Local-event BitQuads scalar attributes require a regular or dual adjacency context."));
        MMCFILTERS_CONTRACT_REQUIRE(
            decreasingAdjacency->isCanonical4Or8Connectivity() && increasingAdjacency->isCanonical4Or8Connectivity(),
            throw std::invalid_argument("Local-event BitQuads scalar projection requires canonical 4- or 8-connectivity."));
        MMCFILTERS_CONTRACT_REQUIRE(
            decreasingAdjacency->is4connectivity() == increasingAdjacency->is4connectivity(),
            throw std::invalid_argument(
                "Local-event BitQuads scalar projection requires an altitude buffer when decreasing and increasing connectivity differ."));

        const AttributeComputeContext<Real> context{tree, buffer, attrNames, requestedAttributes};
        const BitquadRequest request = BitquadRequest::from(requestedAttributes);
        kernel::materializeBitquadAttributes(context, request, familyCounts,
                                             [is4Connectivity = decreasingAdjacency->is4connectivity()](NodeId) { return is4Connectivity; });
    }

    /**
     * @brief Projects family counters using a generic altitude span.
     *
     * @details
     * Trees with a regular adjacency ignore `altitude`. Trees with a dual
     * adjacency policy use the span to select the decreasing or increasing
     * relation for each non-root branch.
     *
     * @param tree Tree topology used by the operation.
     * @param altitude Altitude data indexed by node identifier.
     * @param familyCounts Per-node counters for the attribute family.
     * @param buffer Buffer read or written by the operation.
     * @param attrNames Layout that maps attributes to buffer columns.
     * @param requestedAttributes Attributes requested for materialization.
     *
     * @throws std::invalid_argument If family counts or altitude shape are
     * invalid or if the topology lacks regular and dual adjacency context.
     *
     */
    template <AltitudeValue T, std::floating_point Real>
    static void materializeAttributesFromBitquadFamilyCounts(const MorphologicalTree& tree, std::span<const T> altitude,
                                                             std::span<const BitquadFamilyCounts> familyCounts, std::span<Real> buffer,
                                                             const AttributeNames& attrNames, std::span<const Attribute> requestedAttributes) {
        const std::size_t numNodeSlots = static_cast<std::size_t>(tree.getNumInternalNodeSlots());
        MMCFILTERS_CONTRACT_REQUIRE(familyCounts.size() >= numNodeSlots,
                                    throw std::invalid_argument("Local-event bitquad family counts do not cover all tree node slots."));

        const RegularGridAdjacency2D* adjacency = tree.getUniformGridAdjacency2D();
        if (adjacency != nullptr) {
            MMCFILTERS_CONTRACT_REQUIRE(
                adjacency->isCanonical4Or8Connectivity(),
                throw std::invalid_argument("Local-event BitQuads scalar projection requires canonical 4- or 8-connectivity."));
            const AttributeComputeContext<Real> context{tree, buffer, attrNames, requestedAttributes};
            const BitquadRequest request = BitquadRequest::from(requestedAttributes);
            kernel::materializeBitquadAttributes(
                context, request, familyCounts,
                [is4Connectivity = adjacency->is4connectivity()](NodeId) { return is4Connectivity; });
            return;
        }

        MMCFILTERS_CONTRACT_REQUIRE(tree.hasDirectionalGridAdjacency2D(),
                                    throw std::invalid_argument("Local-event BitQuads scalar attributes require a regular or dual adjacency context."));
        TreeAltitudeAlgorithms::validateAltitudeBufferShape(tree, altitude);
        const RegularGridAdjacency2D* decreasingAdjacency = tree.getDecreasingGridAdjacency2D();
        const RegularGridAdjacency2D* increasingAdjacency = tree.getIncreasingGridAdjacency2D();
        MMCFILTERS_CONTRACT_REQUIRE(decreasingAdjacency != nullptr && increasingAdjacency != nullptr,
                                    throw std::invalid_argument("Local-event BitQuads scalar projection requires decreasing and increasing adjacency relations."));
        MMCFILTERS_CONTRACT_REQUIRE(
            decreasingAdjacency->isCanonical4Or8Connectivity() && increasingAdjacency->isCanonical4Or8Connectivity(),
            throw std::invalid_argument("Local-event BitQuads scalar projection requires canonical 4- or 8-connectivity."));

        const bool decreasingIs4Connectivity = decreasingAdjacency->is4connectivity();
        const bool increasingIs4Connectivity = increasingAdjacency->is4connectivity();
        const NodeId root = tree.getRoot();
        if constexpr (contract::validationsEnabled) {
            if (decreasingIs4Connectivity != increasingIs4Connectivity) {
                for (NodeId node : tree.getAliveNodeIds()) {
                    if (node == root) {
                        continue;
                    }
                    const NodeId parent = ::mmcfilters::detail::CommittedTreeAccess::nodeParent(tree, node);
                    if (altitude[static_cast<std::size_t>(node)] == altitude[static_cast<std::size_t>(parent)]) {
                        throw std::runtime_error("Local-event BitQuads scalar projection cannot select decreasing or increasing connectivity from equal node and parent altitudes.");
                    }
                }
            }
        }

        const AttributeComputeContext<Real> context{tree, buffer, attrNames, requestedAttributes};
        const BitquadRequest request = BitquadRequest::from(requestedAttributes);
        kernel::materializeBitquadAttributes(context, request, familyCounts, [&](NodeId node) {
            if (node == root) {
                return decreasingIs4Connectivity == increasingIs4Connectivity ? decreasingIs4Connectivity : false;
            }
            const NodeId parent = ::mmcfilters::detail::CommittedTreeAccess::nodeParent(tree, node);
            if (altitude[static_cast<std::size_t>(node)] > altitude[static_cast<std::size_t>(parent)]) {
                return increasingIs4Connectivity;
            }
            return decreasingIs4Connectivity;
        });
    }

    /**
     * @brief Projects family counters using a non-owning weighted-tree view.
     *
     * The view must remain valid for the duration of the call and its altitude
     * span must cover the topology's dense internal node-slot domain.
     *
     * @param tree Tree topology used by the operation.
     * @param familyCounts Per-node counters for the attribute family.
     * @param buffer Buffer read or written by the operation.
     * @param attrNames Layout that maps attributes to buffer columns.
     * @param requestedAttributes Attributes requested for materialization.
     */
    template <AltitudeValue T, std::floating_point Real>
    static void materializeAttributesFromBitquadFamilyCounts(const WeightedTreeView<T>& tree, std::span<const BitquadFamilyCounts> familyCounts,
                                                             std::span<Real> buffer, const AttributeNames& attrNames,
                                                             std::span<const Attribute> requestedAttributes) {
        materializeAttributesFromBitquadFamilyCounts(tree.topology(), tree.altitude(), familyCounts, buffer, attrNames, requestedAttributes);
    }

};

} // namespace mmcfilters::attributes::computers::detail
