#pragma once

#include "BitquadAttributeData.hpp"
#include "../../detail/AttributeKernelSupport.hpp"
#include "../../../trees/MorphologicalTree.hpp"
#include "../../../trees/TreeAltitudeAlgorithms.hpp"
#include "../../../trees/WeightedTreeView.hpp"

#include <concepts>
#include <cstddef>
#include <numbers>
#include <span>
#include <stdexcept>

namespace mmcfilters::attributes::computers::detail {

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
        if (familyCounts.size() < numNodeSlots) {
            throw std::invalid_argument("Local-event bitquad family counts do not cover all tree node slots.");
        }

        const RegularGridAdjacency2D* adjacency = tree.getUniformGridAdjacency2D();
        if (adjacency != nullptr) {
            if (!adjacency->isCanonical4Or8Connectivity()) {
                throw std::invalid_argument("Local-event BitQuads scalar projection requires canonical 4- or 8-connectivity.");
            }
            materializeAttributesFromBitquadFamilyCounts(
                tree, familyCounts, [is4Connectivity = adjacency->is4connectivity()](NodeId) { return is4Connectivity; }, buffer, attrNames,
                requestedAttributes);
            return;
        }

        const RegularGridAdjacency2D* decreasingAdjacency = tree.getDecreasingGridAdjacency2D();
        const RegularGridAdjacency2D* increasingAdjacency = tree.getIncreasingGridAdjacency2D();
        if ((decreasingAdjacency != nullptr && !decreasingAdjacency->isCanonical4Or8Connectivity()) ||
            (increasingAdjacency != nullptr && !increasingAdjacency->isCanonical4Or8Connectivity())) {
            throw std::invalid_argument("Local-event BitQuads scalar projection requires canonical 4- or 8-connectivity.");
        }
        if (decreasingAdjacency != nullptr && increasingAdjacency != nullptr &&
            decreasingAdjacency->is4connectivity() == increasingAdjacency->is4connectivity()) {
            materializeAttributesFromBitquadFamilyCounts(
                tree, familyCounts, [is4Connectivity = decreasingAdjacency->is4connectivity()](NodeId) { return is4Connectivity; }, buffer, attrNames,
                requestedAttributes);
            return;
        }

        if (tree.hasDirectionalGridAdjacency2D()) {
            throw std::invalid_argument("Local-event BitQuads scalar projection requires an altitude "
                                        "buffer when decreasing and increasing connectivity differ.");
        }
        throw std::invalid_argument("Local-event BitQuads scalar attributes require a regular or dual "
                                    "adjacency context.");
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
        if (familyCounts.size() < numNodeSlots) {
            throw std::invalid_argument("Local-event bitquad family counts do not cover all tree node slots.");
        }

        const RegularGridAdjacency2D* adjacency = tree.getUniformGridAdjacency2D();
        if (adjacency != nullptr) {
            if (!adjacency->isCanonical4Or8Connectivity()) {
                throw std::invalid_argument("Local-event BitQuads scalar projection requires canonical 4- or 8-connectivity.");
            }
            materializeAttributesFromBitquadFamilyCounts(
                tree, familyCounts, [is4Connectivity = adjacency->is4connectivity()](NodeId) { return is4Connectivity; }, buffer, attrNames,
                requestedAttributes);
            return;
        }

        if (!tree.hasDirectionalGridAdjacency2D()) {
            throw std::invalid_argument("Local-event BitQuads scalar attributes require a regular or "
                                        "dual adjacency context.");
        }

        TreeAltitudeAlgorithms::validateAltitudeBufferShape(tree, altitude);
        materializeAttributesFromBitquadFamilyCounts(
            tree, familyCounts, [&](NodeId nodeId) { return nodeUses4Connectivity(tree, altitude, nodeId); }, buffer, attrNames, requestedAttributes);
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

  private:
    /**
     * @brief Scalar counter layout needed by Duda-style formulas.
     */
    struct ScalarCounters {
        /** @brief Stores the count pattern c1 c4. */
        int countPatternC1C4 = 0;
        /** @brief Stores the count pattern c1. */
        int countPatternC1 = 0;
        /** @brief Stores the count pattern c2. */
        int countPatternC2 = 0;
        /** @brief Stores the count pattern cd. */
        int countPatternCD = 0;
        /** @brief Stores the count pattern c3. */
        int countPatternC3 = 0;
        /** @brief Stores the count pattern c4. */
        int countPatternC4 = 0;
    };

    /**
     * @brief Converts family counters to the scalar counter basis.
     *
     * @param counts Counters used by the operation.
     * @param is4Connectivity Flag controlling is4 connectivity.
     * @return The converted family counters to the scalar counter basis.
     */
    static ScalarCounters scalarCountersFromBitquadFamilyCounts(const BitquadFamilyCounts& counts, bool is4Connectivity) noexcept {
        ScalarCounters counters;
        if (is4Connectivity) {
            counters.countPatternC1C4 = counts.q1 + (2 * counts.qd);
        } else {
            counters.countPatternC1 = counts.q1;
            counters.countPatternCD = counts.qd;
        }

        counters.countPatternC2 = counts.q2;
        counters.countPatternC3 = counts.q3;
        counters.countPatternC4 = counts.q4;
        return counters;
    }

    /**
     * @brief Euler number from projected bitquad counters.
     *
     * @param counters Counters updated by the operation.
     * @return Euler number from projected bitquad counters.
     */
    static int numberEuler(const ScalarCounters& counters) noexcept { return (counters.countPatternC1C4 - counters.countPatternC3) / 4; }

    /**
     * @brief Number of holes from projected bitquad counters.
     *
     * @param counters Counters updated by the operation.
     * @return The number of holes from projected bitquad counters.
     */
    static int numberHoles(const ScalarCounters& counters) noexcept { return 1 - numberEuler(counters); }

    /**
     * @brief Discrete perimeter from projected bitquad counters.
     *
     * @param counters Counters updated by the operation.
     * @return Discrete perimeter from projected bitquad counters.
     */
    static int perimeter(const ScalarCounters& counters) noexcept {
        return counters.countPatternC1 + counters.countPatternC2 + counters.countPatternC3 + (2 * counters.countPatternCD);
    }

    /**
     * @brief Duda area estimator from projected bitquad counters.
     *
     * @param counters Counters updated by the operation.
     * @return Duda area estimator from projected bitquad counters.
     */
    static double areaDuda(const ScalarCounters& counters) noexcept {
        return (1.0 / 4.0 * counters.countPatternC1) + (1.0 / 2.0 * counters.countPatternC2) + (7.0 / 8.0 * counters.countPatternC3) + counters.countPatternC4 +
               (3.0 / 4.0 * counters.countPatternCD);
    }

    /**
     * @brief Continuous perimeter estimator from projected bitquad counters.
     *
     * @param counters Counters updated by the operation.
     * @return Continuous perimeter estimator from projected bitquad counters.
     */
    static double perimeterContinuous(const ScalarCounters& counters) noexcept {
        return counters.countPatternC2 + ((counters.countPatternC1 + counters.countPatternC3) / 1.5);
    }

    /**
     * @brief Circularity computed from Duda area and continuous perimeter.
     *
     * @param counters Counters updated by the operation.
     * @return Circularity computed from Duda area and continuous perimeter.
     */
    static double circularity(const ScalarCounters& counters) noexcept {
        const double area = areaDuda(counters);
        const double per = perimeterContinuous(counters);
        return ::mmcfilters::attributes::numeric::safeDivide(4.0 * std::numbers::pi * area, per * per);
    }

    /**
     * @brief Average area per connected component.
     *
     * @param counters Counters updated by the operation.
     * @return Average area per connected component.
     */
    static double areaAverage(const ScalarCounters& counters) noexcept {
        const int euler = numberEuler(counters);
        return ::mmcfilters::attributes::numeric::safeDivide(areaDuda(counters), static_cast<double>(euler));
    }

    /**
     * @brief Average continuous perimeter per connected component.
     *
     * @param counters Counters updated by the operation.
     * @return Average continuous perimeter per connected component.
     */
    static double perimeterAverage(const ScalarCounters& counters) noexcept {
        const int euler = numberEuler(counters);
        return ::mmcfilters::attributes::numeric::safeDivide(perimeterContinuous(counters), static_cast<double>(euler));
    }

    /**
     * @brief Average length proxy derived from average perimeter.
     *
     * @param counters Counters updated by the operation.
     * @return Average length proxy derived from average perimeter.
     */
    static double lengthAverage(const ScalarCounters& counters) noexcept {
        return ::mmcfilters::attributes::numeric::safeDivide(perimeterAverage(counters), 2.0);
    }

    /**
     * @brief Average width proxy derived from average area and perimeter.
     *
     * @param counters Counters updated by the operation.
     * @return Average width proxy derived from average area and perimeter.
     */
    static double widthAverage(const ScalarCounters& counters) noexcept {
        const double per = perimeterContinuous(counters);
        return ::mmcfilters::attributes::numeric::safeDivide(2.0 * areaDuda(counters), per);
    }

    /**
     * @brief Selects connectivity from a generic dual adjacency policy.
     *
     * @details
     * A non-root branch is classified by comparing its node altitude with the
     * parent altitude:
     *
     * - higher than parent: use the increasing adjacency;
     * - lower than parent: use the decreasing adjacency;
     * - equal to parent: use the common connectivity when both relations agree,
     *   otherwise reject because the branch direction is ambiguous.
     *
     * The root represents the whole image support. It uses the common
     * connectivity when the two relations agree; when they differ,
     * it falls back to the 8-connectivity scalar projection.
     *
     * @param tree Tree topology used by the operation.
     * @param altitude Altitude data indexed by node identifier.
     * @param nodeId Identifier of the node used by the operation.
     * @return The selected connectivity from a generic dual adjacency policy.
     */
    template <AltitudeValue T> static bool nodeUses4Connectivity(const MorphologicalTree& tree, std::span<const T> altitude, NodeId nodeId) {
        if (!tree.hasDirectionalGridAdjacency2D()) {
            throw std::invalid_argument("Local-event BitQuads scalar projection requires a dual "
                                        "adjacency policy.");
        }

        const RegularGridAdjacency2D* decreasingAdjacency = tree.getDecreasingGridAdjacency2D();
        const RegularGridAdjacency2D* increasingAdjacency = tree.getIncreasingGridAdjacency2D();
        if (decreasingAdjacency == nullptr || increasingAdjacency == nullptr) {
            throw std::invalid_argument("Local-event BitQuads scalar projection requires decreasing "
                                        "and increasing adjacency relations.");
        }
        if (!decreasingAdjacency->isCanonical4Or8Connectivity() || !increasingAdjacency->isCanonical4Or8Connectivity()) {
            throw std::invalid_argument("Local-event BitQuads scalar projection requires canonical 4- or 8-connectivity.");
        }

        const bool decreasingIs4Connectivity = decreasingAdjacency->is4connectivity();
        const bool increasingIs4Connectivity = increasingAdjacency->is4connectivity();
        if (tree.isRoot(nodeId)) {
            return decreasingIs4Connectivity == increasingIs4Connectivity ? decreasingIs4Connectivity : false;
        }

        const NodeId parentNodeId = tree.getNodeParent(nodeId);
        if (parentNodeId == InvalidNode || !tree.isAlive(parentNodeId)) {
            throw std::runtime_error("Local-event BitQuads scalar projection requires every "
                                     "non-root node to have an alive parent.");
        }

        const T nodeAltitude = TreeAltitudeAlgorithms::getAltitude(altitude, nodeId);
        const T parentAltitude = TreeAltitudeAlgorithms::getAltitude(altitude, parentNodeId);
        if (nodeAltitude > parentAltitude) {
            return increasingIs4Connectivity;
        }
        if (nodeAltitude < parentAltitude) {
            return decreasingIs4Connectivity;
        }
        if (decreasingIs4Connectivity == increasingIs4Connectivity) {
            return decreasingIs4Connectivity;
        }

        throw std::runtime_error("Local-event BitQuads scalar projection cannot select decreasing "
                                 "or increasing connectivity from equal node and parent altitudes.");
    }

    /**
     * @brief Shared scalar materialization once per-node connectivity is known.
     *
     * @details
     * `nodeUses4Connectivity(node)` selects the projection basis for each node.
     * The method writes only attributes requested by the caller, preserving the
     * standard flat-buffer attribute-kernel contract.
     *
     * @param tree Tree topology used by the operation.
     * @param familyCounts Per-node counters for the attribute family.
     * @param nodeUses4Connectivity Per-node selector for 4-connected topology.
     * @param buffer Buffer read or written by the operation.
     * @param attrNames Layout that maps attributes to buffer columns.
     * @param requestedAttributes Attributes requested for materialization.
     */
    template <std::floating_point Real, typename ConnectivitySelector>
    static void materializeAttributesFromBitquadFamilyCounts(const MorphologicalTree& tree, std::span<const BitquadFamilyCounts> familyCounts,
                                                             ConnectivitySelector nodeUses4Connectivity, std::span<Real> buffer,
                                                             const AttributeNames& attrNames, std::span<const Attribute> requestedAttributes) {
        auto indexOf = [&](int idx, Attribute attr) { return attrNames.linearIndex(idx, attr); };

        const bool computeArea = requestsAttribute(requestedAttributes, BITQUADS_AREA);
        const bool computeNumberEuler = requestsAttribute(requestedAttributes, BITQUADS_NUMBER_EULER);
        const bool computeNumberHoles = requestsAttribute(requestedAttributes, BITQUADS_NUMBER_HOLES);
        const bool computePerimeter = requestsAttribute(requestedAttributes, BITQUADS_PERIMETER);
        const bool computePerimeterCont = requestsAttribute(requestedAttributes, BITQUADS_PERIMETER_CONTINUOUS);
        const bool computeCircularity = requestsAttribute(requestedAttributes, BITQUADS_CIRCULARITY);
        const bool computePerimeterAverage = requestsAttribute(requestedAttributes, BITQUADS_PERIMETER_AVERAGE);
        const bool computeLengthAverage = requestsAttribute(requestedAttributes, BITQUADS_LENGTH_AVERAGE);
        const bool computeWidthAverage = requestsAttribute(requestedAttributes, BITQUADS_WIDTH_AVERAGE);

        for (NodeId nodeId : tree.getAliveNodeIds()) {
            const ScalarCounters counters =
                scalarCountersFromBitquadFamilyCounts(familyCounts[static_cast<std::size_t>(nodeId)], nodeUses4Connectivity(nodeId));
            if (computeArea) {
                buffer[indexOf(nodeId, BITQUADS_AREA)] = static_cast<Real>(areaDuda(counters));
            }
            if (computeNumberEuler) {
                buffer[indexOf(nodeId, BITQUADS_NUMBER_EULER)] = static_cast<Real>(numberEuler(counters));
            }
            if (computeNumberHoles) {
                buffer[indexOf(nodeId, BITQUADS_NUMBER_HOLES)] = static_cast<Real>(numberHoles(counters));
            }
            if (computePerimeter) {
                buffer[indexOf(nodeId, BITQUADS_PERIMETER)] = static_cast<Real>(perimeter(counters));
            }
            if (computePerimeterCont) {
                buffer[indexOf(nodeId, BITQUADS_PERIMETER_CONTINUOUS)] = static_cast<Real>(perimeterContinuous(counters));
            }
            if (computeCircularity) {
                buffer[indexOf(nodeId, BITQUADS_CIRCULARITY)] = static_cast<Real>(circularity(counters));
            }
            if (computePerimeterAverage) {
                buffer[indexOf(nodeId, BITQUADS_PERIMETER_AVERAGE)] = static_cast<Real>(perimeterAverage(counters));
            }
            if (computeLengthAverage) {
                buffer[indexOf(nodeId, BITQUADS_LENGTH_AVERAGE)] = static_cast<Real>(lengthAverage(counters));
            }
            if (computeWidthAverage) {
                buffer[indexOf(nodeId, BITQUADS_WIDTH_AVERAGE)] = static_cast<Real>(widthAverage(counters));
            }
        }
    }
};

} // namespace mmcfilters::attributes::computers::detail
