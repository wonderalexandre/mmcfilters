#pragma once

#include "../AttributeComputer.hpp"
#include "BitquadAttributeData.hpp"
#include "detail/BitquadLocalEventComputation.hpp"
#include "../../trees/TreeAltitudeAlgorithms.hpp"
#include "../../trees/WeightedTreeView.hpp"

#include <algorithm>
#include <limits>
#include <numbers>
#include <span>
#include <stdexcept>
#include <vector>

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
class BitquadAttributeComputer : public AttributeComputer {
public:
    using AttributeComputer::compute;
    using AttributeComputer::computeUnitAttributes;

    /**
     * @brief Projected bitquad-family bucket consumed by this scalar computer.
     */
    using BitquadFamilyCounts = ::mmcfilters::attributes::computers::BitquadFamilyCounts;

    /**
     * @brief Returns the scalar bitquad attributes materialized by this computer.
     */
    [[nodiscard]] std::vector<Attribute> attributes() const override {
        return {
            BITQUADS_AREA,
            BITQUADS_NUMBER_EULER,
            BITQUADS_NUMBER_HOLES,
            BITQUADS_PERIMETER,
            BITQUADS_PERIMETER_CONTINUOUS,
            BITQUADS_CIRCULARITY,
            BITQUADS_PERIMETER_AVERAGE,
            BITQUADS_LENGTH_AVERAGE,
            BITQUADS_WIDTH_AVERAGE,
        };
    }

    /**
     * @brief Computes requested scalar bitquad attributes for every live node.
     *
     * @details
     * This entry point first computes family deltas, aggregates them into final
     * per-node family counts, and then projects only the requested scalar
     * descriptors into the caller-provided flat buffer. `altitude` may be empty
     * for ordinary component trees with an adjacency relation. It is required
     * for Tree of Shapes scalar projection.
     */
    void compute(const MorphologicalTree& tree, AttributeAltitudeView altitude, std::span<float> buffer, const AttributeNames& attrNames, std::span<const Attribute> requestedAttributes, std::span<const DependencySource>) const override {
        requireAttributeBufferShape(tree, buffer, attrNames);
        const auto familyDeltas = detail::BitquadLocalEventComputation::computeBitquadFamilyDeltas(tree);
        const auto familyCounts = detail::BitquadLocalEventComputation::aggregateBitquadFamilyDeltas(tree, familyDeltas);
        materializeAttributesFromBitquadFamilyCounts(tree, altitude, familyCounts, buffer, attrNames, requestedAttributes);
    }

    /**
     * @brief Projects family counters for component-tree inputs.
     *
     * This overload requires the tree to carry a regular adjacency relation.
     * Tree of Shapes scalar projection should use an altitude-aware overload.
     */
    static void materializeAttributesFromBitquadFamilyCounts(const MorphologicalTree& tree, std::span<const BitquadFamilyCounts> familyCounts, std::span<float> buffer, const AttributeNames& attrNames, std::span<const Attribute> requestedAttributes) {
        materializeAttributesFromBitquadFamilyCounts(tree, AttributeAltitudeView{}, familyCounts, buffer, attrNames, requestedAttributes);
    }

    /**
     * @brief Projects family counters into scalar bitquad attributes.
     *
     * @details
     * When the tree has a regular adjacency relation, all nodes use the same
     * 4-connectivity or 8-connectivity projection. Without a regular adjacency
     * relation, only Tree of Shapes inputs are accepted, and `altitude` is used
     * to select the min-tree or max-tree auxiliary adjacency per node.
     *
     * @throws std::invalid_argument If family counts do not cover every
     * internal node slot, if a non-ToS tree lacks adjacency metadata, or if ToS
     * projection is requested without altitude.
     */
    static void materializeAttributesFromBitquadFamilyCounts(const MorphologicalTree& tree, AttributeAltitudeView altitude, std::span<const BitquadFamilyCounts> familyCounts, std::span<float> buffer, const AttributeNames& attrNames, std::span<const Attribute> requestedAttributes) {
        if (altitude.has_value()) {
            materializeAttributesFromBitquadFamilyCounts(tree, *altitude, familyCounts, buffer, attrNames, requestedAttributes);
            return;
        }

        const std::size_t numNodeSlots = static_cast<std::size_t>(tree.getNumInternalNodeSlots());
        if (familyCounts.size() < numNodeSlots) {
            throw std::invalid_argument("Local-event bitquad family counts do not cover all tree node slots.");
        }

        const AdjacencyRelation* adjacency = tree.getAdjacencyRelation();
        if (adjacency != nullptr) {
            materializeAttributesFromBitquadFamilyCounts(tree, familyCounts,
                [is4Connectivity = adjacency->is4connectivity()](NodeId) {
                    return is4Connectivity;
                },
                buffer,
                attrNames,
                requestedAttributes);
            return;
        }

        if (tree.getTreeType() != MorphologicalTreeKind::TREE_OF_SHAPES) {
            throw std::invalid_argument("Local-event BitQuads scalar attributes require an adjacency relation.");
        }

        throw std::invalid_argument("Local-event ToS BitQuads scalar projection requires an altitude buffer.");
    }

    /**
     * @brief Projects family counters using a generic altitude span.
     *
     * @details
     * Regular component trees ignore `altitude` because their scalar projection
     * uses the tree adjacency relation. Tree of Shapes inputs validate and read
     * the span to infer each non-root node polarity by comparing the node altitude
     * with its parent's altitude.
     *
     * @throws std::invalid_argument If family counts or altitude shape are
     * invalid, if a non-ToS tree lacks adjacency metadata, or if ToS auxiliary
     * adjacency metadata is missing.
     */
    template<AltitudeValue T>
    static void materializeAttributesFromBitquadFamilyCounts(const MorphologicalTree& tree, std::span<const T> altitude, std::span<const BitquadFamilyCounts> familyCounts, std::span<float> buffer, const AttributeNames& attrNames, std::span<const Attribute> requestedAttributes) {
        const std::size_t numNodeSlots = static_cast<std::size_t>(tree.getNumInternalNodeSlots());
        if (familyCounts.size() < numNodeSlots) {
            throw std::invalid_argument("Local-event bitquad family counts do not cover all tree node slots.");
        }

        const AdjacencyRelation* adjacency = tree.getAdjacencyRelation();
        if (adjacency != nullptr) {
            materializeAttributesFromBitquadFamilyCounts(
                tree,
                familyCounts,
                [is4Connectivity = adjacency->is4connectivity()](NodeId) {
                    return is4Connectivity;
                },
                buffer,
                attrNames,
                requestedAttributes);
            return;
        }

        if (tree.getTreeType() != MorphologicalTreeKind::TREE_OF_SHAPES) {
            throw std::invalid_argument("Local-event BitQuads scalar attributes require an adjacency relation.");
        }

        TreeAltitudeAlgorithms::validateAltitudeBufferShape(tree, altitude);
        materializeAttributesFromBitquadFamilyCounts(
            tree,
            familyCounts,
            [&](NodeId nodeId) {
                return treeOfShapesNodeUses4Connectivity(tree, altitude, nodeId);
            },
            buffer,
            attrNames,
            requestedAttributes);
    }

    /**
     * @brief Projects family counters using a non-owning weighted-tree view.
     *
     * The view must remain valid for the duration of the call and its altitude
     * span must cover the topology's dense internal node-slot domain.
     */
    template<AltitudeValue T>
    static void materializeAttributesFromBitquadFamilyCounts(
        const WeightedTreeView<T>& tree,
        std::span<const BitquadFamilyCounts> familyCounts,
        std::span<float> buffer,
        const AttributeNames& attrNames,
        std::span<const Attribute> requestedAttributes) {
        materializeAttributesFromBitquadFamilyCounts(tree.topology(), tree.altitude(), familyCounts, buffer, attrNames, requestedAttributes);
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
    void computeUnitAttributes(
        const MorphologicalTree& tree,
        AttributeAltitudeView,
        std::span<const NodeId> unitProperParts,
        std::span<float> buffer,
        const AttributeNames& attrNames,
        std::span<const Attribute> requestedAttributes) const override {
        requireUnitAttributeBufferShape(tree, unitProperParts, buffer, attrNames);
        const AdjacencyRelation* adjacency = tree.getAdjacencyRelation();
        if (adjacency == nullptr) {
            throw std::invalid_argument("Local-event BitQuads attributes require an adjacency relation.");
        }

        const bool is4Connectivity = adjacency->is4connectivity();
        const float nan = std::numeric_limits<float>::quiet_NaN();
        const float inf = std::numeric_limits<float>::infinity();

        auto unitValue = [&](Attribute attribute) -> float {
            if (is4Connectivity) {
                switch (attribute) {
                    case BITQUADS_AREA: return 0.0f;
                    case BITQUADS_NUMBER_EULER: return 1.0f;
                    case BITQUADS_NUMBER_HOLES: return 0.0f;
                    case BITQUADS_PERIMETER: return 0.0f;
                    case BITQUADS_PERIMETER_CONTINUOUS: return 0.0f;
                    case BITQUADS_CIRCULARITY: return nan;
                    case BITQUADS_PERIMETER_AVERAGE: return 0.0f;
                    case BITQUADS_LENGTH_AVERAGE: return 0.0f;
                    case BITQUADS_WIDTH_AVERAGE: return nan;
                    default: break;
                }
            }

            switch (attribute) {
                case BITQUADS_AREA: return 1.0f;
                case BITQUADS_NUMBER_EULER: return 0.0f;
                case BITQUADS_NUMBER_HOLES: return 1.0f;
                case BITQUADS_PERIMETER: return 4.0f;
                case BITQUADS_PERIMETER_CONTINUOUS: return 8.0f / 3.0f;
                case BITQUADS_CIRCULARITY: return 9.0f * std::numbers::pi_v<float> / 16.0f;
                case BITQUADS_PERIMETER_AVERAGE: return inf;
                case BITQUADS_LENGTH_AVERAGE: return inf;
                case BITQUADS_WIDTH_AVERAGE: return nan;
                default: break;
            }
            throw std::runtime_error("Unsupported local-event BitQuads unit attribute.");
        };

        for (Attribute attribute : attributes()) {
            if (!requestsAttribute(requestedAttributes, attribute)) {
                continue;
            }
            const float value = unitValue(attribute);
            for (NodeId leafIndex = 0; leafIndex < static_cast<NodeId>(unitProperParts.size()); ++leafIndex) {
                buffer[attrNames.linearIndex(leafIndex, attribute)] = value;
            }
        }
    }

private:
    /**
     * @brief Scalar counter layout needed by Duda-style formulas.
     */
    struct ScalarCounters {
        int countPatternC1C4 = 0;
        int countPatternC1 = 0;
        int countPatternC2 = 0;
        int countPatternCD = 0;
        int countPatternC3 = 0;
        int countPatternC4 = 0;
    };

    /**
     * @brief Converts family counters to the scalar counter basis.
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
     */
    static int numberEuler(const ScalarCounters& counters) noexcept {
        return (counters.countPatternC1C4 - counters.countPatternC3) / 4;
    }

    /**
     * @brief Number of holes from projected bitquad counters.
     */
    static int numberHoles(const ScalarCounters& counters) noexcept {
        return 1 - numberEuler(counters);
    }

    /**
     * @brief Discrete perimeter from projected bitquad counters.
     */
    static int perimeter(const ScalarCounters& counters) noexcept {
        return counters.countPatternC1 + counters.countPatternC2 + counters.countPatternC3 + (2 * counters.countPatternCD);
    }

    /**
     * @brief Duda area estimator from projected bitquad counters.
     */
    static double areaDuda(const ScalarCounters& counters) noexcept {
        return (1.0 / 4.0 * counters.countPatternC1) +
               (1.0 / 2.0 * counters.countPatternC2) +
               (7.0 / 8.0 * counters.countPatternC3) +
               counters.countPatternC4 +
               (3.0 / 4.0 * counters.countPatternCD);
    }

    /**
     * @brief Continuous perimeter estimator from projected bitquad counters.
     */
    static double perimeterContinuous(const ScalarCounters& counters) noexcept {
        return counters.countPatternC2 + ((counters.countPatternC1 + counters.countPatternC3) / 1.5);
    }

    /**
     * @brief Circularity computed from Duda area and continuous perimeter.
     */
    static double circularity(const ScalarCounters& counters) noexcept {
        const double area = areaDuda(counters);
        const double per = perimeterContinuous(counters);
        return (4.0 * std::numbers::pi * area) / (per * per);
    }

    /**
     * @brief Average area per connected component.
     */
    static double areaAverage(const ScalarCounters& counters) noexcept {
        return areaDuda(counters) / static_cast<double>(numberEuler(counters));
    }

    /**
     * @brief Average continuous perimeter per connected component.
     */
    static double perimeterAverage(const ScalarCounters& counters) noexcept {
        return perimeterContinuous(counters) / static_cast<double>(numberEuler(counters));
    }

    /**
     * @brief Average length proxy derived from average perimeter.
     */
    static double lengthAverage(const ScalarCounters& counters) noexcept {
        return perimeterAverage(counters) / 2.0;
    }

    /**
     * @brief Average width proxy derived from average area and perimeter.
     */
    static double widthAverage(const ScalarCounters& counters) noexcept {
        return (2.0 * areaAverage(counters)) / perimeterAverage(counters);
    }

    /**
     * @brief Selects the ToS node connectivity used for scalar projection.
     *
     * @details
     * The Tree of Shapes topology alternates nodes that correspond to max-tree
     * and min-tree events. A non-root node is classified by comparing its
     * altitude with its parent's altitude:
     *
     * - higher than parent: use the auxiliary max-tree adjacency;
     * - lower than parent: use the auxiliary min-tree adjacency;
     * - equal to parent: reject, because polarity is not defined.
     *
     * The root represents the whole image support. It uses the common
     * connectivity when min/max auxiliary adjacencies agree; when they differ,
     * it falls back to the 8-connectivity scalar projection.
     */
    template<AltitudeValue T>
    static bool treeOfShapesNodeUses4Connectivity(
        const MorphologicalTree& tree,
        std::span<const T> altitude,
        NodeId nodeId) {
        if (!tree.hasTreeOfShapesAdjacencyPolicy()) {
            throw std::invalid_argument("Local-event ToS BitQuads scalar projection requires min/max adjacency metadata.");
        }

        const AdjacencyRelation* minAdjacency = tree.getTreeOfShapesMinTreeAdjacencyRelation();
        const AdjacencyRelation* maxAdjacency = tree.getTreeOfShapesMaxTreeAdjacencyRelation();
        if (minAdjacency == nullptr || maxAdjacency == nullptr) {
            throw std::invalid_argument("Local-event ToS BitQuads scalar projection requires min/max adjacency relations.");
        }

        if (tree.isRoot(nodeId)) {
            const bool minIs4Connectivity = minAdjacency->is4connectivity();
            const bool maxIs4Connectivity = maxAdjacency->is4connectivity();
            return minIs4Connectivity == maxIs4Connectivity ? minIs4Connectivity : false;
        }

        const NodeId parentNodeId = tree.getNodeParent(nodeId);
        if (parentNodeId == InvalidNode || !tree.isAlive(parentNodeId)) {
            throw std::runtime_error("Local-event ToS BitQuads scalar projection requires every non-root node to have an alive parent.");
        }

        const T nodeAltitude = TreeAltitudeAlgorithms::getAltitude(altitude, nodeId);
        const T parentAltitude = TreeAltitudeAlgorithms::getAltitude(altitude, parentNodeId);
        if (nodeAltitude > parentAltitude) {
            return maxAdjacency->is4connectivity();
        }
        if (nodeAltitude < parentAltitude) {
            return minAdjacency->is4connectivity();
        }

        throw std::runtime_error("Local-event ToS BitQuads scalar projection cannot infer node min/max polarity from equal node and parent altitudes.");
    }

    /**
     * @brief Shared scalar materialization once per-node connectivity is known.
     *
     * @details
     * `nodeUses4Connectivity(node)` selects the projection basis for each node.
     * The method writes only attributes requested by the caller, preserving the
     * standard `AttributeComputer` flat-buffer contract.
     */
    template <typename ConnectivitySelector>
    static void materializeAttributesFromBitquadFamilyCounts(
        const MorphologicalTree& tree,
        std::span<const BitquadFamilyCounts> familyCounts,
        ConnectivitySelector nodeUses4Connectivity,
        std::span<float> buffer,
        const AttributeNames& attrNames,
        std::span<const Attribute> requestedAttributes) {
        auto indexOf = [&](int idx, Attribute attr) {
            return attrNames.linearIndex(idx, attr);
        };

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
            const ScalarCounters counters = scalarCountersFromBitquadFamilyCounts(
                familyCounts[static_cast<std::size_t>(nodeId)],
                nodeUses4Connectivity(nodeId));
            if (computeArea) {
                buffer[indexOf(nodeId, BITQUADS_AREA)] = static_cast<float>(areaDuda(counters));
            }
            if (computeNumberEuler) {
                buffer[indexOf(nodeId, BITQUADS_NUMBER_EULER)] = static_cast<float>(numberEuler(counters));
            }
            if (computeNumberHoles) {
                buffer[indexOf(nodeId, BITQUADS_NUMBER_HOLES)] = static_cast<float>(numberHoles(counters));
            }
            if (computePerimeter) {
                buffer[indexOf(nodeId, BITQUADS_PERIMETER)] = static_cast<float>(perimeter(counters));
            }
            if (computePerimeterCont) {
                buffer[indexOf(nodeId, BITQUADS_PERIMETER_CONTINUOUS)] = static_cast<float>(perimeterContinuous(counters));
            }
            if (computeCircularity) {
                buffer[indexOf(nodeId, BITQUADS_CIRCULARITY)] = static_cast<float>(circularity(counters));
            }
            if (computePerimeterAverage) {
                buffer[indexOf(nodeId, BITQUADS_PERIMETER_AVERAGE)] = static_cast<float>(perimeterAverage(counters));
            }
            if (computeLengthAverage) {
                buffer[indexOf(nodeId, BITQUADS_LENGTH_AVERAGE)] = static_cast<float>(lengthAverage(counters));
            }
            if (computeWidthAverage) {
                buffer[indexOf(nodeId, BITQUADS_WIDTH_AVERAGE)] = static_cast<float>(widthAverage(counters));
            }
        }
    }
};

} // namespace mmcfilters::attributes::computers
