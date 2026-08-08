#pragma once

#include "HierarchySaliencyMap.hpp"
#include "../../utils/RegularGridAdjacency2D.hpp"
#include "../../utils/Common.hpp"
#include "../MorphologicalTree.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <sstream>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mmcfilters {

/**
 * @brief Selects whether lower or higher attribute plateaus are extrema.
 */
enum class ShapeSpaceExtremaPolarity { Minima, Maxima };

/**
 * @brief One regional extremum and its extinction interval in attribute space.
 */
template <std::floating_point Real> struct ShapeSpaceExtremum {
    /// Stable representative node of the regional extremum.
    NodeId representative = InvalidNode;
    /// Attribute level at which the extremum appears.
    Real birthLevel{};
    /// Attribute level at which the extremum is absorbed.
    Real deathLevel{};
    /// Absolute difference between death and birth levels.
    Real extinction{};
};

/**
 * @brief Regional extrema and their sparse dense-domain node scores.
 */
template <std::floating_point Real> struct ShapeSpaceExtinctionResult {
    /// Regional extrema ordered by the computation.
    std::vector<ShapeSpaceExtremum<Real>> extrema;
    /// Sparse extinction score stored in the dense node-id domain.
    std::vector<Real> nodeScores;
};

/**
 * @brief Extrema, their sparse node scores, and the projected contour map.
 */
template <std::floating_point Real> struct ShapeSpaceSaliencyResult {
    /// Regional extrema ordered by the computation.
    std::vector<ShapeSpaceExtremum<Real>> extrema;
    /// Sparse extinction score stored in the dense node-id domain.
    std::vector<Real> nodeScores;
    /// Image-adjacency edge map obtained from the node scores.
    EdgeSaliencyMap<Real> edgeMap;
};

/**
 * @brief Computes Xu-style extinction values in the shape space of a tree.
 *
 * The shape-space graph has one vertex per live morphological-tree node and one
 * undirected edge per parent-child relation. Regional extrema are the connected
 * equal-attribute plateaus born during a lower-level (minima) or upper-level
 * (maxima) component sweep. When components merge, only the strongest extremum
 * survives. Equal-strength extrema are resolved by increasing representative
 * NodeId. The dominant extremum dies at the opposite global attribute level.
 *
 * @par Primary reference
 * Yongchao Xu, Edwin Carlinet, Thierry Géraud, and Laurent Najman,
 * "Hierarchical Segmentation Using Tree-Based Shape Spaces," IEEE Transactions
 * on Pattern Analysis and Machine Intelligence, 39(3):457-469, 2017.
 * [DOI 10.1109/TPAMI.2016.2554550](https://doi.org/10.1109/TPAMI.2016.2554550).
 * Section 4.3 defines the extinction-based contour saliency construction. The
 * paper's local-minimum/Khalimsky-grid path is generalized here to minima or
 * maxima and to an edge-indexed regular-grid representation.
 */
class ShapeSpaceSaliency {
  private:
    /**
     * @brief Validates polarity.
     *
     * @param polarity Polarity that determines whether maxima or minima are processed.
     * @param context Operation name used in diagnostics.
     */
    static void validatePolarity(ShapeSpaceExtremaPolarity polarity, const char* context) {
        switch (polarity) {
        case ShapeSpaceExtremaPolarity::Minima:
        case ShapeSpaceExtremaPolarity::Maxima:
            return;
        }
        throw std::invalid_argument(std::string(context) + " received an unknown extrema polarity.");
    }

    /**
     * @brief Validates rooted tree.
     *
     * @param tree Tree topology used by the operation.
     * @param context Operation name used in diagnostics.
     */
    static void validateRootedTree(const MorphologicalTree& tree, const char* context) { detail::requireCommittedRootedHierarchy(tree, context); }

    /**
     * @brief Validates node buffer.
     *
     * @param tree Tree topology used by the operation.
     * @param values Values read or written by the operation.
     * @param requireNonNegative Whether negative values must be rejected.
     * @param valueName Value label included in validation error messages.
     * @param context Operation name used in diagnostics.
     */
    template <std::floating_point Real>
    static void validateNodeBuffer(const MorphologicalTree& tree, std::span<const Real> values, bool requireNonNegative, const char* valueName,
                                   const char* context) {
        const std::size_t expected = static_cast<std::size_t>(tree.getNumInternalNodeSlots());
        if (values.size() != expected) {
            std::ostringstream oss;
            oss << context << " requires one " << valueName << " value per dense internal node slot; expected " << expected << " values but got "
                << values.size() << ".";
            throw std::invalid_argument(oss.str());
        }

        for (NodeId nodeId : tree.getAliveNodeIds()) {
            const Real value = values[static_cast<std::size_t>(nodeId)];
            if (!std::isfinite(value)) {
                std::ostringstream oss;
                oss << context << " requires finite " << valueName << " values; node " << nodeId << " is non-finite.";
                throw std::invalid_argument(oss.str());
            }
            if (requireNonNegative && value < Real{0}) {
                std::ostringstream oss;
                oss << context << " requires non-negative " << valueName << " values; node " << nodeId << " has value " << value << ".";
                throw std::invalid_argument(oss.str());
            }
        }
    }

    /**
     * @brief Validates tree and adjacency.
     *
     * @param tree Tree topology used by the operation.
     * @param adjacency Adjacency relation used by the operation.
     * @param context Operation name used in diagnostics.
     */
    static void validateTreeAndAdjacency(const MorphologicalTree& tree, const RegularGridAdjacency2D& adjacency, const char* context) {
        validateRootedTree(tree, context);

        const int rows = tree.getNumRowsOfGridDomain2D();
        const int cols = tree.getNumColsOfGridDomain2D();
        const int numProperParts = tree.getNumTotalProperParts();
        if (rows <= 0 || cols <= 0 || numProperParts <= 0) {
            throw std::invalid_argument(std::string(context) + " requires a non-empty image/proper-part domain.");
        }

        const long long pixelCount = static_cast<long long>(rows) * static_cast<long long>(cols);
        if (pixelCount != static_cast<long long>(numProperParts)) {
            std::ostringstream oss;
            oss << context << " requires one proper part per image pixel; got " << numProperParts << " proper parts for a " << rows << "x" << cols
                << " image domain.";
            throw std::invalid_argument(oss.str());
        }

        if (adjacency.getNumRows() != rows || adjacency.getNumCols() != cols) {
            std::ostringstream oss;
            oss << context << " adjacency domain must match the tree image domain; got " << adjacency.getNumRows() << "x" << adjacency.getNumCols()
                << " for a tree domain of " << rows << "x" << cols << ".";
            throw std::invalid_argument(oss.str());
        }
    }

    /**
     * @brief Validates stored adjacency.
     *
     * @param tree Tree topology used by the operation.
     * @param context Operation name used in diagnostics.
     * @return Reference to the resulting object.
     */
    static const RegularGridAdjacency2D& requireStoredAdjacency(const MorphologicalTree& tree, const char* context) {
        const RegularGridAdjacency2D* adjacency = tree.getUniformGridAdjacency2D();
        if (adjacency == nullptr) {
            throw std::invalid_argument(std::string(context) + " requires an attached adjacency relation; pass an explicit adjacency relation instead.");
        }
        return *adjacency;
    }

    /**
     * @brief Checks and converts extinction.
     *
     * @param birthLevel Altitude or level represented by `birthLevel`.
     * @param deathLevel Altitude or level represented by `deathLevel`.
     * @param polarity Polarity that determines whether maxima or minima are processed.
     * @param representative Representative node or pixel used for the shape-space mapping.
     * @param context Operation name used in diagnostics.
     * @return Validated non-negative extinction value.
     */
    template <std::floating_point Real>
    static Real checkedExtinction(Real birthLevel, Real deathLevel, ShapeSpaceExtremaPolarity polarity, NodeId representative, const char* context) {
        const long double birth = static_cast<long double>(birthLevel);
        const long double death = static_cast<long double>(deathLevel);
        const long double difference = polarity == ShapeSpaceExtremaPolarity::Minima ? death - birth : birth - death;

        if (difference < 0.0L) {
            std::ostringstream oss;
            oss << context << " produced an invalid extinction interval for representative node " << representative << ".";
            throw std::invalid_argument(oss.str());
        }

        const long double maximum = static_cast<long double>(std::numeric_limits<Real>::max());
        if (!std::isfinite(difference) || difference > maximum) {
            std::ostringstream oss;
            oss << context << " extinction for representative node " << representative << " is not representable by the requested floating-point type.";
            throw std::overflow_error(oss.str());
        }

        const Real result = static_cast<Real>(difference);
        if (!std::isfinite(result)) {
            std::ostringstream oss;
            oss << context << " extinction for representative node " << representative << " overflowed the requested floating-point type.";
            throw std::overflow_error(oss.str());
        }
        return result;
    }

  public:
    /**
     * @brief Computes regional extrema and finite extinction values.
     *
     * `attribute` is indexed by the tree's dense internal NodeId domain. Birth
     * and death levels remain in the original attribute domain for both
     * polarities. Results are sorted by increasing representative NodeId.
     *
     * @param tree Tree topology used by the operation.
     * @param attribute Attribute requested by the operation.
     * @param polarity Extremum polarity selected by the operation.
     * @return The computed regional extrema and finite extinction values.
     */
    template <std::floating_point Real>
    [[nodiscard]] static ShapeSpaceExtinctionResult<Real> computeExtinctionValues(const MorphologicalTree& tree, std::span<const Real> attribute,
                                                                                  ShapeSpaceExtremaPolarity polarity) {
        constexpr const char* context = "ShapeSpaceSaliency::computeExtinctionValues";
        validatePolarity(polarity, context);
        validateRootedTree(tree, context);
        validateNodeBuffer(tree, attribute, false, "attribute", context);

        std::vector<NodeId> nodes;
        nodes.reserve(static_cast<std::size_t>(tree.getNumInternalNodeSlots()));
        for (NodeId nodeId : tree.getAliveNodeIds()) {
            nodes.push_back(nodeId);
        }
        if (nodes.empty()) {
            throw std::invalid_argument(std::string(context) + " requires at least one live node.");
        }

        std::sort(nodes.begin(), nodes.end(), [&](NodeId lhs, NodeId rhs) {
            const Real lhsLevel = attribute[static_cast<std::size_t>(lhs)];
            const Real rhsLevel = attribute[static_cast<std::size_t>(rhs)];
            if (lhsLevel == rhsLevel) {
                return lhs < rhs;
            }
            if (polarity == ShapeSpaceExtremaPolarity::Minima) {
                return lhsLevel < rhsLevel;
            }
            return rhsLevel < lhsLevel;
        });

        Real globalMinimum = attribute[static_cast<std::size_t>(nodes.front())];
        Real globalMaximum = globalMinimum;
        for (NodeId nodeId : nodes) {
            const Real level = attribute[static_cast<std::size_t>(nodeId)];
            globalMinimum = std::min(globalMinimum, level);
            globalMaximum = std::max(globalMaximum, level);
        }

        const int numSlots = tree.getNumInternalNodeSlots();
        const std::size_t slotCount = static_cast<std::size_t>(numSlots);

        std::vector<int> depth(slotCount, -1);
        std::vector<NodeId> stack;
        stack.push_back(tree.getRoot());
        depth[static_cast<std::size_t>(tree.getRoot())] = 0;
        while (!stack.empty()) {
            const NodeId nodeId = stack.back();
            stack.pop_back();
            for (NodeId childId : tree.getChildren(nodeId)) {
                depth[static_cast<std::size_t>(childId)] = depth[static_cast<std::size_t>(nodeId)] + 1;
                stack.push_back(childId);
            }
        }

        std::vector<NodeId> componentParent(slotCount, InvalidNode);
        std::vector<int> componentSize(slotCount, 0);
        std::vector<int> survivor(slotCount, -1);
        std::vector<std::uint8_t> active(slotCount, 0);
        std::vector<NodeId> plateauRepresentative(slotCount, InvalidNode);
        std::vector<std::vector<NodeId>> priorComponents(slotCount);

        auto findComponent = [&](NodeId nodeId) {
            NodeId root = nodeId;
            while (componentParent[static_cast<std::size_t>(root)] != root) {
                root = componentParent[static_cast<std::size_t>(root)];
            }
            while (componentParent[static_cast<std::size_t>(nodeId)] != nodeId) {
                const NodeId next = componentParent[static_cast<std::size_t>(nodeId)];
                componentParent[static_cast<std::size_t>(nodeId)] = root;
                nodeId = next;
            }
            return root;
        };

        auto joinComponents = [&](NodeId lhs, NodeId rhs) {
            lhs = findComponent(lhs);
            rhs = findComponent(rhs);
            if (lhs == rhs) {
                return lhs;
            }
            const int lhsSize = componentSize[static_cast<std::size_t>(lhs)];
            const int rhsSize = componentSize[static_cast<std::size_t>(rhs)];
            if (lhsSize < rhsSize || (lhsSize == rhsSize && rhs < lhs)) {
                std::swap(lhs, rhs);
            }
            componentParent[static_cast<std::size_t>(rhs)] = lhs;
            componentSize[static_cast<std::size_t>(lhs)] += componentSize[static_cast<std::size_t>(rhs)];
            return lhs;
        };

        auto forEachShapeNeighbor = [&](NodeId nodeId, auto&& visitor) {
            if (!tree.isRoot(nodeId)) {
                visitor(tree.getNodeParent(nodeId));
            }
            for (NodeId childId : tree.getChildren(nodeId)) {
                visitor(childId);
            }
        };

        std::vector<ShapeSpaceExtremum<Real>> extrema;
        std::vector<std::uint8_t> finalized;

        auto finishExtremum = [&](int extremumIndex, Real deathLevel) {
            if (extremumIndex < 0 || static_cast<std::size_t>(extremumIndex) >= extrema.size() || finalized[static_cast<std::size_t>(extremumIndex)] != 0) {
                throw std::runtime_error(std::string(context) + " encountered inconsistent component-extremum state.");
            }
            ShapeSpaceExtremum<Real>& extremum = extrema[static_cast<std::size_t>(extremumIndex)];
            extremum.deathLevel = deathLevel;
            extremum.extinction = checkedExtinction(extremum.birthLevel, deathLevel, polarity, extremum.representative, context);
            finalized[static_cast<std::size_t>(extremumIndex)] = 1;
        };

        auto isStronger = [&](int lhsIndex, int rhsIndex) {
            const auto& lhs = extrema[static_cast<std::size_t>(lhsIndex)];
            const auto& rhs = extrema[static_cast<std::size_t>(rhsIndex)];
            if (lhs.birthLevel != rhs.birthLevel) {
                if (polarity == ShapeSpaceExtremaPolarity::Minima) {
                    return lhs.birthLevel < rhs.birthLevel;
                }
                return rhs.birthLevel < lhs.birthLevel;
            }
            return lhs.representative < rhs.representative;
        };

        std::size_t batchBegin = 0;
        while (batchBegin < nodes.size()) {
            const Real level = attribute[static_cast<std::size_t>(nodes[batchBegin])];
            std::size_t batchEnd = batchBegin + 1;
            while (batchEnd < nodes.size() && attribute[static_cast<std::size_t>(nodes[batchEnd])] == level) {
                ++batchEnd;
            }

            for (std::size_t i = batchBegin; i < batchEnd; ++i) {
                const NodeId nodeId = nodes[i];
                const std::size_t index = static_cast<std::size_t>(nodeId);
                active[index] = 1;
                componentParent[index] = nodeId;
                componentSize[index] = 1;
                survivor[index] = -1;
            }

            for (std::size_t i = batchBegin; i < batchEnd; ++i) {
                const NodeId nodeId = nodes[i];
                forEachShapeNeighbor(nodeId, [&](NodeId neighborId) {
                    if (active[static_cast<std::size_t>(neighborId)] != 0 && attribute[static_cast<std::size_t>(neighborId)] == level) {
                        static_cast<void>(joinComponents(nodeId, neighborId));
                    }
                });
            }

            std::vector<NodeId> plateauRoots;
            plateauRoots.reserve(batchEnd - batchBegin);
            for (std::size_t i = batchBegin; i < batchEnd; ++i) {
                const NodeId nodeId = nodes[i];
                const NodeId plateauRoot = findComponent(nodeId);
                plateauRoots.push_back(plateauRoot);

                NodeId& representative = plateauRepresentative[static_cast<std::size_t>(plateauRoot)];
                if (representative == InvalidNode || depth[static_cast<std::size_t>(nodeId)] < depth[static_cast<std::size_t>(representative)] ||
                    (depth[static_cast<std::size_t>(nodeId)] == depth[static_cast<std::size_t>(representative)] && nodeId < representative)) {
                    representative = nodeId;
                }

                forEachShapeNeighbor(nodeId, [&](NodeId neighborId) {
                    if (active[static_cast<std::size_t>(neighborId)] != 0 && attribute[static_cast<std::size_t>(neighborId)] != level) {
                        priorComponents[static_cast<std::size_t>(plateauRoot)].push_back(findComponent(neighborId));
                    }
                });
            }

            std::sort(plateauRoots.begin(), plateauRoots.end());
            plateauRoots.erase(std::unique(plateauRoots.begin(), plateauRoots.end()), plateauRoots.end());

            for (NodeId plateauRoot : plateauRoots) {
                auto& adjacentComponents = priorComponents[static_cast<std::size_t>(plateauRoot)];
                for (NodeId& component : adjacentComponents) {
                    component = findComponent(component);
                }
                std::sort(adjacentComponents.begin(), adjacentComponents.end());
                adjacentComponents.erase(std::unique(adjacentComponents.begin(), adjacentComponents.end()), adjacentComponents.end());

                int winningExtremum = -1;
                if (adjacentComponents.empty()) {
                    const NodeId representative = plateauRepresentative[static_cast<std::size_t>(plateauRoot)];
                    winningExtremum = static_cast<int>(extrema.size());
                    extrema.push_back(ShapeSpaceExtremum<Real>{representative, level, level, Real{0}});
                    finalized.push_back(0);
                } else {
                    for (NodeId component : adjacentComponents) {
                        const int candidate = survivor[static_cast<std::size_t>(component)];
                        if (candidate < 0) {
                            throw std::runtime_error(std::string(context) + " found an active level component without a surviving extremum.");
                        }
                        if (winningExtremum < 0 || isStronger(candidate, winningExtremum)) {
                            winningExtremum = candidate;
                        }
                    }

                    for (NodeId component : adjacentComponents) {
                        const int candidate = survivor[static_cast<std::size_t>(component)];
                        if (candidate != winningExtremum) {
                            finishExtremum(candidate, level);
                        }
                    }
                }

                NodeId combinedRoot = plateauRoot;
                for (NodeId component : adjacentComponents) {
                    combinedRoot = joinComponents(combinedRoot, component);
                }
                survivor[static_cast<std::size_t>(combinedRoot)] = winningExtremum;

                plateauRepresentative[static_cast<std::size_t>(plateauRoot)] = InvalidNode;
                adjacentComponents.clear();
            }

            batchBegin = batchEnd;
        }

        const NodeId finalComponent = findComponent(nodes.front());
        for (NodeId nodeId : nodes) {
            if (findComponent(nodeId) != finalComponent) {
                throw std::runtime_error(std::string(context) + " did not produce one connected final level component.");
            }
        }

        const int dominantExtremum = survivor[static_cast<std::size_t>(finalComponent)];
        if (dominantExtremum < 0) {
            throw std::runtime_error(std::string(context) + " did not retain a dominant extremum.");
        }
        finishExtremum(dominantExtremum, polarity == ShapeSpaceExtremaPolarity::Minima ? globalMaximum : globalMinimum);

        for (std::uint8_t isFinalized : finalized) {
            if (isFinalized == 0) {
                throw std::runtime_error(std::string(context) + " left an extremum without a death level.");
            }
        }

        std::sort(extrema.begin(), extrema.end(), [](const auto& lhs, const auto& rhs) { return lhs.representative < rhs.representative; });

        ShapeSpaceExtinctionResult<Real> result;
        result.extrema = std::move(extrema);
        result.nodeScores.assign(slotCount, Real{0});
        for (const ShapeSpaceExtremum<Real>& extremum : result.extrema) {
            result.nodeScores[static_cast<std::size_t>(extremum.representative)] = extremum.extinction;
        }
        return result;
    }

    /**
     * @brief Projects sparse node scores onto every image-domain adjacency edge.
     *
     * For an edge `(p, q)`, the result is the maximum score on the two paths
     * `owner(p) -> LCA` and `owner(q) -> LCA`, excluding the LCA. These are
     * exactly the nodes whose regions contain that edge in their contour.
     *
     * @param tree Tree topology used by the operation.
     * @param nodeScores Per-node scores used by the operation.
     * @param adjacency Adjacency relation used by the operation.
     * @return The projected sparse node scores onto every image-domain adjacency edge.
     */
    template <std::floating_point Real>
    [[nodiscard]] static EdgeSaliencyMap<Real> projectContourScores(const MorphologicalTree& tree, std::span<const Real> nodeScores,
                                                                    const RegularGridAdjacency2D& adjacency) {
        constexpr const char* context = "ShapeSpaceSaliency::projectContourScores";
        validateTreeAndAdjacency(tree, adjacency, context);
        validateNodeBuffer(tree, nodeScores, true, "node-score", context);

        EdgeSaliencyMap<Real> edgeMap;
        edgeMap.numRows = tree.getNumRowsOfGridDomain2D();
        edgeMap.numCols = tree.getNumColsOfGridDomain2D();
        edgeMap.adjacencyRadius = adjacency.getRadius();

        const std::size_t slotCount = static_cast<std::size_t>(tree.getNumInternalNodeSlots());
        const std::size_t liftingLevelCount = std::max<std::size_t>(1, static_cast<std::size_t>(std::bit_width(slotCount - 1)));

        std::vector<int> depth(slotCount, -1);
        std::vector<std::vector<NodeId>> ancestor(liftingLevelCount, std::vector<NodeId>(slotCount, InvalidNode));
        std::vector<std::vector<Real>> pathMaximum(liftingLevelCount, std::vector<Real>(slotCount, Real{0}));

        std::vector<NodeId> stack{tree.getRoot()};
        depth[static_cast<std::size_t>(tree.getRoot())] = 0;
        while (!stack.empty()) {
            const NodeId nodeId = stack.back();
            stack.pop_back();
            const std::size_t nodeIndex = static_cast<std::size_t>(nodeId);
            ancestor[0][nodeIndex] = tree.getNodeParent(nodeId);
            pathMaximum[0][nodeIndex] = nodeScores[nodeIndex];

            for (NodeId childId : tree.getChildren(nodeId)) {
                depth[static_cast<std::size_t>(childId)] = depth[nodeIndex] + 1;
                stack.push_back(childId);
            }
        }

        for (std::size_t level = 1; level < liftingLevelCount; ++level) {
            for (NodeId nodeId : tree.getAliveNodeIds()) {
                const std::size_t nodeIndex = static_cast<std::size_t>(nodeId);
                const NodeId middle = ancestor[level - 1][nodeIndex];
                const std::size_t middleIndex = static_cast<std::size_t>(middle);
                ancestor[level][nodeIndex] = ancestor[level - 1][middleIndex];
                pathMaximum[level][nodeIndex] = std::max(pathMaximum[level - 1][nodeIndex], pathMaximum[level - 1][middleIndex]);
            }
        }

        auto accumulateBranch = [&](NodeId owner, NodeId lca, Real& edgeScore) {
            const int ownerDepth = depth[static_cast<std::size_t>(owner)];
            const int lcaDepth = depth[static_cast<std::size_t>(lca)];
            if (ownerDepth < lcaDepth) {
                throw std::runtime_error(std::string(context) + " found an LCA below an edge-endpoint owner.");
            }

            NodeId current = owner;
            std::size_t remaining = static_cast<std::size_t>(ownerDepth - lcaDepth);
            std::size_t level = 0;
            while (remaining != 0) {
                if ((remaining & std::size_t{1}) != 0) {
                    const std::size_t currentIndex = static_cast<std::size_t>(current);
                    edgeScore = std::max(edgeScore, pathMaximum[level][currentIndex]);
                    current = ancestor[level][currentIndex];
                }
                remaining >>= 1;
                ++level;
            }
            if (current != lca) {
                throw std::runtime_error(std::string(context) + " encountered a branch that does not reach its LCA.");
            }
        };

        const NodeId numProperParts = tree.getNumTotalProperParts();
        for (NodeId source = 0; source < numProperParts; ++source) {
            const NodeId sourceOwner = tree.getProperPartOwner(source);
            if (!tree.isAlive(sourceOwner)) {
                throw std::runtime_error(std::string(context) + " found a proper part without a live owner.");
            }

            for (int targetValue : adjacency.getForwardNeighborIndices(source)) {
                const NodeId target = static_cast<NodeId>(targetValue);
                if (target < 0 || target >= numProperParts) {
                    throw std::runtime_error(std::string(context) + " adjacency produced an endpoint outside the proper-part domain.");
                }

                const NodeId targetOwner = tree.getProperPartOwner(target);
                if (!tree.isAlive(targetOwner)) {
                    throw std::runtime_error(std::string(context) + " found a neighbour proper part without a live owner.");
                }

                const NodeId lca = tree.getLowestCommonAncestor(sourceOwner, targetOwner);
                if (lca == InvalidNode || !tree.isAlive(lca)) {
                    throw std::runtime_error(std::string(context) + " could not find a live LCA for an adjacency edge.");
                }

                Real edgeScore{0};
                accumulateBranch(sourceOwner, lca, edgeScore);
                accumulateBranch(targetOwner, lca, edgeScore);

                edgeMap.sources.push_back(source);
                edgeMap.targets.push_back(target);
                edgeMap.values.push_back(edgeScore);
            }
        }

        return edgeMap;
    }

    /**
     * @brief Projects node scores using the adjacency stored by the tree.
     *
     * @param tree Tree topology used by the operation.
     * @param nodeScores Per-node scores used by the operation.
     * @return The projected node scores using the adjacency stored by the tree.
     */
    template <std::floating_point Real>
    [[nodiscard]] static EdgeSaliencyMap<Real> projectContourScores(const MorphologicalTree& tree, std::span<const Real> nodeScores) {
        return projectContourScores(tree, nodeScores, requireStoredAdjacency(tree, "ShapeSpaceSaliency::projectContourScores"));
    }

    /**
     * @brief Computes extinction values, sparse representative scores, and contours.
     *
     * @param tree Tree topology used by the operation.
     * @param attribute Attribute requested by the operation.
     * @param polarity Extremum polarity selected by the operation.
     * @param adjacency Adjacency relation used by the operation.
     * @return The computed extinction values, sparse representative scores, and contours.
     */
    template <std::floating_point Real>
    [[nodiscard]] static ShapeSpaceSaliencyResult<Real> compute(const MorphologicalTree& tree, std::span<const Real> attribute,
                                                                ShapeSpaceExtremaPolarity polarity, const RegularGridAdjacency2D& adjacency) {
        ShapeSpaceSaliencyResult<Real> result;
        ShapeSpaceExtinctionResult<Real> extinction = computeExtinctionValues(tree, attribute, polarity);
        result.extrema = std::move(extinction.extrema);
        result.nodeScores = std::move(extinction.nodeScores);
        result.edgeMap = projectContourScores(tree, std::span<const Real>(result.nodeScores), adjacency);
        return result;
    }

    /**
     * @brief Computes the complete result using the adjacency stored by the tree.
     *
     * @param tree Tree topology used by the operation.
     * @param attribute Attribute requested by the operation.
     * @param polarity Extremum polarity selected by the operation.
     * @return The computed complete result using the adjacency stored by the tree.
     */
    template <std::floating_point Real>
    [[nodiscard]] static ShapeSpaceSaliencyResult<Real> compute(const MorphologicalTree& tree, std::span<const Real> attribute,
                                                                ShapeSpaceExtremaPolarity polarity) {
        return compute(tree, attribute, polarity, requireStoredAdjacency(tree, "ShapeSpaceSaliency::compute"));
    }
};

} // namespace mmcfilters
