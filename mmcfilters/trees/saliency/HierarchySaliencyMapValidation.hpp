#pragma once

#include "../../utils/Altitude.hpp"
#include "../../utils/Common.hpp"
#include "../../utils/RegularGridAdjacency2D.hpp"
#include "../MorphologicalTree.hpp"
#include "../detail/TreeKernelAccess.hpp"
#include "../../utils/CommittedGridAccess.hpp"
#include "../ValuedMorphologicalTree.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <sstream>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace mmcfilters {

namespace detail {

/**
 * @brief Requires a committed non-empty hierarchy before projection.
 *
 * Construction and successful edit commits establish the rooted-tree and
 * support invariants. Reject an open editor, which may temporarily expose a
 * forest, without traversing and revalidating the committed topology.
 *
 * @param tree Tree topology.
 * @param context Operation context or diagnostic label.
 */
inline void requireCommittedRootedHierarchy(const MorphologicalTree& tree, const char* context) {
    if (tree.isEditing()) {
        throw std::invalid_argument(std::string(context) + " requires a committed tree; an edit session is still open.");
    }
    if (tree.root() == InvalidNode || !tree.isAlive(tree.root())) {
        throw std::invalid_argument(std::string(context) + " requires a non-empty connected rooted tree.");
    }
}

} // namespace detail

/**
 * @brief Validation policy for hierarchy valuations used as formal saliency levels.
 *
 * A hierarchy valuation assigns one scalar level to each internal tree node. It
 * is compatible with the quasi-flat-zone saliency definition when it is
 * non-decreasing along the ancestry order:
 *
 *     valuation(parent) >= valuation(child)
 *
 * Under this condition, projecting the valuation to adjacency edges through the
 * LCA of the endpoint smallest nodes implements the edge-indexed saliency map of the
 * hierarchy. Equal parent/child values are allowed by `AllowLevelCollapse`; they
 * merge those explicit tree levels in the induced QFZ hierarchy. Use
 * `RequireStrictHierarchy` when the saliency map must preserve every explicit
 * parent-child level of the tree.
 */
enum class HierarchyValuationPolicy {
    AllowLevelCollapse,
    RequireStrictHierarchy,
};

/**
 * @brief Optional value-domain constraint for hierarchy valuations.
 *
 * Cousty et al.'s saliency-map definition is stated for non-negative edge
 * weights. `AllowAnyFinite` lets validation and transformation helpers accept
 * an arbitrary finite ordered valuation so it can be ranked or normalized
 * without changing the induced hierarchy. The formal projection API always
 * applies `RequireNonNegative` because its finest-region base level is `0`.
 */
enum class HierarchyValuationRangePolicy {
    AllowAnyFinite,
    RequireNonNegative,
};

/**
 * @brief Controls spatial-connectivity verification for formal saliency maps.
 *
 * Cousty's correspondence applies to hierarchies whose regions are connected in
 * the graph on which the saliency map is defined. `ValidateConnected` checks
 * this invariant in O(N + (m + p + e) alpha(p)) time and O(N + p + e) storage,
 * where `N` counts internal node slots, `m` live nodes, `p` pixels, `e` graph
 * edges, and `alpha` is the inverse Ackermann function. These bounds exclude
 * the source tree's LCA caches: DFS intervals require O(N) preparation, and
 * the first incomparable pair may add O(N log(N + 1)) time and retained
 * storage for Euler/RMQ. Subsequent queries reuse the caches in O(1) time.
 * `AssumeConnected` is intended for trusted producer-internal hot paths.
 */
enum class HierarchyConnectivityPolicy {
    AssumeConnected,
    ValidateConnected,
};

/**
 * @brief Validates and transforms hierarchy valuations used by saliency maps.
 *
 * A hierarchy valuation assigns one scalar level to each internal tree node. It
 * is compatible with the quasi-flat-zone saliency definition when it is
 * non-decreasing along the ancestry order:
 *
 *     valuation(parent) >= valuation(child)
 *
 * This class centralizes checks and monotone reparameterizations of that
 * valuation. It does not project the hierarchy onto graph edges; use
 * `HierarchySaliencyMap` for the formal edge-indexed saliency map.
 */
class HierarchySaliencyMapValidation {
    friend class HierarchySaliencyMap;

  private:
    /**
     * @brief Validates valuation value.
     *
     * @param value Value.
     * @param nodeId Dense internal node identifier.
     * @param context Operation name used in diagnostics.
     * @param rangePolicy Policy applied to saliency values outside the accepted range.
     */
    template <class Value>
    static void validateValuationValue(const Value& value, NodeId nodeId, const char* context, HierarchyValuationRangePolicy rangePolicy) {
        if constexpr (std::is_floating_point_v<Value>) {
            if (!std::isfinite(value)) {
                std::ostringstream oss;
                oss << context << " requires finite valuation values; node " << nodeId << " has value " << value << ".";
                throw std::invalid_argument(oss.str());
            }
        }
        if (rangePolicy == HierarchyValuationRangePolicy::RequireNonNegative && value < Value{}) {
            std::ostringstream oss;
            oss << context << " requires non-negative valuation values; node " << nodeId << " has value " << value << ".";
            throw std::invalid_argument(oss.str());
        }
    }

  public:
    /**
     * @brief Validates that every hierarchy support is connected in `adjacency`.
     *
     * Each image-domain edge is assigned to the LCA of its endpoint smallest nodes.
     * Nodes are then processed in post-order. At a node, all edges assigned to
     * it are inserted into one disjoint-set forest; the node is connected if its
     * direct proper parts and already-connected child supports have one common
     * representative. This simultaneously validates the completed hierarchy in
     * which every direct proper-part region is a finest graph region.
     *
     * @param tree Hierarchy topology and smallest-node mapping.
     * @param adjacency Graph defining connectedness.
     * @param context Operation name used in diagnostics.
     * @throws std::invalid_argument If the graph domain differs from the tree or
     * any live-node support is disconnected.
     */
    static void validateHierarchyConnectivity(const MorphologicalTree& tree, const RegularGridAdjacency2D& adjacency,
                                              const char* context = "HierarchySaliencyMapValidation::validateHierarchyConnectivity") {
        detail::requireCommittedRootedHierarchy(tree, context);
        const int rows = tree.numRows();
        const int columns = tree.numColumns();
        if (!tree.hasGridDomain2D() || adjacency.getNumRows() != rows || adjacency.getNumColumns() != columns) {
            throw std::invalid_argument(std::string(context) + " requires one graph vertex per 2D proper part and matching adjacency dimensions.");
        }

        validateHierarchyConnectivityEstablished(tree, adjacency, context);
    }

    /**
     * @brief Validates that a node-indexed valuation is compatible with a hierarchy.
     *
     * This is the contract required by the saliency map definition based on
     * quasi-flat zones. The valuation supplies the hierarchy scale; the tree
     * supplies the nested regions. A compatible valuation is defined for every
     * dense internal `NodeId` slot and is monotone along every live parent-child
     * relation. Floating-point valuations must also be finite.
     *
     * `AllowLevelCollapse` accepts `valuation(parent) >= valuation(child)`. The
     * resulting edge map is a formal saliency map for the hierarchy after merging
     * any adjacent levels with equal valuation. `RequireStrictHierarchy` requires
     * `valuation(parent) > valuation(child)` and should be used when the edge map
     * must recover every explicit level of the current tree.
     *
     * @param tree Tree topology.
     * @param valuation Node valuation.
     * @param policy Policy controlling the operation.
     * @param rangePolicy Policy for values outside the supported range.
     * @param context Operation context or diagnostic label.
     */
    template <class Value>
    static void validateHierarchyValuation(const MorphologicalTree& tree, std::span<const Value> valuation,
                                           HierarchyValuationPolicy policy = HierarchyValuationPolicy::AllowLevelCollapse,
                                           HierarchyValuationRangePolicy rangePolicy = HierarchyValuationRangePolicy::AllowAnyFinite,
                                           const char* context = "HierarchySaliencyMapValidation::validateHierarchyValuation") {
        detail::requireCommittedRootedHierarchy(tree, context);
        validateHierarchyValuationEstablished(tree, valuation, policy, rangePolicy, context);
    }

    /**
     * @brief Converts a compatible valuation to dense non-negative integer levels.
     *
     * Cousty et al. state saliency maps on an integer level range. This helper keeps
     * the caller's hierarchy scale order but re-encodes the distinct live-node
     * valuation values as dense ranks `0..k-1`. Equal valuation values receive the
     * same rank, so level collapse is preserved. The input is validated before
     * ranking; floating-point values must be finite and the selected hierarchy
     * monotonicity policy is enforced.
     *
     * @param tree Tree topology.
     * @param valuation Node valuation.
     * @param policy Policy controlling the operation.
     * @return The converted compatible valuation to dense non-negative integer levels.
     */
    template <class Value>
    [[nodiscard]] static std::vector<int> rankHierarchyValuation(const MorphologicalTree& tree, std::span<const Value> valuation,
                                                                 HierarchyValuationPolicy policy = HierarchyValuationPolicy::AllowLevelCollapse) {
        validateHierarchyValuation(tree, valuation, policy, HierarchyValuationRangePolicy::AllowAnyFinite,
                                   "HierarchySaliencyMapValidation::rankHierarchyValuation");

        std::vector<Value> uniqueValues;
        uniqueValues.reserve(static_cast<std::size_t>(tree.numInternalNodeSlots()));
        for (NodeId nodeId : tree.aliveNodeIds()) {
            uniqueValues.push_back(valuation[static_cast<std::size_t>(nodeId)]);
        }

        std::sort(uniqueValues.begin(), uniqueValues.end());
        uniqueValues.erase(std::unique(uniqueValues.begin(), uniqueValues.end()), uniqueValues.end());

        std::vector<int> ranks(static_cast<std::size_t>(tree.numInternalNodeSlots()), 0);
        for (NodeId nodeId : tree.aliveNodeIds()) {
            const Value& value = valuation[static_cast<std::size_t>(nodeId)];
            const auto it = std::lower_bound(uniqueValues.begin(), uniqueValues.end(), value);
            ranks[static_cast<std::size_t>(nodeId)] = static_cast<int>(std::distance(uniqueValues.begin(), it));
        }
        return ranks;
    }

    /**
     * @brief Normalizes a compatible hierarchy valuation to `[0, 1]`.
     *
     * The input valuation is validated first, so it must be finite for
     * floating-point types and monotone along the hierarchy according to `policy`.
     * The returned double buffer preserves the input order with an increasing
     * affine transform over live-node values:
     *
     *     normalized(node) = (valuation(node) - minLiveValue) / range
     *
     * Equal-valued adjacent levels therefore remain collapsed. If every live node
     * has the same valuation, all normalized values are zero. The computation
     * avoids forming an overflowing `maxLiveValue - minLiveValue` when finite
     * floating-point extrema span both signs, and clamps round-off at the output
     * boundaries so every returned live-node score remains in `[0, 1]`.
     *
     * @param tree Tree topology.
     * @param valuation Node valuation.
     * @param policy Policy controlling the operation.
     * @param rangePolicy Policy for values outside the supported range.
     * @return Values produced by the operation.
     */
    template <class Value>
    [[nodiscard]] static std::vector<double>
    computeNormalizedScores(const MorphologicalTree& tree, std::span<const Value> valuation,
                            HierarchyValuationPolicy policy = HierarchyValuationPolicy::AllowLevelCollapse,
                            HierarchyValuationRangePolicy rangePolicy = HierarchyValuationRangePolicy::AllowAnyFinite) {
        using BareValue = std::remove_cv_t<Value>;
        static_assert(std::is_arithmetic_v<BareValue> && !std::is_same_v<BareValue, bool>,
                      "HierarchySaliencyMapValidation::computeNormalizedScores requires a numeric non-bool valuation type.");
        validateHierarchyValuation(tree, valuation, policy, rangePolicy, "HierarchySaliencyMapValidation::computeNormalizedScores");

        return normalizeHierarchyValuation(tree, valuation);
    }

    /**
     * @brief Computes a dense normalized altitude score buffer in `[0, 1]`.
     *
     * The altitude is first converted to a hierarchy-compatible valuation:
     * max-tree altitudes are inverted, while min-tree altitudes keep their natural
     * coarse-to-fine order. The generic `computeNormalizedScores` helper then
     * normalizes that valuation. Trees without a single component-tree polarity are
     * rejected.
     *
     * @param tree Tree topology.
     * @return The computed dense normalized altitude score buffer in [0, 1].
     */
    template <AltitudeValue T> [[nodiscard]] static std::vector<double> computeNormalizedScores(const ValuedMorphologicalTree<T>& tree) {
        const MorphologicalTree& topology = tree.topology();
        detail::requireCommittedRootedHierarchy(topology, "HierarchySaliencyMapValidation::computeNormalizedScores");
        const NodeAltitudeOrder nodeAltitudeOrder = topology.nodeAltitudeOrder();
        if (nodeAltitudeOrder == NodeAltitudeOrder::Unconstrained) {
            throw std::invalid_argument("HierarchySaliencyMapValidation::computeNormalizedScores requires a globally monotone altitude order.");
        }
        // Use at least double precision while retaining long double when it is
        // the input type. An unconditional cast to double can collapse distinct
        // long-double hierarchy levels and change the induced hierarchy.
        using OrientedAltitude = std::common_type_t<T, double>;
        std::vector<OrientedAltitude> orientedAltitude(static_cast<std::size_t>(topology.numInternalNodeSlots()), OrientedAltitude{});
        for (NodeId nodeId : topology.aliveNodeIds()) {
            const OrientedAltitude altitude = static_cast<OrientedAltitude>(tree.nodeAltitudeSpan()[static_cast<std::size_t>(nodeId)]);
            orientedAltitude[static_cast<std::size_t>(nodeId)] = nodeAltitudeOrder == NodeAltitudeOrder::Increasing ? -altitude : altitude;
        }
        return normalizeHierarchyValuation(topology, std::span<const OrientedAltitude>(orientedAltitude));
    }

  private:
    /** @brief Checks spatial support connectivity after the tree and grid domains are established. */
    static void validateHierarchyConnectivityEstablished(const MorphologicalTree& tree, const RegularGridAdjacency2D& adjacency, const char* context) {
        const int numPixels = tree.numPixels();
        struct DomainDisjointSet {
            std::vector<NodeId> parent;
            std::vector<int> size;

            explicit DomainDisjointSet(int count) : parent(static_cast<std::size_t>(count)), size(static_cast<std::size_t>(count), 1) {
                for (NodeId id = 0; id < count; ++id) {
                    parent[static_cast<std::size_t>(id)] = id;
                }
            }

            NodeId find(NodeId id) {
                NodeId root = id;
                while (parent[static_cast<std::size_t>(root)] != root) {
                    root = parent[static_cast<std::size_t>(root)];
                }
                while (parent[static_cast<std::size_t>(id)] != id) {
                    const NodeId next = parent[static_cast<std::size_t>(id)];
                    parent[static_cast<std::size_t>(id)] = root;
                    id = next;
                }
                return root;
            }

            void unite(NodeId lhs, NodeId rhs) {
                lhs = find(lhs);
                rhs = find(rhs);
                if (lhs == rhs) {
                    return;
                }
                if (size[static_cast<std::size_t>(lhs)] < size[static_cast<std::size_t>(rhs)]) {
                    std::swap(lhs, rhs);
                }
                parent[static_cast<std::size_t>(rhs)] = lhs;
                size[static_cast<std::size_t>(lhs)] += size[static_cast<std::size_t>(rhs)];
            }
        };

        using DomainEdge = std::pair<NodeId, NodeId>;
        std::vector<std::vector<DomainEdge>> edgesByLca(static_cast<std::size_t>(tree.numInternalNodeSlots()));
        for (NodeId source = 0; source < numPixels; ++source) {
            const NodeId sourceSmallestNode = detail::TreeKernelAccess::smallestNodeMap(tree, source);
            for (int targetValue : detail::CommittedGridAccess::forwardNeighbors(adjacency, source)) {
                const NodeId target = static_cast<NodeId>(targetValue);
                const NodeId targetSmallestNode = detail::TreeKernelAccess::smallestNodeMap(tree, target);
                const NodeId lca = sourceSmallestNode == targetSmallestNode ? sourceSmallestNode : detail::TreeKernelAccess::lowestCommonAncestor(tree, sourceSmallestNode, targetSmallestNode);
                edgesByLca[static_cast<std::size_t>(lca)].emplace_back(source, target);
            }
        }

        DomainDisjointSet components(numPixels);
        std::vector<NodeId> supportRepresentative(static_cast<std::size_t>(tree.numInternalNodeSlots()), InvalidNode);
        for (NodeId nodeId : tree.postOrder()) {
            for (const DomainEdge& edge : edgesByLca[static_cast<std::size_t>(nodeId)]) {
                components.unite(edge.first, edge.second);
            }

            NodeId representative = InvalidNode;
            auto acceptRepresentative = [&](NodeId candidate) {
                if (representative == InvalidNode) {
                    representative = candidate;
                    return;
                }
                if (components.find(representative) != components.find(candidate)) {
                    std::ostringstream oss;
                    oss << context << " requires every hierarchy region to be connected in the projection graph; node " << nodeId
                        << " has a disconnected support.";
                    throw std::invalid_argument(oss.str());
                }
            };

            for (PixelId pixel : detail::TreeKernelAccess::properParts(tree, nodeId)) {
                acceptRepresentative(pixel);
            }
            for (NodeId childId : detail::TreeKernelAccess::children(tree, nodeId)) {
                const NodeId childRepresentative = supportRepresentative[static_cast<std::size_t>(childId)];
                acceptRepresentative(childRepresentative);
            }
            supportRepresentative[static_cast<std::size_t>(nodeId)] = components.find(representative);
        }
    }

    /** @brief Validates an external valuation over an established committed topology. */
    template <class Value>
    static void validateHierarchyValuationEstablished(const MorphologicalTree& tree, std::span<const Value> valuation,
                                                       HierarchyValuationPolicy policy, HierarchyValuationRangePolicy rangePolicy, const char* context) {
        if (valuation.size() != static_cast<std::size_t>(tree.numInternalNodeSlots())) {
            std::ostringstream oss;
            oss << context << " requires one valuation value per dense internal node slot; expected " << tree.numInternalNodeSlots() << " values but got "
                << valuation.size() << ".";
            throw std::invalid_argument(oss.str());
        }

        for (NodeId nodeId : tree.aliveNodeIds()) {
            validateValuationValue(valuation[static_cast<std::size_t>(nodeId)], nodeId, context, rangePolicy);
        }

        for (NodeId parentId : tree.aliveNodeIds()) {
            const Value& parentValue = valuation[static_cast<std::size_t>(parentId)];
            for (NodeId childId : detail::TreeKernelAccess::children(tree, parentId)) {
                const Value& childValue = valuation[static_cast<std::size_t>(childId)];
                const bool validOrder = policy == HierarchyValuationPolicy::RequireStrictHierarchy ? childValue < parentValue : !(parentValue < childValue);
                if (!validOrder) {
                    std::ostringstream oss;
                    oss << context << " requires "
                        << (policy == HierarchyValuationPolicy::RequireStrictHierarchy ? "valuation(parent) > valuation(child)"
                                                                                       : "valuation(parent) >= valuation(child)")
                        << "; parent node " << parentId << " has value " << parentValue << " and child node " << childId << " has value " << childValue << ".";
                    throw std::invalid_argument(oss.str());
                }
            }
        }
    }

    /** @brief Normalizes an established finite ordered valuation without revalidating it. */
    template <class Value>
    [[nodiscard]] static std::vector<double> normalizeHierarchyValuation(const MorphologicalTree& tree, std::span<const Value> valuation) {
        using BareValue = std::remove_cv_t<Value>;
        bool initialized = false;
        BareValue minValue{};
        BareValue maxValue{};
        for (NodeId nodeId : tree.aliveNodeIds()) {
            const BareValue value = valuation[static_cast<std::size_t>(nodeId)];
            if (!initialized) {
                minValue = value;
                maxValue = value;
                initialized = true;
            } else {
                minValue = std::min(minValue, value);
                maxValue = std::max(maxValue, value);
            }
        }

        std::vector<double> scores(static_cast<std::size_t>(tree.numInternalNodeSlots()), 0.0);
        if (maxValue == minValue) {
            return scores;
        }

        for (NodeId nodeId : tree.aliveNodeIds()) {
            const BareValue value = valuation[static_cast<std::size_t>(nodeId)];
            long double normalized = 0.0L;

            if constexpr (std::is_integral_v<BareValue>) {
                using UnsignedValue = std::make_unsigned_t<BareValue>;
                const UnsignedValue offset = static_cast<UnsignedValue>(value) - static_cast<UnsignedValue>(minValue);
                const UnsignedValue range = static_cast<UnsignedValue>(maxValue) - static_cast<UnsignedValue>(minValue);
                normalized = static_cast<long double>(offset) / static_cast<long double>(range);
            } else {
                const long double low = static_cast<long double>(minValue);
                const long double high = static_cast<long double>(maxValue);
                const long double current = static_cast<long double>(value);

                if (low < 0.0L && high > 0.0L) {
                    const long double scale = std::max(-low, high);
                    const long double scaledLow = low / scale;
                    const long double scaledHigh = high / scale;
                    normalized = (current / scale - scaledLow) / (scaledHigh - scaledLow);
                } else {
                    normalized = (current - low) / (high - low);
                }
            }

            scores[static_cast<std::size_t>(nodeId)] = static_cast<double>(std::clamp(normalized, 0.0L, 1.0L));
        }
        return scores;
    }
};

} // namespace mmcfilters
