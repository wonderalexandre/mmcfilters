#pragma once

#include "HierarchySaliencyMap.hpp"
#include "../MorphologicalTreeFactory.hpp"
#include "../TreeAltitudeAlgorithms.hpp"
#include "../ValuedMorphologicalTreeView.hpp"
#include "../detail/HierarchyCapabilityValidation.hpp"

#include <algorithm>
#include <cmath>
#include <concepts>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mmcfilters {

/**
 * @brief Builds Cousty-style hierarchical-watershed saliency from extinctions.
 *
 * This class implements the constructive route described in Section 8.1 of
 * Cousty et al. A minimum spanning tree of the altitude-ordered component
 * hierarchy is first selected. Extinctions are propagated through its Kruskal
 * binary partition tree; the persistence of a binary merge is the minimum of
 * the maximum descendant extinctions on its two sides. Finally, the persistence
 * valuedTree MST is converted to its QFZ dendrogram and projected back to every
 * edge of the original graph by LCA.
 *
 * The input is a max-tree or min-tree view with one non-negative extinction
 * value per internal node slot. Values are read only at component-tree leaves.
 * Unlike the paper's arbitrary edge-valuedTree input `(G,w)`, this API derives
 * the graph-edge order from the supplied component tree: same-smallest-node edges come
 * first and transition edges are ordered by the altitude of the LCA of their
 * endpoint smallest nodes. This component-tree adapter feeds the Kruskal constructions.
 *
 * @par Primary reference
 * Jean Cousty, Laurent Najman, Yukiko Kenmochi, and Silvio Guimarães,
 * "Hierarchical segmentations with graphs: quasi-flat zones, minimum spanning
 * trees, and saliency maps," Journal of Mathematical Imaging and Vision,
 * 60(4):479-502, 2018.
 * [DOI 10.1007/s10851-017-0768-7](https://doi.org/10.1007/s10851-017-0768-7).
 * The implementation concretizes the BPTAO/extinction/persistence construction
 * sketched in Section 8.1 and makes its otherwise non-unique ties deterministic;
 * the component-tree-derived edge order above is specific to this library.
 */
class HierarchicalWatershedSaliency {
  private:
    /** @brief Small deterministic disjoint-set forest used by both Kruskal passes. */
    class DisjointSet {
        /// Parent slot of each disjoint-set vertex.
        std::vector<NodeId> parent_;
        /// Cardinality of each current disjoint-set root.
        std::vector<int> size_;

      public:
        /**
         * @brief Creates `count` singleton sets.
         *
         * @param count Number of vertices in the disjoint-set domain.
         */
        explicit DisjointSet(int count) : parent_(static_cast<std::size_t>(count)), size_(static_cast<std::size_t>(count), 1) {
            std::iota(parent_.begin(), parent_.end(), NodeId{0});
        }

        /**
         * @brief Finds the representative of one set and compresses its path.
         *
         * @param id Vertex whose representative is requested.
         * @return Current representative of the set containing `id`.
         */
        [[nodiscard]] NodeId find(NodeId id) {
            NodeId root = id;
            while (parent_[static_cast<std::size_t>(root)] != root) {
                root = parent_[static_cast<std::size_t>(root)];
            }
            while (parent_[static_cast<std::size_t>(id)] != id) {
                const NodeId next = parent_[static_cast<std::size_t>(id)];
                parent_[static_cast<std::size_t>(id)] = root;
                id = next;
            }
            return root;
        }

        /**
         * @brief Merges two sets by size with deterministic left-side tie breaking.
         *
         * @param lhs Vertex in the first set.
         * @param rhs Vertex in the second set.
         * @return Representative of the merged set.
         */
        [[nodiscard]] NodeId unite(NodeId lhs, NodeId rhs) {
            lhs = find(lhs);
            rhs = find(rhs);
            if (lhs == rhs) {
                return lhs;
            }
            if (size_[static_cast<std::size_t>(lhs)] < size_[static_cast<std::size_t>(rhs)]) {
                std::swap(lhs, rhs);
            }
            parent_[static_cast<std::size_t>(rhs)] = lhs;
            size_[static_cast<std::size_t>(lhs)] += size_[static_cast<std::size_t>(rhs)];
            return lhs;
        }
    };

    /** @brief One graph edge together with its component-hierarchy merge key. */
    template <AltitudeValue T> struct OrderedGraphEdge {
        /// Row-major source proper-part id.
        NodeId source = InvalidNode;
        /// Row-major target proper-part id.
        NodeId target = InvalidNode;
        /// Lowest common ancestor of the endpoint smallest nodes.
        NodeId lca = InvalidNode;
        /// Altitude of `lca` in the source component hierarchy.
        T altitude{};
        /// Whether both endpoints have the same finest-region smallest node.
        bool finestRegionEdge = false;
        /// Stable forward-adjacency enumeration position.
        std::size_t order = 0;
    };

    /** @brief One selected MST edge and its hierarchical-watershed persistence. */
    template <AltitudeValue T, std::floating_point Real> struct PersistenceEdge {
        /// Original graph edge selected by the first Kruskal pass.
        OrderedGraphEdge<T> graphEdge;
        /// Persistence assigned by the extinction merge rule.
        Real persistence = Real{0};
    };

    /**
     * @brief Compares graph edges for the deterministic altitude-ordered Kruskal pass.
     *
     * Same-smallest-node edges precede transition edges. Remaining ties are resolved by
     * hierarchy altitude and then by row-major endpoint ids.
     *
     * @tparam T Altitude scalar type.
     * @param lhs First edge.
     * @param rhs Second edge.
     * @param nodeAltitudeOrder Monotone orientation of the component-tree altitude.
     * @return `true` when `lhs` must precede `rhs`.
     */
    template <AltitudeValue T> static bool edgePrecedes(const OrderedGraphEdge<T>& lhs, const OrderedGraphEdge<T>& rhs, NodeAltitudeOrder nodeAltitudeOrder) {
        if (lhs.finestRegionEdge != rhs.finestRegionEdge) {
            return lhs.finestRegionEdge;
        }
        if (lhs.altitude != rhs.altitude) {
            if (nodeAltitudeOrder == NodeAltitudeOrder::Increasing) {
                return rhs.altitude < lhs.altitude;
            }
            return lhs.altitude < rhs.altitude;
        }
        if (lhs.source != rhs.source) {
            return lhs.source < rhs.source;
        }
        return lhs.target < rhs.target;
    }

    /**
     * @brief Validates the dense leaf-extinction input.
     *
     * Only live component-tree leaves are interpreted as extinction-bearing
     * extrema; non-leaf slots are ignored by this construction.
     *
     * @tparam Real Floating-point extinction scalar type.
     * @param tree Component-tree topology whose dense node-id domain is used.
     * @param leafExtinction One slot per internal node id.
     * @param context Operation name included in validation errors.
     * @throws std::invalid_argument If the span size is invalid or a leaf value is
     * non-finite or negative.
     */
    template <std::floating_point Real>
    static void validateLeafExtinctions(const MorphologicalTree& tree, std::span<const Real> leafExtinction, const char* context) {
        if (leafExtinction.size() != static_cast<std::size_t>(tree.numInternalNodeSlots())) {
            throw std::invalid_argument(std::string(context) + " requires one extinction slot per internal NodeId.");
        }
        for (NodeId leaf : tree.leaves()) {
            const Real value = leafExtinction[static_cast<std::size_t>(leaf)];
            if (!std::isfinite(value) || value < Real{0}) {
                throw std::invalid_argument(std::string(context) + " requires finite non-negative leaf extinction values.");
            }
        }
    }

    /**
     * @brief Enumerates and deterministically orders all projection-graph edges.
     *
     * @tparam T Altitude scalar type.
     * @param valuedTree Immutable topology/altitude view of the component tree.
     * @param adjacency Connected graph defined on the tree pixel domain.
     * @param context Operation name included in validation errors.
     * @return Forward-enumerated graph edges sorted for the first Kruskal pass.
     * @throws std::runtime_error If an endpoint-smallest-node pair has no live LCA.
     */
    template <AltitudeValue T>
    static std::vector<OrderedGraphEdge<T>> collectOrderedGraphEdges(const ValuedMorphologicalTreeView<T>& valuedTree, const RegularGridAdjacency2D& adjacency,
                                                                     const char* context) {
        const MorphologicalTree& tree = valuedTree.topology();
        std::vector<OrderedGraphEdge<T>> edges;
        std::size_t order = 0;
        for (NodeId source = 0; source < tree.numPixels(); ++source) {
            const NodeId sourceSmallestNode = tree.smallestNode(source);
            for (int targetValue : adjacency.getForwardNeighborIndices(source)) {
                const NodeId target = static_cast<NodeId>(targetValue);
                const NodeId targetSmallestNode = tree.smallestNode(target);
                const bool finestRegionEdge = sourceSmallestNode == targetSmallestNode;
                const NodeId lca = finestRegionEdge ? sourceSmallestNode : tree.lowestCommonAncestor(sourceSmallestNode, targetSmallestNode);
                if (lca == InvalidNode || !tree.isAlive(lca)) {
                    throw std::runtime_error(std::string(context) + " could not find a live LCA for an adjacency edge.");
                }
                edges.push_back(OrderedGraphEdge<T>{source, target, lca, valuedTree.nodeAltitude(lca), finestRegionEdge, order++});
            }
        }
        const NodeAltitudeOrder nodeAltitudeOrder = tree.nodeAltitudeOrder();
        std::stable_sort(edges.begin(), edges.end(), [nodeAltitudeOrder](const auto& lhs, const auto& rhs) { return edgePrecedes(lhs, rhs, nodeAltitudeOrder); });
        return edges;
    }

    /**
     * @brief Selects a deterministic minimum spanning tree by Kruskal's algorithm.
     *
     * @tparam T Altitude scalar type.
     * @param orderedEdges Graph edges in component-hierarchy order.
     * @param numVertices Number of graph vertices/proper parts.
     * @param context Operation name included in validation errors.
     * @return The `numVertices - 1` selected MST edges in Kruskal order.
     * @throws std::invalid_argument If the projection graph is disconnected.
     */
    template <AltitudeValue T>
    static std::vector<OrderedGraphEdge<T>> selectMinimumSpanningTree(std::vector<OrderedGraphEdge<T>> orderedEdges, int numVertices, const char* context) {
        DisjointSet components(numVertices);
        std::vector<OrderedGraphEdge<T>> mst;
        mst.reserve(static_cast<std::size_t>(std::max(0, numVertices - 1)));
        for (const OrderedGraphEdge<T>& edge : orderedEdges) {
            if (components.find(edge.source) == components.find(edge.target)) {
                continue;
            }
            static_cast<void>(components.unite(edge.source, edge.target));
            mst.push_back(edge);
            if (mst.size() == static_cast<std::size_t>(numVertices - 1)) {
                break;
            }
        }
        if (numVertices > 0 && mst.size() != static_cast<std::size_t>(numVertices - 1)) {
            throw std::invalid_argument(std::string(context) + " requires a connected projection graph.");
        }
        return mst;
    }

    /**
     * @brief Assigns Cousty hierarchical-watershed persistence to MST merges.
     *
     * Each merge receives the minimum of the maximum descendant extinctions on
     * its two sides; the merged component retains their maximum.
     *
     * @tparam T Altitude scalar type.
     * @tparam Real Floating-point extinction scalar type.
     * @param tree Component-tree topology that owns the regional-extremum leaves.
     * @param mst MST edges in altitude-ordered Kruskal order.
     * @param leafExtinction Dense extinction slots indexed by internal node id.
     * @return MST edges paired with their persistence values.
     * @throws std::invalid_argument If a component-tree leaf owns no proper part.
     */
    template <AltitudeValue T, std::floating_point Real>
    static std::vector<PersistenceEdge<T, Real>> assignPersistence(const MorphologicalTree& tree, const std::vector<OrderedGraphEdge<T>>& mst,
                                                                   std::span<const Real> leafExtinction) {
        const int numVertices = tree.numPixels();
        std::vector<Real> componentExtinction(static_cast<std::size_t>(numVertices), Real{0});
        for (NodeId leaf : tree.leaves()) {
            const auto properParts = tree.properPart(leaf);
            const auto it = properParts.begin();
            if (it == properParts.end()) {
                throw std::invalid_argument("HierarchicalWatershedSaliency requires every component-tree leaf to own a proper part.");
            }
            componentExtinction[static_cast<std::size_t>(*it)] = leafExtinction[static_cast<std::size_t>(leaf)];
        }

        DisjointSet components(numVertices);
        std::vector<PersistenceEdge<T, Real>> persistenceEdges;
        persistenceEdges.reserve(mst.size());
        for (const OrderedGraphEdge<T>& edge : mst) {
            const NodeId lhsRoot = components.find(edge.source);
            const NodeId rhsRoot = components.find(edge.target);
            const Real lhsExtinction = componentExtinction[static_cast<std::size_t>(lhsRoot)];
            const Real rhsExtinction = componentExtinction[static_cast<std::size_t>(rhsRoot)];
            persistenceEdges.push_back(PersistenceEdge<T, Real>{edge, std::min(lhsExtinction, rhsExtinction)});
            const NodeId mergedRoot = components.unite(lhsRoot, rhsRoot);
            componentExtinction[static_cast<std::size_t>(mergedRoot)] = std::max(lhsExtinction, rhsExtinction);
        }
        return persistenceEdges;
    }

    /**
     * @brief Builds the binary QFZ dendrogram of the persistence-valuedTree MST.
     *
     * Equal persistence values retain the stable order of the first Kruskal pass.
     *
     * @tparam T Source component-tree altitude scalar type.
     * @tparam Real Floating-point persistence scalar type.
     * @param persistenceEdges Persistence-valuedTree MST edges.
     * @param rows Number of rows in the proper-part grid.
     * @param columns Number of columns in the proper-part grid.
     * @param adjacency Projection graph stored on the returned dendrogram.
     * @return Binary native hierarchy with singleton leaves and persistence altitudes.
     */
    template <AltitudeValue T, std::floating_point Real>
    static ValuedMorphologicalTree<Real> buildPersistenceDendrogram(std::vector<PersistenceEdge<T, Real>> persistenceEdges, int rows, int columns,
                                                                      const RegularGridAdjacency2D& adjacency) {
        const int numVertices = rows * columns;
        std::stable_sort(persistenceEdges.begin(), persistenceEdges.end(), [](const auto& lhs, const auto& rhs) {
            if (lhs.persistence != rhs.persistence) {
                return lhs.persistence < rhs.persistence;
            }
            return lhs.graphEdge.order < rhs.graphEdge.order;
        });

        const int numNodes = numVertices == 0 ? 0 : 2 * numVertices - 1;
        std::vector<NodeId> parent(static_cast<std::size_t>(numNodes), InvalidNode);
        std::vector<NodeId> smallestNodeMap(static_cast<std::size_t>(numVertices), InvalidNode);
        std::vector<Real> altitude(static_cast<std::size_t>(numNodes), Real{0});
        for (PixelId pixel = 0; pixel < numVertices; ++pixel) {
            parent[static_cast<std::size_t>(pixel)] = pixel;
            smallestNodeMap[static_cast<std::size_t>(pixel)] = pixel;
        }

        DisjointSet components(numVertices);
        std::vector<NodeId> componentNode(static_cast<std::size_t>(numVertices));
        std::iota(componentNode.begin(), componentNode.end(), NodeId{0});
        NodeId nextNode = numVertices;
        for (const PersistenceEdge<T, Real>& edge : persistenceEdges) {
            const NodeId lhsRoot = components.find(edge.graphEdge.source);
            const NodeId rhsRoot = components.find(edge.graphEdge.target);
            const NodeId lhsNode = componentNode[static_cast<std::size_t>(lhsRoot)];
            const NodeId rhsNode = componentNode[static_cast<std::size_t>(rhsRoot)];
            parent[static_cast<std::size_t>(lhsNode)] = nextNode;
            parent[static_cast<std::size_t>(rhsNode)] = nextNode;
            parent[static_cast<std::size_t>(nextNode)] = nextNode;
            altitude[static_cast<std::size_t>(nextNode)] = edge.persistence;
            const NodeId mergedRoot = components.unite(lhsRoot, rhsRoot);
            componentNode[static_cast<std::size_t>(mergedRoot)] = nextNode;
            ++nextNode;
        }

        const NodeId root = numVertices == 1 ? NodeId{0} : nextNode - 1;
        MorphologicalTreeSemantics semantics{
            MorphologicalTreeKind::Generic, NodeAltitudeOrder::Unconstrained, SharedAdjacencyContext{adjacency}};
        return MorphologicalTreeFactory::createFromNativeTopology(std::span<const NodeId>(parent), std::span<const NodeId>(smallestNodeMap),
                                                                  std::span<const Real>(altitude), root, rows, columns, std::move(semantics));
    }

  public:
    /**
     * @brief Computes the full-graph extinction hierarchical-watershed saliency.
     *
     * @tparam T Component-tree altitude scalar type.
     * @tparam Real Floating-point extinction and output scalar type.
     * @param valuedTree Immutable view of a committed max-tree or min-tree with
     * finite monotone altitudes.
     * @param leafExtinction Dense internal-node-id buffer. Values at live leaves
     * must be finite and non-negative; other slots are ignored.
     * @param adjacency Connected projection graph on exactly the tree proper-part
     * grid domain.
     * @return Full edge-indexed QFZ saliency map with persistence values.
     * @throws std::logic_error If the topology changed after `valuedTree` was made.
     * @throws std::invalid_argument If the hierarchy, adjacency, extinction buffer,
     * leaf smallest-node mapping, or monotone-altitude contract is invalid.
     * @throws std::runtime_error If a live LCA cannot be recovered for a graph edge.
     *
     * @par Complexity
     * Let `m` be the number of tree nodes, `p` the number of proper-part graph
     * vertices, and `e` the number of adjacency edges. The two stable sorts cost
     * `O(e log e + p log p)` time; validation, Kruskal passes, and final LCA
     * projection add linear terms after tree preprocessing. Dominant auxiliary
     * memory is `O(m + e + p)`.
     */
    template <AltitudeValue T, std::floating_point Real>
    [[nodiscard]] static EdgeSaliencyMap<Real> compute(const ValuedMorphologicalTreeView<T>& valuedTree, std::span<const Real> leafExtinction,
                                                       const RegularGridAdjacency2D& adjacency) {
        constexpr const char* context = "HierarchicalWatershedSaliency::compute";
        valuedTree.requireTopologyUnchanged(context);
        const MorphologicalTree& tree = valuedTree.topology();
        detail::validateGlobalMonotoneAltitudeOrder(tree, context);
        TreeAltitudeAlgorithms::validateFiniteAltitudeValues(valuedTree.nodeAltitudes(), context);
        TreeAltitudeAlgorithms::validateMonotoneNodeAltitudes(tree, valuedTree.nodeAltitudes());
        validateLeafExtinctions(tree, leafExtinction, context);
        HierarchySaliencyMapValidation::validateHierarchyConnectivity(tree, adjacency, context);

        const int numVertices = tree.numPixels();
        auto graphEdges = collectOrderedGraphEdges(valuedTree, adjacency, context);
        auto mst = selectMinimumSpanningTree(std::move(graphEdges), numVertices, context);
        auto persistenceEdges = assignPersistence(tree, mst, leafExtinction);
        auto dendrogram = buildPersistenceDendrogram(std::move(persistenceEdges), tree.numRows(), tree.numColumns(), adjacency);
        return HierarchySaliencyMap::computeSaliencyEdgeMap(dendrogram.topology(), adjacency, dendrogram.nodeAltitudeSpan(),
                                                            HierarchyValuationPolicy::AllowLevelCollapse, HierarchyLevelConvention::EdgeSaliencyValue,
                                                            HierarchyConnectivityPolicy::AssumeConnected);
    }

    /**
     * @brief Computes the canonical dense rank scale of `compute`.
     *
     * @tparam T Component-tree altitude scalar type.
     * @tparam Real Floating-point extinction scalar type.
     * @param valuedTree Immutable view of a committed max-tree or min-tree.
     * @param leafExtinction Dense internal-node-id extinction buffer.
     * @param adjacency Connected projection graph on the pixel domain.
     * @return Full edge-indexed saliency map with dense effective-edge ranks.
     * @throws std::logic_error If the topology changed after `valuedTree` was made.
     * @throws std::invalid_argument If an input contract required by `compute` is
     * invalid.
     */
    template <AltitudeValue T, std::floating_point Real>
    [[nodiscard]] static EdgeSaliencyMap<int> computeRanked(const ValuedMorphologicalTreeView<T>& valuedTree, std::span<const Real> leafExtinction,
                                                            const RegularGridAdjacency2D& adjacency) {
        return HierarchySaliencyMap::rankEdgeSaliencyMap(compute(valuedTree, leafExtinction, adjacency));
    }
};

} // namespace mmcfilters
