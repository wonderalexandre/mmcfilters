#pragma once

#include "HierarchySaliencyMap.hpp"
#include "../../utils/Image.hpp"
#include "../detail/MorphologicalTreeConstructionContextQueries.hpp"

#include <cstddef>
#include <cstdint>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace mmcfilters {

/**
 * @brief Edge-indexed contour set induced by thresholding a saliency map.
 *
 * The arrays are parallel: `sources[i]` and `targets[i]` describe one
 * undirected image-domain graph edge selected as a contour edge.
 */
struct EdgeContourMap {
    /// Number of rows in the image grid.
    int numRows = 0;
    /// Number of columns in the image grid.
    int numColumns = 0;
    /// Radius of the adjacency used to enumerate the edges.
    double adjacencyRadius = 0.0;
    /// Source pixel identifier of each selected edge.
    std::vector<NodeId> sources;
    /// Target pixel identifier of each selected edge.
    std::vector<NodeId> targets;

    /**
     * @brief Returns the number of selected edges.
     *
     * @return The number of selected edges.
     */
    [[nodiscard]] std::size_t size() const noexcept { return sources.size(); }

    /**
     * @brief Returns whether no contour edge is selected.
     *
     * @return Whether no contour edge is selected.
     */
    [[nodiscard]] bool empty() const noexcept { return sources.empty(); }
};

/**
 * @brief Pixel aggregation policy for rasterizing edge maps for display.
 */
enum class EdgeToPixelReducer {
    /// Use the maximum value among incident edge values for each pixel.
    Max,
    /// Use the arithmetic mean of incident edge values for each pixel.
    Mean
};

/**
 * @brief Edge-indexed contour set projected onto hierarchy nodes.
 *
 * `nodes[i]` is the hierarchy node whose boundary contains the image-domain
 * transition edge `(sources[i], targets[i])`. It is computed as
 * `LCA(smallestNode(sources[i]), smallestNode(targets[i]))`.
 */
struct NodeContourEdgeMap {
    /// Number of rows in the image grid.
    int numRows = 0;
    /// Number of columns in the image grid.
    int numColumns = 0;
    /// Radius of the adjacency used to enumerate the edges.
    double adjacencyRadius = 0.0;
    /// Source pixel identifier of each contour edge.
    std::vector<NodeId> sources;
    /// Target pixel identifier of each contour edge.
    std::vector<NodeId> targets;
    /// Boundary-owning hierarchy node parallel to each edge.
    std::vector<NodeId> nodes;

    /**
     * @brief Returns the number of projected contour edges.
     *
     * @return The number of projected contour edges.
     */
    [[nodiscard]] std::size_t size() const noexcept { return nodes.size(); }

    /**
     * @brief Returns whether the projection has no contour edges.
     *
     * @return Whether the projection has no contour edges.
     */
    [[nodiscard]] bool empty() const noexcept { return nodes.empty(); }
};

/**
 * @brief Per-node incremental contour edges induced by a hierarchy.
 *
 * For a node `u`, the slice
 * `sources[offsets[u]:offsets[u + 1]]` and
 * `targets[offsets[u]:offsets[u + 1]]` stores exactly the transition adjacency
 * edges whose endpoint smallest nodes are distinct and whose endpoint-smallest-node LCA is `u`.
 * Edges whose endpoints already belong to the same finest region are omitted;
 * in the full formal saliency map they have the implicit base value `0`.
 */
struct IncrementalNodeContourMap {
    /// Number of rows in the image grid.
    int numRows = 0;
    /// Number of columns in the image grid.
    int numColumns = 0;
    /// Number of slots in the dense node-id domain.
    int numNodeSlots = 0;
    /// Radius of the adjacency used to enumerate the edges.
    double adjacencyRadius = 0.0;
    /// CSR-style boundaries of the edge slice owned by each node.
    std::vector<std::size_t> offsets;
    /// Source pixel ids grouped by boundary-owning node.
    std::vector<NodeId> sources;
    /// Target pixel ids parallel to `sources`.
    std::vector<NodeId> targets;

    /**
     * @brief Returns the total number of stored contour edges.
     *
     * @return The total number of stored contour edges.
     */
    [[nodiscard]] std::size_t size() const noexcept { return sources.size(); }

    /**
     * @brief Returns whether no contour edge is stored.
     *
     * @return Whether no contour edge is stored.
     */
    [[nodiscard]] bool empty() const noexcept { return sources.empty(); }

    /**
     * @brief Returns the inclusive edge offset for `nodeId`.
     *
     * @param nodeId Dense internal node identifier.
     * @return The inclusive edge offset for nodeId.
     */
    [[nodiscard]] std::size_t nodeBegin(NodeId nodeId) const { return offsets.at(static_cast<std::size_t>(nodeId)); }

    /**
     * @brief Returns the exclusive edge offset for `nodeId`.
     *
     * @param nodeId Dense internal node identifier.
     * @return The exclusive edge offset for nodeId.
     */
    [[nodiscard]] std::size_t nodeEnd(NodeId nodeId) const { return offsets.at(static_cast<std::size_t>(nodeId) + 1); }

    /**
     * @brief Returns the number of contour edges owned by `nodeId`.
     *
     * @param nodeId Dense internal node identifier.
     * @return The number of contour edges owned by nodeId.
     */
    [[nodiscard]] std::size_t nodeSize(NodeId nodeId) const { return nodeEnd(nodeId) - nodeBegin(nodeId); }
};

/**
 * @brief Derived projections and contour materializations of saliency maps.
 *
 * This class does not construct the formal edge-indexed saliency map. It
 * consumes edge maps, trees, or precomputed transition-contour slices to produce
 * display images, threshold cuts, and sparse transition-edge projections.
 */
class HierarchySaliencyMapProjection {
  private:
    /**
     * @brief Validates edge map.
     *
     * @param edgeMap Saliency values indexed by graph edge.
     * @param context Operation name used in diagnostics.
     */
    template <class Value> static void validateEdgeMap(const EdgeSaliencyMap<Value>& edgeMap, const char* context) {
        if (edgeMap.sources.size() != edgeMap.targets.size() || edgeMap.sources.size() != edgeMap.values.size()) {
            throw std::invalid_argument(std::string(context) + " requires sources, targets, and values to have the same length.");
        }
        if (edgeMap.numRows <= 0 || edgeMap.numColumns <= 0) {
            throw std::invalid_argument(std::string(context) + " requires a non-empty image domain.");
        }
    }

    /**
     * @brief Validates incremental contours.
     *
     * @param contours Contour representation used to project saliency values.
     * @param context Operation name used in diagnostics.
     */
    static void validateIncrementalContours(const IncrementalNodeContourMap& contours, const char* context) {
        if (contours.numRows <= 0 || contours.numColumns <= 0) {
            throw std::invalid_argument(std::string(context) + " requires a non-empty image domain.");
        }
        if (contours.numNodeSlots < 0) {
            throw std::invalid_argument(std::string(context) + " requires a non-negative node-slot count.");
        }
        if (contours.sources.size() != contours.targets.size()) {
            throw std::invalid_argument(std::string(context) + " requires sources and targets to have the same length.");
        }
        if (contours.offsets.size() != static_cast<std::size_t>(contours.numNodeSlots) + 1) {
            throw std::invalid_argument(std::string(context) + " requires offsets to have numNodeSlots + 1 entries.");
        }
        if (contours.offsets.empty() || contours.offsets.front() != 0 || contours.offsets.back() != contours.sources.size()) {
            throw std::invalid_argument(std::string(context) + " requires offsets to span exactly the contour edge arrays.");
        }
        for (std::size_t i = 1; i < contours.offsets.size(); ++i) {
            if (contours.offsets[i] < contours.offsets[i - 1] || contours.offsets[i] > contours.sources.size()) {
                throw std::invalid_argument(std::string(context) + " requires monotonically increasing offsets inside the edge-array bounds.");
            }
        }
    }

    /**
     * @brief Validates edge endpoints.
     *
     * @param edgeMap Saliency values indexed by graph edge.
     * @param context Operation name used in diagnostics.
     */
    template <class Value> static void validateEdgeEndpoints(const EdgeSaliencyMap<Value>& edgeMap, const char* context) {
        const auto pixelCount = static_cast<std::size_t>(edgeMap.numRows) * static_cast<std::size_t>(edgeMap.numColumns);
        for (std::size_t i = 0; i < edgeMap.size(); ++i) {
            const NodeId source = edgeMap.sources[i];
            const NodeId target = edgeMap.targets[i];
            if (source < 0 || static_cast<std::size_t>(source) >= pixelCount || target < 0 || static_cast<std::size_t>(target) >= pixelCount) {
                throw std::invalid_argument(std::string(context) + " requires sources and targets inside the image domain.");
            }
        }
    }

    /**
     * @brief Validates tree and adjacency.
     *
     * @param tree Tree topology.
     * @param adjacency Adjacency relation.
     * @param context Operation name used in diagnostics.
     */
    static void validateTreeAndAdjacency(const MorphologicalTree& tree, const RegularGridAdjacency2D& adjacency, const char* context) {
        if (tree.root() == InvalidNode) {
            throw std::invalid_argument(std::string(context) + " requires a non-empty rooted tree.");
        }
        if (tree.numRows() <= 0 || tree.numColumns() <= 0 || tree.numPixels() <= 0) {
            throw std::invalid_argument(std::string(context) + " requires a non-empty image/pixel domain.");
        }
        if (adjacency.getNumRows() != tree.numRows() || adjacency.getNumColumns() != tree.numColumns()) {
            throw std::invalid_argument(std::string(context) + " adjacency domain must match the tree image domain.");
        }
    }

    /**
     * @brief Groups node contour edges.
     *
     * @param flat Flattened saliency values to reshape into the image domain.
     * @param numNodeSlots Count.
     * @return Contour edges grouped by tree node.
     */
    static IncrementalNodeContourMap groupNodeContourEdges(NodeContourEdgeMap&& flat, int numNodeSlots) {
        IncrementalNodeContourMap contours;
        contours.numRows = flat.numRows;
        contours.numColumns = flat.numColumns;
        contours.numNodeSlots = numNodeSlots;
        contours.adjacencyRadius = flat.adjacencyRadius;
        contours.offsets.assign(static_cast<std::size_t>(contours.numNodeSlots) + 1, 0);

        for (NodeId nodeId : flat.nodes) {
            if (nodeId < 0 || nodeId >= contours.numNodeSlots) {
                throw std::runtime_error("HierarchySaliencyMapProjection::computeIncrementalNodeContours found an LCA outside the dense node-id domain.");
            }
            ++contours.offsets[static_cast<std::size_t>(nodeId) + 1];
        }
        std::partial_sum(contours.offsets.begin(), contours.offsets.end(), contours.offsets.begin());

        contours.sources.resize(flat.size());
        contours.targets.resize(flat.size());
        std::vector<std::size_t> cursor = contours.offsets;
        for (std::size_t i = 0; i < flat.size(); ++i) {
            const NodeId nodeId = flat.nodes[i];
            const std::size_t dst = cursor[static_cast<std::size_t>(nodeId)]++;
            contours.sources[dst] = flat.sources[i];
            contours.targets[dst] = flat.targets[i];
        }

        return contours;
    }

  public:
    /**
     * @brief Rasterizes an edge-indexed map into a pixel image for display.
     *
     * This helper is only a visualization projection. It does not change the
     * formal edge-indexed saliency representation. Each edge value contributes to
     * both endpoint pixels. `Max` writes the maximum incident edge value per pixel,
     * while `Mean` writes the arithmetic mean of incident edge values. Pixels with
     * no incident edges keep value `0`.
     *
     * @param edgeMap Edge-indexed saliency or score map.
     * @param reducer Pixel aggregation policy.
     * @return Row-major `double` image on the same image domain as `edgeMap`.
     */
    template <class Value>
    [[nodiscard]] static ImagePtr<double> edgeMapToPixelImage(const EdgeSaliencyMap<Value>& edgeMap, EdgeToPixelReducer reducer = EdgeToPixelReducer::Max) {
        constexpr const char* context = "HierarchySaliencyMapProjection::edgeMapToPixelImage";
        validateEdgeMap(edgeMap, context);
        validateEdgeEndpoints(edgeMap, context);

        ImagePtr<double> image = Image<double>::create(edgeMap.numRows, edgeMap.numColumns, 0.0);
        double* pixels = image->rawData();
        const auto pixelCount = static_cast<std::size_t>(edgeMap.numRows) * static_cast<std::size_t>(edgeMap.numColumns);

        switch (reducer) {
        case EdgeToPixelReducer::Max: {
            std::vector<uint8_t> initialized(pixelCount, false);
            for (std::size_t i = 0; i < edgeMap.size(); ++i) {
                const double value = static_cast<double>(edgeMap.values[i]);
                const auto source = static_cast<std::size_t>(edgeMap.sources[i]);
                const auto target = static_cast<std::size_t>(edgeMap.targets[i]);
                if (!initialized[source] || pixels[source] < value) {
                    pixels[source] = value;
                    initialized[source] = true;
                }
                if (!initialized[target] || pixels[target] < value) {
                    pixels[target] = value;
                    initialized[target] = true;
                }
            }
            return image;
        }
        case EdgeToPixelReducer::Mean: {
            std::vector<std::size_t> count(pixelCount, 0);
            for (std::size_t i = 0; i < edgeMap.size(); ++i) {
                const double value = static_cast<double>(edgeMap.values[i]);
                const auto source = static_cast<std::size_t>(edgeMap.sources[i]);
                const auto target = static_cast<std::size_t>(edgeMap.targets[i]);
                pixels[source] += value;
                pixels[target] += value;
                ++count[source];
                ++count[target];
            }
            for (std::size_t i = 0; i < pixelCount; ++i) {
                if (count[i] != 0) {
                    pixels[i] /= static_cast<double>(count[i]);
                }
            }
            return image;
        }
        default:
            throw std::invalid_argument(std::string(context) + " received an unknown edge-to-pixel reducer.");
        }
    }

    /**
     * @brief Thresholds an edge saliency map into an edge contour set.
     *
     * A contour edge is selected when `saliency(edge) >= threshold`. This is the
     * cut complement of the quasi-flat-zone convention where graph edges with
     * saliency strictly below the level remain connected.
     *
     * @param edgeMap Edge saliency map to threshold.
     * @param threshold Threshold applied by the operation.
     * @return The thresholded edge saliency map into an edge contour set.
     */
    template <class Value, class Threshold> [[nodiscard]] static EdgeContourMap thresholdCut(const EdgeSaliencyMap<Value>& edgeMap, Threshold threshold) {
        constexpr const char* context = "HierarchySaliencyMapProjection::thresholdCut";
        validateEdgeMap(edgeMap, context);

        EdgeContourMap contours;
        contours.numRows = edgeMap.numRows;
        contours.numColumns = edgeMap.numColumns;
        contours.adjacencyRadius = edgeMap.adjacencyRadius;
        contours.sources.reserve(edgeMap.size());
        contours.targets.reserve(edgeMap.size());

        for (std::size_t i = 0; i < edgeMap.size(); ++i) {
            if (edgeMap.values[i] >= threshold) {
                contours.sources.push_back(edgeMap.sources[i]);
                contours.targets.push_back(edgeMap.targets[i]);
            }
        }
        return contours;
    }

    /**
     * @brief Projects transition adjacency edges onto their hierarchy node.
     *
     * Edges whose endpoints already have the same smallest hierarchy node are
     * omitted. They are internal to the finest region represented by that smallest node
     * and carry the implicit base value `0` in a full formal saliency map.
     *
     * @param tree Tree topology.
     * @param adjacency Adjacency relation.
     * @return The projected transition adjacency edges onto their hierarchy node.
     */
    [[nodiscard]] static NodeContourEdgeMap nodeContourEdges(const MorphologicalTree& tree, const RegularGridAdjacency2D& adjacency) {
        constexpr const char* context = "HierarchySaliencyMapProjection::nodeContourEdges";
        validateTreeAndAdjacency(tree, adjacency, context);

        NodeContourEdgeMap contours;
        contours.numRows = tree.numRows();
        contours.numColumns = tree.numColumns();
        contours.adjacencyRadius = adjacency.getRadius();

        const int numPixels = tree.numPixels();
        for (NodeId source = 0; source < numPixels; ++source) {
            const NodeId sourceSmallestNode = tree.smallestNode(source);
            if (!tree.isAlive(sourceSmallestNode)) {
                throw std::runtime_error(std::string(context) + " found a pixel without a live smallest node.");
            }

            for (int target : adjacency.getForwardNeighborIndices(source)) {
                const NodeId targetSmallestNode = tree.smallestNode(target);
                if (!tree.isAlive(targetSmallestNode)) {
                    throw std::runtime_error(std::string(context) + " found a neighbour pixel without a live smallest node.");
                }
                if (sourceSmallestNode == targetSmallestNode) {
                    continue;
                }

                const NodeId lca = tree.lowestCommonAncestor(sourceSmallestNode, targetSmallestNode);
                if (lca == InvalidNode || !tree.isAlive(lca)) {
                    throw std::runtime_error(std::string(context) + " could not find a live LCA for an adjacency edge.");
                }

                contours.sources.push_back(source);
                contours.targets.push_back(target);
                contours.nodes.push_back(lca);
            }
        }
        return contours;
    }

    /**
     * @brief Projects stored-adjacency transition edges onto hierarchy nodes.
     *
     * @param tree Tree topology.
     * @return The projected stored-adjacency transition edges onto hierarchy nodes.
     */
    [[nodiscard]] static NodeContourEdgeMap nodeContourEdges(const MorphologicalTree& tree) {
        const RegularGridAdjacency2D* adjacency = ::mmcfilters::detail::constructionAdjacency(tree);
        if (adjacency == nullptr) {
            throw std::invalid_argument(
                "HierarchySaliencyMapProjection::nodeContourEdges requires an attached adjacency relation; pass an explicit adjacency relation instead.");
        }
        return nodeContourEdges(tree, *adjacency);
    }

    /**
     * @brief Projects valued-tree adjacency edges onto their hierarchy node.
     *
     * @param tree Tree topology.
     * @param adjacency Adjacency relation.
     * @return The projected valued-tree adjacency edges onto their hierarchy node.
     */
    template <AltitudeValue T>
    [[nodiscard]] static NodeContourEdgeMap nodeContourEdges(const ValuedMorphologicalTree<T>& tree, const RegularGridAdjacency2D& adjacency) {
        return nodeContourEdges(tree.topology(), adjacency);
    }

    /**
     * @brief Projects valued-tree stored-adjacency edges onto hierarchy nodes.
     *
     * @param tree Tree topology.
     * @return The projected valued-tree stored-adjacency edges onto hierarchy nodes.
     */
    template <AltitudeValue T> [[nodiscard]] static NodeContourEdgeMap nodeContourEdges(const ValuedMorphologicalTree<T>& tree) {
        return nodeContourEdges(tree.topology());
    }

    /**
     * @brief Computes per-node incremental contour edges with explicit adjacency.
     *
     * The output stores one contiguous edge slice per dense internal `NodeId`.
     * A transition edge belongs to node `u` iff its endpoint smallest nodes are distinct
     * and `u = LCA(smallestNode(source), smallestNode(target))`.
     *
     * @param tree Tree topology.
     * @param adjacency Adjacency relation.
     * @return The computed per-node incremental contour edges with explicit adjacency.
     */
    [[nodiscard]] static IncrementalNodeContourMap computeIncrementalNodeContours(const MorphologicalTree& tree, const RegularGridAdjacency2D& adjacency) {
        return groupNodeContourEdges(nodeContourEdges(tree, adjacency), tree.numInternalNodeSlots());
    }

    /**
     * @brief Computes per-node incremental contour edges using stored adjacency.
     *
     * @param tree Tree topology.
     * @return The computed per-node incremental contour edges using stored adjacency.
     */
    [[nodiscard]] static IncrementalNodeContourMap computeIncrementalNodeContours(const MorphologicalTree& tree) {
        return groupNodeContourEdges(nodeContourEdges(tree), tree.numInternalNodeSlots());
    }

    /**
     * @brief Computes valued-tree incremental contours with explicit adjacency.
     *
     * @param tree Tree topology.
     * @param adjacency Adjacency relation.
     * @return The computed valued-tree incremental contours with explicit adjacency.
     */
    template <AltitudeValue T>
    [[nodiscard]] static IncrementalNodeContourMap computeIncrementalNodeContours(const ValuedMorphologicalTree<T>& tree,
                                                                                  const RegularGridAdjacency2D& adjacency) {
        return computeIncrementalNodeContours(tree.topology(), adjacency);
    }

    /**
     * @brief Computes valued-tree incremental contours using stored adjacency.
     *
     * @param tree Tree topology.
     * @return The computed valued-tree incremental contours using stored adjacency.
     */
    template <AltitudeValue T> [[nodiscard]] static IncrementalNodeContourMap computeIncrementalNodeContours(const ValuedMorphologicalTree<T>& tree) {
        return computeIncrementalNodeContours(tree.topology());
    }

    /**
     * @brief Projects dense node valuation onto transition contour edges.
     *
     * This is a sparse projection over the transition-edge support returned by
     * `computeIncrementalNodeContours`. Edges internal to the finest represented
     * regions are omitted and have implicit value `0` in the full formal saliency
     * map. Use `HierarchySaliencyMap::computeSaliencyEdgeMap` when the output must
     * explicitly contain every graph edge.
     *
     * @param contours Contour data.
     * @param nodeValuation Dense valuation indexed by node identifier.
     * @return The projected dense node valuation onto transition contour edges.
     */
    template <class Value>
    [[nodiscard]] static EdgeSaliencyMap<Value> projectNodeValuation(const IncrementalNodeContourMap& contours, std::span<const Value> nodeValuation) {
        constexpr const char* context = "HierarchySaliencyMapProjection::projectNodeValuation";
        validateIncrementalContours(contours, context);
        if (nodeValuation.size() != static_cast<std::size_t>(contours.numNodeSlots)) {
            throw std::invalid_argument(std::string(context) + " requires one valuation value per dense node slot.");
        }

        EdgeSaliencyMap<Value> edgeMap;
        edgeMap.numRows = contours.numRows;
        edgeMap.numColumns = contours.numColumns;
        edgeMap.adjacencyRadius = contours.adjacencyRadius;
        edgeMap.sources.reserve(contours.size());
        edgeMap.targets.reserve(contours.size());
        edgeMap.values.reserve(contours.size());

        for (NodeId nodeId = 0; nodeId < contours.numNodeSlots; ++nodeId) {
            const Value value = nodeValuation[static_cast<std::size_t>(nodeId)];
            for (std::size_t i = contours.offsets[static_cast<std::size_t>(nodeId)]; i < contours.offsets[static_cast<std::size_t>(nodeId) + 1]; ++i) {
                edgeMap.sources.push_back(contours.sources[i]);
                edgeMap.targets.push_back(contours.targets[i]);
                edgeMap.values.push_back(value);
            }
        }
        return edgeMap;
    }

    /**
     * @brief Thresholds incremental contour edges by dense node valuation.
     *
     * This is a generic fast cut helper over precomputed transition-contour
     * slices. It does not see the tree topology and therefore cannot validate that
     * `nodeValuation` is a formal hierarchy valuation. Validate the node valuation
     * separately when this cut is intended to represent a quasi-flat-zone saliency
     * level.
     *
     * @param contours Contour data.
     * @param nodeValuation Dense valuation indexed by node identifier.
     * @param threshold Threshold applied by the operation.
     * @return The thresholded incremental contour edges by dense node valuation.
     */
    template <class Value, class Threshold>
    [[nodiscard]] static EdgeContourMap thresholdByNodeValuation(const IncrementalNodeContourMap& contours, std::span<const Value> nodeValuation,
                                                                 Threshold threshold) {
        constexpr const char* context = "HierarchySaliencyMapProjection::thresholdByNodeValuation";
        validateIncrementalContours(contours, context);
        if (nodeValuation.size() != static_cast<std::size_t>(contours.numNodeSlots)) {
            throw std::invalid_argument(std::string(context) + " requires one valuation value per dense node slot.");
        }

        EdgeContourMap cut;
        cut.numRows = contours.numRows;
        cut.numColumns = contours.numColumns;
        cut.adjacencyRadius = contours.adjacencyRadius;
        cut.sources.reserve(contours.size());
        cut.targets.reserve(contours.size());

        for (NodeId nodeId = 0; nodeId < contours.numNodeSlots; ++nodeId) {
            if (nodeValuation[static_cast<std::size_t>(nodeId)] >= threshold) {
                for (std::size_t i = contours.offsets[static_cast<std::size_t>(nodeId)]; i < contours.offsets[static_cast<std::size_t>(nodeId) + 1]; ++i) {
                    cut.sources.push_back(contours.sources[i]);
                    cut.targets.push_back(contours.targets[i]);
                }
            }
        }
        return cut;
    }
};

} // namespace mmcfilters
