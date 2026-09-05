#pragma once

#include "ContourTrace.hpp"
#include "../trees/ValuedMorphologicalTreeView.hpp"
#include "../trees/detail/CommittedTreeAccess.hpp"
#include "../trees/detail/MorphologicalTreeConstructionContextQueries.hpp"
#include "../utils/Image.hpp"
#include "detail/ContourBoundaryTracer.hpp"
#include "detail/ContourEdgeDeltaStore.hpp"
#include "detail/ContourTraceTraversal.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <span>
#include <stdexcept>
#include <utility>
#include <variant>
#include <vector>

namespace mmcfilters {

/**
 * @brief Incremental ordered contour traces on the image domain.
 *
 * Construction stores compact edge changes shared by independent traversals.
 * Iteration emits every live node in post-order and retains only the edge sets
 * needed to build an unvisited parent. A query for one node scans that node's
 * support and returns an independently owned trace. The tree must outlive this
 * object and its iterators and remain unchanged.
 */
class ContourTraceComputation {
  private:
    using EdgeDeltas = contours::detail::ContourEdgeDeltaStore;       ///< Compact edge changes by node.
    using ForegroundConnectivity = contours::detail::ForegroundConnectivity; ///< Foreground connectivity at diagonal contacts.
    using BoundaryTracer = contours::detail::ContourBoundaryTracer;   ///< Ordered boundary tracing implementation.
    using VertexIndex = contours::detail::ContourVertexIndex;         ///< Vertex-index representation.

    /** @brief Immutable buffers prepared during construction. */
    struct ConstructionData {
        EdgeDeltas edgeDeltas;                                   ///< Compact edge changes by node.
        std::vector<ForegroundConnectivity> connectivityByNode; ///< Digital connectivity by node.
    };

    /** @brief Immutable indexes shared with independent iterators. */
    struct SharedIndexes {
        const MorphologicalTree& tree;                           ///< Stable source topology.
        std::size_t mutationVersion = 0;                         ///< Captured tree mutation version.
        EdgeDeltas edgeDeltas;                                   ///< Compact edge changes by node.
        std::vector<ForegroundConnectivity> connectivityByNode; ///< Digital connectivity by node.

        /** @brief Takes ownership of prepared indexes for one stable tree. @param source Stable source tree. @param data Prepared indexes. */
        SharedIndexes(const MorphologicalTree& source, ConstructionData data)
            : tree(source), mutationVersion(source.getMutationVersion()), edgeDeltas(std::move(data.edgeDeltas)),
              connectivityByNode(std::move(data.connectivityByNode)) {}

        /** @brief Rejects access after a topology mutation. */
        void requireStableTree() const { tree.requireMutationVersion(mutationVersion, "ContourTraceComputation"); }
    };

    /** @brief Mutable traversal position shared by copies of one input iterator. */
    struct TraversalState {
        std::shared_ptr<const SharedIndexes> indexes;             ///< Shared immutable trace indexes.
        contours::detail::ContourTraceTraversal traversal;       ///< Incremental traversal at the current node.
        bool hasCurrentTrace = false;                             ///< Whether dereference currently yields a trace.

        /** @brief Starts an independent traversal over shared indexes. @param source Shared trace indexes. */
        explicit TraversalState(std::shared_ptr<const SharedIndexes> source)
            : indexes(std::move(source)), traversal(indexes->tree, indexes->edgeDeltas, indexes->connectivityByNode),
              hasCurrentTrace(traversal.advance()) {}
    };

  public:
    /**
     * @brief Single-pass iterator yielding a node and borrowed ordered trace.
     *
     * The trace view remains valid until this iterator or one of its copies
     * advances. Copies share one traversal position. Each call to `begin()`
     * creates an independent traversal.
     */
    class iterator {
      public:
        /// C++20 iterator concept for a single-pass traversal.
        using iterator_concept = std::input_iterator_tag;
        /// Iterator category used by standard algorithms.
        using iterator_category = std::input_iterator_tag;
        /// Node identifier and borrowed ordered trace yielded by dereference.
        using value_type = std::pair<NodeId, ContourTraceView>;
        /// Signed iterator-distance type.
        using difference_type = std::ptrdiff_t;

        iterator() = default;

        /**
         * @brief Borrows the current node and ordered trace.
         * @return Current node identifier and trace view.
         */
        [[nodiscard]] value_type operator*() const {
            if (!state_ || !state_->hasCurrentTrace) {
                throw std::out_of_range("Contour trace iterator is exhausted.");
            }
            return state_->traversal.current();
        }

        /**
         * @brief Advances to the next node trace.
         * @return This iterator after advancing.
         */
        iterator& operator++() {
            if (!state_ || !state_->hasCurrentTrace) {
                throw std::out_of_range("Contour trace iterator is exhausted.");
            }
            state_->hasCurrentTrace = state_->traversal.advance();
            return *this;
        }

        /** @brief Advances without retaining the preceding borrowed trace. */
        void operator++(int) { ++*this; }

        /** @brief Tests exhaustion against the traversal sentinel. */
        friend bool operator==(const iterator& position, std::default_sentinel_t) noexcept {
            return !position.state_ || !position.state_->hasCurrentTrace;
        }

      private:
        friend class ContourTraceComputation;

        /** @brief Creates the first position of an independent traversal. @param indexes Shared trace indexes. */
        explicit iterator(std::shared_ptr<const SharedIndexes> indexes)
            : state_(std::make_shared<TraversalState>(std::move(indexes))) {}

        std::shared_ptr<TraversalState> state_; ///< Shared mutable traversal position.
    };

    /**
     * @brief Packs one pixel-side edge into a compact integer identifier.
     * @param pixel Support pixel identifier.
     * @param side Side occupied by the contour edge.
     * @return Packed contour edge identifier.
     */
    [[nodiscard]] static int packEdge(PixelId pixel, ContourSide side) {
        return contours::detail::packContourEdge(pixel, side);
    }

    /**
     * @brief Unpacks one compact contour edge identifier.
     * @param packedEdge Packed contour edge identifier.
     * @return Support pixel and side represented by the identifier.
     */
    [[nodiscard]] static ContourEdge unpackEdge(int packedEdge) {
        return contours::detail::unpackContourEdge(packedEdge);
    }

    /**
     * @brief Prepares compact edge changes for a stable tree with a 2D domain.
     * @param tree Source topology, which must outlive the computation.
     */
    explicit ContourTraceComputation(const MorphologicalTree& tree)
        : indexes_(std::make_shared<SharedIndexes>(tree, prepareConstructionData(tree))) {}

    /**
     * @brief Prepares traces using the current valued view's shape connectivity.
     * @param view Current valued view whose topology must outlive the computation.
     */
    template <AltitudeValue T>
    explicit ContourTraceComputation(const ValuedMorphologicalTreeView<T>& view)
        : indexes_(std::make_shared<SharedIndexes>(view.topology(), prepareConstructionData(view))) {}

    /**
     * @brief Computes an independently owned ordered trace for one live node.
     *
     * The support is scanned directly. The result can outlive this computation.
     *
     * @param node Live internal node identifier.
     * @return Owned ordered contour trace.
     */
    [[nodiscard]] ContourTrace trace(NodeId node) const {
        indexes_->requireStableTree();
        if (!indexes_->tree.isAlive(node)) {
            throw std::invalid_argument("ContourTraceComputation::trace requires a live internal NodeId.");
        }

        std::vector<int> packedEdges;
        for (PixelId pixel : indexes_->tree.nodeSupport(node)) {
            for (ContourSide side : contourSides()) {
                const PixelId neighbor = adjacentPixel(indexes_->tree, pixel, side);
                if (neighbor == InvalidPixel || !indexes_->tree.isAncestor(node, indexes_->tree.smallestNode(neighbor))) {
                    packedEdges.push_back(packEdge(pixel, side));
                }
            }
        }
        if (packedEdges.empty()) {
            throw std::logic_error("ContourTraceComputation::trace produced an empty edge set for a live node.");
        }

        std::vector<ContourBoundary> boundaries;
        BoundaryTracer tracer(indexes_->tree.numRows(), indexes_->tree.numColumns(), VertexIndex::Sparse);
        tracer.trace(packedEdges, boundaries, indexes_->connectivityByNode[static_cast<std::size_t>(node)]);
        return ContourTrace(std::move(packedEdges), std::move(boundaries));
    }

    /**
     * @brief Starts an independent incremental post-order traversal.
     * @return Iterator positioned at the first trace, or exhausted.
     */
    [[nodiscard]] iterator begin() const {
        indexes_->requireStableTree();
        return iterator(indexes_);
    }

    /** @brief Returns the exhaustion sentinel shared by all traversals. @return Traversal sentinel. */
    [[nodiscard]] std::default_sentinel_t end() const noexcept { return {}; }

    /**
     * @brief Calls `consumer(node, trace)` once for every live node.
     *
     * The borrowed trace view expires after the callback. Exceptions propagate
     * and release traversal storage.
     *
     * @param consumer Callback accepting a node identifier and trace view.
     */
    template <typename Consumer> void forEachTrace(Consumer&& consumer) const {
        for (auto [node, traceView] : *this) {
            consumer(node, traceView);
            indexes_->requireStableTree();
        }
    }

  private:
    /** @brief Returns the four public contour sides in channel order. @return Ordered contour sides. */
    [[nodiscard]] static constexpr std::array<ContourSide, 4> contourSides() {
        return {ContourSide::North, ContourSide::West, ContourSide::East, ContourSide::South};
    }

    /**
     * @brief Converts canonical grid adjacency to digital connectivity.
     * @param adjacency Grid adjacency retained by tree construction.
     * @return Four, eight, or unknown connectivity.
     */
    [[nodiscard]] static ForegroundConnectivity foregroundConnectivity(const RegularGridAdjacency2D& adjacency) {
        if (adjacency.is4connectivity()) {
            return ForegroundConnectivity::Four;
        }
        if (adjacency.is8connectivity()) {
            return ForegroundConnectivity::Eight;
        }
        return ForegroundConnectivity::Unknown;
    }

    /**
     * @brief Returns lower and upper connectivity from construction semantics.
     * @param tree Source tree.
     * @return Lower-shape and upper-shape foreground connectivity.
     */
    [[nodiscard]] static std::array<ForegroundConnectivity, 2> shapeForegroundConnectivities(const MorphologicalTree& tree) {
        if (const auto* adjacency = detail::constructionAdjacency(tree)) {
            const auto connectivity = foregroundConnectivity(*adjacency);
            return {connectivity, connectivity};
        }
        if (const auto* adjacencies = detail::complementaryAdjacencies(tree)) {
            return {foregroundConnectivity(adjacencies->minAdjacency), foregroundConnectivity(adjacencies->maxAdjacency)};
        }
        if (const auto* convention = tree.topographicConvention();
            convention && std::holds_alternative<SelfDualSpanImmersion>(convention->immersion)) {
            return {ForegroundConnectivity::Four, ForegroundConnectivity::Four};
        }
        return {ForegroundConnectivity::Unknown, ForegroundConnectivity::Unknown};
    }

    /**
     * @brief Captures topology-only connectivity choices for every node slot.
     * @param tree Source tree.
     * @return Foreground connectivity indexed by node slot.
     */
    [[nodiscard]] static std::vector<ForegroundConnectivity> foregroundConnectivityByNode(const MorphologicalTree& tree) {
        const auto [lowerShape, upperShape] = shapeForegroundConnectivities(tree);
        return std::vector<ForegroundConnectivity>(static_cast<std::size_t>(tree.numInternalNodeSlots()),
                                                   lowerShape == upperShape ? lowerShape : ForegroundConnectivity::Unknown);
    }

    /**
     * @brief Captures lower and upper shape connectivity from a valued view.
     * @param view Current valued tree view.
     * @return Foreground connectivity indexed by node slot.
     */
    template <AltitudeValue T>
    [[nodiscard]] static std::vector<ForegroundConnectivity> foregroundConnectivityByNode(const ValuedMorphologicalTreeView<T>& view) {
        const MorphologicalTree& tree = view.topology();
        auto connectivityByNode = foregroundConnectivityByNode(tree);
        const auto [lowerShape, upperShape] = shapeForegroundConnectivities(tree);
        if (lowerShape != upperShape) {
            for (NodeId node : tree.aliveNodeIds()) {
                if (tree.isRoot(node)) {
                    continue;
                }
                const auto altitude = view.nodeAltitude(node);
                const auto parentAltitude = view.nodeAltitude(tree.parent(node));
                if (altitude < parentAltitude) {
                    connectivityByNode[static_cast<std::size_t>(node)] = lowerShape;
                } else if (altitude > parentAltitude) {
                    connectivityByNode[static_cast<std::size_t>(node)] = upperShape;
                }
            }
        }
        return connectivityByNode;
    }

    /**
     * @brief Returns the pixel adjacent to one side, or `InvalidPixel`.
     * @param tree Source tree with a grid domain.
     * @param pixel Source pixel.
     * @param side Side whose neighbor is requested.
     * @return Adjacent pixel, or `InvalidPixel` at the image border.
     */
    [[nodiscard]] static PixelId adjacentPixel(const MorphologicalTree& tree, PixelId pixel, ContourSide side) {
        const int rows = tree.numRows();
        const int columns = tree.numColumns();
        const auto [row, column] = ImageUtils::to2D(pixel, columns);

        switch (side) {
        case ContourSide::North:
            return row == 0 ? InvalidPixel : ImageUtils::to1D(row - 1, column, columns);
        case ContourSide::West:
            return column == 0 ? InvalidPixel : ImageUtils::to1D(row, column - 1, columns);
        case ContourSide::East:
            return column == columns - 1 ? InvalidPixel : ImageUtils::to1D(row, column + 1, columns);
        case ContourSide::South:
            return row == rows - 1 ? InvalidPixel : ImageUtils::to1D(row + 1, column, columns);
        }
        return InvalidPixel;
    }

    /**
     * @brief Builds compact edge changes for every live node.
     * @param tree Source tree.
     * @return Compact edge additions and removals indexed by node.
     */
    [[nodiscard]] static EdgeDeltas prepareEdgeDeltas(const MorphologicalTree& tree) {
        if (tree.numRows() <= 0 || tree.numColumns() <= 0) {
            throw std::invalid_argument("Contour tracing requires a non-empty image domain.");
        }
        if (!tree.isAlive(tree.root())) {
            throw std::invalid_argument("Contour tracing requires a live tree root.");
        }

        const int numNodes = tree.numInternalNodeSlots();
        std::vector<EdgeDeltas::Event> additions;
        std::vector<EdgeDeltas::Event> removals;
        additions.reserve(static_cast<std::size_t>(std::max(tree.numPixels(), 1)));
        removals.reserve(static_cast<std::size_t>(std::max(tree.numPixels(), 1)));
        const int rows = tree.numRows();
        const int columns = tree.numColumns();
        const std::span<const NodeId> smallestNodes = tree.smallestNodeMap();

        const auto addBorderEdge = [&](PixelId pixel, ContourSide side) {
            additions.push_back({smallestNodes[static_cast<std::size_t>(pixel)], packEdge(pixel, side)});
        };
        for (int column = 0; column < columns; ++column) {
            addBorderEdge(column, ContourSide::North);
            addBorderEdge((rows - 1) * columns + column, ContourSide::South);
        }
        for (int row = 0; row < rows; ++row) {
            addBorderEdge(row * columns, ContourSide::West);
            addBorderEdge(row * columns + columns - 1, ContourSide::East);
        }

        std::vector<uint8_t> isRightBorder(static_cast<std::size_t>(tree.numPixels()), uint8_t{0});
        for (int row = 0; row < rows; ++row) {
            isRightBorder[static_cast<std::size_t>(row * columns + columns - 1)] = uint8_t{1};
        }
        const std::size_t numAdjacentQueries = 2 * static_cast<std::size_t>(tree.numPixels());
        const auto adjacentPixels = [&](std::size_t queryIndex) {
            const PixelId firstPixel = static_cast<PixelId>(queryIndex >> 1);
            const bool horizontal = (queryIndex & 1) == 0;
            if (horizontal) {
                const PixelId secondPixel = isRightBorder[static_cast<std::size_t>(firstPixel)] ? firstPixel : firstPixel + 1;
                return std::pair{firstPixel, secondPixel};
            }
            const PixelId secondPixel = firstPixel >= tree.numPixels() - columns ? firstPixel : firstPixel + columns;
            return std::pair{firstPixel, secondPixel};
        };
        const auto lcaQuery = [&](std::size_t queryIndex) {
            const auto [firstPixel, secondPixel] = adjacentPixels(queryIndex);
            return std::pair{smallestNodes[static_cast<std::size_t>(firstPixel)],
                             smallestNodes[static_cast<std::size_t>(secondPixel)]};
        };
        detail::CommittedTreeAccess::forEachLowestCommonAncestor(
            tree, numAdjacentQueries, lcaQuery,
            [&](std::size_t queryIndex, NodeId entryNode) {
                const auto [firstPixel, secondPixel] = adjacentPixels(queryIndex);
                if (firstPixel == secondPixel) {
                    return;
                }
                const bool horizontal = (queryIndex & 1) == 0;
                const ContourSide firstSide = horizontal ? ContourSide::East : ContourSide::South;
                const ContourSide secondSide = horizontal ? ContourSide::West : ContourSide::North;
                const NodeId firstNode = smallestNodes[static_cast<std::size_t>(firstPixel)];
                const NodeId secondNode = smallestNodes[static_cast<std::size_t>(secondPixel)];
                if (firstNode != entryNode) {
                    const int packedEdge = packEdge(firstPixel, firstSide);
                    additions.push_back({firstNode, packedEdge});
                    removals.push_back({entryNode, packedEdge});
                }
                if (secondNode != entryNode) {
                    const int packedEdge = packEdge(secondPixel, secondSide);
                    additions.push_back({secondNode, packedEdge});
                    removals.push_back({entryNode, packedEdge});
                }
            });
        return EdgeDeltas::groupDistinct(numNodes, additions, removals);
    }

    /**
     * @brief Validates a topology and prepares shared trace indexes.
     * @param tree Source tree.
     * @return Prepared edge changes and connectivity values.
     */
    [[nodiscard]] static ConstructionData prepareConstructionData(const MorphologicalTree& tree) {
        tree.requireNotEditing("ContourTraceComputation");
        return {prepareEdgeDeltas(tree), foregroundConnectivityByNode(tree)};
    }

    /**
     * @brief Validates a valued view and prepares shared trace indexes.
     * @param view Current valued tree view.
     * @return Prepared edge changes and connectivity values.
     */
    template <AltitudeValue T>
    [[nodiscard]] static ConstructionData prepareConstructionData(const ValuedMorphologicalTreeView<T>& view) {
        view.requireTopologyUnchanged("ContourTraceComputation");
        view.topology().requireNotEditing("ContourTraceComputation");
        return {prepareEdgeDeltas(view.topology()), foregroundConnectivityByNode(view)};
    }

    std::shared_ptr<const SharedIndexes> indexes_; ///< Immutable data shared with active iterators.
};

} // namespace mmcfilters
