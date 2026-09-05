#pragma once

#include "../ContourTrace.hpp"
#include "../../utils/Image.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace mmcfilters::contours::detail {

/** @brief Digital foreground connectivity used at diagonal boundary contacts. */
enum class ForegroundConnectivity : uint8_t { Unknown, Four, Eight };

/** @brief Vertex-index representation used while ordering one edge set. */
enum class ContourVertexIndex : uint8_t { Dense, Sparse };

/**
 * @brief Orders packed contour edges into closed external and internal boundaries.
 *
 * Scratch buffers are retained across calls. The supplied packed-edge vector is
 * replaced by its boundary-ordered permutation without copying the result.
 */
class ContourBoundaryTracer {
  private:
    /** @brief Cardinal direction of an oriented contour edge. */
    enum class Direction : uint8_t { North = 0, East = 1, South = 2, West = 3 };

    /** @brief Geometry retained while ordering one contour edge. */
    struct TraceEdgeData {
        int endVertex = -1;        ///< Grid vertex where the edge ends.
        int doubledSignedArea = 0; ///< Contribution to twice the signed area.
    };

    /** @brief Grid geometry derived from one packed contour edge. */
    struct EdgeGeometry {
        int startVertex = -1;      ///< Grid vertex where the edge begins.
        int endVertex = -1;        ///< Grid vertex where the edge ends.
        int doubledSignedArea = 0; ///< Contribution to twice the signed area.
    };

  public:
    /**
     * @brief Allocates the reusable direct vertex index for an image domain.
     * @param rows Number of image rows.
     * @param columns Number of image columns.
     * @param vertexIndex Vertex-index representation used by the trace.
     */
    ContourBoundaryTracer(int rows, int columns, ContourVertexIndex vertexIndex = ContourVertexIndex::Dense)
        : columns_(columns), vertexIndex_(vertexIndex) {
        if (vertexIndex_ == ContourVertexIndex::Dense) {
            outgoingEdgeHeads_.assign(static_cast<std::size_t>((rows + 1) * (columns + 1)), -1);
        }
    }

    /**
     * @brief Orders one node's contour edges and records its closed boundaries.
     * @param packedEdges Distinct unordered input edges, replaced by ordered edges.
     * @param boundaries Output descriptors with offsets relative to `packedEdges`.
     * @param connectivity Foreground connectivity used at diagonal contacts.
     */
    void trace(std::vector<int>& packedEdges, std::vector<ContourBoundary>& boundaries, ForegroundConnectivity connectivity) {
        resetVertexIndex();
        edgeData_.clear();
        edgeData_.reserve(packedEdges.size());
        nextOutgoingEdges_.clear();
        nextOutgoingEdges_.reserve(packedEdges.size());
        if (vertexIndex_ == ContourVertexIndex::Sparse) {
            sparseOutgoingEdgeHeads_.reserve(packedEdges.size());
        }

        for (int packedEdge : packedEdges) {
            const EdgeGeometry geometry = edgeGeometry(packedEdge);
            const int edgeIndex = static_cast<int>(edgeData_.size());
            edgeData_.push_back(TraceEdgeData{geometry.endVertex, geometry.doubledSignedArea});
            addOutgoingEdge(geometry.startVertex, edgeIndex);
        }

        visitedEdges_.assign(edgeData_.size(), 0);
        orderedEdges_.clear();
        orderedEdges_.reserve(edgeData_.size());
        boundaries.clear();
        boundaries.reserve((edgeData_.size() / 4) + 1);

        try {
            for (int startEdgeIndex = 0; startEdgeIndex < static_cast<int>(edgeData_.size()); ++startEdgeIndex) {
                if (isVisited(startEdgeIndex)) {
                    continue;
                }

                const std::size_t firstBoundaryEdge = orderedEdges_.size();
                int doubledSignedArea = 0;
                int currentEdgeIndex = startEdgeIndex;

                do {
                    if (isVisited(currentEdgeIndex)) {
                        throw std::runtime_error("Contour successor revisited an edge before closing its cycle.");
                    }
                    const TraceEdgeData& edge = edgeData_[static_cast<std::size_t>(currentEdgeIndex)];
                    visitedEdges_[static_cast<std::size_t>(currentEdgeIndex)] = 1;
                    orderedEdges_.push_back(packedEdges[static_cast<std::size_t>(currentEdgeIndex)]);
                    doubledSignedArea += edge.doubledSignedArea;

                    const int firstSuccessor = firstOutgoingEdge(edge.endVertex);
                    if (firstSuccessor == -1) {
                        throw std::runtime_error("Contour boundary contains an edge without a successor.");
                    }
                    if (nextOutgoingEdges_[static_cast<std::size_t>(firstSuccessor)] == -1) {
                        currentEdgeIndex = firstSuccessor;
                    } else {
                        currentEdgeIndex =
                            selectSuccessor(edgeDirection(packedEdges[static_cast<std::size_t>(currentEdgeIndex)]), firstSuccessor,
                                            connectivity, packedEdges);
                        if (currentEdgeIndex == -1) {
                            throw std::runtime_error("Contour boundary contains an unresolved successor.");
                        }
                    }
                } while (currentEdgeIndex != startEdgeIndex);

                const std::size_t edgeCount = orderedEdges_.size() - firstBoundaryEdge;
                if (edgeCount != 0) {
                    boundaries.push_back(ContourBoundary{
                        doubledSignedArea >= 0 ? ContourBoundaryKind::External : ContourBoundaryKind::Internal,
                        checkedUint32(firstBoundaryEdge, "boundary edge offset"), checkedUint32(edgeCount, "boundary edge count"),
                        doubledSignedArea});
                }
            }
        } catch (...) {
            resetVertexIndex();
            throw;
        }
        resetVertexIndex();

        if (orderedEdges_.size() != packedEdges.size()) {
            throw std::runtime_error("Contour traversal did not cover every edge.");
        }
        packedEdges.swap(orderedEdges_);
    }

  private:
    /**
     * @brief Converts one checked size to the compact descriptor representation.
     * @param value Buffer size or offset.
     * @param context Operation name used in diagnostics.
     * @return `value` represented as `uint32_t`.
     */
    [[nodiscard]] static uint32_t checkedUint32(std::size_t value, const char* context) {
        if (value > static_cast<std::size_t>(std::numeric_limits<uint32_t>::max())) {
            throw std::overflow_error(std::string(context) + " exceeds uint32_t.");
        }
        return static_cast<uint32_t>(value);
    }

    /**
     * @brief Returns a row-major grid-vertex identifier.
     * @param row Vertex row.
     * @param column Vertex column.
     * @param vertexColumns Number of columns in the vertex grid.
     * @return Row-major vertex identifier.
     */
    [[nodiscard]] static int vertexId(int row, int column, int vertexColumns) {
        return (row * vertexColumns) + column;
    }

    /**
     * @brief Converts one packed edge to directed grid geometry.
     * @param packedEdge Packed contour edge.
     * @return Directed grid geometry.
     */
    [[nodiscard]] EdgeGeometry edgeGeometry(int packedEdge) const {
        const ContourEdge edge = unpackContourEdge(packedEdge);
        const int vertexColumns = columns_ + 1;
        const auto [row, column] = ImageUtils::to2D(edge.pixel, columns_);

        switch (edge.side) {
        case ContourSide::North:
            return {vertexId(row, column, vertexColumns), vertexId(row, column + 1, vertexColumns), -row};
        case ContourSide::East:
            return {vertexId(row, column + 1, vertexColumns), vertexId(row + 1, column + 1, vertexColumns), column + 1};
        case ContourSide::South:
            return {vertexId(row + 1, column + 1, vertexColumns), vertexId(row + 1, column, vertexColumns), row + 1};
        case ContourSide::West:
            return {vertexId(row + 1, column, vertexColumns), vertexId(row, column, vertexColumns), -column};
        }
        throw std::runtime_error("Invalid contour side.");
    }

    /**
     * @brief Returns the directed grid orientation encoded by a packed edge.
     * @param packedEdge Packed support-pixel edge.
     * @return Direction followed while traversing the edge.
     */
    [[nodiscard]] static Direction edgeDirection(int packedEdge) {
        switch (unpackContourEdge(packedEdge).side) {
        case ContourSide::North:
            return Direction::East;
        case ContourSide::East:
            return Direction::South;
        case ContourSide::South:
            return Direction::West;
        case ContourSide::West:
            return Direction::North;
        }
        throw std::runtime_error("Invalid contour side.");
    }

    /**
     * @brief Adds one directed edge to its start-vertex list.
     * @param startVertex Starting grid vertex.
     * @param edgeIndex Directed-edge index.
     */
    void addOutgoingEdge(int startVertex, int edgeIndex) {
        if (vertexIndex_ == ContourVertexIndex::Sparse) {
            auto [position, inserted] = sparseOutgoingEdgeHeads_.try_emplace(startVertex, -1);
            (void)inserted;
            nextOutgoingEdges_.push_back(position->second);
            position->second = edgeIndex;
            return;
        }

        const auto denseVertexIndex = static_cast<std::size_t>(startVertex);
        if (outgoingEdgeHeads_[denseVertexIndex] == -1) {
            touchedVertices_.push_back(startVertex);
        }
        nextOutgoingEdges_.push_back(outgoingEdgeHeads_[denseVertexIndex]);
        outgoingEdgeHeads_[denseVertexIndex] = edgeIndex;
    }

    /**
     * @brief Returns the first outgoing edge of one grid vertex.
     * @param vertex Grid vertex identifier.
     * @return First outgoing edge index, or -1 when absent.
     */
    [[nodiscard]] int firstOutgoingEdge(int vertex) const {
        if (vertexIndex_ == ContourVertexIndex::Sparse) {
            const auto position = sparseOutgoingEdgeHeads_.find(vertex);
            return position == sparseOutgoingEdgeHeads_.end() ? -1 : position->second;
        }
        return outgoingEdgeHeads_[static_cast<std::size_t>(vertex)];
    }

    /** @brief Clears the direct index entries touched by the preceding trace. */
    void resetVertexIndex() {
        if (vertexIndex_ == ContourVertexIndex::Sparse) {
            sparseOutgoingEdgeHeads_.clear();
            return;
        }
        for (int vertex : touchedVertices_) {
            outgoingEdgeHeads_[static_cast<std::size_t>(vertex)] = -1;
        }
        touchedVertices_.clear();
    }

    /**
     * @brief Returns whether one directed edge has already been traversed.
     * @param edgeIndex Directed-edge index.
     * @return True when the edge was visited.
     */
    [[nodiscard]] bool isVisited(int edgeIndex) const {
        return visitedEdges_[static_cast<std::size_t>(edgeIndex)] != 0;
    }

    /**
     * @brief Ranks an outgoing direction for the selected connectivity.
     * @param incoming Direction of the current edge.
     * @param outgoing Direction of a successor candidate.
     * @param connectivity Foreground connectivity at diagonal contacts.
     * @return Candidate priority, where a lower value is preferred.
     */
    [[nodiscard]] static int successorPriority(Direction incoming, Direction outgoing, ForegroundConnectivity connectivity) {
        const int quarterTurns = (static_cast<int>(outgoing) - static_cast<int>(incoming) + 4) & 3;
        const std::array<int, 4> priorities =
            connectivity == ForegroundConnectivity::Eight ? std::array<int, 4>{1, 2, 3, 0} : std::array<int, 4>{1, 0, 3, 2};
        return priorities[static_cast<std::size_t>(quarterTurns)];
    }

    /**
     * @brief Selects the next edge at a vertex with multiple outgoing edges.
     * @param incoming Direction of the current edge.
     * @param outgoingHead First successor candidate.
     * @param connectivity Foreground connectivity at diagonal contacts.
     * @param packedEdges Packed edges indexed like the temporary trace data.
     * @return Selected successor edge index.
     */
    [[nodiscard]] int selectSuccessor(Direction incoming, int outgoingHead, ForegroundConnectivity connectivity,
                                      const std::vector<int>& packedEdges) const {
        if (connectivity == ForegroundConnectivity::Unknown) {
            throw std::invalid_argument("Diagonal contour tracing requires canonical 4/8 construction adjacency; "
                                        "for complementary shapes, pass a valued-tree view with distinct node and parent altitudes.");
        }

        int bestEdgeIndex = -1;
        int bestPriority = std::numeric_limits<int>::max();
        for (int candidateIndex = outgoingHead; candidateIndex != -1;
             candidateIndex = nextOutgoingEdges_[static_cast<std::size_t>(candidateIndex)]) {
            const int priority =
                successorPriority(incoming, edgeDirection(packedEdges[static_cast<std::size_t>(candidateIndex)]), connectivity);
            if (priority < bestPriority ||
                (priority == bestPriority && packedEdges[static_cast<std::size_t>(candidateIndex)] <
                                                 packedEdges[static_cast<std::size_t>(bestEdgeIndex)])) {
                bestEdgeIndex = candidateIndex;
                bestPriority = priority;
            }
        }
        return bestEdgeIndex;
    }

    int columns_ = 0;                           ///< Number of image columns.
    ContourVertexIndex vertexIndex_;             ///< Active vertex-index representation.
    std::vector<TraceEdgeData> edgeData_;       ///< Geometry retained for the current contour edges.
    std::vector<int> outgoingEdgeHeads_;        ///< First outgoing edge for each grid vertex.
    std::unordered_map<int, int> sparseOutgoingEdgeHeads_; ///< Sparse first outgoing edge by grid vertex.
    std::vector<int> nextOutgoingEdges_;        ///< Successors in outgoing-edge lists.
    std::vector<int> touchedVertices_;          ///< Direct-index entries requiring reset.
    std::vector<uint8_t> visitedEdges_;          ///< Visited flags for the current trace.
    std::vector<int> orderedEdges_;              ///< Reusable output permutation buffer.
};

} // namespace mmcfilters::contours::detail
