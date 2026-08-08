#pragma once

#include "../utils/Image.hpp"
#include "../utils/Common.hpp"
#include "detail/NativeHierarchyValidationDetail.hpp"

#include <algorithm>
#include <array>
#include <concepts>
#include <cstdint>
#include <deque>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace mmcfilters {

/**
 * @brief Selects the interpolation/connectivity convention used to build a tree of shapes.
 */
enum class ToSInterpolation { SelfDual, Min4cMax8c, Min8cMax4c };

/**
 * @brief Selects whether the transient immersion domain has an exterior ring.
 *
 * This is strictly a producer policy. Both choices project the published
 * hierarchy back onto the original proper-part domain.
 */
enum class ToSPaddingPolicy { Exterior, NoPadding };

inline constexpr int ToSInterpolationScale = 2;
inline constexpr int ToSInterpolationPadding = 1;
inline constexpr int ToSUInt8Depth = 8;
inline constexpr int ToSSelfDualDepth = 9;
inline constexpr int ToSDefaultInfinityRow = 0;
inline constexpr int ToSDefaultInfinityCol = 0;

using ToSGrayLevel = uint16_t;
using ToSFloodDepth = uint32_t;

/**
 * @brief Runtime selection of built-in Tree-of-Shapes producer policies.
 *
 * Built-in interpolation dispatch happens once before the propagation loop.
 * The options are never persisted in `MorphologicalTree`.
 */
struct TreeOfShapesProducerOptions {
    /// Interpolation/connectivity convention used on the transient immersion.
    ToSInterpolation interpolation = ToSInterpolation::SelfDual;
    /// Whether the transient immersion contains an exterior ring.
    ToSPaddingPolicy padding = ToSPaddingPolicy::Exterior;
    /// Row of the propagation seed in the selected transient immersion.
    int infinitySeedRow = ToSDefaultInfinityRow;
    /// Column of the propagation seed in the selected transient immersion.
    int infinitySeedCol = ToSDefaultInfinityCol;
};

/**
 * @brief Compatibility encoding of construction levels in source gray units.
 */
struct ToSQuantizedUInt8AltitudeEncoding {
    /// Scalar type stored by the resulting weighted hierarchy.
    using value_type = std::uint8_t;

    /**
     * @brief Converts one transient construction level to source gray-level units.
     *
     * @param constructionLevel Altitude or level represented by `constructionLevel`.
     * @param interpolation Tree-of-Shapes interpolation policy.
     * @return The converted one transient construction level to source gray-level units.
     */
    static value_type encode(ToSGrayLevel constructionLevel, ToSInterpolation interpolation) noexcept {
        const ToSGrayLevel divisor =
            interpolation == ToSInterpolation::SelfDual ? static_cast<ToSGrayLevel>(ToSInterpolationScale) : static_cast<ToSGrayLevel>(1);
        return static_cast<value_type>(constructionLevel / divisor);
    }
};

/**
 * @brief Exact encoding in doubled source-level units.
 *
 * Original source levels are even. Odd values represent half levels generated
 * by the self-dual span interpolation.
 */
struct ToSExactDoubledAltitudeEncoding {
    /// Scalar type storing exact doubled source-level units.
    using value_type = ToSGrayLevel;

    /**
     * @brief Converts one transient construction level to exact doubled units.
     *
     * @param constructionLevel Altitude or level represented by `constructionLevel`.
     * @param interpolation Tree-of-Shapes interpolation policy.
     * @return The converted one transient construction level to exact doubled units.
     */
    static value_type encode(ToSGrayLevel constructionLevel, ToSInterpolation interpolation) noexcept {
        if (interpolation == ToSInterpolation::SelfDual) {
            return constructionLevel;
        }
        return static_cast<value_type>(ToSInterpolationScale * static_cast<unsigned int>(constructionLevel));
    }
};

template <class Policy>
concept ToSAltitudeEncodingPolicy = AltitudeValue<typename Policy::value_type> && requires(ToSGrayLevel level, ToSInterpolation interpolation) {
    { Policy::encode(level, interpolation) } -> std::same_as<typename Policy::value_type>;
};

/**
 * @brief Native representation produced by a Tree-of-Shapes producer.
 *
 * Unlike a pixel-parent image, this representation can preserve a virtual root
 * (or another branching projected node) whose direct proper part is empty.
 * The altitude type is selected by the producer encoding policy.
 */
template <AltitudeValue T> struct TreeOfShapesBuildResultT {
    /// Parent id for every produced internal node.
    std::vector<NodeId> nodeParent;
    /// Direct owning node for every original-domain proper part.
    std::vector<NodeId> properPartOwner;
    /// Encoded altitude for every produced internal node.
    std::vector<T> altitude;
    /// Root node id in the produced internal-node domain.
    NodeId root = InvalidNode;
    /// Number of rows in the original proper-part domain.
    int numRows = 0;
    /// Number of columns in the original proper-part domain.
    int numCols = 0;
    /// Move-only proof that the producer established the topology invariants.
    detail::NativeTopologyProof topologyProof;

    /**
     * @brief Transfers the producer-owned buffers together with their generic
     * topology proof.
     *
     * @param semantics Hierarchy semantics validated by the operation.
     * @return The transferred producer-owned buffers together with their generic topology proof.
     */
    [[nodiscard]] detail::ValidatedNativeHierarchy<T> takeValidatedHierarchy(HierarchySemantics semantics) && {
        return detail::makeValidatedNativeHierarchy<T>(std::move(nodeParent), std::move(properPartOwner), std::move(altitude), root,
                                                       GridDomain2D{numRows, numCols}, std::move(semantics), std::move(topologyProof));
    }
};

using TreeOfShapesBuildResult = TreeOfShapesBuildResultT<std::uint8_t>;
using ExactTreeOfShapesBuildResult = TreeOfShapesBuildResultT<ToSGrayLevel>;

/************************ Tree of Shapes support ************************/

/// @cond INTERNAL
namespace detail {

/*
 * Adaptive adjacency backend used by the tree-of-shapes construction.
 *
 * Diagonal links are activated on demand so that the interpolated grid can
 * emulate the required 4/8-connectivity behaviour during the union-find pass.
 * Each pixel may carry four diagonal flags: SW, NE, SE, and NW.
 */
enum class DiagonalConnection : uint8_t { None = 0, SW = 1 << 0, NE = 1 << 1, SE = 1 << 2, NW = 1 << 3 };

// Helper operators for diagonal-connection flags.
/**
 * @brief Combines two flag sets.
 *
 * @param a First operand.
 * @param b Second operand.
 * @return Combined flag set.
 */
inline DiagonalConnection operator|(DiagonalConnection a, DiagonalConnection b) {
    return static_cast<DiagonalConnection>(static_cast<uint8_t>(a) | static_cast<uint8_t>(b));
}

/**
 * @brief Adds flags to the left-hand flag set.
 *
 * @param a First operand.
 * @param b Second operand.
 * @return Reference to the updated flag set.
 */
inline DiagonalConnection& operator|=(DiagonalConnection& a, DiagonalConnection b) {
    a = a | b;
    return a;
}

/**
 * @brief Tests whether two flag sets overlap.
 *
 * @param a First operand.
 * @param b Second operand.
 * @return True when two flag sets overlap; otherwise false.
 */
inline bool operator&(DiagonalConnection a, DiagonalConnection b) { return static_cast<uint8_t>(a) & static_cast<uint8_t>(b); }

/**
 * @brief Adaptive adjacency used during 4/8-connected tree-of-shapes construction.
 */
class AdjacencyUC {
  private:
    /** @brief Stores the number of rows in the image domain. */
    int numRows;
    /** @brief Stores the number of columns in the image domain. */
    int numCols;
    /** @brief Stores the enabled diagonal-connection flags for each pixel. */
    std::vector<uint8_t> dconnFlags; // 4-connect.  +  diag. connect.
                                     //  N, W, S, E,   SW, NE, SE, NW
    /** @brief Stores row offsets for cardinal and diagonal neighbours. */
    const std::vector<int> offsetRows = {-1, 0, 1, 0, 1, -1, 1, -1};
    /** @brief Stores column offsets for cardinal and diagonal neighbours. */
    const std::vector<int> offsetCols = {0, -1, 0, 1, -1, 1, 1, -1};
    /** @brief Indicates whether diagonal connections may be traversed. */
    bool enableDiagonalConnection;
    /** @brief Maps each diagonal offset to its required connection flag. */
    const std::vector<DiagonalConnection> requiredDiagonal = {DiagonalConnection::SW, DiagonalConnection::NE, DiagonalConnection::SE, DiagonalConnection::NW};

  public:
    /**
     * @brief Constructs `AdjacencyUC` from the supplied inputs.
     *
     * @param rows Number of rows in the domain.
     * @param cols Number of columns in the domain.
     * @param enableDiagonalConnection Flag controlling enable diagonal connection.
     */
    AdjacencyUC(int rows, int cols, bool enableDiagonalConnection) : numRows(rows), numCols(cols), enableDiagonalConnection(enableDiagonalConnection) {
        if (enableDiagonalConnection)
            dconnFlags.resize(rows * cols, 0);
    }

    /**
     * @brief Destroys `AdjacencyUC`.
     */
    ~AdjacencyUC() {}

    /**
     * @brief Sets diagonal connection.
     *
     * @param row Zero-based row coordinate.
     * @param col Zero-based column coordinate.
     * @param conn Diagonal-connection flag.
     */
    void setDiagonalConnection(int row, int col, DiagonalConnection conn) { dconnFlags[ImageUtils::to1D(row, col, numCols)] |= static_cast<uint8_t>(conn); }

    /**
     * @brief Sets diagonal connection.
     *
     * @param idx Zero-based index used by the operation.
     * @param conn Diagonal-connection flag.
     */
    void setDiagonalConnection(int idx, DiagonalConnection conn) { dconnFlags[idx] |= static_cast<uint8_t>(conn); }

    /**
     * @brief Tests whether connection holds.
     *
     * @param row Zero-based row coordinate.
     * @param col Zero-based column coordinate.
     * @param conn Diagonal-connection flag.
     * @return True when connection; otherwise false.
     */
    bool hasConnection(int row, int col, DiagonalConnection conn) const { return dconnFlags[ImageUtils::to1D(row, col, numCols)] & static_cast<uint8_t>(conn); }

    /**
     * @brief Returns connections.
     *
     * @param row Zero-based row coordinate.
     * @param col Zero-based column coordinate.
     * @return Connections.
     */
    uint8_t getConnections(int row, int col) const { return dconnFlags[ImageUtils::to1D(row, col, numCols)]; }

    /**
     * @brief Iterator over valid neighbours including enabled diagonal links.
     */
    class NeighborIterator {
      private:
        /** @brief References the adaptive adjacency being traversed. */
        AdjacencyUC& instance;
        /** @brief Stores the source pixel row. */
        int row;
        /** @brief Stores the source pixel column. */
        int col;
        /** @brief Stores the current neighbour-offset index. */
        std::size_t id;

        /**
         * @brief Advances to valid.
         */
        void advanceToValid() {
            while (id < instance.offsetRows.size()) {
                int r = row + instance.offsetRows[id];
                int c = col + instance.offsetCols[id];
                if (r >= 0 && c >= 0 && r < instance.numRows && c < instance.numCols) {
                    if (id < 4 || (instance.enableDiagonalConnection && instance.dconnFlags[ImageUtils::to1D(row, col, instance.numCols)] &
                                                                            static_cast<uint8_t>(instance.requiredDiagonal[id - 4]))) {
                        return;
                    }
                }
                ++id;
            }
        }

      public:
        /**
         * @brief Constructs `NeighborIterator` from the supplied inputs.
         *
         * @param adj Adjacency relation used to traverse the image domain.
         * @param row Zero-based row coordinate.
         * @param col Zero-based column coordinate.
         * @param id Identifier used by the operation.
         */
        NeighborIterator(AdjacencyUC& adj, int row, int col, int id) : instance(adj), row(row), col(col), id(id) { advanceToValid(); }

        /**
         * @brief Returns the value at the current iterator position.
         *
         * @return The value at the current iterator position.
         */
        int operator*() const {
            int dr = instance.offsetRows[id];
            int dc = instance.offsetCols[id];
            return ImageUtils::to1D(row + dr, col + dc, instance.numCols);
        }

        /**
         * @brief Advances the iterator.
         *
         * @return Reference or value representing the advanced iterator.
         */
        NeighborIterator& operator++() {
            ++id;
            advanceToValid();
            return *this;
        }

        /**
         * @brief Tests whether two objects are equal.
         *
         * @param other Object compared with or copied from.
         * @return True when two objects are equal; otherwise false.
         */
        bool operator==(const NeighborIterator& other) const { return id == other.id; }

        /**
         * @brief Tests whether two objects differ.
         *
         * @param other Object compared with or copied from.
         * @return True when two objects differ; otherwise false.
         */
        bool operator!=(const NeighborIterator& other) const { return !(*this == other); }
    };

    /**
     * @brief Range helper that produces valid-neighbour iterators.
     */
    class NeighborRange {
      private:
        /** @brief References the adaptive adjacency represented by the range. */
        AdjacencyUC& instance;
        /** @brief Stores the source pixel row. */
        int row;
        /** @brief Stores the source pixel column. */
        int col;

      public:
        /**
         * @brief Constructs `NeighborRange` from the supplied inputs.
         *
         * @param instance Producer instance whose priority state is queried.
         * @param row Zero-based row coordinate.
         * @param col Zero-based column coordinate.
         */
        NeighborRange(AdjacencyUC& instance, int row, int col) : instance(instance), row(row), col(col) {}

        /**
         * @brief Returns an iterator positioned at the beginning.
         *
         * @return An iterator positioned at the beginning.
         */
        NeighborIterator begin() { return NeighborIterator(instance, row, col, 0); }
        /**
         * @brief Returns the exclusive-end iterator.
         *
         * @return The exclusive-end iterator.
         */
        NeighborIterator end() { return NeighborIterator(instance, row, col, 8); }
    };

    /**
     * @brief Returns neighbor indices.
     *
     * @param p Row-major pixel index to transform or inspect.
     * @return Neighbor indices.
     */
    NeighborRange getNeighborIndices(int p) {
        auto [row, col] = ImageUtils::to2D(p, numCols);
        return NeighborRange(*this, row, col);
    }

    /**
     * @brief Returns neighbor indices.
     *
     * @param row Zero-based row coordinate.
     * @param col Zero-based column coordinate.
     * @return Neighbor indices.
     */
    NeighborRange getNeighborIndices(int row, int col) { return NeighborRange(*this, row, col); }
};

/**
 * @brief Discrete priority queue used during tree-of-shapes construction.
 */
class PriorityQueueToS {
  private:
    /** @brief Stores queued pixels grouped by discrete priority. */
    std::vector<std::deque<int>> buckets;
    /** @brief Stores the priority of the next non-empty bucket. */
    int currentPriority;
    /** @brief Stores the total number of queued elements. */
    int numElements;
    /** @brief Stores the number of allocated priority levels. */
    int maxPriorityLevels;

  public:
    /**
     * @brief Constructs `PriorityQueueToS` from the supplied inputs.
     *
     * @param depthOfImage Bit depth used to size the priority buckets.
     */
    PriorityQueueToS(int depthOfImage = ToSUInt8Depth) : currentPriority(0), numElements(0), maxPriorityLevels(1 << depthOfImage) {
        buckets.resize(maxPriorityLevels);
    }

    /**
     * @brief Initializes the priority queue with its first element.
     *
     * @param element Element handled by the container.
     * @param priority Priority associated with the element.
     */
    void initial(int element, int priority) {
        currentPriority = priority;
        buckets[priority].push_back(element);
        numElements++;
    }
    /**
     * @brief Returns current priority.
     *
     * @return Current priority.
     */
    int getCurrentPriority() { return currentPriority; }
    /**
     * @brief Tests whether empty holds.
     *
     * @return True when empty; otherwise false.
     */
    bool isEmpty() { return numElements == 0; }

    /**
     * @brief Pushes an element into the priority bucket selected from its interval.
     *
     * @param element Element handled by the container.
     * @param lower Lower admissible priority.
     * @param upper Upper admissible priority.
     */
    void priorityPush(int element, int lower, int upper) {
        int priority;
        if (lower > currentPriority) {
            priority = lower;
        } else if (upper < currentPriority) {
            priority = upper;
        } else {
            priority = currentPriority;
        }
        numElements++;
        buckets[priority].push_back(element);
    }

    /**
     * @brief Removes the next element from the circular priority buckets.
     *
     * @return Identifier of the removed element.
     */
    int priorityPop() {
        if (buckets[currentPriority].empty()) {
            int nextPriority = -1;
            for (int distance = 1; distance < maxPriorityLevels; ++distance) {
                const int lowerPriority = currentPriority - distance;
                if (lowerPriority >= 0 && !buckets[lowerPriority].empty()) {
                    nextPriority = lowerPriority;
                    break;
                }

                const int upperPriority = currentPriority + distance;
                if (upperPriority < maxPriorityLevels && !buckets[upperPriority].empty()) {
                    nextPriority = upperPriority;
                    break;
                }
            }

            if (nextPriority == -1) {
                throw std::runtime_error("PriorityQueueToS is empty.");
            }
            currentPriority = nextPriority;
        }

        int element = buckets[currentPriority].front();
        buckets[currentPriority].pop_front();

        numElements--;
        return element;
    }
};

} // namespace detail

/// @endcond

/**
 * @brief Builds trees of shapes (ToS) using a union-find construction.
 *
 * The builder interpolates the input image according to `ToSInterpolation`,
 * constructs a max-tree on the interpolated domain, and projects the resulting
 * hierarchy back to the original image domain.
 */
class TreeOfShapesProducer {
  private:
    /** @brief Defines the `AdjacencyUC` alias used by the component. */
    using AdjacencyUC = detail::AdjacencyUC;
    /** @brief Defines the `DiagonalConnection` alias used by the component. */
    using DiagonalConnection = detail::DiagonalConnection;
    /** @brief Defines the `PriorityQueueToS` alias used by the component. */
    using PriorityQueueToS = detail::PriorityQueueToS;

    /** @brief Stores the intermediate topology produced by the tree-of-shapes flood. */
    struct FloodResult {
        /** @brief Stores the tree level. */
        std::vector<ToSFloodDepth> treeLevel;
        /** @brief Stores the gray level. */
        std::vector<ToSGrayLevel> grayLevel;
        /** @brief Stores the order. */
        std::vector<int> order;
        /** @brief Stores the adjacency. */
        AdjacencyUC adjacency;
    };

    /** @brief Stores the canonical flood-tree topology and interpolation metadata. */
    struct CanonicalFloodTree {
        /** @brief Stores the tree level. */
        std::vector<ToSFloodDepth> treeLevel;
        /** @brief Stores the gray level. */
        std::vector<ToSGrayLevel> grayLevel;
        /** @brief Stores the order. */
        std::vector<int> order;
        /** @brief Stores the parent. */
        std::vector<int> parent;
    };

    /** @brief Stores the options. */
    TreeOfShapesProducerOptions options_;

    /**
     * @brief Validates options.
     *
     * @param options Options controlling the operation.
     */
    static void validateOptions(const TreeOfShapesProducerOptions& options) {
        switch (options.interpolation) {
        case ToSInterpolation::SelfDual:
        case ToSInterpolation::Min4cMax8c:
        case ToSInterpolation::Min8cMax4c:
            break;
        default:
            throw std::invalid_argument("Unsupported tree-of-shapes interpolation policy.");
        }
        switch (options.padding) {
        case ToSPaddingPolicy::Exterior:
        case ToSPaddingPolicy::NoPadding:
            break;
        default:
            throw std::invalid_argument("Unsupported tree-of-shapes padding policy.");
        }
    }

    /**
     * @brief Tests whether connectivity map holds.
     *
     * @return True when connectivity map; otherwise false.
     */
    inline bool usesConnectivityMap() const noexcept {
        return options_.interpolation == ToSInterpolation::Min4cMax8c || options_.interpolation == ToSInterpolation::Min8cMax4c;
    }

    /**
     * @brief Tests whether high diagonal at saddle holds.
     *
     * @return True when high diagonal at saddle; otherwise false.
     */
    inline bool usesHighDiagonalAtSaddle() const noexcept { return options_.interpolation == ToSInterpolation::Min4cMax8c; }

    /**
     * @brief Tests whether exterior padding holds.
     *
     * @return True when exterior padding; otherwise false.
     */
    inline bool usesExteriorPadding() const noexcept { return options_.padding == ToSPaddingPolicy::Exterior; }

    /**
     * @brief Returns the padding offset applied to the interpolated domain.
     *
     * @return Padding offset in interpolated-domain coordinates.
     */
    inline int interpolationOffset() const noexcept { return usesExteriorPadding() ? ToSInterpolationPadding : 0; }

    /**
     * @brief Checks and converts interpolated extent.
     *
     * @param extent Number of entries in the interpolated ordering.
     * @return Validated extent of the interpolated dimension.
     */
    inline int checkedInterpolatedExtent(int extent) const {
        if (extent <= 0) {
            throw std::invalid_argument("TreeOfShapesProducer requires positive image dimensions.");
        }
        const std::int64_t interpolated = ToSInterpolationScale * static_cast<std::int64_t>(extent) + (usesExteriorPadding() ? 1 : -1);
        if (interpolated <= 0 || interpolated > std::numeric_limits<int>::max()) {
            throw std::overflow_error("Tree-of-shapes interpolated dimension exceeds int range.");
        }
        return static_cast<int>(interpolated);
    }

    /**
     * @brief Computes the number of rows in the interpolated domain.
     *
     * @param numRows Number of rows in the domain.
     * @return Validated extent of the interpolated dimension.
     */
    inline int interpolatedNumRows(int numRows) const { return checkedInterpolatedExtent(numRows); }

    /**
     * @brief Computes the number of columns in the interpolated domain.
     *
     * @param numCols Number of columns in the domain.
     * @return Validated extent of the interpolated dimension.
     */
    inline int interpolatedNumCols(int numCols) const { return checkedInterpolatedExtent(numCols); }

    /**
     * @brief Checks and converts domain size.
     *
     * @param numRows Number of rows in the domain.
     * @param numCols Number of columns in the domain.
     * @return Validated number of pixels in the original domain.
     */
    inline int checkedDomainSize(int numRows, int numCols) const {
        if (numRows <= 0 || numCols <= 0 || numRows > std::numeric_limits<int>::max() / numCols) {
            throw std::overflow_error("Tree-of-shapes domain size exceeds int range.");
        }
        return numRows * numCols;
    }

    /**
     * @brief Maps point row.
     *
     * @param row Zero-based row coordinate.
     * @return Mapped point row.
     */
    inline int originalPointRow(int row) const noexcept { return ToSInterpolationScale * row + interpolationOffset(); }

    /**
     * @brief Maps point col.
     *
     * @param col Zero-based column coordinate.
     * @return Mapped point col.
     */
    inline int originalPointCol(int col) const noexcept { return ToSInterpolationScale * col + interpolationOffset(); }

    /**
     * @brief Scales original level.
     *
     * @param value Value used by the operation.
     * @return Original gray level expressed in the interpolation scale.
     */
    inline ToSGrayLevel scaledOriginalLevel(uint8_t value) const noexcept { return static_cast<ToSGrayLevel>(ToSInterpolationScale * static_cast<int>(value)); }

    /**
     * @brief Computes the configured infinity-seed index in the interpolated domain.
     *
     * @param interpNumRows Number of rows in the interpolated domain.
     * @param interpNumCols Number of columns in the interpolated domain.
     * @return Row-major index of the configured infinity seed.
     */
    inline int infinitySeedIndex(int interpNumRows, int interpNumCols) const {
        if (options_.infinitySeedRow < 0 || options_.infinitySeedCol < 0 || options_.infinitySeedRow >= interpNumRows ||
            options_.infinitySeedCol >= interpNumCols) {
            throw std::invalid_argument("Tree-of-shapes infinity seed must be inside the interpolated domain.");
        }
        return ImageUtils::to1D(options_.infinitySeedRow, options_.infinitySeedCol, interpNumCols);
    }

    /**
     * @brief Sets diagonal0 connection.
     *
     * @param adj Adjacency relation used to traverse the image domain.
     * @param row Zero-based row coordinate.
     * @param col Zero-based column coordinate.
     */
    inline void setDiagonal0Connection(AdjacencyUC& adj, int row, int col) const {
        adj.setDiagonalConnection(row, col - 1, DiagonalConnection::SE);
        adj.setDiagonalConnection(row + 1, col, DiagonalConnection::NW);

        adj.setDiagonalConnection(row - 1, col - 1, DiagonalConnection::SE);
        adj.setDiagonalConnection(row, col, DiagonalConnection::SE | DiagonalConnection::NW);
        adj.setDiagonalConnection(row + 1, col + 1, DiagonalConnection::NW);

        adj.setDiagonalConnection(row - 1, col, DiagonalConnection::SE);
        adj.setDiagonalConnection(row, col + 1, DiagonalConnection::NW);
    }

    /**
     * @brief Sets diagonal1 connection.
     *
     * @param adj Adjacency relation used to traverse the image domain.
     * @param row Zero-based row coordinate.
     * @param col Zero-based column coordinate.
     */
    inline void setDiagonal1Connection(AdjacencyUC& adj, int row, int col) const {
        adj.setDiagonalConnection(row, col - 1, DiagonalConnection::NE);
        adj.setDiagonalConnection(row - 1, col, DiagonalConnection::SW);

        adj.setDiagonalConnection(row - 1, col + 1, DiagonalConnection::SW);
        adj.setDiagonalConnection(row, col, DiagonalConnection::SW | DiagonalConnection::NE);
        adj.setDiagonalConnection(row + 1, col - 1, DiagonalConnection::NE);

        adj.setDiagonalConnection(row + 1, col, DiagonalConnection::NE);
        adj.setDiagonalConnection(row, col + 1, DiagonalConnection::SW);
    }

    std::tuple<std::vector<ToSGrayLevel>, std::vector<ToSGrayLevel>, AdjacencyUC>
    /**
     * @brief Crops exterior interpolation.
     *
     * @param paddedMin Minimum values in the padded interpolation domain.
     * @param paddedMax Maximum values in the padded interpolation domain.
     * @param paddedAdjacency Adjacency relation defined on the padded domain.
     * @param paddedRows Number of rows in the padded domain.
     * @param paddedCols Number of columns in the padded domain.
     * @param adaptiveDiagonal Whether diagonal connections are chosen adaptively.
     * @return Values produced by the operation.
     */
    cropExteriorInterpolation(std::vector<ToSGrayLevel> paddedMin, std::vector<ToSGrayLevel> paddedMax, AdjacencyUC paddedAdjacency, int paddedRows,
                              int paddedCols, bool adaptiveDiagonal) const {
        if (paddedRows < 3 || paddedCols < 3) {
            throw std::logic_error("Tree-of-shapes padded immersion cannot be cropped.");
        }
        const int rows = paddedRows - 2;
        const int cols = paddedCols - 2;
        const int size = checkedDomainSize(rows, cols);
        std::vector<ToSGrayLevel> croppedMin(static_cast<std::size_t>(size));
        std::vector<ToSGrayLevel> croppedMax(static_cast<std::size_t>(size));
        AdjacencyUC croppedAdjacency(rows, cols, adaptiveDiagonal);

        for (int row = 0; row < rows; ++row) {
            for (int col = 0; col < cols; ++col) {
                const int source = ImageUtils::to1D(row + 1, col + 1, paddedCols);
                const int target = ImageUtils::to1D(row, col, cols);
                croppedMin[static_cast<std::size_t>(target)] = paddedMin[static_cast<std::size_t>(source)];
                croppedMax[static_cast<std::size_t>(target)] = paddedMax[static_cast<std::size_t>(source)];
                if (adaptiveDiagonal) {
                    const std::uint8_t connections = paddedAdjacency.getConnections(row + 1, col + 1);
                    if (connections != 0) {
                        croppedAdjacency.setDiagonalConnection(target, static_cast<DiagonalConnection>(connections));
                    }
                }
            }
        }
        return {std::move(croppedMin), std::move(croppedMax), std::move(croppedAdjacency)};
    }

  public:
    /**
     * @brief Creates a tree-of-shapes builder.
     *
     * @param interpolation Interpolation/connectivity convention.
     * @param infinitySeedRow Row of the propagation seed in the interpolated domain.
     * @param infinitySeedCol Column of the propagation seed in the interpolated domain.
     */
    explicit TreeOfShapesProducer(ToSInterpolation interpolation, int infinitySeedRow = ToSDefaultInfinityRow, int infinitySeedCol = ToSDefaultInfinityCol)
        : options_{interpolation, ToSPaddingPolicy::Exterior, infinitySeedRow, infinitySeedCol} {
        validateOptions(options_);
    }

    /**
     * @brief Creates a producer from the complete runtime policy set.
     *
     * @param options Policy options controlling the operation.
     */
    explicit TreeOfShapesProducer(TreeOfShapesProducerOptions options = {}) : options_(options) { validateOptions(options_); }

    /**
     * @brief Returns the immutable policies selected for this producer.
     *
     * @return The immutable policies selected for this producer.
     */
    [[nodiscard]] const TreeOfShapesProducerOptions& options() const noexcept { return options_; }

    /**
     * @brief Destroys `TreeOfShapesProducer`.
     */
    ~TreeOfShapesProducer() = default;

    /// @cond INTERNAL
    /**
     * @brief Interpolates image.
     *
     * @param imgPtr Image data represented by `imgPtr`.
     * @return Values produced by the operation.
     */
    std::tuple<std::vector<ToSGrayLevel>, std::vector<ToSGrayLevel>, AdjacencyUC> interpolateImage(const ImageUInt8Ptr& imgPtr) const {
        // Implements the self-dual span-based immersion used by Boutry's PhD
        // thesis: ISpan(u) followed by front propagation from a median-valued
        // outer boundary. Internal levels are stored in Z/2 by scaling values
        // by 2, so odd integers represent half gray levels.
        if (!imgPtr) {
            throw std::invalid_argument("TreeOfShapesProducer requires a non-null image.");
        }
        if (!usesExteriorPadding()) {
            TreeOfShapesProducerOptions paddedOptions = options_;
            paddedOptions.padding = ToSPaddingPolicy::Exterior;
            TreeOfShapesProducer paddedProducer(paddedOptions);
            auto [paddedMin, paddedMax, paddedAdjacency] = paddedProducer.interpolateImage(imgPtr);
            const int paddedRows = paddedProducer.interpolatedNumRows(imgPtr->getNumRows());
            const int paddedCols = paddedProducer.interpolatedNumCols(imgPtr->getNumCols());
            return cropExteriorInterpolation(std::move(paddedMin), std::move(paddedMax), std::move(paddedAdjacency), paddedRows, paddedCols, false);
        }
        auto img = imgPtr->rawData();
        int numRows = imgPtr->getNumRows();
        int numCols = imgPtr->getNumCols();
        if (numRows <= 0 || numCols <= 0) {
            throw std::invalid_argument("TreeOfShapesProducer requires a non-empty image.");
        }
        if (numRows == 1 && numCols == 1) {
            const int interpNumRows = interpolatedNumRows(numRows);
            const int interpNumCols = interpolatedNumCols(numCols);
            const int interpSize = checkedDomainSize(interpNumRows, interpNumCols);
            std::vector<ToSGrayLevel> interpolationMin(static_cast<std::size_t>(interpSize), scaledOriginalLevel(img[0]));
            std::vector<ToSGrayLevel> interpolationMax(static_cast<std::size_t>(interpSize), scaledOriginalLevel(img[0]));
            AdjacencyUC adj(interpNumRows, interpNumCols, false);
            return std::make_tuple(std::move(interpolationMin), std::move(interpolationMax), std::move(adj));
        }

        constexpr int adjCircleCol[] = {-1, +1, -1, +1};
        constexpr int adjCircleRow[] = {-1, -1, +1, +1};

        constexpr int adjRetHorCol[] = {0, 0};
        constexpr int adjRetHorRow[] = {-1, +1};

        constexpr int adjRetVerCol[] = {+1, -1};
        constexpr int adjRetVerRow[] = {0, 0};

        int interpNumCols = interpolatedNumCols(numCols);
        int interpNumRows = interpolatedNumRows(numRows);
        int size = checkedDomainSize(interpNumRows, interpNumCols);

        // Allocate interpolation result buffers for minimum and maximum levels.
        std::vector<ToSGrayLevel> interpolationMin(size);
        std::vector<ToSGrayLevel> interpolationMax(size);

        // Collect each boundary pixel exactly once. The closed-form perimeter
        // size 2 * (rows + cols) - 4 is valid only when both dimensions are at
        // least two and previously left zero-filled entries for thin images.
        std::array<int, 256> boundaryHistogram{};
        int numBoundary = 0;

        int pT;

        for (int p = 0; p < numCols * numRows; p++) {
            auto [row, col] = ImageUtils::to2D(p, numCols);

            // Check whether the pixel lies on the image border.
            if (row == 0 || row == numRows - 1 || col == 0 || col == numCols - 1) {
                ++boundaryHistogram[static_cast<std::size_t>(img[p])];
                ++numBoundary;
            }

            // Compute the interpolated-image index.
            pT = ImageUtils::to1D(originalPointRow(row), originalPointCol(col), interpNumCols);

            // Assign interpolation values.
            interpolationMin[pT] = interpolationMax[pT] = scaledOriginalLevel(img[p]);
        }

        auto boundaryValueAtRank = [&](int rank) {
            int cumulative = 0;
            for (int value = 0; value < 256; ++value) {
                cumulative += boundaryHistogram[static_cast<std::size_t>(value)];
                if (rank < cumulative) {
                    return value;
                }
            }
            throw std::runtime_error("Tree-of-shapes boundary histogram is inconsistent.");
        };
        int median;
        if (numBoundary % 2 == 0) {
            median = boundaryValueAtRank(numBoundary / 2 - 1) + boundaryValueAtRank(numBoundary / 2);
        } else {
            median = ToSInterpolationScale * boundaryValueAtRank(numBoundary / 2);
        }
        // std::cout << "Interpolation (Median): " << median << std::endl;

        int qT, qCol, qRow, min, max;
        const int* adjCol = nullptr;
        const int* adjRow = nullptr;
        int adjSize;
        AdjacencyUC adj(interpNumRows, interpNumCols, false);

        for (int row = 0; row < interpNumRows; row++) {
            for (int col = 0; col < interpNumCols; col++) {
                if (col % 2 == 1 && row % 2 == 1)
                    continue;
                pT = ImageUtils::to1D(row, col, interpNumCols);
                if (col == 0 || col == interpNumCols - 1 || row == 0 || row == interpNumRows - 1) {
                    max = median;
                    min = median;
                } else {
                    if (col % 2 == 0 && row % 2 == 0) {
                        adjCol = adjCircleCol;
                        adjRow = adjCircleRow;
                        adjSize = 4;
                    } else if (col % 2 == 0 && row % 2 == 1) {
                        adjCol = adjRetVerCol;
                        adjRow = adjRetVerRow;
                        adjSize = 2;
                    } else if (col % 2 == 1 && row % 2 == 0) {
                        adjCol = adjRetHorCol;
                        adjRow = adjRetHorRow;
                        adjSize = 2;
                    }

                    min = std::numeric_limits<int>::max();
                    max = std::numeric_limits<int>::min();
                    for (int i = 0; i < adjSize; i++) {
                        qRow = row + adjRow[i];
                        qCol = col + adjCol[i];

                        if (qRow >= 0 && qCol >= 0 && qRow < interpNumRows && qCol < interpNumCols) {
                            qT = ImageUtils::to1D(qRow, qCol, interpNumCols);

                            if (interpolationMax[qT] > max) {
                                max = interpolationMax[qT];
                            }
                            if (interpolationMin[qT] < min) {
                                min = interpolationMin[qT];
                            }
                        } else {
                            if (median > max) {
                                max = median;
                            }
                            if (median < min) {
                                min = median;
                            }
                        }
                    }
                }
                interpolationMin[pT] = static_cast<ToSGrayLevel>(min);
                interpolationMax[pT] = static_cast<ToSGrayLevel>(max);
            }
        }
        return std::make_tuple(std::move(interpolationMin), std::move(interpolationMax), std::move(adj));
    }

    /**
     * @brief Interpolates image4c8c.
     *
     * @param imgPtr Image data represented by `imgPtr`.
     * @return Values produced by the operation.
     */
    std::tuple<std::vector<ToSGrayLevel>, std::vector<ToSGrayLevel>, AdjacencyUC> interpolateImage4c8c(const ImageUInt8Ptr& imgPtr) const {
        // Implements the optimized 2D immersion/connectivity-map rules from
        // Carlinet, Crozet, and Geraud, "The Tree of Shapes Turned into a
        // Max-Tree: A Simple and Efficient Linear Algorithm", ICIP 2018.
        if (!imgPtr) {
            throw std::invalid_argument("TreeOfShapesProducer requires a non-null image.");
        }
        if (!usesExteriorPadding()) {
            TreeOfShapesProducerOptions paddedOptions = options_;
            paddedOptions.padding = ToSPaddingPolicy::Exterior;
            TreeOfShapesProducer paddedProducer(paddedOptions);
            auto [paddedMin, paddedMax, paddedAdjacency] = paddedProducer.interpolateImage4c8c(imgPtr);
            const int paddedRows = paddedProducer.interpolatedNumRows(imgPtr->getNumRows());
            const int paddedCols = paddedProducer.interpolatedNumCols(imgPtr->getNumCols());
            return cropExteriorInterpolation(std::move(paddedMin), std::move(paddedMax), std::move(paddedAdjacency), paddedRows, paddedCols, true);
        }
        auto img = imgPtr->rawData();
        int numRows = imgPtr->getNumRows();
        int numCols = imgPtr->getNumCols();
        if (numRows <= 0 || numCols <= 0) {
            throw std::invalid_argument("TreeOfShapesProducer requires a non-empty image.");
        }

        int interpNumCols = interpolatedNumCols(numCols);
        int interpNumRows = interpolatedNumRows(numRows);
        int size = checkedDomainSize(interpNumRows, interpNumCols);
        AdjacencyUC adj(interpNumRows, interpNumCols, true);

        // Allocate interpolation result buffers for minimum and maximum levels.
        std::vector<ToSGrayLevel> interpolationMin(size);
        std::vector<ToSGrayLevel> interpolationMax(size);

        int pT;
        // Compute interval from 2-faces.
        for (int p = 0; p < numCols * numRows; p++) {
            auto [row, col] = ImageUtils::to2D(p, numCols);

            // Compute the interpolated-image index.
            pT = ImageUtils::to1D(originalPointRow(row), originalPointCol(col), interpNumCols);

            // Assign interpolation values.
            interpolationMin[pT] = interpolationMax[pT] = img[p];
        }

        auto getValue = [&](int row, int col) -> int {
            int origRow = (row - 1) / 2;
            int origCol = (col - 1) / 2;
            return img[ImageUtils::to1D(origRow, origCol, numCols)];
        };

        // Borders.
        for (int row = 0; row < interpNumRows; row++) {
            int col;
            if (row % 2 == 1) { // horizontal e vertical
                col = 0;
                int v1 = getValue(row, col + 1);
                interpolationMin[ImageUtils::to1D(row, col, interpNumCols)] = v1;
                interpolationMax[ImageUtils::to1D(row, col, interpNumCols)] = v1;

                col = interpNumCols - 1;
                v1 = getValue(row, col - 1);
                interpolationMin[ImageUtils::to1D(row, col, interpNumCols)] = v1;
                interpolationMax[ImageUtils::to1D(row, col, interpNumCols)] = v1;
            } else { // circulos
                if (row == 0) {
                    col = 0;
                    int v1 = getValue(row + 1, col + 1);
                    interpolationMin[ImageUtils::to1D(row, col, interpNumCols)] = v1;
                    interpolationMax[ImageUtils::to1D(row, col, interpNumCols)] = v1;

                    col = interpNumCols - 1;
                    v1 = getValue(row + 1, col - 1);
                    interpolationMin[ImageUtils::to1D(row, col, interpNumCols)] = v1;
                    interpolationMax[ImageUtils::to1D(row, col, interpNumCols)] = v1;

                } else if (row == interpNumRows - 1) {
                    col = 0;
                    int v1 = getValue(row - 1, 1);
                    interpolationMin[ImageUtils::to1D(row, col, interpNumCols)] = v1;
                    interpolationMax[ImageUtils::to1D(row, col, interpNumCols)] = v1;

                    col = interpNumCols - 1;
                    v1 = getValue(row - 1, col - 1);
                    interpolationMin[ImageUtils::to1D(row, col, interpNumCols)] = v1;
                    interpolationMax[ImageUtils::to1D(row, col, interpNumCols)] = v1;
                } else {
                    col = 0;
                    int v1 = getValue(row - 1, col + 1);
                    int v2 = getValue(row + 1, col + 1);
                    interpolationMin[ImageUtils::to1D(row, 0, interpNumCols)] = std::min(v1, v2);
                    interpolationMax[ImageUtils::to1D(row, 0, interpNumCols)] = std::max(v1, v2);

                    col = interpNumCols - 1;
                    v1 = getValue(row - 1, col - 1);
                    v2 = getValue(row + 1, col - 1);
                    interpolationMin[ImageUtils::to1D(row, col, interpNumCols)] = std::min(v1, v2);
                    interpolationMax[ImageUtils::to1D(row, col, interpNumCols)] = std::max(v1, v2);
                }
            }
        }

        for (int col = 1; col < interpNumCols - 1; col++) {
            int row;
            if (col % 2 == 1) { // horizontal e vertical
                row = 0;
                int v1 = getValue(row + 1, col);
                interpolationMin[ImageUtils::to1D(row, col, interpNumCols)] = v1;
                interpolationMax[ImageUtils::to1D(row, col, interpNumCols)] = v1;

                row = interpNumRows - 1;
                v1 = getValue(row - 1, col);
                interpolationMin[ImageUtils::to1D(row, col, interpNumCols)] = v1;
                interpolationMax[ImageUtils::to1D(row, col, interpNumCols)] = v1;
            } else { // circulos
                row = 0;
                int v1 = getValue(row + 1, col - 1);
                int v2 = getValue(row + 1, col + 1);
                interpolationMin[ImageUtils::to1D(row, col, interpNumCols)] = std::min(v1, v2);
                interpolationMax[ImageUtils::to1D(row, col, interpNumCols)] = std::max(v1, v2);

                row = interpNumRows - 1;
                v1 = getValue(row - 1, col - 1);
                v2 = getValue(row - 1, col + 1);
                interpolationMin[ImageUtils::to1D(row, col, interpNumCols)] = std::min(v1, v2);
                interpolationMax[ImageUtils::to1D(row, col, interpNumCols)] = std::max(v1, v2);
            }
        }

        // Compute interval from 1-faces
        for (int row = 1; row < interpNumRows - 1; row++) {
            for (int col = 1; col < interpNumCols - 1; col++) {
                if (row % 2 == 1 && col % 2 == 1)
                    continue; // Already defined.

                pT = ImageUtils::to1D(row, col, interpNumCols);
                if (col % 2 == 0 && row % 2 == 1) {
                    int v1 = getValue(row, col + 1);
                    int v2 = getValue(row, col - 1);
                    interpolationMin[pT] = std::min(v1, v2);
                    interpolationMax[pT] = std::max(v1, v2);
                } else if (col % 2 == 1 && row % 2 == 0) {
                    int v1 = getValue(row + 1, col);
                    int v2 = getValue(row - 1, col);
                    interpolationMin[pT] = std::min(v1, v2);
                    interpolationMax[pT] = std::max(v1, v2);
                }
            }
        }
        // Compute interval from 0-faces
        for (int row = 1; row < interpNumRows - 1; row++) {
            for (int col = 1; col < interpNumCols - 1; col++) {
                if (row % 2 == 1 && col % 2 == 1)
                    continue; // Already defined.
                pT = ImageUtils::to1D(row, col, interpNumCols);
                if (row % 2 == 0 && col % 2 == 0) {
                    // | v0 | v1 |
                    // | v2 | v3 |
                    int v0 = getValue(row - 1, col - 1);
                    int v1 = getValue(row + 1, col - 1);
                    int v2 = getValue(row - 1, col + 1);
                    int v3 = getValue(row + 1, col + 1);

                    const int min_v0v3 = std::min(v0, v3);
                    const int max_v0v3 = std::max(v0, v3);
                    const int min_v1v2 = std::min(v1, v2);
                    const int max_v1v2 = std::max(v1, v2);
                    const bool diagonal0IsHigh = min_v0v3 > max_v1v2;
                    const bool diagonal1IsHigh = min_v1v2 > max_v0v3;

                    if (diagonal0IsHigh || diagonal1IsHigh) {
                        const bool chooseDiagonal0 = usesHighDiagonalAtSaddle() ? diagonal0IsHigh : !diagonal0IsHigh;
                        if (chooseDiagonal0) {
                            setDiagonal0Connection(adj, row, col);
                            interpolationMin[pT] = min_v0v3;
                            interpolationMax[pT] = max_v0v3;
                        } else {
                            setDiagonal1Connection(adj, row, col);
                            interpolationMin[pT] = min_v1v2;
                            interpolationMax[pT] = max_v1v2;
                        }
                    } else {
                        // Non-critical configuration.
                        interpolationMin[pT] = std::min({v0, v1, v2, v3});
                        interpolationMax[pT] = std::max({v0, v1, v2, v3});
                    }
                }
            }
        }
        return std::make_tuple(std::move(interpolationMin), std::move(interpolationMax), std::move(adj));
    }

  private:
    /**
     * @brief Floods image.
     *
     * @param imgPtr Image data represented by `imgPtr`.
     * @return Number of rows in the interpolated domain.
     */
    FloodResult floodImage(const ImageUInt8Ptr& imgPtr) const {
        if (!imgPtr) {
            throw std::invalid_argument("TreeOfShapesProducer requires a non-null image.");
        }

        int numRows = imgPtr->getNumRows();
        int numCols = imgPtr->getNumCols();

        int interpNumCols = interpolatedNumCols(numCols);
        int interpNumRows = interpolatedNumRows(numRows);
        int size = checkedDomainSize(interpNumRows, interpNumCols);
        const bool connectivityMap = usesConnectivityMap();
        auto [interpolationMin, interpolationMax, adj] = connectivityMap ? interpolateImage4c8c(imgPtr) : interpolateImage(imgPtr);

        std::vector<uint8_t> dejavu(size, 0);
        std::vector<int> imgR(size);               // Ordered pixels.
        std::vector<ToSFloodDepth> imgU(size);     // Monotone max-tree levels.
        std::vector<ToSGrayLevel> grayLevel(size); // Actual propagation levels.

        PriorityQueueToS queue(connectivityMap ? ToSUInt8Depth : ToSSelfDualDepth); // Priority queue.
        int infinityPixel = infinitySeedIndex(interpNumRows, interpNumCols);
        int priorityQueueOld = interpolationMin[infinityPixel];
        queue.initial(infinityPixel, priorityQueueOld);
        dejavu[infinityPixel] = true;

        int order = 0;
        ToSFloodDepth depth = 0;
        while (!queue.isEmpty()) {
            int h = queue.priorityPop();                    // Pop the element with highest priority.
            int priorityQueue = queue.getCurrentPriority(); // Current priority.
            if (connectivityMap) {
                if (priorityQueue != priorityQueueOld)
                    depth++;
                imgU[h] = depth;
            } else {
                imgU[h] = static_cast<ToSFloodDepth>(priorityQueue);
            }
            grayLevel[h] = static_cast<ToSGrayLevel>(priorityQueue);

            // Store h in the correct output order.
            imgR[order++] = h;

            // Adjacencies.
            for (int n : adj.getNeighborIndices(h)) {
                if (!dejavu[n]) {
                    queue.priorityPush(n, interpolationMin[n], interpolationMax[n]);
                    dejavu[n] = true; // Mark as processed.
                }
            }
            priorityQueueOld = priorityQueue;
        }
        if (order != size) {
            throw std::runtime_error("Tree-of-shapes propagation did not visit the complete interpolated domain.");
        }
        return FloodResult{std::move(imgU), std::move(grayLevel), std::move(imgR), std::move(adj)};
    }

    /**
     * @brief Builds canonical flood tree.
     *
     * @param imgPtr Image data represented by `imgPtr`.
     * @return Resulting canonical flood tree.
     */
    CanonicalFloodTree buildCanonicalFloodTree(const ImageUInt8Ptr& imgPtr) const {
        FloodResult flood = floodImage(imgPtr);
        const int numPixelsInterp = static_cast<int>(flood.treeLevel.size());

        std::vector<int> zPar(static_cast<size_t>(numPixelsInterp), InvalidNode);
        std::vector<int> parentInterpolate(static_cast<size_t>(numPixelsInterp), InvalidNode);
        auto findRoot = [&](int pStar) {
            while (zPar[pStar] != pStar) {
                zPar[pStar] = zPar[zPar[pStar]];
                pStar = zPar[pStar];
            }
            return pStar;
        };

        for (int i = numPixelsInterp - 1; i >= 0; --i) {
            const int pStar = flood.order[static_cast<size_t>(i)];
            parentInterpolate[pStar] = pStar;
            zPar[pStar] = pStar;
            for (int qStar : flood.adjacency.getNeighborIndices(pStar)) {
                if (zPar[qStar] != InvalidNode) {
                    const int rStar = findRoot(qStar);
                    if (pStar != rStar) {
                        parentInterpolate[rStar] = pStar;
                        zPar[rStar] = pStar;
                    }
                }
            }
        }

        auto sameLevel = [&](int aStar, int bStar) { return flood.treeLevel[aStar] == flood.treeLevel[bStar]; };
        for (int pStar : flood.order) {
            const int qStar = parentInterpolate[pStar];
            if (sameLevel(parentInterpolate[qStar], qStar)) {
                parentInterpolate[pStar] = parentInterpolate[qStar];
            }
        }

        return CanonicalFloodTree{std::move(flood.treeLevel), std::move(flood.grayLevel), std::move(flood.order), std::move(parentInterpolate)};
    }

  public:
    /**
     * @brief Sorts the interpolated domain for tree-of-shapes flooding.
     *
     * @param imgPtr Image data represented by `imgPtr`.
     * @return Values produced by the operation.
     */
    std::tuple<std::vector<ToSFloodDepth>, std::vector<int>, AdjacencyUC> sort(const ImageUInt8Ptr& imgPtr) const {
        FloodResult flood = floodImage(imgPtr);
        return std::make_tuple(std::move(flood.treeLevel), std::move(flood.order), std::move(flood.adjacency));
    }

    // Tests whether an interpolated pixel corresponds to an original pixel.
    /**
     * @brief Checks whether an interpolated pixel belongs to the original image grid.
     *
     * @param p Row-major pixel index to transform or inspect.
     * @param interpNumCols Number of columns in the interpolated domain.
     * @return True when the interpolated index represents an original image pixel.
     */
    inline bool isOriginal1D(int p, int interpNumCols) const {
        const int row = p / interpNumCols;
        const int col = p - row * interpNumCols;
        const int offset = interpolationOffset();
        return (row & 1) == offset && (col & 1) == offset;
    }

    // Maps an interpolated pixel to the original image domain.
    /**
     * @brief Maps an interpolated-domain pixel index to the original image domain.
     *
     * @param pStar Pixel index in the interpolated domain.
     * @param interNumCols Number of columns in the interpolated domain.
     * @param numCols Number of columns in the domain.
     * @return Padding offset in interpolated-domain coordinates.
     */
    inline int toOriginal1D(int pStar, int interNumCols, int numCols) const {
        int r = pStar / interNumCols;
        int c = pStar - r * interNumCols; // Avoid the modulo operator.
        const int offset = interpolationOffset();
        return ((r - offset) >> 1) * numCols + ((c - offset) >> 1);
    }
    /// @endcond

    /**
     * @brief Builds a native Tree-of-Shapes topology.
     *
     * The projection retains every interpolated plateau that owns an original
     * pixel. A plateau without original pixels is retained exactly when it
     * combines at least two distinct projected child supports. Consequently,
     * unary empty plateaus are contracted while a virtual branching root is
     * represented faithfully with an empty direct proper part.
     *
     * @param imgPtr Image used by the operation.
     * @return The resulting native Tree-of-Shapes topology.
     */
    template <ToSAltitudeEncodingPolicy EncodingPolicy>
    [[nodiscard]] TreeOfShapesBuildResultT<typename EncodingPolicy::value_type> buildWithAltitudeEncoding(const ImageUInt8Ptr& imgPtr,
                                                                                                          EncodingPolicy = {}) const {
        if (!imgPtr) {
            throw std::invalid_argument("TreeOfShapesProducer requires a non-null image.");
        }
        const int numRows = imgPtr->getNumRows();
        const int numCols = imgPtr->getNumCols();
        const int numPixels = checkedDomainSize(numRows, numCols);
        const int interpNumCols = interpolatedNumCols(numCols);

        CanonicalFloodTree canonical = buildCanonicalFloodTree(imgPtr);
        const int numPixelsInterp = static_cast<int>(canonical.parent.size());
        auto sameLevel = [&](int aStar, int bStar) {
            return canonical.treeLevel[static_cast<size_t>(aStar)] == canonical.treeLevel[static_cast<size_t>(bStar)];
        };
        auto repOf = [&](int pStar) {
            const int parentStar = canonical.parent[static_cast<size_t>(pStar)];
            return (parentStar == pStar || sameLevel(parentStar, pStar)) ? parentStar : pStar;
        };

        std::vector<int> canonicalReps;
        canonicalReps.reserve(static_cast<size_t>(numPixels));
        for (int pStar : canonical.order) {
            const int parentStar = canonical.parent[static_cast<size_t>(pStar)];
            if (parentStar == pStar || !sameLevel(parentStar, pStar)) {
                canonicalReps.push_back(pStar);
            }
        }
        if (canonicalReps.empty()) {
            throw std::runtime_error("Tree-of-shapes construction produced no canonical interpolated node.");
        }

        const int rawRoot = canonicalReps.front();
        std::vector<int> rawParent(static_cast<size_t>(numPixelsInterp), InvalidNode);
        std::vector<int> firstRawChild(static_cast<size_t>(numPixelsInterp), InvalidNode);
        std::vector<int> nextRawSibling(static_cast<size_t>(numPixelsInterp), InvalidNode);
        for (int repStar : canonicalReps) {
            if (repStar == rawRoot) {
                rawParent[static_cast<size_t>(repStar)] = repStar;
                continue;
            }
            const int parentRep = repOf(canonical.parent[static_cast<size_t>(repStar)]);
            if (parentRep < 0 || parentRep >= numPixelsInterp || rawParent[static_cast<size_t>(parentRep)] == InvalidNode) {
                throw std::runtime_error("Tree-of-shapes canonical node has no canonical parent.");
            }
            rawParent[static_cast<size_t>(repStar)] = parentRep;
            nextRawSibling[static_cast<size_t>(repStar)] = firstRawChild[static_cast<size_t>(parentRep)];
            firstRawChild[static_cast<size_t>(parentRep)] = repStar;
        }

        std::vector<int> directCount(static_cast<size_t>(numPixelsInterp), 0);
        std::vector<int> ownerRepByElement(static_cast<size_t>(numPixels), InvalidNode);
        for (int pStar = 0; pStar < numPixelsInterp; ++pStar) {
            if (!isOriginal1D(pStar, interpNumCols)) {
                continue;
            }
            const int elementId = toOriginal1D(pStar, interpNumCols, numCols);
            if (elementId < 0 || elementId >= numPixels) {
                throw std::runtime_error("Tree-of-shapes projection produced an element outside the source domain.");
            }
            const int ownerRep = repOf(pStar);
            ownerRepByElement[static_cast<size_t>(elementId)] = ownerRep;
            ++directCount[static_cast<size_t>(ownerRep)];
        }
        if (std::find(ownerRepByElement.begin(), ownerRepByElement.end(), InvalidNode) != ownerRepByElement.end()) {
            throw std::runtime_error("Tree-of-shapes projection did not assign every source-domain element.");
        }

        // The canonical depth, parent, and full interpolated order are no
        // longer needed after the raw hierarchy and source-domain owners have
        // been extracted. Release them before allocating the projection buffer.
        std::vector<ToSFloodDepth>().swap(canonical.treeLevel);
        std::vector<int>().swap(canonical.parent);
        std::vector<int>().swap(canonical.order);

        // Bottom-up support projection. projectedRep[r] is the canonical
        // representative of the distinct projected support rooted at r.
        std::vector<int> projectedRep(static_cast<size_t>(numPixelsInterp), InvalidNode);
        for (auto it = canonicalReps.rbegin(); it != canonicalReps.rend(); ++it) {
            const int repStar = *it;
            int numProjectedChildren = 0;
            int onlyProjectedChild = InvalidNode;
            for (int childStar = firstRawChild[static_cast<size_t>(repStar)]; childStar != InvalidNode;
                 childStar = nextRawSibling[static_cast<size_t>(childStar)]) {
                const int childProjection = projectedRep[static_cast<size_t>(childStar)];
                if (childProjection != InvalidNode) {
                    ++numProjectedChildren;
                    onlyProjectedChild = childProjection;
                }
            }

            if (directCount[static_cast<size_t>(repStar)] > 0 || numProjectedChildren >= 2) {
                projectedRep[static_cast<size_t>(repStar)] = repStar;
            } else if (numProjectedChildren == 1) {
                projectedRep[static_cast<size_t>(repStar)] = onlyProjectedChild;
            }
        }

        const int projectedRootRep = projectedRep[static_cast<size_t>(rawRoot)];
        if (projectedRootRep == InvalidNode || projectedRep[static_cast<size_t>(projectedRootRep)] != projectedRootRep) {
            throw std::runtime_error("Tree-of-shapes projection produced an empty hierarchy.");
        }

        std::vector<int>().swap(firstRawChild);
        std::vector<int>().swap(nextRawSibling);
        std::vector<int>().swap(directCount);

        std::vector<NodeId> nodeIdByRep(static_cast<size_t>(numPixelsInterp), InvalidNode);
        NodeId numNodes = 0;
        for (int repStar : canonicalReps) {
            if (projectedRep[static_cast<size_t>(repStar)] == repStar) {
                nodeIdByRep[static_cast<size_t>(repStar)] = numNodes++;
            }
        }

        using OutputLevel = typename EncodingPolicy::value_type;
        TreeOfShapesBuildResultT<OutputLevel> result;
        result.nodeParent.assign(static_cast<size_t>(numNodes), InvalidNode);
        result.properPartOwner.assign(static_cast<size_t>(numPixels), InvalidNode);
        result.altitude.assign(static_cast<size_t>(numNodes), OutputLevel{});
        result.root = nodeIdByRep[static_cast<size_t>(projectedRootRep)];
        result.numRows = numRows;
        result.numCols = numCols;
        detail::TopologicalNativeHierarchyRecorder proofRecorder(static_cast<std::size_t>(numNodes), static_cast<std::size_t>(numPixels), result.root);

        for (int repStar : canonicalReps) {
            const NodeId nodeId = nodeIdByRep[static_cast<size_t>(repStar)];
            if (nodeId == InvalidNode) {
                continue;
            }

            if (repStar == projectedRootRep) {
                result.nodeParent[static_cast<size_t>(nodeId)] = nodeId;
            } else {
                int ancestorRep = rawParent[static_cast<size_t>(repStar)];
                while (ancestorRep != InvalidNode && nodeIdByRep[static_cast<size_t>(ancestorRep)] == InvalidNode) {
                    if (ancestorRep == rawParent[static_cast<size_t>(ancestorRep)]) {
                        ancestorRep = InvalidNode;
                        break;
                    }
                    ancestorRep = rawParent[static_cast<size_t>(ancestorRep)];
                }
                if (ancestorRep == InvalidNode) {
                    throw std::runtime_error("Tree-of-shapes retained node has no retained parent.");
                }
                result.nodeParent[static_cast<size_t>(nodeId)] = nodeIdByRep[static_cast<size_t>(ancestorRep)];
            }
            proofRecorder.recordSupportedNode(nodeId, result.nodeParent[static_cast<std::size_t>(nodeId)]);

            result.altitude[static_cast<size_t>(nodeId)] = EncodingPolicy::encode(canonical.grayLevel[static_cast<size_t>(repStar)], options_.interpolation);
        }

        for (int elementId = 0; elementId < numPixels; ++elementId) {
            const int ownerRep = ownerRepByElement[static_cast<size_t>(elementId)];
            const NodeId ownerNode = nodeIdByRep[static_cast<size_t>(ownerRep)];
            if (ownerNode == InvalidNode) {
                throw std::runtime_error("Tree-of-shapes proper part belongs to a contracted node.");
            }
            result.properPartOwner[static_cast<size_t>(elementId)] = ownerNode;
            proofRecorder.recordProperPart(elementId, ownerNode);
        }
        result.topologyProof = std::move(proofRecorder).finish();
        return result;
    }

    /**
     * @brief Builds the quantized hierarchy in source uint8 gray units.
     *
     * Half-level structural nodes are quantized downward. Reconstruction is
     * exact because those virtual projected nodes own no source-domain
     * element.
     *
     * @param imgPtr Image used by the operation.
     * @return The resulting quantized hierarchy in source uint8 gray units.
     */
    [[nodiscard]] TreeOfShapesBuildResult build(const ImageUInt8Ptr& imgPtr) const {
        return buildWithAltitudeEncoding(imgPtr, ToSQuantizedUInt8AltitudeEncoding{});
    }

    /**
     * @brief Builds the hierarchy with exact doubled source-level altitudes.
     *
     * Even values represent source gray levels and odd values represent the
     * half levels introduced by self-dual interpolation.
     *
     * @param imgPtr Image used by the operation.
     * @return The resulting hierarchy with exact doubled source-level altitudes.
     */
    [[nodiscard]] ExactTreeOfShapesBuildResult buildExact(const ImageUInt8Ptr& imgPtr) const {
        return buildWithAltitudeEncoding(imgPtr, ToSExactDoubledAltitudeEncoding{});
    }
};

} // namespace mmcfilters
