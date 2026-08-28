#pragma once

#include "../utils/Image.hpp"
#include "../utils/Common.hpp"
#include "MorphologicalTreeSemantics.hpp"
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

namespace detail {
/** @brief Internal dispatch selected from a validated topographic immersion. */
enum class TopographicImmersionMode { SelfDualSpan, Min4Max8, Min8Max4 };
} // namespace detail

inline constexpr int TopographicInterpolationScale = 2;
inline constexpr int TopographicInterpolationPadding = 1;
inline constexpr int ToSUInt8Depth = 8;
inline constexpr int ToSSelfDualDepth = 9;

using ToSGrayLevel = uint16_t;
using ToSFloodDepth = uint32_t;

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
     * @param constructionLevel Altitude or level.
     * @param immersionMode Active topographic immersion mode.
     * @return The converted one transient construction level to exact doubled units.
     */
    static value_type encode(ToSGrayLevel constructionLevel, detail::TopographicImmersionMode immersionMode) noexcept {
        if (immersionMode == detail::TopographicImmersionMode::SelfDualSpan) {
            return constructionLevel;
        }
        return static_cast<value_type>(TopographicInterpolationScale * static_cast<unsigned int>(constructionLevel));
    }
};

/**
 * @brief Exact encoding in 8-bit source-level units.
 *
 * A complementary-grid immersion floods the interpolated domain over the source
 * level set itself, so its construction levels are already source gray levels.
 * A self-dual span immersion floods in doubled units, so its construction
 * levels are halved; they are guaranteed to be even only when no half level can
 * enter the active domain, which the convention validation establishes before
 * this encoding is selected.
 */
struct ToSUInt8AltitudeEncoding {
    /// Scalar type storing source gray levels.
    using value_type = std::uint8_t;

    /**
     * @brief Converts one transient construction level to 8-bit source units.
     *
     * @param constructionLevel Altitude or level.
     * @param immersionMode Active topographic immersion mode.
     * @return The converted one transient construction level to 8-bit source units.
     */
    static value_type encode(ToSGrayLevel constructionLevel, detail::TopographicImmersionMode immersionMode) {
        unsigned int sourceLevel = constructionLevel;
        if (immersionMode == detail::TopographicImmersionMode::SelfDualSpan) {
            if (constructionLevel % TopographicInterpolationScale != 0) {
                throw std::logic_error("Tree-of-shapes 8-bit altitudes cannot represent a self-dual half level.");
            }
            sourceLevel = constructionLevel / TopographicInterpolationScale;
        }
        if (sourceLevel > std::numeric_limits<value_type>::max()) {
            throw std::logic_error("Tree-of-shapes construction level exceeds the 8-bit source level set.");
        }
        return static_cast<value_type>(sourceLevel);
    }
};

/**
 * @brief Selects the altitude encoding policy declared by a topographic convention.
 *
 * @tparam Altitude Published node-altitude type.
 */
template <class Altitude> struct ToSAltitudeEncodingFor;

/** @brief Maps the 8-bit published type to its exact source-level encoding. */
template <> struct ToSAltitudeEncodingFor<std::uint8_t> {
    using type = ToSUInt8AltitudeEncoding; ///< Selected encoding policy.
    static constexpr TopographicAltitudeEncoding declaration = TopographicAltitudeEncoding::UInt8; ///< Matching declaration.
};

/** @brief Maps the doubled published type to its exact doubled-unit encoding. */
template <> struct ToSAltitudeEncodingFor<ToSGrayLevel> {
    using type = ToSExactDoubledAltitudeEncoding; ///< Selected encoding policy.
    static constexpr TopographicAltitudeEncoding declaration = TopographicAltitudeEncoding::ExactDoubled; ///< Matching declaration.
};

/**
 * @brief Native representation produced by a Tree-of-Shapes producer.
 *
 * Unlike a pixel-parent image, this representation can preserve a virtual root
 * (or another branching projected node) whose direct proper part is empty.
 * Node altitudes are published in the scale declared by the topographic
 * convention: unchanged 8-bit source levels, or exact doubled gray-level units
 * where source levels are even and odd values preserve self-dual half levels.
 *
 * @tparam Altitude Published node-altitude type.
 */
template <class Altitude> struct TreeOfShapesBuildResult {
    /// Parent id for every produced internal node.
    std::vector<NodeId> parent;
    /// Direct owning node for every original-domain proper part.
    std::vector<NodeId> smallestNodeMap;
    /// Encoded altitude for every produced internal node.
    std::vector<Altitude> nodeAltitudes;
    /// Root node id in the produced internal-node domain.
    NodeId root = InvalidNode;
    /// Number of rows in the original pixel domain.
    int numRows = 0;
    /// Number of columns in the original pixel domain.
    int numColumns = 0;
    /// Move-only proof that the producer established the topology invariants.
    detail::NativeTopologyProof topologyProof;

    /**
     * @brief Transfers the producer-owned buffers together with their generic
     * topology proof.
     *
     * @param semantics Hierarchy semantics validated by the operation.
     * @return The transferred producer-owned buffers together with their generic topology proof.
     */
    [[nodiscard]] detail::ValidatedNativeHierarchy<Altitude> takeValidatedHierarchy(MorphologicalTreeSemantics semantics) && {
        return detail::makeValidatedNativeHierarchy<Altitude>(std::move(parent), std::move(smallestNodeMap), std::move(nodeAltitudes), root,
                                                              GridDomain2D{numRows, numColumns}, std::move(semantics), std::move(topologyProof));
    }
};

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
enum class DiagonalConnection : uint8_t { None = 0, Sw = 1 << 0, Ne = 1 << 1, Se = 1 << 2, Nw = 1 << 3 };

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
    /** @brief Number of rows in the image domain. */
    int numRows;
    /** @brief Number of columns in the image domain. */
    int numColumns;
    /** @brief Enabled diagonal-connection flags for each pixel buffer. */
    std::vector<uint8_t> dconnFlags; // 4-connect.  +  diag. connect.
                                     //  N, W, S, E,   SW, NE, SE, NW
    /** @brief Stores row offsets for cardinal and diagonal neighbours. */
    const std::vector<int> offsetRows = {-1, 0, 1, 0, 1, -1, 1, -1};
    /** @brief Stores column offsets for cardinal and diagonal neighbours. */
    const std::vector<int> offsetColumns = {0, -1, 0, 1, -1, 1, 1, -1};
    /** @brief Indicates whether diagonal connections may be traversed. */
    bool enableDiagonalConnection;
    /** @brief Maps each diagonal offset to its required connection flag. */
    const std::vector<DiagonalConnection> requiredDiagonal = {DiagonalConnection::Sw, DiagonalConnection::Ne, DiagonalConnection::Se, DiagonalConnection::Nw};

  public:
    /**
     * @brief Constructs `AdjacencyUC` from the supplied inputs.
     *
     * @param rows Number of rows in the domain.
     * @param columns Number of columns in the domain.
     * @param enableDiagonalConnection Flag controlling enable diagonal connection.
     */
    AdjacencyUC(int rows, int columns, bool enableDiagonalConnection) : numRows(rows), numColumns(columns), enableDiagonalConnection(enableDiagonalConnection) {
        if (enableDiagonalConnection)
            dconnFlags.resize(rows * columns, 0);
    }

    /**
     * @brief Destroys `AdjacencyUC`.
     */
    ~AdjacencyUC() {}

    /**
     * @brief Sets diagonal connection.
     *
     * @param row Zero-based row coordinate.
     * @param column Zero-based column coordinate.
     * @param conn Diagonal-connection flag.
     */
    void setDiagonalConnection(int row, int column, DiagonalConnection conn) { dconnFlags[ImageUtils::to1D(row, column, numColumns)] |= static_cast<uint8_t>(conn); }

    /**
     * @brief Sets diagonal connection.
     *
     * @param pixel Row-major pixel index whose diagonal connection is updated.
     * @param conn Diagonal-connection flag.
     */
    void setDiagonalConnection(PixelId pixel, DiagonalConnection conn) { dconnFlags[pixel] |= static_cast<uint8_t>(conn); }

    /**
     * @brief Tests whether connection holds.
     *
     * @param row Zero-based row coordinate.
     * @param column Zero-based column coordinate.
     * @param conn Diagonal-connection flag.
     * @return True when connection; otherwise false.
     */
    bool hasConnection(int row, int column, DiagonalConnection conn) const { return dconnFlags[ImageUtils::to1D(row, column, numColumns)] & static_cast<uint8_t>(conn); }

    /**
     * @brief Returns connections.
     *
     * @param row Zero-based row coordinate.
     * @param column Zero-based column coordinate.
     * @return Connections.
     */
    uint8_t getConnections(int row, int column) const { return dconnFlags[ImageUtils::to1D(row, column, numColumns)]; }

    /**
     * @brief Iterator over valid neighbours including enabled diagonal links.
     */
    class NeighborIterator {
      private:
        /** @brief References the adaptive adjacency being traversed. */
        AdjacencyUC& instance;
        /** @brief Source pixel row. */
        int row;
        /** @brief Source pixel column. */
        int column;
        /** @brief Current neighbour-offset index. */
        std::size_t id;

        /**
         * @brief Advances to valid.
         */
        void advanceToValid() {
            while (id < instance.offsetRows.size()) {
                int r = row + instance.offsetRows[id];
                int c = column + instance.offsetColumns[id];
                if (r >= 0 && c >= 0 && r < instance.numRows && c < instance.numColumns) {
                    if (id < 4 || (instance.enableDiagonalConnection && instance.dconnFlags[ImageUtils::to1D(row, column, instance.numColumns)] &
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
         * @param column Zero-based column coordinate.
         * @param id Identifier.
         */
        NeighborIterator(AdjacencyUC& adj, int row, int column, int id) : instance(adj), row(row), column(column), id(id) { advanceToValid(); }

        /**
         * @brief Returns the value at the current iterator position.
         *
         * @return The value at the current iterator position.
         */
        PixelId operator*() const {
            int dr = instance.offsetRows[id];
            int dc = instance.offsetColumns[id];
            return ImageUtils::to1D(row + dr, column + dc, instance.numColumns);
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
        /** @brief Source pixel row. */
        int row;
        /** @brief Source pixel column. */
        int column;

      public:
        /**
         * @brief Constructs `NeighborRange` from the supplied inputs.
         *
         * @param instance Producer instance whose priority state is queried.
         * @param row Zero-based row coordinate.
         * @param column Zero-based column coordinate.
         */
        NeighborRange(AdjacencyUC& instance, int row, int column) : instance(instance), row(row), column(column) {}

        /**
         * @brief Returns an iterator positioned at the beginning.
         *
         * @return An iterator positioned at the beginning.
         */
        NeighborIterator begin() { return NeighborIterator(instance, row, column, 0); }
        /**
         * @brief Returns the exclusive-end iterator.
         *
         * @return The exclusive-end iterator.
         */
        NeighborIterator end() { return NeighborIterator(instance, row, column, 8); }
    };

    /**
     * @brief Returns neighbor indices.
     *
     * @param pixel Row-major pixel index to transform or inspect.
     * @return Neighbor indices.
     */
    NeighborRange getNeighborIndices(PixelId pixel) {
        auto [row, column] = ImageUtils::to2D(pixel, numColumns);
        return NeighborRange(*this, row, column);
    }

    /**
     * @brief Returns neighbor indices.
     *
     * @param row Zero-based row coordinate.
     * @param column Zero-based column coordinate.
     * @return Neighbor indices.
     */
    NeighborRange getNeighborIndices(int row, int column) { return NeighborRange(*this, row, column); }
};

/**
 * @brief Discrete priority queue used during tree-of-shapes construction.
 */
class PriorityQueueToS {
  private:
    /** @brief Stores queued pixels grouped by discrete priority. */
    std::vector<std::deque<PixelId>> buckets;
    /** @brief Priority of the next non-empty bucket. */
    int currentPriority;
    /** @brief Total number of queued elements. */
    int numElements;
    /** @brief Number of allocated priority levels. */
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
    void initial(PixelId element, int priority) {
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
    void priorityPush(PixelId element, int lower, int upper) {
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
    PixelId priorityPop() {
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

        const PixelId element = buckets[currentPriority].front();
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
 * The builder interpolates the input image according to the retained
 * `TopographicConvention`,
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

    /** @brief Intermediate topology produced by the tree-of-shapes flood. */
    struct FloodResult {
        /** @brief Tree level buffer. */
        std::vector<ToSFloodDepth> treeLevel;
        /** @brief Gray level buffer. */
        std::vector<ToSGrayLevel> grayLevel;
        /** @brief Pixel identifier of the order. */
        std::vector<PixelId> order;
        /** @brief Adjacency. */
        AdjacencyUC adjacency;
    };

    /** @brief Canonical flood-tree topology and interpolation metadata. */
    struct CanonicalFloodTree {
        /** @brief Tree level buffer. */
        std::vector<ToSFloodDepth> treeLevel;
        /** @brief Gray level buffer. */
        std::vector<ToSGrayLevel> grayLevel;
        /** @brief Pixel identifier of the order. */
        std::vector<PixelId> order;
        /** @brief Pixel identifier of the parent. */
        std::vector<PixelId> parent;
    };

    /** @brief Complete topographic convention. */
    TopographicConvention convention_;
    /** @brief Validated internal immersion dispatch. */
    detail::TopographicImmersionMode immersionMode_ = detail::TopographicImmersionMode::SelfDualSpan;

    /**
     * @brief Validates and classifies a topographic convention.
     *
     * @param convention Convention controlling the operation.
     * @return Internal immersion mode corresponding to the convention.
     */
    static detail::TopographicImmersionMode validateConvention(const TopographicConvention& convention) {
        if (convention.infinityPixel < 0) {
            throw std::invalid_argument("Tree-of-shapes infinity pixel must be non-negative.");
        }
        if (std::holds_alternative<SelfDualSpanImmersion>(convention.immersion)) {
            // Only the exterior ring carries the boundary reference level, which
            // is the mean of the two central boundary values on an even boundary
            // and therefore the single source of half levels. Without that ring
            // the reference level is cropped away and never read by an interior
            // cell, so every construction level stays on the source lattice.
            if (convention.altitudeEncoding == TopographicAltitudeEncoding::UInt8 &&
                convention.domainExtension == TopographicDomainExtension::ExteriorRing) {
                throw std::invalid_argument("Self-dual span immersion over an exterior ring cannot publish 8-bit altitudes because its boundary reference "
                                            "level may fall on a half level; declare TopographicAltitudeEncoding::ExactDoubled or "
                                            "TopographicDomainExtension::None.");
            }
            return detail::TopographicImmersionMode::SelfDualSpan;
        }
        if (const auto* canonical = std::get_if<CanonicalComplementaryGridImmersion>(&convention.immersion)) {
            return canonical->pairing == ComplementaryPairing::Min4Max8 ? detail::TopographicImmersionMode::Min4Max8
                                                                       : detail::TopographicImmersionMode::Min8Max4;
        }
        const auto& adjacencies = std::get<ComplementaryGridImmersion>(convention.immersion).complementaryAdjacencies;
        if (adjacencies.minAdjacency.is4connectivity() && adjacencies.maxAdjacency.is8connectivity()) {
            return detail::TopographicImmersionMode::Min4Max8;
        }
        if (adjacencies.minAdjacency.is8connectivity() && adjacencies.maxAdjacency.is4connectivity()) {
            return detail::TopographicImmersionMode::Min8Max4;
        }
        throw std::invalid_argument("Complementary-grid immersion requires minimum/maximum adjacencies in the canonical 4/8 or 8/4 pairing.");
    }

    /** @brief Validates convention data whose domain is known only when an image is supplied. @param image Source image defining the domain. */
    void validateConventionForImage(const ImageUInt8Ptr& image) const {
        if (!image) {
            throw std::invalid_argument("TreeOfShapesProducer requires a non-null image.");
        }
        if (const auto* immersion = std::get_if<ComplementaryGridImmersion>(&convention_.immersion)) {
            const auto& adjacencies = immersion->complementaryAdjacencies;
            const int rows = image->getNumRows();
            const int columns = image->getNumColumns();
            if (adjacencies.minAdjacency.getNumRows() != rows || adjacencies.minAdjacency.getNumColumns() != columns ||
                adjacencies.maxAdjacency.getNumRows() != rows || adjacencies.maxAdjacency.getNumColumns() != columns) {
                throw std::invalid_argument("Tree-of-shapes complementary adjacencies must match the source image domain.");
            }
        }
        static_cast<void>(validatedInfinityPixel(interpolatedNumRows(image->getNumRows()), interpolatedNumColumns(image->getNumColumns())));
    }

    /**
     * @brief Tests whether connectivity map holds.
     *
     * @return True when connectivity map; otherwise false.
     */
    inline bool usesConnectivityMap() const noexcept {
        return immersionMode_ != detail::TopographicImmersionMode::SelfDualSpan;
    }

    /**
     * @brief Tests whether high diagonal at saddle holds.
     *
     * @return True when high diagonal at saddle; otherwise false.
     */
    inline bool usesHighDiagonalAtSaddle() const noexcept { return immersionMode_ == detail::TopographicImmersionMode::Min4Max8; }

    /**
     * @brief Tests whether exterior padding holds.
     *
     * @return True when exterior padding; otherwise false.
     */
    inline bool usesExteriorPadding() const noexcept { return convention_.domainExtension == TopographicDomainExtension::ExteriorRing; }

    /**
     * @brief Returns the padding offset applied to the interpolated domain.
     *
     * @return Padding offset in interpolated-domain coordinates.
     */
    inline int interpolationOffset() const noexcept { return usesExteriorPadding() ? TopographicInterpolationPadding : 0; }

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
        const std::int64_t interpolated = TopographicInterpolationScale * static_cast<std::int64_t>(extent) + (usesExteriorPadding() ? 1 : -1);
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
     * @param numColumns Number of columns in the domain.
     * @return Validated extent of the interpolated dimension.
     */
    inline int interpolatedNumColumns(int numColumns) const { return checkedInterpolatedExtent(numColumns); }

    /**
     * @brief Checks and converts domain size.
     *
     * @param numRows Number of rows in the domain.
     * @param numColumns Number of columns in the domain.
     * @return Validated number of pixels in the original domain.
     */
    inline int checkedDomainSize(int numRows, int numColumns) const {
        if (numRows <= 0 || numColumns <= 0 || numRows > std::numeric_limits<int>::max() / numColumns) {
            throw std::overflow_error("Tree-of-shapes domain size exceeds int range.");
        }
        return numRows * numColumns;
    }

    /**
     * @brief Maps point row.
     *
     * @param row Zero-based row coordinate.
     * @return Mapped point row.
     */
    inline int originalPointRow(int row) const noexcept { return TopographicInterpolationScale * row + interpolationOffset(); }

    /**
     * @brief Maps point column.
     *
     * @param column Zero-based column coordinate.
     * @return Mapped point column.
     */
    inline int originalPointColumn(int column) const noexcept { return TopographicInterpolationScale * column + interpolationOffset(); }

    /**
     * @brief Scales original level.
     *
     * @param value Value.
     * @return Original gray level expressed in the interpolation scale.
     */
    inline ToSGrayLevel scaledOriginalLevel(uint8_t value) const noexcept {
        return static_cast<ToSGrayLevel>(TopographicInterpolationScale * static_cast<int>(value));
    }

    /**
     * @brief Validates and returns the configured infinity pixel in the interpolated domain.
     *
     * @param interpNumRows Number of rows in the interpolated domain.
     * @param interpNumColumns Number of columns in the interpolated domain.
     * @return Validated row-major infinity pixel.
     */
    inline PixelId validatedInfinityPixel(int interpNumRows, int interpNumColumns) const {
        const std::int64_t activeSize = static_cast<std::int64_t>(interpNumRows) * static_cast<std::int64_t>(interpNumColumns);
        if (convention_.infinityPixel < 0 || static_cast<std::int64_t>(convention_.infinityPixel) >= activeSize) {
            throw std::invalid_argument("Tree-of-shapes infinity pixel must belong to the active topographic domain.");
        }
        return convention_.infinityPixel;
    }

    /**
     * @brief Sets diagonal0 connection.
     *
     * @param adj Adjacency relation used to traverse the image domain.
     * @param row Zero-based row coordinate.
     * @param column Zero-based column coordinate.
     */
    inline void setDiagonal0Connection(AdjacencyUC& adj, int row, int column) const {
        adj.setDiagonalConnection(row, column - 1, DiagonalConnection::Se);
        adj.setDiagonalConnection(row + 1, column, DiagonalConnection::Nw);

        adj.setDiagonalConnection(row - 1, column - 1, DiagonalConnection::Se);
        adj.setDiagonalConnection(row, column, DiagonalConnection::Se | DiagonalConnection::Nw);
        adj.setDiagonalConnection(row + 1, column + 1, DiagonalConnection::Nw);

        adj.setDiagonalConnection(row - 1, column, DiagonalConnection::Se);
        adj.setDiagonalConnection(row, column + 1, DiagonalConnection::Nw);
    }

    /**
     * @brief Sets diagonal1 connection.
     *
     * @param adj Adjacency relation used to traverse the image domain.
     * @param row Zero-based row coordinate.
     * @param column Zero-based column coordinate.
     */
    inline void setDiagonal1Connection(AdjacencyUC& adj, int row, int column) const {
        adj.setDiagonalConnection(row, column - 1, DiagonalConnection::Ne);
        adj.setDiagonalConnection(row - 1, column, DiagonalConnection::Sw);

        adj.setDiagonalConnection(row - 1, column + 1, DiagonalConnection::Sw);
        adj.setDiagonalConnection(row, column, DiagonalConnection::Sw | DiagonalConnection::Ne);
        adj.setDiagonalConnection(row + 1, column - 1, DiagonalConnection::Ne);

        adj.setDiagonalConnection(row + 1, column, DiagonalConnection::Ne);
        adj.setDiagonalConnection(row, column + 1, DiagonalConnection::Sw);
    }

    std::tuple<std::vector<ToSGrayLevel>, std::vector<ToSGrayLevel>, AdjacencyUC>
    /**
     * @brief Crops exterior interpolation.
     *
     * @param paddedMin Minimum values in the padded interpolation domain.
     * @param paddedMax Maximum values in the padded interpolation domain.
     * @param paddedAdjacency Adjacency relation defined on the padded domain.
     * @param paddedRows Number of rows in the padded domain.
     * @param paddedColumns Number of columns in the padded domain.
     * @param adaptiveDiagonal Whether diagonal connections are chosen adaptively.
     * @return Values produced by the operation.
     */
    cropExteriorInterpolation(std::vector<ToSGrayLevel> paddedMin, std::vector<ToSGrayLevel> paddedMax, AdjacencyUC paddedAdjacency, int paddedRows,
                              int paddedColumns, bool adaptiveDiagonal) const {
        if (paddedRows < 3 || paddedColumns < 3) {
            throw std::logic_error("Tree-of-shapes padded immersion cannot be cropped.");
        }
        const int rows = paddedRows - 2;
        const int columns = paddedColumns - 2;
        const int size = checkedDomainSize(rows, columns);
        std::vector<ToSGrayLevel> croppedMin(static_cast<std::size_t>(size));
        std::vector<ToSGrayLevel> croppedMax(static_cast<std::size_t>(size));
        AdjacencyUC croppedAdjacency(rows, columns, adaptiveDiagonal);

        for (int row = 0; row < rows; ++row) {
            for (int column = 0; column < columns; ++column) {
                const PixelId source = ImageUtils::to1D(row + 1, column + 1, paddedColumns);
                const PixelId target = ImageUtils::to1D(row, column, columns);
                croppedMin[static_cast<std::size_t>(target)] = paddedMin[static_cast<std::size_t>(source)];
                croppedMax[static_cast<std::size_t>(target)] = paddedMax[static_cast<std::size_t>(source)];
                if (adaptiveDiagonal) {
                    const std::uint8_t connections = paddedAdjacency.getConnections(row + 1, column + 1);
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
     * @param convention Complete topographic convention.
     */
    explicit TreeOfShapesProducer(TopographicConvention convention = {})
        : convention_(std::move(convention)), immersionMode_(validateConvention(convention_)) {}

    /**
     * @brief Returns the immutable topographic convention selected for this producer.
     *
     * @return The immutable topographic convention selected for this producer.
     */
    [[nodiscard]] const TopographicConvention& convention() const noexcept { return convention_; }

    /**
     * @brief Destroys `TreeOfShapesProducer`.
     */
    ~TreeOfShapesProducer() = default;

    /// @cond INTERNAL
    /**
     * @brief Interpolates image.
     *
     * @param imgPtr Image data.
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
            TopographicConvention paddedConvention = convention_;
            paddedConvention.domainExtension = TopographicDomainExtension::ExteriorRing;
            paddedConvention.infinityPixel = 0;
            // The padded producer only supplies interpolation buffers, which are
            // altitude-encoding independent. Declaring the doubled encoding keeps
            // it constructible for a convention whose own encoding is 8-bit.
            paddedConvention.altitudeEncoding = TopographicAltitudeEncoding::ExactDoubled;
            TreeOfShapesProducer paddedProducer(std::move(paddedConvention));
            auto [paddedMin, paddedMax, paddedAdjacency] = paddedProducer.interpolateImage(imgPtr);
            const int paddedRows = paddedProducer.interpolatedNumRows(imgPtr->getNumRows());
            const int paddedColumns = paddedProducer.interpolatedNumColumns(imgPtr->getNumColumns());
            return cropExteriorInterpolation(std::move(paddedMin), std::move(paddedMax), std::move(paddedAdjacency), paddedRows, paddedColumns, false);
        }
        auto img = imgPtr->rawData();
        int numRows = imgPtr->getNumRows();
        int numColumns = imgPtr->getNumColumns();
        if (numRows <= 0 || numColumns <= 0) {
            throw std::invalid_argument("TreeOfShapesProducer requires a non-empty image.");
        }
        if (numRows == 1 && numColumns == 1) {
            const int interpNumRows = interpolatedNumRows(numRows);
            const int interpNumColumns = interpolatedNumColumns(numColumns);
            const int interpSize = checkedDomainSize(interpNumRows, interpNumColumns);
            std::vector<ToSGrayLevel> interpolationMin(static_cast<std::size_t>(interpSize), scaledOriginalLevel(img[0]));
            std::vector<ToSGrayLevel> interpolationMax(static_cast<std::size_t>(interpSize), scaledOriginalLevel(img[0]));
            AdjacencyUC adj(interpNumRows, interpNumColumns, false);
            return std::make_tuple(std::move(interpolationMin), std::move(interpolationMax), std::move(adj));
        }

        constexpr int adjCircleColumn[] = {-1, +1, -1, +1};
        constexpr int adjCircleRow[] = {-1, -1, +1, +1};

        constexpr int adjRetHorColumn[] = {0, 0};
        constexpr int adjRetHorRow[] = {-1, +1};

        constexpr int adjRetVerColumn[] = {+1, -1};
        constexpr int adjRetVerRow[] = {0, 0};

        int interpNumColumns = interpolatedNumColumns(numColumns);
        int interpNumRows = interpolatedNumRows(numRows);
        int size = checkedDomainSize(interpNumRows, interpNumColumns);

        // Allocate interpolation result buffers for minimum and maximum levels.
        std::vector<ToSGrayLevel> interpolationMin(size);
        std::vector<ToSGrayLevel> interpolationMax(size);

        // Collect each boundary pixel exactly once. The closed-form perimeter
        // size 2 * (rows + columns) - 4 is valid only when both dimensions are at
        // least two and previously left zero-filled entries for thin images.
        std::array<int, 256> boundaryHistogram{};
        int numBoundary = 0;

        PixelId pT;

        for (PixelId pixel = 0; pixel < numColumns * numRows; ++pixel) {
            auto [row, column] = ImageUtils::to2D(pixel, numColumns);

            // Check whether the pixel lies on the image border.
            if (row == 0 || row == numRows - 1 || column == 0 || column == numColumns - 1) {
                ++boundaryHistogram[static_cast<std::size_t>(img[pixel])];
                ++numBoundary;
            }

            // Compute the interpolated-image index.
            pT = ImageUtils::to1D(originalPointRow(row), originalPointColumn(column), interpNumColumns);

            // Assign interpolation values.
            interpolationMin[pT] = interpolationMax[pT] = scaledOriginalLevel(img[pixel]);
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
            median = TopographicInterpolationScale * boundaryValueAtRank(numBoundary / 2);
        }
        // std::cout << "Interpolation (Median): " << median << std::endl;

        PixelId qT;
        int qColumn, qRow, min, max;
        const int* adjColumn = nullptr;
        const int* adjRow = nullptr;
        int adjSize;
        AdjacencyUC adj(interpNumRows, interpNumColumns, false);

        for (int row = 0; row < interpNumRows; row++) {
            for (int column = 0; column < interpNumColumns; column++) {
                if (column % 2 == 1 && row % 2 == 1)
                    continue;
                pT = ImageUtils::to1D(row, column, interpNumColumns);
                if (column == 0 || column == interpNumColumns - 1 || row == 0 || row == interpNumRows - 1) {
                    max = median;
                    min = median;
                } else {
                    if (column % 2 == 0 && row % 2 == 0) {
                        adjColumn = adjCircleColumn;
                        adjRow = adjCircleRow;
                        adjSize = 4;
                    } else if (column % 2 == 0 && row % 2 == 1) {
                        adjColumn = adjRetVerColumn;
                        adjRow = adjRetVerRow;
                        adjSize = 2;
                    } else if (column % 2 == 1 && row % 2 == 0) {
                        adjColumn = adjRetHorColumn;
                        adjRow = adjRetHorRow;
                        adjSize = 2;
                    }

                    min = std::numeric_limits<int>::max();
                    max = std::numeric_limits<int>::min();
                    for (int i = 0; i < adjSize; i++) {
                        qRow = row + adjRow[i];
                        qColumn = column + adjColumn[i];

                        if (qRow >= 0 && qColumn >= 0 && qRow < interpNumRows && qColumn < interpNumColumns) {
                            qT = ImageUtils::to1D(qRow, qColumn, interpNumColumns);

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
     * @param imgPtr Image data.
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
            TopographicConvention paddedConvention = convention_;
            paddedConvention.domainExtension = TopographicDomainExtension::ExteriorRing;
            paddedConvention.infinityPixel = 0;
            // The padded producer only supplies interpolation buffers, which are
            // altitude-encoding independent. Declaring the doubled encoding keeps
            // it constructible for a convention whose own encoding is 8-bit.
            paddedConvention.altitudeEncoding = TopographicAltitudeEncoding::ExactDoubled;
            TreeOfShapesProducer paddedProducer(std::move(paddedConvention));
            auto [paddedMin, paddedMax, paddedAdjacency] = paddedProducer.interpolateImage4c8c(imgPtr);
            const int paddedRows = paddedProducer.interpolatedNumRows(imgPtr->getNumRows());
            const int paddedColumns = paddedProducer.interpolatedNumColumns(imgPtr->getNumColumns());
            return cropExteriorInterpolation(std::move(paddedMin), std::move(paddedMax), std::move(paddedAdjacency), paddedRows, paddedColumns, true);
        }
        auto img = imgPtr->rawData();
        int numRows = imgPtr->getNumRows();
        int numColumns = imgPtr->getNumColumns();
        if (numRows <= 0 || numColumns <= 0) {
            throw std::invalid_argument("TreeOfShapesProducer requires a non-empty image.");
        }

        int interpNumColumns = interpolatedNumColumns(numColumns);
        int interpNumRows = interpolatedNumRows(numRows);
        int size = checkedDomainSize(interpNumRows, interpNumColumns);
        AdjacencyUC adj(interpNumRows, interpNumColumns, true);

        // Allocate interpolation result buffers for minimum and maximum levels.
        std::vector<ToSGrayLevel> interpolationMin(size);
        std::vector<ToSGrayLevel> interpolationMax(size);

        PixelId pT;
        // Compute interval from 2-faces.
        for (PixelId pixel = 0; pixel < numColumns * numRows; ++pixel) {
            auto [row, column] = ImageUtils::to2D(pixel, numColumns);

            // Compute the interpolated-image index.
            pT = ImageUtils::to1D(originalPointRow(row), originalPointColumn(column), interpNumColumns);

            // Assign interpolation values.
            interpolationMin[pT] = interpolationMax[pT] = img[pixel];
        }

        auto getValue = [&](int row, int column) -> int {
            int origRow = (row - 1) / 2;
            int originalColumn = (column - 1) / 2;
            return img[ImageUtils::to1D(origRow, originalColumn, numColumns)];
        };

        // Borders.
        for (int row = 0; row < interpNumRows; row++) {
            int column;
            if (row % 2 == 1) { // horizontal e vertical
                column = 0;
                int v1 = getValue(row, column + 1);
                interpolationMin[ImageUtils::to1D(row, column, interpNumColumns)] = v1;
                interpolationMax[ImageUtils::to1D(row, column, interpNumColumns)] = v1;

                column = interpNumColumns - 1;
                v1 = getValue(row, column - 1);
                interpolationMin[ImageUtils::to1D(row, column, interpNumColumns)] = v1;
                interpolationMax[ImageUtils::to1D(row, column, interpNumColumns)] = v1;
            } else { // circulos
                if (row == 0) {
                    column = 0;
                    int v1 = getValue(row + 1, column + 1);
                    interpolationMin[ImageUtils::to1D(row, column, interpNumColumns)] = v1;
                    interpolationMax[ImageUtils::to1D(row, column, interpNumColumns)] = v1;

                    column = interpNumColumns - 1;
                    v1 = getValue(row + 1, column - 1);
                    interpolationMin[ImageUtils::to1D(row, column, interpNumColumns)] = v1;
                    interpolationMax[ImageUtils::to1D(row, column, interpNumColumns)] = v1;

                } else if (row == interpNumRows - 1) {
                    column = 0;
                    int v1 = getValue(row - 1, 1);
                    interpolationMin[ImageUtils::to1D(row, column, interpNumColumns)] = v1;
                    interpolationMax[ImageUtils::to1D(row, column, interpNumColumns)] = v1;

                    column = interpNumColumns - 1;
                    v1 = getValue(row - 1, column - 1);
                    interpolationMin[ImageUtils::to1D(row, column, interpNumColumns)] = v1;
                    interpolationMax[ImageUtils::to1D(row, column, interpNumColumns)] = v1;
                } else {
                    column = 0;
                    int v1 = getValue(row - 1, column + 1);
                    int v2 = getValue(row + 1, column + 1);
                    interpolationMin[ImageUtils::to1D(row, 0, interpNumColumns)] = std::min(v1, v2);
                    interpolationMax[ImageUtils::to1D(row, 0, interpNumColumns)] = std::max(v1, v2);

                    column = interpNumColumns - 1;
                    v1 = getValue(row - 1, column - 1);
                    v2 = getValue(row + 1, column - 1);
                    interpolationMin[ImageUtils::to1D(row, column, interpNumColumns)] = std::min(v1, v2);
                    interpolationMax[ImageUtils::to1D(row, column, interpNumColumns)] = std::max(v1, v2);
                }
            }
        }

        for (int column = 1; column < interpNumColumns - 1; column++) {
            int row;
            if (column % 2 == 1) { // horizontal e vertical
                row = 0;
                int v1 = getValue(row + 1, column);
                interpolationMin[ImageUtils::to1D(row, column, interpNumColumns)] = v1;
                interpolationMax[ImageUtils::to1D(row, column, interpNumColumns)] = v1;

                row = interpNumRows - 1;
                v1 = getValue(row - 1, column);
                interpolationMin[ImageUtils::to1D(row, column, interpNumColumns)] = v1;
                interpolationMax[ImageUtils::to1D(row, column, interpNumColumns)] = v1;
            } else { // circulos
                row = 0;
                int v1 = getValue(row + 1, column - 1);
                int v2 = getValue(row + 1, column + 1);
                interpolationMin[ImageUtils::to1D(row, column, interpNumColumns)] = std::min(v1, v2);
                interpolationMax[ImageUtils::to1D(row, column, interpNumColumns)] = std::max(v1, v2);

                row = interpNumRows - 1;
                v1 = getValue(row - 1, column - 1);
                v2 = getValue(row - 1, column + 1);
                interpolationMin[ImageUtils::to1D(row, column, interpNumColumns)] = std::min(v1, v2);
                interpolationMax[ImageUtils::to1D(row, column, interpNumColumns)] = std::max(v1, v2);
            }
        }

        // Compute interval from 1-faces
        for (int row = 1; row < interpNumRows - 1; row++) {
            for (int column = 1; column < interpNumColumns - 1; column++) {
                if (row % 2 == 1 && column % 2 == 1)
                    continue; // Already defined.

                pT = ImageUtils::to1D(row, column, interpNumColumns);
                if (column % 2 == 0 && row % 2 == 1) {
                    int v1 = getValue(row, column + 1);
                    int v2 = getValue(row, column - 1);
                    interpolationMin[pT] = std::min(v1, v2);
                    interpolationMax[pT] = std::max(v1, v2);
                } else if (column % 2 == 1 && row % 2 == 0) {
                    int v1 = getValue(row + 1, column);
                    int v2 = getValue(row - 1, column);
                    interpolationMin[pT] = std::min(v1, v2);
                    interpolationMax[pT] = std::max(v1, v2);
                }
            }
        }
        // Compute interval from 0-faces
        for (int row = 1; row < interpNumRows - 1; row++) {
            for (int column = 1; column < interpNumColumns - 1; column++) {
                if (row % 2 == 1 && column % 2 == 1)
                    continue; // Already defined.
                pT = ImageUtils::to1D(row, column, interpNumColumns);
                if (row % 2 == 0 && column % 2 == 0) {
                    // | v0 | v1 |
                    // | v2 | v3 |
                    int v0 = getValue(row - 1, column - 1);
                    int v1 = getValue(row + 1, column - 1);
                    int v2 = getValue(row - 1, column + 1);
                    int v3 = getValue(row + 1, column + 1);

                    const int min_v0v3 = std::min(v0, v3);
                    const int max_v0v3 = std::max(v0, v3);
                    const int min_v1v2 = std::min(v1, v2);
                    const int max_v1v2 = std::max(v1, v2);
                    const bool diagonal0IsHigh = min_v0v3 > max_v1v2;
                    const bool diagonal1IsHigh = min_v1v2 > max_v0v3;

                    if (diagonal0IsHigh || diagonal1IsHigh) {
                        const bool chooseDiagonal0 = usesHighDiagonalAtSaddle() ? diagonal0IsHigh : !diagonal0IsHigh;
                        if (chooseDiagonal0) {
                            setDiagonal0Connection(adj, row, column);
                            interpolationMin[pT] = min_v0v3;
                            interpolationMax[pT] = max_v0v3;
                        } else {
                            setDiagonal1Connection(adj, row, column);
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
     * @param imgPtr Image data.
     * @return Number of rows in the interpolated domain.
     */
    FloodResult floodImage(const ImageUInt8Ptr& imgPtr) const {
        if (!imgPtr) {
            throw std::invalid_argument("TreeOfShapesProducer requires a non-null image.");
        }

        int numRows = imgPtr->getNumRows();
        int numColumns = imgPtr->getNumColumns();

        int interpNumColumns = interpolatedNumColumns(numColumns);
        int interpNumRows = interpolatedNumRows(numRows);
        int size = checkedDomainSize(interpNumRows, interpNumColumns);
        const bool connectivityMap = usesConnectivityMap();
        auto [interpolationMin, interpolationMax, adj] = connectivityMap ? interpolateImage4c8c(imgPtr) : interpolateImage(imgPtr);

        std::vector<uint8_t> dejavu(size, 0);
        std::vector<PixelId> imgR(size);           // Ordered pixels.
        std::vector<ToSFloodDepth> imgU(size);     // Monotone max-tree levels.
        std::vector<ToSGrayLevel> grayLevel(size); // Actual propagation levels.

        PriorityQueueToS queue(connectivityMap ? ToSUInt8Depth : ToSSelfDualDepth); // Priority queue.
        const PixelId infinityPixel = validatedInfinityPixel(interpNumRows, interpNumColumns);
        int priorityQueueOld = interpolationMin[infinityPixel];
        queue.initial(infinityPixel, priorityQueueOld);
        dejavu[infinityPixel] = true;

        int order = 0;
        ToSFloodDepth depth = 0;
        while (!queue.isEmpty()) {
            const PixelId pixel = queue.priorityPop();      // Pop the element with highest priority.
            int priorityQueue = queue.getCurrentPriority(); // Current priority.
            if (connectivityMap) {
                if (priorityQueue != priorityQueueOld)
                    depth++;
                imgU[pixel] = depth;
            } else {
                imgU[pixel] = static_cast<ToSFloodDepth>(priorityQueue);
            }
            grayLevel[pixel] = static_cast<ToSGrayLevel>(priorityQueue);

            // Store h in the correct output order.
            imgR[order++] = pixel;

            // Adjacencies.
            for (PixelId neighbor : adj.getNeighborIndices(pixel)) {
                if (!dejavu[neighbor]) {
                    queue.priorityPush(neighbor, interpolationMin[neighbor], interpolationMax[neighbor]);
                    dejavu[neighbor] = true; // Mark as processed.
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
     * @param imgPtr Image data.
     * @return Resulting canonical flood tree.
     */
    CanonicalFloodTree buildCanonicalFloodTree(const ImageUInt8Ptr& imgPtr) const {
        FloodResult flood = floodImage(imgPtr);
        const int numPixelsInterp = static_cast<int>(flood.treeLevel.size());

        std::vector<PixelId> zPar(static_cast<size_t>(numPixelsInterp), InvalidPixel);
        std::vector<PixelId> parentInterpolate(static_cast<size_t>(numPixelsInterp), InvalidPixel);
        auto findRoot = [&](PixelId pStar) {
            while (zPar[pStar] != pStar) {
                zPar[pStar] = zPar[zPar[pStar]];
                pStar = zPar[pStar];
            }
            return pStar;
        };

        for (int i = numPixelsInterp - 1; i >= 0; --i) {
            const PixelId pStar = flood.order[static_cast<size_t>(i)];
            parentInterpolate[pStar] = pStar;
            zPar[pStar] = pStar;
            for (PixelId qStar : flood.adjacency.getNeighborIndices(pStar)) {
                if (zPar[qStar] != InvalidPixel) {
                    const PixelId rStar = findRoot(qStar);
                    if (pStar != rStar) {
                        parentInterpolate[rStar] = pStar;
                        zPar[rStar] = pStar;
                    }
                }
            }
        }

        auto sameLevel = [&](PixelId aStar, PixelId bStar) { return flood.treeLevel[aStar] == flood.treeLevel[bStar]; };
        for (PixelId pStar : flood.order) {
            const PixelId qStar = parentInterpolate[pStar];
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
     * @param imgPtr Image data.
     * @return Values produced by the operation.
     */
    std::tuple<std::vector<ToSFloodDepth>, std::vector<PixelId>, AdjacencyUC> sort(const ImageUInt8Ptr& imgPtr) const {
        FloodResult flood = floodImage(imgPtr);
        return std::make_tuple(std::move(flood.treeLevel), std::move(flood.order), std::move(flood.adjacency));
    }

    // Tests whether an interpolated pixel corresponds to an original pixel.
    /**
     * @brief Checks whether an interpolated pixel belongs to the original image grid.
     *
     * @param pixel Row-major pixel index to transform or inspect.
     * @param interpNumColumns Number of columns in the interpolated domain.
     * @return True when the interpolated index represents an original image pixel.
     */
    inline bool isOriginal1D(PixelId pixel, int interpNumColumns) const {
        const int row = pixel / interpNumColumns;
        const int column = pixel - row * interpNumColumns;
        const int offset = interpolationOffset();
        return (row & 1) == offset && (column & 1) == offset;
    }

    // Maps an interpolated pixel to the original image domain.
    /**
     * @brief Maps an interpolated-domain pixel index to the original image domain.
     *
     * @param pStar Pixel index in the interpolated domain.
     * @param interNumColumns Number of columns in the interpolated domain.
     * @param numColumns Number of columns in the domain.
     * @return Padding offset in interpolated-domain coordinates.
     */
    inline PixelId toOriginal1D(PixelId pStar, int interNumColumns, int numColumns) const {
        int r = pStar / interNumColumns;
        int c = pStar - r * interNumColumns; // Avoid the modulo operator.
        const int offset = interpolationOffset();
        return ((r - offset) >> 1) * numColumns + ((c - offset) >> 1);
    }
    /// @endcond

  private:
    /**
     * @brief Builds a native Tree-of-Shapes topology.
     *
     * The projection retains every interpolated plateau that owns an original
     * pixel. A plateau without original pixels is retained exactly when it
     * combines at least two distinct projected child supports. Consequently,
     * unary empty plateaus are contracted while a virtual branching root is
     * represented faithfully with an empty direct proper part.
     *
     * @tparam Altitude Published node-altitude type.
     * @param imgPtr Image.
     * @return The resulting native Tree-of-Shapes topology.
     */
    template <class Altitude> [[nodiscard]] TreeOfShapesBuildResult<Altitude> buildTreeOfShapes(const ImageUInt8Ptr& imgPtr) const {
        if (!imgPtr) {
            throw std::invalid_argument("TreeOfShapesProducer requires a non-null image.");
        }
        const int numRows = imgPtr->getNumRows();
        const int numColumns = imgPtr->getNumColumns();
        const int numPixels = checkedDomainSize(numRows, numColumns);
        const int interpNumColumns = interpolatedNumColumns(numColumns);

        CanonicalFloodTree canonical = buildCanonicalFloodTree(imgPtr);
        const int numPixelsInterp = static_cast<int>(canonical.parent.size());
        auto sameLevel = [&](PixelId aStar, PixelId bStar) {
            return canonical.treeLevel[static_cast<size_t>(aStar)] == canonical.treeLevel[static_cast<size_t>(bStar)];
        };
        auto repOf = [&](PixelId pStar) {
            const PixelId parentStar = canonical.parent[static_cast<size_t>(pStar)];
            return (parentStar == pStar || sameLevel(parentStar, pStar)) ? parentStar : pStar;
        };

        std::vector<PixelId> canonicalReps;
        canonicalReps.reserve(static_cast<size_t>(numPixels));
        for (PixelId pStar : canonical.order) {
            const PixelId parentStar = canonical.parent[static_cast<size_t>(pStar)];
            if (parentStar == pStar || !sameLevel(parentStar, pStar)) {
                canonicalReps.push_back(pStar);
            }
        }
        if (canonicalReps.empty()) {
            throw std::runtime_error("Tree-of-shapes construction produced no canonical interpolated node.");
        }

        const PixelId rawRoot = canonicalReps.front();
        std::vector<PixelId> rawParent(static_cast<size_t>(numPixelsInterp), InvalidPixel);
        std::vector<PixelId> firstRawChild(static_cast<size_t>(numPixelsInterp), InvalidPixel);
        std::vector<PixelId> nextRawSibling(static_cast<size_t>(numPixelsInterp), InvalidPixel);
        for (PixelId repStar : canonicalReps) {
            if (repStar == rawRoot) {
                rawParent[static_cast<size_t>(repStar)] = repStar;
                continue;
            }
            const PixelId parentRep = repOf(canonical.parent[static_cast<size_t>(repStar)]);
            if (parentRep < 0 || parentRep >= numPixelsInterp || rawParent[static_cast<size_t>(parentRep)] == InvalidPixel) {
                throw std::runtime_error("Tree-of-shapes canonical node has no canonical parent.");
            }
            rawParent[static_cast<size_t>(repStar)] = parentRep;
            nextRawSibling[static_cast<size_t>(repStar)] = firstRawChild[static_cast<size_t>(parentRep)];
            firstRawChild[static_cast<size_t>(parentRep)] = repStar;
        }

        std::vector<int> directCount(static_cast<size_t>(numPixelsInterp), 0);
        std::vector<PixelId> representativePixelByElement(static_cast<size_t>(numPixels), InvalidPixel);
        for (PixelId pStar = 0; pStar < numPixelsInterp; ++pStar) {
            if (!isOriginal1D(pStar, interpNumColumns)) {
                continue;
            }
            const PixelId elementId = toOriginal1D(pStar, interpNumColumns, numColumns);
            if (elementId < 0 || elementId >= numPixels) {
                throw std::runtime_error("Tree-of-shapes projection produced an element outside the source domain.");
            }
            const PixelId representativePixel = repOf(pStar);
            representativePixelByElement[static_cast<size_t>(elementId)] = representativePixel;
            ++directCount[static_cast<size_t>(representativePixel)];
        }
        if (std::find(representativePixelByElement.begin(), representativePixelByElement.end(), InvalidPixel) != representativePixelByElement.end()) {
            throw std::runtime_error("Tree-of-shapes projection did not assign every source-domain element.");
        }

        // The canonical depth, parent, and full interpolated order are no
        // longer needed after the raw hierarchy and source-domain smallest nodes have
        // been extracted. Release them before allocating the projection buffer.
        std::vector<ToSFloodDepth>().swap(canonical.treeLevel);
        std::vector<PixelId>().swap(canonical.parent);
        std::vector<PixelId>().swap(canonical.order);

        // Bottom-up support projection. projectedRep[r] is the canonical
        // representative of the distinct projected support rooted at r.
        std::vector<PixelId> projectedRep(static_cast<size_t>(numPixelsInterp), InvalidPixel);
        for (auto it = canonicalReps.rbegin(); it != canonicalReps.rend(); ++it) {
            const PixelId repStar = *it;
            int numProjectedChildren = 0;
            PixelId onlyProjectedChild = InvalidPixel;
            for (PixelId childStar = firstRawChild[static_cast<size_t>(repStar)]; childStar != InvalidPixel;
                 childStar = nextRawSibling[static_cast<size_t>(childStar)]) {
                const PixelId childProjection = projectedRep[static_cast<size_t>(childStar)];
                if (childProjection != InvalidPixel) {
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

        const PixelId projectedRootRep = projectedRep[static_cast<size_t>(rawRoot)];
        if (projectedRootRep == InvalidPixel || projectedRep[static_cast<size_t>(projectedRootRep)] != projectedRootRep) {
            throw std::runtime_error("Tree-of-shapes projection produced an empty hierarchy.");
        }

        std::vector<PixelId>().swap(firstRawChild);
        std::vector<PixelId>().swap(nextRawSibling);
        std::vector<int>().swap(directCount);

        std::vector<NodeId> nodeIdByRep(static_cast<size_t>(numPixelsInterp), InvalidNode);
        int numNodes = 0;
        for (PixelId repStar : canonicalReps) {
            if (projectedRep[static_cast<size_t>(repStar)] == repStar) {
                nodeIdByRep[static_cast<size_t>(repStar)] = numNodes++;
            }
        }

        TreeOfShapesBuildResult<Altitude> result;
        result.parent.assign(static_cast<size_t>(numNodes), InvalidNode);
        result.smallestNodeMap.assign(static_cast<size_t>(numPixels), InvalidNode);
        result.nodeAltitudes.assign(static_cast<size_t>(numNodes), Altitude{});
        result.root = nodeIdByRep[static_cast<size_t>(projectedRootRep)];
        result.numRows = numRows;
        result.numColumns = numColumns;
        detail::TopologicalNativeHierarchyRecorder proofRecorder(static_cast<std::size_t>(numNodes), static_cast<std::size_t>(numPixels), result.root);

        for (PixelId repStar : canonicalReps) {
            const NodeId nodeId = nodeIdByRep[static_cast<size_t>(repStar)];
            if (nodeId == InvalidNode) {
                continue;
            }

            if (repStar == projectedRootRep) {
                result.parent[static_cast<size_t>(nodeId)] = nodeId;
            } else {
                PixelId ancestorRep = rawParent[static_cast<size_t>(repStar)];
                while (ancestorRep != InvalidPixel && nodeIdByRep[static_cast<size_t>(ancestorRep)] == InvalidNode) {
                    if (ancestorRep == rawParent[static_cast<size_t>(ancestorRep)]) {
                        ancestorRep = InvalidPixel;
                        break;
                    }
                    ancestorRep = rawParent[static_cast<size_t>(ancestorRep)];
                }
                if (ancestorRep == InvalidPixel) {
                    throw std::runtime_error("Tree-of-shapes retained node has no retained parent.");
                }
                result.parent[static_cast<size_t>(nodeId)] = nodeIdByRep[static_cast<size_t>(ancestorRep)];
            }
            proofRecorder.recordSupportedNode(nodeId, result.parent[static_cast<std::size_t>(nodeId)]);

            result.nodeAltitudes[static_cast<size_t>(nodeId)] =
                ToSAltitudeEncodingFor<Altitude>::type::encode(canonical.grayLevel[static_cast<size_t>(repStar)], immersionMode_);
        }

        for (PixelId elementId = 0; elementId < numPixels; ++elementId) {
            const PixelId representativePixel = representativePixelByElement[static_cast<size_t>(elementId)];
            const NodeId smallestNodeId = nodeIdByRep[static_cast<size_t>(representativePixel)];
            if (smallestNodeId == InvalidNode) {
                throw std::runtime_error("Tree-of-shapes proper part belongs to a contracted node.");
            }
            result.smallestNodeMap[static_cast<size_t>(elementId)] = smallestNodeId;
            proofRecorder.recordProperPart(elementId, smallestNodeId);
        }
        for (NodeId nodeId = 0; nodeId < numNodes; ++nodeId) {
            const NodeId parentNodeId = result.parent[static_cast<std::size_t>(nodeId)];
            if (nodeId != result.root &&
                result.nodeAltitudes[static_cast<std::size_t>(nodeId)] == result.nodeAltitudes[static_cast<std::size_t>(parentNodeId)]) {
                throw std::runtime_error("Tree-of-shapes construction produced an equal-altitude parent-child edge.");
            }
        }
        result.topologyProof = std::move(proofRecorder).finish();
        return result;
    }

  public:
    /**
     * @brief Builds the hierarchy in the altitude scale declared by the convention.
     *
     * Under `TopographicAltitudeEncoding::UInt8` the published altitudes are
     * unchanged source gray levels. This encoding is exact for a
     * complementary-grid immersion and for a self-dual span immersion without
     * an exterior ring. Under `TopographicAltitudeEncoding::ExactDoubled`, even
     * values represent source gray levels and odd values represent half levels
     * introduced by self-dual interpolation. Neither encoding loses a
     * parent-child altitude distinction to quantization.
     *
     * @tparam Altitude Published node-altitude type, which must match the
     * altitude encoding declared by the convention.
     * @param imgPtr Image.
     * @return The resulting hierarchy in the declared altitude scale.
     */
    template <class Altitude = std::uint8_t> [[nodiscard]] TreeOfShapesBuildResult<Altitude> build(const ImageUInt8Ptr& imgPtr) const {
        if (convention_.altitudeEncoding != ToSAltitudeEncodingFor<Altitude>::declaration) {
            throw std::invalid_argument("Tree-of-shapes altitude type does not match the altitude encoding declared by the topographic convention.");
        }
        validateConventionForImage(imgPtr);
        return buildTreeOfShapes<Altitude>(imgPtr);
    }
};

} // namespace mmcfilters
