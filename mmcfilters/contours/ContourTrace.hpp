#pragma once

#include "../utils/Common.hpp"

#include <cstddef>
#include <cstdint>
#include <iterator>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcfilters {

class ContourTraceComputation;
namespace contours::detail {
class ContourTraceTraversal;
}

/**
 * @brief Side of a support pixel occupied by a contour edge.
 */
enum class ContourSide : uint8_t { North = 0, West = 1, East = 2, South = 3 };

/**
 * @brief Whether a contour boundary is external or surrounds a hole.
 */
enum class ContourBoundaryKind : uint8_t { External, Internal };

/**
 * @brief One contour edge represented by its support pixel and side.
 */
struct ContourEdge {
    /// Row-major support-pixel index incident to the boundary edge.
    PixelId pixel = InvalidPixel;
    /// Side of the support pixel occupied by the boundary edge.
    ContourSide side = ContourSide::North;

    /// Compares the support pixel and side.
    friend bool operator==(const ContourEdge&, const ContourEdge&) = default;
};

/**
 * @brief Metadata for one ordered external or internal contour boundary.
 */
struct ContourBoundary {
    /// Whether this boundary is external or surrounds a hole.
    ContourBoundaryKind kind = ContourBoundaryKind::External;
    /// First edge in the trace edge buffer.
    uint32_t edgeOffset = 0;
    /// Number of consecutive edges in the boundary.
    uint32_t edgeCount = 0;
    /// Doubled signed area under the contour orientation convention.
    int doubledSignedArea = 0;
};

namespace contours::detail {

/** @brief Packs one pixel-side edge into a compact integer. */
[[nodiscard]] inline int packContourEdge(PixelId pixel, ContourSide side) {
    return (4 * pixel) + static_cast<int>(side);
}

/** @brief Unpacks one compact pixel-side edge. */
[[nodiscard]] inline ContourEdge unpackContourEdge(int packedEdge) {
    if (packedEdge < 0) {
        return {};
    }
    return ContourEdge{packedEdge / 4, static_cast<ContourSide>(packedEdge & 3)};
}

} // namespace contours::detail

/**
 * @brief Immutable range over unpacked contour edges.
 */
class ContourEdgeRange {
  public:
    /**
     * @brief Forward iterator that unpacks contour edges on dereference.
     */
    class iterator {
      public:
        /// Standard category for a multi-pass forward iterator.
        using iterator_category = std::forward_iterator_tag;
        /// Unpacked edge value yielded by dereference.
        using value_type = ContourEdge;
        /// Signed iterator-distance type.
        using difference_type = std::ptrdiff_t;
        /// No pointer type is exposed because dereference returns a value.
        using pointer = void;
        /// Value-returning reference type.
        using reference = value_type;

        iterator() = default;

        /**
         * @brief Creates an iterator at one packed-edge position.
         * @param packedEdge Current packed-edge address.
         */
        explicit iterator(const int* packedEdge) : packedEdge_(packedEdge) {}

        /**
         * @brief Returns the unpacked current edge.
         * @return Current contour edge.
         */
        [[nodiscard]] value_type operator*() const { return contours::detail::unpackContourEdge(*packedEdge_); }

        /**
         * @brief Advances to the next packed edge.
         * @return This iterator after advancing.
         */
        iterator& operator++() {
            ++packedEdge_;
            return *this;
        }

        /**
         * @brief Advances while retaining the previous iterator position.
         * @return Iterator at the position preceding the advance.
         */
        iterator operator++(int) {
            iterator previous(*this);
            ++(*this);
            return previous;
        }

        /// Compares packed-edge positions.
        friend bool operator==(const iterator&, const iterator&) = default;

      private:
        /// Address of the current packed edge.
        const int* packedEdge_ = nullptr;
    };

    ContourEdgeRange() = default;

    /**
     * @brief Creates a view over packed contour edges.
     * @param packedEdges Packed contour-edge span.
     */
    explicit ContourEdgeRange(std::span<const int> packedEdges) : packedEdges_(packedEdges) {}

    /** @brief Returns an iterator to the first contour edge. @return First edge position. */
    [[nodiscard]] iterator begin() const { return iterator(packedEdges_.data()); }

    /** @brief Returns the exclusive end iterator. @return Position following the final edge. */
    [[nodiscard]] iterator end() const {
        return packedEdges_.empty() ? begin() : iterator(packedEdges_.data() + packedEdges_.size());
    }

    /** @brief Returns whether the range contains no edges. @return True when the range is empty. */
    [[nodiscard]] bool empty() const noexcept { return packedEdges_.empty(); }

    /** @brief Returns the number of contour edges. @return Edge count. */
    [[nodiscard]] std::size_t size() const noexcept { return packedEdges_.size(); }

  private:
    std::span<const int> packedEdges_; ///< Packed contour edges viewed by this range.
};

/**
 * @brief Immutable projection of ordered contour edges onto support pixels.
 *
 * One pixel is yielded for each edge. A pixel can therefore occur more than
 * once when different edges occupy different sides of that pixel.
 */
class ContourPixelRange {
  public:
    /** @brief Forward iterator that returns the support pixel of each edge. */
    class iterator {
      public:
        /// Standard category for a multi-pass forward iterator.
        using iterator_category = std::forward_iterator_tag;
        /// Pixel identifier yielded by dereference.
        using value_type = PixelId;
        /// Signed iterator-distance type.
        using difference_type = std::ptrdiff_t;
        /// No pointer type is exposed because dereference returns a value.
        using pointer = void;
        /// Value-returning reference type.
        using reference = value_type;

        iterator() = default;

        /** @brief Creates a projection over an edge iterator. @param edgeIterator Current edge position. */
        explicit iterator(ContourEdgeRange::iterator edgeIterator) : edgeIterator_(edgeIterator) {}

        /** @brief Returns the support pixel of the current edge. @return Current support pixel. */
        [[nodiscard]] value_type operator*() const { return (*edgeIterator_).pixel; }

        /** @brief Advances to the next edge pixel. @return This iterator after advancing. */
        iterator& operator++() {
            ++edgeIterator_;
            return *this;
        }

        /**
         * @brief Advances while retaining the previous iterator position.
         * @return Iterator at the position preceding the advance.
         */
        iterator operator++(int) {
            iterator previous(*this);
            ++(*this);
            return previous;
        }

        /// Compares the underlying edge positions.
        friend bool operator==(const iterator&, const iterator&) = default;

      private:
        ContourEdgeRange::iterator edgeIterator_; ///< Current ordered edge position.
    };

    ContourPixelRange() = default;

    /** @brief Creates a support-pixel projection of an edge range. @param edges Ordered source edges. */
    explicit ContourPixelRange(ContourEdgeRange edges) : edges_(edges) {}

    /** @brief Returns an iterator to the first edge pixel. @return First pixel position. */
    [[nodiscard]] iterator begin() const { return iterator(edges_.begin()); }

    /** @brief Returns the exclusive end iterator. @return Position following the final pixel. */
    [[nodiscard]] iterator end() const { return iterator(edges_.end()); }

    /** @brief Returns whether the range contains no pixels. @return True when the range is empty. */
    [[nodiscard]] bool empty() const noexcept { return edges_.empty(); }

    /** @brief Returns the number of projected edge pixels, including repetitions. @return Projected pixel count. */
    [[nodiscard]] std::size_t size() const noexcept { return edges_.size(); }

  private:
    ContourEdgeRange edges_; ///< Ordered edges projected by this range.
};

class ContourTrace;

/**
 * @brief Borrowed ordered contour trace for one tree node.
 *
 * A view contains all ordered external and internal boundaries without owning
 * their storage. A computation iterator's view remains valid until that
 * iterator or one of its copies advances, or their shared traversal is
 * destroyed. A callback's view is valid only during the callback. Views of
 * an owned ContourTrace require its edge and boundary storage to remain alive
 * and unchanged. The same limits apply to spans and ranges obtained from a view.
 */
class ContourTraceView {
  public:
    ContourTraceView() = default;

    /** @brief Returns all ordered boundary descriptors. @return Borrowed boundary descriptor span. */
    [[nodiscard]] std::span<const ContourBoundary> boundaries() const noexcept { return boundaries_; }

    /**
     * @brief Returns the unique external boundary.
     * @return Unique external boundary descriptor.
     * @throws std::logic_error If zero or multiple external boundaries exist.
     */
    [[nodiscard]] ContourBoundary externalBoundary() const {
        const ContourBoundary* external = nullptr;
        for (const ContourBoundary& boundary : boundaries_) {
            if (boundary.kind != ContourBoundaryKind::External) {
                continue;
            }
            if (external != nullptr) {
                throw std::logic_error("ContourTraceView::externalBoundary requires exactly one external boundary.");
            }
            external = &boundary;
        }
        if (external == nullptr) {
            throw std::logic_error("ContourTraceView::externalBoundary requires exactly one external boundary.");
        }
        return *external;
    }

    /** @brief Returns every contour edge in boundary order. @return Ordered edge range. */
    [[nodiscard]] ContourEdgeRange edges() const noexcept { return ContourEdgeRange(packedEdges_); }

    /**
     * @brief Returns the ordered edges of one boundary.
     * @param boundary Boundary descriptor from this trace.
     * @return Ordered edge range selected by the descriptor.
     */
    [[nodiscard]] ContourEdgeRange boundaryEdges(const ContourBoundary& boundary) const {
        return ContourEdgeRange(boundarySpan(boundary));
    }

    /**
     * @brief Returns one support pixel for each ordered boundary edge.
     * @param boundary Boundary descriptor from this trace.
     * @return Ordered support-pixel projection, including repetitions.
     */
    [[nodiscard]] ContourPixelRange boundaryPixels(const ContourBoundary& boundary) const {
        return ContourPixelRange(boundaryEdges(boundary));
    }

  private:
    friend class ContourTrace;
    friend class ContourTraceComputation;
    friend class contours::detail::ContourTraceTraversal;

    /**
     * @brief Creates a borrowed trace over packed edges and descriptors.
     * @param packedEdges Ordered packed contour edges.
     * @param boundaries Boundary descriptors over `packedEdges`.
     */
    ContourTraceView(std::span<const int> packedEdges, std::span<const ContourBoundary> boundaries)
        : packedEdges_(packedEdges), boundaries_(boundaries) {}

    /**
     * @brief Validates and resolves one boundary edge slice.
     * @param boundary Boundary descriptor to resolve.
     * @return Packed edge interval selected by `boundary`.
     */
    [[nodiscard]] std::span<const int> boundarySpan(const ContourBoundary& boundary) const {
        const std::size_t offset = boundary.edgeOffset;
        const std::size_t count = boundary.edgeCount;
        if (offset > packedEdges_.size() || count > packedEdges_.size() - offset) {
            throw std::invalid_argument("ContourBoundary does not belong to this contour trace.");
        }
        return packedEdges_.subspan(offset, count);
    }

    std::span<const int> packedEdges_;                    ///< Packed ordered boundary edges.
    std::span<const ContourBoundary> boundaries_;         ///< Boundary descriptors over `packedEdges_`.
};

/**
 * @brief Independently owned ordered contour trace for one tree node.
 *
 * Returned views, spans, and ranges borrow this trace's storage. Keep that
 * storage alive and unchanged while using them; replacing it by assignment
 * or destroying it invalidates those borrowed objects.
 */
class ContourTrace {
  public:
    ContourTrace() = default;

    /**
     * @brief Copies a borrowed trace so it can outlive its source traversal.
     * @param trace Borrowed source trace.
     */
    explicit ContourTrace(ContourTraceView trace)
        : packedEdges_(trace.packedEdges_.begin(), trace.packedEdges_.end()),
          boundaries_(trace.boundaries_.begin(), trace.boundaries_.end()) {}

    /** @brief Returns all ordered boundary descriptors. @return Borrowed boundary descriptor span. */
    [[nodiscard]] std::span<const ContourBoundary> boundaries() const noexcept { return view().boundaries(); }

    /**
     * @brief Returns the unique external boundary.
     * @return Unique external boundary descriptor.
     * @throws std::logic_error If zero or multiple external boundaries exist.
     */
    [[nodiscard]] ContourBoundary externalBoundary() const { return view().externalBoundary(); }

    /** @brief Returns every contour edge in boundary order. @return Ordered edge range. */
    [[nodiscard]] ContourEdgeRange edges() const noexcept { return view().edges(); }

    /**
     * @brief Returns the ordered edges of one boundary.
     * @param boundary Boundary descriptor from this trace.
     * @return Ordered edge range selected by the descriptor.
     */
    [[nodiscard]] ContourEdgeRange boundaryEdges(const ContourBoundary& boundary) const { return view().boundaryEdges(boundary); }

    /**
     * @brief Returns one support pixel for each ordered boundary edge.
     * @param boundary Boundary descriptor from this trace.
     * @return Ordered support-pixel projection, including repetitions.
     */
    [[nodiscard]] ContourPixelRange boundaryPixels(const ContourBoundary& boundary) const { return view().boundaryPixels(boundary); }

    /**
     * @brief Returns a borrowed view over this owned trace.
     * @return View borrowing the current edge and boundary storage. Replacing
     * or destroying that storage invalidates the view and its spans and ranges.
     */
    [[nodiscard]] ContourTraceView view() const noexcept { return ContourTraceView(packedEdges_, boundaries_); }

  private:
    friend class ContourTraceComputation;

    /**
     * @brief Takes ownership of packed ordered edges and boundary descriptors.
     * @param packedEdges Ordered packed contour edges.
     * @param boundaries Boundary descriptors over `packedEdges`.
     */
    ContourTrace(std::vector<int> packedEdges, std::vector<ContourBoundary> boundaries)
        : packedEdges_(std::move(packedEdges)), boundaries_(std::move(boundaries)) {}

    std::vector<int> packedEdges_;                ///< Packed ordered edges owned by this trace.
    std::vector<ContourBoundary> boundaries_;     ///< Descriptors over `packedEdges_`.
};

} // namespace mmcfilters
