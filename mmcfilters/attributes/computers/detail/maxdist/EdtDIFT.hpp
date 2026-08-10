#pragma once

#include <iostream>
#include <iomanip>
#include <utility>
#include <vector>
#include <string>

#include "../../../../utils/Image.hpp"
#include "../../../../utils/CommittedGridAccess.hpp"
#include "../../../../utils/CommittedImageAccess.hpp"
#include "../../../../utils/RegularGridAdjacency2D.hpp"
#include "PQueue.hpp"
#include "Geometry.hpp"

namespace mmcfilters::attributes::computers::detail::maxdist {
namespace detail {
/**
 * @brief Squared integer distance helper.
 *
 * @param value Value used by the operation.
 * @return Squared integer distance helper.
 */
inline int square(int value) noexcept { return value * value; }

/**
 * @brief Prints an uint8 image as a small debug table.
 *
 * @param image Image used by the operation.
 */
inline void printUInt8Image(const mmcfilters::ImageUInt8Ptr& image) {
    for (int row = 0; row < image->getNumRows(); ++row) {
        for (int col = 0; col < image->getNumCols(); ++col) {
            std::cout << std::setw(4) << static_cast<int>((*image)[row * image->getNumCols() + col]);
        }
        std::cout << "\n";
    }
}
} // namespace detail

/**
 * @brief Direction-aware adjacency used by the Image Foresting Transform.
 *
 * Each entry stores a geometric offset and the index of the next adaptive
 * adjacency stencil to use if propagation reaches that neighbour. This is a
 * compact way to encode the ordered 8-neighbour propagation rules used by
 * the dynamic Euclidean distance transform while avoiding invalid image
 * border accesses.
 */
class AdaptiveAdj {
  public:
    /**
     * @brief Lightweight iterable view of neighbours around one point.
     *
     * A Neighbors object does not own adjacency data. It points back to the
     * AdaptiveAdj that created it and exposes only the first `end` offsets,
     * which allows propagation to use a prefix of the full stencil.
     */
    class Neighbors {
      public:
        /**
         * @brief Forward iterator yielding `(neighbourPoint, nextAdjIndex)`.
         */
        class Iterator {
          public:
            /**
             * @brief Creates an iterator at a position inside a Neighbors view.
             *
             * @param neighbors Neighbour offsets used by the operation.
             * @param idx Zero-based index used by the operation.
             */
            Iterator(const Neighbors& neighbors, int idx = 0) : idx_{idx}, neighbors_{neighbors} {}

            /**
             * @brief Returns the current neighbour and next-stencil id.
             *
             * @return The current neighbour and next-stencil id.
             */
            std::pair<Point2D, int> operator*() const {
                Point2D p = neighbors_.point(idx_);
                int n = neighbors_.nextAdj(idx_);
                return std::make_pair(p, n);
            }

            /**
             * @brief Advances to the next offset in the view.
             *
             * @return Reference to the resulting object.
             */
            inline Iterator& operator++() noexcept {
                idx_++;
                return *this;
            }

            /**
             * @brief Iterators compare equal when they point to the same offset.
             *
             * @param other Object to compare with or transfer from.
             * @return True when the documented condition holds; otherwise false.
             */
            bool operator==(const Iterator& other) const noexcept { return idx_ == other.idx_; }

            /**
             * @brief Tests whether two objects differ.
             *
             * @param other Object compared with or copied from.
             * @return True when two objects differ; otherwise false.
             */
            inline bool operator!=(Iterator& other) const noexcept { return !(*this == other); }

          private:
            /**
             * @brief Current offset index inside the neighbour view.
             */
            int idx_;

            /**
             * @brief View being iterated.
             */
            const Neighbors& neighbors_;
        };

        /**
         * @brief Creates a neighbour view centered at point `p`.
         *
         * @param p Point used by the operation.
         * @param adj Adaptive adjacency stencil.
         * @param end Exclusive end position.
         */
        Neighbors(const Point2D& p, const AdaptiveAdj& adj, size_t end) : p_{p}, adj_{adj}, end_{end} {}

        /**
         * @brief Returns the neighbour point at offset index `idx`.
         *
         * @param idx Zero-based index used by the operation.
         * @return The neighbour point at offset index idx.
         */
        Point2D point(int idx) const { return p_ + adj_.offset_[idx]; }

        /**
         * @brief Returns the adaptive adjacency id associated with offset `idx`.
         *
         * @param idx Zero-based index used by the operation.
         * @return The adaptive adjacency id associated with offset idx.
         */
        inline int nextAdj(int idx) const { return adj_.nextAdj_[idx]; }

        /**
         * @brief Function-call shorthand for point(idx).
         *
         * @param idx Zero-based index used by the operation.
         * @return Point at the requested offset.
         */
        inline Point2D operator()(int idx) const { return point(idx); }

        /**
         * @brief Number of offsets exposed by this view.
         *
         * @return The number of offsets exposed by this view.
         */
        inline int size() const { return end_; }

        /**
         * @brief Iterator over the first exposed offset.
         *
         * @return Iterator positioned at the first exposed offset.
         */
        Iterator begin() const { return Iterator(*this, 0); }

        /**
         * @brief Sentinel iterator one past the last exposed offset.
         *
         * @return Iterator positioned one past the last exposed offset.
         */
        Iterator end() const { return Iterator(*this, end_); }

      private:
        /**
         * @brief Center point of the neighbour view.
         */
        Point2D p_;

        /**
         * @brief Stencil storage referenced by this view.
         */
        const AdaptiveAdj& adj_;

        /**
         * @brief One-past-last offset index exposed by this view.
         */
        size_t end_;
    };

    /**
     * @brief Constructs a stencil from offsets, next-stencil ids, and prefix size.
     *
     * `npropagation` is the number of offsets used during distance
     * propagation. The full offset vector can be larger because removal
     * traversals need to inspect all neighbours.
     *
     * @param offset Offset into the underlying storage.
     * @param nextAdj Next-stencil identifiers associated with each offset.
     * @param npropagation Number of offsets in the propagation prefix.
     */
    AdaptiveAdj(std::vector<Point2D> offset, std::vector<int> nextAdj, size_t npropagation)
        : offset_(std::move(offset)), nextAdj_(std::move(nextAdj)), npropagation_{npropagation} {}
    /**
     * @brief Constructs `AdaptiveAdj` from the supplied inputs.
     *
     * @param offset Offset used to access the neighboring sample or storage position.
     * @param nextAdj Identifier of the next adjacency relation in the propagation sequence.
     * @param npropagation Number of adjacency relations used during propagation.
     */
    AdaptiveAdj(std::vector<Point2D>&& offset, std::vector<int> nextAdj, size_t npropagation)
        : offset_(std::move(offset)), nextAdj_(std::move(nextAdj)), npropagation_{npropagation} {}
    /**
     * @brief Constructs `AdaptiveAdj` from the supplied inputs.
     *
     * @param offset Offset used to access the neighboring sample or storage position.
     * @param nextAdj Identifier of the next adjacency relation in the propagation sequence.
     * @param npropagation Number of adjacency relations used during propagation.
     */
    AdaptiveAdj(std::vector<Point2D> offset, std::vector<int>&& nextAdj, size_t npropagation)
        : offset_(std::move(offset)), nextAdj_(std::move(nextAdj)), npropagation_{npropagation} {}
    /**
     * @brief Constructs `AdaptiveAdj` from the supplied inputs.
     *
     * @param offset Offset used to access the neighboring sample or storage position.
     * @param nextAdj Identifier of the next adjacency relation in the propagation sequence.
     * @param npropagation Number of adjacency relations used during propagation.
     */
    AdaptiveAdj(std::vector<Point2D>&& offset, std::vector<int>&& nextAdj, size_t npropagation)
        : offset_(std::move(offset)), nextAdj_(std::move(nextAdj)), npropagation_{npropagation} {}

    /**
     * @brief Constructs `AdaptiveAdj` from the supplied inputs.
     *
     * @param offset Offset used to access the neighboring sample or storage position.
     * @param nextAdj Identifier of the next adjacency relation in the propagation sequence.
     * @param npropagation Number of adjacency relations used during propagation.
     */
    AdaptiveAdj(std::vector<Point2D> offset, std::initializer_list<int> nextAdj, size_t npropagation)
        : offset_(std::move(offset)), nextAdj_(nextAdj), npropagation_{npropagation} {}
    /**
     * @brief Constructs `AdaptiveAdj` from the supplied inputs.
     *
     * @param offset Offset used to access the neighboring sample or storage position.
     * @param nextAdj Identifier of the next adjacency relation in the propagation sequence.
     * @param npropagation Number of adjacency relations used during propagation.
     */
    AdaptiveAdj(std::initializer_list<Point2D> offset, std::vector<int> nextAdj, size_t npropagation)
        : offset_(offset), nextAdj_(std::move(nextAdj)), npropagation_{npropagation} {}
    /**
     * @brief Constructs `AdaptiveAdj` from the supplied inputs.
     *
     * @param offset Offset used to access the neighboring sample or storage position.
     * @param nextAdj Identifier of the next adjacency relation in the propagation sequence.
     * @param npropagation Number of adjacency relations used during propagation.
     */
    AdaptiveAdj(std::initializer_list<Point2D> offset, std::initializer_list<int> nextAdj, size_t npropagation)
        : offset_{offset}, nextAdj_{nextAdj}, npropagation_{npropagation} {}

    /**
     * @brief Constructs `AdaptiveAdj` from the supplied inputs.
     *
     * @param offset Offset used to access the neighboring sample or storage position.
     * @param nextAdj Identifier of the next adjacency relation in the propagation sequence.
     * @param npropagation Number of adjacency relations used during propagation.
     */
    AdaptiveAdj(std::vector<Point2D>&& offset, std::initializer_list<int> nextAdj, size_t npropagation)
        : offset_(std::move(offset)), nextAdj_{nextAdj}, npropagation_{npropagation} {}
    /**
     * @brief Constructs `AdaptiveAdj` from the supplied inputs.
     *
     * @param offset Offset used to access the neighboring sample or storage position.
     * @param nextAdj Identifier of the next adjacency relation in the propagation sequence.
     * @param npropagation Number of adjacency relations used during propagation.
     */
    AdaptiveAdj(std::initializer_list<Point2D> offset, std::vector<int>&& nextAdj, size_t npropagation)
        : offset_(offset), nextAdj_(std::move(nextAdj)), npropagation_{npropagation} {}

    /**
     * @brief Returns the complete neighbour view around `p`.
     *
     * @param p Point used by the operation.
     * @return The complete neighbour view around p.
     */
    Neighbors neighbors(const Point2D& p) const { return Neighbors(p, *this, offset_.size()); }

    /**
     * @brief Returns only the propagation prefix around `p`.
     *
     * @param p Point used by the operation.
     * @return Only the propagation prefix around p.
     */
    Neighbors neighborsPropogation(const Point2D& p) const { return Neighbors(p, *this, npropagation_); }

  private:
    /**
     * @brief Neighbour coordinate offsets.
     */
    std::vector<Point2D> offset_;

    /**
     * @brief Next adaptive-adjacency id associated with each offset.
     */
    std::vector<int> nextAdj_;

    /**
     * @brief Number of offsets used during propagation.
     */
    size_t npropagation_; // number of elements for propagation
};

/**
 * @brief Precomputed bank of adaptive adjacency stencils.
 *
 * The distance transform stores one stencil id per pixel in `adjMap_`.
 * Interior pixels use the full 8-neighbour stencil, while border and corner
 * pixels use restricted stencils that do not point outside the image. The
 * `nextAdj` values stored in each AdaptiveAdj entry update that stencil id
 * as the propagation front moves from one pixel to the next.
 */
class AdaptiveAdjBank {
  public:
    /**
     * @brief Builds the fixed stencil table used by EdtDIFT.
     */
    AdaptiveAdjBank() {
        AdaptiveAdj adj1({Point2D(1, -1), Point2D(1, 0), Point2D(1, 1), Point2D(-1, -1), Point2D(-1, 0), Point2D(-1, 1), Point2D(0, -1), Point2D(0, 1)},
                         {5, 1, 6, 7, 2, 8, 3, 4}, 8);
        bank_.push_back(adj1);

        AdaptiveAdj adj2({Point2D(1, -1), Point2D(1, 0), Point2D(1, 1)}, {5, 1, 6}, 3);
        bank_.push_back(adj2);

        AdaptiveAdj adj3({Point2D(-1, -1), Point2D(-1, 0), Point2D(-1, 1)}, {7, 2, 8}, 3);
        bank_.push_back(adj3);

        AdaptiveAdj adj4({Point2D(-1, -1), Point2D(0, -1), Point2D(1, -1)}, {7, 3, 5}, 3);
        bank_.push_back(adj4);

        AdaptiveAdj adj5({Point2D(-1, 1), Point2D(0, 1), Point2D(1, 1)}, {8, 4, 6}, 3);
        bank_.push_back(adj5);

        AdaptiveAdj adj6({Point2D(0, -1), Point2D(1, -1), Point2D(1, 0), Point2D(-1, -1), Point2D(1, 1)}, {3, 5, 1, 7, 6}, 3);
        bank_.push_back(adj6);

        AdaptiveAdj adj7({Point2D(0, 1), Point2D(1, 1), Point2D(1, 0), Point2D(-1, 1), Point2D(1, -1)}, {4, 6, 1, 8, 5}, 3);
        bank_.push_back(adj7);

        AdaptiveAdj adj8({Point2D(0, -1), Point2D(-1, -1), Point2D(-1, 0), Point2D(1, -1), Point2D(-1, 1)}, {3, 7, 2, 5, 8}, 3);
        bank_.push_back(adj8);

        AdaptiveAdj adj9({Point2D(0, 1), Point2D(-1, 1), Point2D(-1, 0), Point2D(1, 1), Point2D(-1, -1)}, {4, 8, 2, 6, 7}, 3);
        bank_.push_back(adj9);
    }

    /**
     * @brief Number of stencils stored in the bank.
     *
     * @return The number of stencils stored in the bank.
     */
    inline size_t size() const noexcept { return bank_.size(); }

    /**
     * @brief Returns a stencil by its adjacency-map id.
     *
     * @param idx Zero-based index used by the operation.
     * @return A stencil by its adjacency-map id.
     */
    const AdaptiveAdj& adj(int idx) const { return bank_[idx]; }

    /**
     * @brief Array-style access to adj(idx).
     *
     * @param idx Zero-based index used by the operation.
     * @return Reference to the resulting object.
     */
    const AdaptiveAdj& operator[](int idx) const { return adj(idx); }

  private:
    /** @brief Stores the bank. */
    std::vector<AdaptiveAdj> bank_;
};

/**
 * @brief Dynamic Image Foresting Transform for MAX_DIST support updates.
 *
 * EdtDIFT maintains a squared Euclidean distance transform on a binary
 * support that changes as tree nodes are processed by altitude level. New
 * contour pixels are inserted as zero-cost seeds, interior pixels are opened
 * with infinite cost, and `run()` propagates the best root labels through a
 * bucket priority queue.
 *
 * Distances are stored as squared integer distances (`dx*dx + dy*dy`). No
 * square root is taken. `Bedt_` accumulates, for each contour root, the
 * largest finalized squared distance reached from that root; MAX_DIST then
 * queries those root values over the active node contour.
 */
class EdtDIFT {
  public:
    /**
     * @brief Invalid pixel index sentinel.
     */
    inline static constexpr int NIL = -1;

    /**
     * @brief Allocates all buffers for an `nrows x ncols` image domain.
     *
     * The queue bucket domain is derived from a conservative squared-distance
     * bound for the image size. All pixels start outside the binary support,
     * with themselves as root representatives.
     *
     * @param nrows Number of rows in the domain.
     * @param ncols Number of columns in the domain.
     */
    EdtDIFT(int nrows, int ncols)
        : bin_{::mmcfilters::detail::CommittedImageAccess::createValue<std::uint8_t>(nrows, ncols)},
          root_{::mmcfilters::detail::CommittedImageAccess::createValue<int>(nrows, ncols)},
          Bedt_{::mmcfilters::detail::CommittedImageAccess::createValue<int>(nrows, ncols)},
          adjMap_{::mmcfilters::detail::CommittedImageAccess::createValue<std::uint8_t>(nrows, ncols)},
          O_{::mmcfilters::detail::CommittedImageAccess::createValue<std::uint8_t>(nrows, ncols)},
          Q_{detail::square(static_cast<int>(std::min(ncols, nrows) / 2.0 + 1)), nrows * ncols},
          adj4_{::mmcfilters::detail::CommittedGridAccess::radiusAdjacency(nrows, ncols, 1.0)}, domain_{ncols, nrows}, stack_(nrows * ncols) {
        bin_.fill(0);
        root_.fill(0);
        Bedt_.fill(0);
        O_.fill(0);
        adjMap_.fill(0);

        for (int pidx = 0; pidx < bin_.getSize(); ++pidx) {
            root_[pidx] = pidx;
        }

        setUpAdjMap();
    }

    /**
     * @brief Propagates all queued distance seeds/updates until convergence.
     *
     * The queue must contain the contour seeds and affected neighbours for
     * the current level before this method is called. Each popped pixel is
     * finalized for the current pass (`O_ = 0`), its root's boundary distance
     * is updated, and admissible neighbours are relaxed through the adaptive
     * stencil associated with the pixel.
     */
    void run() {
        while (!Q_.isEmpty()) {
            int pidx = Q_.popMinFIFO();
            Point2D p = domain_.point(pidx);
            O_[pidx] = 0;

            int ridx = root_[pidx];
            Point2D r = domain_.point(ridx);
            Bedt_[ridx] = std::max(Bedt_[ridx], Q_.cost(pidx));

            const AdaptiveAdj& AA = AAB_[adjMap_[pidx]];
            for (const auto& [q, ai] : AA.neighborsPropogation(p)) {
                int qidx = domain_.index(q);
                if (qidx == NIL)
                    continue;

                int dx = q.x() - r.x();
                int dy = q.y() - r.y();
                int tmp = detail::square(dx) + detail::square(dy);

                if (tmp < Q_.cost(qidx) && O_[qidx] == 1) {
                    if (Q_.state(qidx) != PQueue::State::QUEUED) {
                        Q_.setCost(qidx, tmp);
                        Q_.insert(qidx);
                    } else {
                        Q_.update(qidx, tmp);
                    }

                    root_[qidx] = ridx;
                    adjMap_[qidx] = ai;
                }
            }
        }
    }

    /**
     * @brief Marks a pixel as belonging to the current binary support.
     *
     * @param pidx Row-major pixel index used by the operation.
     */
    inline void addPixelToBinaryImage(int pidx) { bin_[pidx] = 1; }

    /**
     * @brief Requeues finite-distance support neighbours of an opened pixel.
     *
     * When a new interior pixel is opened with infinite cost, adjacent pixels
     * whose labels may propagate into it are inserted back into the queue.
     *
     * @param pidx Row-major pixel index used by the operation.
     */
    void insertNeighborsPQueue(int pidx) {
        for (int qidx : ::mmcfilters::detail::CommittedGridAccess::neighbors(adj4_, pidx)) {
            if (bin_[qidx] > 0 && Q_.cost(qidx) != PQueue::PINF && Q_.state(qidx) != PQueue::State::QUEUED) {
                Q_.setState(qidx, PQueue::State::NOT_PROCESSED);
                Q_.insert(qidx);
            }
        }
    }

    /**
     * @brief Inserts a contour pixel as a zero-cost IFT seed.
     *
     * @param pidx Row-major pixel index used by the operation.
     */
    void seed(int pidx) {
        root_[pidx] = pidx;
        Q_.setCost(pidx, 0);
        Q_.insert(pidx);
    }

    /**
     * @brief Opens an interior pixel for distance propagation.
     *
     * Open pixels (`O_ = 1`) can still receive a better distance label during
     * run(). They start from infinite cost until reached from a contour seed.
     *
     * @param pidx Row-major pixel index used by the operation.
     */
    void open(int pidx) {
        O_[pidx] = 1;
        Q_.setCost(pidx, PQueue::PINF);
    }

    /**
     * @brief Invalidates a set of contour/interior pixels after tree pruning.
     *
     * Removed pixels are reopened with infinite cost. The method then follows
     * the existing root map to reopen any pixels whose closest root was also
     * invalidated, and requeues neighbouring valid support pixels as new
     * propagation sources.
     *
     * @param toRemove Pixel identifiers invalidated by the operation.
     */
    void treeRemoval(const std::vector<int>& toRemove) {
        int top = -1;
        for (int pidx : toRemove) {
            O_[pidx] = 1;
            Q_.setCost(pidx, PQueue::PINF);
            Q_.setState(pidx, PQueue::State::NOT_PROCESSED);
            ++top;
            stack_[top] = pidx;
        }

        while (top > -1) {
            int pidx = stack_[top];
            Point2D p = domain_.point(pidx);
            --top;

            const AdaptiveAdj& AA = AAB_[adjMap_[pidx]];
            for (const auto& [q, ai] : AA.neighbors(p)) {
                int qidx = domain_.index(q);
                if (qidx == NIL)
                    continue;

                if (Q_.cost(root_[qidx]) == PQueue::PINF) {
                    if (O_[qidx] == 0) {
                        O_[qidx] = 1;
                        Q_.setCost(qidx, PQueue::PINF);
                        Q_.setState(qidx, PQueue::State::NOT_PROCESSED);
                        ++top;
                        stack_[top] = qidx;
                    }
                } else if (bin_[qidx] > 0 && Q_.state(qidx) != PQueue::State::QUEUED) {
                    Q_.setState(qidx, PQueue::State::NOT_PROCESSED);
                    Q_.insert(qidx);
                }
            }
        }
    }

    /**
     * @brief Returns the largest root distance attached to a node contour.
     *
     * `Ncontour` is expected to contain contour seed pixel ids for the active
     * node. Each seed indexes `Bedt_`, whose value is the largest squared
     * distance reached by that seed in the current binary support.
     *
     * @param Ncontour Contour pixel identifiers used by the operation.
     * @return The largest root distance attached to a node contour.
     */
    [[nodiscard]] int maxBedt(const std::vector<int>& Ncontour) const {
        int maxValue = 0;
        for (int pidx : Ncontour) {
            int d = Bedt_[pidx];
            if (d > maxValue)
                maxValue = d;
        }

        return maxValue;
    }

    /**
     * @brief Builds an uint8 visualization of the current distance labels.
     *
     * The output is normalized by the largest stored queue cost. It is meant
     * for debugging only; it is not used to compute MAX_DIST attributes.
     *
     * @return The resulting uint8 visualization of the current distance labels.
     */
    [[nodiscard]] ImageUInt8Ptr distanceTransformImage() const {
        const std::vector<int>& cost = Q_.cost();

        ImageUInt8Ptr d = ImageUInt8::create(bin_.getNumRows(), bin_.getNumCols());
        int maxCost = 0;
        for (int pidx = 0; pidx < bin_.getSize(); ++pidx) {
            if (maxCost < cost[pidx]) {
                maxCost = cost[pidx];
            }
        }

        if (maxCost == 0) {
            d->fill(0);
            return d;
        }

        for (int pidx = 0; pidx < bin_.getSize(); ++pidx) {
            (*d)[pidx] = static_cast<uint8_t>((static_cast<float>(cost[pidx]) / static_cast<float>(maxCost) * 255));
        }

        return d;
    }

    /**
     * @brief Builds an uint8 visualization of the current binary support.
     *
     * @return The resulting uint8 visualization of the current binary support.
     */
    [[nodiscard]] ImageUInt8Ptr binaryImageForVisualisation() const {
        ImageUInt8Ptr b = ImageUInt8::create(bin_.getNumRows(), bin_.getNumCols());
        for (int pidx = 0; pidx < b->getSize(); pidx++) {
            (*b)[pidx] = bin_[pidx] == 0 ? 255 : 0;
        }
        return b;
    }

    /**
     * @brief Prints the normalized distance transform for debugging.
     */
    void printDistanceTransform() const {
        detail::printUInt8Image(distanceTransformImage());
        std::cout << "\n ---- Queue Cost --- \n";
    }

    /**
     * @brief Prints the binary support visualization for debugging.
     */
    void printUnderlyingBinaryImage() const { detail::printUInt8Image(binaryImageForVisualisation()); }

    /**
     * @brief Prints the current root map as a table for small images.
     */
    void displayRootMapForSmallImages() const {
        for (int row = 0; row < root_.getNumRows(); row++) {
            std::cout << std::setw(4) << row;
            for (int col = 0; col < root_.getNumCols(); col++) {
                std::cout << std::setw(4) << root_[row * root_.getNumCols() + col];
            }
            std::cout << "\n";
        }
    }

  private:
    /**
     * @brief Initializes per-pixel adaptive-adjacency ids.
     *
     * Interior pixels keep the default full stencil (`0`). Border and corner
     * pixels receive restricted stencil ids so propagation avoids stepping
     * outside the rectangular domain.
     */
    void setUpAdjMap() {
        int lastCol = adjMap_.getNumCols() - 1;
        for (int i = 1; i < adjMap_.getNumRows(); i++) {
            adjMap_[i * adjMap_.getNumCols()] = 1;
            adjMap_[(i * adjMap_.getNumCols()) + lastCol] = 2;
        }

        int lastRow = adjMap_.getNumCols() * (adjMap_.getNumRows() - 1);
        for (int i = 1; i < adjMap_.getNumCols() - 1; i++) {
            adjMap_[i] = 4;
            adjMap_[lastRow + i] = 3;
        }

        adjMap_[0] = 6;
        adjMap_[adjMap_.getNumCols() - 1] = 8;
        adjMap_[adjMap_.getNumCols() * (adjMap_.getNumRows() - 1)] = 5;
        adjMap_[adjMap_.getNumCols() * (adjMap_.getNumRows() - 1) + adjMap_.getNumCols() - 1] = 7;
    }

  private:
    /**
     * @brief Binary support of pixels currently inserted into the level set.
     */
    ImageUInt8 bin_;

    /**
     * @brief Current nearest contour seed/root for each pixel.
     */
    ImageInt32 root_;

    /**
     * @brief Largest finalized squared distance reached from each root seed.
     */
    ImageInt32 Bedt_;

    /**
     * @brief Per-pixel index into AdaptiveAdjBank.
     */
    ImageUInt8 adjMap_;

    /**
     * @brief Open-mask: pixels marked `1` can still be relabelled this pass.
     */
    ImageUInt8 O_;

    /**
     * @brief Bucket queue keyed by squared distance labels.
     */
    PQueue Q_;

    /**
     * @brief Four-neighbour relation used to activate local propagation.
     */
    RegularGridAdjacency2D adj4_;

    /**
     * @brief Fixed table of adaptive 8-neighbour stencils.
     */
    AdaptiveAdjBank AAB_;

    /**
     * @brief Rectangular image-domain mapper between ids and coordinates.
     */
    Box2D domain_;

    /**
     * @brief Scratch stack used by treeRemoval() flood invalidation.
     */
    std::vector<int> stack_;
};
} // namespace mmcfilters::attributes::computers::detail::maxdist
