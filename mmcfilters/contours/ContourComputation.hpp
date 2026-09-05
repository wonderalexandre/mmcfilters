#pragma once

#include "../trees/ValuedMorphologicalTreeView.hpp"
#include "detail/ContourTraversal.hpp"

#include <iterator>
#include <memory>
#include <span>
#include <utility>
#include <vector>

namespace mmcfilters {

/**
 * @brief Incremental foreground A4 contours on the image domain.
 *
 * Stores compact support and boundary-lifetime indexes, without a contour cache.
 * Iteration emits every live node in post-order, reusing child storage. A
 * query for one node scans only the requested support and returns an owned pixel
 * vector.
 * The tree must outlive this object and its iterators, and remain unchanged.
 * Index construction and each traversal use O(P + N) auxiliary storage. The
 * tree resolves comparable pairs of smallest nodes from DFS intervals, then
 * selects the lower estimated storage between RMQ and offline Tarjan for the
 * remaining batch. This bound excludes caller-retained output. Pixel order is
 * unspecified.
 */
class ContourComputation {
  private:
    /** @brief Immutable indexes shared by the computation and its iterators. */
    struct SharedIndexes {
        /// Source tree whose topology defines the contours.
        const MorphologicalTree& tree;
        /// Tree mutation version captured when the indexes were built.
        std::size_t mutationVersion;
        /// Contiguous support pixels for every live node.
        contours::detail::NodeSupportIndex supportIndex;
        /// Node interval over which each pixel belongs to a contour.
        contours::detail::ContourLifetimeIndex contourLifetimes;

        /**
         * @brief Builds the indexes shared by contour traversals.
         * @param source Stable source tree.
         */
        explicit SharedIndexes(const MorphologicalTree& source)
            : tree(source), mutationVersion(source.getMutationVersion()), supportIndex(source), contourLifetimes(source) {}

        /** @brief Rejects access after a topology mutation. */
        void requireStableTree() const { tree.requireMutationVersion(mutationVersion, "ContourComputation"); }
    };

    /** @brief Mutable position shared by copies of one input iterator. */
    struct TraversalState {
        /// Shared immutable indexes and source-tree reference.
        std::shared_ptr<const SharedIndexes> indexes;
        /// Incremental contour traversal at the current position.
        contours::detail::ContourTraversal traversal;
        /// Whether dereferencing the iterator yields a contour.
        bool hasCurrentContour;

        /**
         * @brief Starts a traversal over the supplied indexes.
         * @param source Shared contour indexes.
         */
        explicit TraversalState(std::shared_ptr<const SharedIndexes> source)
            : indexes(std::move(source)), traversal(indexes->tree, indexes->supportIndex, indexes->contourLifetimes),
              hasCurrentContour(traversal.advance()) {}
    };

  public:
    /**
     * @brief Single-pass iterator yielding (node, borrowed pixel span).
     *
     * The span is valid until this iterator, or one of its copies, advances or
     * the shared traversal is destroyed. Copies share a position; each begin()
     * creates an independent traversal. Copy the pixels to retain a contour.
     */
    class iterator {
      public:
        /// C++20 iterator concept for this single-pass traversal.
        using iterator_concept = std::input_iterator_tag;
        /// Iterator category used by standard algorithms.
        using iterator_category = std::input_iterator_tag;
        /// Node identifier and borrowed contour span yielded by dereference.
        using value_type = std::pair<NodeId, std::span<const PixelId>>;
        /// Signed type used to represent iterator distances.
        using difference_type = std::ptrdiff_t;

        iterator() = default;

        /**
         * @brief Borrows the current node contour.
         * @return Current node identifier and contour span.
         */
        [[nodiscard]] value_type operator*() const {
            if (!state_ || !state_->hasCurrentContour) {
                throw std::out_of_range("Contour iterator is exhausted.");
            }
            return state_->traversal.current();
        }

        /**
         * @brief Advances the shared single-pass position.
         * @return This iterator after advancing.
         */
        iterator& operator++() {
            if (!state_ || !state_->hasCurrentContour) {
                throw std::out_of_range("Contour iterator is exhausted.");
            }
            state_->hasCurrentContour = state_->traversal.advance();
            return *this;
        }

        /** @brief Advances without retaining the previous borrowed contour. */
        void operator++(int) { ++*this; }

        /** @brief Tests exhaustion against the range's sentinel. */
        friend bool operator==(const iterator& it, std::default_sentinel_t) noexcept {
            return !it.state_ || !it.state_->hasCurrentContour;
        }

      private:
        friend class ContourComputation;

        /**
         * @brief Creates the first position of an independent traversal.
         * @param indexes Shared contour indexes.
         */
        explicit iterator(std::shared_ptr<const SharedIndexes> indexes) : state_(std::make_shared<TraversalState>(std::move(indexes))) {}
        /// Mutable traversal position shared by iterator copies.
        std::shared_ptr<TraversalState> state_;
    };

    /**
     * @brief Builds the compact indexes for a stable tree with a 2D domain.
     * @param tree Source topology, which must outlive the computation.
     */
    explicit ContourComputation(const MorphologicalTree& tree) : indexes_(std::make_shared<SharedIndexes>(tree)) {}

    /**
     * @brief Builds contours from the topology of a current valued view.
     * @param view Current valued view whose topology must outlive the computation.
     */
    template <AltitudeValue T>
    explicit ContourComputation(const ValuedMorphologicalTreeView<T>& view) : ContourComputation(currentTopology(view)) {}

    /**
     * @brief Returns an owned contour of one live node, without caching it.
     * @param node Live internal node identifier.
     * @return Foreground boundary pixels, obtained in O(|support(node)|) time.
     */
    [[nodiscard]] std::vector<PixelId> contour(NodeId node) const {
        const auto support = indexes_->supportIndex.support(node);
        std::vector<PixelId> contourPixels;
        for (PixelId pixel : support) {
            if (indexes_->contourLifetimes.establishedIsContourPixel(pixel, node)) {
                contourPixels.push_back(pixel);
            }
        }
        return contourPixels;
    }

    /**
     * @brief Starts an independent incremental post-order traversal.
     * @return Single-pass iterator positioned at the first contour, or exhausted.
     */
    [[nodiscard]] iterator begin() const {
        indexes_->requireStableTree();
        return iterator(indexes_);
    }

    /**
     * @brief Returns the exhaustion sentinel.
     * @return Sentinel shared by all contour traversals.
     */
    [[nodiscard]] std::default_sentinel_t end() const noexcept { return {}; }

    /**
     * @brief Calls consumer(node, pixels) once for every live node.
     *
     * Uses the same iterator. The borrowed span expires after the callback;
     * exceptions propagate and release the traversal's temporary storage.
     * @param consumer Callback accepting a node identifier and borrowed contour span.
     */
    template <typename Consumer> void forEachContour(Consumer&& consumer) const {
        for (auto [node, pixels] : *this) {
            consumer(node, pixels);
            indexes_->requireStableTree();
        }
    }

  private:
    /**
     * @brief Returns the stable topology behind a valued view.
     * @param view Current valued tree view.
     * @return Stable topology referenced by `view`.
     */
    template <AltitudeValue T> static const MorphologicalTree& currentTopology(const ValuedMorphologicalTreeView<T>& view) {
        MMCFILTERS_CONTRACT_CHECKED_ONLY(view.requireTopologyUnchanged("ContourComputation"));
        return view.topology();
    }

    /// Immutable indexes shared with active iterators.
    std::shared_ptr<const SharedIndexes> indexes_;
};

} // namespace mmcfilters
