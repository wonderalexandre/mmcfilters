#pragma once

#include "DynamicDistanceFieldMoments2D.hpp"
#include "SquaredDistanceBucketQueue.hpp"
#include "../../../../utils/Common.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace mmcfilters::attributes::computers::detail::distance_transform_approx {

using ApproxSquaredDistance = std::int64_t;

/**
 * @brief Adaptive A8 propagation policy used by the approximate EDT-DIFT.
 *
 * The foreground boundary remains A4; this policy controls only IFT
 * propagation. Each step carries the directional stencil selected after the
 * move, matching the adaptive adjacency proposed for the approximate MAX_DIST
 * implementation.
 */
struct AdaptiveA8SquaredEuclideanPolicy {
    struct Step {
        int rowOffset;
        int columnOffset;
        std::uint8_t nextStencil;
    };

    inline static constexpr std::array<std::array<Step, 8>, 9> stencils{{
        {{{-1, 1, 5}, {0, 1, 1}, {1, 1, 6}, {-1, -1, 7}, {0, -1, 2}, {1, -1, 8}, {-1, 0, 3}, {1, 0, 4}}},
        {{{-1, 1, 5}, {0, 1, 1}, {1, 1, 6}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}},
        {{{-1, -1, 7}, {0, -1, 2}, {1, -1, 8}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}},
        {{{-1, -1, 7}, {-1, 0, 3}, {-1, 1, 5}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}},
        {{{1, -1, 8}, {1, 0, 4}, {1, 1, 6}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}},
        {{{-1, 0, 3}, {-1, 1, 5}, {0, 1, 1}, {-1, -1, 7}, {1, 1, 6}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}},
        {{{1, 0, 4}, {1, 1, 6}, {0, 1, 1}, {1, -1, 8}, {-1, 1, 5}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}},
        {{{-1, 0, 3}, {-1, -1, 7}, {0, -1, 2}, {-1, 1, 5}, {1, -1, 8}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}},
        {{{1, 0, 4}, {1, -1, 8}, {0, -1, 2}, {1, 1, 6}, {-1, -1, 7}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}},
    }};

    inline static constexpr std::array<std::uint8_t, 9> propagationCounts{8, 3, 3, 3, 3, 3, 3, 3, 3};
    inline static constexpr std::array<std::uint8_t, 9> removalCounts{8, 3, 3, 3, 3, 5, 5, 5, 5};

    [[nodiscard]] static constexpr std::uint8_t initialStencil(int row, int column, int rows, int columns) noexcept {
        const bool top = row == 0;
        const bool bottom = row == rows - 1;
        const bool left = column == 0;
        const bool right = column == columns - 1;
        if (top && left)
            return 6;
        if (top && right)
            return 8;
        if (bottom && left)
            return 5;
        if (bottom && right)
            return 7;
        if (left)
            return 1;
        if (right)
            return 2;
        if (bottom)
            return 3;
        if (top)
            return 4;
        return 0;
    }

    [[nodiscard]] static constexpr ApproxSquaredDistance squaredCost(PixelId pixel, PixelId root, int columns) noexcept {
        return squaredCost(pixel / columns, pixel % columns, root / columns, root % columns);
    }

    [[nodiscard]] static constexpr ApproxSquaredDistance squaredCost(int row, int column, int rootRow, int rootColumn) noexcept {
        const std::int64_t rowDelta = static_cast<std::int64_t>(row) - static_cast<std::int64_t>(rootRow);
        const std::int64_t columnDelta = static_cast<std::int64_t>(column) - static_cast<std::int64_t>(rootColumn);
        return rowDelta * rowDelta + columnDelta * columnDelta;
    }
};

/**
 * @brief Initial active-edge state for a standalone or hierarchy-managed DIFT.
 */
enum class PropagationEdgeInitialization : std::uint8_t { AllDomainEdges, NoEdges };

/**
 * @brief Compile-time no-op policy for maximum-location tracking.
 */
class NoopMaximumPixelTracker2D {
  public:
    NoopMaximumPixelTracker2D(int, int) noexcept {}
    void seed(PixelId) noexcept {}
    void observe(PixelId, PixelId, ApproxSquaredDistance, ApproxSquaredDistance) noexcept {}
    [[nodiscard]] PixelId pixelForRoot(PixelId) const noexcept { return InvalidPixel; }
};

/**
 * @brief Per-Bedt-root canonical pixel attaining the recorded maximum.
 */
class MaximumPixelTracker2D {
  public:
    MaximumPixelTracker2D(int numPixels, int) : maximumPixels_(static_cast<std::size_t>(numPixels), InvalidPixel) {}

    void seed(PixelId root) { maximumPixels_[static_cast<std::size_t>(root)] = root; }

    void observe(PixelId root, PixelId pixel, ApproxSquaredDistance cost, ApproxSquaredDistance currentMaximum) {
        PixelId& currentPixel = maximumPixels_[static_cast<std::size_t>(root)];
        const distance_transform::DistanceFieldExtremum current{currentMaximum, currentPixel};
        if (distance_transform::prefersDistanceFieldSample(cost, pixel, current)) {
            currentPixel = pixel;
        }
    }

    [[nodiscard]] PixelId pixelForRoot(PixelId root) const { return maximumPixels_[static_cast<std::size_t>(root)]; }

  private:
    std::vector<PixelId> maximumPixels_;
};

/**
 * @brief Per-Bedt-root geometry of every pixel attaining its recorded maximum.
 *
 * A settled pixel may be queued again only to propagate its unchanged label to
 * newly opened support. Root epochs deduplicate those queue events while still
 * allowing the same pixel to contribute after reassignment or reseeding.
 */
template <bool ValidateInternalOperations = true> class BasicMaximumPlateauTracker2D {
  public:
    BasicMaximumPlateauTracker2D(int numPixels, int columns)
        : columns_(columns), plateaus_(static_cast<std::size_t>(numPixels)), rootEpochs_(static_cast<std::size_t>(numPixels), 0),
          observedEpochs_(static_cast<std::size_t>(numPixels), 0) {}

    void seed(PixelId root) {
        const std::size_t index = static_cast<std::size_t>(root);
        rootEpochs_[index] = ++nextEpoch_;
        plateaus_[index] = {};
    }

    void observe(PixelId root, PixelId pixel, ApproxSquaredDistance cost, ApproxSquaredDistance) {
        const std::size_t pixelIndex = static_cast<std::size_t>(pixel);
        const std::uint64_t epoch = rootEpochs_[static_cast<std::size_t>(root)];
        if constexpr (ValidateInternalOperations) {
            if (epoch == 0) {
                throw std::logic_error("Approximate EDT-DIFT plateau tracking requires a seeded Bedt root.");
            }
        }
        if (observedEpochs_[pixelIndex] == epoch) {
            return;
        }
        observedEpochs_[pixelIndex] = epoch;
        distance_transform::updateDistanceFieldMaximumPlateau(plateaus_[static_cast<std::size_t>(root)], pixel, cost, columns_);
    }

    [[nodiscard]] PixelId pixelForRoot(PixelId root) const { return plateauForRoot(root).pixel; }

    [[nodiscard]] const distance_transform::DistanceFieldMaximumPlateau& plateauForRoot(PixelId root) const {
        return plateaus_[static_cast<std::size_t>(root)];
    }

  private:
    int columns_ = 0;
    std::vector<distance_transform::DistanceFieldMaximumPlateau> plateaus_;
    std::vector<std::uint64_t> rootEpochs_;
    std::vector<std::uint64_t> observedEpochs_;
    std::uint64_t nextEpoch_ = 0;
};

using MaximumPlateauTracker2D = BasicMaximumPlateauTracker2D<true>;
using UncheckedMaximumPlateauTracker2D = BasicMaximumPlateauTracker2D<false>;

/**
 * @brief Dynamic approximate squared Euclidean distance transform.
 *
 * The binary support grows monotonically during one hierarchy sweep. Contour
 * seeds may be inserted or removed. Root, cost and per-root Bedt summaries are
 * maintained differentially. Approximation comes exclusively from the finite
 * adaptive propagation policy; this class never invokes the exact separable
 * EDT.
 *
 * Correspondence with the JMIV 2025 formulation (Algorithms 1-3), within one
 * active binary component:
 *
 * - `support_` represents membership in the paper's set `bin`;
 * - `root_`, `cost_`, and `open_` represent `root`, `cost`, and `open`;
 * - `stencil_` is the indexed adaptive-adjacency map `adjmap`;
 * - `rootMaximum_` is the per-root boundary map `Bedt`;
 * - `queue_` is the priority queue `Q`;
 * - `invalidationStack_` is the stack `T` used by `treeRemoval`.
 *
 * `activeEdges_` has no direct counterpart in the paper. It belongs to the
 * topology-driven generalization and isolates incomparable components that
 * share this global workspace until their owner LCA is processed. The C++
 * member names remain descriptive rather than duplicating paper notation.
 */
template <class PropagationPolicy = AdaptiveA8SquaredEuclideanPolicy, class DistanceFieldObserver = NoopDistanceFieldObserver2D,
          class MaximumPixelTracker = NoopMaximumPixelTracker2D, class Queue = SquaredDistanceBucketQueue, class Cost = ApproxSquaredDistance,
          bool CheckActiveEdges = true, bool AssumeOpenImpliesSupport = false, bool UseLinearRemovalOffsets = false, bool ValidateInternalOperations = true,
          bool EraseInvalidatedPixels = true>
class EdtDIFT2D {
  public:
    static_assert(std::is_integral_v<Cost> && std::is_signed_v<Cost>, "Approximate EDT-DIFT costs must be signed integers.");
    inline static constexpr Cost Infinity = std::numeric_limits<Cost>::max();

    explicit EdtDIFT2D(int rows, int columns, PropagationEdgeInitialization edgeInitialization = PropagationEdgeInitialization::AllDomainEdges)
        : rows_(rows), columns_(columns), numPixels_(checkedNumPixels(rows, columns)), support_(static_cast<std::size_t>(numPixels_), 0),
          root_(static_cast<std::size_t>(numPixels_)), cost_(static_cast<std::size_t>(numPixels_), 0), rootMaximum_(static_cast<std::size_t>(numPixels_), 0),
          stencil_(static_cast<std::size_t>(numPixels_), 0), open_(static_cast<std::size_t>(numPixels_), 0),
          activeEdges_(CheckActiveEdges ? static_cast<std::size_t>(numPixels_) : std::size_t{0}, std::uint8_t{0}),
          queue_(maximumQueueCost(rows, columns), numPixels_), observer_(numPixels_, columns_), maximumPixelTracker_(numPixels_, columns_) {
        for (PixelId pixel = 0; pixel < numPixels_; ++pixel) {
            root_[static_cast<std::size_t>(pixel)] = pixel;
            stencil_[static_cast<std::size_t>(pixel)] = PropagationPolicy::initialStencil(pixel / columns_, pixel % columns_, rows_, columns_);
        }
        if constexpr (CheckActiveEdges) {
            if (edgeInitialization == PropagationEdgeInitialization::AllDomainEdges) {
                activateAllDomainEdges();
            }
        } else if (edgeInitialization != PropagationEdgeInitialization::AllDomainEdges) {
            throw std::invalid_argument("Approximate EDT-DIFT cannot elide edge checks when propagation edges are activated incrementally.");
        }
    }

    EdtDIFT2D(const EdtDIFT2D&) = delete;
    EdtDIFT2D& operator=(const EdtDIFT2D&) = delete;
    EdtDIFT2D(EdtDIFT2D&&) = delete;
    EdtDIFT2D& operator=(EdtDIFT2D&&) = delete;

    [[nodiscard]] int rows() const noexcept { return rows_; }
    [[nodiscard]] int columns() const noexcept { return columns_; }
    [[nodiscard]] int numPixels() const noexcept { return numPixels_; }

    void addPixelToBinaryImage(PixelId pixel) {
        validatePixel(pixel);
        const std::size_t pixelIndex = static_cast<std::size_t>(pixel);
        if constexpr (ValidateInternalOperations) {
            if (support_[pixelIndex] != 0) {
                throw std::logic_error("Approximate EDT-DIFT received a duplicate support pixel insertion.");
            }
        }
        support_[pixelIndex] = 1;
        cost_[pixelIndex] = Infinity;
        if constexpr (!std::is_same_v<DistanceFieldObserver, NoopDistanceFieldObserver2D>) {
            observer_.insertPixel(pixel);
            mergeActiveSupportNeighbours(pixel);
        }
    }

    [[nodiscard]] bool isInBinaryImage(PixelId pixel) const {
        validatePixel(pixel);
        return support_[static_cast<std::size_t>(pixel)] != 0;
    }

    void seed(PixelId pixel) {
        validateSupportPixel(pixel);
        root_[static_cast<std::size_t>(pixel)] = pixel;
        rootMaximum_[static_cast<std::size_t>(pixel)] = 0;
        maximumPixelTracker_.seed(pixel);
        open_[static_cast<std::size_t>(pixel)] = 0;
        setLabelAndQueue(pixel, 0);
    }

    void open(PixelId pixel) {
        validateSupportPixel(pixel);
        open_[static_cast<std::size_t>(pixel)] = 1;
        queue_.erase(pixel);
        assignCost(pixel, Infinity);
    }

    /**
     * @brief Activates one undirected A8 edge between already inserted support pixels.
     */
    void activateEdge(PixelId first, PixelId second) {
        if constexpr (!CheckActiveEdges) {
            throw std::logic_error("Approximate EDT-DIFT with elided edge checks cannot activate individual edges.");
        }
        const int firstRow = first / columns_;
        const int firstColumn = first % columns_;
        const int secondRow = second / columns_;
        const int secondColumn = second % columns_;
        const int rowOffset = secondRow - firstRow;
        const int columnOffset = secondColumn - firstColumn;
        const std::uint8_t forwardBit = validatedDirectionBit(rowOffset, columnOffset);
        const std::uint8_t reverseBit = validatedDirectionBit(-rowOffset, -columnOffset);
        activateEdge(first, second, forwardBit, reverseBit);
    }

    /**
     * @brief Activates a trusted edge whose reciprocal direction bits were precomputed by the topology index.
     */
    void activateEdge(PixelId first, PixelId second, std::uint8_t forwardBit, std::uint8_t reverseBit) {
        if constexpr (!CheckActiveEdges) {
            throw std::logic_error("Approximate EDT-DIFT with elided edge checks cannot activate individual edges.");
        }
        validateSupportPixel(first);
        validateSupportPixel(second);
        if constexpr (ValidateInternalOperations) {
            if (forwardBit == 0 || reverseBit == 0) {
                throw std::logic_error("Approximate EDT-DIFT received an invalid precomputed edge direction.");
            }
        }
        std::uint8_t& firstEdges = activeEdges_[static_cast<std::size_t>(first)];
        if constexpr (ValidateInternalOperations) {
            if ((firstEdges & forwardBit) != 0) {
                return;
            }
        }
        firstEdges = static_cast<std::uint8_t>(firstEdges | forwardBit);
        activeEdges_[static_cast<std::size_t>(second)] = static_cast<std::uint8_t>(activeEdges_[static_cast<std::size_t>(second)] | reverseBit);
        if constexpr (!std::is_same_v<DistanceFieldObserver, NoopDistanceFieldObserver2D>) {
            observer_.mergePixels(first, second);
        }

        if (cost_[static_cast<std::size_t>(first)] != Infinity && open_[static_cast<std::size_t>(second)] != 0) {
            enqueueCurrentLabel(first);
        }
        if (cost_[static_cast<std::size_t>(second)] != Infinity && open_[static_cast<std::size_t>(first)] != 0) {
            enqueueCurrentLabel(second);
        }
    }

    /**
     * @brief Performs the A4 frontier insertion from Algorithm 3, line 27.
     */
    void insertSupportNeighbours(PixelId pixel) {
        validateSupportPixel(pixel);
        const int row = pixel / columns_;
        const int column = pixel % columns_;
        constexpr std::array<std::pair<int, int>, 4> offsets{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};
        for (const auto [rowOffset, columnOffset] : offsets) {
            const PixelId neighbour = index(row + rowOffset, column + columnOffset);
            if (neighbour >= 0 && support_[static_cast<std::size_t>(neighbour)] != 0 && cost_[static_cast<std::size_t>(neighbour)] != Infinity) {
                enqueueCurrentLabel(neighbour);
            }
        }
    }

    /**
     * @brief Implements `treeRemoval` from Algorithm 2 for removed contour seeds.
     */
    void removeSeeds(std::span<const PixelId> removals) {
        invalidationStack_.clear();
        invalidationStack_.reserve(removals.size());
        for (PixelId pixel : removals) {
            validateSupportPixel(pixel);
            if (cost_[static_cast<std::size_t>(pixel)] == Infinity) {
                continue;
            }
            invalidatePixel(pixel);
            invalidationStack_.push_back(pixel);
        }

        while (!invalidationStack_.empty()) {
            const PixelId pixel = invalidationStack_.back();
            invalidationStack_.pop_back();
            forEachRemovalNeighbour(pixel, [&](PixelId neighbour, std::uint8_t, int, int) {
                if (support_[static_cast<std::size_t>(neighbour)] == 0) {
                    return;
                }
                const PixelId neighbourRoot = root_[static_cast<std::size_t>(neighbour)];
                if (cost_[static_cast<std::size_t>(neighbourRoot)] == Infinity) {
                    if (open_[static_cast<std::size_t>(neighbour)] == 0) {
                        invalidatePixel(neighbour);
                        invalidationStack_.push_back(neighbour);
                    }
                } else if (cost_[static_cast<std::size_t>(neighbour)] != Infinity) {
                    enqueueCurrentLabel(neighbour);
                }
            });
        }
    }

    /**
     * @brief Implements `EDTDiff` from Algorithm 1 up to the current fixed point.
     */
    void run() {
        while (!queue_.empty()) {
            const PixelId pixel = queue_.popMinimumFifo();
            open_[static_cast<std::size_t>(pixel)] = 0;
            const PixelId root = root_[static_cast<std::size_t>(pixel)];
            const int rootRow = root / columns_;
            const int rootColumn = root % columns_;
            maximumPixelTracker_.observe(root, pixel, cost_[static_cast<std::size_t>(pixel)], rootMaximum_[static_cast<std::size_t>(root)]);
            rootMaximum_[static_cast<std::size_t>(root)] = std::max(rootMaximum_[static_cast<std::size_t>(root)], cost_[static_cast<std::size_t>(pixel)]);

            forEachNeighbour(pixel, true, [&](PixelId neighbour, std::uint8_t nextStencil, int neighbourRow, int neighbourColumn) {
                if constexpr (AssumeOpenImpliesSupport) {
                    if (open_[static_cast<std::size_t>(neighbour)] == 0) {
                        return;
                    }
                } else {
                    if (support_[static_cast<std::size_t>(neighbour)] == 0 || open_[static_cast<std::size_t>(neighbour)] == 0) {
                        return;
                    }
                }
                const Cost candidate = squaredCost(neighbourRow, neighbourColumn, rootRow, rootColumn);
                if (candidate < cost_[static_cast<std::size_t>(neighbour)]) {
                    root_[static_cast<std::size_t>(neighbour)] = root;
                    stencil_[static_cast<std::size_t>(neighbour)] = nextStencil;
                    setLabelAndQueue(neighbour, candidate);
                }
            });
        }
    }

    /**
     * @brief Implements the `max Bedt` contour reduction from Algorithm 3, lines 31-32.
     */
    [[nodiscard]] ApproxSquaredDistance maximumRootDistance(std::span<const PixelId> contour) const {
        ApproxSquaredDistance maximum = 0;
        for (PixelId pixel : contour) {
            validateSupportPixel(pixel);
            maximum = std::max(maximum, static_cast<ApproxSquaredDistance>(rootMaximum_[static_cast<std::size_t>(pixel)]));
        }
        return maximum;
    }

    /**
     * @brief Reduces Bedt to its maximum and canonical attaining support pixel.
     */
    [[nodiscard]] distance_transform::DistanceFieldExtremum maximumRootExtremum(std::span<const PixelId> contour) const {
        static_assert(!std::is_same_v<MaximumPixelTracker, NoopMaximumPixelTracker2D>, "Maximum-root localization requires MaximumPixelTracker2D.");
        distance_transform::DistanceFieldExtremum extremum;
        for (PixelId root : contour) {
            validateSupportPixel(root);
            const PixelId pixel = maximumPixelTracker_.pixelForRoot(root);
            if constexpr (ValidateInternalOperations) {
                if (pixel == InvalidPixel) {
                    throw std::logic_error("Approximate EDT-DIFT found a Bedt maximum without an attaining pixel.");
                }
            }
            distance_transform::updateDistanceFieldExtremum(extremum, pixel, rootMaximum_[static_cast<std::size_t>(root)]);
        }
        return extremum;
    }

    /**
     * @brief Reduces Bedt to the complete geometry of its maximum plateau.
     */
    [[nodiscard]] distance_transform::DistanceFieldMaximumPlateau maximumRootPlateau(std::span<const PixelId> contour) const {
        static_assert(
            requires(const MaximumPixelTracker& tracker, PixelId root) { tracker.plateauForRoot(root); },
            "Maximum-root plateau geometry requires a plateau tracker.");
        distance_transform::DistanceFieldMaximumPlateau plateau;
        for (PixelId root : contour) {
            validateSupportPixel(root);
            const auto& rootPlateau = maximumPixelTracker_.plateauForRoot(root);
            if constexpr (ValidateInternalOperations) {
                if (rootPlateau.pixel == InvalidPixel || rootPlateau.count == 0) {
                    throw std::logic_error("Approximate EDT-DIFT found a Bedt maximum without a plateau summary.");
                }
                if (rootPlateau.squaredDistance != rootMaximum_[static_cast<std::size_t>(root)]) {
                    throw std::logic_error("Approximate EDT-DIFT plateau summary disagrees with its Bedt maximum.");
                }
            }
            distance_transform::mergeDistanceFieldMaximumPlateau(plateau, rootPlateau);
        }
        return plateau;
    }

    [[nodiscard]] ApproxSquaredDistance squaredDistance(PixelId pixel) const {
        validateSupportPixel(pixel);
        return static_cast<ApproxSquaredDistance>(cost_[static_cast<std::size_t>(pixel)]);
    }

    [[nodiscard]] PixelId assignedRoot(PixelId pixel) const {
        validateSupportPixel(pixel);
        return root_[static_cast<std::size_t>(pixel)];
    }

    /**
     * @brief Returns the compile-time-selected observer that summarizes active DIFT components.
     */
    [[nodiscard]] const DistanceFieldObserver& observer() const noexcept { return observer_; }

  private:
    [[nodiscard]] static int checkedNumPixels(int rows, int columns) {
        if (rows <= 0 || columns <= 0 || columns > std::numeric_limits<int>::max() / rows) {
            throw std::invalid_argument("Approximate EDT-DIFT requires a consistent non-empty 2D domain.");
        }
        return rows * columns;
    }

    void requirePixel(PixelId pixel) const {
        if (pixel < 0 || pixel >= numPixels_) {
            throw std::out_of_range("Approximate EDT-DIFT received an invalid pixel id.");
        }
    }

    void requireSupportPixel(PixelId pixel) const {
        requirePixel(pixel);
        if (support_[static_cast<std::size_t>(pixel)] == 0) {
            throw std::logic_error("Approximate EDT-DIFT operation requires a pixel in the active binary support.");
        }
    }

    void validatePixel(PixelId pixel) const {
        if constexpr (ValidateInternalOperations) {
            requirePixel(pixel);
        }
    }

    void validateSupportPixel(PixelId pixel) const {
        if constexpr (ValidateInternalOperations) {
            requireSupportPixel(pixel);
        }
    }

    [[nodiscard]] PixelId index(int row, int column) const noexcept {
        if (row < 0 || row >= rows_ || column < 0 || column >= columns_)
            return InvalidPixel;
        return row * columns_ + column;
    }

    template <class Consumer> void forEachNeighbour(PixelId pixel, bool propagationOnly, Consumer&& consumer) const {
        const int row = pixel / columns_;
        const int column = pixel % columns_;
        const std::uint8_t stencil = stencil_[static_cast<std::size_t>(pixel)];
        const std::uint8_t count = propagationOnly ? PropagationPolicy::propagationCounts[stencil] : PropagationPolicy::removalCounts[stencil];
        for (std::uint8_t indexValue = 0; indexValue < count; ++indexValue) {
            const auto& step = PropagationPolicy::stencils[stencil][indexValue];
            if constexpr (CheckActiveEdges) {
                if ((activeEdges_[static_cast<std::size_t>(pixel)] & uncheckedDirectionBit(step.rowOffset, step.columnOffset)) == 0) {
                    continue;
                }
            }
            const int neighbourRow = row + step.rowOffset;
            const int neighbourColumn = column + step.columnOffset;
            const PixelId neighbour = neighbourRow * columns_ + neighbourColumn;
            consumer(neighbour, step.nextStencil, neighbourRow, neighbourColumn);
        }
    }

    template <class Consumer> void forEachRemovalNeighbour(PixelId pixel, Consumer&& consumer) const {
        if constexpr (!UseLinearRemovalOffsets) {
            forEachNeighbour(pixel, false, std::forward<Consumer>(consumer));
            return;
        }

        const std::uint8_t stencil = stencil_[static_cast<std::size_t>(pixel)];
        const std::uint8_t count = PropagationPolicy::removalCounts[stencil];
        for (std::uint8_t indexValue = 0; indexValue < count; ++indexValue) {
            const auto& step = PropagationPolicy::stencils[stencil][indexValue];
            if constexpr (CheckActiveEdges) {
                if ((activeEdges_[static_cast<std::size_t>(pixel)] & uncheckedDirectionBit(step.rowOffset, step.columnOffset)) == 0) {
                    continue;
                }
            }
            const PixelId neighbour = pixel + step.rowOffset * columns_ + step.columnOffset;
            consumer(neighbour, step.nextStencil, 0, 0);
        }
    }

    void setLabelAndQueue(PixelId pixel, Cost cost) {
        assignCost(pixel, cost);
        queue_.update(pixel, cost);
    }

    void assignCost(PixelId pixel, Cost cost) {
        const std::size_t indexValue = static_cast<std::size_t>(pixel);
        if constexpr (std::is_same_v<DistanceFieldObserver, NoopDistanceFieldObserver2D>) {
            cost_[indexValue] = cost;
            return;
        }
        const Cost previous = cost_[indexValue];
        if (previous == cost) {
            return;
        }
        if (previous != Infinity) {
            observer_.removeFiniteCost(pixel, previous);
        }
        cost_[indexValue] = cost;
        if (cost != Infinity) {
            observer_.addFiniteCost(pixel, cost);
        }
    }

    void enqueueCurrentLabel(PixelId pixel) {
        const std::size_t indexValue = static_cast<std::size_t>(pixel);
        if (queue_.contains(pixel) || cost_[indexValue] == Infinity)
            return;
        queue_.insert(pixel, cost_[indexValue]);
    }

    void invalidatePixel(PixelId pixel) {
        if constexpr (ValidateInternalOperations) {
            requireSupportPixel(pixel);
            if constexpr (!EraseInvalidatedPixels) {
                if (queue_.contains(pixel)) {
                    throw std::logic_error("Approximate EDT-DIFT removal invariant forbids invalidating a queued pixel.");
                }
            }
        }
        open_[static_cast<std::size_t>(pixel)] = 1;
        if constexpr (EraseInvalidatedPixels) {
            queue_.erase(pixel);
        }
        assignCost(pixel, Infinity);
    }

    [[nodiscard]] static ApproxSquaredDistance maximumQueueCost(int rows, int columns) {
        const ApproxSquaredDistance radius = static_cast<ApproxSquaredDistance>(std::min(rows, columns) / 2 + 1);
        return radius * radius;
    }

    [[nodiscard]] static Cost squaredCost(int row, int column, int rootRow, int rootColumn) noexcept {
        if constexpr (std::is_same_v<Cost, ApproxSquaredDistance>) {
            return PropagationPolicy::squaredCost(row, column, rootRow, rootColumn);
        }
        const Cost rowDelta = static_cast<Cost>(row - rootRow);
        const Cost columnDelta = static_cast<Cost>(column - rootColumn);
        return static_cast<Cost>(rowDelta * rowDelta + columnDelta * columnDelta);
    }

    [[nodiscard]] static constexpr std::uint8_t directionBit(int rowOffset, int columnOffset) {
        if (rowOffset < -1 || rowOffset > 1 || columnOffset < -1 || columnOffset > 1 || (rowOffset == 0 && columnOffset == 0)) {
            throw std::logic_error("Approximate EDT-DIFT requires an A8 propagation edge.");
        }
        return uncheckedDirectionBit(rowOffset, columnOffset);
    }

    [[nodiscard]] static constexpr std::uint8_t validatedDirectionBit(int rowOffset, int columnOffset) {
        if constexpr (ValidateInternalOperations) {
            return directionBit(rowOffset, columnOffset);
        }
        return uncheckedDirectionBit(rowOffset, columnOffset);
    }

    [[nodiscard]] static constexpr std::uint8_t uncheckedDirectionBit(int rowOffset, int columnOffset) noexcept {
        constexpr std::array<std::uint8_t, 9> bits{std::uint8_t{1},  std::uint8_t{2},  std::uint8_t{4},  std::uint8_t{8},  std::uint8_t{0},
                                                   std::uint8_t{16}, std::uint8_t{32}, std::uint8_t{64}, std::uint8_t{128}};
        return bits[static_cast<std::size_t>((rowOffset + 1) * 3 + (columnOffset + 1))];
    }

    void activateAllDomainEdges() {
        constexpr std::array<std::pair<int, int>, 4> forwardOffsets{{{0, 1}, {1, -1}, {1, 0}, {1, 1}}};
        for (PixelId pixel = 0; pixel < numPixels_; ++pixel) {
            const int row = pixel / columns_;
            const int column = pixel % columns_;
            for (const auto [rowOffset, columnOffset] : forwardOffsets) {
                const PixelId neighbour = index(row + rowOffset, column + columnOffset);
                if (neighbour == InvalidPixel) {
                    continue;
                }
                activeEdges_[static_cast<std::size_t>(pixel)] =
                    static_cast<std::uint8_t>(activeEdges_[static_cast<std::size_t>(pixel)] | directionBit(rowOffset, columnOffset));
                activeEdges_[static_cast<std::size_t>(neighbour)] =
                    static_cast<std::uint8_t>(activeEdges_[static_cast<std::size_t>(neighbour)] | directionBit(-rowOffset, -columnOffset));
            }
        }
    }

    void mergeActiveSupportNeighbours(PixelId pixel) {
        const int row = pixel / columns_;
        const int column = pixel % columns_;
        for (int rowOffset = -1; rowOffset <= 1; ++rowOffset) {
            for (int columnOffset = -1; columnOffset <= 1; ++columnOffset) {
                if (rowOffset == 0 && columnOffset == 0) {
                    continue;
                }
                if constexpr (CheckActiveEdges) {
                    if ((activeEdges_[static_cast<std::size_t>(pixel)] & uncheckedDirectionBit(rowOffset, columnOffset)) == 0) {
                        continue;
                    }
                }
                const PixelId neighbour = index(row + rowOffset, column + columnOffset);
                if (neighbour != InvalidPixel && support_[static_cast<std::size_t>(neighbour)] != 0) {
                    observer_.mergePixels(pixel, neighbour);
                }
            }
        }
    }

    int rows_ = 0;
    int columns_ = 0;
    int numPixels_ = 0;
    std::vector<std::uint8_t> support_;
    std::vector<PixelId> root_;
    std::vector<Cost> cost_;
    std::vector<Cost> rootMaximum_;
    std::vector<std::uint8_t> stencil_;
    std::vector<std::uint8_t> open_;
    std::vector<std::uint8_t> activeEdges_;
    Queue queue_;
    std::vector<PixelId> invalidationStack_;
    DistanceFieldObserver observer_;
    MaximumPixelTracker maximumPixelTracker_;
};

} // namespace mmcfilters::attributes::computers::detail::distance_transform_approx
