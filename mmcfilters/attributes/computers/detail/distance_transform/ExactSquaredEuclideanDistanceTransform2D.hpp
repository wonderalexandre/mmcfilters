#pragma once

#include "NodeDistanceFieldProvider.hpp"
#include "../../../../utils/Common.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace mmcfilters::attributes::computers::detail::distance_transform {

/**
 * @brief Reusable exact squared Euclidean distance-transform workspace.
 *
 * The transform is separable. It first computes the one-dimensional lower
 * envelope of squared-distance parabolas down every column, then applies the
 * same transform across every row. Boundary pixels are zero-cost sites and no
 * graph-path or tree-kind approximation is involved. After a full site
 * assignment, exact insertions/removals reuse unchanged vertical and horizontal
 * lines.
 */
class ExactSquaredEuclideanDistanceTransform2D {
  public:
    ExactSquaredEuclideanDistanceTransform2D(const ExactSquaredEuclideanDistanceTransform2D&) = delete;
    ExactSquaredEuclideanDistanceTransform2D& operator=(const ExactSquaredEuclideanDistanceTransform2D&) = delete;
    ExactSquaredEuclideanDistanceTransform2D(ExactSquaredEuclideanDistanceTransform2D&&) = delete;
    ExactSquaredEuclideanDistanceTransform2D& operator=(ExactSquaredEuclideanDistanceTransform2D&&) = delete;

    /**
     * @brief Allocates one workspace for a non-empty row-major 2D domain.
     *
     * @param rows Number of rows.
     * @param columns Number of columns.
     */
    ExactSquaredEuclideanDistanceTransform2D(int rows, int columns) { resetDomain(rows, columns); }

    /**
     * @brief Selects a new non-empty active domain while retaining allocations.
     *
     * The previous result is invalidated. Vector capacities are reused, so one
     * workspace can evaluate differently sized translated node boxes without
     * allocating a fresh transform for every node.
     */
    void resetDomain(int rows, int columns) {
        const int numPixels = validatedPixelCount(rows, columns);
        const std::size_t pixelCount = static_cast<std::size_t>(numPixels);
        const std::size_t lineLength = static_cast<std::size_t>(std::max(rows, columns));

        // Reserve every required capacity before changing any logical size.
        // If an allocation fails, the previous dimensions, result, and vector
        // sizes remain usable; successful earlier reserves only add capacity.
        seedMask_.reserve(pixelCount);
        verticalCost_.reserve(pixelCount);
        squaredDistance_.reserve(pixelCount);
        pendingDelta_.reserve(pixelCount);
        lineInput_.reserve(lineLength);
        lineOutput_.reserve(lineLength);
        envelopeSites_.reserve(lineLength);
        envelopeStarts_.reserve(lineLength);
        dirtyColumnMark_.reserve(static_cast<std::size_t>(columns));
        dirtyRowMark_.reserve(static_cast<std::size_t>(rows));
        touchedDeltaPixels_.reserve(pixelCount);
        dirtyColumns_.reserve(static_cast<std::size_t>(columns));
        dirtyRows_.reserve(static_cast<std::size_t>(rows));

        seedMask_.resize(pixelCount);
        verticalCost_.resize(pixelCount);
        squaredDistance_.resize(pixelCount);
        pendingDelta_.resize(pixelCount);
        lineInput_.resize(lineLength);
        lineOutput_.resize(lineLength);
        envelopeSites_.resize(lineLength);
        envelopeStarts_.resize(lineLength);
        dirtyColumnMark_.resize(static_cast<std::size_t>(columns));
        dirtyRowMark_.resize(static_cast<std::size_t>(rows));

        rows_ = rows;
        columns_ = columns;
        numPixels_ = numPixels;
        hasResult_ = false;
        numActiveSites_ = 0;
        std::fill(pendingDelta_.begin(), pendingDelta_.end(), std::int8_t{0});
        std::fill(dirtyColumnMark_.begin(), dirtyColumnMark_.end(), std::uint8_t{0});
        std::fill(dirtyRowMark_.begin(), dirtyRowMark_.end(), std::uint8_t{0});
        touchedDeltaPixels_.clear();
        dirtyColumns_.clear();
        dirtyRows_.clear();
    }

    /**
     * @brief Returns the active-domain row count.
     */
    [[nodiscard]] int rows() const noexcept { return rows_; }

    /**
     * @brief Returns the active-domain column count.
     */
    [[nodiscard]] int columns() const noexcept { return columns_; }

    /**
     * @brief Computes distances to the supplied non-empty set of boundary sites.
     *
     * Results remain valid until the next call to `compute`, `applySiteDelta`,
     * `resetDomain`, or destruction of this workspace.
     *
     * @param boundaryPixels Row-major boundary-site identifiers.
     */
    void compute(std::span<const PixelId> boundaryPixels) {
        if (boundaryPixels.empty()) {
            throw std::invalid_argument("Exact squared Euclidean distance transform requires at least one boundary pixel.");
        }

        for (PixelId pixel : boundaryPixels) {
            if (pixel < 0 || pixel >= numPixels_) {
                throw std::out_of_range("Exact squared Euclidean distance transform received an invalid boundary pixel.");
            }
        }

        computeEstablished(boundaryPixels);
    }

    /**
     * @brief Computes a field from caller-established, non-empty boundary sites.
     */
    void computeEstablished(std::span<const PixelId> boundaryPixels) {
        hasResult_ = false;
        std::fill(seedMask_.begin(), seedMask_.end(), std::uint8_t{0});
        for (PixelId pixel : boundaryPixels) {
            seedMask_[static_cast<std::size_t>(pixel)] = 1;
        }
        numActiveSites_ = static_cast<std::size_t>(std::count(seedMask_.begin(), seedMask_.end(), std::uint8_t{1}));

        for (int column = 0; column < columns_; ++column) {
            for (int row = 0; row < rows_; ++row) {
                const PixelId pixel = row * columns_ + column;
                lineInput_[static_cast<std::size_t>(row)] = seedMask_[static_cast<std::size_t>(pixel)] != 0 ? SquaredDistance{0} : infinity();
            }
            transformLine(rows_);
            for (int row = 0; row < rows_; ++row) {
                verticalCost_[static_cast<std::size_t>(row * columns_ + column)] = lineOutput_[static_cast<std::size_t>(row)];
            }
        }

        for (int row = 0; row < rows_; ++row) {
            for (int column = 0; column < columns_; ++column) {
                lineInput_[static_cast<std::size_t>(column)] = verticalCost_[static_cast<std::size_t>(row * columns_ + column)];
            }
            transformLine(columns_);
            for (int column = 0; column < columns_; ++column) {
                const PixelId pixel = row * columns_ + column;
                squaredDistance_[static_cast<std::size_t>(pixel)] = lineOutput_[static_cast<std::size_t>(column)];
            }
        }
        hasResult_ = true;
    }

    /**
     * @brief Applies exact boundary-site insertions/removals on the active domain.
     *
     * Only changed columns are transformed vertically. A horizontal row is
     * transformed exactly when at least one of its vertical costs changed.
     * Additions followed by removals of the same previously absent site cancel
     * within the batch.
     */
    void applySiteDelta(std::span<const PixelId> additions, std::span<const PixelId> removals) { applySiteDeltaImpl<true>(additions, removals); }

    /**
     * @brief Applies caller-established insertions and removals to an existing field.
     */
    void applyEstablishedSiteDelta(std::span<const PixelId> additions, std::span<const PixelId> removals) { applySiteDeltaImpl<false>(additions, removals); }

    template <bool ValidateInternalOperations> void applySiteDeltaImpl(std::span<const PixelId> additions, std::span<const PixelId> removals) {
        if constexpr (ValidateInternalOperations) {
            if (!hasResult_) {
                throw std::logic_error("Incremental exact distance transform requires an existing computed result.");
            }
        }

        touchedDeltaPixels_.clear();
        touchedDeltaPixels_.reserve(additions.size() + removals.size());
        try {
            for (PixelId pixel : additions) {
                if constexpr (ValidateInternalOperations) {
                    requireValidPixel(pixel, "addition");
                }
                std::int8_t& delta = pendingDelta_[static_cast<std::size_t>(pixel)];
                if constexpr (ValidateInternalOperations) {
                    if (seedMask_[static_cast<std::size_t>(pixel)] != 0 || delta != 0) {
                        throw std::invalid_argument("Incremental exact distance transform received a duplicate or already-active site addition.");
                    }
                }
                delta = 1;
                touchedDeltaPixels_.push_back(pixel);
            }
            for (PixelId pixel : removals) {
                if constexpr (ValidateInternalOperations) {
                    requireValidPixel(pixel, "removal");
                }
                std::int8_t& delta = pendingDelta_[static_cast<std::size_t>(pixel)];
                if (delta == 1) {
                    delta = 0;
                } else if constexpr (ValidateInternalOperations) {
                    if (delta == 0 && seedMask_[static_cast<std::size_t>(pixel)] != 0) {
                        delta = -1;
                        touchedDeltaPixels_.push_back(pixel);
                    } else {
                        throw std::invalid_argument("Incremental exact distance transform received a duplicate or inactive site removal.");
                    }
                } else {
                    delta = -1;
                    touchedDeltaPixels_.push_back(pixel);
                }
            }
        } catch (...) {
            clearPendingDeltas();
            throw;
        }

        std::size_t effectiveInsertions = 0;
        std::size_t effectiveRemovals = 0;
        for (PixelId pixel : touchedDeltaPixels_) {
            const std::int8_t delta = pendingDelta_[static_cast<std::size_t>(pixel)];
            effectiveInsertions += delta > 0 ? std::size_t{1} : std::size_t{0};
            effectiveRemovals += delta < 0 ? std::size_t{1} : std::size_t{0};
        }
        if constexpr (ValidateInternalOperations) {
            if (effectiveRemovals >= numActiveSites_ + effectiveInsertions) {
                clearPendingDeltas();
                throw std::invalid_argument("Incremental exact distance transform requires at least one active boundary site.");
            }
        }

        dirtyColumns_.clear();
        for (PixelId pixel : touchedDeltaPixels_) {
            const std::int8_t delta = pendingDelta_[static_cast<std::size_t>(pixel)];
            if (delta == 0) {
                continue;
            }
            seedMask_[static_cast<std::size_t>(pixel)] = delta > 0 ? std::uint8_t{1} : std::uint8_t{0};
            markDirtyColumn(pixel % columns_);
        }
        clearPendingDeltas();
        numActiveSites_ += effectiveInsertions;
        numActiveSites_ -= effectiveRemovals;

        dirtyRows_.clear();
        for (int column : dirtyColumns_) {
            for (int row = 0; row < rows_; ++row) {
                const PixelId pixel = row * columns_ + column;
                lineInput_[static_cast<std::size_t>(row)] = seedMask_[static_cast<std::size_t>(pixel)] != 0 ? SquaredDistance{0} : infinity();
            }
            transformLine(rows_);
            for (int row = 0; row < rows_; ++row) {
                const std::size_t pixelIndex = static_cast<std::size_t>(row * columns_ + column);
                const SquaredDistance newCost = lineOutput_[static_cast<std::size_t>(row)];
                if (newCost != verticalCost_[pixelIndex]) {
                    verticalCost_[pixelIndex] = newCost;
                    markDirtyRow(row);
                }
            }
            dirtyColumnMark_[static_cast<std::size_t>(column)] = 0;
        }
        dirtyColumns_.clear();

        for (int row : dirtyRows_) {
            for (int column = 0; column < columns_; ++column) {
                lineInput_[static_cast<std::size_t>(column)] = verticalCost_[static_cast<std::size_t>(row * columns_ + column)];
            }
            transformLine(columns_);
            for (int column = 0; column < columns_; ++column) {
                const PixelId pixel = row * columns_ + column;
                squaredDistance_[static_cast<std::size_t>(pixel)] = lineOutput_[static_cast<std::size_t>(column)];
            }
            dirtyRowMark_[static_cast<std::size_t>(row)] = 0;
        }
        dirtyRows_.clear();
    }

    /**
     * @brief Returns the exact squared distance of one valid domain pixel.
     */
    [[nodiscard]] SquaredDistance squaredDistance(PixelId pixel) const {
        if (!hasResult_) {
            throw std::logic_error("Exact squared Euclidean distance transform has no computed result.");
        }
        if (pixel < 0 || pixel >= numPixels_) {
            throw std::out_of_range("Exact squared Euclidean distance transform received an invalid query pixel.");
        }
        return squaredDistance_[static_cast<std::size_t>(pixel)];
    }

    /**
     * @brief Returns one caller-established computed squared-distance sample.
     */
    [[nodiscard]] SquaredDistance establishedSquaredDistance(PixelId pixel) const noexcept { return squaredDistance_[static_cast<std::size_t>(pixel)]; }

  private:
    [[nodiscard]] static constexpr SquaredDistance infinity() noexcept { return std::numeric_limits<SquaredDistance>::max(); }

    [[nodiscard]] static int validatedPixelCount(int rows, int columns) {
        if (rows <= 0 || columns <= 0) {
            throw std::invalid_argument("Exact squared Euclidean distance transform requires a non-empty 2D domain.");
        }
        if (columns > std::numeric_limits<int>::max() / rows) {
            throw std::overflow_error("Exact squared Euclidean distance-transform pixel count exceeds the supported integer domain.");
        }
        return rows * columns;
    }

    void requireValidPixel(PixelId pixel, const char* operation) const {
        if (pixel < 0 || pixel >= numPixels_) {
            throw std::out_of_range(std::string("Incremental exact distance-transform ") + operation + " received an invalid pixel.");
        }
    }

    void clearPendingDeltas() noexcept {
        for (PixelId pixel : touchedDeltaPixels_) {
            pendingDelta_[static_cast<std::size_t>(pixel)] = 0;
        }
        touchedDeltaPixels_.clear();
    }

    void markDirtyColumn(int column) {
        if (dirtyColumnMark_[static_cast<std::size_t>(column)] == 0) {
            dirtyColumnMark_[static_cast<std::size_t>(column)] = 1;
            dirtyColumns_.push_back(column);
        }
    }

    void markDirtyRow(int row) {
        if (dirtyRowMark_[static_cast<std::size_t>(row)] == 0) {
            dirtyRowMark_[static_cast<std::size_t>(row)] = 1;
            dirtyRows_.push_back(row);
        }
    }

    [[nodiscard]] static SquaredDistance ceilDivide(SquaredDistance numerator, SquaredDistance positiveDenominator) {
        SquaredDistance quotient = numerator / positiveDenominator;
        const SquaredDistance remainder = numerator % positiveDenominator;
        if (remainder > 0) {
            ++quotient;
        }
        return quotient;
    }

    [[nodiscard]] SquaredDistance envelopeSeparation(int previousSite, int newSite) const {
        const SquaredDistance previousCoordinate = static_cast<SquaredDistance>(previousSite);
        const SquaredDistance newCoordinate = static_cast<SquaredDistance>(newSite);
        const SquaredDistance previousValue = lineInput_[static_cast<std::size_t>(previousSite)] + previousCoordinate * previousCoordinate;
        const SquaredDistance newValue = lineInput_[static_cast<std::size_t>(newSite)] + newCoordinate * newCoordinate;
        return ceilDivide(newValue - previousValue, SquaredDistance{2} * (newCoordinate - previousCoordinate));
    }

    void transformLine(int length) {
        int envelopeSize = 0;
        for (int site = 0; site < length; ++site) {
            if (lineInput_[static_cast<std::size_t>(site)] == infinity()) {
                continue;
            }

            SquaredDistance start = std::numeric_limits<SquaredDistance>::min();
            while (envelopeSize > 0) {
                start = envelopeSeparation(envelopeSites_[static_cast<std::size_t>(envelopeSize - 1)], site);
                if (start > envelopeStarts_[static_cast<std::size_t>(envelopeSize - 1)]) {
                    break;
                }
                --envelopeSize;
            }

            envelopeSites_[static_cast<std::size_t>(envelopeSize)] = site;
            envelopeStarts_[static_cast<std::size_t>(envelopeSize)] = envelopeSize == 0 ? std::numeric_limits<SquaredDistance>::min() : start;
            ++envelopeSize;
        }

        if (envelopeSize == 0) {
            std::fill_n(lineOutput_.begin(), length, infinity());
            return;
        }

        int activeEnvelope = 0;
        for (int coordinate = 0; coordinate < length; ++coordinate) {
            while (activeEnvelope + 1 < envelopeSize &&
                   envelopeStarts_[static_cast<std::size_t>(activeEnvelope + 1)] <= static_cast<SquaredDistance>(coordinate)) {
                ++activeEnvelope;
            }
            const int site = envelopeSites_[static_cast<std::size_t>(activeEnvelope)];
            const SquaredDistance delta = static_cast<SquaredDistance>(coordinate - site);
            lineOutput_[static_cast<std::size_t>(coordinate)] = lineInput_[static_cast<std::size_t>(site)] + delta * delta;
        }
    }

    int rows_ = 0;
    int columns_ = 0;
    int numPixels_ = 0;
    std::size_t numActiveSites_ = 0;
    std::vector<std::uint8_t> seedMask_;
    std::vector<SquaredDistance> verticalCost_;
    std::vector<SquaredDistance> squaredDistance_;
    std::vector<std::int8_t> pendingDelta_;
    std::vector<SquaredDistance> lineInput_;
    std::vector<SquaredDistance> lineOutput_;
    std::vector<int> envelopeSites_;
    std::vector<SquaredDistance> envelopeStarts_;
    std::vector<std::uint8_t> dirtyColumnMark_;
    std::vector<std::uint8_t> dirtyRowMark_;
    std::vector<PixelId> touchedDeltaPixels_;
    std::vector<int> dirtyColumns_;
    std::vector<int> dirtyRows_;
    bool hasResult_ = false;
};

} // namespace mmcfilters::attributes::computers::detail::distance_transform
