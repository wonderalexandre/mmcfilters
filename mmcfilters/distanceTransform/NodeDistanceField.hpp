#pragma once

#include "DistanceTransformTypes.hpp"

#include <cstdint>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcfilters {

/** @brief Self-contained compact per-pixel distance field for one live node. */
class NodeDistanceField {
  public:
    /**
     * @brief Constructs a self-contained compact node field.
     * @param node Live internal node represented by the samples.
     * @param rows Number of rows in the original image domain.
     * @param columns Number of columns in the original image domain.
     * @param accuracy Numerical backend that produced the samples.
     * @param supportPixels Row-major support pixel ids.
     * @param squaredDistances Squared distances aligned with `supportPixels`.
     */
    NodeDistanceField(NodeId node, int rows, int columns, DistanceTransformAccuracy accuracy, std::vector<PixelId> supportPixels,
                      std::vector<std::int64_t> squaredDistances)
        : node_(node), rows_(rows), columns_(columns), accuracy_(accuracy), supportPixels_(std::move(supportPixels)),
          squaredDistances_(std::move(squaredDistances)) {
        if (node < 0 || rows <= 0 || columns <= 0 || (accuracy != DistanceTransformAccuracy::Exact && accuracy != DistanceTransformAccuracy::Approximate) ||
            supportPixels_.empty() || supportPixels_.size() != squaredDistances_.size()) {
            throw std::invalid_argument("NodeDistanceField requires a live node and matching non-empty sample arrays.");
        }
        const std::size_t domainSize = static_cast<std::size_t>(rows_) * static_cast<std::size_t>(columns_);
        for (std::size_t index = 0; index < supportPixels_.size(); ++index) {
            if (supportPixels_[index] < 0 || static_cast<std::size_t>(supportPixels_[index]) >= domainSize || squaredDistances_[index] < 0) {
                throw std::invalid_argument("NodeDistanceField samples must belong to the domain and have non-negative squared distances.");
            }
        }
    }

    /** @brief Returns the represented live internal node id. @return Node id. */
    [[nodiscard]] NodeId node() const noexcept { return node_; }
    /** @brief Returns the original domain height. @return Number of rows. */
    [[nodiscard]] int rows() const noexcept { return rows_; }
    /** @brief Returns the original domain width. @return Number of columns. */
    [[nodiscard]] int columns() const noexcept { return columns_; }
    /** @brief Returns the numerical backend used to produce the field. @return Accuracy selection. */
    [[nodiscard]] DistanceTransformAccuracy accuracy() const noexcept { return accuracy_; }
    /** @brief Returns the number of compact support samples. @return Sample count. */
    [[nodiscard]] std::size_t size() const noexcept { return supportPixels_.size(); }
    /** @brief Returns the row-major support pixel ids. @return Borrowed support span. */
    [[nodiscard]] std::span<const PixelId> supportPixels() const noexcept { return supportPixels_; }
    /** @brief Returns squared distances aligned with `supportPixels()`. @return Borrowed squared-distance span. */
    [[nodiscard]] std::span<const std::int64_t> squaredDistances() const noexcept { return squaredDistances_; }

    /**
     * @brief Returns a full-domain row-major field with `outsideValue` off support.
     * @param outsideValue Value assigned to pixels outside the represented support.
     * @return Dense row-major squared-distance field.
     */
    [[nodiscard]] std::vector<std::int64_t> denseSquaredDistances(std::int64_t outsideValue = -1) const {
        const std::size_t domainSize = static_cast<std::size_t>(rows_) * static_cast<std::size_t>(columns_);
        std::vector<std::int64_t> dense(domainSize, outsideValue);
        for (std::size_t index = 0; index < supportPixels_.size(); ++index) {
            dense[static_cast<std::size_t>(supportPixels_[index])] = squaredDistances_[index];
        }
        return dense;
    }

    /** @cond */
  private:
    NodeId node_ = InvalidNode;
    int rows_ = 0;
    int columns_ = 0;
    DistanceTransformAccuracy accuracy_ = DistanceTransformAccuracy::Exact;
    std::vector<PixelId> supportPixels_;
    std::vector<std::int64_t> squaredDistances_;
    /** @endcond */
};

} // namespace mmcfilters
