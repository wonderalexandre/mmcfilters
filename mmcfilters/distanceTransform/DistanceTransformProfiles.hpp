#pragma once

#include "DistanceTransformTypes.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcfilters {

/** @brief Compact CSR-like sparse distance histograms for every live tree node. */
class DistanceTransformProfiles {
  public:
    /** @brief Offset and length of one node histogram in the flat bin array. */
    struct Span {
        /// First bin belonging to the node in the flat bin buffer.
        std::size_t offset = 0;
        /// Number of consecutive bins belonging to the node.
        std::size_t size = 0;
    };

    /**
     * @brief Constructs validated CSR-like profiles.
     * @param accuracy Numerical backend that produced the profiles.
     * @param bins Flat sorted-bin storage.
     * @param spans Per-node offsets and lengths into `bins`.
     * @param alive Per-node liveness flags aligned with `spans`.
     */
    DistanceTransformProfiles(DistanceTransformAccuracy accuracy, std::vector<DistanceTransformHistogramBin> bins, std::vector<Span> spans,
                              std::vector<std::uint8_t> alive)
        : accuracy_(accuracy), bins_(std::move(bins)), spans_(std::move(spans)), alive_(std::move(alive)) {
        if ((accuracy_ != DistanceTransformAccuracy::Exact && accuracy_ != DistanceTransformAccuracy::Approximate) || spans_.size() != alive_.size()) {
            throw std::invalid_argument("DistanceTransformProfiles requires aligned node spans and liveness flags.");
        }
        for (std::size_t node = 0; node < spans_.size(); ++node) {
            const Span span = spans_[node];
            if (span.offset > bins_.size() || span.size > bins_.size() - span.offset || (alive_[node] != 0 && span.size == 0) ||
                (alive_[node] == 0 && span.size != 0)) {
                throw std::invalid_argument("DistanceTransformProfiles received an invalid node histogram span.");
            }
            std::int64_t previousDistance = -1;
            for (std::size_t binIndex = span.offset; binIndex < span.offset + span.size; ++binIndex) {
                const DistanceTransformHistogramBin& bin = bins_[binIndex];
                if (bin.squaredDistance < 0 || bin.squaredDistance <= previousDistance || bin.count == 0) {
                    throw std::invalid_argument("DistanceTransformProfiles requires positive counts at strictly increasing non-negative levels.");
                }
                previousDistance = bin.squaredDistance;
            }
        }
    }

    /** @brief Returns the numerical backend used to produce the profiles. @return Accuracy selection. */
    [[nodiscard]] DistanceTransformAccuracy accuracy() const noexcept { return accuracy_; }
    /** @brief Returns the number of addressable internal node slots. @return Node-slot count. */
    [[nodiscard]] std::size_t numNodeSlots() const noexcept { return spans_.size(); }
    /** @brief Returns the number of bins stored across all live nodes. @return Flat bin count. */
    [[nodiscard]] std::size_t numBins() const noexcept { return bins_.size(); }

    /**
     * @brief Returns the sorted sparse histogram for one live internal node.
     * @param node Live internal node id.
     * @return Borrowed span of strictly increasing squared-distance bins.
     */
    [[nodiscard]] std::span<const DistanceTransformHistogramBin> histogram(NodeId node) const {
        const Span span = checkedSpan(node);
        return std::span<const DistanceTransformHistogramBin>(bins_.data() + span.offset, span.size);
    }

    /**
     * @brief Computes counts at or above every represented squared-distance level.
     * @param node Live internal node id.
     * @return Owning survival-profile bins in increasing threshold order.
     */
    [[nodiscard]] std::vector<DistanceTransformSurvivalBin> survivalProfile(NodeId node) const {
        const std::span<const DistanceTransformHistogramBin> nodeHistogram = histogram(node);
        std::uint64_t total = 0;
        for (const DistanceTransformHistogramBin& bin : nodeHistogram) {
            if (bin.count > std::numeric_limits<std::uint64_t>::max() - total) {
                throw std::overflow_error("Distance-transform survival profile exceeds the 64-bit count domain.");
            }
            total += bin.count;
        }
        std::vector<DistanceTransformSurvivalBin> result;
        result.reserve(nodeHistogram.size());
        std::uint64_t below = 0;
        for (const DistanceTransformHistogramBin& bin : nodeHistogram) {
            result.push_back({bin.squaredDistance, total - below});
            below += bin.count;
        }
        return result;
    }

    /** @cond */
  private:
    [[nodiscard]] Span checkedSpan(NodeId node) const {
        if (node < 0 || static_cast<std::size_t>(node) >= spans_.size() || alive_[static_cast<std::size_t>(node)] == 0) {
            throw std::out_of_range("Distance-transform profile access requires a live internal node id.");
        }
        return spans_[static_cast<std::size_t>(node)];
    }

    DistanceTransformAccuracy accuracy_ = DistanceTransformAccuracy::Exact;
    std::vector<DistanceTransformHistogramBin> bins_;
    std::vector<Span> spans_;
    std::vector<std::uint8_t> alive_;
    /** @endcond */
};

} // namespace mmcfilters
