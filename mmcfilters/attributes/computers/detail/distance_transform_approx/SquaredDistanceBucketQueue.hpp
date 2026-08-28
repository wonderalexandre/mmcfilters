#pragma once

#include "../../../../utils/Common.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace mmcfilters::attributes::computers::detail::distance_transform_approx {

/**
 * @brief GFT-style fixed-range FIFO bucket queue for the approximate EDT hot path.
 *
 * This production queue mirrors the `Fast*` operations of the PQueue32 used
 * by the JMIV 2025 implementation: callers establish valid pixel ids and costs,
 * each element has compact intrusive links, and decrease-key unlinks and
 * reinserts directly. The queued bucket doubles as the membership state, so no
 * separate 64-bit per-element bucket map is needed.
 */
class SquaredDistanceBucketQueue {
  public:
    using Cost = std::int64_t;

    SquaredDistanceBucketQueue(Cost maximumCost, int numElements)
        : maximumCost_(requireMaximumCost(maximumCost)), numElements_(requireNumElements(numElements)), first_(checkedBucketCount(maximumCost_), InvalidPixel),
          last_(checkedBucketCount(maximumCost_), InvalidPixel), nodes_(static_cast<std::size_t>(numElements_)) {}

    SquaredDistanceBucketQueue(const SquaredDistanceBucketQueue&) = delete;
    SquaredDistanceBucketQueue& operator=(const SquaredDistanceBucketQueue&) = delete;
    SquaredDistanceBucketQueue(SquaredDistanceBucketQueue&&) = delete;
    SquaredDistanceBucketQueue& operator=(SquaredDistanceBucketQueue&&) = delete;

    [[nodiscard]] bool empty() const noexcept { return size_ == 0; }
    [[nodiscard]] std::size_t size() const noexcept { return size_; }
    [[nodiscard]] Cost maximumCost() const noexcept { return maximumCost_; }
    [[nodiscard]] bool contains(PixelId element) const noexcept { return nodes_[static_cast<std::size_t>(element)].bucket != InvalidBucket; }

    void insert(PixelId element, Cost cost) noexcept { link(element, static_cast<Bucket>(cost)); }

    void update(PixelId element, Cost cost) noexcept {
        Node& node = nodes_[static_cast<std::size_t>(element)];
        const Bucket nextBucket = static_cast<Bucket>(cost);
        if (node.bucket == InvalidBucket) {
            link(element, nextBucket);
            return;
        }
        if (node.bucket == nextBucket) {
            return;
        }
        unlink(element, node.bucket);
        link(element, nextBucket);
    }

    void erase(PixelId element) noexcept {
        Node& node = nodes_[static_cast<std::size_t>(element)];
        if (node.bucket == InvalidBucket) {
            return;
        }
        unlink(element, node.bucket);
    }

    [[nodiscard]] PixelId popMinimumFifo() noexcept {
        advanceMinimumBucket();
        const Bucket bucket = minimumBucket_;
        const PixelId element = first_[static_cast<std::size_t>(bucket)];
        unlink(element, bucket);
        return element;
    }

  private:
    using Bucket = std::int32_t;
    inline static constexpr Bucket InvalidBucket = Bucket{-1};

    struct Node {
        PixelId next = InvalidPixel;
        PixelId previous = InvalidPixel;
        Bucket bucket = InvalidBucket;
    };

    [[nodiscard]] static Cost requireMaximumCost(Cost value) {
        if (value < 0 || value > static_cast<Cost>(std::numeric_limits<Bucket>::max())) {
            throw std::invalid_argument("Squared-distance queue requires a non-negative 32-bit bucket domain.");
        }
        return value;
    }

    [[nodiscard]] static int requireNumElements(int value) {
        if (value <= 0) {
            throw std::invalid_argument("Squared-distance queue requires a non-empty element domain.");
        }
        return value;
    }

    [[nodiscard]] static std::size_t checkedBucketCount(Cost maximumCost) { return static_cast<std::size_t>(maximumCost) + 1; }

    void link(PixelId element, Bucket bucket) noexcept {
        Node& node = nodes_[static_cast<std::size_t>(element)];
        const std::size_t bucketIndex = static_cast<std::size_t>(bucket);
        const PixelId tail = last_[bucketIndex];
        node.previous = tail;
        node.next = InvalidPixel;
        node.bucket = bucket;
        if (tail == InvalidPixel) {
            first_[bucketIndex] = element;
        } else {
            nodes_[static_cast<std::size_t>(tail)].next = element;
        }
        last_[bucketIndex] = element;
        ++size_;
        if (minimumBucket_ == InvalidBucket || bucket < minimumBucket_) {
            minimumBucket_ = bucket;
        }
    }

    void unlink(PixelId element, Bucket bucket) noexcept {
        Node& node = nodes_[static_cast<std::size_t>(element)];
        const std::size_t bucketIndex = static_cast<std::size_t>(bucket);
        if (node.previous == InvalidPixel) {
            first_[bucketIndex] = node.next;
        } else {
            nodes_[static_cast<std::size_t>(node.previous)].next = node.next;
        }
        if (node.next == InvalidPixel) {
            last_[bucketIndex] = node.previous;
        } else {
            nodes_[static_cast<std::size_t>(node.next)].previous = node.previous;
        }
        node.bucket = InvalidBucket;
        --size_;
        if (size_ == 0) {
            minimumBucket_ = InvalidBucket;
        }
    }

    void advanceMinimumBucket() noexcept {
        while (first_[static_cast<std::size_t>(minimumBucket_)] == InvalidPixel) {
            ++minimumBucket_;
        }
    }

    Cost maximumCost_ = 0;
    int numElements_ = 0;
    std::vector<PixelId> first_;
    std::vector<PixelId> last_;
    std::vector<Node> nodes_;
    Bucket minimumBucket_ = InvalidBucket;
    std::size_t size_ = 0;
};

} // namespace mmcfilters::attributes::computers::detail::distance_transform_approx
