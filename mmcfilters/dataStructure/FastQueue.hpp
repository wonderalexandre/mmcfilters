#pragma once

#include <vector>
#include <cstddef>
#include <utility>

namespace mmcfilters {
/**
 * @brief Lightweight FIFO queue backed by a contiguous `std::vector`.
 *
 * The queue avoids the container indirection of `std::queue` by storing
 * elements in a single vector and advancing a read head during `pop()`.
 * Clearing the queue reuses the allocated storage.
 */
template <typename T> struct FastQueue {
  private:
    /** @brief Data buffer. */
    std::vector<T> data_;
    /** @brief Head. */
    size_t head_ = 0;

  public:
    /**
     * @brief Constructs an empty reusable queue.
     */
    FastQueue() = default;

    /**
     * @brief Creates an empty queue and reserves capacity for `n` elements.
     *
     * @param n Requested element count or capacity.
     */
    FastQueue(size_t n) { data_.reserve(n); }

    /**
     * @brief Reserves storage to avoid future reallocations.
     *
     * @param n Requested element count or capacity.
     */
    void reserve(size_t n) { data_.reserve(n); }

    /**
     * @brief Clears the queue and resets the read head.
     */
    void clear() {
        data_.clear();
        head_ = 0;
    }

    /**
     * @brief Returns whether the queue is empty.
     *
     * @return Whether the queue is empty.
     */
    bool empty() const { return head_ >= data_.size(); }

    /**
     * @brief Returns the number of unread elements.
     *
     * @return The number of unread elements.
     */
    size_t size() const { return data_.size() - head_; }

    /**
     * @brief Appends an element to the tail.
     *
     * @param value Value.
     */
    void push(const T& value) { data_.push_back(value); }

    /**
     * @brief Appends an element to the tail by moving it into storage.
     *
     * @param value Value.
     */
    void push(T&& value) { data_.push_back(std::move(value)); }

    /**
     * @brief Removes and returns the next element.
     *
     * @return The removed next element.
     */
    T pop() { return std::move(data_[head_++]); }

    /**
     * @brief Returns the next element without removing it.
     *
     * @return The next element without removing it.
     */
    T& front() { return data_[head_]; }

    /**
     * @brief Returns the next element without removing it.
     *
     * @return The next element without removing it.
     */
    const T& front() const { return data_[head_]; }
};
} // namespace mmcfilters
