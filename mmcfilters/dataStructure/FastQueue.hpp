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
template <typename T>
struct FastQueue {
private:
    std::vector<T> data_;
    size_t head_ = 0;

public:
    FastQueue() = default;

    FastQueue(size_t n){
        data_.reserve(n); 
    } 

    /// Reserves storage to avoid future reallocations.
    void reserve(size_t n) { data_.reserve(n); }

    /// Clears the queue and resets the read head.
    void clear() { data_.clear(); head_ = 0; }

    /// Returns whether the queue is empty.
    bool empty() const { return head_ >= data_.size(); }

    /// Returns the number of unread elements.
    size_t size() const { return data_.size() - head_; }

    /// Appends an element to the tail.
    void push(const T& value) { data_.push_back(value); }

    void push(T&& value) { data_.push_back(std::move(value)); }

    /// Removes and returns the next element.
    T pop() { return std::move(data_[head_++]); }

    /// Returns the next element without removing it.
    T& front() { return data_[head_]; }
    const T& front() const { return data_[head_]; }
};
}
