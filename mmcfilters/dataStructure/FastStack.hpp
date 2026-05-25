#pragma once

#include <vector>
#include <cstddef>
#include <utility>

namespace mmcfilters {
/**
 * @brief Lightweight LIFO stack backed by a contiguous `std::vector`.
 *
 * `FastStack<T>` provides the core stack operations needed by DFS-like
 * traversals while keeping memory layout simple and allocation behaviour
 * predictable.
 */
template <typename T>
struct FastStack {
private:
    std::vector<T> data_;

public:
    FastStack() = default;

    /**
     * @brief Creates an empty stack and reserves capacity for `n` elements.
     */
    explicit FastStack(size_t n) {
        data_.reserve(n);
    }

    /// Reserves storage to avoid future reallocations.
    void reserve(size_t n) { data_.reserve(n); }

    /// Clears the stack.
    void clear() { data_.clear(); }

    /// Returns whether the stack is empty.
    bool empty() const { return data_.empty(); }

    /// Returns the number of stored elements.
    size_t size() const { return data_.size(); }

    /// Pushes an element on top of the stack.
    void push(const T& value) { data_.push_back(value); }

    /// Pushes an element on top of the stack by moving it into storage.
    void push(T&& value) { data_.push_back(std::move(value)); }

    /// Removes and returns the top element.
    T pop() {
        T value = std::move(data_.back());
        data_.pop_back();
        return value;
    }

    /// Returns the top element without removing it.
    T& top() { return data_.back(); }

    /// Returns the top element without removing it.
    const T& top() const { return data_.back(); }
};
}
