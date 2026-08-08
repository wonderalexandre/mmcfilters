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
template <typename T> struct FastStack {
  private:
    /** @brief Stores the data. */
    std::vector<T> data_;

  public:
    /**
     * @brief Constructs an empty reusable stack.
     */
    FastStack() = default;

    /**
     * @brief Creates an empty stack and reserves capacity for `n` elements.
     *
     * @param n Requested element count or capacity.
     */
    explicit FastStack(size_t n) { data_.reserve(n); }

    /**
     * @brief Reserves storage to avoid future reallocations.
     *
     * @param n Requested element count or capacity.
     */
    void reserve(size_t n) { data_.reserve(n); }

    /**
     * @brief Clears the stack.
     */
    void clear() { data_.clear(); }

    /**
     * @brief Returns whether the stack is empty.
     *
     * @return Whether the stack is empty.
     */
    bool empty() const { return data_.empty(); }

    /**
     * @brief Returns the number of stored elements.
     *
     * @return The number of stored elements.
     */
    size_t size() const { return data_.size(); }

    /**
     * @brief Pushes an element on top of the stack.
     *
     * @param value Value used by the operation.
     */
    void push(const T& value) { data_.push_back(value); }

    /**
     * @brief Pushes an element on top of the stack by moving it into storage.
     *
     * @param value Value used by the operation.
     */
    void push(T&& value) { data_.push_back(std::move(value)); }

    /**
     * @brief Removes and returns the top element.
     *
     * @return The removed top element.
     */
    T pop() {
        T value = std::move(data_.back());
        data_.pop_back();
        return value;
    }

    /**
     * @brief Returns the top element without removing it.
     *
     * @return The top element without removing it.
     */
    T& top() { return data_.back(); }

    /**
     * @brief Returns the top element without removing it.
     *
     * @return The top element without removing it.
     */
    const T& top() const { return data_.back(); }
};
} // namespace mmcfilters
