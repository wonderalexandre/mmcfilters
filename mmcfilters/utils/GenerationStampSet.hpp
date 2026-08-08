#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>

namespace mmcfilters {

/**
 * @brief Efficient visited-set implementation based on generation stamps.
 *
 * `GenerationStampSet` stores one integer stamp per index. Instead of clearing
 * the whole buffer between traversals, the current generation counter is
 * incremented and used as the logical mark. This keeps the common reset path
 * at O(1) while still allowing a full clear when the counter wraps around.
 *
 * @code
 * GenerationStampSet visited(numNodes);
 * visited.mark(nodeIdx);
 *
 * if (!visited.isMarked(otherIdx)) {
 *     // process one still-unvisited node
 * }
 *
 * visited.resetAll();  // O(1) logical reset for the next traversal
 * @endcode
 */
struct GenerationStampSet {
    /// Integer type used for individual generation stamps.
    using gen_t = uint32_t;

    /// Stamp buffer with one entry per logical index.
    std::unique_ptr<gen_t[]> stamp;

    /// Number of logical entries in `stamp`.
    size_t n{0};

    /// Current generation; zero is reserved for physically cleared slots.
    gen_t cur{1};

    /**
     * @brief Constructs a default `GenerationStampSet`.
     */
    GenerationStampSet() = default;

    /**
     * @brief Creates a stamp set with `n` logical entries.
     *
     * @param n Requested element count or capacity.
     */
    explicit GenerationStampSet(size_t n) { resize(n); }

    /**
     * @brief Resizes the stamp buffer and physically clears all marks.
     *
     * @param newN New number of logical entries.
     */
    void resize(size_t newN) {
        n = newN;
        stamp = std::make_unique<gen_t[]>(n);
        std::fill_n(stamp.get(), n, 0);
        cur = 1;
    }

    /**
     * @brief Marks `idx` in the current generation.
     *
     * @param idx Zero-based index used by the operation.
     */
    inline void mark(size_t idx) noexcept { stamp[idx] = cur; }

    /**
     * @brief Returns true when `idx` is marked in the current generation.
     *
     * @param idx Zero-based index used by the operation.
     * @return True when idx is marked in the current generation.
     */
    inline bool isMarked(size_t idx) const noexcept { return stamp[idx] == cur; }

    /**
     * @brief Removes one mark from the current logical generation.
     *
     * @param idx Zero-based index used by the operation.
     */
    inline void unmark(size_t idx) noexcept {
        if (stamp[idx] == cur) {
            stamp[idx] = 0;
        }
    }

    /**
     * @brief Performs an O(1) logical reset by advancing the generation counter.
     */
    void resetAll() {
        if (++cur == 0) {
            std::fill_n(stamp.get(), n, 0);
            cur = 1;
        }
    }

    /**
     * @brief Performs an O(N) physical clear of the whole stamp buffer.
     */
    void clearAll() {
        std::fill_n(stamp.get(), n, 0);
        cur = 1;
    }

    /**
     * @brief Returns the current generation counter.
     *
     * @return The current generation counter.
     */
    gen_t generation() const noexcept { return cur; }
};

} // namespace mmcfilters
