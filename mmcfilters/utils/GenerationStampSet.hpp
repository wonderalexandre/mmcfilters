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
    using gen_t = uint32_t;

    std::unique_ptr<gen_t[]> stamp; // stamp buffer
    size_t n{0};                    // logical size
    gen_t cur{1};                   // current generation (0 means cleared)

    GenerationStampSet() = default;
    explicit GenerationStampSet(size_t n) { resize(n); }

    void resize(size_t newN) {
        n = newN;
        stamp = std::make_unique<gen_t[]>(n);
        std::fill_n(stamp.get(), n, 0);
        cur = 1;
    }

    inline void mark(size_t idx) noexcept {
        stamp[idx] = cur;
    }

    inline bool isMarked(size_t idx) const noexcept {
        return stamp[idx] == cur;
    }

    /// Removes one mark from the current logical generation.
    inline void unmark(size_t idx) noexcept {
        if (stamp[idx] == cur) {
            stamp[idx] = 0;
        }
    }

    /// Performs an O(1) logical reset by advancing the generation counter.
    void resetAll() {
        if (++cur == 0) {
            std::fill_n(stamp.get(), n, 0);
            cur = 1;
        }
    }

    /// Performs an O(N) physical clear of the whole stamp buffer.
    void clearAll() {
        std::fill_n(stamp.get(), n, 0);
        cur = 1;
    }

    gen_t generation() const noexcept { return cur; }
};

} // namespace mmcfilters
