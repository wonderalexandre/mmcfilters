#pragma once

/**
 * @file ResidualTreeBuildStatistics.hpp
 * @brief Minimal correctness-oriented diagnostics for residual-tree construction.
 */

#include <cstddef>

namespace mmcfilters::sdrt {

/**
 * @brief Diagnostics retained by the synchronized min-tree/max-tree builders.
 *
 * These counters participate in consistency checks and regression tests.
 * Detailed profiling remains the responsibility of external benchmarks.
 */
struct ResidualTreeBuildStatistics {
    std::size_t residualEvents = 0;                  ///< Emitted non-root residual nodes.
    std::size_t rejectedExtrema = 0;                 ///< Candidates rejected by saturation.
    std::size_t complementTraversalCertificates = 0; ///< Exact fallback certificates.
};

} // namespace mmcfilters::sdrt
