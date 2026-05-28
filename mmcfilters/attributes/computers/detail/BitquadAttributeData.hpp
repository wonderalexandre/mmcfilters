#pragma once

namespace mmcfilters::attributes::computers::detail {

/**
 * @brief Compact per-node counters grouped by bitquad family.
 *
 * @details
 * `empty` corresponds to state `0000`, `qd` to diagonal two-pixel states, and
 * the other fields are the usual Q1/Q2/Q3/Q4 families. These counters are the
 * internal bucket consumed by the bitquad local-event implementation and its
 * scalar materialization helper.
 */
struct BitquadFamilyCounts {
    /// Count of all-background 2x2 configurations.
    int empty = 0;

    /// Count of Q1 configurations with one foreground sample.
    int q1 = 0;

    /// Count of Q2 configurations with adjacent foreground samples.
    int q2 = 0;

    /// Count of diagonal Q2 configurations.
    int qd = 0;

    /// Count of Q3 configurations with three foreground samples.
    int q3 = 0;

    /// Count of Q4 configurations with four foreground samples.
    int q4 = 0;
};

} // namespace mmcfilters::attributes::computers::detail
