#pragma once

/**
 * @file Contract.hpp
 * @brief Compile-time policy for defensive checks at untrusted API boundaries.
 *
 * Public signatures are independent of the selected mode.  Checked builds
 * reject invalid caller input; unchecked builds discard those boundary checks
 * for controlled scientific benchmarks.  Internal kernels must not use these
 * macros: they are entered only after their caller established the required
 * domains.
 */

/** @brief Compile-time value selecting a build without defensive boundary validation. */
#define MMCFILTERS_CONTRACT_UNCHECKED 0

/** @brief Compile-time value selecting a build with defensive boundary validation. */
#define MMCFILTERS_CONTRACT_CHECKED 1

#ifndef MMCFILTERS_CONTRACT_MODE
/** @brief Contract mode selected for the current translation unit. */
#define MMCFILTERS_CONTRACT_MODE MMCFILTERS_CONTRACT_CHECKED
#endif

#if MMCFILTERS_CONTRACT_MODE != MMCFILTERS_CONTRACT_UNCHECKED && MMCFILTERS_CONTRACT_MODE != MMCFILTERS_CONTRACT_CHECKED
#error "MMCFILTERS_CONTRACT_MODE must be MMCFILTERS_CONTRACT_CHECKED or MMCFILTERS_CONTRACT_UNCHECKED."
#endif

namespace mmcfilters::contract {

/** @brief Available compile-time policies for defensive API-boundary validation. */
enum class Mode {
    /** Caller-established preconditions are trusted without defensive checks. */
    Unchecked = MMCFILTERS_CONTRACT_UNCHECKED,
    /** Public boundaries validate caller-provided preconditions. */
    Checked = MMCFILTERS_CONTRACT_CHECKED
};

/** @brief Contract policy selected when this translation unit was compiled. */
inline constexpr Mode mode = static_cast<Mode>(MMCFILTERS_CONTRACT_MODE);

/** @brief Whether defensive API-boundary validation is enabled in this build. */
inline constexpr bool validationsEnabled = mode == Mode::Checked;

} // namespace mmcfilters::contract

/**
 * @brief Evaluates a caller precondition and its failure action only in checked builds.
 *
 * @param condition Boolean caller precondition established at the public API boundary.
 * @param ... Statements that report or handle a failed precondition.
 */
#define MMCFILTERS_CONTRACT_REQUIRE(condition, ...)                                                                                                            \
    do {                                                                                                                                                       \
        if constexpr (::mmcfilters::contract::validationsEnabled) {                                                                                            \
            if (!(condition)) {                                                                                                                                \
                __VA_ARGS__;                                                                                                                                   \
            }                                                                                                                                                  \
        }                                                                                                                                                      \
    } while (false)

/**
 * @brief Executes validation statements only when defensive checks are enabled.
 *
 * @param ... Validation statements omitted from unchecked builds.
 */
#define MMCFILTERS_CONTRACT_CHECKED_ONLY(...)                                                                                                                  \
    do {                                                                                                                                                       \
        if constexpr (::mmcfilters::contract::validationsEnabled) {                                                                                            \
            __VA_ARGS__;                                                                                                                                       \
        }                                                                                                                                                      \
    } while (false)
