#pragma once

/**
 * @file Assert.hpp
 * @brief Project assertion configuration wrapper.
 *
 * Defining `MMCFILTERS_ENABLE_ASSERTS` forces standard `assert` checks on even
 * when the including translation unit was compiled with `NDEBUG`. This is a
 * development and validation aid; public runtime precondition checks should use
 * explicit exceptions instead of relying on assertions.
 */

#if defined(MMCFILTERS_ENABLE_ASSERTS)
#ifdef NDEBUG
#undef NDEBUG
#endif
#endif

#include <cassert>
