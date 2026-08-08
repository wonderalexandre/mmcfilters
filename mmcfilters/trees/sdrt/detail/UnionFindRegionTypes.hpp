#pragma once

/**
 * @file UnionFindRegionTypes.hpp
 * @brief Dense region identifiers shared by residual-tree assembly.
 */

namespace mmcfilters::sdrt::detail {

using RegionId = int;                         ///< Dense initial flat-zone id and union-find slot index.
inline constexpr RegionId InvalidRegion = -1; ///< Sentinel outside the valid id domain.

} // namespace mmcfilters::sdrt::detail
