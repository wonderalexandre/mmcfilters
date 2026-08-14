#pragma once

/**
 * @file Common.hpp
 * @brief Shared scalar conventions used by the C++ API.
 */

namespace mmcfilters {

/**
 * @brief Node identifier type used throughout the project.
 *
 * Node ids index this repository's dense internal node-slot domain. They are
 * intentionally signed so `InvalidNode` can be represented without reserving a
 * valid dense slot.
 */
using NodeId = int; // keep signed to preserve InvalidNode semantics

/**
 * @brief Pixel identifier type used by source and active construction domains.
 *
 * Pixel ids are dense zero-based indices in the domain declared by the
 * surrounding API. In particular, a topographic convention's infinity pixel
 * belongs to that convention's active domain.
 */
using PixelId = int;

/**
 * @brief Sentinel value used to denote an invalid node identifier.
 *
 * Public APIs return or accept this value when no node exists, for example for
 * a missing parent, failed lookup, or invalid smallest node.
 */
constexpr NodeId InvalidNode = -1;

/**
 * @brief Sentinel value used to denote an invalid pixel identifier.
 *
 * Pixel and node identifiers currently share the same signed representation,
 * but their sentinels remain separately named so code records the domain to
 * which an invalid identifier belongs.
 */
constexpr PixelId InvalidPixel = -1;

/**
 * @brief Returns true when `id` is not the invalid-node sentinel.
 *
 * @param id Identifier.
 * @return True when id is not the invalid-node sentinel.
 */
inline bool isValidNode(NodeId id) noexcept { return id != InvalidNode; }

/**
 * @brief Returns true when `id` is the invalid-node sentinel.
 *
 * @param id Node identifier.
 * @return True when id is the invalid-node sentinel.
 */
inline bool isInvalidNode(NodeId id) noexcept { return id == InvalidNode; }

/**
 * @brief Returns true when `id` is not the invalid-pixel sentinel.
 *
 * @param id Pixel identifier.
 * @return True when id is not the invalid-pixel sentinel.
 */
inline bool isValidPixel(PixelId id) noexcept { return id != InvalidPixel; }

/**
 * @brief Returns true when `id` is the invalid-pixel sentinel.
 *
 * @param id Pixel identifier.
 * @return True when id is the invalid-pixel sentinel.
 */
inline bool isInvalidPixel(PixelId id) noexcept { return id == InvalidPixel; }

/**
 * @brief Compile-time switch for optional diagnostic logging.
 *
 * The default public build keeps this disabled so library code remains quiet
 * unless a local development build changes the constant.
 */
constexpr bool PRINT_LOG = false;

} // namespace mmcfilters
