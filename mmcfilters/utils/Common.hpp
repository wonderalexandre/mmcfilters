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
 * @brief Sentinel value used to denote an invalid node identifier.
 *
 * Public APIs return or accept this value when no node exists, for example for
 * a missing parent, failed lookup, or invalid proper-part owner.
 */
constexpr NodeId InvalidNode = -1;

/**
 * @brief Returns true when `id` is not the invalid-node sentinel.
 *
 * @param id Identifier used by the operation.
 * @return True when id is not the invalid-node sentinel.
 */
inline bool isValidNode(NodeId id) noexcept { return id != InvalidNode; }

/**
 * @brief Returns true when `id` is the invalid-node sentinel.
 *
 * @param id Identifier used by the operation.
 * @return True when id is the invalid-node sentinel.
 */
inline bool isInvalid(NodeId id) noexcept { return id == InvalidNode; }

/**
 * @brief Compile-time switch for optional diagnostic logging.
 *
 * The default public build keeps this disabled so library code remains quiet
 * unless a local development build changes the constant.
 */
constexpr bool PRINT_LOG = false;

} // namespace mmcfilters
