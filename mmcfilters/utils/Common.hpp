#pragma once


// ---------------------------------------------------------------------------
// Assertion control
// ---------------------------------------------------------------------------
#if defined(MMCFILTERS_ENABLE_ASSERTS)
#  ifdef NDEBUG
#    undef NDEBUG
#  endif
#endif
#include <cassert>      // assert()

// ---------------------------------------------------------------------------
// Standard containers and algorithms
// ---------------------------------------------------------------------------
#include <list>          // Doubly-linked list
#include <vector>        // Resizable contiguous container
#include <array>         // Fixed-size array
#include <deque>         // Double-ended queue
#include <stack>         // Stack adaptor
#include <unordered_set> // Hash set
#include <unordered_map> // Hash map
#include <typeindex>     // std::type_index for type-aware maps and sets
#include <set>           // Ordered set / multiset
#include <map>           // Ordered map / multimap
#include <span>          // Non-owning contiguous view
#include <tuple>         // Fixed heterogeneous tuples
#include <algorithm>     // Generic algorithms
#include <iterator>      // Iterator utilities
#include <utility>       // std::pair, std::move, std::swap

// ---------------------------------------------------------------------------
// General utilities
// ---------------------------------------------------------------------------
#include <cstdint>   // Fixed-width integer types
#include <limits>    // Numeric limits

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <cmath>     // Standard maths functions
#include <iostream>  // Standard streams
#include <string>    // std::string
#include <iomanip>   // Stream manipulators
#include <numeric>   // Numeric algorithms
#include <stdexcept> // Standard exceptions
#include <sstream>   // String-based streams
#include <numbers>   // C++20 mathematical constants


// ---------------------------------------------------------------------------
// Memory management, callables, and metaprogramming
// ---------------------------------------------------------------------------
#include <memory>       // Smart pointers
#include <variant>      // Type-safe tagged union
#include <optional>     // Optional values
#include <functional>   // std::function, std::bind, std::hash
#include <type_traits>  // Traits and compile-time utilities

// ---------------------------------------------------------------------------
// mmcfilters common project-wide includes
// ---------------------------------------------------------------------------
#include "Image.hpp"
#include "../dataStructure/FastStack.hpp"
#include "../dataStructure/FastQueue.hpp"


namespace mmcfilters {

/// Compile-time switches for optional logging and debug traces.
constexpr bool PRINT_LOG   = false;
constexpr bool PRINT_DEBUG = false;

/// Node identifier type used throughout the project.
using NodeId = int; // keep signed to preserve InvalidNode semantics
/// Sentinel value used to denote an invalid node identifier.
constexpr NodeId InvalidNode = -1;
inline bool isValidNode(NodeId id) noexcept { return id != InvalidNode;}
inline bool isInvalid(NodeId id) noexcept { return id == InvalidNode; }

/// Canonical value type stored for node altitudes / gray levels in component trees.
using AltitudeType = int;

/// Dense altitude buffer indexed by internal node id.
using AltitudeBuffer = std::vector<AltitudeType>;

/// Signed/arithmetic type used for altitude differences such as residues.
using AltitudeDiffType = decltype(std::declval<AltitudeType>() - std::declval<AltitudeType>());

// Optional deprecation hook for the legacy weighted API that still lives on
// MorphologicalTree during the transition to the split between topology/
// ownership and external altitude buffers.
#if defined(MMCFILTERS_ENABLE_LEGACY_WEIGHTED_TREE_API_DEPRECATION)
#define MMCFILTERS_LEGACY_WEIGHTED_TREE_API [[deprecated("Use WeightedMorphologicalTree or tree_altitude_ops with an explicit altitude buffer.")]]
#else
#define MMCFILTERS_LEGACY_WEIGHTED_TREE_API
#endif


} // namespace mmcfilters
