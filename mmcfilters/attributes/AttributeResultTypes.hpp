#pragma once

#include "AttributeNames.hpp"
#include "../trees/MorphologicalTree.hpp"
#include "../utils/Common.hpp"

#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace mmcfilters {

/**
 * @brief Owning result for one computed scalar attribute layout and buffer.
 *
 * @details
 * A computed attribute result is defined by three pieces of information:
 * - the attribute layout used to interpret the flat buffer;
 * - the owned buffer storing the per-node values;
 * - the `NodeIdSpace` in which the buffer indices are expressed.
 *
 * The result is move-only on purpose: the public API treats each computed buffer
 * as having a single owner. Internal dependency reuse is implemented by
 * `detail/AttributeDependencyCache.hpp`, so ordinary callers do not need to
 * traffic in borrowed cache views.
 *
 * The public fields remain named `first` and `second` so existing direct uses
 * and structured bindings keep working naturally.
 */
struct [[nodiscard]] ComputedAttributeData {
    /// Layout used to interpret `second`; kept public for tuple-like access.
    AttributeNames first;

    /// Flat per-node attribute buffer indexed through `first`.
    std::vector<float> second;

    /// Node-id domain used by the rows of `second`.
    NodeIdSpace nodeIdSpace = NodeIdSpace::MORPHOLOGICAL_TREE;

    /**
     * @brief Takes ownership of a scalar-attribute layout and its flat buffer.
     *
     * @param attrNames Layout that maps attributes to per-node offsets.
     * @param buffer Attribute values stored in node-major order.
     * @param outputSpace Node-id domain used by the buffer rows.
     */
    ComputedAttributeData(AttributeNames attrNames, std::vector<float> buffer, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE)
        : first(std::move(attrNames)), second(std::move(buffer)), nodeIdSpace(outputSpace) {}

    ComputedAttributeData(const ComputedAttributeData&) = delete;
    ComputedAttributeData& operator=(const ComputedAttributeData&) = delete;

    /// Transfers ownership of the layout, buffer, and node-id-space marker.
    ComputedAttributeData(ComputedAttributeData&&) noexcept = default;

    ComputedAttributeData& operator=(ComputedAttributeData&&) = delete;

    /// Returns the mutable attribute layout.
    AttributeNames& attributeNames() noexcept { return this->first; }

    /// Returns the immutable attribute layout.
    const AttributeNames& attributeNames() const noexcept { return this->first; }

    /// Returns the mutable flat attribute buffer.
    std::vector<float>& values() noexcept { return this->second; }

    /// Returns the immutable flat attribute buffer.
    const std::vector<float>& values() const noexcept { return this->second; }

};

/**
 * @brief Owning result for one delta-augmented attribute layout and buffer.
 *
 * @details
 * This is the delta-aware counterpart of `ComputedAttributeData`. It is used
 * when one logical attribute is sampled at several ancestor/descendant offsets
 * around each node.
 */
struct [[nodiscard]] ComputedAttributeDataWithDelta {
    /// Delta-aware layout used to interpret `second`.
    AttributeNamesWithDelta first;

    /// Flat per-node, per-delta attribute buffer indexed through `first`.
    std::vector<float> second;

    /// Node-id domain used by the rows of `second`.
    NodeIdSpace nodeIdSpace = NodeIdSpace::MORPHOLOGICAL_TREE;

    /**
     * @brief Takes ownership of a delta-augmented layout and its flat buffer.
     *
     * @param attrNames Layout that maps `(attribute, delta)` keys to offsets.
     * @param buffer Attribute values stored in node-major order.
     * @param outputSpace Node-id domain used by the buffer rows.
     */
    ComputedAttributeDataWithDelta(AttributeNamesWithDelta attrNames, std::vector<float> buffer, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE)
        : first(std::move(attrNames)), second(std::move(buffer)), nodeIdSpace(outputSpace) {}

    ComputedAttributeDataWithDelta(const ComputedAttributeDataWithDelta&) = delete;
    ComputedAttributeDataWithDelta& operator=(const ComputedAttributeDataWithDelta&) = delete;

    /// Transfers ownership of the layout, buffer, and node-id-space marker.
    ComputedAttributeDataWithDelta(ComputedAttributeDataWithDelta&&) noexcept = default;

    ComputedAttributeDataWithDelta& operator=(ComputedAttributeDataWithDelta&&) = delete;

    /// Returns the mutable delta-aware attribute layout.
    AttributeNamesWithDelta& attributeNames() noexcept { return this->first; }

    /// Returns the immutable delta-aware attribute layout.
    const AttributeNamesWithDelta& attributeNames() const noexcept { return this->first; }

    /// Returns the mutable flat attribute buffer.
    std::vector<float>& values() noexcept { return this->second; }

    /// Returns the immutable flat attribute buffer.
    const std::vector<float>& values() const noexcept { return this->second; }
};

} // namespace mmcfilters

namespace std {

template <>
struct tuple_size<mmcfilters::ComputedAttributeData> : integral_constant<std::size_t, 2> {};

template <>
struct tuple_element<0, mmcfilters::ComputedAttributeData> {
    using type = mmcfilters::AttributeNames;
};

template <>
struct tuple_element<1, mmcfilters::ComputedAttributeData> {
    using type = std::vector<float>;
};

template <>
struct tuple_size<mmcfilters::ComputedAttributeDataWithDelta> : integral_constant<std::size_t, 2> {};

template <>
struct tuple_element<0, mmcfilters::ComputedAttributeDataWithDelta> {
    using type = mmcfilters::AttributeNamesWithDelta;
};

template <>
struct tuple_element<1, mmcfilters::ComputedAttributeDataWithDelta> {
    using type = std::vector<float>;
};

} // namespace std

namespace mmcfilters {

/**
 * @brief Tuple-like `std::get` overloads for computed attribute results.
 */
template <std::size_t I>
decltype(auto) get(ComputedAttributeData& computed) noexcept {
    if constexpr (I == 0) {
        return (computed.first);
    } else {
        return (computed.second);
    }
}

template <std::size_t I>
decltype(auto) get(const ComputedAttributeData& computed) noexcept {
    if constexpr (I == 0) {
        return (computed.first);
    } else {
        return (computed.second);
    }
}

template <std::size_t I>
decltype(auto) get(ComputedAttributeData&& computed) noexcept {
    if constexpr (I == 0) {
        return std::move(computed.first);
    } else {
        return std::move(computed.second);
    }
}

template <std::size_t I>
decltype(auto) get(ComputedAttributeDataWithDelta& computed) noexcept {
    if constexpr (I == 0) {
        return (computed.first);
    } else {
        return (computed.second);
    }
}

template <std::size_t I>
decltype(auto) get(const ComputedAttributeDataWithDelta& computed) noexcept {
    if constexpr (I == 0) {
        return (computed.first);
    } else {
        return (computed.second);
    }
}

template <std::size_t I>
decltype(auto) get(ComputedAttributeDataWithDelta&& computed) noexcept {
    if constexpr (I == 0) {
        return std::move(computed.first);
    } else {
        return std::move(computed.second);
    }
}

} // namespace mmcfilters
