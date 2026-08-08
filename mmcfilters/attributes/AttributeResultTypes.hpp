#pragma once

#include "AttributeNames.hpp"
#include "../trees/MorphologicalTree.hpp"
#include "../utils/Common.hpp"

#include <concepts>
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
template <std::floating_point Real = float> struct [[nodiscard]] ComputedAttributeData {
    /// Layout used to interpret `second`; kept public for tuple-like access.
    AttributeNames first;

    /// Flat per-node attribute buffer indexed through `first`.
    std::vector<Real> second;

    /// Node-id domain used by the rows of `second`.
    NodeIdSpace nodeIdSpace = NodeIdSpace::MORPHOLOGICAL_TREE;

    /**
     * @brief Takes ownership of a scalar-attribute layout and its flat buffer.
     *
     * @param attrNames Layout that maps attributes to per-node offsets.
     * @param buffer Attribute values stored in node-major order.
     * @param outputSpace Node-id domain used by the buffer rows.
     */
    ComputedAttributeData(AttributeNames attrNames, std::vector<Real> buffer, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE)
        : first(std::move(attrNames)), second(std::move(buffer)), nodeIdSpace(outputSpace) {}

    /**
     * @brief Disables copy construction.
     */
    ComputedAttributeData(const ComputedAttributeData&) = delete;
    /**
     * @brief Disables copy assignment.
     */
    ComputedAttributeData& operator=(const ComputedAttributeData&) = delete;

    /**
     * @brief Transfers ownership of the layout, buffer, and node-id-space marker.
     */
    ComputedAttributeData(ComputedAttributeData&&) noexcept = default;

    /**
     * @brief Disables move assignment.
     */
    ComputedAttributeData& operator=(ComputedAttributeData&&) = delete;

    /**
     * @brief Returns the mutable attribute layout.
     *
     * @return The mutable attribute layout.
     */
    AttributeNames& attributeNames() noexcept { return this->first; }

    /**
     * @brief Returns the immutable attribute layout.
     *
     * @return The immutable attribute layout.
     */
    const AttributeNames& attributeNames() const noexcept { return this->first; }

    /**
     * @brief Returns the mutable flat attribute buffer.
     *
     * @return The mutable flat attribute buffer.
     */
    std::vector<Real>& values() noexcept { return this->second; }

    /**
     * @brief Returns the immutable flat attribute buffer.
     *
     * @return The immutable flat attribute buffer.
     */
    const std::vector<Real>& values() const noexcept { return this->second; }
};

/**
 * @brief Owning result for one delta-augmented attribute layout and buffer.
 *
 * @details
 * This is the delta-aware counterpart of `ComputedAttributeData<Real>`. It is
 * used when one logical attribute is sampled at several ancestor/descendant
 * offsets around each node.
 */
template <std::floating_point Real = float> struct [[nodiscard]] ComputedAttributeDataWithDelta {
    /// Delta-aware layout used to interpret `second`.
    AttributeNamesWithDelta first;

    /// Flat per-node, per-delta attribute buffer indexed through `first`.
    std::vector<Real> second;

    /// Node-id domain used by the rows of `second`.
    NodeIdSpace nodeIdSpace = NodeIdSpace::MORPHOLOGICAL_TREE;

    /**
     * @brief Takes ownership of a delta-augmented layout and its flat buffer.
     *
     * @param attrNames Layout that maps `(attribute, delta)` keys to offsets.
     * @param buffer Attribute values stored in node-major order.
     * @param outputSpace Node-id domain used by the buffer rows.
     */
    ComputedAttributeDataWithDelta(AttributeNamesWithDelta attrNames, std::vector<Real> buffer, NodeIdSpace outputSpace = NodeIdSpace::MORPHOLOGICAL_TREE)
        : first(std::move(attrNames)), second(std::move(buffer)), nodeIdSpace(outputSpace) {}

    /**
     * @brief Disables copy construction.
     */
    ComputedAttributeDataWithDelta(const ComputedAttributeDataWithDelta&) = delete;

    /**
     * @brief Disables copy assignment.
     */
    ComputedAttributeDataWithDelta& operator=(const ComputedAttributeDataWithDelta&) = delete;

    /**
     * @brief Transfers ownership of the layout, buffer, and node-id-space marker.
     */
    ComputedAttributeDataWithDelta(ComputedAttributeDataWithDelta&&) noexcept = default;

    /**
     * @brief Disables move assignment.
     */
    ComputedAttributeDataWithDelta& operator=(ComputedAttributeDataWithDelta&&) = delete;

    /**
     * @brief Returns the mutable delta-aware attribute layout.
     *
     * @return The mutable delta-aware attribute layout.
     */
    AttributeNamesWithDelta& attributeNames() noexcept { return this->first; }

    /**
     * @brief Returns the immutable delta-aware attribute layout.
     *
     * @return The immutable delta-aware attribute layout.
     */
    const AttributeNamesWithDelta& attributeNames() const noexcept { return this->first; }

    /**
     * @brief Returns the mutable flat attribute buffer.
     *
     * @return The mutable flat attribute buffer.
     */
    std::vector<Real>& values() noexcept { return this->second; }

    /**
     * @brief Returns the immutable flat attribute buffer.
     *
     * @return The immutable flat attribute buffer.
     */
    const std::vector<Real>& values() const noexcept { return this->second; }
};

} // namespace mmcfilters

namespace std {

template <std::floating_point Real> struct tuple_size<mmcfilters::ComputedAttributeData<Real>> : integral_constant<std::size_t, 2> {};

/**
 * @brief Provides tuple-element metadata for structured binding support.
 *
 * @tparam Real Floating-point type stored in the attribute buffer.
 */
template <std::floating_point Real> struct tuple_element<0, mmcfilters::ComputedAttributeData<Real>> {
    /** @brief Defines the `type` alias used by the component. */
    using type = mmcfilters::AttributeNames;
};

/**
 * @brief Provides tuple-element metadata for structured binding support.
 *
 * @tparam Real Floating-point type stored in the attribute buffer.
 */
template <std::floating_point Real> struct tuple_element<1, mmcfilters::ComputedAttributeData<Real>> {
    /** @brief Defines the `type` alias used by the component. */
    using type = std::vector<Real>;
};

template <std::floating_point Real> struct tuple_size<mmcfilters::ComputedAttributeDataWithDelta<Real>> : integral_constant<std::size_t, 2> {};

/**
 * @brief Provides tuple-element metadata for structured binding support.
 *
 * @tparam Real Floating-point type stored in the attribute and delta buffers.
 */
template <std::floating_point Real> struct tuple_element<0, mmcfilters::ComputedAttributeDataWithDelta<Real>> {
    /** @brief Defines the `type` alias used by the component. */
    using type = mmcfilters::AttributeNamesWithDelta;
};

/**
 * @brief Provides tuple-element metadata for structured binding support.
 *
 * @tparam Real Floating-point type stored in the attribute and delta buffers.
 */
template <std::floating_point Real> struct tuple_element<1, mmcfilters::ComputedAttributeDataWithDelta<Real>> {
    /** @brief Defines the `type` alias used by the component. */
    using type = std::vector<Real>;
};

} // namespace std

namespace mmcfilters {

/**
 * @brief Tuple-like `std::get` overloads for computed attribute results.
 *
 * @param computed Flag controlling computed.
 * @return Requested tuple element.
 */
template <std::size_t I, std::floating_point Real> decltype(auto) get(ComputedAttributeData<Real>& computed) noexcept {
    if constexpr (I == 0) {
        return (computed.first);
    } else {
        return (computed.second);
    }
}

/**
 * @brief Returns the requested tuple element.
 *
 * @param computed Computed attribute result to transform.
 * @return The requested tuple element.
 */
template <std::size_t I, std::floating_point Real> decltype(auto) get(const ComputedAttributeData<Real>& computed) noexcept {
    if constexpr (I == 0) {
        return (computed.first);
    } else {
        return (computed.second);
    }
}

/**
 * @brief Returns the requested tuple element.
 *
 * @param computed Computed attribute result to transform.
 * @return The requested tuple element.
 */
template <std::size_t I, std::floating_point Real> decltype(auto) get(ComputedAttributeData<Real>&& computed) noexcept {
    if constexpr (I == 0) {
        return std::move(computed.first);
    } else {
        return std::move(computed.second);
    }
}

/**
 * @brief Returns the requested tuple element.
 *
 * @param computed Computed attribute result to transform.
 * @return The requested tuple element.
 */
template <std::size_t I, std::floating_point Real> decltype(auto) get(ComputedAttributeDataWithDelta<Real>& computed) noexcept {
    if constexpr (I == 0) {
        return (computed.first);
    } else {
        return (computed.second);
    }
}

/**
 * @brief Returns the requested tuple element.
 *
 * @param computed Computed attribute result to transform.
 * @return The requested tuple element.
 */
template <std::size_t I, std::floating_point Real> decltype(auto) get(const ComputedAttributeDataWithDelta<Real>& computed) noexcept {
    if constexpr (I == 0) {
        return (computed.first);
    } else {
        return (computed.second);
    }
}

/**
 * @brief Returns the requested tuple element.
 *
 * @param computed Computed attribute result to transform.
 * @return The requested tuple element.
 */
template <std::size_t I, std::floating_point Real> decltype(auto) get(ComputedAttributeDataWithDelta<Real>&& computed) noexcept {
    if constexpr (I == 0) {
        return std::move(computed.first);
    } else {
        return std::move(computed.second);
    }
}

} // namespace mmcfilters
