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
    NodeIdSpace nodeIdSpace = NodeIdSpace::MorphologicalTree;

    /**
     * @brief Takes ownership of a scalar-attribute layout and its flat buffer.
     *
     * @param attrNames Layout that maps attributes to per-node offsets.
     * @param buffer Attribute values stored in node-major order.
     * @param outputSpace Node-id domain used by the buffer rows.
     */
    ComputedAttributeData(AttributeNames attrNames, std::vector<Real> buffer, NodeIdSpace outputSpace = NodeIdSpace::MorphologicalTree)
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
 * @brief Owning result for one sampled node-attribute layout and buffer.
 *
 * @details
 * This result is used when one logical node attribute is sampled at several
 * ancestor/current/representative-descendant offsets.
 */
template <std::floating_point Real = float> struct [[nodiscard]] SampledNodeAttributeData {
    /// Sample layout used to interpret `second`.
    NodeAttributeSampleLayout first;

    /// Flat per-node, per-sample attribute buffer indexed through `first`.
    std::vector<Real> second;

    /// Node-id domain used by the rows of `second`.
    NodeIdSpace nodeIdSpace = NodeIdSpace::MorphologicalTree;

    /**
     * @brief Takes ownership of a sample layout and its flat buffer.
     *
     * @param attrNames Layout that maps `(attribute, sampleOffset)` keys to offsets.
     * @param buffer Attribute values stored in node-major order.
     * @param outputSpace Node-id domain used by the buffer rows.
     */
    SampledNodeAttributeData(NodeAttributeSampleLayout attrNames, std::vector<Real> buffer, NodeIdSpace outputSpace = NodeIdSpace::MorphologicalTree)
        : first(std::move(attrNames)), second(std::move(buffer)), nodeIdSpace(outputSpace) {}

    /**
     * @brief Disables copy construction.
     */
    SampledNodeAttributeData(const SampledNodeAttributeData&) = delete;

    /**
     * @brief Disables copy assignment.
     */
    SampledNodeAttributeData& operator=(const SampledNodeAttributeData&) = delete;

    /**
     * @brief Transfers ownership of the layout, buffer, and node-id-space marker.
     */
    SampledNodeAttributeData(SampledNodeAttributeData&&) noexcept = default;

    /**
     * @brief Disables move assignment.
     */
    SampledNodeAttributeData& operator=(SampledNodeAttributeData&&) = delete;

    /**
     * @brief Returns the mutable sampled-attribute layout.
     *
     * @return The mutable sampled-attribute layout.
     */
    NodeAttributeSampleLayout& attributeNames() noexcept { return this->first; }

    /**
     * @brief Returns the immutable sampled-attribute layout.
     *
     * @return The immutable sampled-attribute layout.
     */
    const NodeAttributeSampleLayout& attributeNames() const noexcept { return this->first; }

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

template <std::floating_point Real> struct tuple_size<mmcfilters::SampledNodeAttributeData<Real>> : integral_constant<std::size_t, 2> {};

/**
 * @brief Provides tuple-element metadata for structured binding support.
 *
 * @tparam Real Floating-point type stored in the sampled attribute buffer.
 */
template <std::floating_point Real> struct tuple_element<0, mmcfilters::SampledNodeAttributeData<Real>> {
    /** @brief Defines the `type` alias used by the component. */
    using type = mmcfilters::NodeAttributeSampleLayout;
};

/**
 * @brief Provides tuple-element metadata for structured binding support.
 *
 * @tparam Real Floating-point type stored in the sampled attribute buffer.
 */
template <std::floating_point Real> struct tuple_element<1, mmcfilters::SampledNodeAttributeData<Real>> {
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
template <std::size_t I, std::floating_point Real> decltype(auto) get(SampledNodeAttributeData<Real>& computed) noexcept {
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
template <std::size_t I, std::floating_point Real> decltype(auto) get(const SampledNodeAttributeData<Real>& computed) noexcept {
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
template <std::size_t I, std::floating_point Real> decltype(auto) get(SampledNodeAttributeData<Real>&& computed) noexcept {
    if constexpr (I == 0) {
        return std::move(computed.first);
    } else {
        return std::move(computed.second);
    }
}

} // namespace mmcfilters
