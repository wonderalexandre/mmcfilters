#pragma once

#include "AttributeComputerDomain.hpp"
#include "AttributeComputerFamily.hpp"
#include "AreaComputer.hpp"
#include "BitquadAttributeComputer.hpp"
#include "BoundingBoxComputer.hpp"
#include "ContourSideAttributeComputer.hpp"
#include "GrayLevelStatsComputer.hpp"
#include "MaxDistComputer.hpp"
#include "MomentBasedAttributeComputer.hpp"
#include "TreeTopologyComputer.hpp"
#include "VolumeComputer.hpp"

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <string_view>
#include <type_traits>
#include <tuple>
#include <vector>

namespace mmcfilters::attributes::computers {

/**
 * @file AttributeComputerRegistry.hpp
 * @brief Concepts and canonical family lists for concrete attribute computers.
 */

namespace detail {

/**
 * @brief Fallback trait for types that are not canonical attribute arrays.
 */
template <class T> struct IsAttributeArray : std::false_type {};

/**
 * @brief Matches the fixed-size attribute arrays used by computer metadata.
 */
template <std::size_t N> struct IsAttributeArray<std::array<Attribute, N>> : std::true_type {};

/**
 * @brief Convenience value for checking whether a type is a canonical attribute array.
 */
template <class T> inline constexpr bool IsAttributeArrayV = IsAttributeArray<std::remove_cvref_t<T>>::value;

} // namespace detail

/**
 * @brief Verifies that a computer declares its canonical output list.
 */
template <class Computer>
concept AttributeComputerProducesAttributes = requires { Computer::producedAttributes; } && detail::IsAttributeArrayV<decltype(Computer::producedAttributes)>;

/**
 * @brief Verifies the complete static protocol of one attribute computer.
 */
template <class Computer>
concept AttributeComputer = AttributeComputerProducesAttributes<Computer> && requires {
    { Computer::familyName } -> std::convertible_to<std::string_view>;
    { Computer::family } -> std::convertible_to<AttributeComputerFamily>;
    { Computer::domain } -> std::convertible_to<AttributeComputerDomain>;
};

/**
 * @brief Computer concept for families that do not read node altitudes.
 */
template <class Computer, class Real = float>
concept TopologyAttributeComputer = std::floating_point<Real> && AttributeComputer<Computer> && Computer::domain == AttributeComputerDomain::Topology &&
                                    requires(const AttributeComputeContext<Real>& context, const UnitAttributeComputeContext<Real>& unitContext) {
                                        { Computer::compute(context) } -> std::same_as<void>;
                                        { Computer::computeUnitRows(unitContext) } -> std::same_as<void>;
                                    };

/**
 * @brief Computer concept for families that require a typed altitude span.
 */
template <class Computer, class Real = float, class T = std::uint8_t>
concept AltitudeAttributeComputer =
    std::floating_point<Real> && AltitudeValue<T> && AttributeComputer<Computer> && Computer::domain == AttributeComputerDomain::Altitude &&
    requires(const AltitudeAttributeComputeContext<Real, T>& context, const AltitudeUnitAttributeComputeContext<Real, T>& unitContext) {
        { Computer::compute(context) } -> std::same_as<void>;
        { Computer::computeUnitRows(unitContext) } -> std::same_as<void>;
    };

/**
 * @brief Tests whether `Computer` is the declared producer of `attribute`.
 *
 * @param attribute Attribute requested by the operation.
 * @return True if Computer is the declared producer of attribute; otherwise false.
 */
template <AttributeComputer Computer> [[nodiscard]] constexpr bool producesAttribute(Attribute attribute) noexcept {
    const auto& attributes = Computer::producedAttributes;
    return std::find(attributes.begin(), attributes.end(), attribute) != attributes.end();
}

/**
 * @brief Number of scalar descriptors declared by the family.
 *
 * @return The number of scalar descriptors declared by the family.
 */
template <AttributeComputer Computer> [[nodiscard]] constexpr std::size_t numProducedAttributes() noexcept { return Computer::producedAttributes.size(); }

/**
 * @brief Materializes a computer's canonical output list as a runtime vector.
 *
 * @return The materialized computer's canonical output list as a runtime vector.
 */
template <AttributeComputer Computer> [[nodiscard]] std::vector<Attribute> runtimeProducedAttributes() {
    const auto& attributes = Computer::producedAttributes;
    return {attributes.begin(), attributes.end()};
}

/**
 * @brief Canonical list of topology/support families known to the backend.
 */
using TopologyAttributeComputers = std::tuple<AreaComputer, BoundingBoxComputer, TreeTopologyComputer, CentralMomentsComputer, HuMomentsComputer,
                                              MomentBasedAttributeComputer, BitquadAttributeComputer, ContourSideAttributeComputer>;

/**
 * @brief Canonical list of altitude-aware families known to the pipeline.
 */
using AltitudeAttributeComputers = std::tuple<VolumeComputer, GrayLevelStatsComputer, MaxDistComputer>;

/**
 * @brief Canonical list used by contract tests to cover every public family.
 */
using RegisteredAttributeComputers =
    std::tuple<AreaComputer, BoundingBoxComputer, TreeTopologyComputer, CentralMomentsComputer, HuMomentsComputer, MomentBasedAttributeComputer,
               BitquadAttributeComputer, ContourSideAttributeComputer, VolumeComputer, GrayLevelStatsComputer, MaxDistComputer>;

} // namespace mmcfilters::attributes::computers
