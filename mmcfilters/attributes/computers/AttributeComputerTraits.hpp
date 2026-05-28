#pragma once

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
 * @brief Stable execution domain declared by an attribute computer family.
 *
 * @details
 * The domain is part of the orchestration contract, not a mathematical
 * classification of the descriptor itself. `Topology` families can run from a
 * `MorphologicalTree` alone. `Altitude` families need a dense altitude span in
 * addition to topology.
 */
enum class AttributeComputerDomain {
    /// Computes from support, image geometry, or parent/child topology only.
    Topology,

    /// Computes from topology plus a typed altitude span.
    Altitude
};

template <class T>
struct IsAttributeArray : std::false_type {};

template <std::size_t N>
struct IsAttributeArray<std::array<Attribute, N>> : std::true_type {};

template <class T>
inline constexpr bool IsAttributeArrayV = IsAttributeArray<std::remove_cvref_t<T>>::value;

/**
 * @brief Compile-time metadata for an attribute computer family.
 *
 * @details
 * A specialization declares metadata that is not already owned by the
 * computer class: its diagnostic family name and execution domain. Produced
 * attributes live on `Computer::producedAttributes`; attribute-level
 * dependency rules live in the scheduler.
 */
template <class Computer>
struct AttributeComputerTraits;

/**
 * @brief Verifies that a computer declares its canonical output list.
 */
template <class Computer>
concept AttributeComputerProducesAttributes =
    requires {
        Computer::producedAttributes;
    } &&
    IsAttributeArrayV<decltype(Computer::producedAttributes)>;

/**
 * @brief Verifies that a concrete computer has a complete trait contract.
 */
template <class Computer>
concept AttributeComputerWithTraits =
    AttributeComputerProducesAttributes<Computer> &&
    requires {
        { AttributeComputerTraits<Computer>::familyName } -> std::convertible_to<std::string_view>;
        { AttributeComputerTraits<Computer>::domain } -> std::convertible_to<AttributeComputerDomain>;
    };

/**
 * @brief Computer concept for families that do not read node altitudes.
 */
template <class Computer, class Real = float>
concept TopologyAttributeComputer =
    std::floating_point<Real> &&
    AttributeComputerWithTraits<Computer> &&
    AttributeComputerTraits<Computer>::domain == AttributeComputerDomain::Topology &&
    requires(
        const AttributeComputeContext<Real>& context,
        const UnitAttributeComputeContext<Real>& unitContext) {
        { Computer::compute(context) } -> std::same_as<void>;
        { Computer::computeUnitRows(unitContext) } -> std::same_as<void>;
    };

/**
 * @brief Computer concept for families that require a typed altitude span.
 */
template <class Computer, class Real = float, class T = std::uint8_t>
concept AltitudeAttributeComputer =
    std::floating_point<Real> &&
    AltitudeValue<T> &&
    AttributeComputerWithTraits<Computer> &&
    AttributeComputerTraits<Computer>::domain == AttributeComputerDomain::Altitude &&
    requires(
        const AltitudeAttributeComputeContext<Real, T>& context,
        const AltitudeUnitAttributeComputeContext<Real, T>& unitContext) {
        { Computer::compute(context) } -> std::same_as<void>;
        { Computer::computeUnitRows(unitContext) } -> std::same_as<void>;
    };

/**
 * @brief Tests whether `Computer` is the declared producer of `attribute`.
 */
template <AttributeComputerWithTraits Computer>
[[nodiscard]] constexpr bool producesAttribute(Attribute attribute) noexcept
{
    const auto& attributes = Computer::producedAttributes;
    return std::find(attributes.begin(), attributes.end(), attribute) != attributes.end();
}

/**
 * @brief Number of scalar descriptors declared by the family.
 */
template <AttributeComputerWithTraits Computer>
[[nodiscard]] constexpr std::size_t numProducedAttributes() noexcept
{
    return Computer::producedAttributes.size();
}

/**
 * @brief Materializes a computer's canonical output list as a runtime vector.
 */
template <AttributeComputerWithTraits Computer>
[[nodiscard]] std::vector<Attribute> runtimeProducedAttributes()
{
    const auto& attributes = Computer::producedAttributes;
    return {attributes.begin(), attributes.end()};
}

/**
 * @brief Scheduler metadata for support-area attributes.
 */
template <>
struct AttributeComputerTraits<AreaComputer> {
    /// Family name used in dependency-plan diagnostics.
    static constexpr std::string_view familyName = "area";
    /// Execution domain required by the computer.
    static constexpr AttributeComputerDomain domain = AttributeComputerDomain::Topology;
};

/**
 * @brief Scheduler metadata for bounding-box shape attributes.
 */
template <>
struct AttributeComputerTraits<BoundingBoxComputer> {
    /// Family name used in dependency-plan diagnostics.
    static constexpr std::string_view familyName = "bounding-box";
    /// Execution domain required by the computer.
    static constexpr AttributeComputerDomain domain = AttributeComputerDomain::Topology;
};

/**
 * @brief Scheduler metadata for pure tree-topology attributes.
 */
template <>
struct AttributeComputerTraits<TreeTopologyComputer> {
    /// Family name used in dependency-plan diagnostics.
    static constexpr std::string_view familyName = "tree-topology";
    /// Execution domain required by the computer.
    static constexpr AttributeComputerDomain domain = AttributeComputerDomain::Topology;
};

/**
 * @brief Scheduler metadata for central-moment attributes.
 */
template <>
struct AttributeComputerTraits<CentralMomentsComputer> {
    /// Family name used in dependency-plan diagnostics.
    static constexpr std::string_view familyName = "central-moments";
    /// Execution domain required by the computer.
    static constexpr AttributeComputerDomain domain = AttributeComputerDomain::Topology;
};

/**
 * @brief Scheduler metadata for Hu invariant moment attributes.
 */
template <>
struct AttributeComputerTraits<HuMomentsComputer> {
    /// Family name used in dependency-plan diagnostics.
    static constexpr std::string_view familyName = "hu-moments";
    /// Execution domain required by the computer.
    static constexpr AttributeComputerDomain domain = AttributeComputerDomain::Topology;
};

/**
 * @brief Scheduler metadata for descriptors derived from second-order moments.
 */
template <>
struct AttributeComputerTraits<MomentBasedAttributeComputer> {
    /// Family name used in dependency-plan diagnostics.
    static constexpr std::string_view familyName = "moment-derived";
    /// Execution domain required by the computer.
    static constexpr AttributeComputerDomain domain = AttributeComputerDomain::Topology;
};

/**
 * @brief Scheduler metadata for bitquad-based shape attributes.
 */
template <>
struct AttributeComputerTraits<BitquadAttributeComputer> {
    /// Family name used in dependency-plan diagnostics.
    static constexpr std::string_view familyName = "bitquad";
    /// Execution domain required by the computer.
    static constexpr AttributeComputerDomain domain = AttributeComputerDomain::Topology;
};

/**
 * @brief Scheduler metadata for contour-side counting attributes.
 */
template <>
struct AttributeComputerTraits<ContourSideAttributeComputer> {
    /// Family name used in dependency-plan diagnostics.
    static constexpr std::string_view familyName = "contour-side";
    /// Execution domain required by the computer.
    static constexpr AttributeComputerDomain domain = AttributeComputerDomain::Topology;
};

/**
 * @brief Scheduler metadata for altitude-volume attributes.
 */
template <>
struct AttributeComputerTraits<VolumeComputer> {
    /// Family name used in dependency-plan diagnostics.
    static constexpr std::string_view familyName = "volume";
    /// Execution domain required by the computer.
    static constexpr AttributeComputerDomain domain = AttributeComputerDomain::Altitude;
};

/**
 * @brief Scheduler metadata for gray-level statistics.
 */
template <>
struct AttributeComputerTraits<GrayLevelStatsComputer> {
    /// Family name used in dependency-plan diagnostics.
    static constexpr std::string_view familyName = "gray-level-stats";
    /// Execution domain required by the computer.
    static constexpr AttributeComputerDomain domain = AttributeComputerDomain::Altitude;
};

/**
 * @brief Scheduler metadata for maximum-distance attributes.
 */
template <>
struct AttributeComputerTraits<MaxDistComputer> {
    /// Family name used in dependency-plan diagnostics.
    static constexpr std::string_view familyName = "max-dist";
    /// Execution domain required by the computer.
    static constexpr AttributeComputerDomain domain = AttributeComputerDomain::Altitude;
};

/**
 * @brief Canonical list of topology/support families known to the backend.
 */
using TopologyAttributeComputers = std::tuple<
    AreaComputer,
    BoundingBoxComputer,
    TreeTopologyComputer,
    CentralMomentsComputer,
    HuMomentsComputer,
    MomentBasedAttributeComputer,
    BitquadAttributeComputer,
    ContourSideAttributeComputer>;

/**
 * @brief Canonical list of altitude-aware families known to the pipeline.
 */
using AltitudeAttributeComputers = std::tuple<
    VolumeComputer,
    GrayLevelStatsComputer,
    MaxDistComputer>;

/**
 * @brief Canonical list used by contract tests to cover every public family.
 */
using RegisteredAttributeComputers = std::tuple<
    AreaComputer,
    BoundingBoxComputer,
    TreeTopologyComputer,
    CentralMomentsComputer,
    HuMomentsComputer,
    MomentBasedAttributeComputer,
    BitquadAttributeComputer,
    ContourSideAttributeComputer,
    VolumeComputer,
    GrayLevelStatsComputer,
    MaxDistComputer>;

} // namespace mmcfilters::attributes::computers
