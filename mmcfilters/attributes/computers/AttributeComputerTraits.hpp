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
 * A specialization is the declarative contract between one concrete computer
 * and the attribute scheduler. It states which scalar descriptors the family
 * can produce and which other descriptors may be needed as dependencies. The
 * execution domain determines whether a computer receives an altitude-aware
 * context.
 *
 * `requiredAttributes` is a family-level superset. The scheduler refines this
 * through per-attribute dependency rules, because only some descriptors in a
 * family may consume a dependency.
 */
template <class Computer>
struct AttributeComputerTraits;

/**
 * @brief Verifies that a concrete computer has a complete trait contract.
 */
template <class Computer>
concept AttributeComputerWithTraits =
    requires {
        { AttributeComputerTraits<Computer>::familyName } -> std::convertible_to<std::string_view>;
        { AttributeComputerTraits<Computer>::domain } -> std::convertible_to<AttributeComputerDomain>;
        AttributeComputerTraits<Computer>::producedAttributes;
        AttributeComputerTraits<Computer>::requiredAttributes;
    } &&
    IsAttributeArrayV<decltype(AttributeComputerTraits<Computer>::producedAttributes)> &&
    IsAttributeArrayV<decltype(AttributeComputerTraits<Computer>::requiredAttributes)>;

/**
 * @brief Verifies the runtime metadata method exposed by a computer instance.
 */
template <class Computer>
concept AttributeComputerInstance =
    std::default_initializable<Computer> &&
    requires(const Computer& computer) {
        { computer.attributes() } -> std::same_as<std::vector<Attribute>>;
    };

/**
 * @brief Computer concept for families that do not read node altitudes.
 */
template <class Computer, class Real = float>
concept TopologyAttributeComputer =
    std::floating_point<Real> &&
    AttributeComputerWithTraits<Computer> &&
    AttributeComputerInstance<Computer> &&
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
    AttributeComputerInstance<Computer> &&
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
    const auto& attributes = AttributeComputerTraits<Computer>::producedAttributes;
    return std::find(attributes.begin(), attributes.end(), attribute) != attributes.end();
}

/**
 * @brief Tests whether `attribute` belongs to the family dependency superset.
 */
template <AttributeComputerWithTraits Computer>
[[nodiscard]] constexpr bool requiresAttribute(Attribute attribute) noexcept
{
    const auto& attributes = AttributeComputerTraits<Computer>::requiredAttributes;
    return std::find(attributes.begin(), attributes.end(), attribute) != attributes.end();
}

/**
 * @brief Number of scalar descriptors declared by the family.
 */
template <AttributeComputerWithTraits Computer>
[[nodiscard]] constexpr std::size_t numProducedAttributes() noexcept
{
    return AttributeComputerTraits<Computer>::producedAttributes.size();
}

/**
 * @brief Number of dependency descriptors declared by the family.
 */
template <AttributeComputerWithTraits Computer>
[[nodiscard]] constexpr std::size_t numRequiredAttributes() noexcept
{
    return AttributeComputerTraits<Computer>::requiredAttributes.size();
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
    /// Attributes produced by the family.
    inline static constexpr std::array<Attribute, 1> producedAttributes{AREA};
    /// Attributes that may be required by the family.
    inline static constexpr std::array<Attribute, 0> requiredAttributes{};
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
    /// Attributes produced by the family.
    inline static constexpr std::array<Attribute, 9> producedAttributes{
        BOX_WIDTH,
        BOX_HEIGHT,
        DIAGONAL_LENGTH,
        RECTANGULARITY,
        RATIO_WH,
        BOX_COL_MIN,
        BOX_COL_MAX,
        BOX_ROW_MIN,
        BOX_ROW_MAX};
    /// Attributes that may be required by the family.
    inline static constexpr std::array<Attribute, 1> requiredAttributes{AREA};
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
    /// Attributes produced by the family.
    inline static constexpr std::array<Attribute, 11> producedAttributes{
        HEIGHT_NODE,
        DEPTH_NODE,
        IS_LEAF_NODE,
        IS_ROOT_NODE,
        NUM_CHILDREN_NODE,
        NUM_SIBLINGS_NODE,
        NUM_DESCENDANTS_NODE,
        NUM_LEAF_DESCENDANTS_NODE,
        LEAF_RATIO_NODE,
        BALANCE_NODE,
        AVG_CHILD_HEIGHT_NODE};
    /// Attributes that may be required by the family.
    inline static constexpr std::array<Attribute, 0> requiredAttributes{};
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
    /// Attributes produced by the family.
    inline static constexpr std::array<Attribute, 7> producedAttributes{
        CENTRAL_MOMENT_20,
        CENTRAL_MOMENT_02,
        CENTRAL_MOMENT_11,
        CENTRAL_MOMENT_30,
        CENTRAL_MOMENT_03,
        CENTRAL_MOMENT_21,
        CENTRAL_MOMENT_12};
    /// Attributes that may be required by the family.
    inline static constexpr std::array<Attribute, 0> requiredAttributes{};
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
    /// Attributes produced by the family.
    inline static constexpr std::array<Attribute, 7> producedAttributes{
        HU_MOMENT_1,
        HU_MOMENT_2,
        HU_MOMENT_3,
        HU_MOMENT_4,
        HU_MOMENT_5,
        HU_MOMENT_6,
        HU_MOMENT_7};
    /// Attributes that may be required by the family.
    inline static constexpr std::array<Attribute, 8> requiredAttributes{
        AREA,
        CENTRAL_MOMENT_20,
        CENTRAL_MOMENT_02,
        CENTRAL_MOMENT_11,
        CENTRAL_MOMENT_30,
        CENTRAL_MOMENT_03,
        CENTRAL_MOMENT_21,
        CENTRAL_MOMENT_12};
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
    /// Attributes produced by the family.
    inline static constexpr std::array<Attribute, 7> producedAttributes{
        INERTIA,
        COMPACTNESS,
        ECCENTRICITY,
        LENGTH_MAJOR_AXIS,
        LENGTH_MINOR_AXIS,
        AXIS_ORIENTATION,
        CIRCULARITY};
    /// Attributes that may be required by the family.
    inline static constexpr std::array<Attribute, 4> requiredAttributes{
        AREA,
        CENTRAL_MOMENT_20,
        CENTRAL_MOMENT_02,
        CENTRAL_MOMENT_11};
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
    /// Attributes produced by the family.
    inline static constexpr std::array<Attribute, 9> producedAttributes{
        BITQUADS_AREA,
        BITQUADS_NUMBER_EULER,
        BITQUADS_NUMBER_HOLES,
        BITQUADS_PERIMETER,
        BITQUADS_PERIMETER_CONTINUOUS,
        BITQUADS_CIRCULARITY,
        BITQUADS_PERIMETER_AVERAGE,
        BITQUADS_LENGTH_AVERAGE,
        BITQUADS_WIDTH_AVERAGE};
    /// Attributes that may be required by the family.
    inline static constexpr std::array<Attribute, 0> requiredAttributes{};
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
    /// Attributes produced by the family.
    inline static constexpr std::array<Attribute, 6> producedAttributes{
        CONTOUR_PIXELS,
        CONTOUR_PERIMETER,
        CONTOUR_SIDE_NORTH,
        CONTOUR_SIDE_WEST,
        CONTOUR_SIDE_EAST,
        CONTOUR_SIDE_SOUTH};
    /// Attributes that may be required by the family.
    inline static constexpr std::array<Attribute, 0> requiredAttributes{};
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
    /// Attributes produced by the family.
    inline static constexpr std::array<Attribute, 2> producedAttributes{
        VOLUME,
        RELATIVE_VOLUME};
    /// Attributes that may be required by the family.
    inline static constexpr std::array<Attribute, 1> requiredAttributes{AREA};
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
    /// Attributes produced by the family.
    inline static constexpr std::array<Attribute, 4> producedAttributes{
        LEVEL,
        MEAN_LEVEL,
        VARIANCE_LEVEL,
        GRAY_HEIGHT};
    /// Attributes that may be required by the family.
    inline static constexpr std::array<Attribute, 2> requiredAttributes{
        AREA,
        VOLUME};
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
    /// Attributes produced by the family.
    inline static constexpr std::array<Attribute, 1> producedAttributes{MAX_DIST};
    /// Attributes that may be required by the family.
    inline static constexpr std::array<Attribute, 0> requiredAttributes{};
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
