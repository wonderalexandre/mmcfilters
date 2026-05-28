#pragma once

#include "../AttributeTypes.hpp"
#include "../computers/AttributeComputerProtocol.hpp"

#include <algorithm>
#include <span>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace mmcfilters::detail {

/**
 * @brief Attribute-computer family identifiers used by orchestration code.
 *
 * @details
 * These ids are the runtime counterpart of the registered computer list. They are
 * intentionally small and local to the scheduler/backend; public users should
 * select attributes, not families.
 */
enum class AttributeFamily {
    Area,
    Volume,
    GrayLevelStats,
    MaxDist,
    BoundingBox,
    TreeTopology,
    CentralMoments,
    HuMoments,
    MomentDerived,
    Bitquad,
    ContourSide,
    Unsupported
};

template <class Computer>
struct AttributeComputerFamily;

template <>
struct AttributeComputerFamily<attributes::computers::AreaComputer> {
    static constexpr AttributeFamily value = AttributeFamily::Area;
};

template <>
struct AttributeComputerFamily<attributes::computers::VolumeComputer> {
    static constexpr AttributeFamily value = AttributeFamily::Volume;
};

template <>
struct AttributeComputerFamily<attributes::computers::GrayLevelStatsComputer> {
    static constexpr AttributeFamily value = AttributeFamily::GrayLevelStats;
};

template <>
struct AttributeComputerFamily<attributes::computers::MaxDistComputer> {
    static constexpr AttributeFamily value = AttributeFamily::MaxDist;
};

template <>
struct AttributeComputerFamily<attributes::computers::BoundingBoxComputer> {
    static constexpr AttributeFamily value = AttributeFamily::BoundingBox;
};

template <>
struct AttributeComputerFamily<attributes::computers::TreeTopologyComputer> {
    static constexpr AttributeFamily value = AttributeFamily::TreeTopology;
};

template <>
struct AttributeComputerFamily<attributes::computers::CentralMomentsComputer> {
    static constexpr AttributeFamily value = AttributeFamily::CentralMoments;
};

template <>
struct AttributeComputerFamily<attributes::computers::HuMomentsComputer> {
    static constexpr AttributeFamily value = AttributeFamily::HuMoments;
};

template <>
struct AttributeComputerFamily<attributes::computers::MomentBasedAttributeComputer> {
    static constexpr AttributeFamily value = AttributeFamily::MomentDerived;
};

template <>
struct AttributeComputerFamily<attributes::computers::BitquadAttributeComputer> {
    static constexpr AttributeFamily value = AttributeFamily::Bitquad;
};

template <>
struct AttributeComputerFamily<attributes::computers::ContourSideAttributeComputer> {
    static constexpr AttributeFamily value = AttributeFamily::ContourSide;
};

template <class Computer>
inline constexpr AttributeFamily attributeComputerFamilyV = AttributeComputerFamily<Computer>::value;

/**
 * @brief Tests membership in a scheduler-owned attribute list.
 */
inline bool containsScheduledAttribute(std::span<const Attribute> attributes, Attribute attribute) noexcept
{
    return std::find(attributes.begin(), attributes.end(), attribute) != attributes.end();
}

inline bool containsScheduledAttribute(const std::vector<Attribute>& attributes, Attribute attribute) noexcept
{
    return containsScheduledAttribute(std::span<const Attribute>(attributes), attribute);
}

/**
 * @brief Appends `attribute` while preserving first-seen request order.
 */
inline void appendScheduledAttributeOnce(std::vector<Attribute>& attributes, Attribute attribute)
{
    if (!containsScheduledAttribute(attributes, attribute)) {
        attributes.push_back(attribute);
    }
}

template <class Computer>
inline AttributeFamily familyForAttributeWithComputer(Attribute attribute) noexcept
{
    return attributes::computers::producesAttribute<Computer>(attribute)
        ? attributeComputerFamilyV<Computer>
        : AttributeFamily::Unsupported;
}

template <class Computer, class... Rest>
inline AttributeFamily familyForAttributeInComputers(Attribute attribute) noexcept
{
    const AttributeFamily family = familyForAttributeWithComputer<Computer>(attribute);
    if (family != AttributeFamily::Unsupported) {
        return family;
    }
    if constexpr (sizeof...(Rest) > 0) {
        return familyForAttributeInComputers<Rest...>(attribute);
    } else {
        return AttributeFamily::Unsupported;
    }
}

template <class Tuple>
struct AttributeFamilyLookup;

template <class... Computers>
struct AttributeFamilyLookup<std::tuple<Computers...>> {
    static AttributeFamily familyForAttribute(Attribute attribute) noexcept
    {
        return familyForAttributeInComputers<Computers...>(attribute);
    }
};

/**
 * @brief Returns the unique family declared as producer of `attribute`.
 */
inline AttributeFamily familyForAttribute(Attribute attribute) noexcept
{
    return AttributeFamilyLookup<attributes::computers::RegisteredAttributeComputers>::familyForAttribute(attribute);
}

template <attributes::computers::AttributeComputerDomain Domain, class Computer>
inline bool computerProducesAttributeInDomain(Attribute attribute) noexcept
{
    if constexpr (Computer::domain == Domain) {
        return attributes::computers::producesAttribute<Computer>(attribute);
    } else {
        return false;
    }
}

template <attributes::computers::AttributeComputerDomain Domain, class Computer, class... Rest>
inline bool attributeHasComputerDomainInComputers(Attribute attribute) noexcept
{
    if (computerProducesAttributeInDomain<Domain, Computer>(attribute)) {
        return true;
    }
    if constexpr (sizeof...(Rest) > 0) {
        return attributeHasComputerDomainInComputers<Domain, Rest...>(attribute);
    } else {
        return false;
    }
}

template <attributes::computers::AttributeComputerDomain Domain, class Tuple>
struct AttributeComputerDomainLookup;

template <attributes::computers::AttributeComputerDomain Domain, class... Computers>
struct AttributeComputerDomainLookup<Domain, std::tuple<Computers...>> {
    static bool contains(Attribute attribute) noexcept
    {
        return attributeHasComputerDomainInComputers<Domain, Computers...>(attribute);
    }
};

template <attributes::computers::AttributeComputerDomain Domain>
inline bool attributeHasComputerDomain(Attribute attribute) noexcept
{
    return AttributeComputerDomainLookup<Domain, attributes::computers::RegisteredAttributeComputers>::contains(attribute);
}

/**
 * @brief Attribute-level dependencies required to materialize one descriptor.
 *
 * @details
 * This is intentionally more precise than family-level traits. For example,
 * `BoundingBoxComputer` can consume `AREA`, but only `RECTANGULARITY` needs it.
 *
 * The returned order is dependency order for the current descriptor only. The
 * recursive closure builder expands nested dependencies before the consumer is
 * appended to the materialization list.
 */
inline std::vector<Attribute> dependenciesForAttribute(Attribute attribute)
{
    switch (attribute) {
        case RECTANGULARITY:
        case RELATIVE_VOLUME:
            return {AREA};

        case MEAN_LEVEL:
        case VARIANCE_LEVEL:
            return {VOLUME, AREA};

        case HU_MOMENT_1:
        case HU_MOMENT_2:
        case HU_MOMENT_3:
        case HU_MOMENT_4:
        case HU_MOMENT_5:
        case HU_MOMENT_6:
        case HU_MOMENT_7:
            return {
                AREA,
                CENTRAL_MOMENT_20,
                CENTRAL_MOMENT_02,
                CENTRAL_MOMENT_11,
                CENTRAL_MOMENT_30,
                CENTRAL_MOMENT_03,
                CENTRAL_MOMENT_21,
                CENTRAL_MOMENT_12};

        case INERTIA:
        case COMPACTNESS:
        case ECCENTRICITY:
        case LENGTH_MAJOR_AXIS:
        case LENGTH_MINOR_AXIS:
        case AXIS_ORIENTATION:
        case CIRCULARITY:
            return {
                AREA,
                CENTRAL_MOMENT_20,
                CENTRAL_MOMENT_02,
                CENTRAL_MOMENT_11};

        default:
            return {};
    }
}

/**
 * @brief Tests whether one descriptor directly consumes another descriptor.
 */
inline bool attributeRequiresDependency(Attribute attribute, Attribute dependency)
{
    const std::vector<Attribute> dependencies = dependenciesForAttribute(attribute);
    return containsScheduledAttribute(dependencies, dependency);
}

/**
 * @brief Tests whether any descriptor in `attributes` consumes `dependency`.
 */
inline bool anyAttributeRequiresDependency(std::span<const Attribute> attributes, Attribute dependency)
{
    return std::any_of(
        attributes.begin(),
        attributes.end(),
        [&](Attribute attribute) {
            return attributeRequiresDependency(attribute, dependency);
        });
}

/**
 * @brief Adds `attribute` and its recursive dependencies to `closure`.
 *
 * @details
 * Dependencies are appended before their consumer, so the final closure can be
 * traversed in materialization order. `visiting` is a recursion stack used only
 * for cycle detection.
 */
inline void appendDependencyClosure(
    Attribute attribute,
    std::vector<Attribute>& closure,
    std::vector<Attribute>& visiting)
{
    if (containsScheduledAttribute(closure, attribute)) {
        return;
    }
    if (containsScheduledAttribute(visiting, attribute)) {
        throw std::logic_error("Attribute dependency cycle detected while building the computation plan.");
    }

    visiting.push_back(attribute);
    for (Attribute dependency : dependenciesForAttribute(attribute)) {
        appendDependencyClosure(dependency, closure, visiting);
    }
    visiting.pop_back();
    appendScheduledAttributeOnce(closure, attribute);
}

/**
 * @brief Builds the dependency-closed materialization order for a request.
 */
inline std::vector<Attribute> expandDependencyClosure(std::span<const Attribute> requestedAttributes)
{
    std::vector<Attribute> closure;
    std::vector<Attribute> visiting;
    for (Attribute attribute : requestedAttributes) {
        appendDependencyClosure(attribute, closure, visiting);
    }
    return closure;
}

/**
 * @brief Dependency-expanded request plan grouped by computer family.
 *
 * @details
 * `requestedAttributes` is the public scalar request after group expansion and
 * deduplication. `materializedAttributes` also contains hidden dependencies in
 * dependency-first order. `hiddenDependencyAttributes` is the subset that must
 * be computed internally but omitted from the public result layout.
 */
struct AttributeComputationPlan {
    /// Public requested attributes, without hidden dependencies.
    std::vector<Attribute> requestedAttributes;

    /// Requested attributes plus recursive dependencies in compute order.
    std::vector<Attribute> materializedAttributes;

    /// Materialized attributes that are not part of the public request.
    std::vector<Attribute> hiddenDependencyAttributes;

    /**
     * @brief Tests whether `attribute` belongs to the public request.
     */
    [[nodiscard]] bool requests(Attribute attribute) const noexcept
    {
        return containsScheduledAttribute(requestedAttributes, attribute);
    }

    /**
     * @brief Tests whether `attribute` will be computed by the plan.
     */
    [[nodiscard]] bool materializes(Attribute attribute) const noexcept
    {
        return containsScheduledAttribute(materializedAttributes, attribute);
    }

    /**
     * @brief Tests whether `attribute` is an internal-only dependency.
     */
    [[nodiscard]] bool hides(Attribute attribute) const noexcept
    {
        return containsScheduledAttribute(hiddenDependencyAttributes, attribute);
    }

    /**
     * @brief Returns public requested attributes produced by `family`.
     */
    [[nodiscard]] std::vector<Attribute> requestedForFamily(AttributeFamily family) const
    {
        return attributesForFamily(requestedAttributes, family);
    }

    /**
     * @brief Returns requested and hidden attributes produced by `family`.
     */
    [[nodiscard]] std::vector<Attribute> materializedForFamily(AttributeFamily family) const
    {
        return attributesForFamily(materializedAttributes, family);
    }

private:
    [[nodiscard]] static std::vector<Attribute> attributesForFamily(
        const std::vector<Attribute>& attributes,
        AttributeFamily family)
    {
        std::vector<Attribute> filtered;
        for (Attribute attribute : attributes) {
            if (familyForAttribute(attribute) == family) {
                filtered.push_back(attribute);
            }
        }
        return filtered;
    }
};

/**
 * @brief Creates a computation plan from scalar attributes already expanded from groups.
 *
 * @details
 * Duplicate public requests are ignored after the first occurrence. Dependencies
 * are added recursively and marked hidden unless they were requested directly.
 */
inline AttributeComputationPlan makeAttributeComputationPlan(std::span<const Attribute> requestedAttributes)
{
    AttributeComputationPlan plan;
    for (Attribute attribute : requestedAttributes) {
        appendScheduledAttributeOnce(plan.requestedAttributes, attribute);
    }
    plan.materializedAttributes = expandDependencyClosure(std::span<const Attribute>(plan.requestedAttributes));

    for (Attribute attribute : plan.materializedAttributes) {
        if (!plan.requests(attribute)) {
            plan.hiddenDependencyAttributes.push_back(attribute);
        }
    }

    return plan;
}

} // namespace mmcfilters::detail
