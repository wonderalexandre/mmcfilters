#pragma once

#include "../AttributeTypes.hpp"
#include "../computers/AttributeComputerRegistry.hpp"

#include <algorithm>
#include <span>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace mmcfilters::detail {

/**
 * @brief Tests membership in a scheduler-owned attribute list.
 *
 * @param attributes Attributes requested by the operation.
 * @param attribute Attribute requested by the operation.
 * @return True when the documented condition holds; otherwise false.
 */
inline bool containsScheduledAttribute(std::span<const Attribute> attributes, Attribute attribute) noexcept {
    return std::find(attributes.begin(), attributes.end(), attribute) != attributes.end();
}

/**
 * @brief Tests whether scheduled attribute holds.
 *
 * @param attributes Attributes requested by the operation.
 * @param attribute Attribute requested by the operation.
 * @return True when scheduled attribute; otherwise false.
 */
inline bool containsScheduledAttribute(const std::vector<Attribute>& attributes, Attribute attribute) noexcept {
    return containsScheduledAttribute(std::span<const Attribute>(attributes), attribute);
}

/**
 * @brief Appends `attribute` while preserving first-seen request order.
 *
 * @param attributes Attributes requested by the operation.
 * @param attribute Attribute requested by the operation.
 */
inline void appendScheduledAttributeOnce(std::vector<Attribute>& attributes, Attribute attribute) {
    if (!containsScheduledAttribute(attributes, attribute)) {
        attributes.push_back(attribute);
    }
}

/**
 * @brief Returns the family associated with an attribute for one computer type.
 *
 * @param attribute Attribute requested by the operation.
 * @return Matching computer family, or `None` when unsupported.
 */
template <class Computer> inline attributes::computers::AttributeComputerFamily familyForAttributeWithComputer(Attribute attribute) noexcept {
    return attributes::computers::producesAttribute<Computer>(attribute) ? Computer::family : attributes::computers::AttributeComputerFamily::Unsupported;
}

/**
 * @brief Finds the family associated with an attribute across computer types.
 *
 * @param attribute Attribute requested by the operation.
 * @return First matching computer family, or `None` when unsupported.
 */
template <class Computer, class... Rest> inline attributes::computers::AttributeComputerFamily familyForAttributeInComputers(Attribute attribute) noexcept {
    const attributes::computers::AttributeComputerFamily family = familyForAttributeWithComputer<Computer>(attribute);
    if (family != attributes::computers::AttributeComputerFamily::Unsupported) {
        return family;
    }
    if constexpr (sizeof...(Rest) > 0) {
        return familyForAttributeInComputers<Rest...>(attribute);
    } else {
        return attributes::computers::AttributeComputerFamily::Unsupported;
    }
}

template <class Tuple> struct AttributeFamilyLookup;

/**
 * @brief Resolves attribute-computer families from a registered computer tuple.
 *
 * @tparam Computers Registered attribute-computer types searched by the lookup.
 */
template <class... Computers> struct AttributeFamilyLookup<std::tuple<Computers...>> {
    /**
     * @brief Returns the registered computer family for an attribute.
     *
     * @param attribute Attribute requested by the operation.
     * @return Registered computer family, or `None` when unsupported.
     */
    static attributes::computers::AttributeComputerFamily familyForAttribute(Attribute attribute) noexcept {
        return familyForAttributeInComputers<Computers...>(attribute);
    }
};

/**
 * @brief Returns the unique family declared as producer of `attribute`.
 *
 * @param attribute Attribute requested by the operation.
 * @return The unique family declared as producer of attribute.
 */
inline attributes::computers::AttributeComputerFamily familyForAttribute(Attribute attribute) noexcept {
    return AttributeFamilyLookup<attributes::computers::RegisteredAttributeComputers>::familyForAttribute(attribute);
}

/**
 * @brief Checks whether a computer produces an attribute in the requested domain.
 *
 * @param attribute Attribute requested by the operation.
 * @return True when the computer produces the attribute in the requested domain.
 */
template <attributes::computers::AttributeComputerDomain Domain, class Computer> inline bool computerProducesAttributeInDomain(Attribute attribute) noexcept {
    if constexpr (Computer::domain == Domain) {
        return attributes::computers::producesAttribute<Computer>(attribute);
    } else {
        return false;
    }
}

/**
 * @brief Tests whether has computer domain in computers holds.
 *
 * @param attribute Attribute requested by the operation.
 * @return True when has computer domain in computers; otherwise false.
 */
template <attributes::computers::AttributeComputerDomain Domain, class Computer, class... Rest>
inline bool attributeHasComputerDomainInComputers(Attribute attribute) noexcept {
    if (computerProducesAttributeInDomain<Domain, Computer>(attribute)) {
        return true;
    }
    if constexpr (sizeof...(Rest) > 0) {
        return attributeHasComputerDomainInComputers<Domain, Rest...>(attribute);
    } else {
        return false;
    }
}

template <attributes::computers::AttributeComputerDomain Domain, class Tuple> struct AttributeComputerDomainLookup;

/**
 * @brief Tests a registered computer tuple for one attribute-computation domain.
 *
 * @tparam Domain Attribute-computation domain required by the lookup.
 * @tparam Computers Registered attribute-computer types searched by the lookup.
 */
template <attributes::computers::AttributeComputerDomain Domain, class... Computers> struct AttributeComputerDomainLookup<Domain, std::tuple<Computers...>> {
    /**
     * @brief Tests whether contains holds.
     *
     * @param attribute Attribute requested by the operation.
     * @return True when contains; otherwise false.
     */
    static bool contains(Attribute attribute) noexcept { return attributeHasComputerDomainInComputers<Domain, Computers...>(attribute); }
};

/**
 * @brief Tests whether has computer domain holds.
 *
 * @param attribute Attribute requested by the operation.
 * @return True when has computer domain; otherwise false.
 */
template <attributes::computers::AttributeComputerDomain Domain> inline bool attributeHasComputerDomain(Attribute attribute) noexcept {
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
 *
 * @param attribute Attribute requested by the operation.
 * @return Values produced by the operation.
 */
inline std::vector<Attribute> dependenciesForAttribute(Attribute attribute) {
    switch (attribute) {
    case Rectangularity:
    case RelativeVolume:
        return {Area};

    case MeanGrayLevel:
    case GrayLevelVariance:
        return {Volume, Area};

    case HuMoment1:
    case HuMoment2:
    case HuMoment3:
    case HuMoment4:
    case HuMoment5:
    case HuMoment6:
    case HuMoment7:
        return {Area, CentralMoment20, CentralMoment02, CentralMoment11, CentralMoment30, CentralMoment03, CentralMoment21, CentralMoment12};

    case Inertia:
    case Compactness:
    case Eccentricity:
    case LengthMajorAxis:
    case LengthMinorAxis:
    case AxisOrientation:
    case Circularity:
        return {Area, CentralMoment20, CentralMoment02, CentralMoment11};

    default:
        return {};
    }
}

/**
 * @brief Tests whether one descriptor directly consumes another descriptor.
 *
 * @param attribute Attribute requested by the operation.
 * @param dependency Required dependent attribute.
 * @return True if one descriptor directly consumes another descriptor; otherwise false.
 */
inline bool attributeRequiresDependency(Attribute attribute, Attribute dependency) {
    const std::vector<Attribute> dependencies = dependenciesForAttribute(attribute);
    return containsScheduledAttribute(dependencies, dependency);
}

/**
 * @brief Tests whether any descriptor in `attributes` consumes `dependency`.
 *
 * @param attributes Attributes requested by the operation.
 * @param dependency Required dependent attribute.
 * @return True if any descriptor in attributes consumes dependency; otherwise false.
 */
inline bool anyAttributeRequiresDependency(std::span<const Attribute> attributes, Attribute dependency) {
    return std::any_of(attributes.begin(), attributes.end(), [&](Attribute attribute) { return attributeRequiresDependency(attribute, dependency); });
}

/**
 * @brief Adds `attribute` and its recursive dependencies to `closure`.
 *
 * @details
 * Dependencies are appended before their consumer, so the final closure can be
 * traversed in materialization order. `visiting` is a recursion stack used only
 * for cycle detection.
 *
 * @param attribute Attribute requested by the operation.
 * @param closure Dependency closure accumulated by the traversal.
 * @param visiting Attributes on the current dependency traversal path.
 */
inline void appendDependencyClosure(Attribute attribute, std::vector<Attribute>& closure, std::vector<Attribute>& visiting) {
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
 *
 * @param requestedAttributes Attributes requested for materialization.
 * @return The resulting dependency-closed materialization order for a request.
 */
inline std::vector<Attribute> expandDependencyClosure(std::span<const Attribute> requestedAttributes) {
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
     *
     * @param attribute Attribute requested by the operation.
     * @return True if attribute belongs to the public request; otherwise false.
     */
    [[nodiscard]] bool requests(Attribute attribute) const noexcept { return containsScheduledAttribute(requestedAttributes, attribute); }

    /**
     * @brief Tests whether `attribute` will be computed by the plan.
     *
     * @param attribute Attribute requested by the operation.
     * @return True if attribute will be computed by the plan; otherwise false.
     */
    [[nodiscard]] bool materializes(Attribute attribute) const noexcept { return containsScheduledAttribute(materializedAttributes, attribute); }

    /**
     * @brief Tests whether `attribute` is an internal-only dependency.
     *
     * @param attribute Attribute requested by the operation.
     * @return True if attribute is an internal-only dependency; otherwise false.
     */
    [[nodiscard]] bool hides(Attribute attribute) const noexcept { return containsScheduledAttribute(hiddenDependencyAttributes, attribute); }

    /**
     * @brief Returns public requested attributes produced by `family`.
     *
     * @param family Attribute-computer family.
     * @return Public requested attributes produced by family.
     */
    [[nodiscard]] std::vector<Attribute> requestedForFamily(attributes::computers::AttributeComputerFamily family) const {
        return attributesForFamily(requestedAttributes, family);
    }

    /**
     * @brief Returns requested and hidden attributes produced by `family`.
     *
     * @param family Attribute-computer family.
     * @return Requested and hidden attributes produced by family.
     */
    [[nodiscard]] std::vector<Attribute> materializedForFamily(attributes::computers::AttributeComputerFamily family) const {
        return attributesForFamily(materializedAttributes, family);
    }

  private:
    /**
     * @brief Checks whether any requested attribute belongs to a computer family.
     *
     * @param attributes Attributes requested by the operation.
     * @param family Attribute family to schedule.
     * @return True when at least one requested attribute belongs to `family`.
     */
    [[nodiscard]] static std::vector<Attribute> attributesForFamily(const std::vector<Attribute>& attributes,
                                                                    attributes::computers::AttributeComputerFamily family) {
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
 *
 * @param requestedAttributes Attributes requested for materialization.
 * @return The created computation plan from scalar attributes already expanded from groups.
 */
inline AttributeComputationPlan makeAttributeComputationPlan(std::span<const Attribute> requestedAttributes) {
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
