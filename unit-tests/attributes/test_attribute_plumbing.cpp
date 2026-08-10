#include "support/TestSupport.hpp"

#include "mmcfilters/attributes/computers/AttributeComputerRegistry.hpp"
#include "mmcfilters/attributes/computers/BoundingBoxComputer.hpp"
#include "mmcfilters/attributes/computers/GrayLevelStatsComputer.hpp"
#include "mmcfilters/attributes/computers/MomentBasedAttributeComputer.hpp"
#include "mmcfilters/attributes/computers/VolumeComputer.hpp"
#include "mmcfilters/attributes/detail/AttributeFamilyScheduler.hpp"
#include "mmcfilters/attributes/detail/TopologyAttributeBackend.hpp"
#include "mmcfilters/utils/Contract.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <initializer_list>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

using namespace mmcfilters;
using namespace mmcfilters::attributes::computers;
using namespace mmcfilters::unit_tests;

static_assert(AttributeComputer<AreaComputer>);
static_assert(AttributeComputer<BoundingBoxComputer>);
static_assert(AttributeComputer<TreeTopologyComputer>);
static_assert(AttributeComputer<CentralMomentsComputer>);
static_assert(AttributeComputer<HuMomentsComputer>);
static_assert(AttributeComputer<MomentBasedAttributeComputer>);
static_assert(AttributeComputer<BitquadAttributeComputer>);
static_assert(AttributeComputer<ContourSideAttributeComputer>);
static_assert(AttributeComputer<VolumeComputer>);
static_assert(AttributeComputer<GrayLevelStatsComputer>);
static_assert(AttributeComputer<MaxDistComputer>);

static_assert(TopologyAttributeComputer<AreaComputer>);
static_assert(TopologyAttributeComputer<BoundingBoxComputer>);
static_assert(TopologyAttributeComputer<TreeTopologyComputer>);
static_assert(TopologyAttributeComputer<CentralMomentsComputer>);
static_assert(TopologyAttributeComputer<HuMomentsComputer>);
static_assert(TopologyAttributeComputer<MomentBasedAttributeComputer>);
static_assert(TopologyAttributeComputer<BitquadAttributeComputer>);
static_assert(TopologyAttributeComputer<ContourSideAttributeComputer>);

static_assert(AltitudeAttributeComputer<VolumeComputer>);
static_assert(AltitudeAttributeComputer<GrayLevelStatsComputer>);
static_assert(AltitudeAttributeComputer<MaxDistComputer>);
static_assert(std::tuple_size_v<TopologyAttributeComputers> == 8);
static_assert(std::tuple_size_v<AltitudeAttributeComputers> == 3);
static_assert(std::tuple_size_v<RegisteredAttributeComputers> == 11);

template <class Computer>
void requireComputerContract(std::initializer_list<Attribute> producedAttributes, AttributeComputerDomain domain, AttributeComputerFamily family,
                             std::string_view label) {
    const auto& produced = Computer::producedAttributes;
    const std::string labelText(label);

    require(!std::string_view(Computer::familyName).empty(), labelText + " family name");
    requireEqual(produced.size(), producedAttributes.size(), labelText + " produced attribute count");
    requireEqual(static_cast<int>(Computer::domain), static_cast<int>(domain), labelText + " execution domain");
    requireEqual(static_cast<int>(Computer::family), static_cast<int>(family), labelText + " scheduler family");
    requireEqual(numProducedAttributes<Computer>(), producedAttributes.size(), labelText + " produced helper count");

    for (Attribute attribute : producedAttributes) {
        require(producesAttribute<Computer>(attribute), labelText + " produces expected attribute");
    }
    for (Attribute attribute : produced) {
        require(attributes::registry::metadata(attribute) != nullptr, labelText + " produced attribute is registered");
        require(std::find(producedAttributes.begin(), producedAttributes.end(), attribute) != producedAttributes.end(),
                labelText + " has no unexpected produced attribute");
    }
}

template <class Computer> void requireRuntimeProducedAttributesMatchCanonical(std::string_view label) {
    const auto& canonicalAttributes = Computer::producedAttributes;
    const std::vector<Attribute> runtimeAttributes = runtimeProducedAttributes<Computer>();
    const std::string labelText(label);

    requireEqual(runtimeAttributes.size(), canonicalAttributes.size(), labelText + " runtime attribute count");
    for (std::size_t i = 0; i < canonicalAttributes.size(); ++i) {
        requireEqual(static_cast<int>(runtimeAttributes[i]), static_cast<int>(canonicalAttributes[i]),
                     labelText + " runtime vector uses canonical produced attribute order " + std::to_string(i));
    }
}

template <class Computer> void requireRegisteredComputerFamily(AttributeComputerFamily expectedFamily, std::string_view label) {
    const std::string labelText(label);
    for (Attribute attribute : Computer::producedAttributes) {
        requireEqual(static_cast<int>(mmcfilters::detail::familyForAttribute(attribute)), static_cast<int>(expectedFamily),
                     labelText + " registered family for " + AttributeNames::toString(attribute));
        require(mmcfilters::detail::attributeHasComputerDomain<Computer::domain>(attribute),
                labelText + " registered domain for " + AttributeNames::toString(attribute));
    }
}

template <class Computer> void countProducedAttributes(std::array<int, static_cast<std::size_t>(CONTOUR_SIDE_SOUTH) + 1>& counts) {
    for (Attribute attribute : Computer::producedAttributes) {
        const auto index = static_cast<std::size_t>(attribute);
        require(index < counts.size(), "computer produced attribute must be in registry range");
        counts[index] += 1;
    }
}

void requireGlobalAttributeRegistryContracts() {
    std::array<int, static_cast<std::size_t>(CONTOUR_SIDE_SOUTH) + 1> producedCounts{};
    countProducedAttributes<AreaComputer>(producedCounts);
    countProducedAttributes<BoundingBoxComputer>(producedCounts);
    countProducedAttributes<TreeTopologyComputer>(producedCounts);
    countProducedAttributes<CentralMomentsComputer>(producedCounts);
    countProducedAttributes<HuMomentsComputer>(producedCounts);
    countProducedAttributes<MomentBasedAttributeComputer>(producedCounts);
    countProducedAttributes<BitquadAttributeComputer>(producedCounts);
    countProducedAttributes<ContourSideAttributeComputer>(producedCounts);
    countProducedAttributes<VolumeComputer>(producedCounts);
    countProducedAttributes<GrayLevelStatsComputer>(producedCounts);
    countProducedAttributes<MaxDistComputer>(producedCounts);

    const std::vector<Attribute>& allAttributes = ATTRIBUTE_GROUPS.at(AttributeGroup::ALL);
    requireEqual(allAttributes.size(), producedCounts.size(), "ALL group and registered computer trait range must have the same size");

    for (Attribute attribute : allAttributes) {
        const std::string attributeName(AttributeNames::toString(attribute));
        const auto index = static_cast<std::size_t>(attribute);
        const auto* metadata = attributes::registry::metadata(attribute);
        require(metadata != nullptr, attributeName + " must have registry metadata");
        requireEqual(producedCounts[index], 1, attributeName + " must be produced by exactly one registered computer family");
        require(mmcfilters::detail::familyForAttribute(attribute) != AttributeComputerFamily::Unsupported,
                attributeName + " must resolve to a scheduler family");
        require(attributes::registry::isPipelineComputed(attribute), attributeName + " must be computable by the ordinary attribute pipeline");

        const bool topologyFamily =
            mmcfilters::detail::familyForAttribute(attribute) == AttributeComputerFamily::Area || attributes::registry::isTopologyOnly(attribute);
        if (attributes::registry::requiresAltitude(attribute)) {
            require(!topologyFamily, attributeName + " altitude attribute must not be topology-only");
        }
    }
}

void requireGlobalDependencySchedulerContracts() {
    const std::vector<Attribute>& allAttributes = ATTRIBUTE_GROUPS.at(AttributeGroup::ALL);
    const auto allPlan = mmcfilters::detail::makeAttributeComputationPlan(std::span<const Attribute>(allAttributes));

    requireEqual(allPlan.requestedAttributes.size(), allAttributes.size(), "ALL plan requested attribute count");
    requireEqual(allPlan.materializedAttributes.size(), allAttributes.size(), "ALL plan materialized attribute count");
    require(allPlan.hiddenDependencyAttributes.empty(), "ALL plan must not create hidden dependencies because every dependency is requested");

    for (Attribute attribute : allAttributes) {
        const std::string attributeName(AttributeNames::toString(attribute));
        require(allPlan.requests(attribute), attributeName + " must be requested in ALL plan");
        require(allPlan.materializes(attribute), attributeName + " must be materialized in ALL plan");
        require(!allPlan.hides(attribute), attributeName + " must not be hidden in ALL plan");

        const std::array<Attribute, 1> singleRequest{attribute};
        const auto singlePlan = mmcfilters::detail::makeAttributeComputationPlan(std::span<const Attribute>(singleRequest));
        require(singlePlan.requests(attribute), attributeName + " single plan must preserve the request");
        require(singlePlan.materializes(attribute), attributeName + " single plan must materialize the request");
        require(!singlePlan.hides(attribute), attributeName + " single plan must not hide the requested attribute");

        const auto attributePosition = std::find(singlePlan.materializedAttributes.begin(), singlePlan.materializedAttributes.end(), attribute);
        require(attributePosition != singlePlan.materializedAttributes.end(), attributeName + " single plan must contain the requested attribute");

        const std::vector<Attribute> dependencies = mmcfilters::detail::dependenciesForAttribute(attribute);
        for (Attribute dependency : dependencies) {
            const std::string dependencyName(AttributeNames::toString(dependency));
            require(attributes::registry::metadata(dependency) != nullptr, attributeName + " dependency " + dependencyName + " must be registered");
            require(mmcfilters::detail::familyForAttribute(dependency) != AttributeComputerFamily::Unsupported,
                    attributeName + " dependency " + dependencyName + " must resolve to a scheduler family");
            require(singlePlan.materializes(dependency), attributeName + " single plan must materialize dependency " + dependencyName);
            require(singlePlan.hides(dependency), attributeName + " single plan must hide dependency " + dependencyName);

            const auto dependencyPosition = std::find(singlePlan.materializedAttributes.begin(), singlePlan.materializedAttributes.end(), dependency);
            require(dependencyPosition != singlePlan.materializedAttributes.end(),
                    attributeName + " dependency " + dependencyName + " must be present in closure");
            require(dependencyPosition < attributePosition, attributeName + " dependency " + dependencyName + " must be ordered before the consumer");
        }
    }
}

void requireCapabilityRegistryContracts() {
    using mmcfilters::attributes::registry::AttributeAdjacencyRequirement;
    using mmcfilters::attributes::registry::capabilityRequirements;

    const auto area = capabilityRequirements(AREA);
    require(!area.altitude && !area.gridDomain2D && area.adjacency == AttributeAdjacencyRequirement::NONE && !area.monotoneAltitudeOrder,
            "AREA capability contract");

    const auto box = capabilityRequirements(BOX_WIDTH);
    require(!box.altitude && box.gridDomain2D && box.adjacency == AttributeAdjacencyRequirement::NONE, "BOX_WIDTH capability contract");

    const auto bitquad = capabilityRequirements(BITQUADS_AREA);
    require(!bitquad.altitude && bitquad.gridDomain2D && bitquad.adjacency == AttributeAdjacencyRequirement::UNIFORM_OR_DIRECTIONAL &&
                bitquad.altitudeForDirectionalAdjacency && bitquad.canonical4Or8Adjacency,
            "BITQUADS_AREA capability contract");

    const auto maxDist = capabilityRequirements(MAX_DIST);
    require(maxDist.altitude && maxDist.gridDomain2D && maxDist.adjacency == AttributeAdjacencyRequirement::UNIFORM && maxDist.monotoneAltitudeOrder,
            "MAX_DIST capability contract");

    for (Attribute attribute : ATTRIBUTE_GROUPS.at(AttributeGroup::ALL)) {
        const auto requirements = capabilityRequirements(attribute);
        const auto* item = attributes::registry::metadata(attribute);
        require(item != nullptr && item->requirements == requirements, std::string(AttributeNames::toString(attribute)) + " capability metadata agreement");
        require(requirements.altitude == attributes::registry::requiresAltitude(attribute),
                std::string(AttributeNames::toString(attribute)) + " altitude metadata agreement");
    }

    const auto unknown = capabilityRequirements(static_cast<Attribute>(-1));
    require(unknown == attributes::registry::AttributeCapabilityRequirements{}, "unknown attribute must have an empty capability contract");
}

int main() {
    {
        auto names = AttributeNames::fromList({AREA, VOLUME, MAX_DIST});
        requireEqual(names.NUM_ATTRIBUTES, 3, "AttributeNames count");
        requireEqual(names.getIndex(AREA), 0, "AREA index in AttributeNames");
        requireEqual(names.offset(VOLUME), 1, "VOLUME offset in AttributeNames");
        require(names.contains(AREA), "AttributeNames::contains must find requested attribute");
        require(!names.contains(LEVEL), "AttributeNames::contains must reject absent attribute");
        requireEqual(names.getIndex(VOLUME), 1, "VOLUME index in AttributeNames");
        requireEqual(names.getIndex(MAX_DIST), 2, "MAX_DIST index in AttributeNames");
        requireEqual(names.linearIndex(2, MAX_DIST), 8, "MAX_DIST linear index");
        for (int nodeIndex = 0; nodeIndex < 6; ++nodeIndex) {
            for (Attribute attribute : {AREA, VOLUME, MAX_DIST}) {
                require(names.linearIndex(nodeIndex, attribute) < 6 * names.NUM_ATTRIBUTES, "AttributeNames dense linear index must stay within buffer bounds");
            }
        }
        if constexpr (contract::validationsEnabled) {
            requireThrows<std::invalid_argument>([]() { static_cast<void>(AttributeNames::fromList({AREA, AREA, VOLUME})); },
                                                 "AttributeNames::fromList must reject duplicate attributes");
        }
        auto grayLevelNames = AttributeNames::fromGroup(AttributeGroup::GRAY_LEVEL);
        requireEqual(grayLevelNames.NUM_ATTRIBUTES, 6, "GRAY_LEVEL AttributeNames count");
        requireEqual(grayLevelNames.getIndex(VOLUME), 0, "GRAY_LEVEL VOLUME index");
        requireEqual(grayLevelNames.getIndex(VARIANCE_LEVEL), 5, "GRAY_LEVEL VARIANCE_LEVEL index");
        requireEqual(grayLevelNames.linearIndex(5, VARIANCE_LEVEL), 35, "GRAY_LEVEL dense linear index");
        auto shapeNames = AttributeNames::fromGroup(AttributeGroup::SHAPE);
        requireEqual(shapeNames.NUM_ATTRIBUTES, 47, "SHAPE AttributeNames count");
        requireEqual(shapeNames.getIndex(AREA), 0, "SHAPE AREA index");
        require(shapeNames.contains(MAX_DIST), "SHAPE AttributeNames must include MAX_DIST");
        require(shapeNames.contains(CIRCULARITY), "SHAPE AttributeNames must include CIRCULARITY");
        requireEqual(shapeNames.getIndex(CONTOUR_SIDE_SOUTH), 46, "SHAPE CONTOUR_SIDE_SOUTH index");
        auto momentNames = AttributeNames::fromGroup(AttributeGroup::MOMENTS);
        requireEqual(momentNames.NUM_ATTRIBUTES, 21, "MOMENTS AttributeNames count");
        requireEqual(momentNames.getIndex(CENTRAL_MOMENT_20), 0, "MOMENTS CENTRAL_MOMENT_20 index");
        requireEqual(momentNames.getIndex(CIRCULARITY), 20, "MOMENTS CIRCULARITY index");
        auto boundaryNames = AttributeNames::fromGroup(AttributeGroup::BOUNDARY);
        requireEqual(boundaryNames.NUM_ATTRIBUTES, 15, "BOUNDARY AttributeNames count");
        requireEqual(boundaryNames.getIndex(BITQUADS_AREA), 0, "BOUNDARY BITQUADS_AREA index");
        requireEqual(boundaryNames.getIndex(CONTOUR_SIDE_SOUTH), 14, "BOUNDARY CONTOUR_SIDE_SOUTH index");
        auto allNames = AttributeNames::fromGroup(AttributeGroup::ALL);
        requireEqual(allNames.NUM_ATTRIBUTES, static_cast<int>(CONTOUR_SIDE_SOUTH) + 1, "ALL AttributeNames count");
        requireEqual(allNames.getIndex(CONTOUR_PIXELS), static_cast<int>(CONTOUR_PIXELS), "ALL CONTOUR_PIXELS index");
        requireEqual(allNames.getIndex(CONTOUR_SIDE_SOUTH), static_cast<int>(CONTOUR_SIDE_SOUTH), "ALL CONTOUR_SIDE_SOUTH index");
        const std::vector<Attribute>& allAttributes = ATTRIBUTE_GROUPS.at(AttributeGroup::ALL);
        for (int index = 0; index < static_cast<int>(allAttributes.size()); ++index) {
            const Attribute attribute = allAttributes[static_cast<std::size_t>(index)];
            requireEqual(static_cast<int>(attribute), index, "ALL group preserves enum order");
            const auto parsed = AttributeNames::parse(AttributeNames::toString(attribute));
            require(parsed.has_value(), "ALL group attribute string must parse");
            requireEqual(static_cast<int>(*parsed), index, "ALL group parse round-trip");
        }
        requireEqual(AttributeNames::toString(MAX_DIST), std::string("MAX_DIST"), "AttributeNames::toString");
        require(AttributeNames::describe(AREA).starts_with("Area:"), "AttributeNames::describe AREA");

        auto parsedArea = AttributeNames::parse("AREA");
        require(parsedArea.has_value(), "AttributeNames::parse AREA must succeed");
        requireEqual(static_cast<int>(*parsedArea), static_cast<int>(AREA), "AttributeNames::parse AREA value");
        auto parsedContour = AttributeNames::parse("CONTOUR_PIXELS");
        require(parsedContour.has_value(), "AttributeNames::parse CONTOUR_PIXELS must succeed");
        requireEqual(static_cast<int>(*parsedContour), static_cast<int>(CONTOUR_PIXELS), "AttributeNames::parse CONTOUR_PIXELS value");
        require(!AttributeNames::parse("NOT_AN_ATTRIBUTE").has_value(), "AttributeNames::parse invalid must fail");

        auto deltaNames = AttributeNamesWithDelta::create(2, {AREA, LEVEL});
        requireEqual(deltaNames.NUM_ATTRIBUTES, 10, "AttributeNamesWithDelta count");
        requireEqual(deltaNames.getIndex(AREA, -2), 0, "AREA asc2 offset");
        requireEqual(deltaNames.getIndex(LEVEL, 0), 5, "LEVEL center offset");
        requireEqual(deltaNames.linearIndex(1, AREA, -2), 10, "AREA asc2 linear index");
        requireEqual(deltaNames.linearIndex(1, LEVEL, 2), 19, "LEVEL desc2 linear index");
        requireEqual(AttributeNamesWithDelta::toString(AREA, -2), std::string("AREA_ASC_2"), "AttributeNamesWithDelta asc label");
        requireEqual(AttributeNamesWithDelta::toString(LEVEL, 2), std::string("LEVEL_DESC_2"), "AttributeNamesWithDelta desc label");
        if constexpr (contract::validationsEnabled) {
            requireThrows<std::invalid_argument>([]() { static_cast<void>(AttributeNamesWithDelta::create(1, {AREA, AREA})); },
                                                 "AttributeNamesWithDelta::create must reject duplicate attributes");
        }
    }

    {
        requireComputerContract<AreaComputer>({AREA}, AttributeComputerDomain::Topology, AttributeComputerFamily::Area, "AreaComputer");
        requireComputerContract<BoundingBoxComputer>(
            {BOX_WIDTH, BOX_HEIGHT, DIAGONAL_LENGTH, RECTANGULARITY, RATIO_WH, BOX_COL_MIN, BOX_COL_MAX, BOX_ROW_MIN, BOX_ROW_MAX},
            AttributeComputerDomain::Topology, AttributeComputerFamily::BoundingBox, "BoundingBoxComputer");
        requireComputerContract<TreeTopologyComputer>({HEIGHT_NODE, DEPTH_NODE, IS_LEAF_NODE, IS_ROOT_NODE, NUM_CHILDREN_NODE, NUM_SIBLINGS_NODE,
                                                       NUM_DESCENDANTS_NODE, NUM_LEAF_DESCENDANTS_NODE, LEAF_RATIO_NODE, BALANCE_NODE, AVG_CHILD_HEIGHT_NODE},
                                                      AttributeComputerDomain::Topology, AttributeComputerFamily::TreeTopology, "TreeTopologyComputer");
        requireComputerContract<CentralMomentsComputer>(
            {CENTRAL_MOMENT_20, CENTRAL_MOMENT_02, CENTRAL_MOMENT_11, CENTRAL_MOMENT_30, CENTRAL_MOMENT_03, CENTRAL_MOMENT_21, CENTRAL_MOMENT_12},
            AttributeComputerDomain::Topology, AttributeComputerFamily::CentralMoments, "CentralMomentsComputer");
        requireComputerContract<HuMomentsComputer>({HU_MOMENT_1, HU_MOMENT_2, HU_MOMENT_3, HU_MOMENT_4, HU_MOMENT_5, HU_MOMENT_6, HU_MOMENT_7},
                                                   AttributeComputerDomain::Topology, AttributeComputerFamily::HuMoments, "HuMomentsComputer");
        requireComputerContract<MomentBasedAttributeComputer>(
            {INERTIA, COMPACTNESS, ECCENTRICITY, LENGTH_MAJOR_AXIS, LENGTH_MINOR_AXIS, AXIS_ORIENTATION, CIRCULARITY}, AttributeComputerDomain::Topology,
            AttributeComputerFamily::MomentDerived, "MomentBasedAttributeComputer");
        requireComputerContract<BitquadAttributeComputer>({BITQUADS_AREA, BITQUADS_NUMBER_EULER, BITQUADS_NUMBER_HOLES, BITQUADS_PERIMETER,
                                                           BITQUADS_PERIMETER_CONTINUOUS, BITQUADS_CIRCULARITY, BITQUADS_PERIMETER_AVERAGE,
                                                           BITQUADS_LENGTH_AVERAGE, BITQUADS_WIDTH_AVERAGE},
                                                          AttributeComputerDomain::Topology, AttributeComputerFamily::Bitquad, "BitquadAttributeComputer");
        requireComputerContract<ContourSideAttributeComputer>(
            {CONTOUR_PIXELS, CONTOUR_PERIMETER, CONTOUR_SIDE_NORTH, CONTOUR_SIDE_WEST, CONTOUR_SIDE_EAST, CONTOUR_SIDE_SOUTH},
            AttributeComputerDomain::Topology, AttributeComputerFamily::ContourSide, "ContourSideAttributeComputer");
        requireComputerContract<VolumeComputer>({VOLUME, RELATIVE_VOLUME}, AttributeComputerDomain::Altitude, AttributeComputerFamily::Volume,
                                                "VolumeComputer");
        requireComputerContract<GrayLevelStatsComputer>({LEVEL, MEAN_LEVEL, VARIANCE_LEVEL, GRAY_HEIGHT}, AttributeComputerDomain::Altitude,
                                                        AttributeComputerFamily::GrayLevelStats, "GrayLevelStatsComputer");
        requireComputerContract<MaxDistComputer>({MAX_DIST}, AttributeComputerDomain::Altitude, AttributeComputerFamily::MaxDist, "MaxDistComputer");

        requireRuntimeProducedAttributesMatchCanonical<AreaComputer>("AreaComputer");
        requireRuntimeProducedAttributesMatchCanonical<BoundingBoxComputer>("BoundingBoxComputer");
        requireRuntimeProducedAttributesMatchCanonical<TreeTopologyComputer>("TreeTopologyComputer");
        requireRuntimeProducedAttributesMatchCanonical<CentralMomentsComputer>("CentralMomentsComputer");
        requireRuntimeProducedAttributesMatchCanonical<HuMomentsComputer>("HuMomentsComputer");
        requireRuntimeProducedAttributesMatchCanonical<MomentBasedAttributeComputer>("MomentBasedAttributeComputer");
        requireRuntimeProducedAttributesMatchCanonical<BitquadAttributeComputer>("BitquadAttributeComputer");
        requireRuntimeProducedAttributesMatchCanonical<ContourSideAttributeComputer>("ContourSideAttributeComputer");
        requireRuntimeProducedAttributesMatchCanonical<VolumeComputer>("VolumeComputer");
        requireRuntimeProducedAttributesMatchCanonical<GrayLevelStatsComputer>("GrayLevelStatsComputer");
        requireRuntimeProducedAttributesMatchCanonical<MaxDistComputer>("MaxDistComputer");

        requireRegisteredComputerFamily<AreaComputer>(AttributeComputerFamily::Area, "AreaComputer");
        requireRegisteredComputerFamily<BoundingBoxComputer>(AttributeComputerFamily::BoundingBox, "BoundingBoxComputer");
        requireRegisteredComputerFamily<TreeTopologyComputer>(AttributeComputerFamily::TreeTopology, "TreeTopologyComputer");
        requireRegisteredComputerFamily<CentralMomentsComputer>(AttributeComputerFamily::CentralMoments, "CentralMomentsComputer");
        requireRegisteredComputerFamily<HuMomentsComputer>(AttributeComputerFamily::HuMoments, "HuMomentsComputer");
        requireRegisteredComputerFamily<MomentBasedAttributeComputer>(AttributeComputerFamily::MomentDerived, "MomentBasedAttributeComputer");
        requireRegisteredComputerFamily<BitquadAttributeComputer>(AttributeComputerFamily::Bitquad, "BitquadAttributeComputer");
        requireRegisteredComputerFamily<ContourSideAttributeComputer>(AttributeComputerFamily::ContourSide, "ContourSideAttributeComputer");
        requireRegisteredComputerFamily<VolumeComputer>(AttributeComputerFamily::Volume, "VolumeComputer");
        requireRegisteredComputerFamily<GrayLevelStatsComputer>(AttributeComputerFamily::GrayLevelStats, "GrayLevelStatsComputer");
        requireRegisteredComputerFamily<MaxDistComputer>(AttributeComputerFamily::MaxDist, "MaxDistComputer");

        requireEqual(static_cast<int>(mmcfilters::detail::familyForAttribute(RECTANGULARITY)), static_cast<int>(AttributeComputerFamily::BoundingBox),
                     "scheduler family lookup for bounding boxes");
        requireEqual(static_cast<int>(mmcfilters::detail::familyForAttribute(AREA)), static_cast<int>(AttributeComputerFamily::Area),
                     "scheduler family lookup for AREA");
        requireEqual(static_cast<int>(mmcfilters::detail::familyForAttribute(AVG_CHILD_HEIGHT_NODE)), static_cast<int>(AttributeComputerFamily::TreeTopology),
                     "scheduler family lookup for tree topology");
        requireEqual(static_cast<int>(mmcfilters::detail::familyForAttribute(CENTRAL_MOMENT_20)), static_cast<int>(AttributeComputerFamily::CentralMoments),
                     "scheduler family lookup for central moments");
        requireEqual(static_cast<int>(mmcfilters::detail::familyForAttribute(HU_MOMENT_7)), static_cast<int>(AttributeComputerFamily::HuMoments),
                     "scheduler family lookup for Hu moments");
        requireEqual(static_cast<int>(mmcfilters::detail::familyForAttribute(ECCENTRICITY)), static_cast<int>(AttributeComputerFamily::MomentDerived),
                     "scheduler family lookup for moment-derived attributes");
        requireEqual(static_cast<int>(mmcfilters::detail::familyForAttribute(BITQUADS_CIRCULARITY)), static_cast<int>(AttributeComputerFamily::Bitquad),
                     "scheduler family lookup for bitquads");
        requireEqual(static_cast<int>(mmcfilters::detail::familyForAttribute(CONTOUR_SIDE_SOUTH)), static_cast<int>(AttributeComputerFamily::ContourSide),
                     "scheduler family lookup for contour sides");

        requireGlobalAttributeRegistryContracts();
    }

    {
        const std::array<Attribute, 3> scheduledRequest{MEAN_LEVEL, ECCENTRICITY, RECTANGULARITY};
        const mmcfilters::detail::AttributeComputationPlan schedulerPlan =
            mmcfilters::detail::makeAttributeComputationPlan(std::span<const Attribute>(scheduledRequest));

        require(schedulerPlan.requests(MEAN_LEVEL), "scheduler preserves requested MEAN_LEVEL");
        require(schedulerPlan.requests(ECCENTRICITY), "scheduler preserves requested ECCENTRICITY");
        require(schedulerPlan.requests(RECTANGULARITY), "scheduler preserves requested RECTANGULARITY");
        require(schedulerPlan.materializes(VOLUME), "scheduler adds hidden VOLUME dependency");
        require(schedulerPlan.materializes(AREA), "scheduler adds hidden AREA dependency");
        require(schedulerPlan.materializes(CENTRAL_MOMENT_20), "scheduler adds hidden central moment dependency");
        require(schedulerPlan.materializes(CENTRAL_MOMENT_02), "scheduler adds hidden central moment dependency");
        require(schedulerPlan.materializes(CENTRAL_MOMENT_11), "scheduler adds hidden central moment dependency");
        require(schedulerPlan.hides(VOLUME), "scheduler marks VOLUME as hidden");
        require(schedulerPlan.hides(AREA), "scheduler marks AREA as hidden");
        require(!schedulerPlan.hides(MEAN_LEVEL), "scheduler does not hide requested attributes");

        const std::vector<Attribute> scheduledVolume = schedulerPlan.materializedForFamily(AttributeComputerFamily::Volume);
        requireEqual(scheduledVolume.size(), static_cast<std::size_t>(1), "scheduler volume family count");
        requireEqual(static_cast<int>(scheduledVolume.front()), static_cast<int>(VOLUME), "scheduler volume family attribute");

        const std::vector<Attribute> scheduledGray = schedulerPlan.requestedForFamily(AttributeComputerFamily::GrayLevelStats);
        requireEqual(scheduledGray.size(), static_cast<std::size_t>(1), "scheduler gray family count");
        requireEqual(static_cast<int>(scheduledGray.front()), static_cast<int>(MEAN_LEVEL), "scheduler gray family attribute");

        const std::array<Attribute, 1> rectangularityOnly{RECTANGULARITY};
        const std::array<Attribute, 1> widthOnly{BOX_WIDTH};
        require(mmcfilters::detail::anyAttributeRequiresDependency(std::span<const Attribute>(rectangularityOnly), AREA),
                "scheduler identifies RECTANGULARITY AREA dependency");
        require(!mmcfilters::detail::anyAttributeRequiresDependency(std::span<const Attribute>(widthOnly), AREA),
                "scheduler does not add AREA to independent bounding-box attributes");
        requireEqual(static_cast<int>(mmcfilters::detail::familyForAttribute(ECCENTRICITY)), static_cast<int>(AttributeComputerFamily::MomentDerived),
                     "scheduler family lookup for ECCENTRICITY");
        requireEqual(static_cast<int>(mmcfilters::detail::familyForAttribute(CONTOUR_SIDE_SOUTH)), static_cast<int>(AttributeComputerFamily::ContourSide),
                     "scheduler family lookup for contour sides");

        requireGlobalDependencySchedulerContracts();
        requireCapabilityRegistryContracts();
    }

    {
        const AttributeNames areaNames = AttributeNames::fromList({AREA});
        const AttributeNames volumeNames = AttributeNames::fromList({VOLUME, RELATIVE_VOLUME});
        const std::array<float, 2> areaValues{4.0f, 8.0f};
        const std::array<float, 4> volumeValues{12.0f, 1.5f, 30.0f, 3.75f};
        const std::array<DependencySourceT<float>, 2> sources{{
            DependencySourceT<float>{&volumeNames, volumeValues.data()},
            DependencySourceT<float>{&areaNames, areaValues.data()},
        }};
        const DependencyResolver<float> resolver{std::span<const DependencySourceT<float>>(sources)};

        require(&resolver.require(AREA) == &sources[1], "DependencyResolver resolves AREA by name");
        require(&resolver.require(VOLUME) == &sources[0], "DependencyResolver resolves VOLUME by name");
        require(&resolver.requireAll({VOLUME, RELATIVE_VOLUME}) == &sources[0], "DependencyResolver resolves multi-attribute source");
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(resolver.require(MEAN_LEVEL)); }, "DependencyResolver rejects missing named dependency");

        const std::array<DependencySourceT<float>, 1> invalidSources{{DependencySourceT<float>{nullptr, areaValues.data()}}};
        const DependencyResolver<float> invalidResolver{std::span<const DependencySourceT<float>>(invalidSources)};
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(invalidResolver.require(AREA)); },
                                             "DependencyResolver rejects invalid dependency source");

        requireNear(::mmcfilters::attributes::numeric::safeDivide(6.0f, 3.0f), 2.0f, 0.0f, "safeDivide normal division");
        requireEqual(::mmcfilters::attributes::numeric::safeDivide(6.0f, 0.0f, -1.0f), -1.0f, "safeDivide returns finite fallback");
        requireEqual(::mmcfilters::attributes::numeric::safeSqrt(-4.0f), 0.0f, "safeSqrt clamps negative input");
        requireEqual(::mmcfilters::attributes::numeric::clampUpper(5.0f, 3.0f), 3.0f, "clampUpper caps values");
    }

    {
        auto tree = makeComponentTree(makeComponentTreeFixture(), true);

        const std::array<Attribute, 1> widthRequest{BOX_WIDTH};
        const AttributeNames widthNames = AttributeNames::fromList({BOX_WIDTH});
        std::vector<float> widthBuffer(static_cast<std::size_t>(tree->getNumInternalNodeSlots()) * static_cast<std::size_t>(widthNames.NUM_ATTRIBUTES), 0.0f);
        BoundingBoxComputer::compute(
            AttributeComputeContext<float>{*tree, std::span<float>(widthBuffer), widthNames, std::span<const Attribute>(widthRequest)});
        requireEqual(widthBuffer[widthNames.linearIndex(tree->getRoot(), BOX_WIDTH)], 4.0f, "BOX_WIDTH context compute does not require AREA dependency");

        if constexpr (contract::validationsEnabled) {
            const std::array<Attribute, 1> rectangularityRequest{RECTANGULARITY};
            const AttributeNames rectangularityNames = AttributeNames::fromList({RECTANGULARITY});
            std::vector<float> rectangularityBuffer(
                static_cast<std::size_t>(tree->getNumInternalNodeSlots()) * static_cast<std::size_t>(rectangularityNames.NUM_ATTRIBUTES), 0.0f);
            requireThrows<std::invalid_argument>(
                [&]() {
                    BoundingBoxComputer::compute(AttributeComputeContext<float>{*tree, std::span<float>(rectangularityBuffer), rectangularityNames,
                                                                                std::span<const Attribute>(rectangularityRequest)});
                },
                "RECTANGULARITY context compute must reject missing AREA dependency");
        }
    }

    {
        auto weighted = makeWeightedComponentTree(makeComponentTreeFixture(), true);
        const MorphologicalTree& tree = weighted->topology();
        const std::span<const std::uint8_t> altitude = weighted->altitudeSpan();

        const std::array<Attribute, 1> volumeRequest{VOLUME};
        const AttributeNames volumeNames = AttributeNames::fromList({VOLUME});
        std::vector<float> volumeBuffer(static_cast<std::size_t>(tree.getNumInternalNodeSlots()) * static_cast<std::size_t>(volumeNames.NUM_ATTRIBUTES), 0.0f);
        VolumeComputer::compute(AltitudeAttributeComputeContext<float, std::uint8_t>{tree, altitude, std::span<float>(volumeBuffer), volumeNames,
                                                                                     std::span<const Attribute>(volumeRequest)});
        requireEqual(volumeBuffer[volumeNames.linearIndex(tree.getRoot(), VOLUME)], 42.0f, "VOLUME context compute does not require AREA dependency");

        if constexpr (contract::validationsEnabled) {
            const std::array<Attribute, 1> relativeVolumeRequest{RELATIVE_VOLUME};
            const AttributeNames relativeVolumeNames = AttributeNames::fromList({RELATIVE_VOLUME});
            std::vector<float> relativeVolumeBuffer(
                static_cast<std::size_t>(tree.getNumInternalNodeSlots()) * static_cast<std::size_t>(relativeVolumeNames.NUM_ATTRIBUTES), 0.0f);
            requireThrows<std::invalid_argument>(
                [&]() {
                    VolumeComputer::compute(AltitudeAttributeComputeContext<float, std::uint8_t>{
                        tree, altitude, std::span<float>(relativeVolumeBuffer), relativeVolumeNames, std::span<const Attribute>(relativeVolumeRequest)});
                },
                "RELATIVE_VOLUME context compute must reject missing AREA dependency");
        }

        const std::array<Attribute, 1> levelRequest{LEVEL};
        const AttributeNames levelNames = AttributeNames::fromList({LEVEL});
        std::vector<float> levelBuffer(static_cast<std::size_t>(tree.getNumInternalNodeSlots()) * static_cast<std::size_t>(levelNames.NUM_ATTRIBUTES), 0.0f);
        GrayLevelStatsComputer::compute(AltitudeAttributeComputeContext<float, std::uint8_t>{tree, altitude, std::span<float>(levelBuffer), levelNames,
                                                                                             std::span<const Attribute>(levelRequest)});
        requireEqual(levelBuffer[levelNames.linearIndex(tree.getRoot(), LEVEL)], static_cast<float>(weighted->getAltitude(tree.getRoot())),
                     "LEVEL context compute does not require aggregate dependencies");

        const std::array<Attribute, 1> grayHeightRequest{GRAY_HEIGHT};
        const AttributeNames grayHeightNames = AttributeNames::fromList({GRAY_HEIGHT});
        std::vector<float> grayHeightBuffer(static_cast<std::size_t>(tree.getNumInternalNodeSlots()) * static_cast<std::size_t>(grayHeightNames.NUM_ATTRIBUTES),
                                            0.0f);
        GrayLevelStatsComputer::compute(AltitudeAttributeComputeContext<float, std::uint8_t>{tree, altitude, std::span<float>(grayHeightBuffer),
                                                                                             grayHeightNames, std::span<const Attribute>(grayHeightRequest)});
        require(grayHeightBuffer[grayHeightNames.linearIndex(tree.getRoot(), GRAY_HEIGHT)] >= 0.0f,
                "GRAY_HEIGHT context compute does not require aggregate dependencies");

        if constexpr (contract::validationsEnabled) {
            const std::array<Attribute, 1> meanRequest{MEAN_LEVEL};
            const AttributeNames meanNames = AttributeNames::fromList({MEAN_LEVEL});
            std::vector<float> meanBuffer(static_cast<std::size_t>(tree.getNumInternalNodeSlots()) * static_cast<std::size_t>(meanNames.NUM_ATTRIBUTES), 0.0f);
            requireThrows<std::invalid_argument>(
                [&]() {
                    GrayLevelStatsComputer::compute(AltitudeAttributeComputeContext<float, std::uint8_t>{tree, altitude, std::span<float>(meanBuffer), meanNames,
                                                                                                         std::span<const Attribute>(meanRequest)});
                },
                "MEAN_LEVEL context compute must reject missing VOLUME and AREA dependencies");

            const std::array<Attribute, 1> huRequest{HU_MOMENT_1};
            const AttributeNames huNames = AttributeNames::fromList({HU_MOMENT_1});
            std::vector<float> huBuffer(static_cast<std::size_t>(tree.getNumInternalNodeSlots()) * static_cast<std::size_t>(huNames.NUM_ATTRIBUTES), 0.0f);
            requireThrows<std::invalid_argument>(
                [&]() {
                    HuMomentsComputer::compute(AttributeComputeContext<float>{tree, std::span<float>(huBuffer), huNames, std::span<const Attribute>(huRequest)});
                },
                "HU_MOMENT context compute must reject missing central moment and AREA dependencies");

            const std::array<Attribute, 1> inertiaRequest{INERTIA};
            const AttributeNames inertiaNames = AttributeNames::fromList({INERTIA});
            std::vector<float> inertiaBuffer(
                static_cast<std::size_t>(tree.getNumInternalNodeSlots()) * static_cast<std::size_t>(inertiaNames.NUM_ATTRIBUTES), 0.0f);
            requireThrows<std::invalid_argument>(
                [&]() {
                    MomentBasedAttributeComputer::compute(
                        AttributeComputeContext<float>{tree, std::span<float>(inertiaBuffer), inertiaNames, std::span<const Attribute>(inertiaRequest)});
                },
                "moment-based context compute must reject missing central moment and AREA dependencies");
        }
    }

    return 0;
}
