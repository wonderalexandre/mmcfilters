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

template <class Computer> void countProducedAttributes(std::array<int, static_cast<std::size_t>(ContourSideSouth) + 1>& counts) {
    for (Attribute attribute : Computer::producedAttributes) {
        const auto index = static_cast<std::size_t>(attribute);
        require(index < counts.size(), "computer produced attribute must be in registry range");
        counts[index] += 1;
    }
}

void requireGlobalAttributeRegistryContracts() {
    std::array<int, static_cast<std::size_t>(ContourSideSouth) + 1> producedCounts{};
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

    const std::vector<Attribute>& allAttributes = ATTRIBUTE_GROUPS.at(AttributeGroup::All);
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
    const std::vector<Attribute>& allAttributes = ATTRIBUTE_GROUPS.at(AttributeGroup::All);
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

    const auto area = capabilityRequirements(Area);
    require(!area.altitude && !area.gridDomain2D && area.adjacency == AttributeAdjacencyRequirement::None && !area.monotoneAltitudeOrder,
            "AREA capability contract");

    const auto box = capabilityRequirements(BoxWidth);
    require(!box.altitude && box.gridDomain2D && box.adjacency == AttributeAdjacencyRequirement::None, "BOX_WIDTH capability contract");

    const auto bitquad = capabilityRequirements(BitquadArea);
    require(!bitquad.altitude && bitquad.gridDomain2D && bitquad.adjacency == AttributeAdjacencyRequirement::UniformOrDirectional &&
                bitquad.altitudeForDirectionalAdjacency && bitquad.canonical4Or8Adjacency,
            "BITQUAD_AREA capability contract");

    const auto maxDist = capabilityRequirements(MaxDist);
    require(maxDist.altitude && maxDist.gridDomain2D && maxDist.adjacency == AttributeAdjacencyRequirement::Uniform && maxDist.monotoneAltitudeOrder,
            "MAX_DIST capability contract");

    for (Attribute attribute : ATTRIBUTE_GROUPS.at(AttributeGroup::All)) {
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
        auto names = AttributeNames::fromList({Area, Volume, MaxDist});
        requireEqual(names.NUM_ATTRIBUTES, 3, "AttributeNames count");
        requireEqual(names.getIndex(Area), 0, "AREA index in AttributeNames");
        requireEqual(names.offset(Volume), 1, "VOLUME offset in AttributeNames");
        require(names.contains(Area), "AttributeNames::contains must find requested attribute");
        require(!names.contains(GrayLevelHeight), "AttributeNames::contains must reject absent attribute");
        requireEqual(names.getIndex(Volume), 1, "VOLUME index in AttributeNames");
        requireEqual(names.getIndex(MaxDist), 2, "MAX_DIST index in AttributeNames");
        requireEqual(names.linearIndex(2, MaxDist), 8, "MAX_DIST linear index");
        for (int nodeIndex = 0; nodeIndex < 6; ++nodeIndex) {
            for (Attribute attribute : {Area, Volume, MaxDist}) {
                require(names.linearIndex(nodeIndex, attribute) < 6 * names.NUM_ATTRIBUTES, "AttributeNames dense linear index must stay within buffer bounds");
            }
        }
        if constexpr (contract::validationsEnabled) {
            requireThrows<std::invalid_argument>([]() { static_cast<void>(AttributeNames::fromList({Area, Area, Volume})); },
                                                 "AttributeNames::fromList must reject duplicate attributes");
        }
        auto grayLevelNames = AttributeNames::fromGroup(AttributeGroup::GrayLevel);
        requireEqual(grayLevelNames.NUM_ATTRIBUTES, 5, "GRAY_LEVEL AttributeNames count");
        requireEqual(grayLevelNames.getIndex(Volume), 0, "GRAY_LEVEL VOLUME index");
        requireEqual(grayLevelNames.getIndex(GrayLevelVariance), 4, "GRAY_LEVEL GrayLevelVariance index");
        requireEqual(grayLevelNames.linearIndex(5, GrayLevelVariance), 29, "GRAY_LEVEL dense linear index");
        auto shapeNames = AttributeNames::fromGroup(AttributeGroup::Shape);
        requireEqual(shapeNames.NUM_ATTRIBUTES, 47, "SHAPE AttributeNames count");
        requireEqual(shapeNames.getIndex(Area), 0, "SHAPE AREA index");
        require(shapeNames.contains(MaxDist), "SHAPE AttributeNames must include MAX_DIST");
        require(shapeNames.contains(Circularity), "SHAPE AttributeNames must include CIRCULARITY");
        requireEqual(shapeNames.getIndex(ContourSideSouth), 46, "SHAPE CONTOUR_SIDE_SOUTH index");
        auto momentNames = AttributeNames::fromGroup(AttributeGroup::Moments);
        requireEqual(momentNames.NUM_ATTRIBUTES, 21, "MOMENTS AttributeNames count");
        requireEqual(momentNames.getIndex(CentralMoment20), 0, "MOMENTS CENTRAL_MOMENT_20 index");
        requireEqual(momentNames.getIndex(Circularity), 20, "MOMENTS CIRCULARITY index");
        auto boundaryNames = AttributeNames::fromGroup(AttributeGroup::Boundary);
        requireEqual(boundaryNames.NUM_ATTRIBUTES, 15, "BOUNDARY AttributeNames count");
        requireEqual(boundaryNames.getIndex(BitquadArea), 0, "BOUNDARY BITQUAD_AREA index");
        requireEqual(boundaryNames.getIndex(ContourSideSouth), 14, "BOUNDARY CONTOUR_SIDE_SOUTH index");
        auto allNames = AttributeNames::fromGroup(AttributeGroup::All);
        requireEqual(allNames.NUM_ATTRIBUTES, static_cast<int>(ContourSideSouth) + 1, "ALL AttributeNames count");
        requireEqual(allNames.getIndex(ContourPixels), static_cast<int>(ContourPixels), "ALL CONTOUR_PIXELS index");
        requireEqual(allNames.getIndex(ContourSideSouth), static_cast<int>(ContourSideSouth), "ALL CONTOUR_SIDE_SOUTH index");
        const std::vector<Attribute>& allAttributes = ATTRIBUTE_GROUPS.at(AttributeGroup::All);
        for (int index = 0; index < static_cast<int>(allAttributes.size()); ++index) {
            const Attribute attribute = allAttributes[static_cast<std::size_t>(index)];
            requireEqual(static_cast<int>(attribute), index, "ALL group preserves enum order");
            const auto parsed = AttributeNames::parse(AttributeNames::toString(attribute));
            require(parsed.has_value(), "ALL group attribute string must parse");
            requireEqual(static_cast<int>(*parsed), index, "ALL group parse round-trip");
        }
        requireEqual(AttributeNames::toString(MaxDist), std::string("MAX_DIST"), "AttributeNames::toString");
        require(AttributeNames::describe(Area).starts_with("Area:"), "AttributeNames::describe AREA");

        auto parsedArea = AttributeNames::parse("AREA");
        require(parsedArea.has_value(), "AttributeNames::parse AREA must succeed");
        requireEqual(static_cast<int>(*parsedArea), static_cast<int>(Area), "AttributeNames::parse AREA value");
        auto parsedContour = AttributeNames::parse("CONTOUR_PIXELS");
        require(parsedContour.has_value(), "AttributeNames::parse CONTOUR_PIXELS must succeed");
        requireEqual(static_cast<int>(*parsedContour), static_cast<int>(ContourPixels), "AttributeNames::parse CONTOUR_PIXELS value");
        require(!AttributeNames::parse("NOT_AN_ATTRIBUTE").has_value(), "AttributeNames::parse invalid must fail");

        auto sampleLayout = NodeAttributeSampleLayout::create(2, {Area, MeanGrayLevel});
        requireEqual(sampleLayout.NUM_ATTRIBUTES, 10, "NodeAttributeSampleLayout count");
        requireEqual(sampleLayout.getIndex(Area, -2), 0, "AREA ancestor-2 offset");
        requireEqual(sampleLayout.getIndex(MeanGrayLevel, 0), 5, "MeanGrayLevel current-node offset");
        requireEqual(sampleLayout.linearIndex(1, Area, -2), 10, "AREA ancestor-2 linear index");
        requireEqual(sampleLayout.linearIndex(1, MeanGrayLevel, 2), 19, "MeanGrayLevel descendant-2 linear index");
        requireEqual(NodeAttributeSampleLayout::toString(Area, -2), std::string("AREA_ANCESTOR_2"),
                     "NodeAttributeSampleLayout ancestor label");
        requireEqual(NodeAttributeSampleLayout::toString(MeanGrayLevel, 2), std::string("MEAN_GRAY_LEVEL_DESCENDANT_2"),
                     "NodeAttributeSampleLayout descendant label");
        if constexpr (contract::validationsEnabled) {
            requireThrows<std::invalid_argument>([]() { static_cast<void>(NodeAttributeSampleLayout::create(1, {Area, Area})); },
                                                 "NodeAttributeSampleLayout::create must reject duplicate attributes");
        }
    }

    {
        requireComputerContract<AreaComputer>({Area}, AttributeComputerDomain::Topology, AttributeComputerFamily::Area, "AreaComputer");
        requireComputerContract<BoundingBoxComputer>(
            {BoxWidth, BoundingBoxHeight, DiagonalLength, Rectangularity, RatioWh, BoxColumnMin, BoxColumnMax, BoxRowMin, BoxRowMax},
            AttributeComputerDomain::Topology, AttributeComputerFamily::BoundingBox, "BoundingBoxComputer");
        requireComputerContract<TreeTopologyComputer>({SubtreeHeight, DepthNode, IsLeafNode, IsRootNode, NumChildrenNode, NumSiblingsNode,
                                                       NumDescendantsNode, NumLeafDescendantsNode, LeafRatioNode, BalanceNode, AvgChildHeightNode},
                                                      AttributeComputerDomain::Topology, AttributeComputerFamily::TreeTopology, "TreeTopologyComputer");
        requireComputerContract<CentralMomentsComputer>(
            {CentralMoment20, CentralMoment02, CentralMoment11, CentralMoment30, CentralMoment03, CentralMoment21, CentralMoment12},
            AttributeComputerDomain::Topology, AttributeComputerFamily::CentralMoments, "CentralMomentsComputer");
        requireComputerContract<HuMomentsComputer>({HuMoment1, HuMoment2, HuMoment3, HuMoment4, HuMoment5, HuMoment6, HuMoment7},
                                                   AttributeComputerDomain::Topology, AttributeComputerFamily::HuMoments, "HuMomentsComputer");
        requireComputerContract<MomentBasedAttributeComputer>(
            {Inertia, Compactness, Eccentricity, LengthMajorAxis, LengthMinorAxis, AxisOrientation, Circularity}, AttributeComputerDomain::Topology,
            AttributeComputerFamily::MomentDerived, "MomentBasedAttributeComputer");
        requireComputerContract<BitquadAttributeComputer>({BitquadArea, BitquadNumberEuler, BitquadNumberHoles, BitquadPerimeter,
                                                           BitquadPerimeterContinuous, BitquadCircularity, BitquadPerimeterAverage,
                                                           BitquadLengthAverage, BitquadWidthAverage},
                                                          AttributeComputerDomain::Topology, AttributeComputerFamily::Bitquad, "BitquadAttributeComputer");
        requireComputerContract<ContourSideAttributeComputer>(
            {ContourPixels, ContourPerimeter, ContourSideNorth, ContourSideWest, ContourSideEast, ContourSideSouth},
            AttributeComputerDomain::Topology, AttributeComputerFamily::ContourSide, "ContourSideAttributeComputer");
        requireComputerContract<VolumeComputer>({Volume, RelativeVolume}, AttributeComputerDomain::Altitude, AttributeComputerFamily::Volume,
                                                "VolumeComputer");
        requireComputerContract<GrayLevelStatsComputer>({MeanGrayLevel, GrayLevelVariance, GrayLevelHeight}, AttributeComputerDomain::Altitude,
                                                        AttributeComputerFamily::GrayLevelStats, "GrayLevelStatsComputer");
        requireComputerContract<MaxDistComputer>({MaxDist}, AttributeComputerDomain::Altitude, AttributeComputerFamily::MaxDist, "MaxDistComputer");

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

        requireEqual(static_cast<int>(mmcfilters::detail::familyForAttribute(Rectangularity)), static_cast<int>(AttributeComputerFamily::BoundingBox),
                     "scheduler family lookup for bounding boxes");
        requireEqual(static_cast<int>(mmcfilters::detail::familyForAttribute(Area)), static_cast<int>(AttributeComputerFamily::Area),
                     "scheduler family lookup for AREA");
        requireEqual(static_cast<int>(mmcfilters::detail::familyForAttribute(AvgChildHeightNode)), static_cast<int>(AttributeComputerFamily::TreeTopology),
                     "scheduler family lookup for tree topology");
        requireEqual(static_cast<int>(mmcfilters::detail::familyForAttribute(CentralMoment20)), static_cast<int>(AttributeComputerFamily::CentralMoments),
                     "scheduler family lookup for central moments");
        requireEqual(static_cast<int>(mmcfilters::detail::familyForAttribute(HuMoment7)), static_cast<int>(AttributeComputerFamily::HuMoments),
                     "scheduler family lookup for Hu moments");
        requireEqual(static_cast<int>(mmcfilters::detail::familyForAttribute(Eccentricity)), static_cast<int>(AttributeComputerFamily::MomentDerived),
                     "scheduler family lookup for moment-derived attributes");
        requireEqual(static_cast<int>(mmcfilters::detail::familyForAttribute(BitquadCircularity)), static_cast<int>(AttributeComputerFamily::Bitquad),
                     "scheduler family lookup for bitquads");
        requireEqual(static_cast<int>(mmcfilters::detail::familyForAttribute(ContourSideSouth)), static_cast<int>(AttributeComputerFamily::ContourSide),
                     "scheduler family lookup for contour sides");

        requireGlobalAttributeRegistryContracts();
    }

    {
        const std::array<Attribute, 3> scheduledRequest{MeanGrayLevel, Eccentricity, Rectangularity};
        const mmcfilters::detail::AttributeComputationPlan schedulerPlan =
            mmcfilters::detail::makeAttributeComputationPlan(std::span<const Attribute>(scheduledRequest));

        require(schedulerPlan.requests(MeanGrayLevel), "scheduler preserves requested MeanGrayLevel");
        require(schedulerPlan.requests(Eccentricity), "scheduler preserves requested ECCENTRICITY");
        require(schedulerPlan.requests(Rectangularity), "scheduler preserves requested RECTANGULARITY");
        require(schedulerPlan.materializes(Volume), "scheduler adds hidden VOLUME dependency");
        require(schedulerPlan.materializes(Area), "scheduler adds hidden AREA dependency");
        require(schedulerPlan.materializes(CentralMoment20), "scheduler adds hidden central moment dependency");
        require(schedulerPlan.materializes(CentralMoment02), "scheduler adds hidden central moment dependency");
        require(schedulerPlan.materializes(CentralMoment11), "scheduler adds hidden central moment dependency");
        require(schedulerPlan.hides(Volume), "scheduler marks VOLUME as hidden");
        require(schedulerPlan.hides(Area), "scheduler marks AREA as hidden");
        require(!schedulerPlan.hides(MeanGrayLevel), "scheduler does not hide requested attributes");

        const std::vector<Attribute> scheduledVolume = schedulerPlan.materializedForFamily(AttributeComputerFamily::Volume);
        requireEqual(scheduledVolume.size(), static_cast<std::size_t>(1), "scheduler volume family count");
        requireEqual(static_cast<int>(scheduledVolume.front()), static_cast<int>(Volume), "scheduler volume family attribute");

        const std::vector<Attribute> scheduledGray = schedulerPlan.requestedForFamily(AttributeComputerFamily::GrayLevelStats);
        requireEqual(scheduledGray.size(), static_cast<std::size_t>(1), "scheduler gray family count");
        requireEqual(static_cast<int>(scheduledGray.front()), static_cast<int>(MeanGrayLevel), "scheduler gray family attribute");

        const std::array<Attribute, 1> rectangularityOnly{Rectangularity};
        const std::array<Attribute, 1> widthOnly{BoxWidth};
        require(mmcfilters::detail::anyAttributeRequiresDependency(std::span<const Attribute>(rectangularityOnly), Area),
                "scheduler identifies RECTANGULARITY AREA dependency");
        require(!mmcfilters::detail::anyAttributeRequiresDependency(std::span<const Attribute>(widthOnly), Area),
                "scheduler does not add AREA to independent bounding-box attributes");
        requireEqual(static_cast<int>(mmcfilters::detail::familyForAttribute(Eccentricity)), static_cast<int>(AttributeComputerFamily::MomentDerived),
                     "scheduler family lookup for ECCENTRICITY");
        requireEqual(static_cast<int>(mmcfilters::detail::familyForAttribute(ContourSideSouth)), static_cast<int>(AttributeComputerFamily::ContourSide),
                     "scheduler family lookup for contour sides");

        requireGlobalDependencySchedulerContracts();
        requireCapabilityRegistryContracts();
    }

    {
        const AttributeNames areaNames = AttributeNames::fromList({Area});
        const AttributeNames volumeNames = AttributeNames::fromList({Volume, RelativeVolume});
        const std::array<float, 2> areaValues{4.0f, 8.0f};
        const std::array<float, 4> volumeValues{12.0f, 1.5f, 30.0f, 3.75f};
        const std::array<DependencySourceT<float>, 2> sources{{
            DependencySourceT<float>{&volumeNames, volumeValues.data()},
            DependencySourceT<float>{&areaNames, areaValues.data()},
        }};
        const DependencyResolver<float> resolver{std::span<const DependencySourceT<float>>(sources)};

        require(&resolver.require(Area) == &sources[1], "DependencyResolver resolves AREA by name");
        require(&resolver.require(Volume) == &sources[0], "DependencyResolver resolves VOLUME by name");
        require(&resolver.requireAll({Volume, RelativeVolume}) == &sources[0], "DependencyResolver resolves multi-attribute source");
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(resolver.require(MeanGrayLevel)); }, "DependencyResolver rejects missing named dependency");

        const std::array<DependencySourceT<float>, 1> invalidSources{{DependencySourceT<float>{nullptr, areaValues.data()}}};
        const DependencyResolver<float> invalidResolver{std::span<const DependencySourceT<float>>(invalidSources)};
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(invalidResolver.require(Area)); },
                                             "DependencyResolver rejects invalid dependency source");

        requireNear(::mmcfilters::attributes::numeric::safeDivide(6.0f, 3.0f), 2.0f, 0.0f, "safeDivide normal division");
        requireEqual(::mmcfilters::attributes::numeric::safeDivide(6.0f, 0.0f, -1.0f), -1.0f, "safeDivide returns finite fallback");
        requireEqual(::mmcfilters::attributes::numeric::safeSqrt(-4.0f), 0.0f, "safeSqrt clamps negative input");
        requireEqual(::mmcfilters::attributes::numeric::clampUpper(5.0f, 3.0f), 3.0f, "clampUpper caps values");
    }

    {
        auto tree = makeComponentTree(makeComponentTreeFixture(), true);

        const std::array<Attribute, 1> widthRequest{BoxWidth};
        const AttributeNames widthNames = AttributeNames::fromList({BoxWidth});
        std::vector<float> widthBuffer(static_cast<std::size_t>(tree->numInternalNodeSlots()) * static_cast<std::size_t>(widthNames.NUM_ATTRIBUTES), 0.0f);
        BoundingBoxComputer::compute(
            AttributeComputeContext<float>{*tree, std::span<float>(widthBuffer), widthNames, std::span<const Attribute>(widthRequest)});
        requireEqual(widthBuffer[widthNames.linearIndex(tree->root(), BoxWidth)], 4.0f, "BOX_WIDTH context compute does not require AREA dependency");

        if constexpr (contract::validationsEnabled) {
            const std::array<Attribute, 1> rectangularityRequest{Rectangularity};
            const AttributeNames rectangularityNames = AttributeNames::fromList({Rectangularity});
            std::vector<float> rectangularityBuffer(
                static_cast<std::size_t>(tree->numInternalNodeSlots()) * static_cast<std::size_t>(rectangularityNames.NUM_ATTRIBUTES), 0.0f);
            requireThrows<std::invalid_argument>(
                [&]() {
                    BoundingBoxComputer::compute(AttributeComputeContext<float>{*tree, std::span<float>(rectangularityBuffer), rectangularityNames,
                                                                                std::span<const Attribute>(rectangularityRequest)});
                },
                "RECTANGULARITY context compute must reject missing AREA dependency");
        }
    }

    {
        auto valuedTree = makeValuedComponentTree(makeComponentTreeFixture(), true);
        const MorphologicalTree& tree = valuedTree->topology();
        const std::span<const std::uint8_t> altitude = valuedTree->nodeAltitudeSpan();

        const std::array<Attribute, 1> volumeRequest{Volume};
        const AttributeNames volumeNames = AttributeNames::fromList({Volume});
        std::vector<float> volumeBuffer(static_cast<std::size_t>(tree.numInternalNodeSlots()) * static_cast<std::size_t>(volumeNames.NUM_ATTRIBUTES), 0.0f);
        VolumeComputer::compute(AltitudeAttributeComputeContext<float, std::uint8_t>{tree, altitude, std::span<float>(volumeBuffer), volumeNames,
                                                                                     std::span<const Attribute>(volumeRequest)});
        requireEqual(volumeBuffer[volumeNames.linearIndex(tree.root(), Volume)], 42.0f, "VOLUME context compute does not require AREA dependency");

        if constexpr (contract::validationsEnabled) {
            const std::array<Attribute, 1> relativeVolumeRequest{RelativeVolume};
            const AttributeNames relativeVolumeNames = AttributeNames::fromList({RelativeVolume});
            std::vector<float> relativeVolumeBuffer(
                static_cast<std::size_t>(tree.numInternalNodeSlots()) * static_cast<std::size_t>(relativeVolumeNames.NUM_ATTRIBUTES), 0.0f);
            requireThrows<std::invalid_argument>(
                [&]() {
                    VolumeComputer::compute(AltitudeAttributeComputeContext<float, std::uint8_t>{
                        tree, altitude, std::span<float>(relativeVolumeBuffer), relativeVolumeNames, std::span<const Attribute>(relativeVolumeRequest)});
                },
                "RELATIVE_VOLUME context compute must reject missing AREA dependency");
        }

        const std::array<Attribute, 1> grayHeightRequest{GrayLevelHeight};
        const AttributeNames grayHeightNames = AttributeNames::fromList({GrayLevelHeight});
        std::vector<float> grayHeightBuffer(static_cast<std::size_t>(tree.numInternalNodeSlots()) * static_cast<std::size_t>(grayHeightNames.NUM_ATTRIBUTES),
                                            0.0f);
        GrayLevelStatsComputer::compute(AltitudeAttributeComputeContext<float, std::uint8_t>{tree, altitude, std::span<float>(grayHeightBuffer),
                                                                                             grayHeightNames, std::span<const Attribute>(grayHeightRequest)});
        require(grayHeightBuffer[grayHeightNames.linearIndex(tree.root(), GrayLevelHeight)] >= 0.0f,
                "GrayLevelHeight context compute does not require aggregate dependencies");

        if constexpr (contract::validationsEnabled) {
            const std::array<Attribute, 1> meanRequest{MeanGrayLevel};
            const AttributeNames meanNames = AttributeNames::fromList({MeanGrayLevel});
            std::vector<float> meanBuffer(static_cast<std::size_t>(tree.numInternalNodeSlots()) * static_cast<std::size_t>(meanNames.NUM_ATTRIBUTES), 0.0f);
            requireThrows<std::invalid_argument>(
                [&]() {
                    GrayLevelStatsComputer::compute(AltitudeAttributeComputeContext<float, std::uint8_t>{tree, altitude, std::span<float>(meanBuffer), meanNames,
                                                                                                         std::span<const Attribute>(meanRequest)});
                },
                "MeanGrayLevel context compute must reject missing VOLUME and AREA dependencies");

            const std::array<Attribute, 1> huRequest{HuMoment1};
            const AttributeNames huNames = AttributeNames::fromList({HuMoment1});
            std::vector<float> huBuffer(static_cast<std::size_t>(tree.numInternalNodeSlots()) * static_cast<std::size_t>(huNames.NUM_ATTRIBUTES), 0.0f);
            requireThrows<std::invalid_argument>(
                [&]() {
                    HuMomentsComputer::compute(AttributeComputeContext<float>{tree, std::span<float>(huBuffer), huNames, std::span<const Attribute>(huRequest)});
                },
                "HU_MOMENT context compute must reject missing central moment and AREA dependencies");

            const std::array<Attribute, 1> inertiaRequest{Inertia};
            const AttributeNames inertiaNames = AttributeNames::fromList({Inertia});
            std::vector<float> inertiaBuffer(
                static_cast<std::size_t>(tree.numInternalNodeSlots()) * static_cast<std::size_t>(inertiaNames.NUM_ATTRIBUTES), 0.0f);
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
