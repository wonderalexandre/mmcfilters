#include "support/TestSupport.hpp"

#include "mmcfilters/attributes/computers/BoundingBoxComputer.hpp"
#include "mmcfilters/attributes/computers/GrayLevelStatsComputer.hpp"
#include "mmcfilters/attributes/computers/MomentBasedAttributeComputer.hpp"
#include "mmcfilters/attributes/computers/VolumeComputer.hpp"

#include <array>
#include <optional>
#include <span>

using namespace mmcfilters;
using namespace mmcfilters::attributes::computers;
using namespace mmcfilters::unit_tests;

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
                require(
                    names.linearIndex(nodeIndex, attribute) < 6 * names.NUM_ATTRIBUTES,
                    "AttributeNames dense linear index must stay within buffer bounds");
            }
        }
        requireThrows<std::invalid_argument>(
            []() { static_cast<void>(AttributeNames::fromList({AREA, AREA, VOLUME})); },
            "AttributeNames::fromList must reject duplicate attributes");
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
        requireThrows<std::invalid_argument>(
            []() { static_cast<void>(AttributeNamesWithDelta::create(1, {AREA, AREA})); },
            "AttributeNamesWithDelta::create must reject duplicate attributes");
    }

    {
        auto tree = makeComponentTree(makeComponentTreeFixture(), true);
        BoundingBoxComputer boxComputer;

        const std::array<Attribute, 1> widthRequest{BOX_WIDTH};
        const AttributeNames widthNames = AttributeNames::fromList({BOX_WIDTH});
        std::vector<float> widthBuffer(
            static_cast<std::size_t>(tree->getNumInternalNodeSlots()) *
                static_cast<std::size_t>(widthNames.NUM_ATTRIBUTES),
            0.0f);
        boxComputer.compute(
            *tree,
            AttributeAltitudeView{},
            widthBuffer,
            widthNames,
            std::span<const Attribute>(widthRequest),
            {});
        requireEqual(widthBuffer[widthNames.linearIndex(tree->getRoot(), BOX_WIDTH)], 4.0f, "BOX_WIDTH direct compute does not require AREA dependency");

        const std::array<Attribute, 1> rectangularityRequest{RECTANGULARITY};
        const AttributeNames rectangularityNames = AttributeNames::fromList({RECTANGULARITY});
        std::vector<float> rectangularityBuffer(
            static_cast<std::size_t>(tree->getNumInternalNodeSlots()) *
                static_cast<std::size_t>(rectangularityNames.NUM_ATTRIBUTES),
            0.0f);
        requireThrows<std::invalid_argument>(
            [&]() {
                boxComputer.compute(
                    *tree,
                    AttributeAltitudeView{},
                    rectangularityBuffer,
                    rectangularityNames,
                    std::span<const Attribute>(rectangularityRequest),
                    {});
            },
            "RECTANGULARITY direct compute must reject missing AREA dependency");
    }

    {
        auto weighted = makeWeightedComponentTree(makeComponentTreeFixture(), true);
        const MorphologicalTree& tree = weighted->topology();
        const AttributeAltitudeView altitudeView = makeAttributeAltitudeView(weighted->getAltitudeBuffer());

        VolumeComputer volumeComputer;
        const std::array<Attribute, 1> volumeRequest{VOLUME};
        const AttributeNames volumeNames = AttributeNames::fromList({VOLUME});
        std::vector<float> volumeBuffer(
            static_cast<std::size_t>(tree.getNumInternalNodeSlots()) *
                static_cast<std::size_t>(volumeNames.NUM_ATTRIBUTES),
            0.0f);
        volumeComputer.compute(
            tree,
            altitudeView,
            volumeBuffer,
            volumeNames,
            std::span<const Attribute>(volumeRequest),
            {});
        requireEqual(volumeBuffer[volumeNames.linearIndex(tree.getRoot(), VOLUME)], 42.0f, "VOLUME direct compute does not require AREA dependency");

        const std::array<Attribute, 1> relativeVolumeRequest{RELATIVE_VOLUME};
        const AttributeNames relativeVolumeNames = AttributeNames::fromList({RELATIVE_VOLUME});
        std::vector<float> relativeVolumeBuffer(
            static_cast<std::size_t>(tree.getNumInternalNodeSlots()) *
                static_cast<std::size_t>(relativeVolumeNames.NUM_ATTRIBUTES),
            0.0f);
        requireThrows<std::invalid_argument>(
            [&]() {
                volumeComputer.compute(
                    tree,
                    altitudeView,
                    relativeVolumeBuffer,
                    relativeVolumeNames,
                    std::span<const Attribute>(relativeVolumeRequest),
                    {});
            },
            "RELATIVE_VOLUME direct compute must reject missing AREA dependency");

        GrayLevelStatsComputer grayComputer;
        const std::array<Attribute, 1> levelRequest{LEVEL};
        const AttributeNames levelNames = AttributeNames::fromList({LEVEL});
        std::vector<float> levelBuffer(
            static_cast<std::size_t>(tree.getNumInternalNodeSlots()) *
                static_cast<std::size_t>(levelNames.NUM_ATTRIBUTES),
            0.0f);
        grayComputer.compute(
            tree,
            altitudeView,
            levelBuffer,
            levelNames,
            std::span<const Attribute>(levelRequest),
            {});
        requireEqual(
            levelBuffer[levelNames.linearIndex(tree.getRoot(), LEVEL)],
            static_cast<float>(weighted->getAltitude(tree.getRoot())),
            "LEVEL direct compute does not require aggregate dependencies");

        const std::array<Attribute, 1> grayHeightRequest{GRAY_HEIGHT};
        const AttributeNames grayHeightNames = AttributeNames::fromList({GRAY_HEIGHT});
        std::vector<float> grayHeightBuffer(
            static_cast<std::size_t>(tree.getNumInternalNodeSlots()) *
                static_cast<std::size_t>(grayHeightNames.NUM_ATTRIBUTES),
            0.0f);
        grayComputer.compute(
            tree,
            altitudeView,
            grayHeightBuffer,
            grayHeightNames,
            std::span<const Attribute>(grayHeightRequest),
            {});
        require(
            grayHeightBuffer[grayHeightNames.linearIndex(tree.getRoot(), GRAY_HEIGHT)] >= 0.0f,
            "GRAY_HEIGHT direct compute does not require aggregate dependencies");

        const std::array<Attribute, 1> meanRequest{MEAN_LEVEL};
        const AttributeNames meanNames = AttributeNames::fromList({MEAN_LEVEL});
        std::vector<float> meanBuffer(
            static_cast<std::size_t>(tree.getNumInternalNodeSlots()) *
                static_cast<std::size_t>(meanNames.NUM_ATTRIBUTES),
            0.0f);
        requireThrows<std::invalid_argument>(
            [&]() {
                grayComputer.compute(
                    tree,
                    altitudeView,
                    meanBuffer,
                    meanNames,
                    std::span<const Attribute>(meanRequest),
                    {});
            },
            "MEAN_LEVEL direct compute must reject missing VOLUME and AREA dependencies");

        HuMomentsComputer huComputer;
        const std::array<Attribute, 1> huRequest{HU_MOMENT_1};
        const AttributeNames huNames = AttributeNames::fromList({HU_MOMENT_1});
        std::vector<float> huBuffer(
            static_cast<std::size_t>(tree.getNumInternalNodeSlots()) *
                static_cast<std::size_t>(huNames.NUM_ATTRIBUTES),
            0.0f);
        requireThrows<std::invalid_argument>(
            [&]() {
                huComputer.compute(
                    tree,
                    AttributeAltitudeView{},
                    huBuffer,
                    huNames,
                    std::span<const Attribute>(huRequest),
                    {});
            },
            "HU_MOMENT direct compute must reject missing central moment and AREA dependencies");

        MomentBasedAttributeComputer momentBasedComputer;
        const std::array<Attribute, 1> inertiaRequest{INERTIA};
        const AttributeNames inertiaNames = AttributeNames::fromList({INERTIA});
        std::vector<float> inertiaBuffer(
            static_cast<std::size_t>(tree.getNumInternalNodeSlots()) *
                static_cast<std::size_t>(inertiaNames.NUM_ATTRIBUTES),
            0.0f);
        requireThrows<std::invalid_argument>(
            [&]() {
                momentBasedComputer.compute(
                    tree,
                    AttributeAltitudeView{},
                    inertiaBuffer,
                    inertiaNames,
                    std::span<const Attribute>(inertiaRequest),
                    {});
            },
            "moment-based direct compute must reject missing central moment and AREA dependencies");
    }

    return 0;
}
