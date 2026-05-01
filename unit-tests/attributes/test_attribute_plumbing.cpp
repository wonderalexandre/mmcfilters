#include "support/TestSupport.hpp"

#include "mmcfilters/attributes/AreaComputer.hpp"
#include "mmcfilters/attributes/AttributeFactory.hpp"
#include "mmcfilters/attributes/BitquadsComputer.hpp"
#include "mmcfilters/attributes/BoundingBoxComputer.hpp"
#include "mmcfilters/attributes/GrayLevelStatsComputer.hpp"
#include "mmcfilters/attributes/MaxDistComputer.hpp"
#include "mmcfilters/attributes/MomentBasedAttributeComputer.hpp"
#include "mmcfilters/attributes/TreeTopologyComputer.hpp"
#include "mmcfilters/attributes/VolumeComputer.hpp"

#include <optional>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

namespace {

template <class T>
void requireInstanceOf(const AttributeComputer& computer, const std::string& label) {
    require(dynamic_cast<const T*>(&computer) != nullptr, label);
}

} // namespace

int main() {
    {
        const auto& areaComputer = AttributeFactory::create(AREA);
        requireInstanceOf<AreaComputer>(areaComputer, "AREA must map to AreaComputer");
        requireEqual(static_cast<int>(areaComputer.attributes().size()), 1, "AreaComputer attribute count");
        requireEqual(static_cast<int>(areaComputer.attributes().front()), static_cast<int>(AREA), "AreaComputer AREA attribute");
        require(areaComputer.requiredAttributes().empty(), "AreaComputer has no dependencies");

        const auto& volumeComputer = AttributeFactory::create(VOLUME);
        requireInstanceOf<VolumeComputer>(volumeComputer, "VOLUME must map to VolumeComputer");
        requireEqual(static_cast<int>(volumeComputer.attributes().size()), 2, "VolumeComputer attribute count");
        requireEqual(static_cast<int>(volumeComputer.attributes()[0]), static_cast<int>(VOLUME), "VolumeComputer first attribute");
        requireEqual(static_cast<int>(volumeComputer.attributes()[1]), static_cast<int>(RELATIVE_VOLUME), "VolumeComputer second attribute");
        requireEqual(static_cast<int>(volumeComputer.requiredAttributes().size()), 1, "VolumeComputer dependency count");

        const auto& grayComputer = AttributeFactory::create(LEVEL);
        requireInstanceOf<GrayLevelStatsComputer>(grayComputer, "LEVEL must map to GrayLevelStatsComputer");

        const auto& boxComputer = AttributeFactory::create(AttributeGroup::BOUNDING_BOX);
        requireInstanceOf<BoundingBoxComputer>(boxComputer, "BOUNDING_BOX group must map to BoundingBoxComputer");

        const auto& momentsComputer = AttributeFactory::create(HU_MOMENT_5);
        requireInstanceOf<HuMomentsComputer>(momentsComputer, "Hu moments must map to HuMomentsComputer");

        const auto& momentBasedComputer = AttributeFactory::create(COMPACTNESS);
        requireInstanceOf<MomentBasedAttributeComputer>(momentBasedComputer, "moment-based attributes must map to MomentBasedAttributeComputer");

        const auto& topologyComputer = AttributeFactory::create(AttributeGroup::TREE_TOPOLOGY);
        requireInstanceOf<TreeTopologyComputer>(topologyComputer, "TREE_TOPOLOGY group must map to TreeTopologyComputer");

        const auto& bitquadsComputer = AttributeFactory::create(BITQUADS_PERIMETER);
        requireInstanceOf<BitquadsComputer>(bitquadsComputer, "BitQuads attributes must map to BitquadsComputer");

        const auto& maxDistComputer = AttributeFactory::create(MAX_DIST);
        requireInstanceOf<MaxDistComputer>(maxDistComputer, "MAX_DIST must map to MaxDistComputer");

        bool allGroupRejected = false;
        try {
            (void)AttributeFactory::create(AttributeGroup::ALL);
        } catch (const std::runtime_error&) {
            allGroupRejected = true;
        }
        require(allGroupRejected, "AttributeFactory must reject ALL group as a direct computer target");
    }

    {
        auto names = AttributeNames::fromList(6, {AREA, VOLUME, MAX_DIST});
        requireEqual(names.NUM_ATTRIBUTES, 3, "AttributeNames count");
        requireEqual(names.getIndex(AREA), 0, "AREA index in AttributeNames");
        requireEqual(names.getIndex(VOLUME), 6, "VOLUME index in AttributeNames");
        requireEqual(names.linearIndex(2, MAX_DIST), 18, "MAX_DIST linear index");
        requireEqual(AttributeNames::toString(MAX_DIST), std::string("MAX_DIST"), "AttributeNames::toString");
        require(AttributeNames::describe(AREA).starts_with("Area:"), "AttributeNames::describe AREA");

        auto parsedArea = AttributeNames::parse("AREA");
        require(parsedArea.has_value(), "AttributeNames::parse AREA must succeed");
        requireEqual(static_cast<int>(*parsedArea), static_cast<int>(AREA), "AttributeNames::parse AREA value");
        require(!AttributeNames::parse("NOT_AN_ATTRIBUTE").has_value(), "AttributeNames::parse invalid must fail");

        auto deltaNames = AttributeNamesWithDelta::create(2, {AREA, LEVEL});
        requireEqual(deltaNames.NUM_ATTRIBUTES, 10, "AttributeNamesWithDelta count");
        requireEqual(deltaNames.getIndex(AREA, -2), 0, "AREA asc2 offset");
        requireEqual(deltaNames.getIndex(LEVEL, 0), 5, "LEVEL center offset");
        requireEqual(deltaNames.linearIndex(1, AREA, -2), 10, "AREA asc2 linear index");
        requireEqual(deltaNames.linearIndex(1, LEVEL, 2), 19, "LEVEL desc2 linear index");
        requireEqual(AttributeNamesWithDelta::toString(AREA, -2), std::string("AREA_ASC_2"), "AttributeNamesWithDelta asc label");
        requireEqual(AttributeNamesWithDelta::toString(LEVEL, 2), std::string("LEVEL_DESC_2"), "AttributeNamesWithDelta desc label");
    }

    return 0;
}
