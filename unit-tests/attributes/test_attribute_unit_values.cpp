#include "support/TestSupport.hpp"

#include "mmcfilters/attributes/computers/BoundingBoxComputer.hpp"
#include "mmcfilters/attributes/computers/GrayLevelStatsComputer.hpp"
#include "mmcfilters/attributes/computers/VolumeComputer.hpp"

#include <array>
#include <cmath>
#include <span>
#include <sstream>

using namespace mmcfilters;
using namespace mmcfilters::attributes::computers;
using namespace mmcfilters::unit_tests;

namespace {

AttributeNames makeDenseAttributeNames(const std::vector<Attribute>& attributes) {
    std::unordered_map<Attribute, int> offsets;
    for (int index = 0; index < static_cast<int>(attributes.size()); ++index) {
        offsets[attributes[static_cast<size_t>(index)]] = index;
    }
    return AttributeNames(std::move(offsets));
}

void requireFloatEquivalent(float actual, float expected, const std::string& label) {
    if (std::isnan(expected)) {
        require(std::isnan(actual), label + " expected NaN");
        return;
    }
    if (std::isinf(expected)) {
        require(std::isinf(actual) && std::signbit(actual) == std::signbit(expected), label + " expected infinity");
        return;
    }
    requireNear(actual, expected, 1.0e-5f, label);
}

void requireExportedHigraUnitProjectionMatchesOnePixelOracle(double adjacencyRadius) {
    auto image = makeImage(1, 1, {7});
    auto weighted = makeWeightedComponentTree(image, true, adjacencyRadius);
    const MorphologicalTree& topology = weighted->topology();

    for (Attribute attribute : ATTRIBUTE_GROUPS.at(AttributeGroup::ALL)) {
        const auto computed =
            AttributeComputation::computeSingleAttribute(*weighted, attribute);
        const auto projected = AttributeComputation::projectNodeValuesToExportedHigra(
            *weighted,
            computed.attributeNames(),
            computed.values());
        const float expected = computed.values()[computed.attributeNames().linearIndex(topology.getRoot(), attribute)];

        std::ostringstream label;
        label << "unit value for " << AttributeNames::toString(attribute)
              << " at radius " << adjacencyRadius;
        requireFloatEquivalent(projected[computed.attributeNames().linearIndex(0, attribute)], expected, label.str());
    }
}

void requireUnitComputersUseProvidedLeafOrder() {
    auto weighted = makeWeightedComponentTree(makeComponentTreeFixture(), true);
    const MorphologicalTree& topology = weighted->topology();
    const std::vector<NodeId> exportedProperParts{10, 0};
    const std::vector<Attribute> attributes{LEVEL, VOLUME, BOX_COL_MIN, BOX_ROW_MIN};
    const AttributeNames unitNames = makeDenseAttributeNames(attributes);
    std::vector<float> unitValues(
        exportedProperParts.size() * static_cast<size_t>(unitNames.NUM_ATTRIBUTES),
        std::numeric_limits<float>::quiet_NaN());

    GrayLevelStatsComputer grayComputer;
    VolumeComputer volumeComputer;
    BoundingBoxComputer boundingBoxComputer;
    const AttributeAltitudeView altitude = makeAttributeAltitudeView(weighted->getAltitudeBuffer());

    const std::array<Attribute, 1> levelRequest{LEVEL};
    grayComputer.computeUnitAttributes(
        topology,
        altitude,
        exportedProperParts,
        unitValues,
        unitNames,
        std::span<const Attribute>(levelRequest));

    const std::array<Attribute, 1> volumeRequest{VOLUME};
    volumeComputer.computeUnitAttributes(
        topology,
        altitude,
        exportedProperParts,
        unitValues,
        unitNames,
        std::span<const Attribute>(volumeRequest));

    const std::array<Attribute, 2> boundingBoxRequest{BOX_COL_MIN, BOX_ROW_MIN};
    boundingBoxComputer.computeUnitAttributes(
        topology,
        AttributeAltitudeView{},
        exportedProperParts,
        unitValues,
        unitNames,
        std::span<const Attribute>(boundingBoxRequest));

    for (NodeId leafIndex = 0; leafIndex < static_cast<NodeId>(exportedProperParts.size()); ++leafIndex) {
        const NodeId properPart = exportedProperParts[static_cast<size_t>(leafIndex)];
        const NodeId owner = topology.getProperPartOwner(properPart);
        const auto [row, col] = ImageUtils::to2D(properPart, topology.getNumColsOfImage());
        requireEqual(
            unitValues[unitNames.linearIndex(leafIndex, LEVEL)],
            static_cast<float>(weighted->getAltitude(owner)),
            "unit LEVEL must follow the exported leaf order");
        requireEqual(
            unitValues[unitNames.linearIndex(leafIndex, VOLUME)],
            static_cast<float>(weighted->getAltitude(owner)),
            "unit VOLUME must follow the exported leaf order");
        requireEqual(
            unitValues[unitNames.linearIndex(leafIndex, BOX_COL_MIN)],
            static_cast<float>(col),
            "unit BOX_COL_MIN must follow the exported leaf order");
        requireEqual(
            unitValues[unitNames.linearIndex(leafIndex, BOX_ROW_MIN)],
            static_cast<float>(row),
            "unit BOX_ROW_MIN must follow the exported leaf order");
    }
}

} // namespace

int main() {
    requireExportedHigraUnitProjectionMatchesOnePixelOracle(1.0);
    requireExportedHigraUnitProjectionMatchesOnePixelOracle(1.5);
    requireUnitComputersUseProvidedLeafOrder();
    return 0;
}
