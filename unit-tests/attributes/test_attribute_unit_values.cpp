#include "support/TestSupport.hpp"

#include "mmcfilters/attributes/AttributeFactory.hpp"

#include <cmath>
#include <sstream>

using namespace mmcfilters;
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

void requireUnitComputerMatchesSingletonOracle(double adjacencyRadius) {
    auto image = makeImage(1, 1, {7});
    auto weighted = makeWeightedComponentTree(image, true, adjacencyRadius);
    const MorphologicalTree& topology = weighted->topology();
    const std::vector<NodeId> unitProperParts{0};

    for (Attribute attribute : ATTRIBUTE_GROUPS.at(AttributeGroup::ALL)) {
        const AttributeComputer& computer = AttributeFactory::create(attribute);
        const auto [oracleNames, oracleValues] =
            AttributeComputedIncrementally::computeSingleAttribute(*weighted, attribute);
        const float expected = oracleValues[oracleNames.linearIndex(topology.getRoot(), attribute)];

        const AttributeNames unitNames = makeDenseAttributeNames({attribute});
        std::vector<float> unitValues(
            unitProperParts.size() * static_cast<size_t>(unitNames.NUM_ATTRIBUTES),
            std::numeric_limits<float>::quiet_NaN());
        const std::array<Attribute, 1> requested{attribute};

        computer.computeUnitAttributes(
            topology,
            &weighted->getAltitudeBuffer(),
            unitProperParts,
            unitValues,
            unitNames,
            std::span<const Attribute>(requested));

        std::ostringstream label;
        label << "unit value for " << AttributeNames::toString(attribute)
              << " at radius " << adjacencyRadius;
        requireFloatEquivalent(unitValues[unitNames.linearIndex(0, attribute)], expected, label.str());
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

    for (const Attribute attribute : attributes) {
        const AttributeComputer& computer = AttributeFactory::create(attribute);
        const std::array<Attribute, 1> requested{attribute};
        computer.computeUnitAttributes(
            topology,
            &weighted->getAltitudeBuffer(),
            exportedProperParts,
            unitValues,
            unitNames,
            std::span<const Attribute>(requested));
    }

    for (NodeId leafIndex = 0; leafIndex < static_cast<NodeId>(exportedProperParts.size()); ++leafIndex) {
        const NodeId properPart = exportedProperParts[static_cast<size_t>(leafIndex)];
        const NodeId owner = topology.getSmallestComponent(properPart);
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
    requireUnitComputerMatchesSingletonOracle(1.0);
    requireUnitComputerMatchesSingletonOracle(1.5);
    requireUnitComputersUseProvidedLeafOrder();
    return 0;
}
