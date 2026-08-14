#include "support/TestSupport.hpp"

#include "mmcfilters/attributes/computers/BoundingBoxComputer.hpp"
#include "mmcfilters/attributes/computers/GrayLevelStatsComputer.hpp"
#include "mmcfilters/attributes/computers/VolumeComputer.hpp"

#include <array>
#include <cmath>
#include <cstdint>
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
    auto valuedTree = makeValuedComponentTree(image, true, adjacencyRadius);
    const MorphologicalTree& topology = valuedTree->topology();

    for (Attribute attribute : ATTRIBUTE_GROUPS.at(AttributeGroup::All)) {
        const auto computed = AttributeComputation::computeSingleAttribute(*valuedTree, attribute);
        const auto projected = AttributeComputation::projectNodeValuesToExportedHigra(*valuedTree, computed.attributeNames(), computed.values());
        const float expected = computed.values()[computed.attributeNames().linearIndex(topology.root(), attribute)];

        std::ostringstream label;
        label << "unit value for " << AttributeNames::toString(attribute) << " at radius " << adjacencyRadius;
        requireFloatEquivalent(projected[computed.attributeNames().linearIndex(0, attribute)], expected, label.str());
    }
}

void requireUnitComputersUseProvidedLeafOrder() {
    auto valuedTree = makeValuedComponentTree(makeComponentTreeFixture(), true);
    const MorphologicalTree& topology = valuedTree->topology();
    const std::vector<PixelId> exportedProperParts{10, 0};
    const std::vector<Attribute> attributes{MeanGrayLevel, Volume, BoxColumnMin, BoxRowMin};
    const AttributeNames unitNames = makeDenseAttributeNames(attributes);
    std::vector<float> unitValues(exportedProperParts.size() * static_cast<size_t>(unitNames.NUM_ATTRIBUTES), std::numeric_limits<float>::quiet_NaN());

    const std::span<const std::uint8_t> altitude = valuedTree->nodeAltitudeSpan();

    const std::array<Attribute, 1> meanGrayLevelRequest{MeanGrayLevel};
    GrayLevelStatsComputer::computeUnitRows(AltitudeUnitAttributeComputeContext<float, std::uint8_t>{
        topology, altitude, std::span<const PixelId>(exportedProperParts), std::span<float>(unitValues), unitNames,
        std::span<const Attribute>(meanGrayLevelRequest)});

    const std::array<Attribute, 1> volumeRequest{Volume};
    VolumeComputer::computeUnitRows(AltitudeUnitAttributeComputeContext<float, std::uint8_t>{
        topology, altitude, std::span<const PixelId>(exportedProperParts), std::span<float>(unitValues), unitNames, std::span<const Attribute>(volumeRequest)});

    const std::array<Attribute, 2> boundingBoxRequest{BoxColumnMin, BoxRowMin};
    BoundingBoxComputer::computeUnitRows(UnitAttributeComputeContext<float>{
        topology, std::span<const PixelId>(exportedProperParts), std::span<float>(unitValues), unitNames, std::span<const Attribute>(boundingBoxRequest)});

    for (NodeId leafIndex = 0; leafIndex < static_cast<NodeId>(exportedProperParts.size()); ++leafIndex) {
        const PixelId pixel = exportedProperParts[static_cast<size_t>(leafIndex)];
        const NodeId smallestNodeId = topology.smallestNode(pixel);
        const auto [row, column] = ImageUtils::to2D(pixel, topology.numColumns());
        requireEqual(unitValues[unitNames.linearIndex(leafIndex, MeanGrayLevel)], static_cast<float>(valuedTree->nodeAltitude(smallestNodeId)),
                     "unit MEAN_GRAY_LEVEL must follow the exported leaf order");
        requireEqual(unitValues[unitNames.linearIndex(leafIndex, Volume)], static_cast<float>(valuedTree->nodeAltitude(smallestNodeId)),
                     "unit VOLUME must follow the exported leaf order");
        requireEqual(unitValues[unitNames.linearIndex(leafIndex, BoxColumnMin)], static_cast<float>(column),
                     "unit BOX_COLUMN_MIN must follow the exported leaf order");
        requireEqual(unitValues[unitNames.linearIndex(leafIndex, BoxRowMin)], static_cast<float>(row),
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
