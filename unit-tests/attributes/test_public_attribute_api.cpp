#include <mmcfilters/attributes/Attributes.hpp>
#include <mmcfilters/trees/MorphologicalTreeFactory.hpp>
#include <mmcfilters/utils/Image.hpp>

#include <array>
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

template<class T, class U>
void requireEqual(const T& actual, const U& expected, const std::string& message) {
    if (!(actual == expected)) {
        throw std::runtime_error(message);
    }
}

void requireHiddenDependenciesAreAbsent(
    const mmcfilters::AttributeNames& names,
    std::initializer_list<mmcfilters::Attribute> hiddenDependencies,
    const std::string& label)
{
    for (mmcfilters::Attribute dependency : hiddenDependencies) {
        require(!names.contains(dependency), label + " must not expose hidden dependency");
    }
}

void requireFiniteLiveNodeValues(
    const mmcfilters::MorphologicalTree& tree,
    const mmcfilters::AttributeNames& names,
    const std::vector<float>& values,
    mmcfilters::Attribute attribute,
    const std::string& label)
{
    for (mmcfilters::NodeId nodeId : tree.getAliveNodeIds()) {
        const float value = values[names.linearIndex(nodeId, attribute)];
        require(std::isfinite(value), label + " must be finite for every live node");
    }
}

} // namespace

int main() {
    using namespace mmcfilters;

    std::array<std::uint8_t, 16> pixels{
        3, 3, 2, 2,
        3, 4, 4, 2,
        1, 4, 5, 2,
        1, 1, 5, 0};
    auto image = ImageUInt8::fromExternal(pixels.data(), 4, 4);
    auto weighted = MorphologicalTreeFactory::createMaxTree(image, 1.5);

    auto weightedAttributes = AttributeComputation::computeAttributes(
        weighted,
        std::vector<AttributeOrGroup>{AREA, LEVEL, MAX_DIST});
    requireEqual(weightedAttributes.first.NUM_ATTRIBUTES, 3, "weighted public attribute stride");
    require(
        weightedAttributes.second.size() ==
            static_cast<std::size_t>(weighted.topology().getNumInternalNodeSlots()) *
                static_cast<std::size_t>(weightedAttributes.first.NUM_ATTRIBUTES),
        "weighted public attribute buffer shape");
    requireEqual(
        weightedAttributes.second[weightedAttributes.first.linearIndex(weighted.topology().getRoot(), AREA)],
        16.0f,
        "root AREA through public weighted facade");

    auto topologyAttributes = AttributeComputation::computeTopologyAttributes(
        weighted.topology(),
        std::vector<AttributeOrGroup>{AREA, BOX_WIDTH, BALANCE_NODE});
    requireEqual(topologyAttributes.first.NUM_ATTRIBUTES, 3, "topology public attribute stride");
    requireEqual(
        topologyAttributes.second[topologyAttributes.first.linearIndex(weighted.topology().getRoot(), AREA)],
        16.0f,
        "root AREA through explicit topology facade");

    auto deltaLevel = AttributeComputation::computeSingleAttributeWithDelta(
        weighted,
        LEVEL,
        AltitudeDiff<std::uint8_t>{1},
        1);
    requireEqual(deltaLevel.first.NUM_ATTRIBUTES, 3, "delta public attribute stride");
    require(
        deltaLevel.second.size() ==
            static_cast<std::size_t>(weighted.topology().getNumInternalNodeSlots()) *
                static_cast<std::size_t>(deltaLevel.first.NUM_ATTRIBUTES),
        "delta public attribute buffer shape");

    const NodeId root = weighted.topology().getRoot();
    requireEqual(
        deltaLevel.second[deltaLevel.first.linearIndex(root, LEVEL, 0)],
        static_cast<float>(weighted.getAltitude(root)),
        "delta center LEVEL through public weighted facade");

    auto meanOnly = AttributeComputation::computeAttributes(
        weighted,
        std::vector<AttributeOrGroup>{MEAN_LEVEL});
    requireEqual(meanOnly.first.NUM_ATTRIBUTES, 1, "MEAN_LEVEL request exposes one public attribute");
    require(meanOnly.first.contains(MEAN_LEVEL), "MEAN_LEVEL request exposes requested attribute");
    requireHiddenDependenciesAreAbsent(meanOnly.first, {AREA, VOLUME}, "MEAN_LEVEL request");
    requireFiniteLiveNodeValues(weighted.topology(), meanOnly.first, meanOnly.second, MEAN_LEVEL, "MEAN_LEVEL request");

    auto eccentricityOnly = AttributeComputation::computeAttributes(
        weighted,
        std::vector<AttributeOrGroup>{ECCENTRICITY});
    requireEqual(eccentricityOnly.first.NUM_ATTRIBUTES, 1, "ECCENTRICITY weighted request exposes one public attribute");
    require(eccentricityOnly.first.contains(ECCENTRICITY), "ECCENTRICITY weighted request exposes requested attribute");
    requireHiddenDependenciesAreAbsent(
        eccentricityOnly.first,
        {AREA, CENTRAL_MOMENT_20, CENTRAL_MOMENT_02, CENTRAL_MOMENT_11},
        "ECCENTRICITY weighted request");
    requireFiniteLiveNodeValues(weighted.topology(), eccentricityOnly.first, eccentricityOnly.second, ECCENTRICITY, "ECCENTRICITY weighted request");

    auto topologyEccentricityOnly = AttributeComputation::computeTopologyAttributes(
        weighted.topology(),
        std::vector<AttributeOrGroup>{ECCENTRICITY});
    requireEqual(topologyEccentricityOnly.first.NUM_ATTRIBUTES, 1, "ECCENTRICITY topology request exposes one public attribute");
    require(topologyEccentricityOnly.first.contains(ECCENTRICITY), "ECCENTRICITY topology request exposes requested attribute");
    requireHiddenDependenciesAreAbsent(
        topologyEccentricityOnly.first,
        {AREA, CENTRAL_MOMENT_20, CENTRAL_MOMENT_02, CENTRAL_MOMENT_11},
        "ECCENTRICITY topology request");
    requireFiniteLiveNodeValues(
        weighted.topology(),
        topologyEccentricityOnly.first,
        topologyEccentricityOnly.second,
        ECCENTRICITY,
        "ECCENTRICITY topology request");

    auto huOnly = AttributeComputation::computeAttributes(
        weighted,
        std::vector<AttributeOrGroup>{HU_MOMENT_1});
    requireEqual(huOnly.first.NUM_ATTRIBUTES, 1, "HU_MOMENT_1 request exposes one public attribute");
    require(huOnly.first.contains(HU_MOMENT_1), "HU_MOMENT_1 request exposes requested attribute");
    requireHiddenDependenciesAreAbsent(
        huOnly.first,
        {AREA, CENTRAL_MOMENT_20, CENTRAL_MOMENT_02, CENTRAL_MOMENT_11, CENTRAL_MOMENT_30, CENTRAL_MOMENT_03, CENTRAL_MOMENT_21, CENTRAL_MOMENT_12},
        "HU_MOMENT_1 request");
    requireFiniteLiveNodeValues(weighted.topology(), huOnly.first, huOnly.second, HU_MOMENT_1, "HU_MOMENT_1 request");

    auto rectangularityOnly = AttributeComputation::computeTopologyAttributes(
        weighted.topology(),
        std::vector<AttributeOrGroup>{RECTANGULARITY});
    requireEqual(rectangularityOnly.first.NUM_ATTRIBUTES, 1, "RECTANGULARITY topology request exposes one public attribute");
    require(rectangularityOnly.first.contains(RECTANGULARITY), "RECTANGULARITY topology request exposes requested attribute");
    requireHiddenDependenciesAreAbsent(rectangularityOnly.first, {AREA}, "RECTANGULARITY topology request");
    requireFiniteLiveNodeValues(
        weighted.topology(),
        rectangularityOnly.first,
        rectangularityOnly.second,
        RECTANGULARITY,
        "RECTANGULARITY topology request");

    auto grayGroup = AttributeComputation::computeAttributes(
        weighted,
        std::vector<AttributeOrGroup>{AttributeGroup::GRAY_LEVEL});
    requireEqual(grayGroup.first.NUM_ATTRIBUTES, 6, "GRAY_LEVEL group public attribute stride");
    require(!grayGroup.first.contains(AREA), "GRAY_LEVEL group must not expose hidden AREA dependency");

    return 0;
}
