#include <mmcfilters/attributes/Attributes.hpp>
#include <mmcfilters/trees/MorphologicalTreeFactory.hpp>
#include <mmcfilters/utils/Image.hpp>

#include <array>
#include <cstdint>
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

    return 0;
}
