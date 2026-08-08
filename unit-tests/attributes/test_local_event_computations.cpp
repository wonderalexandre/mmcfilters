#include "support/TestSupport.hpp"

#include "mmcfilters/attributes/AttributeComputation.hpp"
#include "mmcfilters/attributes/computers/AttributeComputerRegistry.hpp"
#include "mmcfilters/attributes/computers/BitquadAttributeComputer.hpp"
#include "mmcfilters/attributes/computers/ContourSideAttributeComputer.hpp"
#include "mmcfilters/attributes/computers/detail/BitquadAttributeMaterialization.hpp"
#include "mmcfilters/attributes/computers/detail/BitquadLocalEventComputation.hpp"
#include "mmcfilters/attributes/computers/detail/ContourSideAttributeMaterialization.hpp"
#include "mmcfilters/attributes/computers/detail/ContourSideLocalEventComputation.hpp"
#include "mmcfilters/localEvents/EventEngine.hpp"
#include "mmcfilters/trees/detail/ProperPartEntryNode.hpp"
#include "mmcfilters/trees/WeightedTreeView.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <span>
#include <unordered_map>

using namespace mmcfilters;
using namespace mmcfilters::attributes::computers;
using namespace mmcfilters::attributes::computers::detail;
using namespace mmcfilters::local_events;
using namespace mmcfilters::unit_tests;

namespace {

using Offset = mmcfilters::local_events::WindowOffset;
using BitquadFamilyCounts = mmcfilters::attributes::computers::detail::BitquadFamilyCounts;
using ContourSideCounts = mmcfilters::attributes::computers::detail::ContourSideCounts;

AttributeNames makeDenseAttributeNames(const std::vector<Attribute>& attributes) {
    std::unordered_map<Attribute, int> offsets;
    for (int i = 0; i < static_cast<int>(attributes.size()); ++i) {
        offsets[attributes[static_cast<std::size_t>(i)]] = i;
    }
    return AttributeNames(std::move(offsets));
}

std::vector<uint8_t> supportMask(const MorphologicalTree& tree, NodeId nodeId) {
    std::vector<uint8_t> mask(static_cast<std::size_t>(tree.getNumTotalProperParts()), 0);
    for (NodeId subtreeNodeId : tree.getNodeSubtree(nodeId)) {
        for (int pixel : tree.getProperParts(subtreeNodeId)) {
            mask[static_cast<std::size_t>(pixel)] = 1;
        }
    }
    return mask;
}

ImageUInt8Ptr makeRampImage(int rows, int cols) {
    auto image = ImageUInt8::create(rows, cols);
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            (*image)[row * cols + col] = static_cast<uint8_t>((row * 7 + col * 11) & 0xff);
        }
    }
    return image;
}

std::vector<ContourSideCounts> expectedContourSideCounts(const MorphologicalTree& tree) {
    const int rows = tree.getNumRowsOfGridDomain2D();
    const int cols = tree.getNumColsOfGridDomain2D();

    std::vector<ContourSideCounts> expected(static_cast<std::size_t>(tree.getNumInternalNodeSlots()));
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        const std::vector<uint8_t> mask = supportMask(tree, nodeId);
        ContourSideCounts counts;
        for (int p = 0; p < rows * cols; ++p) {
            if (!mask[static_cast<std::size_t>(p)]) {
                continue;
            }

            const auto [row, col] = ImageUtils::to2D(p, cols);
            auto sideIsExposed = [&](int qRow, int qCol) {
                if (qRow < 0 || qRow >= rows || qCol < 0 || qCol >= cols) {
                    return true;
                }
                const int q = ImageUtils::to1D(qRow, qCol, cols);
                return !mask[static_cast<std::size_t>(q)];
            };

            const int north = sideIsExposed(row - 1, col) ? 1 : 0;
            const int west = sideIsExposed(row, col - 1) ? 1 : 0;
            const int east = sideIsExposed(row, col + 1) ? 1 : 0;
            const int south = sideIsExposed(row + 1, col) ? 1 : 0;
            const int exposedSides = north + west + east + south;
            if (exposedSides > 0) {
                counts.contourPixels += 1;
            }
            counts.exposedSides += exposedSides;
            counts.north += north;
            counts.west += west;
            counts.east += east;
            counts.south += south;
        }
        expected[static_cast<std::size_t>(nodeId)] = counts;
    }
    return expected;
}

std::vector<int> expectedContourCounts(const MorphologicalTree& tree) {
    return ContourSideLocalEventComputation::projectContourPixels(expectedContourSideCounts(tree));
}

std::vector<std::array<int, 16>> expectedBitquadHistograms(const MorphologicalTree& tree) {
    const int rows = tree.getNumRowsOfGridDomain2D();
    const int cols = tree.getNumColsOfGridDomain2D();
    const std::array<Offset, 4> cellOffsets = {{
        {0, 0},
        {1, 0},
        {0, 1},
        {1, 1},
    }};

    std::vector<std::array<int, 16>> expected(static_cast<std::size_t>(tree.getNumInternalNodeSlots()));
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        const std::vector<uint8_t> mask = supportMask(tree, nodeId);
        auto& hist = expected[static_cast<std::size_t>(nodeId)];
        for (int row = -1; row < rows; ++row) {
            for (int col = -1; col < cols; ++col) {
                uint32_t code = 0;
                for (std::size_t bit = 0; bit < cellOffsets.size(); ++bit) {
                    const int qRow = row + cellOffsets[bit].rowOffset;
                    const int qCol = col + cellOffsets[bit].colOffset;
                    if (qRow < 0 || qRow >= rows || qCol < 0 || qCol >= cols) {
                        continue;
                    }
                    const int q = ImageUtils::to1D(qRow, qCol, cols);
                    if (mask[static_cast<std::size_t>(q)]) {
                        code |= uint32_t{1} << bit;
                    }
                }
                hist[static_cast<std::size_t>(code)] += 1;
            }
        }
    }
    return expected;
}

void requireArrayEqual(const std::array<int, 16>& actual, const std::array<int, 16>& expected, const std::string& label) {
    requireVectorEqual(std::vector<int>(actual.begin(), actual.end()), std::vector<int>(expected.begin(), expected.end()), label);
}

void requireFamily(BitquadLocalEventComputation::StateFamily actual, BitquadLocalEventComputation::StateFamily expected, const std::string& label) {
    require(actual == expected, label);
}

void requireBitquadFamilyCountsEqual(const BitquadFamilyCounts& actual, const BitquadFamilyCounts& expected, const std::string& label) {
    requireEqual(actual.empty, expected.empty, label + " empty");
    requireEqual(actual.q1, expected.q1, label + " q1");
    requireEqual(actual.q2, expected.q2, label + " q2");
    requireEqual(actual.qd, expected.qd, label + " qd");
    requireEqual(actual.q3, expected.q3, label + " q3");
    requireEqual(actual.q4, expected.q4, label + " q4");
}

void requireFloatEquivalent(float actual, float expected, const std::string& label) {
    if (std::isnan(expected)) {
        require(std::isnan(actual), label + ": expected NaN");
        return;
    }
    if (std::isinf(expected)) {
        require(std::isinf(actual) && std::signbit(actual) == std::signbit(expected), label + ": expected infinity");
        return;
    }
    requireNear(actual, expected, 1.0e-5f, label);
}

float projectedContourSideScalarValue(const ContourSideCounts& counts, Attribute attribute) {
    switch (attribute) {
    case CONTOUR_PIXELS:
        return static_cast<float>(counts.contourPixels);
    case CONTOUR_PERIMETER:
        return static_cast<float>(counts.exposedSides);
    case CONTOUR_SIDE_NORTH:
        return static_cast<float>(counts.north);
    case CONTOUR_SIDE_WEST:
        return static_cast<float>(counts.west);
    case CONTOUR_SIDE_EAST:
        return static_cast<float>(counts.east);
    case CONTOUR_SIDE_SOUTH:
        return static_cast<float>(counts.south);
    default:
        throw std::runtime_error("Unsupported projected contour-side scalar attribute.");
    }
}

std::vector<BitquadFamilyCounts> makeQ1ConnectivityProbeCounts(const MorphologicalTree& tree) {
    std::vector<BitquadFamilyCounts> counts(static_cast<std::size_t>(tree.getNumInternalNodeSlots()));
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        counts[static_cast<std::size_t>(nodeId)].q1 = 4;
    }
    return counts;
}

std::vector<float> materializeTreeOfShapesQ1ProbeEuler(const WeightedMorphologicalTree<std::uint8_t>& weighted) {
    const MorphologicalTree& tree = weighted.topology();
    const std::vector<Attribute> requestedAttributes = {BITQUADS_NUMBER_EULER};
    const AttributeNames names = makeDenseAttributeNames(requestedAttributes);
    std::vector<float> buffer(static_cast<std::size_t>(tree.getNumInternalNodeSlots()) * static_cast<std::size_t>(names.NUM_ATTRIBUTES), 0.0f);

    const auto counts = makeQ1ConnectivityProbeCounts(tree);
    BitquadAttributeMaterialization::materializeAttributesFromBitquadFamilyCounts(tree, weighted.altitudeSpan(), std::span<const BitquadFamilyCounts>(counts),
                                                                                  std::span<float>(buffer), names, requestedAttributes);
    return buffer;
}

template <class T> std::vector<T> copyAltitudeAs(const WeightedMorphologicalTree<std::uint8_t>& weighted) {
    const AltitudeBuffer<std::uint8_t>& altitude = weighted.getAltitudeBuffer();
    std::vector<T> converted;
    converted.reserve(altitude.size());
    for (std::uint8_t value : altitude) {
        converted.push_back(static_cast<T>(value));
    }
    return converted;
}

void verifyTreeOfShapesScalarConnectivityPolicy(const ImageUInt8Ptr& image, ToSInterpolation interpolation, const std::string& label,
                                                bool requireBothPolarities) {
    const auto weighted = makeWeightedTreeOfShapes(image, interpolation);
    const MorphologicalTree& tree = weighted->topology();
    const AttributeNames names = makeDenseAttributeNames({BITQUADS_NUMBER_EULER});
    const auto buffer = materializeTreeOfShapesQ1ProbeEuler(*weighted);

    const bool minRootIs4Connectivity = tree.getDecreasingGridAdjacency2D()->is4connectivity();
    const bool maxRootIs4Connectivity = tree.getIncreasingGridAdjacency2D()->is4connectivity();
    const bool rootUses4Connectivity = minRootIs4Connectivity == maxRootIs4Connectivity ? minRootIs4Connectivity : false;
    requireFloatEquivalent(buffer[static_cast<std::size_t>(names.linearIndex(tree.getRoot(), BITQUADS_NUMBER_EULER))], rootUses4Connectivity ? 1.0f : 0.0f,
                           label + " ToS root scalar projection connectivity");

    bool sawMinTreeNode = false;
    bool sawMaxTreeNode = false;
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        if (tree.isRoot(nodeId)) {
            continue;
        }

        const NodeId parentNodeId = tree.getNodeParent(nodeId);
        const std::uint8_t nodeAltitude = weighted->getAltitude(nodeId);
        const std::uint8_t parentAltitude = weighted->getAltitude(parentNodeId);
        if (nodeAltitude == parentAltitude) {
            continue;
        }

        const bool isMaxTreeNode = nodeAltitude > parentAltitude;
        const bool expected4Connectivity =
            isMaxTreeNode ? tree.getIncreasingGridAdjacency2D()->is4connectivity() : tree.getDecreasingGridAdjacency2D()->is4connectivity();
        sawMaxTreeNode = sawMaxTreeNode || isMaxTreeNode;
        sawMinTreeNode = sawMinTreeNode || !isMaxTreeNode;

        requireFloatEquivalent(buffer[static_cast<std::size_t>(names.linearIndex(nodeId, BITQUADS_NUMBER_EULER))], expected4Connectivity ? 1.0f : 0.0f,
                               label + " ToS scalar projection must choose min/max adjacency by node polarity " + std::to_string(nodeId));
    }

    if (requireBothPolarities) {
        require(sawMinTreeNode, label + " ToS scalar connectivity policy fixture must contain a min-tree node");
        require(sawMaxTreeNode, label + " ToS scalar connectivity policy fixture must contain a max-tree node");
    }
}

void verifyTreeOfShapesScalarConnectivityPolicy(const ImageUInt8Ptr& image) {
    verifyTreeOfShapesScalarConnectivityPolicy(image, ToSInterpolation::SelfDual, "SelfDual", false);
    verifyTreeOfShapesScalarConnectivityPolicy(image, ToSInterpolation::Min4cMax8c, "Min4cMax8c", true);
    verifyTreeOfShapesScalarConnectivityPolicy(image, ToSInterpolation::Min8cMax4c, "Min8cMax4c", true);
}

void verifyBitquadStateFamilyTable() {
    using Family = BitquadLocalEventComputation::StateFamily;
    const std::array<Family, 16> expected = {{
        Family::Empty,
        Family::Q1,
        Family::Q1,
        Family::Q2,
        Family::Q1,
        Family::Q2,
        Family::QD,
        Family::Q3,
        Family::Q1,
        Family::QD,
        Family::Q2,
        Family::Q3,
        Family::Q2,
        Family::Q3,
        Family::Q3,
        Family::Q4,
    }};

    const auto actual = BitquadLocalEventComputation::stateFamilyTable();
    for (std::size_t state = 0; state < expected.size(); ++state) {
        requireFamily(actual[state], expected[state], "bitquad 2x2 state family " + std::to_string(state));

        std::array<int, 16> oneHot{};
        oneHot[state] = 1;
        const auto counts = BitquadLocalEventComputation::projectBitquadFamilyCounts(oneHot);
        BitquadFamilyCounts expectedCounts;
        switch (expected[state]) {
        case Family::Empty:
            expectedCounts.empty = 1;
            break;
        case Family::Q1:
            expectedCounts.q1 = 1;
            break;
        case Family::Q2:
            expectedCounts.q2 = 1;
            break;
        case Family::QD:
            expectedCounts.qd = 1;
            break;
        case Family::Q3:
            expectedCounts.q3 = 1;
            break;
        case Family::Q4:
            expectedCounts.q4 = 1;
            break;
        }
        requireBitquadFamilyCountsEqual(counts, expectedCounts, "bitquad one-hot family projection " + std::to_string(state));
    }

    requireThrows<std::out_of_range>([]() { static_cast<void>(BitquadLocalEventComputation::stateFamily(16)); },
                                     "bitquad state family must reject states outside 0..15");
}

void verifyLocalEventBitquadScalarComputer(const MorphologicalTree& tree, const std::vector<Attribute>& requestedAttributes, const std::string& label) {
    const auto scalarAttributes = runtimeProducedAttributes<BitquadAttributeComputer>();
    requireEqual(scalarAttributes.size(), static_cast<std::size_t>(9), label + " attribute count");

    const AttributeNames names = makeDenseAttributeNames(requestedAttributes);
    const std::size_t bufferSize = static_cast<std::size_t>(tree.getNumInternalNodeSlots()) * static_cast<std::size_t>(names.NUM_ATTRIBUTES);
    std::vector<float> localBuffer(bufferSize, 0.0f);
    std::vector<float> directFamilyBuffer(bufferSize, 0.0f);

    const auto directBitquadFamilyDeltas = BitquadLocalEventComputation::computeBitquadFamilyDeltas(tree);
    const auto directBitquadFamilyCounts = BitquadLocalEventComputation::aggregateBitquadFamilyDeltas(tree, directBitquadFamilyDeltas);

    BitquadAttributeComputer::compute(
        AttributeComputeContext<float>{tree, std::span<float>(localBuffer), names, std::span<const Attribute>(requestedAttributes)});
    BitquadAttributeMaterialization::materializeAttributesFromBitquadFamilyCounts(tree, std::span<const BitquadFamilyCounts>(directBitquadFamilyCounts),
                                                                                  std::span<float>(directFamilyBuffer), names, requestedAttributes);

    for (NodeId nodeId : tree.getAliveNodeIds()) {
        for (Attribute attribute : requestedAttributes) {
            const int index = names.linearIndex(nodeId, attribute);
            requireFloatEquivalent(localBuffer[static_cast<std::size_t>(index)], directFamilyBuffer[static_cast<std::size_t>(index)],
                                   label + " local compute uses direct family scalar " + AttributeNames::toString(attribute) + " node " +
                                       std::to_string(nodeId));
        }
    }
}

void verifyLocalEventBitquadScalarComputer(const MorphologicalTree& tree, const std::string& label) {
    verifyLocalEventBitquadScalarComputer(tree, runtimeProducedAttributes<BitquadAttributeComputer>(), label + " full");
    verifyLocalEventBitquadScalarComputer(tree, {BITQUADS_AREA, BITQUADS_PERIMETER, BITQUADS_CIRCULARITY}, label + " subset");
}

void verifyLocalEventTreeOfShapesBitquadScalarComputer(const WeightedMorphologicalTree<std::uint8_t>& weighted, const std::string& label) {
    const MorphologicalTree& tree = weighted.topology();
    const auto requestedAttributes = runtimeProducedAttributes<BitquadAttributeComputer>();
    const AttributeNames names = makeDenseAttributeNames(requestedAttributes);
    const std::size_t bufferSize = static_cast<std::size_t>(tree.getNumInternalNodeSlots()) * static_cast<std::size_t>(names.NUM_ATTRIBUTES);
    std::vector<float> computeBuffer(bufferSize, 0.0f);
    std::vector<float> directBuffer(bufferSize, 0.0f);
    std::vector<float> genericSpanBuffer(bufferSize, 0.0f);
    std::vector<float> genericViewBuffer(bufferSize, 0.0f);

    BitquadAttributeComputer::compute(AltitudeAttributeComputeContext<float, std::uint8_t>{tree, weighted.altitudeSpan(), std::span<float>(computeBuffer),
                                                                                           names, std::span<const Attribute>(requestedAttributes)});

    const auto directBitquadFamilyCounts = BitquadLocalEventComputation::computeBitquadFamilyCounts(tree);
    BitquadAttributeMaterialization::materializeAttributesFromBitquadFamilyCounts(tree, weighted.altitudeSpan(),
                                                                                  std::span<const BitquadFamilyCounts>(directBitquadFamilyCounts),
                                                                                  std::span<float>(directBuffer), names, requestedAttributes);

    const std::vector<std::int32_t> externalIntAltitude = copyAltitudeAs<std::int32_t>(weighted);
    BitquadAttributeMaterialization::materializeAttributesFromBitquadFamilyCounts(tree, std::span<const std::int32_t>(externalIntAltitude),
                                                                                  std::span<const BitquadFamilyCounts>(directBitquadFamilyCounts),
                                                                                  std::span<float>(genericSpanBuffer), names, requestedAttributes);

    const std::vector<float> externalFloatAltitude = copyAltitudeAs<float>(weighted);
    const WeightedTreeView<float> externalFloatView(tree, std::span<const float>(externalFloatAltitude));
    BitquadAttributeMaterialization::materializeAttributesFromBitquadFamilyCounts(
        externalFloatView, std::span<const BitquadFamilyCounts>(directBitquadFamilyCounts), std::span<float>(genericViewBuffer), names, requestedAttributes);

    for (NodeId nodeId : tree.getAliveNodeIds()) {
        for (Attribute attribute : requestedAttributes) {
            const int index = names.linearIndex(nodeId, attribute);
            requireFloatEquivalent(computeBuffer[static_cast<std::size_t>(index)], directBuffer[static_cast<std::size_t>(index)],
                                   label + " ToS direct scalar " + AttributeNames::toString(attribute) + " node " + std::to_string(nodeId));
            requireFloatEquivalent(genericSpanBuffer[static_cast<std::size_t>(index)], directBuffer[static_cast<std::size_t>(index)],
                                   label + " ToS generic int32 span scalar " + AttributeNames::toString(attribute) + " node " + std::to_string(nodeId));
            requireFloatEquivalent(genericViewBuffer[static_cast<std::size_t>(index)], directBuffer[static_cast<std::size_t>(index)],
                                   label + " ToS generic float view scalar " + AttributeNames::toString(attribute) + " node " + std::to_string(nodeId));
        }
    }

    auto [pipelineNames, pipelineBuffer] = AttributeComputation::computeSingleAttribute(weighted, AttributeGroup::BOUNDARY);
    auto [weightedTopologyNames, weightedTopologyBuffer] = AttributeComputation::computeTopologyAttributes(weighted, {AttributeGroup::BOUNDARY});
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        for (Attribute attribute : requestedAttributes) {
            const int directIndex = names.linearIndex(nodeId, attribute);
            requireFloatEquivalent(pipelineBuffer[static_cast<std::size_t>(pipelineNames.linearIndex(nodeId, attribute))],
                                   directBuffer[static_cast<std::size_t>(directIndex)],
                                   label + " public weighted pipeline " + AttributeNames::toString(attribute) + " node " + std::to_string(nodeId));
            requireFloatEquivalent(weightedTopologyBuffer[static_cast<std::size_t>(weightedTopologyNames.linearIndex(nodeId, attribute))],
                                   directBuffer[static_cast<std::size_t>(directIndex)],
                                   label + " public weighted topology pipeline " + AttributeNames::toString(attribute) + " node " + std::to_string(nodeId));
        }
    }

    requireThrows<std::invalid_argument>([&]() { static_cast<void>(AttributeComputation::computeTopologyAttributes(tree, {AttributeGroup::BOUNDARY})); },
                                         label + " topology-only ToS boundary request must require altitude");

    requireThrows<std::invalid_argument>(
        [&]() {
            BitquadAttributeComputer::compute(
                AttributeComputeContext<float>{tree, std::span<float>(computeBuffer), names, std::span<const Attribute>(requestedAttributes)});
        },
        label + " ToS scalar projection without altitude must throw");

    std::vector<float> shortAltitude = externalFloatAltitude;
    if (!shortAltitude.empty()) {
        shortAltitude.pop_back();
    }
    requireThrows<std::runtime_error>(
        [&]() {
            std::vector<float> throwBuffer(bufferSize, 0.0f);
            BitquadAttributeMaterialization::materializeAttributesFromBitquadFamilyCounts(tree, std::span<const float>(shortAltitude),
                                                                                          std::span<const BitquadFamilyCounts>(directBitquadFamilyCounts),
                                                                                          std::span<float>(throwBuffer), names, requestedAttributes);
        },
        label + " ToS generic scalar projection with wrong altitude size must throw");

    NodeId nonRootNode = InvalidNode;
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        if (!tree.isRoot(nodeId)) {
            nonRootNode = nodeId;
            break;
        }
    }
    if (nonRootNode != InvalidNode) {
        std::vector<float> ambiguousAltitude = externalFloatAltitude;
        ambiguousAltitude[static_cast<std::size_t>(nonRootNode)] = ambiguousAltitude[static_cast<std::size_t>(tree.getNodeParent(nonRootNode))];
        requireThrows<std::runtime_error>(
            [&]() {
                std::vector<float> ambiguousBuffer(bufferSize, 0.0f);
                BitquadAttributeMaterialization::materializeAttributesFromBitquadFamilyCounts(tree, std::span<const float>(ambiguousAltitude),
                                                                                              std::span<const BitquadFamilyCounts>(directBitquadFamilyCounts),
                                                                                              std::span<float>(ambiguousBuffer), names, requestedAttributes);
            },
            label + " dual adjacency must reject an altitude-equal ambiguous "
                    "branch");
    }
}

void verifyGenericContourCounts(const MorphologicalTree& tree, const std::string& label) {
    const auto actualSideCounts = ContourSideLocalEventComputation::computeContourSideCounts(tree);
    const auto expectedSideCounts = expectedContourSideCounts(tree);
    const auto actualProjected = ContourSideLocalEventComputation::projectContourPixels(actualSideCounts);
    const auto actualExposedSides = ContourSideLocalEventComputation::projectExposedSides(actualSideCounts);
    const auto expected = expectedContourCounts(tree);
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        const std::size_t index = static_cast<std::size_t>(nodeId);
        requireEqual(actualProjected[index], expected[index], label + " contour count node " + std::to_string(nodeId));
        requireEqual(actualSideCounts[index].contourPixels, expectedSideCounts[index].contourPixels,
                     label + " side contour pixels node " + std::to_string(nodeId));
        requireEqual(actualSideCounts[index].exposedSides, expectedSideCounts[index].exposedSides, label + " exposed sides node " + std::to_string(nodeId));
        requireEqual(actualSideCounts[index].north, expectedSideCounts[index].north, label + " north exposed sides node " + std::to_string(nodeId));
        requireEqual(actualSideCounts[index].west, expectedSideCounts[index].west, label + " west exposed sides node " + std::to_string(nodeId));
        requireEqual(actualSideCounts[index].east, expectedSideCounts[index].east, label + " east exposed sides node " + std::to_string(nodeId));
        requireEqual(actualSideCounts[index].south, expectedSideCounts[index].south, label + " south exposed sides node " + std::to_string(nodeId));
        requireEqual(actualSideCounts[index].exposedSides,
                     actualSideCounts[index].north + actualSideCounts[index].west + actualSideCounts[index].east + actualSideCounts[index].south,
                     label + " exposed sides must equal directional sum node " + std::to_string(nodeId));
        requireEqual(actualExposedSides[index], actualSideCounts[index].exposedSides, label + " exposed-side projection node " + std::to_string(nodeId));
    }
}

void verifyPublicContourSideCountsComputer(const MorphologicalTree& tree, const std::vector<Attribute>& requestedAttributes, const std::string& label) {
    const auto allContourAttributes = runtimeProducedAttributes<ContourSideAttributeComputer>();
    const std::vector<Attribute> expectedAttributes = {
        CONTOUR_PIXELS, CONTOUR_PERIMETER, CONTOUR_SIDE_NORTH, CONTOUR_SIDE_WEST, CONTOUR_SIDE_EAST, CONTOUR_SIDE_SOUTH,
    };
    requireEqual(allContourAttributes.size(), expectedAttributes.size(), label + " contour attribute count");
    for (std::size_t i = 0; i < expectedAttributes.size(); ++i) {
        requireEqual(static_cast<int>(allContourAttributes[i]), static_cast<int>(expectedAttributes[i]),
                     label + " contour attribute order " + std::to_string(i));
    }

    const AttributeNames names = makeDenseAttributeNames(requestedAttributes);
    const std::size_t bufferSize = static_cast<std::size_t>(tree.getNumInternalNodeSlots()) * static_cast<std::size_t>(names.NUM_ATTRIBUTES);
    std::vector<float> computeBuffer(bufferSize, 0.0f);
    std::vector<float> directBuffer(bufferSize, 0.0f);

    ContourSideAttributeComputer::compute(
        AttributeComputeContext<float>{tree, std::span<float>(computeBuffer), names, std::span<const Attribute>(requestedAttributes)});

    const auto expectedSideCounts = expectedContourSideCounts(tree);
    ContourSideAttributeMaterialization::materializeAttributesFromContourSideCounts(tree, std::span<const ContourSideCounts>(expectedSideCounts),
                                                                                    std::span<float>(directBuffer), names, requestedAttributes);

    std::vector<AttributeOrGroup> publicRequests;
    publicRequests.reserve(requestedAttributes.size());
    for (Attribute attribute : requestedAttributes) {
        publicRequests.emplace_back(attribute);
    }
    auto publicComputed = AttributeComputation::computeTopologyAttributes(tree, publicRequests);

    for (NodeId nodeId : tree.getAliveNodeIds()) {
        const auto& expectedCounts = expectedSideCounts[static_cast<std::size_t>(nodeId)];
        for (Attribute attribute : requestedAttributes) {
            const int index = names.linearIndex(nodeId, attribute);
            requireFloatEquivalent(computeBuffer[static_cast<std::size_t>(index)], projectedContourSideScalarValue(expectedCounts, attribute),
                                   label + " scalar " + AttributeNames::toString(attribute) + " node " + std::to_string(nodeId));
            requireFloatEquivalent(computeBuffer[static_cast<std::size_t>(index)], directBuffer[static_cast<std::size_t>(index)],
                                   label + " direct contour projection " + AttributeNames::toString(attribute) + " node " + std::to_string(nodeId));
            requireFloatEquivalent(publicComputed.values()[static_cast<std::size_t>(publicComputed.attributeNames().linearIndex(nodeId, attribute))],
                                   computeBuffer[static_cast<std::size_t>(index)],
                                   label + " public contour request " + AttributeNames::toString(attribute) + " node " + std::to_string(nodeId));
        }
    }
}

void verifyPublicContourSideCountsComputer(const MorphologicalTree& tree, const std::string& label) {
    verifyPublicContourSideCountsComputer(tree, runtimeProducedAttributes<ContourSideAttributeComputer>(), label + " full");
    verifyPublicContourSideCountsComputer(tree, {CONTOUR_PIXELS, CONTOUR_PERIMETER, CONTOUR_SIDE_EAST}, label + " subset");
}

void verifyPublicContourSideCountUnitAttributes(const MorphologicalTree& tree, const std::string& label) {
    const auto attributes = runtimeProducedAttributes<ContourSideAttributeComputer>();
    const AttributeNames names = makeDenseAttributeNames(attributes);
    const std::vector<NodeId> unitProperParts = {
        0,
        static_cast<NodeId>(tree.getNumTotalProperParts() - 1),
    };
    std::vector<float> buffer(unitProperParts.size() * static_cast<std::size_t>(names.NUM_ATTRIBUTES), 0.0f);

    ContourSideAttributeComputer::computeUnitRows(UnitAttributeComputeContext<float>{tree, std::span<const NodeId>(unitProperParts), std::span<float>(buffer),
                                                                                     names, std::span<const Attribute>(attributes)});

    ContourSideCounts expected;
    expected.contourPixels = 1;
    expected.exposedSides = 4;
    expected.north = 1;
    expected.west = 1;
    expected.east = 1;
    expected.south = 1;
    for (NodeId leafIndex = 0; leafIndex < static_cast<NodeId>(unitProperParts.size()); ++leafIndex) {
        for (Attribute attribute : attributes) {
            const int index = names.linearIndex(leafIndex, attribute);
            requireFloatEquivalent(buffer[static_cast<std::size_t>(index)], projectedContourSideScalarValue(expected, attribute),
                                   label + " unit " + AttributeNames::toString(attribute) + " leaf " + std::to_string(leafIndex));
        }
    }
}

void verifyPublicContourAttributeNames() {
    requireEqual(AttributeNames::toString(CONTOUR_PIXELS), std::string("CONTOUR_PIXELS"), "contour pixels name");
    requireEqual(AttributeNames::toString(CONTOUR_PERIMETER), std::string("CONTOUR_PERIMETER"), "contour perimeter name");
    require(AttributeNames::describe(CONTOUR_SIDE_NORTH).starts_with("Contour north sides:"), "contour north description");
    auto parsed = AttributeNames::parse("CONTOUR_PIXELS");
    require(parsed.has_value(), "public contour attributes must be parsed");
    requireEqual(static_cast<int>(*parsed), static_cast<int>(CONTOUR_PIXELS), "public contour parse value");
}

void verifyBitquadDeltaProjection(const MorphologicalTree& tree, const std::string& label) {
    const auto stateDeltas = BitquadLocalEventComputation::computeBitquadStateHistogramDeltas(tree);
    const auto familyDeltas = BitquadLocalEventComputation::computeBitquadFamilyDeltas(tree);
    const auto familyDeltasFromStates = BitquadLocalEventComputation::computeBitquadFamilyCounts(stateDeltas);

    const auto aggregatedStates = BitquadLocalEventComputation::aggregateBitquadStateHistogramDeltas(tree, stateDeltas);
    const auto aggregatedFamilies = BitquadLocalEventComputation::aggregateBitquadFamilyDeltas(tree, familyDeltas);
    const auto directStates = BitquadLocalEventComputation::computeBitquadStateHistograms(tree);
    const auto directFamilies = BitquadLocalEventComputation::computeBitquadFamilyCounts(tree);

    for (NodeId nodeId : tree.getAliveNodeIds()) {
        const std::size_t nodeIndex = static_cast<std::size_t>(nodeId);
        requireArrayEqual(aggregatedStates[nodeIndex], directStates[nodeIndex],
                          label + " bitquad state deltas aggregate to direct histogram node " + std::to_string(nodeId));
        requireBitquadFamilyCountsEqual(aggregatedFamilies[nodeIndex], directFamilies[nodeIndex],
                                        label + " bitquad family deltas aggregate to direct counts node " + std::to_string(nodeId));
        requireBitquadFamilyCountsEqual(familyDeltasFromStates[nodeIndex], familyDeltas[nodeIndex],
                                        label + " bitquad family deltas match projected state deltas node " + std::to_string(nodeId));
    }

    const auto properPartStates = BitquadLocalEventComputation::projectBitquadStateHistogramsToProperParts(tree, directStates);
    const auto properPartFamilies = BitquadLocalEventComputation::projectBitquadFamilyCountsToProperParts(tree, directFamilies);
    const auto properPartFamilyDeltas = BitquadLocalEventComputation::projectBitquadFamilyCountsToProperParts(tree, familyDeltas);
    for (NodeId properPart = 0; properPart < tree.getNumTotalProperParts(); ++properPart) {
        const NodeId owner = tree.getProperPartOwner(properPart);
        const std::size_t properPartIndex = static_cast<std::size_t>(properPart);
        const std::size_t ownerIndex = static_cast<std::size_t>(owner);
        requireArrayEqual(properPartStates[properPartIndex], directStates[ownerIndex], label + " proper-part state projection " + std::to_string(properPart));
        requireBitquadFamilyCountsEqual(properPartFamilies[properPartIndex], directFamilies[ownerIndex],
                                        label + " proper-part family projection " + std::to_string(properPart));
        requireBitquadFamilyCountsEqual(properPartFamilyDeltas[properPartIndex], familyDeltas[ownerIndex],
                                        label + " proper-part family-delta projection " + std::to_string(properPart));
    }

    std::vector<BitquadFamilyCounts> shortFamilyDeltas = familyDeltas;
    if (!shortFamilyDeltas.empty()) {
        shortFamilyDeltas.pop_back();
    }
    requireThrows<std::invalid_argument>([&]() { static_cast<void>(BitquadLocalEventComputation::aggregateBitquadFamilyDeltas(tree, shortFamilyDeltas)); },
                                         label + " bitquad family delta aggregation must validate node-slot coverage");
}

void verifyGenericBitquadHistograms(const MorphologicalTree& tree, const std::string& label) {
    BitquadLocalEventComputation computer(tree);
    const auto& actual = computer.getBitquadStateHistograms();
    const auto& actualBitquadFamilyCounts = computer.getBitquadFamilyCounts();
    const auto directBitquadFamilyCounts = BitquadLocalEventComputation::computeBitquadFamilyCounts(tree);
    const auto expected = expectedBitquadHistograms(tree);
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        requireArrayEqual(actual[static_cast<std::size_t>(nodeId)], expected[static_cast<std::size_t>(nodeId)],
                          label + " bitquad histogram node " + std::to_string(nodeId));
        requireBitquadFamilyCountsEqual(actualBitquadFamilyCounts[static_cast<std::size_t>(nodeId)],
                                        BitquadLocalEventComputation::projectBitquadFamilyCounts(expected[static_cast<std::size_t>(nodeId)]),
                                        label + " bitquad family count node " + std::to_string(nodeId));
        requireBitquadFamilyCountsEqual(directBitquadFamilyCounts[static_cast<std::size_t>(nodeId)],
                                        actualBitquadFamilyCounts[static_cast<std::size_t>(nodeId)],
                                        label + " direct bitquad family count node " + std::to_string(nodeId));
    }
    verifyBitquadDeltaProjection(tree, label);
}

MorphologicalTree makeTwoBranchTreeOfShapes() {
    std::vector<NodeId> parent = {
        4, 4, 5, 5, 6, 6, 6,
    };
    std::vector<std::uint8_t> altitude(parent.size(), std::uint8_t{});
    auto weighted = MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(parent), std::span<const std::uint8_t>(altitude), 2, 2,
                                                                    MorphologicalTreeKind::TREE_OF_SHAPES);
    return weighted.topology().clone();
}

} // namespace

int main() {
    verifyBitquadStateFamilyTable();
    verifyPublicContourAttributeNames();

    auto image = makeComponentTreeFixture();
    auto rampImage = makeRampImage(16, 16);

    {
        auto tree = makeComponentTree(image, true);
        verifyGenericContourCounts(*tree, "max-tree");
        verifyPublicContourSideCountsComputer(*tree, "max-tree");
        verifyPublicContourSideCountUnitAttributes(*tree, "max-tree");
        verifyGenericBitquadHistograms(*tree, "max-tree");
        verifyLocalEventBitquadScalarComputer(*tree, "max-tree");
    }

    {
        auto tree = makeComponentTree(image, false);
        verifyGenericContourCounts(*tree, "min-tree");
        verifyPublicContourSideCountsComputer(*tree, "min-tree");
        verifyGenericBitquadHistograms(*tree, "min-tree");
        verifyLocalEventBitquadScalarComputer(*tree, "min-tree");
    }

    {
        auto tree = makeComponentTree(image, true, 1.0);
        verifyGenericContourCounts(*tree, "max-tree 4-connectivity");
        verifyPublicContourSideCountsComputer(*tree, "max-tree 4-connectivity");
        verifyGenericBitquadHistograms(*tree, "max-tree 4-connectivity");
        verifyLocalEventBitquadScalarComputer(*tree, "max-tree 4-connectivity");
    }

    {
        auto tree = makeComponentTree(image, false, 1.0);
        verifyGenericContourCounts(*tree, "min-tree 4-connectivity");
        verifyPublicContourSideCountsComputer(*tree, "min-tree 4-connectivity");
        verifyGenericBitquadHistograms(*tree, "min-tree 4-connectivity");
        verifyLocalEventBitquadScalarComputer(*tree, "min-tree 4-connectivity");
    }

    {
        auto tree = makeComponentTree(rampImage, true, 1.0);
        verifyGenericBitquadHistograms(*tree, "ramp max-tree 4-connectivity");
        verifyLocalEventBitquadScalarComputer(*tree, "ramp max-tree 4-connectivity");
    }

    {
        auto tree = makeComponentTree(rampImage, false, 1.0);
        verifyGenericBitquadHistograms(*tree, "ramp min-tree 4-connectivity");
        verifyLocalEventBitquadScalarComputer(*tree, "ramp min-tree 4-connectivity");
    }

    {
        auto tree = makeTreeOfShapes(image, ToSInterpolation::SelfDual);
        verifyGenericContourCounts(*tree, "tree of shapes");
        verifyPublicContourSideCountsComputer(*tree, "tree of shapes");
        verifyGenericBitquadHistograms(*tree, "tree of shapes");
    }

    {
        verifyTreeOfShapesScalarConnectivityPolicy(image);
        auto weightedMin4Max8 = makeWeightedTreeOfShapes(image, ToSInterpolation::Min4cMax8c);
        verifyGenericBitquadHistograms(weightedMin4Max8->topology(), "tree of shapes Min4cMax8c");
        verifyLocalEventTreeOfShapesBitquadScalarComputer(*weightedMin4Max8, "weighted tree of shapes Min4cMax8c");
        auto weightedMin8Max4 = makeWeightedTreeOfShapes(image, ToSInterpolation::Min8cMax4c);
        verifyGenericBitquadHistograms(weightedMin8Max4->topology(), "tree of shapes Min8cMax4c");
        verifyLocalEventTreeOfShapesBitquadScalarComputer(*weightedMin8Max4, "weighted tree of shapes Min8cMax4c");
    }

    {
        auto tree = makeTwoBranchTreeOfShapes();
        requireEqual(tree.getRoot(), 2, "two-branch ToS root slot");
        requireEqual(mmcfilters::detail::properPartEntryNode(tree, 0, 1), 0, "same-branch sample enters at anchor owner");
        requireEqual(mmcfilters::detail::properPartEntryNode(tree, 0, 2), 2, "cross-branch sample enters at LCA");
        requireEqual(mmcfilters::local_events::EventEngine::entryNode(tree, 0, Offset{-1, 0}), InvalidNode, "out-of-domain sample has no entry");

        verifyGenericContourCounts(tree, "two-branch tree of shapes");
        verifyPublicContourSideCountsComputer(tree, "two-branch tree of shapes");
        verifyGenericBitquadHistograms(tree, "two-branch tree of shapes");
    }

    return 0;
}
