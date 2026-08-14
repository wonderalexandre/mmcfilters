#include "support/TestSupport.hpp"

#include "mmcfilters/attributes/AttributeComputation.hpp"
#include "mmcfilters/attributes/computers/AttributeComputerRegistry.hpp"
#include "mmcfilters/attributes/computers/BitquadAttributeComputer.hpp"
#include "mmcfilters/attributes/computers/ContourSideAttributeComputer.hpp"
#include "mmcfilters/attributes/computers/detail/BitquadAttributeProjection.hpp"
#include "mmcfilters/attributes/computers/detail/BitquadFiniteWindowComputation.hpp"
#include "mmcfilters/attributes/computers/detail/ContourSideAttributeMaterialization.hpp"
#include "mmcfilters/attributes/computers/detail/ContourSideFiniteWindowComputation.hpp"
#include "mmcfilters/localAttributes/FiniteWindowLocalAttributeComputer.hpp"
#include "mmcfilters/trees/detail/MorphologicalTreeConstructionContextQueries.hpp"
#include "mmcfilters/trees/ValuedMorphologicalTreeView.hpp"
#include "mmcfilters/utils/Contract.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <optional>
#include <random>
#include <span>
#include <type_traits>
#include <unordered_map>
#include <vector>

using namespace mmcfilters;
using namespace mmcfilters::attributes::computers;
using namespace mmcfilters::attributes::computers::detail;
using namespace mmcfilters::local_attributes;
using namespace mmcfilters::unit_tests;

namespace {

using Offset = mmcfilters::local_attributes::WindowOffset;
using BitquadCode = mmcfilters::attributes::computers::detail::BitquadCode;
using BitquadFamily = mmcfilters::attributes::computers::detail::BitquadFamily;
using BitquadFamilyIncrement = mmcfilters::attributes::computers::detail::BitquadFamilyIncrement;
using BitquadFamilyCounts = mmcfilters::attributes::computers::detail::BitquadFamilyCounts;
using BitquadState = mmcfilters::attributes::computers::detail::BitquadState;
using NonemptyBitquadStateHistogramIncrement = mmcfilters::attributes::computers::detail::NonemptyBitquadStateHistogramIncrement;
using NonemptyBitquadStateHistogram = mmcfilters::attributes::computers::detail::NonemptyBitquadStateHistogram;
using BitquadStateHistogram = mmcfilters::attributes::computers::detail::BitquadStateHistogram;
using ContourSideCounts = mmcfilters::attributes::computers::detail::ContourSideCounts;

static_assert(!std::is_same_v<BitquadFamilyIncrement, BitquadFamilyCounts>);
static_assert(!std::is_same_v<NonemptyBitquadStateHistogramIncrement, NonemptyBitquadStateHistogram>);
static_assert(!std::is_same_v<NonemptyBitquadStateHistogram, BitquadStateHistogram>);

AttributeNames makeDenseAttributeNames(const std::vector<Attribute>& attributes) {
    std::unordered_map<Attribute, int> offsets;
    for (int i = 0; i < static_cast<int>(attributes.size()); ++i) {
        offsets[attributes[static_cast<std::size_t>(i)]] = i;
    }
    return AttributeNames(std::move(offsets));
}

std::vector<uint8_t> supportMask(const MorphologicalTree& tree, NodeId nodeId) {
    std::vector<uint8_t> mask(static_cast<std::size_t>(tree.numPixels()), 0);
    for (NodeId subtreeNodeId : tree.subtreeNodes(nodeId)) {
        for (PixelId pixel : tree.properPart(subtreeNodeId)) {
            mask[static_cast<std::size_t>(pixel)] = 1;
        }
    }
    return mask;
}

ImageUInt8Ptr makeRampImage(int rows, int columns) {
    auto image = ImageUInt8::create(rows, columns);
    for (int row = 0; row < rows; ++row) {
        for (int column = 0; column < columns; ++column) {
            (*image)[row * columns + column] = static_cast<uint8_t>((row * 7 + column * 11) & 0xff);
        }
    }
    return image;
}

std::vector<ContourSideCounts> expectedContourSideCounts(const MorphologicalTree& tree) {
    const int rows = tree.numRows();
    const int columns = tree.numColumns();

    std::vector<ContourSideCounts> expected(static_cast<std::size_t>(tree.numInternalNodeSlots()));
    for (NodeId nodeId : tree.aliveNodeIds()) {
        const std::vector<uint8_t> mask = supportMask(tree, nodeId);
        ContourSideCounts counts;
        for (int p = 0; p < rows * columns; ++p) {
            if (!mask[static_cast<std::size_t>(p)]) {
                continue;
            }

            const auto [row, column] = ImageUtils::to2D(p, columns);
            auto sideIsExposed = [&](int qRow, int qColumn) {
                if (qRow < 0 || qRow >= rows || qColumn < 0 || qColumn >= columns) {
                    return true;
                }
                const int q = ImageUtils::to1D(qRow, qColumn, columns);
                return !mask[static_cast<std::size_t>(q)];
            };

            const int north = sideIsExposed(row - 1, column) ? 1 : 0;
            const int west = sideIsExposed(row, column - 1) ? 1 : 0;
            const int east = sideIsExposed(row, column + 1) ? 1 : 0;
            const int south = sideIsExposed(row + 1, column) ? 1 : 0;
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
    return ContourSideFiniteWindowComputation::projectContourPixels(expectedContourSideCounts(tree));
}

std::vector<BitquadStateHistogram> expectedBitquadHistograms(const MorphologicalTree& tree) {
    const int rows = tree.numRows();
    const int columns = tree.numColumns();
    const std::array<Offset, 4> cellOffsets = {{
        {0, 0},
        {0, 1},
        {1, 0},
        {1, 1},
    }};

    std::vector<BitquadStateHistogram> expected(static_cast<std::size_t>(tree.numInternalNodeSlots()));
    for (NodeId nodeId : tree.aliveNodeIds()) {
        const std::vector<uint8_t> mask = supportMask(tree, nodeId);
        auto& hist = expected[static_cast<std::size_t>(nodeId)];
        for (int row = -1; row < rows; ++row) {
            for (int column = -1; column < columns; ++column) {
                BitquadCode code = 0;
                for (std::size_t bit = 0; bit < cellOffsets.size(); ++bit) {
                    const int qRow = row + cellOffsets[bit].rowOffset;
                    const int qColumn = column + cellOffsets[bit].columnOffset;
                    if (qRow < 0 || qRow >= rows || qColumn < 0 || qColumn >= columns) {
                        continue;
                    }
                    const int q = ImageUtils::to1D(qRow, qColumn, columns);
                    if (mask[static_cast<std::size_t>(q)]) {
                        code = static_cast<BitquadCode>(code | (BitquadCode{1} << bit));
                    }
                }
                hist.count(code) += 1;
            }
        }
    }
    return expected;
}

void requireBitquadStateHistogramEqual(const BitquadStateHistogram& actual, const BitquadStateHistogram& expected, const std::string& label) {
    requireVectorEqual(std::vector<int>(actual.bins.begin(), actual.bins.end()), std::vector<int>(expected.bins.begin(), expected.bins.end()), label);
}

void requireBitquadFamilyCountsEqual(const BitquadFamilyCounts& actual, const BitquadFamilyCounts& expected, const std::string& label) {
    requireEqual(actual.q1, expected.q1, label + " q1");
    requireEqual(actual.q2, expected.q2, label + " q2");
    requireEqual(actual.qd, expected.qd, label + " qd");
    requireEqual(actual.q3, expected.q3, label + " q3");
    requireEqual(actual.q4, expected.q4, label + " q4");
}

void requireBitquadFamilyIncrementEqual(const BitquadFamilyIncrement& actual, const BitquadFamilyIncrement& expected, const std::string& label) {
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
    case ContourPixels:
        return static_cast<float>(counts.contourPixels);
    case ContourPerimeter:
        return static_cast<float>(counts.exposedSides);
    case ContourSideNorth:
        return static_cast<float>(counts.north);
    case ContourSideWest:
        return static_cast<float>(counts.west);
    case ContourSideEast:
        return static_cast<float>(counts.east);
    case ContourSideSouth:
        return static_cast<float>(counts.south);
    default:
        throw std::runtime_error("Unsupported projected contour-side scalar attribute.");
    }
}

std::vector<BitquadFamilyCounts> makeQ1ConnectivityProbeCounts(const MorphologicalTree& tree) {
    std::vector<BitquadFamilyCounts> counts(static_cast<std::size_t>(tree.numInternalNodeSlots()));
    for (NodeId nodeId : tree.aliveNodeIds()) {
        counts[static_cast<std::size_t>(nodeId)].q1 = 4;
    }
    return counts;
}

std::vector<float> materializeTreeOfShapesQ1ProbeEuler(const ValuedMorphologicalTree<ToSGrayLevel>& valuedTree) {
    const MorphologicalTree& tree = valuedTree.topology();
    const std::vector<Attribute> requestedAttributes = {BitquadNumberEuler};
    const AttributeNames names = makeDenseAttributeNames(requestedAttributes);
    std::vector<float> buffer(static_cast<std::size_t>(tree.numInternalNodeSlots()) * static_cast<std::size_t>(names.NUM_ATTRIBUTES), 0.0f);

    const auto counts = makeQ1ConnectivityProbeCounts(tree);
    const BitquadConnectivityPolicy connectivityPolicy = makeBitquadConnectivityPolicy(tree, true);
    BitquadAttributeProjection::materializeBitquadAttributes(tree, valuedTree.nodeAltitudeSpan(), std::span<const BitquadFamilyCounts>(counts),
                                                             connectivityPolicy, std::span<float>(buffer), names, requestedAttributes);
    return buffer;
}

template <class T> std::vector<T> copyAltitudeAs(const ValuedMorphologicalTree<ToSGrayLevel>& valuedTree) {
    const NodeAltitudeBuffer<ToSGrayLevel>& altitude = valuedTree.nodeAltitudes();
    std::vector<T> converted;
    converted.reserve(altitude.size());
    for (ToSGrayLevel value : altitude) {
        converted.push_back(static_cast<T>(value));
    }
    return converted;
}

void verifyTreeOfShapesScalarConnectivityPolicy(const ImageUInt8Ptr& image, TestTopographicImmersion interpolation, const std::string& label,
                                                bool requireBothPolarities) {
    const auto valuedTree = makeValuedTreeOfShapes(image, interpolation);
    const MorphologicalTree& tree = valuedTree->topology();
    const AttributeNames names = makeDenseAttributeNames({BitquadNumberEuler});
    const auto buffer = materializeTreeOfShapesQ1ProbeEuler(*valuedTree);

    const auto projectionAdjacencies = mmcfilters::detail::currentBitquadProjectionAdjacencies(tree);
    require(projectionAdjacencies.has_value(), label + " ToS retains a bitquad projection convention");
    const BitquadConnectivityPolicy connectivityPolicy = makeBitquadConnectivityPolicy(tree, true);
    const BitquadConnectivity expectedLowerConnectivity =
        projectionAdjacencies->minAdjacency.is4connectivity() ? BitquadConnectivity::Four : BitquadConnectivity::Eight;
    const BitquadConnectivity expectedUpperConnectivity =
        projectionAdjacencies->maxAdjacency.is4connectivity() ? BitquadConnectivity::Four : BitquadConnectivity::Eight;
    const BitquadConnectivity expectedRootConnectivity =
        expectedLowerConnectivity == expectedUpperConnectivity ? expectedLowerConnectivity : BitquadConnectivity::Eight;
    require(connectivityPolicy.lowerShapeConnectivity() == expectedLowerConnectivity, label + " ToS policy must retain lower-shape connectivity");
    require(connectivityPolicy.upperShapeConnectivity() == expectedUpperConnectivity, label + " ToS policy must retain upper-shape connectivity");
    require(connectivityPolicy.rootConnectivity() == expectedRootConnectivity,
            label + " ToS policy must represent root connectivity independently of shape polarity");
    const std::optional<ShapePolarity> rootPolarity = shapePolarity(tree, valuedTree->nodeAltitudeSpan(), tree.root());
    require(!rootPolarity.has_value(), label + " ToS root must have no shape polarity");
    requireFloatEquivalent(buffer[static_cast<std::size_t>(names.linearIndex(tree.root(), BitquadNumberEuler))],
                           connectivityPolicy.uses4Connectivity(rootPolarity) ? 1.0f : 0.0f, label + " ToS root scalar projection connectivity");

    bool sawLowerShape = false;
    bool sawUpperShape = false;
    for (NodeId nodeId : tree.aliveNodeIds()) {
        if (tree.isRoot(nodeId)) {
            continue;
        }

        const NodeId parentNodeId = tree.parent(nodeId);
        const ToSGrayLevel nodeAltitude = valuedTree->nodeAltitude(nodeId);
        const ToSGrayLevel parentAltitude = valuedTree->nodeAltitude(parentNodeId);
        require(nodeAltitude != parentAltitude, label + " exact ToS must not contain equal-altitude parent-child edges");

        const std::optional<ShapePolarity> polarity = shapePolarity(tree, valuedTree->nodeAltitudeSpan(), nodeId);
        require(polarity.has_value(), label + " ToS non-root node must have shape polarity");
        const bool isUpperShape = *polarity == ShapePolarity::Upper;
        require(isUpperShape == (nodeAltitude > parentAltitude), label + " ToS shape polarity must be derived from exact altitude");
        const bool expected4Connectivity =
            isUpperShape ? projectionAdjacencies->maxAdjacency.is4connectivity() : projectionAdjacencies->minAdjacency.is4connectivity();
        sawUpperShape = sawUpperShape || isUpperShape;
        sawLowerShape = sawLowerShape || !isUpperShape;

        requireEqual(connectivityPolicy.uses4Connectivity(polarity), expected4Connectivity,
                     label + " ToS policy must map lower/upper shape polarity to min/max adjacency " + std::to_string(nodeId));
        requireFloatEquivalent(buffer[static_cast<std::size_t>(names.linearIndex(nodeId, BitquadNumberEuler))], expected4Connectivity ? 1.0f : 0.0f,
                               label + " ToS scalar projection must use explicit lower/upper shape connectivity " + std::to_string(nodeId));
    }

    if (requireBothPolarities) {
        require(sawLowerShape, label + " ToS scalar connectivity policy fixture must contain a lower shape");
        require(sawUpperShape, label + " ToS scalar connectivity policy fixture must contain an upper shape");
    }
}

void verifyTreeOfShapesScalarConnectivityPolicy(const ImageUInt8Ptr& image) {
    verifyTreeOfShapesScalarConnectivityPolicy(image, TestTopographicImmersion::SelfDualSpan, "SelfDual", false);
    verifyTreeOfShapesScalarConnectivityPolicy(image, TestTopographicImmersion::Min4Max8, "Min4cMax8c", true);
    verifyTreeOfShapesScalarConnectivityPolicy(image, TestTopographicImmersion::Min8Max4, "Min8cMax4c", true);
}

void verifyCanonicalBitquadContract() {
    using Family = BitquadFamily;
    const std::array<std::optional<Family>, 16> expected = {{
        std::nullopt,
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

    const auto actual = BitquadFiniteWindowComputation::bitquadFamilyTable();
    for (BitquadCode code = 0; code <= 15; ++code) {
        const std::size_t index = static_cast<std::size_t>(code);
        require(actual[index] == expected[index], "bitquad state family " + std::to_string(code));
        require(BitquadFiniteWindowComputation::bitquadFamily(code) == expected[index], "bitquad family lookup " + std::to_string(code));

        const BitquadState spatialState{(code & (BitquadCode{1} << topLeftBit)) != 0, (code & (BitquadCode{1} << topRightBit)) != 0,
                                        (code & (BitquadCode{1} << bottomLeftBit)) != 0, (code & (BitquadCode{1} << bottomRightBit)) != 0};
        requireEqual(BitquadFiniteWindowComputation::bitquadCode(spatialState), code, "canonical row-major bitquad code " + std::to_string(code));

        const auto lowestVisibleBit = BitquadFiniteWindowComputation::lowestVisibleBitIndex(code);
        if (code == 0) {
            require(!lowestVisibleBit.has_value(), "empty state has no owning bitquad position");
        } else {
            require(lowestVisibleBit.has_value() && *lowestVisibleBit == static_cast<std::size_t>(std::countr_zero(static_cast<unsigned int>(code))),
                    "lowest-index visible bit owns state " + std::to_string(code));
        }

        BitquadStateHistogram oneHot;
        oneHot.count(code) = 1;
        const auto counts = BitquadFiniteWindowComputation::projectBitquadFamilyCounts(oneHot);
        BitquadFamilyCounts expectedCounts;
        if (expected[index].has_value()) {
            switch (*expected[index]) {
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
        }
        requireBitquadFamilyCountsEqual(counts, expectedCounts, "bitquad one-hot family projection " + std::to_string(code));
    }

    requireThrows<std::out_of_range>([]() { static_cast<void>(BitquadFiniteWindowComputation::bitquadFamily(BitquadCode{16})); },
                                     "bitquad family must reject codes outside 0..15");
    requireThrows<std::out_of_range>(
        []() {
            NonemptyBitquadStateHistogramIncrement increment;
            static_cast<void>(increment.count(0));
        },
        "nonempty bitquad state increments must not expose code zero");
    requireThrows<std::out_of_range>(
        []() {
            NonemptyBitquadStateHistogram histogram;
            static_cast<void>(histogram.count(0));
        },
        "nonempty aggregated bitquad histograms must not expose code zero");

    const std::array<std::array<Offset, 4>, 4> expectedWindows{{
        {{{0, 0}, {0, 1}, {1, 0}, {1, 1}}},
        {{{0, -1}, {0, 0}, {1, -1}, {1, 0}}},
        {{{-1, 0}, {-1, 1}, {0, 0}, {0, 1}}},
        {{{-1, -1}, {-1, 0}, {0, -1}, {0, 0}}},
    }};
    for (std::size_t anchorPosition = 0; anchorPosition < expectedWindows.size(); ++anchorPosition) {
        const auto window = BitquadFiniteWindowComputation::bitquadObservationWindow(anchorPosition);
        for (std::size_t coordinate = 0; coordinate < expectedWindows[anchorPosition].size(); ++coordinate) {
            require(window[coordinate] == expectedWindows[anchorPosition][coordinate],
                    "canonical bitquad observation window position " + std::to_string(anchorPosition) + " coordinate " + std::to_string(coordinate));
        }
    }
    requireThrows<std::out_of_range>([]() { static_cast<void>(BitquadFiniteWindowComputation::bitquadObservationWindow(4)); },
                                     "bitquad observation window must reject anchor positions outside 0..3");
}

void verifyExhaustiveBitquadStateHistograms() {
    for (BitquadCode code = 1; code <= 15; ++code) {
        auto image = ImageUInt8::create(2, 2);
        for (std::size_t position = 0; position < 4; ++position) {
            (*image)[static_cast<int>(position)] = (code & (BitquadCode{1} << position)) != 0 ? std::uint8_t{1} : std::uint8_t{0};
        }
        const auto tree = makeComponentTree(image, true);
        const std::size_t firstForegroundPosition = static_cast<std::size_t>(std::countr_zero(static_cast<unsigned int>(code)));
        const NodeId foregroundNode = tree->smallestNode(static_cast<int>(firstForegroundPosition));
        const auto actual = BitquadFiniteWindowComputation::computeBitquadStateHistograms(*tree);
        const auto expected = expectedBitquadHistograms(*tree);

        requireBitquadStateHistogramEqual(actual[static_cast<std::size_t>(foregroundNode)], expected[static_cast<std::size_t>(foregroundNode)],
                                          "exhaustive canonical state histogram code " + std::to_string(code));
        require(actual[static_cast<std::size_t>(foregroundNode)].count(code) > 0, "center framed cell must retain canonical code " + std::to_string(code));
    }

    auto onePixelImage = ImageUInt8::create(1, 1);
    (*onePixelImage)[0] = 1;
    const auto onePixelTree = makeComponentTree(onePixelImage, true);
    const auto histogram = BitquadFiniteWindowComputation::computeBitquadStateHistograms(*onePixelTree)[static_cast<std::size_t>(onePixelTree->root())];
    for (BitquadCode code = 0; code <= 15; ++code) {
        const int expectedCount = code == 1 || code == 2 || code == 4 || code == 8 ? 1 : 0;
        requireEqual(histogram.count(code), expectedCount, "one-pixel framed-boundary canonical code " + std::to_string(code));
    }
}

void verifyFiniteWindowBitquadScalarComputer(const MorphologicalTree& tree, const std::vector<Attribute>& requestedAttributes, const std::string& label) {
    const auto scalarAttributes = runtimeProducedAttributes<BitquadAttributeComputer>();
    requireEqual(scalarAttributes.size(), static_cast<std::size_t>(9), label + " attribute count");

    const AttributeNames names = makeDenseAttributeNames(requestedAttributes);
    const std::size_t bufferSize = static_cast<std::size_t>(tree.numInternalNodeSlots()) * static_cast<std::size_t>(names.NUM_ATTRIBUTES);
    std::vector<float> localBuffer(bufferSize, 0.0f);
    std::vector<float> directFamilyBuffer(bufferSize, 0.0f);

    const auto directBitquadFamilyIncrements = BitquadFiniteWindowComputation::computeBitquadFamilyIncrements(tree);
    const auto directBitquadFamilyCounts = BitquadFiniteWindowComputation::aggregateBitquadFamilyIncrements(tree, directBitquadFamilyIncrements);

    BitquadAttributeComputer::compute(
        AttributeComputeContext<float>{tree, std::span<float>(localBuffer), names, std::span<const Attribute>(requestedAttributes)});
    const BitquadConnectivityPolicy connectivityPolicy = makeBitquadConnectivityPolicy(tree, false);
    BitquadAttributeProjection::materializeBitquadAttributes(tree, std::span<const BitquadFamilyCounts>(directBitquadFamilyCounts), connectivityPolicy,
                                                             std::span<float>(directFamilyBuffer), names, requestedAttributes);

    for (NodeId nodeId : tree.aliveNodeIds()) {
        for (Attribute attribute : requestedAttributes) {
            const int index = names.linearIndex(nodeId, attribute);
            requireFloatEquivalent(localBuffer[static_cast<std::size_t>(index)], directFamilyBuffer[static_cast<std::size_t>(index)],
                                   label + " local compute uses direct family scalar " + AttributeNames::toString(attribute) + " node " +
                                       std::to_string(nodeId));
        }
    }
}

void verifyFiniteWindowBitquadScalarComputer(const MorphologicalTree& tree, const std::string& label) {
    verifyFiniteWindowBitquadScalarComputer(tree, runtimeProducedAttributes<BitquadAttributeComputer>(), label + " full");
    verifyFiniteWindowBitquadScalarComputer(tree, {BitquadArea, BitquadPerimeter, BitquadCircularity}, label + " subset");
}

void verifyFiniteWindowTreeOfShapesBitquadScalarComputer(const ValuedMorphologicalTree<ToSGrayLevel>& valuedTree, const std::string& label) {
    const MorphologicalTree& tree = valuedTree.topology();
    const auto requestedAttributes = runtimeProducedAttributes<BitquadAttributeComputer>();
    const AttributeNames names = makeDenseAttributeNames(requestedAttributes);
    const std::size_t bufferSize = static_cast<std::size_t>(tree.numInternalNodeSlots()) * static_cast<std::size_t>(names.NUM_ATTRIBUTES);
    std::vector<float> computeBuffer(bufferSize, 0.0f);
    std::vector<float> directBuffer(bufferSize, 0.0f);
    std::vector<float> genericSpanBuffer(bufferSize, 0.0f);
    std::vector<float> genericViewBuffer(bufferSize, 0.0f);

    BitquadAttributeComputer::compute(AltitudeAttributeComputeContext<float, ToSGrayLevel>{tree, valuedTree.nodeAltitudeSpan(), std::span<float>(computeBuffer),
                                                                                           names, std::span<const Attribute>(requestedAttributes)});

    const auto directBitquadFamilyCounts = BitquadFiniteWindowComputation::computeBitquadFamilyCounts(tree);
    const BitquadConnectivityPolicy connectivityPolicy = makeBitquadConnectivityPolicy(tree, true);
    BitquadAttributeProjection::materializeBitquadAttributes(tree, valuedTree.nodeAltitudeSpan(), std::span<const BitquadFamilyCounts>(directBitquadFamilyCounts),
                                                             connectivityPolicy, std::span<float>(directBuffer), names, requestedAttributes);

    const std::vector<std::int32_t> externalIntAltitude = copyAltitudeAs<std::int32_t>(valuedTree);
    BitquadAttributeProjection::materializeBitquadAttributes(tree, std::span<const std::int32_t>(externalIntAltitude),
                                                             std::span<const BitquadFamilyCounts>(directBitquadFamilyCounts), connectivityPolicy,
                                                             std::span<float>(genericSpanBuffer), names, requestedAttributes);

    const std::vector<float> externalFloatAltitude = copyAltitudeAs<float>(valuedTree);
    const ValuedMorphologicalTreeView<float> externalFloatView(tree, std::span<const float>(externalFloatAltitude));
    BitquadAttributeProjection::materializeBitquadAttributes(externalFloatView, std::span<const BitquadFamilyCounts>(directBitquadFamilyCounts),
                                                             connectivityPolicy, std::span<float>(genericViewBuffer), names, requestedAttributes);

    for (NodeId nodeId : tree.aliveNodeIds()) {
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

    auto [pipelineNames, pipelineBuffer] = AttributeComputation::computeSingleAttribute(valuedTree, AttributeGroup::Boundary);
    auto [valuedTreeTopologyNames, valuedTreeTopologyBuffer] = AttributeComputation::computeTopologyAttributes(valuedTree, {AttributeGroup::Boundary});
    for (NodeId nodeId : tree.aliveNodeIds()) {
        for (Attribute attribute : requestedAttributes) {
            const int directIndex = names.linearIndex(nodeId, attribute);
            requireFloatEquivalent(pipelineBuffer[static_cast<std::size_t>(pipelineNames.linearIndex(nodeId, attribute))],
                                   directBuffer[static_cast<std::size_t>(directIndex)],
                                   label + " public valuedTree pipeline " + AttributeNames::toString(attribute) + " node " + std::to_string(nodeId));
            requireFloatEquivalent(valuedTreeTopologyBuffer[static_cast<std::size_t>(valuedTreeTopologyNames.linearIndex(nodeId, attribute))],
                                   directBuffer[static_cast<std::size_t>(directIndex)],
                                   label + " public valuedTree topology pipeline " + AttributeNames::toString(attribute) + " node " + std::to_string(nodeId));
        }
    }

    if constexpr (contract::validationsEnabled) {
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(AttributeComputation::computeTopologyAttributes(tree, {AttributeGroup::Boundary})); },
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
                BitquadAttributeProjection::materializeBitquadAttributes(tree, std::span<const float>(shortAltitude),
                                                                         std::span<const BitquadFamilyCounts>(directBitquadFamilyCounts), connectivityPolicy,
                                                                         std::span<float>(throwBuffer), names, requestedAttributes);
            },
            label + " ToS generic scalar projection with wrong altitude size must throw");

        NodeId nonRootNode = InvalidNode;
        for (NodeId nodeId : tree.aliveNodeIds()) {
            if (!tree.isRoot(nodeId)) {
                nonRootNode = nodeId;
                break;
            }
        }
        if (nonRootNode != InvalidNode) {
            std::vector<float> ambiguousAltitude = externalFloatAltitude;
            ambiguousAltitude[static_cast<std::size_t>(nonRootNode)] = ambiguousAltitude[static_cast<std::size_t>(tree.parent(nonRootNode))];
            requireThrows<std::runtime_error>(
                [&]() {
                    std::vector<float> ambiguousBuffer(bufferSize, 0.0f);
                    BitquadAttributeProjection::materializeBitquadAttributes(tree, std::span<const float>(ambiguousAltitude),
                                                                             std::span<const BitquadFamilyCounts>(directBitquadFamilyCounts),
                                                                             connectivityPolicy, std::span<float>(ambiguousBuffer), names, requestedAttributes);
                },
                label + " dual adjacency must reject an altitude-equal ambiguous "
                        "branch");
        }
    }
}

void verifyGenericContourCounts(const MorphologicalTree& tree, const std::string& label) {
    const auto actualSideCounts = ContourSideFiniteWindowComputation::computeContourSideCounts(tree);
    const auto expectedSideCounts = expectedContourSideCounts(tree);
    const auto actualProjected = ContourSideFiniteWindowComputation::projectContourPixels(actualSideCounts);
    const auto actualExposedSides = ContourSideFiniteWindowComputation::projectExposedSides(actualSideCounts);
    const auto expected = expectedContourCounts(tree);
    for (NodeId nodeId : tree.aliveNodeIds()) {
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
        ContourPixels, ContourPerimeter, ContourSideNorth, ContourSideWest, ContourSideEast, ContourSideSouth,
    };
    requireEqual(allContourAttributes.size(), expectedAttributes.size(), label + " contour attribute count");
    for (std::size_t i = 0; i < expectedAttributes.size(); ++i) {
        requireEqual(static_cast<int>(allContourAttributes[i]), static_cast<int>(expectedAttributes[i]),
                     label + " contour attribute order " + std::to_string(i));
    }

    const AttributeNames names = makeDenseAttributeNames(requestedAttributes);
    const std::size_t bufferSize = static_cast<std::size_t>(tree.numInternalNodeSlots()) * static_cast<std::size_t>(names.NUM_ATTRIBUTES);
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

    for (NodeId nodeId : tree.aliveNodeIds()) {
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
    verifyPublicContourSideCountsComputer(tree, {ContourPixels, ContourPerimeter, ContourSideEast}, label + " subset");
}

void verifyPublicContourSideCountUnitAttributes(const MorphologicalTree& tree, const std::string& label) {
    const auto attributes = runtimeProducedAttributes<ContourSideAttributeComputer>();
    const AttributeNames names = makeDenseAttributeNames(attributes);
    const std::vector<PixelId> unitPixels = {
        0,
        static_cast<PixelId>(tree.numPixels() - 1),
    };
    std::vector<float> buffer(unitPixels.size() * static_cast<std::size_t>(names.NUM_ATTRIBUTES), 0.0f);

    ContourSideAttributeComputer::computeUnitRows(UnitAttributeComputeContext<float>{tree, std::span<const PixelId>(unitPixels), std::span<float>(buffer),
                                                                                     names, std::span<const Attribute>(attributes)});

    ContourSideCounts expected;
    expected.contourPixels = 1;
    expected.exposedSides = 4;
    expected.north = 1;
    expected.west = 1;
    expected.east = 1;
    expected.south = 1;
    for (NodeId leafIndex = 0; leafIndex < static_cast<NodeId>(unitPixels.size()); ++leafIndex) {
        for (Attribute attribute : attributes) {
            const int index = names.linearIndex(leafIndex, attribute);
            requireFloatEquivalent(buffer[static_cast<std::size_t>(index)], projectedContourSideScalarValue(expected, attribute),
                                   label + " unit " + AttributeNames::toString(attribute) + " leaf " + std::to_string(leafIndex));
        }
    }
}

void verifyPublicContourAttributeNames() {
    requireEqual(AttributeNames::toString(ContourPixels), std::string("CONTOUR_PIXELS"), "contour pixels name");
    requireEqual(AttributeNames::toString(ContourPerimeter), std::string("CONTOUR_PERIMETER"), "contour perimeter name");
    require(AttributeNames::describe(ContourSideNorth).starts_with("Contour north sides:"), "contour north description");
    auto parsed = AttributeNames::parse("CONTOUR_PIXELS");
    require(parsed.has_value(), "public contour attributes must be parsed");
    requireEqual(static_cast<int>(*parsed), static_cast<int>(ContourPixels), "public contour parse value");
}

void verifyBitquadIncrementPipeline(const MorphologicalTree& tree, const std::string& label) {
    const auto stateIncrements = BitquadFiniteWindowComputation::computeNonemptyBitquadStateHistogramIncrements(tree);
    const auto familyIncrements = BitquadFiniteWindowComputation::computeBitquadFamilyIncrements(tree);
    std::vector<BitquadFamilyIncrement> familyIncrementsFromStates;
    familyIncrementsFromStates.reserve(stateIncrements.size());
    for (const auto& stateIncrement : stateIncrements) {
        familyIncrementsFromStates.push_back(BitquadFiniteWindowComputation::projectBitquadFamilyIncrement(stateIncrement));
    }

    const auto nonemptyStates = BitquadFiniteWindowComputation::aggregateNonemptyBitquadStateHistogramIncrements(tree, stateIncrements);
    const auto aggregatedStates = BitquadFiniteWindowComputation::materializeEmptyBitquadCount(tree, nonemptyStates);
    const auto aggregatedFamilies = BitquadFiniteWindowComputation::aggregateBitquadFamilyIncrements(tree, familyIncrements);
    const auto directStates = BitquadFiniteWindowComputation::computeBitquadStateHistograms(tree);
    const auto directFamilies = BitquadFiniteWindowComputation::computeBitquadFamilyCounts(tree);

    for (NodeId nodeId : tree.aliveNodeIds()) {
        const std::size_t nodeIndex = static_cast<std::size_t>(nodeId);
        requireBitquadStateHistogramEqual(aggregatedStates[nodeIndex], directStates[nodeIndex],
                                          label + " bitquad state increments aggregate before empty materialization node " + std::to_string(nodeId));
        requireBitquadFamilyCountsEqual(aggregatedFamilies[nodeIndex], directFamilies[nodeIndex],
                                        label + " bitquad family increments aggregate to direct counts node " + std::to_string(nodeId));
        requireBitquadFamilyIncrementEqual(familyIncrementsFromStates[nodeIndex], familyIncrements[nodeIndex],
                                           label + " bitquad family increments match projected state increments node " + std::to_string(nodeId));

        int nonemptyCount = 0;
        for (BitquadCode code = 1; code <= 15; ++code) {
            nonemptyCount += nonemptyStates[nodeIndex].count(code);
        }
        const int expectedEmptyCount = (tree.numRows() + 1) * (tree.numColumns() + 1) - nonemptyCount;
        requireEqual(aggregatedStates[nodeIndex].count(0), expectedEmptyCount,
                     label + " empty state is materialized by complement node " + std::to_string(nodeId));
    }

    const auto properPartStates = BitquadFiniteWindowComputation::projectBitquadStateHistogramsToProperParts(tree, directStates);
    const auto properPartFamilies = BitquadFiniteWindowComputation::projectBitquadFamilyCountsToProperParts(tree, directFamilies);
    const auto properPartFamilyIncrements = BitquadFiniteWindowComputation::projectBitquadFamilyIncrementsToProperParts(tree, familyIncrements);
    for (PixelId pixel = 0; pixel < tree.numPixels(); ++pixel) {
        const NodeId smallestNodeId = tree.smallestNode(pixel);
        const std::size_t pixelIndex = static_cast<std::size_t>(pixel);
        const std::size_t smallestNodeIndex = static_cast<std::size_t>(smallestNodeId);
        requireBitquadStateHistogramEqual(properPartStates[pixelIndex], directStates[smallestNodeIndex],
                                          label + " proper-part state projection " + std::to_string(pixel));
        requireBitquadFamilyCountsEqual(properPartFamilies[pixelIndex], directFamilies[smallestNodeIndex],
                                        label + " proper-part family projection " + std::to_string(pixel));
        requireBitquadFamilyIncrementEqual(properPartFamilyIncrements[pixelIndex], familyIncrements[smallestNodeIndex],
                                           label + " proper-part family-increment projection " + std::to_string(pixel));
    }

    std::vector<BitquadFamilyIncrement> shortFamilyIncrements = familyIncrements;
    if (!shortFamilyIncrements.empty()) {
        shortFamilyIncrements.pop_back();
    }
    requireThrows<std::invalid_argument>(
        [&]() { static_cast<void>(BitquadFiniteWindowComputation::aggregateBitquadFamilyIncrements(tree, shortFamilyIncrements)); },
        label + " bitquad family increment aggregation must validate exact node-slot coverage");

    std::vector<NonemptyBitquadStateHistogramIncrement> shortStateIncrements = stateIncrements;
    if (!shortStateIncrements.empty()) {
        shortStateIncrements.pop_back();
    }
    requireThrows<std::invalid_argument>(
        [&]() { static_cast<void>(BitquadFiniteWindowComputation::aggregateNonemptyBitquadStateHistogramIncrements(tree, shortStateIncrements)); },
        label + " nonempty bitquad state increment aggregation must validate exact node-slot coverage");

    std::vector<NonemptyBitquadStateHistogram> shortNonemptyStates = nonemptyStates;
    if (!shortNonemptyStates.empty()) {
        shortNonemptyStates.pop_back();
    }
    requireThrows<std::invalid_argument>([&]() { static_cast<void>(BitquadFiniteWindowComputation::materializeEmptyBitquadCount(tree, shortNonemptyStates)); },
                                         label + " empty bitquad count materialization must validate exact node-slot coverage");
}

void verifyGenericBitquadHistograms(const MorphologicalTree& tree, const std::string& label) {
    BitquadFiniteWindowComputation computer(tree);
    const auto& actual = computer.bitquadStateHistograms();
    const auto& actualBitquadFamilyCounts = computer.bitquadFamilyCounts();
    const auto directBitquadFamilyCounts = BitquadFiniteWindowComputation::computeBitquadFamilyCounts(tree);
    const auto expected = expectedBitquadHistograms(tree);
    for (NodeId nodeId : tree.aliveNodeIds()) {
        requireBitquadStateHistogramEqual(actual[static_cast<std::size_t>(nodeId)], expected[static_cast<std::size_t>(nodeId)],
                                          label + " bitquad histogram node " + std::to_string(nodeId));
        requireBitquadFamilyCountsEqual(actualBitquadFamilyCounts[static_cast<std::size_t>(nodeId)],
                                        BitquadFiniteWindowComputation::projectBitquadFamilyCounts(expected[static_cast<std::size_t>(nodeId)]),
                                        label + " bitquad family count node " + std::to_string(nodeId));
        requireBitquadFamilyCountsEqual(directBitquadFamilyCounts[static_cast<std::size_t>(nodeId)],
                                        actualBitquadFamilyCounts[static_cast<std::size_t>(nodeId)],
                                        label + " direct bitquad family count node " + std::to_string(nodeId));
    }
    verifyBitquadIncrementPipeline(tree, label);
}

MorphologicalTree makeTwoBranchTreeOfShapes() {
    std::vector<NodeId> parent = {
        4, 4, 5, 5, 6, 6, 6,
    };
    std::vector<std::uint8_t> altitude(parent.size(), std::uint8_t{});
    auto valuedTree = MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(parent), std::span<const std::uint8_t>(altitude), 2, 2,
                                                                    MorphologicalTreeKind::TreeOfShapes);
    return valuedTree.topology().clone();
}

MorphologicalTree makeThreeNodeChain() {
    const std::vector<NodeId> parent = {3, 4, 5, 4, 5, 5};
    const std::vector<std::uint8_t> altitude = {2, 1, 0, 2, 1, 0};
    auto valuedTree = MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(parent), std::span<const std::uint8_t>(altitude), 1, 3,
                                                                    MorphologicalTreeKind::MaxTree, RegularGridAdjacency2D(1, 3, 1.0));
    return valuedTree.topology().clone();
}

struct VisibleSampleCountLocalRule {
    using Value = int;

    [[nodiscard]] Value additiveIdentity() const { return 0; }
    [[nodiscard]] Value evaluateLocalRule(BinaryVisibilityState state) const { return static_cast<int>(std::popcount(state.bits())); }
    void addAssign(Value& target, const Value& source) const { target += source; }
    void subtractAssign(Value& target, const Value& source) const { target -= source; }
};

static_assert(LocalRule<VisibleSampleCountLocalRule>);

struct SignedVisibilityLocalRule {
    using Value = int;

    [[nodiscard]] Value additiveIdentity() const { return 0; }
    [[nodiscard]] Value evaluateLocalRule(BinaryVisibilityState state) const {
        const int numVisibleSamples = static_cast<int>(std::popcount(state.bits()));
        return (numVisibleSamples % 2 == 0) ? -numVisibleSamples : 2 * numVisibleSamples + 1;
    }
    void addAssign(Value& target, const Value& source) const { target += source; }
    void subtractAssign(Value& target, const Value& source) const { target -= source; }
};

static_assert(LocalRule<SignedVisibilityLocalRule>);

template <class Value> std::vector<Value> nodeAttributeValues(std::span<const NodeAttribute<Value>> nodeAttributes) {
    std::vector<Value> values;
    values.reserve(nodeAttributes.size());
    for (std::size_t slot = 0; slot < nodeAttributes.size(); ++slot) {
        requireEqual(nodeAttributes[slot].node, static_cast<NodeId>(slot), "node attributes use dense node order");
        values.push_back(nodeAttributes[slot].value);
    }
    return values;
}

MorphologicalTree makeRandomSmallTree(std::mt19937& generator) {
    constexpr int numPixels = 6;
    std::vector<NodeId> parent(static_cast<std::size_t>(2 * numPixels - 1), InvalidNode);
    std::vector<NodeId> activeNodes;
    activeNodes.reserve(numPixels);
    for (PixelId pixel = 0; pixel < numPixels; ++pixel) {
        activeNodes.push_back(pixel);
    }

    NodeId nextNode = numPixels;
    while (activeNodes.size() > 1) {
        std::uniform_int_distribution<std::size_t> firstDistribution(0, activeNodes.size() - 1);
        const std::size_t firstIndex = firstDistribution(generator);
        NodeId first = activeNodes[firstIndex];
        activeNodes.erase(activeNodes.begin() + static_cast<std::ptrdiff_t>(firstIndex));

        std::uniform_int_distribution<std::size_t> secondDistribution(0, activeNodes.size() - 1);
        const std::size_t secondIndex = secondDistribution(generator);
        NodeId second = activeNodes[secondIndex];
        activeNodes.erase(activeNodes.begin() + static_cast<std::ptrdiff_t>(secondIndex));

        parent[static_cast<std::size_t>(first)] = nextNode;
        parent[static_cast<std::size_t>(second)] = nextNode;
        activeNodes.push_back(nextNode);
        ++nextNode;
    }
    parent[static_cast<std::size_t>(activeNodes.front())] = activeNodes.front();

    std::vector<std::uint8_t> altitude(parent.size(), std::uint8_t{});
    auto valuedTree = MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(parent), std::span<const std::uint8_t>(altitude), 2, 3,
                                                                    MorphologicalTreeKind::TreeOfShapes);
    return valuedTree.topology().clone();
}

std::vector<int> directFiniteWindowAttribute(const MorphologicalTree& tree, const ObservationWindow& observationWindow,
                                             const SignedVisibilityLocalRule& localRule) {
    std::vector<int> expected(static_cast<std::size_t>(tree.numInternalNodeSlots()), localRule.additiveIdentity());
    for (NodeId node : tree.aliveNodeIds()) {
        const std::vector<uint8_t> mask = supportMask(tree, node);
        for (PixelId anchorPixel = 0; anchorPixel < tree.numPixels(); ++anchorPixel) {
            if (mask[static_cast<std::size_t>(anchorPixel)] == 0) {
                continue;
            }

            const int anchorRow = anchorPixel / tree.numColumns();
            const int anchorColumn = anchorPixel % tree.numColumns();
            std::uint32_t bits = 0;
            for (std::size_t coordinate = 0; coordinate < observationWindow.size(); ++coordinate) {
                const Offset offset = observationWindow[coordinate];
                const int sampleRow = anchorRow + offset.rowOffset;
                const int sampleColumn = anchorColumn + offset.columnOffset;
                if (sampleRow < 0 || sampleRow >= tree.numRows() || sampleColumn < 0 || sampleColumn >= tree.numColumns()) {
                    continue;
                }
                const PixelId samplePixel = sampleRow * tree.numColumns() + sampleColumn;
                if (mask[static_cast<std::size_t>(samplePixel)] != 0) {
                    bits |= std::uint32_t{1} << coordinate;
                }
            }
            expected[static_cast<std::size_t>(node)] += localRule.evaluateLocalRule(BinaryVisibilityState(bits, observationWindow.size()));
        }
    }
    return expected;
}

void verifyObservationWindowContract() {
    const ObservationWindow window{{0, 0}, {-1, 2}, {3, -4}};
    requireEqual(window.size(), std::size_t{3}, "observation window size");
    requireEqual(window[1].rowOffset, -1, "observation window row offset");
    requireEqual(window[1].columnOffset, 2, "observation window column offset");

    requireThrows<std::invalid_argument>([]() { static_cast<void>(ObservationWindow(std::vector<Offset>{})); }, "empty observation window must be rejected");
    requireThrows<std::invalid_argument>([]() { static_cast<void>(ObservationWindow{{0, 1}}); }, "observation window without zero offset must be rejected");
    requireThrows<std::invalid_argument>([]() { static_cast<void>(ObservationWindow{{0, 0}, {0, 1}, {0, 1}}); },
                                         "duplicate observation offset must be rejected");
    requireThrows<std::invalid_argument>([]() { static_cast<void>(ObservationWindow{{0, 0}, {0, 0}}); }, "duplicate zero offset must be rejected");

    std::vector<Offset> oversized{{0, 0}};
    for (int column = 1; column <= 32; ++column) {
        oversized.push_back({0, column});
    }
    requireThrows<std::invalid_argument>([&]() { static_cast<void>(ObservationWindow(oversized)); },
                                         "observation window with more than 32 offsets must be rejected");
}

void verifyMaximumObservationWindowComputation() {
    const auto tree = makeComponentTree(makeRampImage(1, 32), true, 1.0);
    std::vector<Offset> offsets;
    offsets.reserve(ObservationWindow::maxNumOffsets);
    for (int column = 0; column < static_cast<int>(ObservationWindow::maxNumOffsets); ++column) {
        offsets.push_back({0, column});
    }
    const ObservationWindow observationWindow(std::move(offsets));
    const SignedVisibilityLocalRule localRule;
    const std::vector<int> expected = directFiniteWindowAttribute(*tree, observationWindow, localRule);
    requireVectorEqual(nodeAttributeValues<int>(FiniteWindowLocalAttributeComputer::compute(*tree, observationWindow, localRule)), expected,
                       "maximum-width observation window equals direct definition");
}

void verifyFiniteWindowStructuralContract() {
    using Computer = FiniteWindowLocalAttributeComputer;

    const MorphologicalTree chain = makeThreeNodeChain();
    requireVectorEqual(Computer::anchorBranch(chain, 0), std::vector<NodeId>{0, 1, 2}, "anchor branch follows increasing inclusion");

    const auto equalEntry = Computer::anchoredEntry(chain, 0, 0);
    require(equalEntry.has_value() && *equalEntry == 0, "equal smallest nodes enter at the anchor smallest node");
    const auto ancestorEntry = Computer::anchoredEntry(chain, 0, 1);
    require(ancestorEntry.has_value() && *ancestorEntry == 1, "ancestor sample enters at its smallest node");
    const auto descendantEntry = Computer::anchoredEntry(chain, 1, 0);
    require(descendantEntry.has_value() && *descendantEntry == 1, "descendant sample enters at the anchor smallest node");
    require(!Computer::anchoredEntry(chain, 0, Offset{0, -1}).has_value(), "out-of-domain sample has no anchored entry");

    const ObservationWindow chainWindow{{0, 0}, {0, 1}, {0, 2}, {0, -1}};
    const AnchoredEntryMap chainEntries = Computer::anchoredEntryMap(chain, 0, chainWindow);
    require(chainEntries[0].has_value() && *chainEntries[0] == 0, "zero-offset anchored entry");
    require(chainEntries[1].has_value() && *chainEntries[1] == 1, "middle-chain anchored entry");
    require(chainEntries[2].has_value() && *chainEntries[2] == 2, "root anchored entry");
    require(!chainEntries[3].has_value(), "out-of-domain coordinate remains missing in anchored-entry map");

    const AnchoredEntrySet chainEntrySet = Computer::anchoredEntrySet(chain, 0, chainWindow);
    require(chainEntrySet == AnchoredEntrySet{0, 1, 2}, "anchored-entry set contains distinct live entries");
    const OrderedAnchoredEntries chainOrderedEntries = Computer::orderAnchoredEntriesByInclusion(chain, chainEntries);
    requireEqual(chainOrderedEntries.size(), std::size_t{3}, "ordered anchored-entry count");
    require(chainOrderedEntries[0] == AnchoredEntryMask{0, 1}, "first anchored-entry mask");
    require(chainOrderedEntries[1] == AnchoredEntryMask{1, 2}, "second anchored-entry mask");
    require(chainOrderedEntries[2] == AnchoredEntryMask{2, 4}, "third anchored-entry mask");

    requireEqual(Computer::binaryVisibilityState(chain, chainEntries, 1).bits(), std::uint32_t{0b0011}, "visibility state on middle chain node");
    requireEqual(Computer::binaryVisibilityState(chain, chainEntries, 2).bits(), std::uint32_t{0b0111}, "visibility state on root");

    const MorphologicalTree twoBranchTree = makeTwoBranchTreeOfShapes();
    const auto incomparableEntry = Computer::anchoredEntry(twoBranchTree, 0, 2);
    require(incomparableEntry.has_value() && *incomparableEntry == twoBranchTree.root(), "incomparable sample enters at the LCA");

    const ObservationWindow groupedWindow{{0, 0}, {0, 1}, {1, 0}};
    const OrderedAnchoredEntries groupedEntries = Computer::orderedAnchoredEntries(twoBranchTree, 0, groupedWindow);
    requireEqual(groupedEntries.size(), std::size_t{2}, "equal anchored entries are grouped");
    require(groupedEntries[0] == AnchoredEntryMask{0, std::uint32_t{0b011}}, "equal entry coordinates share one mask");
    require(groupedEntries[1] == AnchoredEntryMask{twoBranchTree.root(), std::uint32_t{0b100}}, "LCA entry follows anchor entry");
}

void verifyAdditiveLocalRuleComputation() {
    const MorphologicalTree tree = makeThreeNodeChain();
    const ObservationWindow window{{0, 0}, {0, 1}};
    const auto computedNodeAttributes = FiniteWindowLocalAttributeComputer::compute(tree, window, VisibleSampleCountLocalRule{});
    const std::vector<int> actual = nodeAttributeValues<int>(computedNodeAttributes);

    std::vector<int> expected(static_cast<std::size_t>(tree.numInternalNodeSlots()), 0);
    for (NodeId node : tree.aliveNodeIds()) {
        const std::vector<uint8_t> mask = supportMask(tree, node);
        for (PixelId anchorPixel = 0; anchorPixel < tree.numPixels(); ++anchorPixel) {
            if (mask[static_cast<std::size_t>(anchorPixel)] == 0) {
                continue;
            }
            for (Offset offset : window) {
                const int sampleColumn = anchorPixel + offset.columnOffset;
                if (sampleColumn >= 0 && sampleColumn < tree.numColumns() && mask[static_cast<std::size_t>(sampleColumn)] != 0) {
                    expected[static_cast<std::size_t>(node)] += 1;
                }
            }
        }
    }
    requireVectorEqual(actual, expected, "additive local rule equals direct finite-window definition");
}

void verifyExplicitFiniteWindowPipeline() {
    using Computer = FiniteWindowLocalAttributeComputer;
    const MorphologicalTree tree = makeThreeNodeChain();
    const ObservationWindow window{{0, 0}, {0, 1}};
    const SignedVisibilityLocalRule localRule;

    const auto firstAnchorEventDeltas = Computer::computeEventDeltas(tree, 0, window, localRule);
    requireEqual(firstAnchorEventDeltas.size(), std::size_t{2}, "first anchor event-delta count");
    requireEqual(firstAnchorEventDeltas[0].anchorPixel, 0, "first event-delta anchor");
    requireEqual(firstAnchorEventDeltas[0].anchoredEntry, NodeId{0}, "first event-delta entry");
    requireEqual(firstAnchorEventDeltas[0].value, 3, "first event delta is measured from the additive identity");
    requireEqual(firstAnchorEventDeltas[1].anchoredEntry, NodeId{1}, "second event-delta entry");
    requireEqual(firstAnchorEventDeltas[1].value, -5, "later event delta is a consecutive signed rule difference");

    const auto localAttributeIncrements = Computer::computeLocalAttributeIncrements(tree, window, localRule);
    requireEqual(localAttributeIncrements.size(), static_cast<std::size_t>(tree.numInternalNodeSlots()), "local-attribute increments cover dense node slots");
    requireEqual(localAttributeIncrements[0].value, 3, "leaf local-attribute increment");
    requireEqual(localAttributeIncrements[1].value, -2, "middle signed local-attribute increment");
    requireEqual(localAttributeIncrements[2].value, -2, "root signed local-attribute increment");

    const auto aggregated =
        Computer::aggregateLocalAttributeIncrements(tree, std::span<const LocalAttributeIncrement<int>>(localAttributeIncrements), localRule);
    requireVectorEqual(nodeAttributeValues<int>(aggregated), std::vector<int>{3, 1, -1}, "bottom-up aggregation materializes final node attributes");
    requireVectorEqual(nodeAttributeValues<int>(Computer::compute(tree, window, localRule)), nodeAttributeValues<int>(aggregated),
                       "complete finite-window computation follows the explicit pipeline");

    requireThrows<std::invalid_argument>([&]() { static_cast<void>(Computer::computeEventDeltas(tree, -1, window, localRule)); },
                                         "event-delta computation rejects an invalid anchor");
    auto malformedIncrements = localAttributeIncrements;
    malformedIncrements[0].node = 1;
    requireThrows<std::invalid_argument>(
        [&]() {
            static_cast<void>(Computer::aggregateLocalAttributeIncrements(tree, std::span<const LocalAttributeIncrement<int>>(malformedIncrements), localRule));
        },
        "bottom-up aggregation rejects increments outside dense node order");
}

void verifyRandomizedFiniteWindowReference() {
    std::mt19937 generator(0x5a17u);
    const std::array<Offset, 8> candidateOffsets = {{{-1, -1}, {-1, 0}, {-1, 1}, {0, -1}, {0, 1}, {1, -1}, {1, 0}, {1, 1}}};
    const SignedVisibilityLocalRule localRule;

    for (int iteration = 0; iteration < 64; ++iteration) {
        const MorphologicalTree tree = makeRandomSmallTree(generator);
        std::vector<Offset> shuffledOffsets(candidateOffsets.begin(), candidateOffsets.end());
        std::shuffle(shuffledOffsets.begin(), shuffledOffsets.end(), generator);
        const std::size_t numAdditionalOffsets = 1 + static_cast<std::size_t>(generator() % 5);
        shuffledOffsets.resize(numAdditionalOffsets);
        shuffledOffsets.insert(shuffledOffsets.begin() + static_cast<std::ptrdiff_t>(generator() % (numAdditionalOffsets + 1)), Offset{0, 0});
        const ObservationWindow observationWindow(std::move(shuffledOffsets));

        const auto localAttributeIncrements = FiniteWindowLocalAttributeComputer::computeLocalAttributeIncrements(tree, observationWindow, localRule);
        std::vector<int> independentlySummedIncrements(static_cast<std::size_t>(tree.numInternalNodeSlots()), 0);
        for (PixelId anchorPixel = 0; anchorPixel < tree.numPixels(); ++anchorPixel) {
            for (const EventDelta<int>& eventDelta : FiniteWindowLocalAttributeComputer::computeEventDeltas(tree, anchorPixel, observationWindow, localRule)) {
                independentlySummedIncrements[static_cast<std::size_t>(eventDelta.anchoredEntry)] += eventDelta.value;
            }
        }
        for (std::size_t node = 0; node < localAttributeIncrements.size(); ++node) {
            requireEqual(localAttributeIncrements[node].value, independentlySummedIncrements[node], "randomized event deltas sum to local-attribute increment");
        }

        const auto nodeAttributes = FiniteWindowLocalAttributeComputer::aggregateLocalAttributeIncrements(
            tree, std::span<const LocalAttributeIncrement<int>>(localAttributeIncrements), localRule);
        const std::vector<int> expected = directFiniteWindowAttribute(tree, observationWindow, localRule);
        requireVectorEqual(nodeAttributeValues<int>(nodeAttributes), expected, "randomized explicit pipeline equals direct finite-window definition");
        requireVectorEqual(nodeAttributeValues<int>(FiniteWindowLocalAttributeComputer::compute(tree, observationWindow, localRule)), expected,
                           "randomized complete compute equals direct finite-window definition");
    }
}

} // namespace

int main() {
    verifyObservationWindowContract();
    verifyMaximumObservationWindowComputation();
    verifyFiniteWindowStructuralContract();
    verifyAdditiveLocalRuleComputation();
    verifyExplicitFiniteWindowPipeline();
    verifyRandomizedFiniteWindowReference();
    verifyCanonicalBitquadContract();
    verifyExhaustiveBitquadStateHistograms();
    verifyPublicContourAttributeNames();

    auto image = makeComponentTreeFixture();
    auto rampImage = makeRampImage(16, 16);

    {
        auto tree = makeComponentTree(image, true);
        verifyGenericContourCounts(*tree, "max-tree");
        verifyPublicContourSideCountsComputer(*tree, "max-tree");
        verifyPublicContourSideCountUnitAttributes(*tree, "max-tree");
        verifyGenericBitquadHistograms(*tree, "max-tree");
        verifyFiniteWindowBitquadScalarComputer(*tree, "max-tree");
    }

    {
        auto tree = makeComponentTree(image, false);
        verifyGenericContourCounts(*tree, "min-tree");
        verifyPublicContourSideCountsComputer(*tree, "min-tree");
        verifyGenericBitquadHistograms(*tree, "min-tree");
        verifyFiniteWindowBitquadScalarComputer(*tree, "min-tree");
    }

    {
        auto tree = makeComponentTree(image, true, 1.0);
        verifyGenericContourCounts(*tree, "max-tree 4-connectivity");
        verifyPublicContourSideCountsComputer(*tree, "max-tree 4-connectivity");
        verifyGenericBitquadHistograms(*tree, "max-tree 4-connectivity");
        verifyFiniteWindowBitquadScalarComputer(*tree, "max-tree 4-connectivity");
    }

    {
        auto tree = makeComponentTree(image, false, 1.0);
        verifyGenericContourCounts(*tree, "min-tree 4-connectivity");
        verifyPublicContourSideCountsComputer(*tree, "min-tree 4-connectivity");
        verifyGenericBitquadHistograms(*tree, "min-tree 4-connectivity");
        verifyFiniteWindowBitquadScalarComputer(*tree, "min-tree 4-connectivity");
    }

    {
        auto tree = makeComponentTree(rampImage, true, 1.0);
        verifyGenericBitquadHistograms(*tree, "ramp max-tree 4-connectivity");
        verifyFiniteWindowBitquadScalarComputer(*tree, "ramp max-tree 4-connectivity");
    }

    {
        auto tree = makeComponentTree(rampImage, false, 1.0);
        verifyGenericBitquadHistograms(*tree, "ramp min-tree 4-connectivity");
        verifyFiniteWindowBitquadScalarComputer(*tree, "ramp min-tree 4-connectivity");
    }

    {
        auto tree = makeTreeOfShapes(image, TestTopographicImmersion::SelfDualSpan);
        verifyGenericContourCounts(*tree, "tree of shapes");
        verifyPublicContourSideCountsComputer(*tree, "tree of shapes");
        verifyGenericBitquadHistograms(*tree, "tree of shapes");
    }

    {
        verifyTreeOfShapesScalarConnectivityPolicy(image);
        auto valuedTreeMin4Max8 = makeValuedTreeOfShapes(image, TestTopographicImmersion::Min4Max8);
        verifyGenericBitquadHistograms(valuedTreeMin4Max8->topology(), "tree of shapes Min4cMax8c");
        verifyFiniteWindowTreeOfShapesBitquadScalarComputer(*valuedTreeMin4Max8, "valued tree of shapes Min4cMax8c");
        auto valuedTreeMin8Max4 = makeValuedTreeOfShapes(image, TestTopographicImmersion::Min8Max4);
        verifyGenericBitquadHistograms(valuedTreeMin8Max4->topology(), "tree of shapes Min8cMax4c");
        verifyFiniteWindowTreeOfShapesBitquadScalarComputer(*valuedTreeMin8Max4, "valued tree of shapes Min8cMax4c");
    }

    {
        auto tree = makeTwoBranchTreeOfShapes();
        requireEqual(tree.root(), 2, "two-branch ToS root slot");
        verifyGenericContourCounts(tree, "two-branch tree of shapes");
        verifyPublicContourSideCountsComputer(tree, "two-branch tree of shapes");
        verifyGenericBitquadHistograms(tree, "two-branch tree of shapes");
    }

    return 0;
}
