#include "support/TestSupport.hpp"

#include "mmcfilters/attributes/AttributeComputation.hpp"
#include "mmcfilters/utils/Contract.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <numbers>
#include <vector>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

namespace {

struct DirectCentralMoments {
    double area = 0.0;
    double mu20 = 0.0;
    double mu02 = 0.0;
    double mu11 = 0.0;
    double mu30 = 0.0;
    double mu03 = 0.0;
    double mu21 = 0.0;
    double mu12 = 0.0;
};

constexpr double MAX_FINITE_ECCENTRICITY = 1.0e6;

DirectCentralMoments computeDirectCentralMoments(const MorphologicalTree& tree, NodeId nodeId) {
    struct Pixel {
        double x;
        double y;
    };

    std::vector<Pixel> support;
    const int numColumns = tree.numColumns();
    for (NodeId subtreeNodeId : tree.subtreeNodes(nodeId)) {
        for (PixelId pixel : tree.properPart(subtreeNodeId)) {
            const auto [row, column] = ImageUtils::to2D(pixel, numColumns);
            support.push_back({static_cast<double>(column), static_cast<double>(row)});
        }
    }
    require(!support.empty(), "direct moment oracle requires a non-empty support");

    DirectCentralMoments moments;
    moments.area = static_cast<double>(support.size());
    double sumX = 0.0;
    double sumY = 0.0;
    for (const Pixel& pixel : support) {
        sumX += pixel.x;
        sumY += pixel.y;
    }

    const double centroidX = sumX / moments.area;
    const double centroidY = sumY / moments.area;
    for (const Pixel& pixel : support) {
        const double dx = pixel.x - centroidX;
        const double dy = pixel.y - centroidY;
        moments.mu20 += dx * dx;
        moments.mu02 += dy * dy;
        moments.mu11 += dx * dy;
        moments.mu30 += dx * dx * dx;
        moments.mu03 += dy * dy * dy;
        moments.mu21 += dx * dx * dy;
        moments.mu12 += dx * dy * dy;
    }
    return moments;
}

double centralMomentValue(const DirectCentralMoments& moments, Attribute attribute) {
    switch (attribute) {
    case CentralMoment20:
        return moments.mu20;
    case CentralMoment02:
        return moments.mu02;
    case CentralMoment11:
        return moments.mu11;
    case CentralMoment30:
        return moments.mu30;
    case CentralMoment03:
        return moments.mu03;
    case CentralMoment21:
        return moments.mu21;
    case CentralMoment12:
        return moments.mu12;
    default:
        throw std::runtime_error("Unsupported central moment oracle attribute.");
    }
}

double huMomentValue(const DirectCentralMoments& moments, Attribute attribute) {
    auto normMoment = [&](double moment, int p, int q) { return moment / std::pow(moments.area, (p + q + 2.0) / 2.0); };

    const double eta20 = normMoment(moments.mu20, 2, 0);
    const double eta02 = normMoment(moments.mu02, 0, 2);
    const double eta11 = normMoment(moments.mu11, 1, 1);
    const double eta30 = normMoment(moments.mu30, 3, 0);
    const double eta03 = normMoment(moments.mu03, 0, 3);
    const double eta21 = normMoment(moments.mu21, 2, 1);
    const double eta12 = normMoment(moments.mu12, 1, 2);

    switch (attribute) {
    case HuMoment1:
        return eta20 + eta02;
    case HuMoment2:
        return std::pow(eta20 - eta02, 2.0) + 4.0 * std::pow(eta11, 2.0);
    case HuMoment3:
        return std::pow(eta30 - 3.0 * eta12, 2.0) + std::pow(3.0 * eta21 - eta03, 2.0);
    case HuMoment4:
        return std::pow(eta30 + eta12, 2.0) + std::pow(eta21 + eta03, 2.0);
    case HuMoment5:
        return (eta30 - 3.0 * eta12) * (eta30 + eta12) * (std::pow(eta30 + eta12, 2.0) - 3.0 * std::pow(eta21 + eta03, 2.0)) +
               (3.0 * eta21 - eta03) * (eta21 + eta03) * (3.0 * std::pow(eta30 + eta12, 2.0) - std::pow(eta21 + eta03, 2.0));
    case HuMoment6:
        return (eta20 - eta02) * (std::pow(eta30 + eta12, 2.0) - std::pow(eta21 + eta03, 2.0)) + 4.0 * eta11 * (eta30 + eta12) * (eta21 + eta03);
    case HuMoment7:
        return (3.0 * eta21 - eta03) * (eta30 + eta12) * (std::pow(eta30 + eta12, 2.0) - 3.0 * std::pow(eta21 + eta03, 2.0)) -
               (eta30 - 3.0 * eta12) * (eta21 + eta03) * (3.0 * std::pow(eta30 + eta12, 2.0) - std::pow(eta21 + eta03, 2.0));
    default:
        throw std::runtime_error("Unsupported Hu moment oracle attribute.");
    }
}

double momentDerivedValue(const DirectCentralMoments& moments, Attribute attribute) {
    double discriminant = std::pow(moments.mu20 - moments.mu02, 2.0) + 4.0 * std::pow(moments.mu11, 2.0);
    discriminant = std::max(discriminant, 0.0);
    const double lambda1 = moments.mu20 + moments.mu02 + std::sqrt(discriminant);
    const double lambda2 = moments.mu20 + moments.mu02 - std::sqrt(discriminant);

    switch (attribute) {
    case LengthMajorAxis:
        return moments.area > 0.0 && lambda1 > 0.0 ? std::sqrt((2.0 * lambda1) / moments.area) : 0.0;
    case LengthMinorAxis:
        return moments.area > 0.0 && lambda2 > 0.0 ? std::sqrt((2.0 * lambda2) / moments.area) : 0.0;
    case Eccentricity:
        if (lambda1 <= std::numeric_limits<float>::epsilon() && std::abs(lambda2) <= std::numeric_limits<float>::epsilon()) {
            return 1.0;
        }
        return lambda2 <= std::numeric_limits<float>::epsilon() ? MAX_FINITE_ECCENTRICITY : std::min(lambda1 / lambda2, MAX_FINITE_ECCENTRICITY);
    case Compactness: {
        const double denom = moments.mu20 + moments.mu02;
        return denom > std::numeric_limits<float>::epsilon() ? (1.0 / (2.0 * std::numbers::pi_v<double>)) * (moments.area / denom) : 0.0;
    }
    case AxisOrientation:
        if (moments.mu20 != moments.mu02 || moments.mu11 != 0.0) {
            const double radians = 0.5 * std::atan2(2.0 * moments.mu11, moments.mu20 - moments.mu02);
            const double degrees = radians * (180.0 / std::numbers::pi_v<double>);
            return std::fmod(std::abs(degrees), 360.0);
        }
        return 0.0;
    case Inertia:
        return moments.mu20 / std::pow(moments.area, 2.0) + moments.mu02 / std::pow(moments.area, 2.0);
    case Circularity:
        if (lambda1 <= std::numeric_limits<float>::epsilon() && std::abs(lambda2) <= std::numeric_limits<float>::epsilon()) {
            return 1.0;
        }
        return lambda1 <= std::numeric_limits<float>::epsilon() || lambda2 <= std::numeric_limits<float>::epsilon() ? 0.0 : lambda2 / lambda1;
    default:
        throw std::runtime_error("Unsupported moment-derived oracle attribute.");
    }
}

void requireAttributesMatchPixelOracle(const MorphologicalTree& tree, const std::vector<Attribute>& attributes,
                                       double (*oracle)(const DirectCentralMoments&, Attribute), const std::string& label, double tolerance = 1.0e-5) {
    std::vector<AttributeOrGroup> request;
    for (Attribute attribute : attributes) {
        request.emplace_back(attribute);
    }

    auto [names, buffer] = AttributeComputation::computeTopologyAttributes(tree, request);
    for (NodeId nodeId : tree.aliveNodeIds()) {
        const DirectCentralMoments moments = computeDirectCentralMoments(tree, nodeId);
        for (Attribute attribute : attributes) {
            const double actual = static_cast<double>(buffer[names.linearIndex(nodeId, attribute)]);
            const double expected = oracle(moments, attribute);
            if (std::isinf(expected)) {
                require(std::isinf(actual) && std::signbit(actual) == std::signbit(expected),
                        label + " " + AttributeNames::toString(attribute) + " node " + std::to_string(nodeId) + " expected infinity");
                continue;
            }
            if (std::isnan(expected)) {
                require(std::isnan(actual), label + " " + AttributeNames::toString(attribute) + " node " + std::to_string(nodeId) + " expected NaN");
                continue;
            }
            const double roundingTolerance = std::abs(expected) * 512.0 * static_cast<double>(std::numeric_limits<float>::epsilon());
            const double attrTolerance = std::max(attribute == AxisOrientation ? 1.0e-4 : tolerance, roundingTolerance);
            requireNear(actual, expected, attrTolerance, label + " " + AttributeNames::toString(attribute) + " node " + std::to_string(nodeId));
        }
    }
}

void requireMomentFamiliesMatchPixelOracle(const MorphologicalTree& tree, const std::string& label) {
    requireAttributesMatchPixelOracle(
        tree, {CentralMoment20, CentralMoment02, CentralMoment11, CentralMoment30, CentralMoment03, CentralMoment21, CentralMoment12},
        centralMomentValue, label + " central moment pixel oracle", 1.0e-4);
    requireAttributesMatchPixelOracle(tree, {HuMoment1, HuMoment2, HuMoment3, HuMoment4, HuMoment5, HuMoment6, HuMoment7}, huMomentValue,
                                      label + " Hu moment pixel oracle", 1.0e-6);
    requireAttributesMatchPixelOracle(tree, {Compactness, Eccentricity, LengthMajorAxis, LengthMinorAxis, AxisOrientation, Inertia, Circularity},
                                      momentDerivedValue, label + " moment-derived pixel oracle", 1.0e-4);
}

void requireMomentOutputDtypeOnlyAffectsStorage() {
    constexpr int columns = 8193;
    auto image = ImageUInt8::create(1, columns, static_cast<std::uint8_t>(7));
    auto valuedTree = MorphologicalTreeFactory::createMaxTree(image, 1.0);
    const MorphologicalTree& tree = valuedTree.topology();
    const NodeId root = tree.root();

    auto [floatNames, floatValues] = AttributeComputation::computeSingleTopologyAttribute(tree, CentralMoment20);
    auto [doubleNames, doubleValues] = AttributeComputation::computeSingleTopologyAttribute<double>(tree, CentralMoment20);

    const float floatValue = floatValues[floatNames.linearIndex(root, CentralMoment20)];
    const double doubleValue = doubleValues[doubleNames.linearIndex(root, CentralMoment20)];
    const long double n = static_cast<long double>(columns);
    const double expected = static_cast<double>(n * (n * n - 1.0L) / 12.0L);

    requireEqual(doubleValue, expected, "double CENTRAL_MOMENT_20 must preserve exact integer moment on large line support");
    requireEqual(floatValue, static_cast<float>(doubleValue), "float CENTRAL_MOMENT_20 must be the rounded double-facade output");
}

int directTopologyHeight(const MorphologicalTree& tree, NodeId nodeId) {
    if (tree.isLeaf(nodeId)) {
        return 0;
    }

    int height = 0;
    for (NodeId childId : tree.children(nodeId)) {
        height = std::max(height, directTopologyHeight(tree, childId) + 1);
    }
    return height;
}

int directChildHeightBalance(const MorphologicalTree& tree, NodeId nodeId) {
    if (tree.isLeaf(nodeId)) {
        return 0;
    }

    int minChildHeight = std::numeric_limits<int>::max();
    int maxChildHeight = 0;
    for (NodeId childId : tree.children(nodeId)) {
        const int childHeight = directTopologyHeight(tree, childId);
        minChildHeight = std::min(minChildHeight, childHeight);
        maxChildHeight = std::max(maxChildHeight, childHeight);
    }
    return maxChildHeight - minChildHeight;
}

void requireBalanceMatchesTopologyOracle(const MorphologicalTree& tree, const std::string& label) {
    auto [names, buffer] = AttributeComputation::computeSingleTopologyAttribute(tree, BalanceNode);
    for (NodeId nodeId : tree.aliveNodeIds()) {
        requireEqual(buffer[names.linearIndex(nodeId, BalanceNode)], static_cast<float>(directChildHeightBalance(tree, nodeId)),
                     label + " node " + std::to_string(nodeId));
    }
}

template <class ValuedTree> float canonicalGrayHeightOracle(const ValuedTree& valuedTree, NodeId nodeId) {
    const MorphologicalTree& tree = valuedTree.topology();
    const float nodeAltitude = static_cast<float>(valuedTree.nodeAltitude(nodeId));
    float maximumAbsoluteExcursion = 0.0f;
    for (NodeId subtreeNodeId : tree.subtreeNodes(nodeId)) {
        const float altitude = static_cast<float>(valuedTree.nodeAltitude(subtreeNodeId));
        maximumAbsoluteExcursion = std::max(maximumAbsoluteExcursion, std::abs(altitude - nodeAltitude));
    }
    return maximumAbsoluteExcursion;
}

template <class ValuedTree> void requireGrayHeightMatchesCanonicalSpanSemantics(const ValuedTree& valuedTree, const std::string& label) {
    const MorphologicalTree& tree = valuedTree.topology();
    auto [names, buffer] = AttributeComputation::computeSingleAttribute(valuedTree, GrayLevelHeight);
    for (NodeId nodeId : tree.aliveNodeIds()) {
        requireNear(buffer[names.linearIndex(nodeId, GrayLevelHeight)], canonicalGrayHeightOracle(valuedTree, nodeId), 1.0e-6f,
                    label + " canonical GrayLevelHeight span node " + std::to_string(nodeId));
    }
}

template <class ValuedTree> void requireGrayLevelStatisticsMatchSupportOracle(const ValuedTree& valuedTree, const std::string& label) {
    const MorphologicalTree& tree = valuedTree.topology();
    auto [names, buffer] = AttributeComputation::computeAttributes(valuedTree, {MeanGrayLevel, GrayLevelVariance});

    for (NodeId nodeId : tree.aliveNodeIds()) {
        std::vector<double> supportValues;
        for (PixelId pixel : tree.nodeSupport(nodeId)) {
            supportValues.push_back(static_cast<double>(valuedTree.nodeAltitude(tree.smallestNode(pixel))));
        }
        require(!supportValues.empty(), label + " gray-level support oracle requires a non-empty support");

        double sum = 0.0;
        for (double value : supportValues) {
            sum += value;
        }
        const double mean = sum / static_cast<double>(supportValues.size());

        double squaredDeviationSum = 0.0;
        for (double value : supportValues) {
            const double deviation = value - mean;
            squaredDeviationSum += deviation * deviation;
        }
        const double populationVariance = squaredDeviationSum / static_cast<double>(supportValues.size());

        requireNear(buffer[names.linearIndex(nodeId, MeanGrayLevel)], static_cast<float>(mean), 1.0e-5f,
                    label + " MEAN_GRAY_LEVEL support oracle node " + std::to_string(nodeId));
        requireNear(buffer[names.linearIndex(nodeId, GrayLevelVariance)], static_cast<float>(populationVariance), 1.0e-5f,
                    label + " GRAY_LEVEL_VARIANCE population oracle node " + std::to_string(nodeId));
    }
}

std::shared_ptr<MorphologicalTree> makeBranchingTopologyFixture() {
    constexpr int numPixels = 2;
    constexpr int numNodeSlots = 6;
    auto higraNodeId = [](NodeId slotId) { return numPixels + slotId; };

    std::vector<NodeId> parent(static_cast<std::size_t>(numPixels + numNodeSlots), InvalidNode);
    parent[0] = higraNodeId(2);
    parent[1] = higraNodeId(5);
    parent[static_cast<std::size_t>(higraNodeId(0))] = higraNodeId(0);
    parent[static_cast<std::size_t>(higraNodeId(1))] = higraNodeId(0);
    parent[static_cast<std::size_t>(higraNodeId(2))] = higraNodeId(1);
    parent[static_cast<std::size_t>(higraNodeId(3))] = higraNodeId(0);
    parent[static_cast<std::size_t>(higraNodeId(4))] = higraNodeId(3);
    parent[static_cast<std::size_t>(higraNodeId(5))] = higraNodeId(4);

    return makeTreeFromHigraParent(parent, 1, 2, true);
}

} // namespace

int main() {
    auto image = makeComponentTreeFixture();

    requireMomentOutputDtypeOnlyAffectsStorage();

    for (bool isMaxtree : {true, false}) {
        auto valuedTree = makeValuedComponentTree(image, isMaxtree);
        const auto& tree = valuedTree->topology();

        const auto reconstruction = valuedTree->reconstructFromNodeAltitudes();
        for (PixelId pixel = 0; pixel < image->getSize(); ++pixel) {
            requireEqual((*reconstruction)[pixel], (*image)[pixel],
                         isMaxtree ? "max-tree node-altitude reconstruction" : "min-tree node-altitude reconstruction");
        }

        auto [grayNames, grayBuffer] = AttributeComputation::computeSingleAttribute(*valuedTree, GrayLevelHeight);
        for (NodeId leafId : tree.aliveNodeIds()) {
            if (!tree.isLeaf(leafId)) {
                continue;
            }
            requireEqual(grayBuffer[grayNames.linearIndex(leafId, GrayLevelHeight)], 0.0f, isMaxtree ? "max-tree leaf gray height" : "min-tree leaf gray height");
        }
        requireGrayHeightMatchesCanonicalSpanSemantics(*valuedTree, isMaxtree ? "max-tree" : "min-tree");
        requireGrayLevelStatisticsMatchSupportOracle(*valuedTree, isMaxtree ? "max-tree" : "min-tree");

        auto [areaNames, areaBuffer] = AttributeComputation::computeSingleAttribute(*valuedTree, Area);
        requireEqual(areaBuffer[areaNames.linearIndex(tree.root(), Area)], 16.0f,
                     isMaxtree ? "max-tree root area attribute" : "min-tree root area attribute");
        for (NodeId nodeId : tree.aliveNodeIds()) {
            require(areaBuffer[areaNames.linearIndex(nodeId, Area)] >= static_cast<float>(tree.properPartCardinality(nodeId)),
                    isMaxtree ? "max-tree area must dominate direct proper parts" : "min-tree area must dominate direct proper parts");
        }
        if (isMaxtree) {
            requireEqual(areaBuffer[areaNames.linearIndex(0, Area)], 16.0f, "max-tree exact area root");
            requireEqual(areaBuffer[areaNames.linearIndex(1, Area)], 15.0f, "max-tree exact area node 1");
            requireEqual(areaBuffer[areaNames.linearIndex(2, Area)], 12.0f, "max-tree exact area node 2");
            requireEqual(areaBuffer[areaNames.linearIndex(3, Area)], 8.0f, "max-tree exact area node 3");
            requireEqual(areaBuffer[areaNames.linearIndex(4, Area)], 5.0f, "max-tree exact area node 4");
            requireEqual(areaBuffer[areaNames.linearIndex(5, Area)], 2.0f, "max-tree exact area node 5");
        }

        auto [maxDistNames, maxDistBuffer] = AttributeComputation::computeSingleAttribute(*valuedTree, MaxDist);
        const std::vector<float> expectedMaxDist =
            isMaxtree ? std::vector<float>{1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f} : std::vector<float>{1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f};
        for (NodeId nodeId : tree.aliveNodeIds()) {
            requireEqual(maxDistBuffer[maxDistNames.linearIndex(nodeId, MaxDist)], expectedMaxDist[static_cast<std::size_t>(nodeId)],
                         isMaxtree ? "max-tree MAX_DIST regression" : "min-tree MAX_DIST regression");
        }

        if constexpr (contract::validationsEnabled) {
            auto topologyOnly = makeComponentTree(image, isMaxtree);
            requireThrows<std::invalid_argument>([&]() { (void)AttributeComputation::computeSingleTopologyAttribute(*topologyOnly, MaxDist); },
                                                 isMaxtree ? "max-tree MAX_DIST requires explicit altitude" : "min-tree MAX_DIST requires explicit altitude");
        }

        std::vector<int> shiftedAltitude;
        shiftedAltitude.reserve(valuedTree->nodeAltitudes().size());
        const int offset = isMaxtree ? 300 : -300;
        for (std::uint8_t value : valuedTree->nodeAltitudes()) {
            shiftedAltitude.push_back(static_cast<int>(value) + offset);
        }
        const ValuedMorphologicalTreeView<int> shiftedView(valuedTree->topology(), std::span<const int>(shiftedAltitude));
        auto [shiftedMaxDistNames, shiftedMaxDistBuffer] = AttributeComputation::computeSingleAttribute(shiftedView, MaxDist);
        for (NodeId nodeId : tree.aliveNodeIds()) {
            requireEqual(shiftedMaxDistBuffer[shiftedMaxDistNames.linearIndex(nodeId, MaxDist)], maxDistBuffer[maxDistNames.linearIndex(nodeId, MaxDist)],
                         isMaxtree ? "max-tree MAX_DIST supports altitudes above 255" : "min-tree MAX_DIST supports negative altitudes");
        }

        requireMomentFamiliesMatchPixelOracle(tree, isMaxtree ? "max-tree" : "min-tree");
        requireBalanceMatchesTopologyOracle(tree, isMaxtree ? "max-tree balance" : "min-tree balance");
    }

    if constexpr (contract::validationsEnabled) {
        auto treeOfShapes = makeValuedTreeOfShapes(image, TestTopographicImmersion::SelfDualSpan);
        requireThrowsContaining<std::invalid_argument>([&]() { (void)AttributeComputation::computeSingleAttribute(*treeOfShapes, MaxDist); },
                                                       "globally monotone altitude order",
                                                       "MAX_DIST must reject a hierarchy without monotone altitude capability");
    }

    {
        auto edited = makeValuedComponentTree(image, true);
        edited->mergeNodeIntoParent(4);
        requireMomentFamiliesMatchPixelOracle(edited->topology(), "max-tree after mergeNodeIntoParent");
        requireBalanceMatchesTopologyOracle(edited->topology(), "max-tree after mergeNodeIntoParent balance");
    }

    {
        auto branching = makeBranchingTopologyFixture();
        requireBalanceMatchesTopologyOracle(*branching, "branching topology balance");

        auto [branchNames, branchBuffer] =
            AttributeComputation::computeTopologyAttributes(*branching, {SubtreeHeight, BalanceNode, AvgChildHeightNode, NumChildrenNode});
        requireEqual(branchBuffer[branchNames.linearIndex(0, NumChildrenNode)], 2.0f, "branching fixture root child count");
        requireEqual(branchBuffer[branchNames.linearIndex(0, SubtreeHeight)], 3.0f, "branching fixture root height");
        requireEqual(branchBuffer[branchNames.linearIndex(0, BalanceNode)], 1.0f, "branching fixture root balance");
        requireEqual(branchBuffer[branchNames.linearIndex(1, BalanceNode)], 0.0f, "branching fixture unary branch balance");
        requireEqual(branchBuffer[branchNames.linearIndex(3, BalanceNode)], 0.0f, "branching fixture deeper unary branch balance");
        requireNear(branchBuffer[branchNames.linearIndex(0, AvgChildHeightNode)], 1.5f, 1.0e-6f, "branching fixture average child height");
    }

    {
        auto tree = makeComponentTree(image, true);

        auto requireAttributeValues = [&](Attribute attr, const std::vector<float>& expected, const std::string& label, float tolerance = 1e-5f) {
            auto [names, buffer] = AttributeComputation::computeSingleTopologyAttribute(*tree, attr);
            for (NodeId nodeId : tree->aliveNodeIds()) {
                requireNear(buffer[names.linearIndex(nodeId, attr)], expected[static_cast<std::size_t>(nodeId)], tolerance,
                            label + " node " + std::to_string(nodeId));
            }
        };

        requireAttributeValues(RatioWh, {1.0f, 1.0f, 1.0f, 1.3333334f, 1.5f, 2.0f}, "RATIO_WH");
        requireAttributeValues(BoxColumnMin, {0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 2.0f}, "BOX_COLUMN_MIN");
        requireAttributeValues(BoxColumnMax, {3.0f, 3.0f, 3.0f, 2.0f, 2.0f, 2.0f}, "BOX_COLUMN_MAX");
        requireAttributeValues(BoxRowMin, {0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 2.0f}, "BOX_ROW_MIN");
        requireAttributeValues(BoxRowMax, {3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f}, "BOX_ROW_MAX");

        requireAttributeValues(NumSiblingsNode, {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}, "NUM_SIBLINGS_NODE");
        requireAttributeValues(NumLeafDescendantsNode, {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f}, "NUM_LEAF_DESCENDANTS_NODE");
        requireAttributeValues(LeafRatioNode, {0.16666667f, 0.2f, 0.25f, 0.33333334f, 0.5f, 1.0f}, "LEAF_RATIO_NODE");
        requireAttributeValues(BalanceNode, {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}, "BALANCE_NODE");
        requireAttributeValues(AvgChildHeightNode, {4.0f, 3.0f, 2.0f, 1.0f, 0.0f, 0.0f}, "AVG_CHILD_HEIGHT_NODE");

        requireAttributeValues(BitquadNumberEuler, {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}, "BITQUAD_NUMBER_EULER");
        requireAttributeValues(BitquadNumberHoles, {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f}, "BITQUAD_NUMBER_HOLES");
        requireAttributeValues(BitquadPerimeter, {16.0f, 16.0f, 16.0f, 14.0f, 10.0f, 6.0f}, "BITQUAD_PERIMETER");
        requireAttributeValues(BitquadPerimeterContinuous, {14.666667f, 14.0f, 12.666667f, 10.666667f, 8.0f, 4.6666665f}, "BITQUAD_PERIMETER_CONTINUOUS");
        requireAttributeValues(BitquadCircularity, {0.93468875f, 0.96972632f, 0.96923792f, 0.92499042f, 1.0062914f, 1.1540544f}, "BITQUAD_CIRCULARITY");

        auto [bitquadsNames, bitquadsBuffer] =
            AttributeComputation::computeTopologyAttributes(*tree, {BitquadPerimeterAverage, BitquadLengthAverage, BitquadWidthAverage});
        for (NodeId nodeId : tree->aliveNodeIds()) {
            requireEqual(bitquadsBuffer[bitquadsNames.linearIndex(nodeId, BitquadPerimeterAverage)], 0.0f,
                         "BITQUAD_PERIMETER_AVERAGE finite zero fallback");
            requireEqual(bitquadsBuffer[bitquadsNames.linearIndex(nodeId, BitquadLengthAverage)], 0.0f, "BITQUAD_LENGTH_AVERAGE finite zero fallback");
            require(std::isfinite(bitquadsBuffer[bitquadsNames.linearIndex(nodeId, BitquadWidthAverage)]), "BITQUAD_WIDTH_AVERAGE finite fallback");
        }

        auto [groupMomentNames, groupMomentBuffer] = AttributeComputation::computeSingleTopologyAttribute(*tree, AttributeGroup::Moments);
        requireNear(static_cast<double>(groupMomentBuffer[groupMomentNames.linearIndex(3, HuMoment4)]),
                    huMomentValue(computeDirectCentralMoments(*tree, 3), HuMoment4), 1.0e-6, "MOMENTS group Hu path");
        auto [groupTopoNames, groupTopoBuffer] = AttributeComputation::computeSingleTopologyAttribute(*tree, AttributeGroup::TreeTopology);
        requireEqual(groupTopoBuffer[groupTopoNames.linearIndex(4, BalanceNode)], 0.0f, "TREE_TOPOLOGY group path");
        require(groupMomentNames.contains(Circularity), "MOMENTS group includes moment-derived descriptors");
        const double groupEccentricity = static_cast<double>(groupMomentBuffer[groupMomentNames.linearIndex(5, Eccentricity)]);
        const double expectedGroupEccentricity = momentDerivedValue(computeDirectCentralMoments(*tree, 5), Eccentricity);
        requireNear(groupEccentricity, expectedGroupEccentricity, 1.0e-5, "MOMENTS group moment-derived path");
        auto [groupBoundaryNames, groupBoundaryBuffer] = AttributeComputation::computeSingleTopologyAttribute(*tree, AttributeGroup::Boundary);
        requireEqual(groupBoundaryBuffer[groupBoundaryNames.linearIndex(5, BitquadPerimeter)], 6.0f, "BOUNDARY group Bitquad path");
        requireEqual(groupBoundaryNames.NUM_ATTRIBUTES, 15, "BOUNDARY group count");
        require(groupBoundaryNames.contains(ContourPixels), "BOUNDARY group includes contour pixels");
        require(groupBoundaryNames.contains(ContourSideSouth), "BOUNDARY group includes south sides");
        require(groupBoundaryBuffer.size() ==
                    static_cast<std::size_t>(tree->numInternalNodeSlots()) * static_cast<std::size_t>(groupBoundaryNames.NUM_ATTRIBUTES),
                "BOUNDARY group path");

        auto valuedTreeForSampling = makeValuedComponentTree(image, true);
        auto [groupShapeNames, groupShapeBuffer] = AttributeComputation::computeSingleAttribute(*valuedTreeForSampling, AttributeGroup::Shape);
        requireEqual(groupShapeNames.NUM_ATTRIBUTES, 47, "SHAPE group count");
        requireEqual(groupShapeBuffer[groupShapeNames.linearIndex(5, BoxColumnMin)], 2.0f, "SHAPE group bounding-box path");
        require(groupShapeNames.contains(MaxDist), "SHAPE group includes MAX_DIST");
        require(groupShapeNames.contains(ContourSideSouth), "SHAPE group includes contour sides");
        auto [sampleLayout, sampledValues] =
            AttributeComputation::computeSampledNodeAttribute(*valuedTreeForSampling, Area, AltitudeDifference<std::uint8_t>{1}, 2);
        requireEqual(sampleLayout.NUM_ATTRIBUTES, 5, "sampled attribute count");
        requireNear(sampledValues[sampleLayout.linearIndex(0, Area, -2)], 16.0f, 1e-6f, "root ancestor-2 sample");
        requireNear(sampledValues[sampleLayout.linearIndex(3, Area, -1)], 12.0f, 1e-6f, "node 3 ancestor-1 sample");
        requireNear(sampledValues[sampleLayout.linearIndex(3, Area, 0)], 8.0f, 1e-6f, "node 3 current sample");
        requireNear(sampledValues[sampleLayout.linearIndex(3, Area, 1)], 5.0f, 1e-6f, "node 3 descendant-1 sample");
        requireNear(sampledValues[sampleLayout.linearIndex(5, Area, 2)], 2.0f, 1e-6f, "leaf repeated descendant-2 sample");

        auto [centerOnlyLayout, centerOnlyValues] =
            AttributeComputation::computeSampledNodeAttribute(*valuedTreeForSampling, Area, AltitudeDifference<std::uint8_t>{1}, 0);
        requireEqual(centerOnlyLayout.NUM_ATTRIBUTES, 1, "center-only sampled attribute count");
        requireNear(centerOnlyValues[centerOnlyLayout.linearIndex(3, Area, 0)], 8.0f, 1e-6f, "center-only node sample");

        auto [nanSampleLayout, nanSampleValues] =
            AttributeComputation::computeSampledNodeAttribute(
                *valuedTreeForSampling, Area, AltitudeDifference<std::uint8_t>{1}, 1,
                NodeAttributeSamplingPolicy::LargestSupportDescendant, MissingNodeAttributeSamplePolicy::NotANumber);
        require(std::isnan(nanSampleValues[nanSampleLayout.linearIndex(0, Area, -1)]), "root missing ancestor sample");
        requireNear(nanSampleValues[nanSampleLayout.linearIndex(0, Area, 0)], 16.0f, 1e-6f, "root current sample");
        requireNear(nanSampleValues[nanSampleLayout.linearIndex(0, Area, 1)], 15.0f, 1e-6f, "root descendant sample");
        requireNear(nanSampleValues[nanSampleLayout.linearIndex(5, Area, -1)], 5.0f, 1e-6f, "leaf ancestor sample");
        requireNear(nanSampleValues[nanSampleLayout.linearIndex(5, Area, 0)], 2.0f, 1e-6f, "leaf current sample");
        require(std::isnan(nanSampleValues[nanSampleLayout.linearIndex(5, Area, 1)]), "leaf missing descendant sample");

        auto [repeatNearestLayout, repeatNearestValues] = AttributeComputation::computeSampledNodeAttribute(
            *valuedTreeForSampling, Area, AltitudeDifference<std::uint8_t>{100}, 2,
            NodeAttributeSamplingPolicy::LargestSupportDescendant, MissingNodeAttributeSamplePolicy::RepeatNearest);
        requireNear(repeatNearestValues[repeatNearestLayout.linearIndex(3, Area, -1)], 8.0f, 1e-6f,
                    "unavailable ancestor must repeat the current-node sample");
        requireNear(repeatNearestValues[repeatNearestLayout.linearIndex(3, Area, -2)], 8.0f, 1e-6f,
                    "consecutive unavailable ancestors must repeat the nearest sample");
        requireNear(repeatNearestValues[repeatNearestLayout.linearIndex(3, Area, 1)], 8.0f, 1e-6f,
                    "unavailable descendant must repeat the current-node sample");

        auto [notANumberLayout, notANumberValues] = AttributeComputation::computeSampledNodeAttribute(
            *valuedTreeForSampling, Area, AltitudeDifference<std::uint8_t>{100}, 1,
            NodeAttributeSamplingPolicy::LargestSupportDescendant, MissingNodeAttributeSamplePolicy::NotANumber);
        require(std::isnan(notANumberValues[notANumberLayout.linearIndex(3, Area, -1)]),
                "unavailable ancestor must materialize as NaN");
        require(std::isnan(notANumberValues[notANumberLayout.linearIndex(3, Area, 1)]),
                "unavailable descendant must materialize as NaN");

        auto [zeroLayout, zeroValues] = AttributeComputation::computeSampledNodeAttribute(
            *valuedTreeForSampling, Area, AltitudeDifference<std::uint8_t>{100}, 1,
            NodeAttributeSamplingPolicy::LargestSupportDescendant, MissingNodeAttributeSamplePolicy::Zero);
        requireNear(zeroValues[zeroLayout.linearIndex(3, Area, -1)], 0.0f, 1e-6f,
                    "unavailable ancestor must materialize as zero");
        requireNear(zeroValues[zeroLayout.linearIndex(3, Area, 1)], 0.0f, 1e-6f,
                    "unavailable descendant must materialize as zero");

        if constexpr (contract::validationsEnabled) {
            bool invalidMissingSamplePolicyRejected = false;
            try {
                (void)AttributeComputation::computeSampledNodeAttribute(
                    *valuedTreeForSampling, Area, AltitudeDifference<std::uint8_t>{1}, 1,
                    NodeAttributeSamplingPolicy::LargestSupportDescendant, static_cast<MissingNodeAttributeSamplePolicy>(999));
            } catch (const std::invalid_argument&) {
                invalidMissingSamplePolicyRejected = true;
            }
            require(invalidMissingSamplePolicyRejected, "invalid missing-sample policy must throw");

            requireThrows<std::invalid_argument>(
                [&]() {
                    static_cast<void>(
                        AttributeComputation::computeSampledNodeAttribute(*valuedTreeForSampling, Area, AltitudeDifference<std::uint8_t>{1}, -1));
                },
                "negative sampling radius must throw");
        }
    }

    {
        const std::vector<NodeId> parent{5, 5, 4, 4, 6, 6, 6};
        const std::vector<std::uint8_t> altitude{2, 2, 2, 2, 2, 2, 0};
        auto branchingTree = MorphologicalTreeFactory::createFromHigraParent(
            std::span<const NodeId>(parent), std::span<const std::uint8_t>(altitude), 1, 4, MorphologicalTreeKind::MaxTree,
            RegularGridAdjacency2D(1, 4, 1.0));

        auto [sampleLayout, sampledValues] = AttributeComputation::computeSampledNodeAttribute(
            branchingTree, BoxColumnMin, AltitudeDifference<std::uint8_t>{2}, 1,
            NodeAttributeSamplingPolicy::LargestSupportDescendant, MissingNodeAttributeSamplePolicy::NotANumber);
        const NodeId root = branchingTree.topology().root();
        requireNear(sampledValues[sampleLayout.linearIndex(root, BoxColumnMin, 1)], 0.0f, 1e-6f,
                    "equal-support descendant tie must select the support with the smallest row-major pixel");
    }

    {
        auto valuedTree = makeValuedComponentTree(image, true);
        valuedTree->setNodeAltitude(5, 7);

        requireEqual(valuedTree->nodeAltitude(5), 7, "valuedTree wrapper must read the external altitude buffer");

        auto [volumeNames, volumeBuffer] = AttributeComputation::computeSingleAttribute(*valuedTree, Volume);
        requireEqual(volumeBuffer[volumeNames.linearIndex(5, Volume)], 14.0f, "valuedTree VOLUME leaf must use external altitude buffer");
        requireEqual(volumeBuffer[volumeNames.linearIndex(4, Volume)], 26.0f, "valuedTree VOLUME ancestor must aggregate external altitude buffer");

        const auto reconstructed = valuedTree->reconstructFromNodeAltitudes();
        requireEqual((*reconstructed)[10], std::uint8_t{7}, "node-altitude reconstruction must project the external altitude on the first leaf pixel");
        requireEqual((*reconstructed)[14], std::uint8_t{7}, "node-altitude reconstruction must project the external altitude on the second leaf pixel");

        auto [valuedTreeParent, valuedTreeAltitude] = valuedTree->exportHigraHierarchy();
        auto valuedTreeRoundtrip = MorphologicalTreeFactory::createFromHigraParent(
            std::span<const NodeId>(valuedTreeParent), std::span<const std::uint8_t>(valuedTreeAltitude), valuedTree->topology().numRows(),
            valuedTree->topology().numColumns(), MorphologicalTreeKind::MaxTree,
            RegularGridAdjacency2D(valuedTree->topology().numRows(), valuedTree->topology().numColumns(), 1.5));

        auto requireMappedAttributeMatch = [&](Attribute attr, const std::string& label) {
            auto valuedTreeMapping = AttributeComputation::computeAttributeMapping(*valuedTree, attr);
            auto roundtripMapping = AttributeComputation::computeAttributeMapping(valuedTreeRoundtrip, attr);
            requireVectorEqual(collectImageValues(valuedTreeMapping), collectImageValues(roundtripMapping), label);
        };

        requireMappedAttributeMatch(MeanGrayLevel, "valuedTree MeanGrayLevel mapping must match equivalent valuedTree round-trip hierarchy");
        requireMappedAttributeMatch(GrayLevelVariance, "valuedTree GrayLevelVariance mapping must match equivalent valuedTree round-trip hierarchy");
        requireMappedAttributeMatch(GrayLevelHeight, "valuedTree GrayLevelHeight mapping must match equivalent valuedTree round-trip hierarchy");
        requireMappedAttributeMatch(Volume, "valuedTree VOLUME mapping must match equivalent valuedTree round-trip hierarchy");
        requireMappedAttributeMatch(RelativeVolume, "valuedTree RELATIVE_VOLUME mapping must match equivalent valuedTree round-trip hierarchy");
        requireMappedAttributeMatch(MaxDist, "valuedTree MAX_DIST mapping must match equivalent valuedTree round-trip hierarchy");
    }

    return 0;
}
