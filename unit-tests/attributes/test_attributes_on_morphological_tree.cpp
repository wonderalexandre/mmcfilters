#include "support/TestSupport.hpp"

#include "mmcfilters/attributes/AttributeComputation.hpp"

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

DirectCentralMoments computeDirectCentralMoments(const MorphologicalTree& tree, NodeId nodeId) {
    struct Pixel {
        double x;
        double y;
    };

    std::vector<Pixel> support;
    const int numCols = tree.getNumColsOfImage();
    for (NodeId subtreeNodeId : tree.getNodeSubtree(nodeId)) {
        for (int properPart : tree.getProperParts(subtreeNodeId)) {
            const auto [row, col] = ImageUtils::to2D(properPart, numCols);
            support.push_back({static_cast<double>(col), static_cast<double>(row)});
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
        case CENTRAL_MOMENT_20: return moments.mu20;
        case CENTRAL_MOMENT_02: return moments.mu02;
        case CENTRAL_MOMENT_11: return moments.mu11;
        case CENTRAL_MOMENT_30: return moments.mu30;
        case CENTRAL_MOMENT_03: return moments.mu03;
        case CENTRAL_MOMENT_21: return moments.mu21;
        case CENTRAL_MOMENT_12: return moments.mu12;
        default: throw std::runtime_error("Unsupported central moment oracle attribute.");
    }
}

double huMomentValue(const DirectCentralMoments& moments, Attribute attribute) {
    auto normMoment = [&](double moment, int p, int q) {
        return moment / std::pow(moments.area, (p + q + 2.0) / 2.0);
    };

    const double eta20 = normMoment(moments.mu20, 2, 0);
    const double eta02 = normMoment(moments.mu02, 0, 2);
    const double eta11 = normMoment(moments.mu11, 1, 1);
    const double eta30 = normMoment(moments.mu30, 3, 0);
    const double eta03 = normMoment(moments.mu03, 0, 3);
    const double eta21 = normMoment(moments.mu21, 2, 1);
    const double eta12 = normMoment(moments.mu12, 1, 2);

    switch (attribute) {
        case HU_MOMENT_1:
            return eta20 + eta02;
        case HU_MOMENT_2:
            return std::pow(eta20 - eta02, 2.0) + 4.0 * std::pow(eta11, 2.0);
        case HU_MOMENT_3:
            return std::pow(eta30 - 3.0 * eta12, 2.0) + std::pow(3.0 * eta21 - eta03, 2.0);
        case HU_MOMENT_4:
            return std::pow(eta30 + eta12, 2.0) + std::pow(eta21 + eta03, 2.0);
        case HU_MOMENT_5:
            return (eta30 - 3.0 * eta12) * (eta30 + eta12) *
                       (std::pow(eta30 + eta12, 2.0) - 3.0 * std::pow(eta21 + eta03, 2.0)) +
                   (3.0 * eta21 - eta03) * (eta21 + eta03) *
                       (3.0 * std::pow(eta30 + eta12, 2.0) - std::pow(eta21 + eta03, 2.0));
        case HU_MOMENT_6:
            return (eta20 - eta02) * (std::pow(eta30 + eta12, 2.0) - std::pow(eta21 + eta03, 2.0)) +
                   4.0 * eta11 * (eta30 + eta12) * (eta21 + eta03);
        case HU_MOMENT_7:
            return (3.0 * eta21 - eta03) * (eta30 + eta12) *
                       (std::pow(eta30 + eta12, 2.0) - 3.0 * std::pow(eta21 + eta03, 2.0)) -
                   (eta30 - 3.0 * eta12) * (eta21 + eta03) *
                       (3.0 * std::pow(eta30 + eta12, 2.0) - std::pow(eta21 + eta03, 2.0));
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
        case LENGTH_MAJOR_AXIS:
            return moments.area > 0.0 && lambda1 > 0.0 ? std::sqrt((2.0 * lambda1) / moments.area) : 0.0;
        case LENGTH_MINOR_AXIS:
            return moments.area > 0.0 && lambda2 > 0.0 ? std::sqrt((2.0 * lambda2) / moments.area) : 0.0;
        case ECCENTRICITY:
            if (lambda1 <= std::numeric_limits<float>::epsilon() &&
                std::abs(lambda2) <= std::numeric_limits<float>::epsilon()) {
                return 1.0;
            }
            return lambda2 <= std::numeric_limits<float>::epsilon()
                ? std::numeric_limits<double>::infinity()
                : lambda1 / lambda2;
        case COMPACTNESS: {
            const double denom = moments.mu20 + moments.mu02;
            return denom > std::numeric_limits<float>::epsilon()
                ? (1.0 / (2.0 * std::numbers::pi_v<double>)) * (moments.area / denom)
                : 0.0;
        }
        case AXIS_ORIENTATION:
            if (moments.mu20 != moments.mu02 || moments.mu11 != 0.0) {
                const double radians = 0.5 * std::atan2(2.0 * moments.mu11, moments.mu20 - moments.mu02);
                const double degrees = radians * (180.0 / std::numbers::pi_v<double>);
                return std::fmod(std::abs(degrees), 360.0);
            }
            return 0.0;
        case INERTIA:
            return moments.mu20 / std::pow(moments.area, 2.0) +
                   moments.mu02 / std::pow(moments.area, 2.0);
        case CIRCULARITY:
            if (lambda1 <= std::numeric_limits<float>::epsilon() &&
                std::abs(lambda2) <= std::numeric_limits<float>::epsilon()) {
                return 1.0;
            }
            return lambda1 <= std::numeric_limits<float>::epsilon() ||
                   lambda2 <= std::numeric_limits<float>::epsilon()
                ? 0.0
                : lambda2 / lambda1;
        default:
            throw std::runtime_error("Unsupported moment-derived oracle attribute.");
    }
}

void requireAttributesMatchPixelOracle(
    const MorphologicalTree& tree,
    const std::vector<Attribute>& attributes,
    double (*oracle)(const DirectCentralMoments&, Attribute),
    const std::string& label,
    double tolerance = 1.0e-5)
{
    std::vector<AttributeOrGroup> request;
    for (Attribute attribute : attributes) {
        request.emplace_back(attribute);
    }

    auto [names, buffer] = AttributeComputation::computeTopologyAttributes(tree, request);
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        const DirectCentralMoments moments = computeDirectCentralMoments(tree, nodeId);
        for (Attribute attribute : attributes) {
            const double actual = static_cast<double>(buffer[names.linearIndex(nodeId, attribute)]);
            const double expected = oracle(moments, attribute);
            if (std::isinf(expected)) {
                require(
                    std::isinf(actual) && std::signbit(actual) == std::signbit(expected),
                    label + " " + AttributeNames::toString(attribute) + " node " + std::to_string(nodeId) + " expected infinity");
                continue;
            }
            if (std::isnan(expected)) {
                require(
                    std::isnan(actual),
                    label + " " + AttributeNames::toString(attribute) + " node " + std::to_string(nodeId) + " expected NaN");
                continue;
            }
            const double attrTolerance = attribute == AXIS_ORIENTATION ? 1.0e-4 : tolerance;
            requireNear(
                actual,
                expected,
                attrTolerance,
                label + " " + AttributeNames::toString(attribute) + " node " + std::to_string(nodeId));
        }
    }
}

void requireMomentFamiliesMatchPixelOracle(const MorphologicalTree& tree, const std::string& label) {
    requireAttributesMatchPixelOracle(
        tree,
        {CENTRAL_MOMENT_20, CENTRAL_MOMENT_02, CENTRAL_MOMENT_11, CENTRAL_MOMENT_30, CENTRAL_MOMENT_03, CENTRAL_MOMENT_21, CENTRAL_MOMENT_12},
        centralMomentValue,
        label + " central moment pixel oracle");
    requireAttributesMatchPixelOracle(
        tree,
        {HU_MOMENT_1, HU_MOMENT_2, HU_MOMENT_3, HU_MOMENT_4, HU_MOMENT_5, HU_MOMENT_6, HU_MOMENT_7},
        huMomentValue,
        label + " Hu moment pixel oracle",
        1.0e-6);
    requireAttributesMatchPixelOracle(
        tree,
        {COMPACTNESS, ECCENTRICITY, LENGTH_MAJOR_AXIS, LENGTH_MINOR_AXIS, AXIS_ORIENTATION, INERTIA, CIRCULARITY},
        momentDerivedValue,
        label + " moment-derived pixel oracle");
}

int directTopologyHeight(const MorphologicalTree& tree, NodeId nodeId) {
    if (tree.isLeaf(nodeId)) {
        return 0;
    }

    int height = 0;
    for (NodeId childId : tree.getChildren(nodeId)) {
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
    for (NodeId childId : tree.getChildren(nodeId)) {
        const int childHeight = directTopologyHeight(tree, childId);
        minChildHeight = std::min(minChildHeight, childHeight);
        maxChildHeight = std::max(maxChildHeight, childHeight);
    }
    return maxChildHeight - minChildHeight;
}

void requireBalanceMatchesTopologyOracle(const MorphologicalTree& tree, const std::string& label) {
    auto [names, buffer] = AttributeComputation::computeSingleTopologyAttribute(tree, BALANCE_NODE);
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        requireEqual(
            buffer[names.linearIndex(nodeId, BALANCE_NODE)],
            static_cast<float>(directChildHeightBalance(tree, nodeId)),
            label + " node " + std::to_string(nodeId));
    }
}

template<class WeightedTree>
float canonicalGrayHeightOracle(const WeightedTree& weighted, NodeId nodeId) {
    const MorphologicalTree& tree = weighted.topology();
    if (tree.isLeaf(nodeId)) {
        return 0.0f;
    }

    float extreme = static_cast<float>(weighted.getAltitude(nodeId));
    for (NodeId subtreeNodeId : tree.getNodeSubtree(nodeId)) {
        const float altitude = static_cast<float>(weighted.getAltitude(subtreeNodeId));
        if (tree.getTreeType() == MorphologicalTreeKind::MAX_TREE) {
            extreme = std::max(extreme, altitude);
        } else {
            extreme = std::min(extreme, altitude);
        }
    }

    return std::abs(static_cast<float>(weighted.getAltitude(nodeId)) - extreme);
}

template<class WeightedTree>
void requireGrayHeightMatchesCanonicalSpanSemantics(
    const WeightedTree& weighted,
    const std::string& label)
{
    const MorphologicalTree& tree = weighted.topology();
    auto [names, buffer] = AttributeComputation::computeSingleAttribute(weighted, GRAY_HEIGHT);
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        requireNear(
            buffer[names.linearIndex(nodeId, GRAY_HEIGHT)],
            canonicalGrayHeightOracle(weighted, nodeId),
            1.0e-6f,
            label + " canonical GRAY_HEIGHT span node " + std::to_string(nodeId));
    }
}

std::shared_ptr<MorphologicalTree> makeBranchingTopologyFixture() {
    constexpr NodeId numProperParts = 1;
    constexpr NodeId numNodeSlots = 6;
    auto higraNodeId = [](NodeId slotId) {
        return numProperParts + slotId;
    };

    std::vector<NodeId> parent(static_cast<std::size_t>(numProperParts + numNodeSlots), InvalidNode);
    parent[0] = higraNodeId(2);
    parent[static_cast<std::size_t>(higraNodeId(0))] = higraNodeId(0);
    parent[static_cast<std::size_t>(higraNodeId(1))] = higraNodeId(0);
    parent[static_cast<std::size_t>(higraNodeId(2))] = higraNodeId(1);
    parent[static_cast<std::size_t>(higraNodeId(3))] = higraNodeId(0);
    parent[static_cast<std::size_t>(higraNodeId(4))] = higraNodeId(3);
    parent[static_cast<std::size_t>(higraNodeId(5))] = higraNodeId(4);

    return makeTreeFromHigraParent(parent, 1, 1, true);
}

} // namespace

int main() {
    auto image = makeComponentTreeFixture();

    for (bool isMaxtree : {true, false}) {
        auto weighted = makeWeightedComponentTree(image, isMaxtree);
        const auto& tree = weighted->topology();

        auto levelMapping = AttributeComputation::computeAttributeMapping(*weighted, LEVEL);
        std::vector<float> expectedMapping;
        expectedMapping.reserve(static_cast<std::size_t>(image->getSize()));
        for (int p = 0; p < image->getSize(); ++p) {
            expectedMapping.push_back(static_cast<float>((*image)[p]));
        }
        requireVectorEqual(collectImageValues(levelMapping), expectedMapping, isMaxtree ? "max-tree level mapping" : "min-tree level mapping");

        auto [levelNames, levelBuffer] = AttributeComputation::computeSingleAttribute(*weighted, LEVEL);
        for (NodeId nodeId : tree.getAliveNodeIds()) {
            requireEqual(
                levelBuffer[levelNames.linearIndex(nodeId, LEVEL)],
                static_cast<float>(weighted->getAltitude(nodeId)),
                isMaxtree ? "max-tree node level attribute" : "min-tree node level attribute"
            );
        }

        auto [grayNames, grayBuffer] = AttributeComputation::computeSingleAttribute(*weighted, GRAY_HEIGHT);
        for (NodeId leafId : tree.getAliveNodeIds()) {
            if (!tree.isLeaf(leafId)) {
                continue;
            }
            requireEqual(
                grayBuffer[grayNames.linearIndex(leafId, GRAY_HEIGHT)],
                0.0f,
                isMaxtree ? "max-tree leaf gray height" : "min-tree leaf gray height"
            );
        }
        requireGrayHeightMatchesCanonicalSpanSemantics(
            *weighted,
            isMaxtree ? "max-tree" : "min-tree");

        auto [areaNames, areaBuffer] = AttributeComputation::computeSingleAttribute(*weighted, AREA);
        requireEqual(
            areaBuffer[areaNames.linearIndex(tree.getRoot(), AREA)],
            16.0f,
            isMaxtree ? "max-tree root area attribute" : "min-tree root area attribute"
        );
        for (NodeId nodeId : tree.getAliveNodeIds()) {
            require(
                areaBuffer[areaNames.linearIndex(nodeId, AREA)] >= static_cast<float>(tree.getNumProperParts(nodeId)),
                isMaxtree ? "max-tree area must dominate direct proper parts" : "min-tree area must dominate direct proper parts"
            );
        }
        if (isMaxtree) {
            requireEqual(areaBuffer[areaNames.linearIndex(0, AREA)], 16.0f, "max-tree exact area root");
            requireEqual(areaBuffer[areaNames.linearIndex(1, AREA)], 15.0f, "max-tree exact area node 1");
            requireEqual(areaBuffer[areaNames.linearIndex(2, AREA)], 12.0f, "max-tree exact area node 2");
            requireEqual(areaBuffer[areaNames.linearIndex(3, AREA)], 8.0f, "max-tree exact area node 3");
            requireEqual(areaBuffer[areaNames.linearIndex(4, AREA)], 5.0f, "max-tree exact area node 4");
            requireEqual(areaBuffer[areaNames.linearIndex(5, AREA)], 2.0f, "max-tree exact area node 5");
        }

        auto [maxDistNames, maxDistBuffer] = AttributeComputation::computeSingleAttribute(*weighted, MAX_DIST);
        const std::vector<float> expectedMaxDist = isMaxtree
            ? std::vector<float>{1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f}
            : std::vector<float>{1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f};
        for (NodeId nodeId : tree.getAliveNodeIds()) {
            requireEqual(
                maxDistBuffer[maxDistNames.linearIndex(nodeId, MAX_DIST)],
                expectedMaxDist[static_cast<std::size_t>(nodeId)],
                isMaxtree ? "max-tree MAX_DIST regression" : "min-tree MAX_DIST regression"
            );
        }

        auto unweighted = makeComponentTree(image, isMaxtree);
        requireThrows<std::invalid_argument>(
            [&]() { (void)AttributeComputation::computeSingleTopologyAttribute(*unweighted, MAX_DIST); },
            isMaxtree ? "max-tree MAX_DIST requires explicit altitude" : "min-tree MAX_DIST requires explicit altitude"
        );

        std::vector<int> shiftedAltitude;
        shiftedAltitude.reserve(weighted->getAltitudeBuffer().size());
        const int offset = isMaxtree ? 300 : -300;
        for (std::uint8_t value : weighted->getAltitudeBuffer()) {
            shiftedAltitude.push_back(static_cast<int>(value) + offset);
        }
        const WeightedTreeView<int> shiftedView(
            weighted->topology(),
            std::span<const int>(shiftedAltitude));
        auto [shiftedMaxDistNames, shiftedMaxDistBuffer] =
            AttributeComputation::computeSingleAttribute(shiftedView, MAX_DIST);
        for (NodeId nodeId : tree.getAliveNodeIds()) {
            requireEqual(
                shiftedMaxDistBuffer[shiftedMaxDistNames.linearIndex(nodeId, MAX_DIST)],
                maxDistBuffer[maxDistNames.linearIndex(nodeId, MAX_DIST)],
                isMaxtree
                    ? "max-tree MAX_DIST supports altitudes above 255"
                    : "min-tree MAX_DIST supports negative altitudes"
            );
        }

        requireMomentFamiliesMatchPixelOracle(
            tree,
            isMaxtree ? "max-tree" : "min-tree");
        requireBalanceMatchesTopologyOracle(
            tree,
            isMaxtree ? "max-tree balance" : "min-tree balance");
    }

    {
        auto treeOfShapes = makeWeightedTreeOfShapes(image, ToSInterpolation::SelfDual);
        requireThrowsContaining<std::invalid_argument>(
            [&]() { (void)AttributeComputation::computeSingleAttribute(*treeOfShapes, MAX_DIST); },
            "TREE_OF_SHAPES",
            "MAX_DIST must explicitly reject tree of shapes");

        auto selfDualResidualTree = MorphologicalTreeFactory::createSelfDualResidualTree(image, 1.5);
        requireThrowsContaining<std::invalid_argument>(
            [&]() { (void)AttributeComputation::computeSingleAttribute(selfDualResidualTree, MAX_DIST); },
            "SELF_DUAL_RESIDUAL_TREE",
            "MAX_DIST must explicitly reject self-dual residual trees");
    }

    {
        auto edited = makeWeightedComponentTree(image, true);
        edited->mergeNodeIntoParent(4);
        requireMomentFamiliesMatchPixelOracle(
            edited->topology(),
            "max-tree after mergeNodeIntoParent");
        requireBalanceMatchesTopologyOracle(
            edited->topology(),
            "max-tree after mergeNodeIntoParent balance");
    }

    {
        auto branching = makeBranchingTopologyFixture();
        requireBalanceMatchesTopologyOracle(*branching, "branching topology balance");

        auto [branchNames, branchBuffer] = AttributeComputation::computeTopologyAttributes(
            *branching,
            {HEIGHT_NODE, BALANCE_NODE, AVG_CHILD_HEIGHT_NODE, NUM_CHILDREN_NODE});
        requireEqual(branchBuffer[branchNames.linearIndex(0, NUM_CHILDREN_NODE)], 2.0f, "branching fixture root child count");
        requireEqual(branchBuffer[branchNames.linearIndex(0, HEIGHT_NODE)], 3.0f, "branching fixture root height");
        requireEqual(branchBuffer[branchNames.linearIndex(0, BALANCE_NODE)], 1.0f, "branching fixture root balance");
        requireEqual(branchBuffer[branchNames.linearIndex(1, BALANCE_NODE)], 0.0f, "branching fixture unary branch balance");
        requireEqual(branchBuffer[branchNames.linearIndex(3, BALANCE_NODE)], 0.0f, "branching fixture deeper unary branch balance");
        requireNear(branchBuffer[branchNames.linearIndex(0, AVG_CHILD_HEIGHT_NODE)], 1.5f, 1.0e-6f, "branching fixture average child height");
    }

    {
        auto tree = makeComponentTree(image, true);

        auto requireAttributeValues = [&](Attribute attr, const std::vector<float>& expected, const std::string& label, float tolerance = 1e-5f) {
            auto [names, buffer] = AttributeComputation::computeSingleTopologyAttribute(*tree, attr);
            for (NodeId nodeId : tree->getAliveNodeIds()) {
                requireNear(
                    buffer[names.linearIndex(nodeId, attr)],
                    expected[static_cast<std::size_t>(nodeId)],
                    tolerance,
                    label + " node " + std::to_string(nodeId)
                );
            }
        };

        requireAttributeValues(RATIO_WH, {1.0f, 1.0f, 1.0f, 1.3333334f, 1.5f, 2.0f}, "RATIO_WH");
        requireAttributeValues(BOX_COL_MIN, {0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 2.0f}, "BOX_COL_MIN");
        requireAttributeValues(BOX_COL_MAX, {3.0f, 3.0f, 3.0f, 2.0f, 2.0f, 2.0f}, "BOX_COL_MAX");
        requireAttributeValues(BOX_ROW_MIN, {0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 2.0f}, "BOX_ROW_MIN");
        requireAttributeValues(BOX_ROW_MAX, {3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f}, "BOX_ROW_MAX");

        requireAttributeValues(NUM_SIBLINGS_NODE, {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}, "NUM_SIBLINGS_NODE");
        requireAttributeValues(NUM_LEAF_DESCENDANTS_NODE, {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f}, "NUM_LEAF_DESCENDANTS_NODE");
        requireAttributeValues(LEAF_RATIO_NODE, {0.16666667f, 0.2f, 0.25f, 0.33333334f, 0.5f, 1.0f}, "LEAF_RATIO_NODE");
        requireAttributeValues(BALANCE_NODE, {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}, "BALANCE_NODE");
        requireAttributeValues(AVG_CHILD_HEIGHT_NODE, {4.0f, 3.0f, 2.0f, 1.0f, 0.0f, 0.0f}, "AVG_CHILD_HEIGHT_NODE");

        requireAttributeValues(BITQUADS_NUMBER_EULER, {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}, "BITQUADS_NUMBER_EULER");
        requireAttributeValues(BITQUADS_NUMBER_HOLES, {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f}, "BITQUADS_NUMBER_HOLES");
        requireAttributeValues(BITQUADS_PERIMETER, {16.0f, 16.0f, 16.0f, 14.0f, 10.0f, 6.0f}, "BITQUADS_PERIMETER");
        requireAttributeValues(BITQUADS_PERIMETER_CONTINUOUS, {14.666667f, 14.0f, 12.666667f, 10.666667f, 8.0f, 4.6666665f}, "BITQUADS_PERIMETER_CONTINUOUS");
        requireAttributeValues(BITQUADS_CIRCULARITY, {0.93468875f, 0.96972632f, 0.96923792f, 0.92499042f, 1.0062914f, 1.1540544f}, "BITQUADS_CIRCULARITY");

        auto [bitquadsNames, bitquadsBuffer] = AttributeComputation::computeTopologyAttributes(
            *tree,
            {BITQUADS_PERIMETER_AVERAGE, BITQUADS_LENGTH_AVERAGE, BITQUADS_WIDTH_AVERAGE}
        );
        for (NodeId nodeId : tree->getAliveNodeIds()) {
            require(std::isinf(bitquadsBuffer[bitquadsNames.linearIndex(nodeId, BITQUADS_PERIMETER_AVERAGE)]), "BITQUADS_PERIMETER_AVERAGE should be inf on fixture");
            require(std::isinf(bitquadsBuffer[bitquadsNames.linearIndex(nodeId, BITQUADS_LENGTH_AVERAGE)]), "BITQUADS_LENGTH_AVERAGE should be inf on fixture");
            require(std::isnan(bitquadsBuffer[bitquadsNames.linearIndex(nodeId, BITQUADS_WIDTH_AVERAGE)]), "BITQUADS_WIDTH_AVERAGE should be nan on fixture");
        }

        auto [groupMomentNames, groupMomentBuffer] = AttributeComputation::computeSingleTopologyAttribute(*tree, AttributeGroup::MOMENTS);
        requireNear(
            static_cast<double>(groupMomentBuffer[groupMomentNames.linearIndex(3, HU_MOMENT_4)]),
            huMomentValue(computeDirectCentralMoments(*tree, 3), HU_MOMENT_4),
            1.0e-6,
            "MOMENTS group Hu path");
        auto [groupTopoNames, groupTopoBuffer] = AttributeComputation::computeSingleTopologyAttribute(*tree, AttributeGroup::TREE_TOPOLOGY);
        requireEqual(groupTopoBuffer[groupTopoNames.linearIndex(4, BALANCE_NODE)], 0.0f, "TREE_TOPOLOGY group path");
        require(groupMomentNames.contains(CIRCULARITY), "MOMENTS group includes moment-derived descriptors");
        const double groupEccentricity = static_cast<double>(groupMomentBuffer[groupMomentNames.linearIndex(5, ECCENTRICITY)]);
        const double expectedGroupEccentricity = momentDerivedValue(computeDirectCentralMoments(*tree, 5), ECCENTRICITY);
        require(
            std::isinf(groupEccentricity) && std::isinf(expectedGroupEccentricity),
            "MOMENTS group moment-derived path");
        auto [groupBoundaryNames, groupBoundaryBuffer] = AttributeComputation::computeSingleTopologyAttribute(*tree, AttributeGroup::BOUNDARY);
        requireEqual(groupBoundaryBuffer[groupBoundaryNames.linearIndex(5, BITQUADS_PERIMETER)], 6.0f, "BOUNDARY group BitQuads path");
        requireEqual(groupBoundaryNames.NUM_ATTRIBUTES, 15, "BOUNDARY group count");
        require(groupBoundaryNames.contains(CONTOUR_PIXELS), "BOUNDARY group includes contour pixels");
        require(groupBoundaryNames.contains(CONTOUR_SIDE_SOUTH), "BOUNDARY group includes south sides");
        require(
            groupBoundaryBuffer.size() ==
                static_cast<std::size_t>(tree->getNumInternalNodeSlots()) *
                static_cast<std::size_t>(groupBoundaryNames.NUM_ATTRIBUTES),
            "BOUNDARY group path");

        auto weightedForDelta = makeWeightedComponentTree(image, true);
        auto [groupShapeNames, groupShapeBuffer] = AttributeComputation::computeSingleAttribute(*weightedForDelta, AttributeGroup::SHAPE);
        requireEqual(groupShapeNames.NUM_ATTRIBUTES, 47, "SHAPE group count");
        requireEqual(groupShapeBuffer[groupShapeNames.linearIndex(5, BOX_COL_MIN)], 2.0f, "SHAPE group bounding-box path");
        require(groupShapeNames.contains(MAX_DIST), "SHAPE group includes MAX_DIST");
        require(groupShapeNames.contains(CONTOUR_SIDE_SOUTH), "SHAPE group includes contour sides");
        auto [deltaNames, deltaBuffer] = AttributeComputation::computeSingleAttributeWithDelta(*weightedForDelta, AREA, AltitudeDiff<std::uint8_t>{1}, 2, "last-padding");
        requireEqual(deltaNames.NUM_ATTRIBUTES, 5, "delta attribute count");
        requireNear(deltaBuffer[deltaNames.linearIndex(0, AREA, -2)], 16.0f, 1e-6f, "delta asc2 node0");
        requireNear(deltaBuffer[deltaNames.linearIndex(3, AREA, -1)], 12.0f, 1e-6f, "delta asc1 node3");
        requireNear(deltaBuffer[deltaNames.linearIndex(3, AREA, 0)], 8.0f, 1e-6f, "delta center node3");
        requireNear(deltaBuffer[deltaNames.linearIndex(3, AREA, 1)], 5.0f, 1e-6f, "delta desc1 node3");
        requireNear(deltaBuffer[deltaNames.linearIndex(5, AREA, 2)], 2.0f, 1e-6f, "delta last-padding leaf");

        auto [zeroDeltaNames, zeroDeltaBuffer] = AttributeComputation::computeSingleAttributeWithDelta(*weightedForDelta, AREA, AltitudeDiff<std::uint8_t>{1}, 0, "last-padding");
        requireEqual(zeroDeltaNames.NUM_ATTRIBUTES, 1, "zero delta topology-only attribute count");
        requireNear(zeroDeltaBuffer[zeroDeltaNames.linearIndex(3, AREA, 0)], 8.0f, 1e-6f, "zero delta topology-only center");

        auto [deltaNullNames, deltaNullBuffer] = AttributeComputation::computeSingleAttributeWithDelta(*weightedForDelta, AREA, AltitudeDiff<std::uint8_t>{1}, 1, "null-padding");
        require(std::isnan(deltaNullBuffer[deltaNullNames.linearIndex(0, AREA, -1)]), "delta null-padding root missing asc");
        requireNear(deltaNullBuffer[deltaNullNames.linearIndex(0, AREA, 0)], 16.0f, 1e-6f, "delta null-padding root center");
        requireNear(deltaNullBuffer[deltaNullNames.linearIndex(0, AREA, 1)], 15.0f, 1e-6f, "delta null-padding root desc");
        requireNear(deltaNullBuffer[deltaNullNames.linearIndex(5, AREA, -1)], 5.0f, 1e-6f, "delta null-padding leaf asc");
        requireNear(deltaNullBuffer[deltaNullNames.linearIndex(5, AREA, 0)], 2.0f, 1e-6f, "delta null-padding leaf center");
        require(std::isnan(deltaNullBuffer[deltaNullNames.linearIndex(5, AREA, 1)]), "delta null-padding leaf missing desc");

        bool invalidPaddingRejected = false;
        try {
            (void)AttributeComputation::computeSingleAttributeWithDelta(*weightedForDelta, AREA, AltitudeDiff<std::uint8_t>{1}, 1, "unsupported-padding");
        } catch (const std::invalid_argument&) {
            invalidPaddingRejected = true;
        }
        require(invalidPaddingRejected, "delta invalid padding must throw");

        requireThrows<std::invalid_argument>(
            [&]() {
                static_cast<void>(AttributeComputation::computeSingleAttributeWithDelta(
                    *weightedForDelta,
                    AREA,
                    AltitudeDiff<std::uint8_t>{1},
                    -1,
                    "last-padding"));
            },
            "delta negative radius must throw");
    }

    {
        auto weighted = makeWeightedComponentTree(image, true);
        weighted->setAltitude(5, 7);

        requireEqual(weighted->getAltitude(5), 7, "weighted wrapper must read the external altitude buffer");

        auto [levelNames, levelBuffer] = AttributeComputation::computeSingleAttribute(*weighted, LEVEL);
        requireEqual(levelBuffer[levelNames.linearIndex(5, LEVEL)], 7.0f, "weighted LEVEL must use external altitude buffer");

        auto [volumeNames, volumeBuffer] = AttributeComputation::computeSingleAttribute(*weighted, VOLUME);
        requireEqual(volumeBuffer[volumeNames.linearIndex(5, VOLUME)], 14.0f, "weighted VOLUME leaf must use external altitude buffer");
        requireEqual(volumeBuffer[volumeNames.linearIndex(4, VOLUME)], 26.0f, "weighted VOLUME ancestor must aggregate external altitude buffer");

        auto levelMapping = AttributeComputation::computeAttributeMapping(*weighted, LEVEL);
        auto levelValues = collectImageValues(levelMapping);
        requireEqual(levelValues[10], 7.0f, "weighted LEVEL mapping must project external altitude buffer on first leaf pixel");
        requireEqual(levelValues[14], 7.0f, "weighted LEVEL mapping must project external altitude buffer on second leaf pixel");

        auto [weightedParent, weightedAltitude] = weighted->exportHigraHierarchy();
        auto weightedRoundtrip = MorphologicalTreeFactory::createFromHigraParent(
            std::span<const NodeId>(weightedParent),
            std::span<const std::uint8_t>(weightedAltitude),
            weighted->topology().getNumRowsOfImage(),
            weighted->topology().getNumColsOfImage(),
            MorphologicalTreeKind::MAX_TREE,
            AdjacencyRelation(weighted->topology().getNumRowsOfImage(), weighted->topology().getNumColsOfImage(), 1.5));

        auto requireMappedAttributeMatch = [&](Attribute attr, const std::string& label) {
            auto weightedMapping = AttributeComputation::computeAttributeMapping(*weighted, attr);
            auto roundtripMapping = AttributeComputation::computeAttributeMapping(weightedRoundtrip, attr);
            requireVectorEqual(
                collectImageValues(weightedMapping),
                collectImageValues(roundtripMapping),
                label);
        };

        requireMappedAttributeMatch(LEVEL, "weighted LEVEL mapping must match equivalent weighted round-trip hierarchy");
        requireMappedAttributeMatch(MEAN_LEVEL, "weighted MEAN_LEVEL mapping must match equivalent weighted round-trip hierarchy");
        requireMappedAttributeMatch(VARIANCE_LEVEL, "weighted VARIANCE_LEVEL mapping must match equivalent weighted round-trip hierarchy");
        requireMappedAttributeMatch(GRAY_HEIGHT, "weighted GRAY_HEIGHT mapping must match equivalent weighted round-trip hierarchy");
        requireMappedAttributeMatch(VOLUME, "weighted VOLUME mapping must match equivalent weighted round-trip hierarchy");
        requireMappedAttributeMatch(RELATIVE_VOLUME, "weighted RELATIVE_VOLUME mapping must match equivalent weighted round-trip hierarchy");
        requireMappedAttributeMatch(MAX_DIST, "weighted MAX_DIST mapping must match equivalent weighted round-trip hierarchy");
    }

    return 0;
}
