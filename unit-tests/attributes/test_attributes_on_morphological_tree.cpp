#include "support/TestSupport.hpp"

#include "mmcfilters/attributes/AttributeComputedIncrementally.hpp"

#include <cmath>
#include <memory>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

int main() {
    auto image = makeComponentTreeFixture();

    for (bool isMaxtree : {true, false}) {
        auto weighted = makeWeightedComponentTree(image, isMaxtree);
        const auto& tree = weighted->topology();

        auto levelMapping = AttributeComputedIncrementally::computeAttributeMapping(*weighted, LEVEL);
        std::vector<float> expectedMapping;
        expectedMapping.reserve(static_cast<std::size_t>(image->getSize()));
        for (int p = 0; p < image->getSize(); ++p) {
            expectedMapping.push_back(static_cast<float>((*image)[p]));
        }
        requireVectorEqual(collectImageValues(levelMapping), expectedMapping, isMaxtree ? "max-tree level mapping" : "min-tree level mapping");

        auto [levelNames, levelBuffer] = AttributeComputedIncrementally::computeSingleAttribute(*weighted, LEVEL);
        for (NodeId nodeId : tree.getAliveNodeIds()) {
            requireEqual(
                levelBuffer[levelNames.linearIndex(nodeId, LEVEL)],
                static_cast<float>(weighted->getAltitude(nodeId)),
                isMaxtree ? "max-tree node level attribute" : "min-tree node level attribute"
            );
        }

        auto [grayNames, grayBuffer] = AttributeComputedIncrementally::computeSingleAttribute(*weighted, GRAY_HEIGHT);
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

        auto [areaNames, areaBuffer] = AttributeComputedIncrementally::computeSingleAttribute(*weighted, AREA);
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

        auto [maxDistNames, maxDistBuffer] = AttributeComputedIncrementally::computeSingleAttribute(*weighted, MAX_DIST);
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
            [&]() { (void)AttributeComputedIncrementally::computeSingleAttribute(*unweighted, MAX_DIST); },
            isMaxtree ? "max-tree MAX_DIST requires explicit altitude" : "min-tree MAX_DIST requires explicit altitude"
        );
    }

    {
        auto tree = makeComponentTree(image, true);

        auto requireAttributeValues = [&](Attribute attr, const std::vector<float>& expected, const std::string& label, float tolerance = 1e-5f) {
            auto [names, buffer] = AttributeComputedIncrementally::computeSingleAttribute(*tree, attr);
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

        requireAttributeValues(CENTRAL_MOMENT_11, {-0.28069448f, -2.5306945f, 1.1893053f, 1.88375f, 0.03999999f, 0.0f}, "CENTRAL_MOMENT_11");
        requireAttributeValues(CENTRAL_MOMENT_30, {1.7535391f, -1.6214609f, 3.9305391f, -3.2176094f, -0.36800003f, 0.0f}, "CENTRAL_MOMENT_30");
        requireAttributeValues(CENTRAL_MOMENT_03, {5.0719776f, 1.6969776f, -6.7110233f, -4.937875f, -1.0159998f, 0.0f}, "CENTRAL_MOMENT_03");
        requireAttributeValues(CENTRAL_MOMENT_21, {5.1162157f, 1.7412160f, -2.8267840f, -2.2619686f, -0.34399995f, 0.0f}, "CENTRAL_MOMENT_21");
        requireAttributeValues(CENTRAL_MOMENT_12, {-0.82678986f, -4.20179f, 0.91021085f, -2.1754375f, -0.15200004f, 0.0f}, "CENTRAL_MOMENT_12");

        requireAttributeValues(HU_MOMENT_1, {0.12161431f, 0.11837006f, 0.11856434f, 0.1317871f, 0.107999995f, 0.125f}, "HU_MOMENT_1");
        requireAttributeValues(HU_MOMENT_2, {1.2954243e-05f, 5.1657221e-04f, 2.9445661e-04f, 4.0753004e-03f, 1.4239997e-03f, 1.5625e-02f}, "HU_MOMENT_2");
        requireAttributeValues(HU_MOMENT_3, {1.1781303e-04f, 1.7525420e-04f, 1.8367014e-05f, 4.3831591e-04f, 2.5600054e-06f, 0.0f}, "HU_MOMENT_3");
        requireAttributeValues(HU_MOMENT_4, {9.9809789e-05f, 6.0222454e-05f, 4.5975851e-04f, 2.4695650e-03f, 6.7839975e-04f, 0.0f}, "HU_MOMENT_4");
        requireAttributeValues(HU_MOMENT_5, {-1.0747108e-08f, 2.0628408e-09f, -2.9037215e-08f, 2.5399299e-06f, 2.2419876e-08f, 0.0f}, "HU_MOMENT_5", 1e-8f);
        requireAttributeValues(HU_MOMENT_6, {2.4068953e-07f, 1.0917399e-06f, -7.3914539e-06f, 1.5666048e-04f, 2.0449266e-05f, 0.0f}, "HU_MOMENT_6");
        requireAttributeValues(HU_MOMENT_7, {1.2812862e-09f, -5.8328555e-09f, 3.0688721e-08f, -3.8774493e-07f, -1.7222872e-08f, 0.0f}, "HU_MOMENT_7", 1e-8f);

        requireAttributeValues(COMPACTNESS, {0.08179287f, 0.08963694f, 0.11186257f, 0.15095837f, 0.29473138f, 0.63661975f}, "COMPACTNESS");
        requireAttributeValues(ECCENTRICITY, {1.0609956f, 1.4752778f, 1.3384410f, 2.8789990f, 2.0741172f, 10.0f}, "ECCENTRICITY");
        requireAttributeValues(LENGTH_MAJOR_AXIS, {2.0017073f, 2.0574131f, 1.8048208f, 1.7691815f, 1.2072113f, 1.0f}, "LENGTH_MAJOR_AXIS");
        requireAttributeValues(LENGTH_MINOR_AXIS, {1.9433177f, 1.6938876f, 1.5600353f, 1.0426813f, 0.8382367f, 0.0f}, "LENGTH_MINOR_AXIS");
        requireAttributeValues(AXIS_ORIENTATION, {71.231224f, 49.107025f, 37.141315f, 56.380051f, 87.567741f, 89.999992f}, "AXIS_ORIENTATION", 1e-4f);

        requireAttributeValues(NUM_SIBLINGS_NODE, {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}, "NUM_SIBLINGS_NODE");
        requireAttributeValues(NUM_LEAF_DESCENDANTS_NODE, {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f}, "NUM_LEAF_DESCENDANTS_NODE");
        requireAttributeValues(LEAF_RATIO_NODE, {0.16666667f, 0.2f, 0.25f, 0.33333334f, 0.5f, 1.0f}, "LEAF_RATIO_NODE");
        requireAttributeValues(BALANCE_NODE, {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f}, "BALANCE_NODE");
        requireAttributeValues(AVG_CHILD_HEIGHT_NODE, {4.0f, 3.0f, 2.0f, 1.0f, 0.0f, 0.0f}, "AVG_CHILD_HEIGHT_NODE");

        requireAttributeValues(BITQUADS_NUMBER_EULER, {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}, "BITQUADS_NUMBER_EULER");
        requireAttributeValues(BITQUADS_NUMBER_HOLES, {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f}, "BITQUADS_NUMBER_HOLES");
        requireAttributeValues(BITQUADS_PERIMETER, {16.0f, 16.0f, 16.0f, 14.0f, 10.0f, 6.0f}, "BITQUADS_PERIMETER");
        requireAttributeValues(BITQUADS_PERIMETER_CONTINUOUS, {14.666667f, 14.0f, 12.666667f, 10.666667f, 8.0f, 4.6666665f}, "BITQUADS_PERIMETER_CONTINUOUS");
        requireAttributeValues(BITQUADS_CIRCULARITY, {0.93468875f, 0.96972632f, 0.96923792f, 0.92499042f, 1.0062914f, 1.1540544f}, "BITQUADS_CIRCULARITY");

        auto [bitquadsNames, bitquadsBuffer] = AttributeComputedIncrementally::computeAttributes(
            *tree,
            {BITQUADS_PERIMETER_AVERAGE, BITQUADS_LENGTH_AVERAGE, BITQUADS_WIDTH_AVERAGE}
        );
        for (NodeId nodeId : tree->getAliveNodeIds()) {
            require(std::isinf(bitquadsBuffer[bitquadsNames.linearIndex(nodeId, BITQUADS_PERIMETER_AVERAGE)]), "BITQUADS_PERIMETER_AVERAGE should be inf on fixture");
            require(std::isinf(bitquadsBuffer[bitquadsNames.linearIndex(nodeId, BITQUADS_LENGTH_AVERAGE)]), "BITQUADS_LENGTH_AVERAGE should be inf on fixture");
            require(std::isnan(bitquadsBuffer[bitquadsNames.linearIndex(nodeId, BITQUADS_WIDTH_AVERAGE)]), "BITQUADS_WIDTH_AVERAGE should be nan on fixture");
        }

        auto [groupHuNames, groupHuBuffer] = AttributeComputedIncrementally::computeSingleAttribute(*tree, AttributeGroup::HU_MOMENTS);
        requireNear(groupHuBuffer[groupHuNames.linearIndex(3, HU_MOMENT_4)], 0.0024695650f, 1e-6f, "HU group path");
        auto [groupTopoNames, groupTopoBuffer] = AttributeComputedIncrementally::computeSingleAttribute(*tree, AttributeGroup::TREE_TOPOLOGY);
        requireEqual(groupTopoBuffer[groupTopoNames.linearIndex(4, BALANCE_NODE)], 1.0f, "TREE_TOPOLOGY group path");
        auto [groupBoxNames, groupBoxBuffer] = AttributeComputedIncrementally::computeSingleAttribute(*tree, AttributeGroup::BOUNDING_BOX);
        requireEqual(groupBoxBuffer[groupBoxNames.linearIndex(5, BOX_COL_MIN)], 2.0f, "BOUNDING_BOX group path");
        auto [groupMomNames, groupMomBuffer] = AttributeComputedIncrementally::computeSingleAttribute(*tree, AttributeGroup::MOMENT_BASED);
        requireNear(groupMomBuffer[groupMomNames.linearIndex(5, ECCENTRICITY)], 10.0f, 1e-6f, "MOMENT_BASED group path");
        auto [groupBitquadsNames, groupBitquadsBuffer] = AttributeComputedIncrementally::computeSingleAttribute(*tree, AttributeGroup::BITQUADS);
        requireEqual(groupBitquadsBuffer[groupBitquadsNames.linearIndex(5, BITQUADS_PERIMETER)], 6.0f, "BITQUADS group path");

        auto weightedForDelta = makeWeightedComponentTree(image, true);
        auto [deltaNames, deltaBuffer] = AttributeComputedIncrementally::computeSingleAttributeWithDelta(*weightedForDelta, AREA, 2, "last-padding");
        requireEqual(deltaNames.NUM_ATTRIBUTES, 5, "delta attribute count");
        requireNear(deltaBuffer[deltaNames.linearIndex(0, AREA, -2)], 16.0f, 1e-6f, "delta asc2 node0");
        requireNear(deltaBuffer[deltaNames.linearIndex(3, AREA, -1)], 12.0f, 1e-6f, "delta asc1 node3");
        requireNear(deltaBuffer[deltaNames.linearIndex(3, AREA, 0)], 8.0f, 1e-6f, "delta center node3");
        requireNear(deltaBuffer[deltaNames.linearIndex(3, AREA, 1)], 5.0f, 1e-6f, "delta desc1 node3");
        requireNear(deltaBuffer[deltaNames.linearIndex(5, AREA, 2)], 2.0f, 1e-6f, "delta last-padding leaf");

        auto [deltaNullNames, deltaNullBuffer] = AttributeComputedIncrementally::computeSingleAttributeWithDelta(*weightedForDelta, AREA, 1, "null-padding");
        require(std::isnan(deltaNullBuffer[deltaNullNames.linearIndex(0, AREA, -1)]), "delta null-padding root missing asc");
        requireNear(deltaNullBuffer[deltaNullNames.linearIndex(0, AREA, 0)], 16.0f, 1e-6f, "delta null-padding root center");
        requireNear(deltaNullBuffer[deltaNullNames.linearIndex(0, AREA, 1)], 15.0f, 1e-6f, "delta null-padding root desc");
        requireNear(deltaNullBuffer[deltaNullNames.linearIndex(5, AREA, -1)], 5.0f, 1e-6f, "delta null-padding leaf asc");
        requireNear(deltaNullBuffer[deltaNullNames.linearIndex(5, AREA, 0)], 2.0f, 1e-6f, "delta null-padding leaf center");
        require(std::isnan(deltaNullBuffer[deltaNullNames.linearIndex(5, AREA, 1)]), "delta null-padding leaf missing desc");

        bool invalidPaddingRejected = false;
        try {
            (void)AttributeComputedIncrementally::computeSingleAttributeWithDelta(*weightedForDelta, AREA, 1, "unsupported-padding");
        } catch (const std::invalid_argument&) {
            invalidPaddingRejected = true;
        }
        require(invalidPaddingRejected, "delta invalid padding must throw");
    }

    {
        auto weighted = makeWeightedComponentTree(image, true);
        weighted->setAltitude(5, 7);

        requireEqual(weighted->getAltitude(5), 7, "weighted wrapper must read the external altitude buffer");

        auto [levelNames, levelBuffer] = AttributeComputedIncrementally::computeSingleAttribute(*weighted, LEVEL);
        requireEqual(levelBuffer[levelNames.linearIndex(5, LEVEL)], 7.0f, "weighted LEVEL must use external altitude buffer");

        auto [volumeNames, volumeBuffer] = AttributeComputedIncrementally::computeSingleAttribute(*weighted, VOLUME);
        requireEqual(volumeBuffer[volumeNames.linearIndex(5, VOLUME)], 14.0f, "weighted VOLUME leaf must use external altitude buffer");
        requireEqual(volumeBuffer[volumeNames.linearIndex(4, VOLUME)], 26.0f, "weighted VOLUME ancestor must aggregate external altitude buffer");

        auto levelMapping = AttributeComputedIncrementally::computeAttributeMapping(*weighted, LEVEL);
        auto levelValues = collectImageValues(levelMapping);
        requireEqual(levelValues[10], 7.0f, "weighted LEVEL mapping must project external altitude buffer on first leaf pixel");
        requireEqual(levelValues[14], 7.0f, "weighted LEVEL mapping must project external altitude buffer on second leaf pixel");

        auto [weightedParent, weightedAltitude] = weighted->exportHigraHierarchy();
        auto weightedRoundtrip = WeightedMorphologicalTree::createFromHigraParent(
            weightedParent,
            weightedAltitude,
            weighted->topology().getNumRowsOfImage(),
            weighted->topology().getNumColsOfImage(),
            MorphologicalTree::MAX_TREE,
            AdjacencyRelation(weighted->topology().getNumRowsOfImage(), weighted->topology().getNumColsOfImage(), 1.5));

        auto requireMappedAttributeMatch = [&](Attribute attr, const std::string& label) {
            auto weightedMapping = AttributeComputedIncrementally::computeAttributeMapping(*weighted, attr);
            auto roundtripMapping = AttributeComputedIncrementally::computeAttributeMapping(weightedRoundtrip, attr);
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
