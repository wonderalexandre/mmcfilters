#include "support/TestSupport.hpp"

#include <cmath>
#include <memory>
#include <stdexcept>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

void requireCompactHigraHierarchy(
    const std::vector<NodeId>& parent,
    const std::vector<AltitudeType>& altitude,
    int numProperParts,
    const std::string& label) {
    requireEqual(parent.size(), altitude.size(), label + " parent/altitude size agreement");

    const NodeId totalNodes = static_cast<NodeId>(parent.size());
    int numRoots = 0;

    for (NodeId properPartId = 0; properPartId < numProperParts; ++properPartId) {
        const NodeId ownerNodeId = parent[static_cast<size_t>(properPartId)];
        require(ownerNodeId >= numProperParts && ownerNodeId < totalNodes, label + " leaf parent must be an internal Higra node");
        requireEqual(
            altitude[static_cast<size_t>(properPartId)],
            altitude[static_cast<size_t>(ownerNodeId)],
            label + " leaf altitude must match its owning internal node");
    }

    for (NodeId nodeId = numProperParts; nodeId < totalNodes; ++nodeId) {
        const NodeId parentNodeId = parent[static_cast<size_t>(nodeId)];
        require(parentNodeId >= numProperParts && parentNodeId < totalNodes, label + " internal parent must stay in the internal Higra domain");
        if (parentNodeId == nodeId) {
            ++numRoots;
        } else {
            require(parentNodeId > nodeId, label + " internal parent must appear after its child");
        }
    }

    requireEqual(numRoots, 1, label + " must encode exactly one self-parented root");
}

void requireFloatVectorEqualAllowingNaN(
    const std::vector<float>& actual,
    const std::vector<float>& expected,
    const std::string& label) {
    requireEqual(actual.size(), expected.size(), label + " size");
    for (size_t i = 0; i < actual.size(); ++i) {
        if (std::isnan(actual[i]) && std::isnan(expected[i])) {
            continue;
        }
        requireNear(actual[i], expected[i], 1.0e-6f, label + " value");
    }
}

AttributeNames makeDenseAttributeNames(const std::vector<Attribute>& attributes) {
    std::unordered_map<Attribute, int> offsets;
    for (int index = 0; index < static_cast<int>(attributes.size()); ++index) {
        offsets[attributes[static_cast<size_t>(index)]] = index;
    }
    return AttributeNames(std::move(offsets));
}

int main() {
    {
        auto image = makeComponentTreeFixture();
        auto weighted = makeWeightedComponentTree(image, true);

        weighted->validateAltitudeBufferShape();
        weighted->validateMonotoneAltitude();

        requireEqual(
            static_cast<int>(weighted->getAltitudeBuffer().size()),
            weighted->topology().getNumInternalNodeSlots(),
            "weighted altitude buffer must match the dense internal-node domain");
        NodeId builtSampleNodeId = InvalidNode;
        AltitudeType builtSampleAltitude = AltitudeType{};
        for (NodeId nodeId : weighted->topology().getAliveNodeIds()) {
            if (weighted->getAltitude(nodeId) != 0) {
                builtSampleNodeId = nodeId;
                builtSampleAltitude = weighted->getAltitude(nodeId);
                break;
            }
        }
        require(builtSampleNodeId != InvalidNode, "weighted image build must expose at least one non-zero internal altitude");
        requireEqual(weighted->getAltitude(builtSampleNodeId), builtSampleAltitude, "weighted image build external altitude sample");
        requireEqual(weighted->getAltitude(weighted->topology().getRoot()), 0, "weighted root altitude");
        requireEqual(weighted->getNodeResidue(5), 1, "weighted node residue");

        auto reconstruction = weighted->reconstructionImage();
        requireImageShape(reconstruction, 4, 4);
        requireVectorEqual(
            collectImageValues(reconstruction),
            collectImageValues(image),
            "weighted reconstruction must match the original image");

        const auto [higraParent, higraAltitude] = weighted->exportHigraHierarchy();
        requireEqual(
            static_cast<int>(higraParent.size()),
            weighted->topology().getNumTotalProperParts() + weighted->topology().getNumNodes(),
            "weighted Higra export parent size");
        requireEqual(
            static_cast<int>(higraAltitude.size()),
            weighted->topology().getNumTotalProperParts() + weighted->topology().getNumNodes(),
            "weighted Higra export altitude size");
        requireCompactHigraHierarchy(
            higraParent,
            higraAltitude,
            weighted->topology().getNumTotalProperParts(),
            "weighted Higra export");

        auto roundtrip = WeightedMorphologicalTree::createFromHigraParent(
            higraParent,
            higraAltitude,
            4,
            4,
            MorphologicalTree::MAX_TREE,
            AdjacencyRelation(4, 4, 1.5));
        roundtrip.validateAltitudeBufferShape();
        roundtrip.validateMonotoneAltitude();
        NodeId importedSampleNodeId = InvalidNode;
        AltitudeType importedSampleAltitude = AltitudeType{};
        for (NodeId nodeId : roundtrip.topology().getAliveNodeIds()) {
            if (roundtrip.getAltitude(nodeId) != 0) {
                importedSampleNodeId = nodeId;
                importedSampleAltitude = roundtrip.getAltitude(nodeId);
                break;
            }
        }
        require(importedSampleNodeId != InvalidNode, "weighted Higra import must expose at least one non-zero internal altitude");
        requireVectorEqual(
            collectImageValues(roundtrip.reconstructionImage()),
            collectImageValues(image),
            "weighted Higra round-trip reconstruction");
        const auto [reexportedParent, reexportedAltitude] = roundtrip.exportHigraHierarchy();
        requireVectorEqual(reexportedParent, higraParent, "weighted Higra parent round-trip");
        requireVectorEqual(reexportedAltitude, higraAltitude, "weighted Higra altitude round-trip");

        auto [maxDistNames, maxDistByNode] = AttributeComputedIncrementally::computeSingleAttribute(*weighted, MAX_DIST);
        auto maxDistByExportedHigra = AttributeComputedIncrementally::projectNodeValuesToExportedHigra(
            *weighted,
            maxDistNames,
            maxDistByNode);
        auto [maxDistHigraNames, maxDistByImportedHigra] = AttributeComputedIncrementally::computeSingleAttribute(
            roundtrip,
            MAX_DIST,
            {},
            NodeIdSpace::HIGRA);
        for (NodeId properPart = 0; properPart < weighted->topology().getNumTotalProperParts(); ++properPart) {
            maxDistByImportedHigra[static_cast<size_t>(properPart)] = 0.0f;
        }
        requireFloatVectorEqualAllowingNaN(
            maxDistByExportedHigra,
            maxDistByImportedHigra,
            "weighted exported-Higra MAX_DIST projection");

        auto [areaNames, areaByNode] = AttributeComputedIncrementally::computeSingleAttribute(*weighted, AREA);
        std::vector<float> areaAndMaxDist(static_cast<size_t>(weighted->topology().getNumInternalNodeSlots()) * 2, 0.0f);
        for (NodeId nodeId = 0; nodeId < weighted->topology().getNumInternalNodeSlots(); ++nodeId) {
            areaAndMaxDist[static_cast<size_t>(nodeId) * 2] = areaByNode[areaNames.linearIndex(nodeId, AREA)];
            areaAndMaxDist[static_cast<size_t>(nodeId) * 2 + 1] = maxDistByNode[maxDistNames.linearIndex(nodeId, MAX_DIST)];
        }
        AttributeNames areaAndMaxDistNames = makeDenseAttributeNames({AREA, MAX_DIST});
        const auto projectedPair = AttributeComputedIncrementally::projectNodeValuesToExportedHigra(
            *weighted,
            areaAndMaxDistNames,
            areaAndMaxDist);
        requireEqual(projectedPair.size(), maxDistByExportedHigra.size() * 2, "2D exported-Higra projection size");
        requireEqual(projectedPair[0], 1.0f, "2D exported-Higra unit AREA value");
        requireEqual(projectedPair[1], 0.0f, "2D exported-Higra unit MAX_DIST value");

        const std::vector<Attribute> unitProjectionAttributes{LEVEL, VOLUME, BOX_COL_MIN, BOX_ROW_MIN};
        const AttributeNames unitProjectionNames = makeDenseAttributeNames(unitProjectionAttributes);
        std::vector<float> unitProjectionNodeValues(
            static_cast<size_t>(weighted->topology().getNumInternalNodeSlots()) *
            static_cast<size_t>(unitProjectionNames.NUM_ATTRIBUTES),
            0.0f);
        for (Attribute attribute : unitProjectionAttributes) {
            auto [singleNames, singleValues] = AttributeComputedIncrementally::computeSingleAttribute(*weighted, attribute);
            for (NodeId nodeId = 0; nodeId < weighted->topology().getNumInternalNodeSlots(); ++nodeId) {
                unitProjectionNodeValues[unitProjectionNames.linearIndex(nodeId, attribute)] =
                    singleValues[singleNames.linearIndex(nodeId, attribute)];
            }
        }

        const auto projectedUnitValues = AttributeComputedIncrementally::projectNodeValuesToExportedHigra(
            *weighted,
            unitProjectionNames,
            unitProjectionNodeValues);
        const NodeId sampleProperPart = 10;
        const NodeId sampleOwner = weighted->topology().getSmallestComponent(sampleProperPart);
        const auto [sampleRow, sampleCol] = ImageUtils::to2D(sampleProperPart, weighted->topology().getNumColsOfImage());
        requireEqual(
            projectedUnitValues[unitProjectionNames.linearIndex(sampleProperPart, LEVEL)],
            static_cast<float>(weighted->getAltitude(sampleOwner)),
            "exported-Higra unit LEVEL must use the proper-part altitude");
        requireEqual(
            projectedUnitValues[unitProjectionNames.linearIndex(sampleProperPart, VOLUME)],
            static_cast<float>(weighted->getAltitude(sampleOwner)),
            "exported-Higra unit VOLUME must use one pixel at the proper-part altitude");
        requireEqual(
            projectedUnitValues[unitProjectionNames.linearIndex(sampleProperPart, BOX_COL_MIN)],
            static_cast<float>(sampleCol),
            "exported-Higra unit BOX_COL_MIN must use the proper-part column");
        requireEqual(
            projectedUnitValues[unitProjectionNames.linearIndex(sampleProperPart, BOX_ROW_MIN)],
            static_cast<float>(sampleRow),
            "exported-Higra unit BOX_ROW_MIN must use the proper-part row");

        requireThrows<std::invalid_argument>(
            [&]() {
                const std::vector<float> invalidValues{1.0f};
                static_cast<void>(AttributeComputedIncrementally::projectNodeValuesToExportedHigra(*weighted, maxDistNames, invalidValues));
            },
            "exported-Higra projection must reject invalid node-value size");

        roundtrip.setAltitude(importedSampleNodeId, static_cast<AltitudeType>(importedSampleAltitude + 4));
        auto reimported = WeightedMorphologicalTree::createFromHigraParent(
            higraParent,
            higraAltitude,
            4,
            4,
            MorphologicalTree::MAX_TREE,
            AdjacencyRelation(4, 4, 1.5));
        requireEqual(reimported.getAltitude(importedSampleNodeId), importedSampleAltitude, "weighted Higra import must repopulate the external altitude buffer");
    }

    {
        auto image = makeComponentTreeFixture();
        auto topology = MorphologicalTree::createMaxTree(image);
        const auto originalAliveNodeIds = collectNodeIds(topology.getAliveNodeIds());

        auto weightedFromTopology = WeightedMorphologicalTree::createFromTopology(topology, image);
        auto rebuiltWeighted = WeightedMorphologicalTree::createComponentTree(image, true);

        weightedFromTopology.validateAltitudeBufferShape();
        weightedFromTopology.validateMonotoneAltitude();
        requireEqual(weightedFromTopology.topology().getRoot(), topology.getRoot(), "weighted-from-topology root");
        require(weightedFromTopology.topology().isMaxtree(), "weighted-from-topology tree type");
        requireVectorEqual(weightedFromTopology.getAltitudeBuffer(), rebuiltWeighted.getAltitudeBuffer(), "weighted-from-topology altitude buffer");
        requireVectorEqual(
            collectImageValues(weightedFromTopology.reconstructionImage()),
            collectImageValues(image),
            "weighted-from-topology reconstruction");

        weightedFromTopology.mergeNodeIntoParent(5);
        requireVectorEqual(
            collectNodeIds(topology.getAliveNodeIds()),
            originalAliveNodeIds,
            "weighted-from-topology must clone the caller topology");
    }

    {
        auto image = makeComponentTreeFixture();
        auto topology = MorphologicalTree::createMinTree(image);
        auto weightedFromTopology = WeightedMorphologicalTree::createFromTopology(topology, image);
        auto rebuiltWeighted = WeightedMorphologicalTree::createComponentTree(image, false);

        weightedFromTopology.validateAltitudeBufferShape();
        weightedFromTopology.validateMonotoneAltitude();
        require(!weightedFromTopology.topology().isMaxtree(), "min weighted-from-topology tree type");
        requireVectorEqual(weightedFromTopology.getAltitudeBuffer(), rebuiltWeighted.getAltitudeBuffer(), "min weighted-from-topology altitude buffer");
    }

    {
        auto topology = MorphologicalTree::createMaxTree(makeComponentTreeFixture());
        auto mismatchedImage = makeImage(
            2,
            2,
            {
                1, 2,
                3, 4,
            });

        requireThrows<std::invalid_argument>(
            [&]() { static_cast<void>(WeightedMorphologicalTree::createFromTopology(topology, mismatchedImage)); },
            "weighted-from-topology must reject mismatched image domains");
    }

    {
        auto weighted = makeWeightedComponentTree(makeComponentTreeFixture(), true);
        weighted->setAltitudeBuffer(AltitudeBuffer(static_cast<size_t>(weighted->topology().getNumInternalNodeSlots()), AltitudeType{7}));

        const auto [higraParent, higraAltitude] = weighted->exportHigraHierarchy();
        requireCompactHigraHierarchy(
            higraParent,
            higraAltitude,
            weighted->topology().getNumTotalProperParts(),
            "equal-altitude Higra export");

        auto roundtrip = WeightedMorphologicalTree::createFromHigraParent(
            higraParent,
            higraAltitude,
            4,
            4,
            MorphologicalTree::MAX_TREE,
            AdjacencyRelation(4, 4, 1.5));
        const auto [reexportedParent, reexportedAltitude] = roundtrip.exportHigraHierarchy();
        requireVectorEqual(reexportedParent, higraParent, "equal-altitude Higra parent round-trip");
        requireVectorEqual(reexportedAltitude, higraAltitude, "equal-altitude Higra altitude round-trip");
    }

    {
        auto weighted = makeWeightedComponentTree(makeComponentTreeFixture(), true);

        const auto [weightedAsc, weightedDesc] = WeightedMorphologicalTree::computeAscendantsAndDescendants(weighted->topology(), &weighted->getAltitudeBuffer(), 2);
        requireVectorEqual(weightedAsc, {InvalidNode, 0, 0, 1, 2, 3}, "weighted level-based ascendants must use the external altitude buffer");
        requireVectorEqual(weightedDesc, {1, 3, 4, 5, InvalidNode, InvalidNode}, "weighted level-based descendants must use the external altitude buffer");
    }

    {
        auto weighted = makeWeightedComponentTree(makeComponentTreeFixture(), true);

        weighted->setAltitude(5, 9);
        requireEqual(weighted->getAltitude(5), 9, "weighted setAltitude must update the external buffer");

        auto altitude = weighted->getAltitudeBuffer();
        altitude[4] = 11;
        weighted->setAltitudeBuffer(altitude);
        requireEqual(weighted->getAltitude(4), 11, "weighted setAltitudeBuffer must replace the external buffer");
    }

    {
        auto weighted = makeWeightedComponentTree(makeComponentTreeFixture(), true);
        weighted->setAltitude(4, 10);

        weighted->mergeNodeIntoParent(5);
        requireEqual(weighted->getAltitude(4), 10, "weighted merge must preserve external altitude values");
    }

    {
        auto weighted = makeWeightedComponentTree(makeComponentTreeFixture(), true);
        weighted->setAltitude(3, 12);

        weighted->pruneNode(5);
        requireEqual(weighted->getAltitude(3), 12, "weighted prune must preserve external altitude values");
    }

    {
        auto weighted = makeWeightedComponentTree(makeComponentTreeFixture(), true);
        const auto rootChildren = collectNodeIds(weighted->topology().getChildren(weighted->topology().getRoot()));
        require(!rootChildren.empty(), "fixture must expose at least one root child");

        weighted->setAltitude(rootChildren.front(), weighted->getAltitude(weighted->topology().getRoot()) - 1);

        bool threw = false;
        try {
            weighted->validateMonotoneAltitude();
        } catch (const std::runtime_error&) {
            threw = true;
        }
        require(threw, "weighted monotone validation must reject a max-tree child below its parent");
    }

    {
        auto tosImage = makeImage(
            3,
            3,
            {
                5, 5, 4,
                5, 3, 3,
                2, 2, 1,
            });
        auto weighted = makeWeightedTreeOfShapes(tosImage, ToSInterpolation::SelfDual);
        weighted->validateAltitudeBufferShape();

        auto topology = MorphologicalTree::createTreeOfShapes(tosImage, ToSInterpolation::SelfDual);
        auto weightedFromTopology = WeightedMorphologicalTree::createFromTopology(topology, tosImage);
        weightedFromTopology.validateAltitudeBufferShape();
        requireVectorEqual(
            weightedFromTopology.getAltitudeBuffer(),
            weighted->getAltitudeBuffer(),
            "tree-of-shapes weighted-from-topology altitude buffer");
        requireVectorEqual(
            collectImageValues(weightedFromTopology.reconstructionImage()),
            collectImageValues(tosImage),
            "tree-of-shapes weighted-from-topology reconstruction");

        auto tosAltitude = weighted->getAltitudeBuffer();
        if (!tosAltitude.empty()) {
            tosAltitude.front() += 100;
            weighted->setAltitudeBuffer(std::move(tosAltitude));
        }
        weighted->validateMonotoneAltitude();
    }

    {
        auto weighted = makeWeightedComponentTree(makeComponentTreeFixture(), false);
        auto editor = weighted->edit();

        auto [areaNamesBeforeEdit, areaBeforeEdit] = AttributeComputedIncrementally::computeSingleAttribute(*weighted, AREA);
        const int expectedInsertedArea =
            static_cast<int>(areaBeforeEdit[areaNamesBeforeEdit.linearIndex(3, AREA)]) +
            static_cast<int>(areaBeforeEdit[areaNamesBeforeEdit.linearIndex(4, AREA)]);
        const AltitudeType insertedAltitude = std::min(
            weighted->getAltitude(2),
            std::max(weighted->getAltitude(3), weighted->getAltitude(4)));
        const NodeId insertedNode = editor.createDetachedNode();

        requireEqual(insertedNode, 6, "weighted editor must append a fresh slot when none is free");
        requireEqual(
            static_cast<int>(weighted->getAltitudeBuffer().size()),
            weighted->topology().getNumInternalNodeSlots(),
            "weighted editor must resize the altitude buffer after node creation");

        editor.setNodeAltitude(insertedNode, insertedAltitude);
        requireEqual(weighted->getAltitude(insertedNode), insertedAltitude, "weighted editor setNodeAltitude");

        editor.reparent(3, insertedNode);
        editor.reparent(4, insertedNode);
        editor.attach(2, insertedNode);
        editor.commit();

        requireEqual(weighted->topology().getNodeParent(insertedNode), 2, "weighted editor inserted node parent after commit");
        requireEqual(weighted->getAltitude(insertedNode), insertedAltitude, "weighted editor inserted node altitude after commit");
        auto [areaNamesAfterEdit, areaAfterEdit] = AttributeComputedIncrementally::computeSingleAttribute(*weighted, AREA);
        requireEqual(
            static_cast<int>(areaAfterEdit[areaNamesAfterEdit.linearIndex(insertedNode, AREA)]),
            expectedInsertedArea,
            "weighted editor inserted node area after commit");

        const auto [higraParent, higraAltitude] = weighted->exportHigraHierarchy();
        requireEqual(
            static_cast<int>(higraParent.size()),
            weighted->topology().getNumTotalProperParts() + weighted->topology().getNumNodes(),
            "weighted editor Higra export parent size after commit");
        requireEqual(
            static_cast<int>(higraAltitude.size()),
            weighted->topology().getNumTotalProperParts() + weighted->topology().getNumNodes(),
            "weighted editor Higra export altitude size after commit");
    }

    return 0;
}
