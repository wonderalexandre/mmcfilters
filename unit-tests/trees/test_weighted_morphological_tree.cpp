#include "support/TestSupport.hpp"
#include "mmcfilters/trees/detail/TreeAltitudeDeltaNeighborhood.hpp"

#include <cmath>
#include <memory>
#include <span>
#include <stdexcept>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

void requireCompactHigraHierarchy(
    const std::vector<NodeId>& parent,
    const std::vector<std::uint8_t>& altitude,
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
        std::uint8_t builtSampleAltitude = std::uint8_t{};
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

        auto roundtrip = MorphologicalTreeFactory::createFromHigraParent(
            std::span<const NodeId>(higraParent),
            std::span<const std::uint8_t>(higraAltitude),
            4,
            4,
            MorphologicalTreeKind::MAX_TREE,
            AdjacencyRelation(4, 4, 1.5));
        roundtrip.validateAltitudeBufferShape();
        roundtrip.validateMonotoneAltitude();
        NodeId importedSampleNodeId = InvalidNode;
        std::uint8_t importedSampleAltitude = std::uint8_t{};
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

        auto [maxDistNames, maxDistByNode] = AttributeComputation::computeSingleAttribute(*weighted, MAX_DIST);
        auto maxDistByExportedHigra = AttributeComputation::projectNodeValuesToExportedHigra(
            *weighted,
            maxDistNames,
            maxDistByNode);
        requireEqual(
            maxDistByExportedHigra.size(),
            higraParent.size(),
            "weighted exported-Higra MAX_DIST projection size");
        auto [maxDistHigraNames, maxDistByImportedHigra] = AttributeComputation::computeSingleAttribute(
            roundtrip,
            MAX_DIST,
            NodeIdSpace::HIGRA);
        requireFloatVectorEqualAllowingNaN(
            maxDistByExportedHigra,
            maxDistByImportedHigra,
            "weighted imported/exported Higra MAX_DIST projection");

        auto [areaNames, areaByNode] = AttributeComputation::computeSingleAttribute(*weighted, AREA);
        std::vector<float> areaAndMaxDist(static_cast<size_t>(weighted->topology().getNumInternalNodeSlots()) * 2, 0.0f);
        for (NodeId nodeId = 0; nodeId < weighted->topology().getNumInternalNodeSlots(); ++nodeId) {
            areaAndMaxDist[static_cast<size_t>(nodeId) * 2] = areaByNode[areaNames.linearIndex(nodeId, AREA)];
            areaAndMaxDist[static_cast<size_t>(nodeId) * 2 + 1] = maxDistByNode[maxDistNames.linearIndex(nodeId, MAX_DIST)];
        }
        AttributeNames areaAndMaxDistNames = makeDenseAttributeNames({AREA, MAX_DIST});
        const auto projectedPair = AttributeComputation::projectNodeValuesToExportedHigra(
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
            auto [singleNames, singleValues] = AttributeComputation::computeSingleAttribute(*weighted, attribute);
            for (NodeId nodeId = 0; nodeId < weighted->topology().getNumInternalNodeSlots(); ++nodeId) {
                unitProjectionNodeValues[unitProjectionNames.linearIndex(nodeId, attribute)] =
                    singleValues[singleNames.linearIndex(nodeId, attribute)];
            }
        }

        const auto projectedUnitValues = AttributeComputation::projectNodeValuesToExportedHigra(
            *weighted,
            unitProjectionNames,
            unitProjectionNodeValues);
        const NodeId sampleProperPart = 10;
        const NodeId sampleOwner = weighted->topology().getProperPartOwner(sampleProperPart);
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
                static_cast<void>(AttributeComputation::projectNodeValuesToExportedHigra(*weighted, maxDistNames, invalidValues));
            },
            "exported-Higra projection must reject invalid node-value size");

        roundtrip.setAltitude(importedSampleNodeId, static_cast<std::uint8_t>(importedSampleAltitude + 4));
        auto reimported = MorphologicalTreeFactory::createFromHigraParent(
            std::span<const NodeId>(higraParent),
            std::span<const std::uint8_t>(higraAltitude),
            4,
            4,
            MorphologicalTreeKind::MAX_TREE,
            AdjacencyRelation(4, 4, 1.5));
        requireEqual(reimported.getAltitude(importedSampleNodeId), importedSampleAltitude, "weighted Higra import must repopulate the external altitude buffer");
    }

    {
        auto weighted = makeWeightedComponentTree(makeComponentTreeFixture(), true);
        weighted->setAltitudeBuffer(AltitudeBuffer<std::uint8_t>(static_cast<size_t>(weighted->topology().getNumInternalNodeSlots()), std::uint8_t{7}));

        const auto [higraParent, higraAltitude] = weighted->exportHigraHierarchy();
        requireCompactHigraHierarchy(
            higraParent,
            higraAltitude,
            weighted->topology().getNumTotalProperParts(),
            "equal-altitude Higra export");

        auto roundtrip = MorphologicalTreeFactory::createFromHigraParent(
            std::span<const NodeId>(higraParent),
            std::span<const std::uint8_t>(higraAltitude),
            4,
            4,
            MorphologicalTreeKind::MAX_TREE,
            AdjacencyRelation(4, 4, 1.5));
        const auto [reexportedParent, reexportedAltitude] = roundtrip.exportHigraHierarchy();
        requireVectorEqual(reexportedParent, higraParent, "equal-altitude Higra parent round-trip");
        requireVectorEqual(reexportedAltitude, higraAltitude, "equal-altitude Higra altitude round-trip");
    }

    {
        auto weighted = makeWeightedComponentTree(makeComponentTreeFixture(), true);

        const auto [weightedAsc, weightedDesc] = detail::computeAscendantsAndDescendantsByAltitude(
            weighted->topology(),
            std::span<const std::uint8_t>(weighted->getAltitudeBuffer()),
            static_cast<AltitudeDiff<std::uint8_t>>(2));
        requireVectorEqual(weightedAsc, {InvalidNode, 0, 0, 1, 2, 3}, "weighted level-based ascendants must use the external altitude buffer");
        requireVectorEqual(weightedDesc, {1, 3, 4, 5, InvalidNode, InvalidNode}, "weighted level-based descendants must use the external altitude buffer");
    }

    {
        auto weighted = makeWeightedComponentTree(makeComponentTreeFixture(), true);

        weighted->setAltitude(5, 9);
        requireEqual(weighted->getAltitude(5), 9, "weighted setAltitude must update the external buffer");

        auto altitude = weighted->getAltitudeBuffer();
        altitude.assign(altitude.size(), std::uint8_t{11});
        weighted->setAltitudeBuffer(altitude);
        requireEqual(weighted->getAltitude(weighted->topology().getRoot()), 11, "weighted setAltitudeBuffer must replace the external buffer");
    }

    {
        auto weighted = makeWeightedComponentTree(makeComponentTreeFixture(), true);
        weighted->setAltitude(4, 5);

        weighted->mergeNodeIntoParent(5);
        requireEqual(weighted->getAltitude(4), 5, "weighted merge must preserve external altitude values");
    }

    {
        auto weighted = makeWeightedComponentTree(makeComponentTreeFixture(), true);
        weighted->setAltitude(4, 5);

        weighted->pruneNode(5);
        requireEqual(weighted->getAltitude(4), 5, "weighted prune must preserve external altitude values");
    }

    {
        auto weighted = makeWeightedComponentTree(makeComponentTreeFixture(), true);
        weighted->setAltitude(4, 5);
        auto editor = weighted->edit();

        editor.mergeNodeIntoParent(5);
        editor.commit();

        require(!weighted->topology().isAlive(5), "weighted editor merge must release the merged node");
        requireEqual(weighted->getAltitude(4), 5, "weighted editor merge must preserve external altitude values");
    }

    {
        auto weighted = makeWeightedComponentTree(makeComponentTreeFixture(), true);
        weighted->setAltitude(4, 5);
        auto editor = weighted->edit();

        editor.pruneNode(5);
        editor.commit();

        require(!weighted->topology().isAlive(5), "weighted editor prune must release the pruned node");
        requireEqual(weighted->getAltitude(4), 5, "weighted editor prune must preserve external altitude values");
    }

    {
        auto weighted = makeWeightedComponentTree(makeComponentTreeFixture(), true);
        const auto rootChildren = collectNodeIds(weighted->topology().getChildren(weighted->topology().getRoot()));
        require(!rootChildren.empty(), "fixture must expose at least one root child");

        const std::uint8_t invalidChildAltitude = weighted->getAltitude(weighted->topology().getRoot()) - 1;
        requireThrows<std::runtime_error>(
            [&]() {
                weighted->setAltitude(rootChildren.front(), invalidChildAltitude);
            },
            "weighted setAltitude must reject a max-tree child below its parent");

        auto invalidAltitude = weighted->getAltitudeBuffer();
        invalidAltitude[static_cast<size_t>(rootChildren.front())] = invalidChildAltitude;
        requireThrows<std::runtime_error>(
            [&]() {
                weighted->setAltitudeBuffer(invalidAltitude);
            },
            "weighted setAltitudeBuffer must reject a max-tree child below its parent");

        weighted->setAltitudeUnchecked(rootChildren.front(), invalidChildAltitude);
        requireThrows<std::runtime_error>(
            [&]() {
                weighted->validateMonotoneAltitude();
            },
            "weighted monotone validation must reject unchecked max-tree altitude edits");
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

        auto tosAltitude = weighted->getAltitudeBuffer();
        if (!tosAltitude.empty()) {
            tosAltitude.front() += 100;
            weighted->setAltitudeBuffer(std::move(tosAltitude));
        }
        weighted->validateMonotoneAltitude();
    }

    {
        auto weighted = makeWeightedComponentTree(makeComponentTreeFixture(), false);
        auto [areaNamesBeforeEdit, areaBeforeEdit] = AttributeComputation::computeSingleAttribute(*weighted, AREA);
        const int expectedInsertedArea =
            static_cast<int>(areaBeforeEdit[areaNamesBeforeEdit.linearIndex(3, AREA)]) +
            static_cast<int>(areaBeforeEdit[areaNamesBeforeEdit.linearIndex(4, AREA)]);
        const std::uint8_t insertedAltitude = std::min(
            weighted->getAltitude(2),
            std::max(weighted->getAltitude(3), weighted->getAltitude(4)));
        auto editor = weighted->edit();
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
        auto [areaNamesAfterEdit, areaAfterEdit] = AttributeComputation::computeSingleAttribute(*weighted, AREA);
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
