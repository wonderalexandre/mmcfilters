#include "support/TestSupport.hpp"
#include "mmcfilters/trees/detail/TreeAttributeSamplingNeighborhood.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

void requireCompactHigraHierarchy(const std::vector<NodeId>& parent, const std::vector<std::uint8_t>& altitude, int numPixels, const std::string& label) {
    requireEqual(parent.size(), altitude.size(), label + " parent/altitude size agreement");

    const NodeId totalNodes = static_cast<NodeId>(parent.size());
    int numRoots = 0;

    for (PixelId pixel = 0; pixel < numPixels; ++pixel) {
        const NodeId smallestNodeId = parent[static_cast<size_t>(pixel)];
        require(smallestNodeId >= numPixels && smallestNodeId < totalNodes, label + " leaf parent must be an internal Higra node");
        requireEqual(altitude[static_cast<size_t>(pixel)], altitude[static_cast<size_t>(smallestNodeId)],
                     label + " leaf altitude must match its owning internal node");
    }

    for (NodeId nodeId = numPixels; nodeId < totalNodes; ++nodeId) {
        const NodeId parentNodeId = parent[static_cast<size_t>(nodeId)];
        require(parentNodeId >= numPixels && parentNodeId < totalNodes, label + " internal parent must stay in the internal Higra domain");
        if (parentNodeId == nodeId) {
            ++numRoots;
        } else {
            require(parentNodeId > nodeId, label + " internal parent must appear after its child");
        }
    }

    requireEqual(numRoots, 1, label + " must encode exactly one self-parented root");
}

void requireFloatVectorEqualAllowingNaN(const std::vector<float>& actual, const std::vector<float>& expected, const std::string& label) {
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
        auto valuedTree = makeValuedComponentTree(image, true);

        valuedTree->validateNodeAltitudeBufferShape();
        valuedTree->validateMonotoneNodeAltitudes();

        requireEqual(static_cast<int>(valuedTree->nodeAltitudes().size()), valuedTree->topology().numInternalNodeSlots(),
                     "valuedTree altitude buffer must match the dense internal-node domain");
        NodeId builtSampleNodeId = InvalidNode;
        std::uint8_t builtSampleAltitude = std::uint8_t{};
        for (NodeId nodeId : valuedTree->topology().aliveNodeIds()) {
            if (valuedTree->nodeAltitude(nodeId) != 0) {
                builtSampleNodeId = nodeId;
                builtSampleAltitude = valuedTree->nodeAltitude(nodeId);
                break;
            }
        }
        require(builtSampleNodeId != InvalidNode, "valuedTree image build must expose at least one non-zero internal altitude");
        requireEqual(valuedTree->nodeAltitude(builtSampleNodeId), builtSampleAltitude, "valuedTree image build external altitude sample");
        requireEqual(valuedTree->nodeAltitude(valuedTree->topology().root()), 0, "valuedTree root altitude");
        requireEqual(valuedTree->nodeResidue(5), 1, "valuedTree node residue");

        auto reconstruction = valuedTree->reconstructFromNodeAltitudes();
        requireImageShape(reconstruction, 4, 4);
        requireVectorEqual(collectImageValues(reconstruction), collectImageValues(image), "valuedTree reconstruction must match the original image");

        const auto [higraParent, higraAltitude] = valuedTree->exportHigraHierarchy();
        requireEqual(static_cast<int>(higraParent.size()), valuedTree->topology().numPixels() + valuedTree->topology().numNodes(),
                     "valuedTree Higra export parent size");
        requireEqual(static_cast<int>(higraAltitude.size()), valuedTree->topology().numPixels() + valuedTree->topology().numNodes(),
                     "valuedTree Higra export altitude size");
        requireCompactHigraHierarchy(higraParent, higraAltitude, valuedTree->topology().numPixels(), "valuedTree Higra export");

        auto roundtrip = MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(higraParent), std::span<const std::uint8_t>(higraAltitude), 4,
                                                                         4, MorphologicalTreeKind::MaxTree, RegularGridAdjacency2D(4, 4, 1.5));
        roundtrip.validateNodeAltitudeBufferShape();
        roundtrip.validateMonotoneNodeAltitudes();
        NodeId importedSampleNodeId = InvalidNode;
        std::uint8_t importedSampleAltitude = std::uint8_t{};
        for (NodeId nodeId : roundtrip.topology().aliveNodeIds()) {
            if (roundtrip.nodeAltitude(nodeId) != 0) {
                importedSampleNodeId = nodeId;
                importedSampleAltitude = roundtrip.nodeAltitude(nodeId);
                break;
            }
        }
        require(importedSampleNodeId != InvalidNode, "valuedTree Higra import must expose at least one non-zero internal altitude");
        requireVectorEqual(collectImageValues(roundtrip.reconstructFromNodeAltitudes()), collectImageValues(image), "valuedTree Higra round-trip reconstruction");
        const auto [reexportedParent, reexportedAltitude] = roundtrip.exportHigraHierarchy();
        requireVectorEqual(reexportedParent, higraParent, "valuedTree Higra parent round-trip");
        requireVectorEqual(reexportedAltitude, higraAltitude, "valuedTree Higra altitude round-trip");

        auto [maxDistNames, maxDistByNode] = AttributeComputation::computeSingleAttribute(*valuedTree, MaxDist);
        auto maxDistByExportedHigra = AttributeComputation::projectNodeValuesToExportedHigra(*valuedTree, maxDistNames, maxDistByNode);
        requireEqual(maxDistByExportedHigra.size(), higraParent.size(), "valuedTree exported-Higra MAX_DIST projection size");
        auto [maxDistHigraNames, maxDistByImportedHigra] = AttributeComputation::computeSingleAttribute(roundtrip, MaxDist, NodeIdSpace::Higra);
        requireFloatVectorEqualAllowingNaN(maxDistByExportedHigra, maxDistByImportedHigra, "valuedTree imported/exported Higra MAX_DIST projection");

        auto [areaNames, areaByNode] = AttributeComputation::computeSingleAttribute(*valuedTree, Area);
        std::vector<float> areaAndMaxDist(static_cast<size_t>(valuedTree->topology().numInternalNodeSlots()) * 2, 0.0f);
        for (NodeId nodeId = 0; nodeId < valuedTree->topology().numInternalNodeSlots(); ++nodeId) {
            areaAndMaxDist[static_cast<size_t>(nodeId) * 2] = areaByNode[areaNames.linearIndex(nodeId, Area)];
            areaAndMaxDist[static_cast<size_t>(nodeId) * 2 + 1] = maxDistByNode[maxDistNames.linearIndex(nodeId, MaxDist)];
        }
        AttributeNames areaAndMaxDistNames = makeDenseAttributeNames({Area, MaxDist});
        const auto projectedPair = AttributeComputation::projectNodeValuesToExportedHigra(*valuedTree, areaAndMaxDistNames, areaAndMaxDist);
        requireEqual(projectedPair.size(), maxDistByExportedHigra.size() * 2, "2D exported-Higra projection size");
        requireEqual(projectedPair[0], 1.0f, "2D exported-Higra unit AREA value");
        requireEqual(projectedPair[1], 0.0f, "2D exported-Higra unit MAX_DIST value");

        const std::vector<Attribute> unitProjectionAttributes{MeanGrayLevel, Volume, BoxColumnMin, BoxRowMin};
        const AttributeNames unitProjectionNames = makeDenseAttributeNames(unitProjectionAttributes);
        std::vector<float> unitProjectionNodeValues(
            static_cast<size_t>(valuedTree->topology().numInternalNodeSlots()) * static_cast<size_t>(unitProjectionNames.NUM_ATTRIBUTES), 0.0f);
        for (Attribute attribute : unitProjectionAttributes) {
            auto [singleNames, singleValues] = AttributeComputation::computeSingleAttribute(*valuedTree, attribute);
            for (NodeId nodeId = 0; nodeId < valuedTree->topology().numInternalNodeSlots(); ++nodeId) {
                unitProjectionNodeValues[unitProjectionNames.linearIndex(nodeId, attribute)] = singleValues[singleNames.linearIndex(nodeId, attribute)];
            }
        }

        const auto projectedUnitValues = AttributeComputation::projectNodeValuesToExportedHigra(*valuedTree, unitProjectionNames, unitProjectionNodeValues);
        const NodeId sampleProperPart = 10;
        const NodeId sampleSmallestNode = valuedTree->topology().smallestNode(sampleProperPart);
        const auto [sampleRow, sampleColumn] = ImageUtils::to2D(sampleProperPart, valuedTree->topology().numColumns());
        requireEqual(projectedUnitValues[unitProjectionNames.linearIndex(sampleProperPart, MeanGrayLevel)],
                     static_cast<float>(valuedTree->nodeAltitude(sampleSmallestNode)),
                     "exported-Higra unit MEAN_GRAY_LEVEL must use the proper-part altitude");
        requireEqual(projectedUnitValues[unitProjectionNames.linearIndex(sampleProperPart, Volume)], static_cast<float>(valuedTree->nodeAltitude(sampleSmallestNode)),
                     "exported-Higra unit VOLUME must use one pixel at the proper-part altitude");
        requireEqual(projectedUnitValues[unitProjectionNames.linearIndex(sampleProperPart, BoxColumnMin)], static_cast<float>(sampleColumn),
                     "exported-Higra unit BOX_COLUMN_MIN must use the proper-part column");
        requireEqual(projectedUnitValues[unitProjectionNames.linearIndex(sampleProperPart, BoxRowMin)], static_cast<float>(sampleRow),
                     "exported-Higra unit BOX_ROW_MIN must use the proper-part row");

        if constexpr (contract::validationsEnabled) {
            requireThrows<std::invalid_argument>(
                [&]() {
                    const std::vector<float> invalidValues{1.0f};
                    static_cast<void>(AttributeComputation::projectNodeValuesToExportedHigra(*valuedTree, maxDistNames, invalidValues));
                },
                "exported-Higra projection must reject invalid node-value size");
        }

        roundtrip.setNodeAltitude(importedSampleNodeId, static_cast<std::uint8_t>(importedSampleAltitude + 4));
        auto reimported = MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(higraParent), std::span<const std::uint8_t>(higraAltitude), 4,
                                                                          4, MorphologicalTreeKind::MaxTree, RegularGridAdjacency2D(4, 4, 1.5));
        requireEqual(reimported.nodeAltitude(importedSampleNodeId), importedSampleAltitude,
                     "valuedTree Higra import must repopulate the external altitude buffer");
    }

    {
        auto valuedTree = makeValuedComponentTree(makeComponentTreeFixture(), true);
        const NodeAltitudeBuffer<std::uint8_t> equalAltitude(static_cast<size_t>(valuedTree->topology().numInternalNodeSlots()), std::uint8_t{7});
        requireThrows<std::runtime_error>([&]() { valuedTree->setNodeAltitudes(equalAltitude); },
                                          "ordered component tree must reject equal parent-child altitudes");

        auto [higraParent, higraAltitude] = valuedTree->exportHigraHierarchy();
        std::fill(higraAltitude.begin(), higraAltitude.end(), std::uint8_t{7});
        requireThrows<std::runtime_error>(
            [&]() {
                static_cast<void>(MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(higraParent),
                                                                                  std::span<const std::uint8_t>(higraAltitude), 4, 4,
                                                                                  MorphologicalTreeKind::MaxTree, RegularGridAdjacency2D(4, 4, 1.5)));
            },
            "ordered Higra import must reject equal parent-child altitudes");
    }

    {
        auto valuedTree = makeValuedComponentTree(makeComponentTreeFixture(), true);

        const auto neighborhood = detail::computeNodeAttributeSamplingNeighborhood(
            valuedTree->topology(), std::span<const std::uint8_t>(valuedTree->nodeAltitudes()), static_cast<AltitudeDifference<std::uint8_t>>(2),
            NodeAttributeSamplingPolicy::LargestSupportDescendant);
        requireVectorEqual(neighborhood.ancestors, {InvalidNode, InvalidNode, 0, 1, 2, 3},
                           "valuedTree level-based ancestors must use the external altitude buffer");
        requireVectorEqual(neighborhood.representativeDescendants, {2, 3, 4, 5, InvalidNode, InvalidNode},
                           "valuedTree representative descendants must use the external altitude buffer");
    }

    {
        auto valuedTree = makeValuedComponentTree(makeComponentTreeFixture(), true);

        valuedTree->setNodeAltitude(5, 9);
        requireEqual(valuedTree->nodeAltitude(5), 9, "valuedTree setNodeAltitude must update the external buffer");

        valuedTree->setNodeAltitudes(NodeAltitudeBuffer<std::uint8_t>{0, 1, 2, 3, 4, 5});
        requireEqual(valuedTree->nodeAltitude(5), 5, "valuedTree setNodeAltitudes must replace the external buffer");
    }

    {
        auto valuedTree = makeValuedComponentTree(makeComponentTreeFixture(), true);
        const auto preservedAltitude = valuedTree->nodeAltitude(4);

        valuedTree->mergeNodeIntoParent(5);
        requireEqual(valuedTree->nodeAltitude(4), preservedAltitude, "valuedTree merge must preserve external altitude values");
    }

    {
        auto valuedTree = makeValuedComponentTree(makeComponentTreeFixture(), true);
        const auto preservedAltitude = valuedTree->nodeAltitude(4);

        valuedTree->pruneNode(5);
        requireEqual(valuedTree->nodeAltitude(4), preservedAltitude, "valuedTree prune must preserve external altitude values");
    }

    {
        auto valuedTree = makeValuedComponentTree(makeComponentTreeFixture(), true);
        const auto preservedAltitude = valuedTree->nodeAltitude(4);
        auto editor = valuedTree->edit();

        editor.mergeNodeIntoParent(5);
        editor.commit();

        require(!valuedTree->topology().isAlive(5), "valuedTree editor merge must release the merged node");
        requireEqual(valuedTree->nodeAltitude(4), preservedAltitude, "valuedTree editor merge must preserve external altitude values");
    }

    {
        auto valuedTree = makeValuedComponentTree(makeComponentTreeFixture(), true);
        const auto preservedAltitude = valuedTree->nodeAltitude(4);
        auto editor = valuedTree->edit();

        editor.pruneNode(5);
        editor.commit();

        require(!valuedTree->topology().isAlive(5), "valuedTree editor prune must release the pruned node");
        requireEqual(valuedTree->nodeAltitude(4), preservedAltitude, "valuedTree editor prune must preserve external altitude values");
    }

    {
        auto valuedTree = makeValuedComponentTree(makeComponentTreeFixture(), true);
        const auto rootChildren = collectNodeIds(valuedTree->topology().children(valuedTree->topology().root()));
        require(!rootChildren.empty(), "fixture must expose at least one root child");

        const std::uint8_t equalChildAltitude = valuedTree->nodeAltitude(valuedTree->topology().root());
        requireThrows<std::runtime_error>([&]() { valuedTree->setNodeAltitude(rootChildren.front(), equalChildAltitude); },
                                          "valuedTree setNodeAltitude must reject a max-tree child equal to its parent");

        auto invalidAltitude = valuedTree->nodeAltitudes();
        invalidAltitude[static_cast<size_t>(rootChildren.front())] = equalChildAltitude;
        requireThrows<std::runtime_error>([&]() { valuedTree->setNodeAltitudes(invalidAltitude); },
                                          "valuedTree setNodeAltitudes must reject a max-tree child equal to its parent");

        const std::uint8_t originalChildAltitude = valuedTree->nodeAltitude(rootChildren.front());
        {
            auto editor = valuedTree->edit();
            require(editor.canRollback(), "ordinary valuedTree edits must be recoverable");
            editor.setNodeAltitude(rootChildren.front(), equalChildAltitude);
            requireThrows<std::runtime_error>([&]() { editor.commit(); }, "valuedTree commit must reject a staged non-strict altitude");
        }
        requireEqual(valuedTree->nodeAltitude(rootChildren.front()), originalChildAltitude, "abandoned invalid valuedTree edit must restore the altitude");
        valuedTree->validateMonotoneNodeAltitudes();
    }

    {
        auto tosImage = makeImage(3, 3,
                                  {
                                      5,
                                      5,
                                      4,
                                      5,
                                      3,
                                      3,
                                      2,
                                      2,
                                      1,
                                  });
        auto valuedTree = makeValuedTreeOfShapes(tosImage, TestTopographicImmersion::SelfDualSpan);
        valuedTree->validateNodeAltitudeBufferShape();

        auto tosAltitude = valuedTree->nodeAltitudes();
        if (!tosAltitude.empty()) {
            tosAltitude.front() += 100;
            valuedTree->setNodeAltitudes(std::move(tosAltitude));
        }
        valuedTree->validateMonotoneNodeAltitudes();
    }

    {
        auto valuedTree = makeValuedComponentTree(makeComponentTreeFixture(), false);
        auto spacedAltitude = valuedTree->nodeAltitudes();
        for (auto& value : spacedAltitude) {
            value = static_cast<std::uint8_t>(value * 2);
        }
        valuedTree->setNodeAltitudes(std::move(spacedAltitude));
        auto [areaNamesBeforeEdit, areaBeforeEdit] = AttributeComputation::computeSingleAttribute(*valuedTree, Area);
        const int expectedInsertedArea = static_cast<int>(areaBeforeEdit[areaNamesBeforeEdit.linearIndex(3, Area)]) +
                                         static_cast<int>(areaBeforeEdit[areaNamesBeforeEdit.linearIndex(4, Area)]);
        const std::uint8_t insertedAltitude = static_cast<std::uint8_t>(std::max(valuedTree->nodeAltitude(3), valuedTree->nodeAltitude(4)) + 1);
        auto editor = valuedTree->edit();
        const NodeId insertedNode = editor.createDetachedNode();

        requireEqual(insertedNode, 6, "valuedTree editor must append a fresh slot when none is free");
        requireEqual(static_cast<int>(valuedTree->nodeAltitudes().size()), valuedTree->topology().numInternalNodeSlots(),
                     "valuedTree editor must resize the altitude buffer after node creation");

        editor.setNodeAltitude(insertedNode, insertedAltitude);
        requireEqual(valuedTree->nodeAltitude(insertedNode), insertedAltitude, "valuedTree editor setNodeAltitude");

        editor.reparent(3, insertedNode);
        editor.reparent(4, insertedNode);
        editor.attach(2, insertedNode);
        editor.commit();

        requireEqual(valuedTree->topology().parent(insertedNode), 2, "valuedTree editor inserted node parent after commit");
        requireEqual(valuedTree->nodeAltitude(insertedNode), insertedAltitude, "valuedTree editor inserted node altitude after commit");
        auto [areaNamesAfterEdit, areaAfterEdit] = AttributeComputation::computeSingleAttribute(*valuedTree, Area);
        requireEqual(static_cast<int>(areaAfterEdit[areaNamesAfterEdit.linearIndex(insertedNode, Area)]), expectedInsertedArea,
                     "valuedTree editor inserted node area after commit");

        const auto [higraParent, higraAltitude] = valuedTree->exportHigraHierarchy();
        requireEqual(static_cast<int>(higraParent.size()), valuedTree->topology().numPixels() + valuedTree->topology().numNodes(),
                     "valuedTree editor Higra export parent size after commit");
        requireEqual(static_cast<int>(higraAltitude.size()), valuedTree->topology().numPixels() + valuedTree->topology().numNodes(),
                     "valuedTree editor Higra export altitude size after commit");
    }

    {
        auto valuedTree = makeValuedComponentTree(makeComponentTreeFixture(), true);
        const auto originalAltitude = valuedTree->nodeAltitudes();
        const NodeId root = valuedTree->topology().root();

        {
            auto editor = valuedTree->edit();
            require(editor.canRollback(), "ordinary valuedTree editor must advertise rollback");
            editor.setNodeAltitude(root, std::numeric_limits<std::uint8_t>::max());
            const TreeValidationResult result = editor.validateAndCommit();
            require(!result.ok, "recoverable valuedTree commit must reject invalid altitude order");
        }

        require(!valuedTree->topology().isEditing(), "valuedTree editor destructor must close through rollback");
        requireVectorEqual(valuedTree->nodeAltitudes(), originalAltitude, "valuedTree delta rollback must restore altitude");
        valuedTree->validateMonotoneNodeAltitudes();
    }

    {
        auto valuedTree = makeValuedComponentTree(makeComponentTreeFixture(), true);
        const auto originalAltitude = valuedTree->nodeAltitudes();
        const int originalSlots = valuedTree->topology().numInternalNodeSlots();

        auto editor = valuedTree->edit();
        static_cast<void>(editor.createDetachedNode(std::uint8_t{17}));
        editor.rollback();

        requireEqual(valuedTree->topology().numInternalNodeSlots(), originalSlots, "valuedTree explicit rollback must restore topology slots");
        requireVectorEqual(valuedTree->nodeAltitudes(), originalAltitude, "valuedTree explicit rollback must restore altitude size and values");
    }

    {
        const std::vector<NodeId> nodeParent{2, 2, 2};
        const std::vector<NodeId> smallestNodeMap{0, 1};
        const std::vector<std::uint8_t> altitude{0, 20, 10};
        auto mixedPolarity = MorphologicalTreeFactory::createFromNativeTopology(
            std::span<const NodeId>(nodeParent), std::span<const NodeId>(smallestNodeMap), std::span<const std::uint8_t>(altitude), 2, 1, 2,
            MorphologicalTreeSemantics{MorphologicalTreeKind::TreeOfShapes, NodeAltitudeOrder::Unconstrained, NoConstructionContext{}});

        const auto [parent, exportedAltitude] = mixedPolarity.exportHigraHierarchy();
        requireCompactHigraHierarchy(parent, exportedAltitude, 2, "mixed-polarity Higra export");

        auto roundtrip = MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(parent), std::span<const std::uint8_t>(exportedAltitude), 1, 2,
                                                                         MorphologicalTreeKind::TreeOfShapes);
        const auto [roundtripParent, roundtripAltitude] = roundtrip.exportHigraHierarchy();
        requireVectorEqual(roundtripParent, parent, "mixed-polarity Higra parent round-trip");
        requireVectorEqual(roundtripAltitude, exportedAltitude, "mixed-polarity Higra altitude round-trip");
    }

    {
        auto valuedTree = makeValuedComponentTree(makeComponentTreeFixture(), true);
        const auto originalAltitude = valuedTree->nodeAltitudes();
        auto editor = valuedTree->edit();
        editor.detach(3);

        if constexpr (contract::validationsEnabled) {
            requireThrows<std::logic_error>([&] { static_cast<void>(valuedTree->exportHigraHierarchy()); }, "Higra export must reject a staged topology");
        }
        if constexpr (contract::validationsEnabled) {
            requireThrows<std::logic_error>(
                [&] {
                    ValuedMorphologicalTree<std::uint8_t> moved(std::move(*valuedTree));
                    static_cast<void>(moved);
                },
                "valuedTree move construction must reject an active editor");
            requireVectorEqual(valuedTree->nodeAltitudes(), originalAltitude, "rejected valuedTree move must preserve altitude ownership");
            require(valuedTree->topology().isEditing(), "rejected valuedTree move must preserve the edit session");

            auto assignmentTarget = makeValuedComponentTree(makeComponentTreeFixture(), false);
            const auto targetAltitude = assignmentTarget->nodeAltitudes();
            requireThrows<std::logic_error>([&] { *assignmentTarget = std::move(*valuedTree); },
                                            "valuedTree move assignment must reject an active source editor");
            requireVectorEqual(assignmentTarget->nodeAltitudes(), targetAltitude,
                               "rejected valuedTree move assignment must preserve destination altitude");
            assignmentTarget->topology().validateConnectedRootedTree();
            require(valuedTree->topology().isEditing(), "rejected valuedTree move assignment must preserve the source editor");
        }

        editor.rollback();
        ValuedMorphologicalTree<std::uint8_t> moved(std::move(*valuedTree));
        moved.validateNodeAltitudeBufferShape();
        moved.validateMonotoneNodeAltitudes();
    }

    return 0;
}
