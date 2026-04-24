#include "support/TestSupport.hpp"
#include "../../mmcfilters/trees/TreeAltitudeOps.hpp"

#include <memory>
#include <stdexcept>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

int main() {
    {
        auto image = makeComponentTreeFixture();
        auto weighted = std::make_shared<WeightedMorphologicalTree>(image, true);

        weighted->validateAltitudeBufferShape();
        weighted->validateMonotoneAltitude();

        requireEqual(
            static_cast<int>(weighted->altitude.size()),
            weighted->tree.getNumInternalNodeSlots(),
            "weighted altitude buffer must match the dense internal-node domain");
        NodeId builtSampleNodeId = InvalidNode;
        AltitudeType builtSampleAltitude = AltitudeType{};
        for (NodeId nodeId : weighted->tree.getAliveNodeIds()) {
            if (weighted->getAltitude(nodeId) != 0) {
                builtSampleNodeId = nodeId;
                builtSampleAltitude = weighted->getAltitude(nodeId);
                break;
            }
        }
        require(builtSampleNodeId != InvalidNode, "weighted image build must expose at least one non-zero internal altitude");
        requireEqual(weighted->getAltitude(builtSampleNodeId), builtSampleAltitude, "weighted image build external altitude sample");
        requireEqual(weighted->getAltitude(weighted->tree.getRoot()), 0, "weighted root altitude");
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
            weighted->tree.getNumTotalProperParts() + weighted->tree.getNumNodes(),
            "weighted Higra export parent size");
        requireEqual(
            static_cast<int>(higraAltitude.size()),
            weighted->tree.getNumTotalProperParts() + weighted->tree.getNumNodes(),
            "weighted Higra export altitude size");

        auto roundtrip = WeightedMorphologicalTree::createFromHigra(
            higraParent,
            higraAltitude,
            16,
            4,
            4,
            true,
            AdjacencyRelation(4, 4, 1.5));
        roundtrip.validateAltitudeBufferShape();
        roundtrip.validateMonotoneAltitude();
        NodeId importedSampleNodeId = InvalidNode;
        AltitudeType importedSampleAltitude = AltitudeType{};
        for (NodeId nodeId : roundtrip.tree.getAliveNodeIds()) {
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

        roundtrip.setAltitude(importedSampleNodeId, static_cast<AltitudeType>(importedSampleAltitude + 4));
        roundtrip.resetFromHigra(higraParent, higraAltitude, 16);
        requireEqual(roundtrip.getAltitude(importedSampleNodeId), importedSampleAltitude, "weighted Higra reset must repopulate the external altitude buffer");
    }

    {
        auto weighted = std::make_shared<WeightedMorphologicalTree>(makeComponentTreeFixture(), true);

        const auto [weightedAsc, weightedDesc] = tree_altitude_ops::computeAscendantsAndDescendants(weighted->tree, &weighted->altitude, 2, true);
        requireVectorEqual(weightedAsc, {0, 0, 0, 1, 2, 3}, "weighted level-based ascendants must use the external altitude buffer");
        requireVectorEqual(weightedDesc, {0, 3, 4, 5, InvalidNode, InvalidNode}, "weighted level-based descendants must use the external altitude buffer");
    }

    {
        auto weighted = std::make_shared<WeightedMorphologicalTree>(makeComponentTreeFixture(), true);

        weighted->setAltitude(5, 9);
        requireEqual(weighted->getAltitude(5), 9, "weighted setAltitude must update the external buffer");

        auto altitude = weighted->getAltitudeBuffer();
        altitude[4] = 11;
        weighted->setAltitudeBuffer(altitude);
        requireEqual(weighted->getAltitude(4), 11, "weighted setAltitudeBuffer must replace the external buffer");
    }

    {
        auto weighted = std::make_shared<WeightedMorphologicalTree>(makeComponentTreeFixture(), true);
        weighted->setAltitude(4, 10);

        weighted->mergeNodeIntoParent(5);
        requireEqual(weighted->getAltitude(4), 10, "weighted merge must preserve external altitude values");
    }

    {
        auto weighted = std::make_shared<WeightedMorphologicalTree>(makeComponentTreeFixture(), true);
        weighted->setAltitude(3, 12);

        weighted->pruneNode(5);
        requireEqual(weighted->getAltitude(3), 12, "weighted prune must preserve external altitude values");
    }

    {
        auto weighted = std::make_shared<WeightedMorphologicalTree>(makeComponentTreeFixture(), true);
        const auto rootChildren = collectNodeIds(weighted->tree.getChildren(weighted->tree.getRoot()));
        require(!rootChildren.empty(), "fixture must expose at least one root child");

        weighted->altitude[static_cast<size_t>(rootChildren.front())] = weighted->getAltitude(weighted->tree.getRoot()) - 1;

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
        auto weighted = std::make_shared<WeightedMorphologicalTree>(tosImage, ToSInterpolation::SelfDual);
        weighted->validateAltitudeBufferShape();
        if (!weighted->altitude.empty()) {
            weighted->altitude.front() += 100;
        }
        weighted->validateMonotoneAltitude();
    }

    {
        auto weighted = std::make_shared<WeightedMorphologicalTree>(makeComponentTreeFixture(), false);
        auto editor = weighted->edit();

        const int expectedInsertedArea =
            computeAreaAttribute(weighted->tree, 3) + computeAreaAttribute(weighted->tree, 4);
        const AltitudeType insertedAltitude = std::min(
            weighted->getAltitude(2),
            std::max(weighted->getAltitude(3), weighted->getAltitude(4)));
        const NodeId insertedNode = editor.createDetachedNode();

        requireEqual(insertedNode, 6, "weighted editor must append a fresh slot when none is free");
        requireEqual(
            static_cast<int>(weighted->altitude.size()),
            weighted->tree.getNumInternalNodeSlots(),
            "weighted editor must resize the altitude buffer after node creation");

        editor.setNodeAltitude(insertedNode, insertedAltitude);
        requireEqual(weighted->getAltitude(insertedNode), insertedAltitude, "weighted editor setNodeAltitude");

        editor.reparent(3, insertedNode);
        editor.reparent(4, insertedNode);
        editor.attach(2, insertedNode);
        editor.commit();

        requireEqual(weighted->tree.getNodeParent(insertedNode), 2, "weighted editor inserted node parent after commit");
        requireEqual(weighted->getAltitude(insertedNode), insertedAltitude, "weighted editor inserted node altitude after commit");
        requireEqual(
            computeAreaAttribute(weighted->tree, insertedNode),
            expectedInsertedArea,
            "weighted editor inserted node area after commit");

        const auto [higraParent, higraAltitude] = weighted->exportHigraHierarchy();
        requireEqual(
            static_cast<int>(higraParent.size()),
            weighted->tree.getNumTotalProperParts() + weighted->tree.getNumNodes(),
            "weighted editor Higra export parent size after commit");
        requireEqual(
            static_cast<int>(higraAltitude.size()),
            weighted->tree.getNumTotalProperParts() + weighted->tree.getNumNodes(),
            "weighted editor Higra export altitude size after commit");
    }

    return 0;
}
