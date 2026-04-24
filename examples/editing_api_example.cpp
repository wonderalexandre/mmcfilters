#include "../mmcfilters/trees/TreeEditor.hpp"
#include "../mmcfilters/trees/WeightedMorphologicalTree.hpp"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

mmcfilters::ImageUInt8Ptr makeFixtureImage()
{
    const std::vector<uint8_t> values = {
        3, 3, 2, 2,
        3, 4, 4, 2,
        1, 4, 5, 2,
        1, 1, 5, 0,
    };

    auto image = mmcfilters::ImageUInt8::create(4, 4);
    for (std::size_t i = 0; i < values.size(); ++i) {
        (*image)[static_cast<int>(i)] = values[i];
    }
    return image;
}

void require(bool condition, const char* message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

mmcfilters::NodeId firstNonRootLeaf(const mmcfilters::MorphologicalTree& tree)
{
    for (mmcfilters::NodeId nodeId : tree.getLeaves()) {
        if (!tree.isRoot(nodeId)) {
            return nodeId;
        }
    }
    throw std::runtime_error("fixture must expose at least one non-root leaf");
}

void safeMutatorExample()
{
    auto tree = mmcfilters::MorphologicalTree::createComponentTree(makeFixtureImage(), true);
    const int nodesBeforePrune = tree.getNumNodes();

    tree.pruneNode(firstNonRootLeaf(tree));
    require(tree.getNumNodes() < nodesBeforePrune, "pruneNode should remove a live subtree");

    auto weighted = mmcfilters::WeightedMorphologicalTree::createComponentTree(makeFixtureImage(), true);
    const int nodesBeforeMerge = weighted.topology().getNumNodes();

    weighted.mergeNodeIntoParent(firstNonRootLeaf(weighted.topology()));
    require(weighted.topology().getNumNodes() < nodesBeforeMerge, "weighted merge should update the owned topology");
}

void treeEditorExample()
{
    auto tree = mmcfilters::MorphologicalTree::createComponentTree(makeFixtureImage(), false);
    auto editor = tree.edit();

    const mmcfilters::NodeId insertedNode = editor.createDetachedNode();
    editor.reparent(3, insertedNode);
    editor.reparent(4, insertedNode);
    editor.attach(2, insertedNode);
    editor.commit();

    require(tree.getNodeParent(insertedNode) == 2, "TreeEditor should validate and commit the staged topology");
}

void weightedTreeEditorExample()
{
    auto weighted = mmcfilters::WeightedMorphologicalTree::createComponentTree(makeFixtureImage(), false);
    auto editor = weighted.edit();

    const mmcfilters::AltitudeType insertedAltitude = std::min(
        weighted.getAltitude(2),
        std::max(weighted.getAltitude(3), weighted.getAltitude(4)));
    const mmcfilters::NodeId insertedNode = editor.createDetachedNode(insertedAltitude);

    editor.reparent(3, insertedNode);
    editor.reparent(4, insertedNode);
    editor.attach(2, insertedNode);
    editor.commit();

    require(weighted.getAltitude(insertedNode) == insertedAltitude, "WeightedTreeEditor should preserve inserted altitude");
}

void rejectedWeightedCommitExample()
{
    auto weighted = mmcfilters::WeightedMorphologicalTree::createComponentTree(makeFixtureImage(), true);
    auto editor = weighted.edit();

    const mmcfilters::NodeId root = weighted.topology().getRoot();
    const mmcfilters::NodeId insertedNode = editor.createDetachedNode(weighted.getAltitude(root) - 1);
    editor.attach(root, insertedNode);

    bool rejected = false;
    try {
        editor.commit();
    } catch (const std::runtime_error&) {
        rejected = true;
    }
    require(rejected, "WeightedTreeEditor::commit should reject non-monotone altitude");
}

} // namespace

int main()
{
    safeMutatorExample();
    treeEditorExample();
    weightedTreeEditorExample();
    rejectedWeightedCommitExample();

    std::cout << "Editing API example completed.\n";
    return 0;
}
