#include "support/TestSupport.hpp"
#include "mmcfilters/trees/TreeEditor.hpp"

#include <memory>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

int main() {
    {
        // Multiple siblings and grandchildren exercise stack growth and child order.
        const std::vector<NodeId> parents{0, 0, 0, 0, 1, 1, 2, 2, 3, 3};
        const std::vector<NodeId> smallestNodes{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
        const std::vector<int> altitudes(10, 0);
        auto valued = MorphologicalTreeFactory::createFromNativeTopology(
            std::span<const NodeId>(parents), std::span<const NodeId>(smallestNodes), std::span<const int>(altitudes), 0,
            MorphologicalTreeSemantics{});
        const auto& tree = valued.topology();
        requireVectorEqual(collectNodeIds(tree.subtreeNodes(tree.root())), {0, 1, 4, 5, 2, 6, 7, 3, 8, 9}, "branching pre-order preserves sibling order");
        requireVectorEqual(collectNodeIds(tree.postOrder()), {4, 5, 1, 6, 7, 2, 8, 9, 3, 0}, "branching post-order preserves sibling order");
        requireVectorEqual(collectNodeIds(tree.subtreeNodes(2)), {2, 6, 7}, "subtree traversal stops at its root boundary");

        const std::vector<std::pair<NodeId, NodeId>> offlineQueries{{4, 5}, {4, 6}, {6, 7}, {0, 9}, {8, 8}, {InvalidNode, 4}, {4, 20}};
        requireVectorEqual(tree.lowestCommonAncestors(offlineQueries), {1, 0, 2, 0, 8, InvalidNode, InvalidNode},
                           "batched LCA uses offline Tarjan before the RMQ cache exists");

        require(tree.lowestCommonAncestor(4, 5) == 1, "scalar LCA builds the lazy RMQ cache");
        const std::vector<std::pair<NodeId, NodeId>> cachedQueries{{4, 5}, {4, 6}, {6, 7}, {0, 9}, {8, 8}};
        requireVectorEqual(tree.lowestCommonAncestors(cachedQueries), {1, 0, 2, 0, 8}, "batched LCA reuses an existing RMQ cache");
    }
    auto image = makeComponentTreeFixture();
    auto maxTree = makeComponentTree(image, true);
    auto minTree = makeComponentTree(image, false);

    auto maxBfsIds = collectNodeIds(maxTree->breadthFirstTraversal());
    auto maxPostOrderIds = collectNodeIds(maxTree->postOrder());
    auto maxSubtreeIds = collectNodeIds(maxTree->subtreeNodes(maxTree->root()));
    auto maxPathIds = collectNodeIds(maxTree->ancestors(5));
    auto maxNodePathIds = collectNodeIds(maxTree->getPathBetweenNodes(5, 2));
    auto maxDescendantIds = collectNodeIds(maxTree->descendants(maxTree->root()));

    requireVectorEqual(maxBfsIds, {0, 1, 2, 3, 4, 5}, "max-tree breadth-first traversal");
    requireVectorEqual(maxPostOrderIds, {5, 4, 3, 2, 1, 0}, "max-tree post-order traversal");
    requireVectorEqual(maxSubtreeIds, {0, 1, 2, 3, 4, 5}, "max-tree subtree traversal");
    requireVectorEqual(maxPathIds, {5, 4, 3, 2, 1, 0}, "max-tree path to root");
    requireVectorEqual(maxNodePathIds, {5, 4, 3, 2}, "max-tree path between nodes");
    requireVectorEqual(maxDescendantIds, {1, 2, 3, 4, 5}, "max-tree descendants");

    requireVectorEqual(collectNodeIds(maxTree->aliveNodeIds()), {0, 1, 2, 3, 4, 5}, "max-tree valid node ids iterator");
    requireVectorEqual(collectNodeIds(maxTree->postOrder()), {5, 4, 3, 2, 1, 0}, "max-tree post-order iterator");
    requireVectorEqual(collectNodeIds(maxTree->breadthFirstTraversal()), {0, 1, 2, 3, 4, 5}, "max-tree breadth-first iterator");
    requireVectorEqual(collectNodeIds(maxTree->ancestors(5)), {5, 4, 3, 2, 1, 0}, "max-tree path-to-root iterator");
    requireVectorEqual(collectNodeIds(maxTree->descendants(maxTree->root())), {1, 2, 3, 4, 5}, "max-tree descendants iterator");
    requireVectorEqual(collectNodeIds(maxTree->subtreeNodes(maxTree->root())), {0, 1, 2, 3, 4, 5}, "max-tree subtree iterator");
    requireVectorEqual(collectNodeIds(maxTree->nodeSupport(3)), {0, 1, 4, 5, 6, 9, 10, 14}, "max-tree connected component iterator");
    requireVectorEqual(collectNodeIds(maxTree->nodeSupport(maxTree->root())), {15, 8, 12, 13, 2, 3, 7, 11, 0, 1, 4, 5, 6, 9, 10, 14},
                       "max-tree root connected component iterator");
    require(maxTree->lowestCommonAncestor(5, 2) == 2, "max-tree LCA");
    require(maxTree->lowestCommonAncestor(5, 2) == 2, "max-tree lazy cached LCA first query");
    require(maxTree->lowestCommonAncestor(4, 5) == 4, "max-tree lazy cached LCA second query");

    auto minBfsIds = collectNodeIds(minTree->breadthFirstTraversal());
    auto minSubtreeIds = collectNodeIds(minTree->subtreeNodes(minTree->root()));
    auto minPathIds = collectNodeIds(minTree->ancestors(5));
    auto minNodePathIds = collectNodeIds(minTree->getPathBetweenNodes(4, 5));
    auto minDescendantIds = collectNodeIds(minTree->descendants(minTree->root()));

    requireVectorEqual(minBfsIds, {0, 1, 2, 3, 4, 5}, "min-tree breadth-first traversal");
    requireVectorEqual(minSubtreeIds, {0, 1, 2, 3, 5, 4}, "min-tree subtree traversal");
    requireVectorEqual(minPathIds, {5, 3, 2, 1, 0}, "min-tree path to root");
    requireVectorEqual(minNodePathIds, {4, 2, 3, 5}, "min-tree path between nodes");
    requireVectorEqual(minDescendantIds, {1, 2, 3, 5, 4}, "min-tree descendants");
    requireVectorEqual(collectNodeIds(minTree->nodeSupport(2)), {0, 1, 4, 2, 3, 7, 11, 15, 8, 12, 13}, "min-tree connected component iterator");

    maxTree->mergeNodeIntoParent(5);
    requireVectorEqual(collectNodeIds(maxTree->aliveNodeIds()), {0, 1, 2, 3, 4}, "max-tree valid node ids iterator after merge");
    require(maxTree->lowestCommonAncestor(4, 2) == 2, "max-tree lazy cached LCA after topology change");
    auto editor = maxTree->edit();
    auto detachedNodeId = editor.createDetachedNode();
    require(detachedNodeId != InvalidNode, "detached node allocation for path traversal");
    requireVectorEqual(collectNodeIds(maxTree->getPathBetweenNodes(detachedNodeId, detachedNodeId)), {detachedNodeId}, "detached node path to itself");
    requireVectorEqual(collectNodeIds(maxTree->getPathBetweenNodes(detachedNodeId, maxTree->root())), {},
                       "detached node path to connected root should be empty");
    const std::vector<std::pair<NodeId, NodeId>> detachedQueries{{detachedNodeId, detachedNodeId}, {detachedNodeId, maxTree->root()}};
    requireVectorEqual(maxTree->lowestCommonAncestors(detachedQueries), {detachedNodeId, InvalidNode},
                       "batched LCA handles a detached live node consistently");

    return 0;
}
