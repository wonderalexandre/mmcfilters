#include "support/TestSupport.hpp"
#include "mmcfilters/trees/TreeEditor.hpp"

#include <memory>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

int main() {
    auto image = makeComponentTreeFixture();
    auto maxTree = makeComponentTree(image, true);
    auto minTree = makeComponentTree(image, false);
    TreeEditor editor(*maxTree);

    auto maxBfsIds = collectNodeIds(maxTree->getIteratorBreadthFirstTraversal());
    auto maxPostOrderIds = collectNodeIds(maxTree->getPostOrderNodes());
    auto maxSubtreeIds = collectNodeIds(maxTree->getNodeSubtree(maxTree->getRoot()));
    auto maxPathIds = collectNodeIds(maxTree->getPathToRootNodes(5));
    auto maxNodePathIds = collectNodeIds(maxTree->getPathBetweenNodes(5, 2));
    auto maxDescendantIds = collectNodeIds(maxTree->getDescendants(maxTree->getRoot()));

    requireVectorEqual(maxBfsIds, {0, 1, 2, 3, 4, 5}, "max-tree breadth-first traversal");
    requireVectorEqual(maxPostOrderIds, {5, 4, 3, 2, 1, 0}, "max-tree post-order traversal");
    requireVectorEqual(maxSubtreeIds, {0, 1, 2, 3, 4, 5}, "max-tree subtree traversal");
    requireVectorEqual(maxPathIds, {5, 4, 3, 2, 1, 0}, "max-tree path to root");
    requireVectorEqual(maxNodePathIds, {5, 4, 3, 2}, "max-tree path between nodes");
    requireVectorEqual(maxDescendantIds, {1, 2, 3, 4, 5}, "max-tree descendants");

    requireVectorEqual(collectNodeIds(maxTree->getAliveNodeIds()), {0, 1, 2, 3, 4, 5}, "max-tree valid node ids iterator");
    requireVectorEqual(collectNodeIds(maxTree->getPostOrderNodes()), {5, 4, 3, 2, 1, 0}, "max-tree post-order iterator");
    requireVectorEqual(collectNodeIds(maxTree->getIteratorBreadthFirstTraversal()), {0, 1, 2, 3, 4, 5}, "max-tree breadth-first iterator");
    requireVectorEqual(collectNodeIds(maxTree->getPathToRootNodes(5)), {5, 4, 3, 2, 1, 0}, "max-tree path-to-root iterator");
    requireVectorEqual(collectNodeIds(maxTree->getDescendants(maxTree->getRoot())), {1, 2, 3, 4, 5}, "max-tree descendants iterator");
    requireVectorEqual(collectNodeIds(maxTree->getNodeSubtree(maxTree->getRoot())), {0, 1, 2, 3, 4, 5}, "max-tree subtree iterator");
    require(maxTree->getLowestCommonAncestor(5, 2) == 2, "max-tree LCA");
    require(maxTree->getLowestCommonAncestor(5, 2) == 2, "max-tree lazy cached LCA first query");
    require(maxTree->getLowestCommonAncestor(4, 5) == 4, "max-tree lazy cached LCA second query");

    auto minBfsIds = collectNodeIds(minTree->getIteratorBreadthFirstTraversal());
    auto minSubtreeIds = collectNodeIds(minTree->getNodeSubtree(minTree->getRoot()));
    auto minPathIds = collectNodeIds(minTree->getPathToRootNodes(5));
    auto minNodePathIds = collectNodeIds(minTree->getPathBetweenNodes(4, 5));
    auto minDescendantIds = collectNodeIds(minTree->getDescendants(minTree->getRoot()));

    requireVectorEqual(minBfsIds, {0, 1, 2, 3, 4, 5}, "min-tree breadth-first traversal");
    requireVectorEqual(minSubtreeIds, {0, 1, 2, 3, 5, 4}, "min-tree subtree traversal");
    requireVectorEqual(minPathIds, {5, 3, 2, 1, 0}, "min-tree path to root");
    requireVectorEqual(minNodePathIds, {4, 2, 3, 5}, "min-tree path between nodes");
    requireVectorEqual(minDescendantIds, {1, 2, 3, 5, 4}, "min-tree descendants");

    maxTree->mergeNodeIntoParent(5);
    requireVectorEqual(collectNodeIds(maxTree->getAliveNodeIds()), {0, 1, 2, 3, 4}, "max-tree valid node ids iterator after merge");
    require(maxTree->getLowestCommonAncestor(4, 2) == 2, "max-tree lazy cached LCA after topology change");
    auto detachedNodeId = editor.createDetachedNode();
    require(detachedNodeId != InvalidNode, "detached node allocation for path traversal");
    requireVectorEqual(collectNodeIds(maxTree->getPathBetweenNodes(detachedNodeId, detachedNodeId)), {detachedNodeId}, "detached node path to itself");
    requireVectorEqual(collectNodeIds(maxTree->getPathBetweenNodes(detachedNodeId, maxTree->getRoot())), {}, "detached node path to connected root should be empty");

    return 0;
}
