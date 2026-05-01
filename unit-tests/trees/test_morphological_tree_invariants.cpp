#include "support/TestSupport.hpp"
#include "mmcfilters/trees/TreeEditor.hpp"

#include <algorithm>
#include <numeric>
#include <random>
#include <set>
#include <unordered_set>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

namespace {

template <class T>
std::vector<T> toSorted(std::vector<T> values) {
    std::sort(values.begin(), values.end());
    return values;
}

bool containsNode(const std::vector<NodeId>& nodes, NodeId target) {
    return std::find(nodes.begin(), nodes.end(), target) != nodes.end();
}

std::vector<NodeId> componentRoots(const MorphologicalTree& tree) {
    std::vector<NodeId> roots{tree.getRoot()};
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        if (nodeId != tree.getRoot() && tree.getNodeParent(nodeId) == nodeId) {
            roots.push_back(nodeId);
        }
    }
    return roots;
}

void requireTreeInvariantSnapshot(const MorphologicalTree& tree, const std::string& label) {
    // These checks intentionally play the role of post-mutation validation for
    // low-level tree editing operations. Runtime mutators only assert the local
    // preconditions needed to avoid obvious cycles in debug builds.
    const auto aliveNodes = collectNodeIds(tree.getAliveNodeIds());
    require(!aliveNodes.empty(), label + ": tree must have alive nodes");
    require(tree.isAlive(tree.getRoot()), label + ": root must be alive");
    requireEqual(tree.getNodeParent(tree.getRoot()), tree.getRoot(), label + ": root parent");

    for (NodeId startNodeId : aliveNodes) {
        std::unordered_set<NodeId> parentChain;
        NodeId currentNodeId = startNodeId;
        while (currentNodeId != InvalidNode) {
            require(parentChain.insert(currentNodeId).second, label + ": parent chains must be acyclic");
            const NodeId parentNodeId = tree.getNodeParent(currentNodeId);
            if (parentNodeId == currentNodeId) {
                break;
            }
            currentNodeId = parentNodeId;
        }
    }

    int selfParentRoots = 0;
    std::vector<uint8_t> seen(static_cast<size_t>(tree.getNumInternalNodeSlots()), 0);
    auto forestRoots = componentRoots(tree);

    for (NodeId rootId : forestRoots) {
        const NodeId parentId = tree.getNodeParent(rootId);
        if (rootId == tree.getRoot()) {
            requireEqual(parentId, rootId, label + ": main root parent must point to itself");
            ++selfParentRoots;
        } else {
            requireEqual(parentId, rootId, label + ": detached root must point to itself");
        }

        for (NodeId nodeId : tree.getNodeSubtree(rootId)) {
            require(tree.isAlive(nodeId), label + ": subtree must contain only alive nodes");
            require(!seen[static_cast<size_t>(nodeId)], label + ": forest components must be disjoint");
            seen[static_cast<size_t>(nodeId)] = 1;

            auto children = collectNodeIds(tree.getChildren(nodeId));
            requireEqual(static_cast<int>(children.size()), tree.getNumChildren(nodeId), label + ": child count");
            requireEqual(static_cast<int>(children.empty()), static_cast<int>(tree.isLeaf(nodeId)), label + ": leaf marker");

            for (NodeId childId : children) {
                requireEqual(tree.getNodeParent(childId), nodeId, label + ": child parent consistency");
                require(tree.hasChild(nodeId, childId), label + ": hasChild consistency");
            }

            auto directProperParts = collectNodeIds(tree.getProperParts(nodeId));
            requireEqual(static_cast<int>(directProperParts.size()), tree.getNumProperParts(nodeId), label + ": direct proper-part count");
            requireEqual(static_cast<int>(directProperParts.size()), tree.getNumProperParts(nodeId), label + ": node proper-part count");

            std::vector<int> expectedPixelsOfCC;
            for (NodeId subtreeNodeId : tree.getNodeSubtree(nodeId)) {
                for (int pixelId : tree.getProperParts(subtreeNodeId)) {
                    expectedPixelsOfCC.push_back(pixelId);
                }
            }
            std::vector<int> actualPixelsOfCC;
            for (int pixelId = 0; pixelId < tree.getNumTotalProperParts(); ++pixelId) {
                const NodeId ownerNodeId = tree.getSmallestComponent(pixelId);
                if (ownerNodeId == InvalidNode) {
                    continue;
                }
                NodeId currentNodeId = ownerNodeId;
                while (currentNodeId != InvalidNode) {
                    if (currentNodeId == nodeId) {
                        actualPixelsOfCC.push_back(pixelId);
                        break;
                    }
                    const NodeId parentNodeId = tree.getNodeParent(currentNodeId);
                    if (parentNodeId == InvalidNode || parentNodeId == currentNodeId) {
                        break;
                    }
                    currentNodeId = parentNodeId;
                }
            }
            requireVectorEqual(
                toSorted(std::move(actualPixelsOfCC)),
                toSorted(std::move(expectedPixelsOfCC)),
                label + ": pixels of connected component"
            );
        }
    }

    requireEqual(selfParentRoots, 1, label + ": exactly one main self-parent root");
    for (NodeId nodeId : aliveNodes) {
        require(seen[static_cast<size_t>(nodeId)], label + ": every alive node must belong to a forest component");
    }

    std::vector<int> properPartOwners(tree.getNumTotalProperParts(), 0);
    for (NodeId nodeId : aliveNodes) {
        for (int properPart : tree.getProperParts(nodeId)) {
            requireEqual(tree.getSmallestComponent(properPart), nodeId, label + ": smallest component ownership");
            properPartOwners[static_cast<size_t>(properPart)] += 1;
        }
    }
    for (int pixelId = 0; pixelId < tree.getNumTotalProperParts(); ++pixelId) {
        requireEqual(properPartOwners[static_cast<size_t>(pixelId)], 1, label + ": every proper part must have exactly one owner");
    }

    const auto mainSubtree = collectNodeIds(tree.getNodeSubtree(tree.getRoot()));
    requireVectorEqual(toSorted(mainSubtree), toSorted(collectNodeIds(tree.getIteratorBreadthFirstTraversal())), label + ": BFS coverage");
    requireVectorEqual(toSorted(mainSubtree), toSorted(collectNodeIds(tree.getPostOrderNodes())), label + ": post-order coverage");

    for (NodeId nodeId : mainSubtree) {
        auto descendants = collectNodeIds(tree.getDescendants(nodeId));
        requireEqual(static_cast<int>(descendants.size()), tree.getNodeNumDescendants(nodeId), label + ": descendant count");

        auto path = collectNodeIds(tree.getPathToRootNodes(nodeId));
        require(!path.empty(), label + ": path to root must not be empty");
        requireEqual(path.front(), nodeId, label + ": path to root start");
        requireEqual(path.back(), tree.getRoot(), label + ": path to root end");
    }

    if (mainSubtree.size() <= 12) {
        for (NodeId lhs : mainSubtree) {
            for (NodeId rhs : mainSubtree) {
                const NodeId lca = tree.getLowestCommonAncestor(lhs, rhs);
                require(containsNode(mainSubtree, lca), label + ": LCA must belong to main component");
                require(tree.isAncestor(lca, lhs), label + ": LCA must be ancestor of lhs");
                require(tree.isAncestor(lca, rhs), label + ": LCA must be ancestor of rhs");
            }
        }
    }
}

template <class URBG>
NodeId randomElement(std::vector<NodeId> values, URBG& rng) {
    std::uniform_int_distribution<size_t> dist(0, values.size() - 1);
    return values[dist(rng)];
}

} // namespace

int main() {
    {
        auto tree = makeComponentTree(makeComponentTreeFixture(), true);
        auto editor = tree->edit();
        requireTreeInvariantSnapshot(*tree, "initial max-tree forest invariants");

        tree->mergeNodeIntoParent(5);
        requireTreeInvariantSnapshot(*tree, "after merge forest invariants");

        const NodeId reused = editor.createDetachedNode();
        requireEqual(reused, 5, "reused node id");
        editor.attach(4, reused);
        editor.moveProperPart(5, 4, 5);
        editor.moveProperParts(5, 4);
        requireTreeInvariantSnapshot(*tree, "after attach and proper-part moves");

        editor.detach(5);
        requireTreeInvariantSnapshot(*tree, "after detach forest invariants");

        editor.attach(3, 5);
        editor.reparent(5, 2);
        requireTreeInvariantSnapshot(*tree, "after moveNode forest invariants");

        editor.setRoot(2);
        requireTreeInvariantSnapshot(*tree, "after setRoot forest invariants");
    }

    for (bool isMaxtree : {true, false}) {
        auto tree = makeComponentTree(makeComponentTreeFixture(), isMaxtree);
        auto editor = tree->edit();
        std::mt19937 rng(isMaxtree ? 1337u : 4242u);

        for (int step = 0; step < 20; ++step) {
            const auto attachedNodes = collectNodeIds(tree->getNodeSubtree(tree->getRoot()));
            std::vector<int> opCodes;

            if (tree->getNumFreeNodeSlots() > 0) {
                opCodes.push_back(0); // allocate + attach
            }

            {
                bool hasSource = false;
                bool hasTarget = attachedNodes.size() > 1;
                for (NodeId nodeId : attachedNodes) {
                    if (!collectNodeIds(tree->getProperParts(nodeId)).empty()) {
                        hasSource = true;
                        break;
                    }
                }
                if (hasSource && hasTarget) {
                    opCodes.push_back(1); // move one proper part
                    opCodes.push_back(2); // move all direct proper parts
                }
            }

            {
                std::vector<NodeId> movableNodes;
                for (NodeId nodeId : attachedNodes) {
                    if (nodeId != tree->getRoot()) {
                        movableNodes.push_back(nodeId);
                    }
                }
                if (!movableNodes.empty()) {
                    opCodes.push_back(3); // move node
                }
            }

            {
                bool hasSourceWithChildren = false;
                for (NodeId nodeId : attachedNodes) {
                    if (!collectNodeIds(tree->getChildren(nodeId)).empty()) {
                        hasSourceWithChildren = true;
                        break;
                    }
                }
                if (hasSourceWithChildren && attachedNodes.size() > 1) {
                    opCodes.push_back(4); // move children
                }
            }

            {
                std::vector<NodeId> mergeableLeaves;
                for (NodeId nodeId : attachedNodes) {
                    if (nodeId != tree->getRoot() && tree->isLeaf(nodeId)) {
                        mergeableLeaves.push_back(nodeId);
                    }
                }
                if (!mergeableLeaves.empty()) {
                    opCodes.push_back(5); // merge leaf
                }
            }

            {
                std::vector<NodeId> prunableNodes;
                for (NodeId nodeId : attachedNodes) {
                    if (nodeId != tree->getRoot()) {
                        prunableNodes.push_back(nodeId);
                    }
                }
                if (!prunableNodes.empty()) {
                    opCodes.push_back(6); // prune subtree
                }
            }

            require(!opCodes.empty(), "connected mutation round must always have an available operation");
            const int op = opCodes[std::uniform_int_distribution<size_t>(0, opCodes.size() - 1)(rng)];

            switch (op) {
                case 0: {
                    const NodeId newNodeId = editor.createDetachedNode();
                    require(newNodeId != InvalidNode, "createDetachedNode during invariant test");
                    editor.attach(randomElement(attachedNodes, rng), newNodeId);
                    break;
                }
                case 1: {
                    std::vector<NodeId> sources;
                    for (NodeId nodeId : attachedNodes) {
                        if (!collectNodeIds(tree->getProperParts(nodeId)).empty()) {
                            sources.push_back(nodeId);
                        }
                    }
                    const NodeId sourceId = randomElement(sources, rng);
                    std::vector<NodeId> targets = attachedNodes;
                    targets.erase(std::remove(targets.begin(), targets.end(), sourceId), targets.end());
                    const NodeId targetId = randomElement(targets, rng);
                    auto directProperParts = collectNodeIds(tree->getProperParts(sourceId));
                    const NodeId pixelId = randomElement(directProperParts, rng);
                    editor.moveProperPart(targetId, sourceId, pixelId);
                    break;
                }
                case 2: {
                    std::vector<NodeId> sources;
                    for (NodeId nodeId : attachedNodes) {
                        if (!collectNodeIds(tree->getProperParts(nodeId)).empty()) {
                            sources.push_back(nodeId);
                        }
                    }
                    const NodeId sourceId = randomElement(sources, rng);
                    std::vector<NodeId> targets = attachedNodes;
                    targets.erase(std::remove(targets.begin(), targets.end(), sourceId), targets.end());
                    const NodeId targetId = randomElement(targets, rng);
                    editor.moveProperParts(targetId, sourceId);
                    break;
                }
                case 3: {
                    std::vector<NodeId> movableNodes;
                    for (NodeId nodeId : attachedNodes) {
                        if (nodeId != tree->getRoot()) {
                            movableNodes.push_back(nodeId);
                        }
                    }
                    const NodeId nodeId = randomElement(movableNodes, rng);
                    auto subtree = collectNodeIds(tree->getNodeSubtree(nodeId));
                    std::vector<NodeId> possibleParents;
                    for (NodeId candidateId : attachedNodes) {
                        if (!containsNode(subtree, candidateId)) {
                            possibleParents.push_back(candidateId);
                        }
                    }
                    if (!possibleParents.empty()) {
                        editor.reparent(nodeId, randomElement(possibleParents, rng));
                    }
                    break;
                }
                case 4: {
                    std::vector<NodeId> sources;
                    for (NodeId nodeId : attachedNodes) {
                        if (!collectNodeIds(tree->getChildren(nodeId)).empty()) {
                            sources.push_back(nodeId);
                        }
                    }
                    const NodeId sourceId = randomElement(sources, rng);
                    auto subtree = collectNodeIds(tree->getNodeSubtree(sourceId));
                    std::vector<NodeId> possibleTargets;
                    for (NodeId candidateId : attachedNodes) {
                        if (candidateId != sourceId && !containsNode(subtree, candidateId)) {
                            possibleTargets.push_back(candidateId);
                        }
                    }
                    if (!possibleTargets.empty()) {
                        editor.moveChildren(randomElement(possibleTargets, rng), sourceId);
                    }
                    break;
                }
                case 5: {
                    std::vector<NodeId> mergeableLeaves;
                    for (NodeId nodeId : attachedNodes) {
                        if (nodeId != tree->getRoot() && tree->isLeaf(nodeId)) {
                            mergeableLeaves.push_back(nodeId);
                        }
                    }
                    tree->mergeNodeIntoParent(randomElement(mergeableLeaves, rng));
                    break;
                }
                case 6: {
                    std::vector<NodeId> prunableNodes;
                    for (NodeId nodeId : attachedNodes) {
                        if (nodeId != tree->getRoot()) {
                            prunableNodes.push_back(nodeId);
                        }
                    }
                    tree->pruneNode(randomElement(prunableNodes, rng));
                    break;
                }
                default:
                    require(false, "unexpected invariant-test operation");
            }

            requireTreeInvariantSnapshot(*tree, (isMaxtree ? "max-tree" : "min-tree") + std::string(" connected mutation step ") + std::to_string(step));

            const auto exportedHigra = exportFlatHigraHierarchy(*tree);
            auto rebuilt = makeTreeFromHigraParent(
                exportedHigra.first,
                tree->getNumRowsOfImage(),
                tree->getNumColsOfImage(),
                isMaxtree);
            requireTreeInvariantSnapshot(*rebuilt, (isMaxtree ? "max-tree" : "min-tree") + std::string(" rebuilt mutation step ") + std::to_string(step));
            const auto reexportedHigra = exportFlatHigraHierarchy(*rebuilt);
            requireVectorEqual(reexportedHigra.first, exportedHigra.first, (isMaxtree ? "max-tree" : "min-tree") + std::string(" Higra parent round-trip step ") + std::to_string(step));
            requireVectorEqual(reexportedHigra.second, exportedHigra.second, (isMaxtree ? "max-tree" : "min-tree") + std::string(" Higra altitude round-trip step ") + std::to_string(step));
        }
    }

    return 0;
}
