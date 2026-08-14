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

template <class T> std::vector<T> toSorted(std::vector<T> values) {
    std::sort(values.begin(), values.end());
    return values;
}

bool containsNode(const std::vector<NodeId>& nodes, NodeId target) { return std::find(nodes.begin(), nodes.end(), target) != nodes.end(); }

std::vector<NodeId> componentRoots(const MorphologicalTree& tree) {
    std::vector<NodeId> roots{tree.root()};
    for (NodeId nodeId : tree.aliveNodeIds()) {
        if (nodeId != tree.root() && tree.parent(nodeId) == nodeId) {
            roots.push_back(nodeId);
        }
    }
    return roots;
}

void requireTreeInvariantSnapshot(const MorphologicalTree& tree, const std::string& label) {
    // These checks intentionally play the role of post-mutation validation for
    // low-level tree editing operations. Runtime mutators only assert the local
    // preconditions needed to avoid obvious cycles in debug builds.
    const auto aliveNodes = collectNodeIds(tree.aliveNodeIds());
    require(!aliveNodes.empty(), label + ": tree must have alive nodes");
    require(tree.isAlive(tree.root()), label + ": root must be alive");
    requireEqual(tree.parent(tree.root()), tree.root(), label + ": root parent");

    for (NodeId startNodeId : aliveNodes) {
        std::unordered_set<NodeId> parentChain;
        NodeId currentNodeId = startNodeId;
        while (currentNodeId != InvalidNode) {
            require(parentChain.insert(currentNodeId).second, label + ": parent chains must be acyclic");
            const NodeId parentNodeId = tree.parent(currentNodeId);
            if (parentNodeId == currentNodeId) {
                break;
            }
            currentNodeId = parentNodeId;
        }
    }

    int selfParentRoots = 0;
    std::vector<uint8_t> seen(static_cast<size_t>(tree.numInternalNodeSlots()), 0);
    auto forestRoots = componentRoots(tree);

    for (NodeId rootId : forestRoots) {
        const NodeId parentId = tree.parent(rootId);
        if (rootId == tree.root()) {
            requireEqual(parentId, rootId, label + ": main root parent must point to itself");
            ++selfParentRoots;
        } else {
            requireEqual(parentId, rootId, label + ": detached root must point to itself");
        }

        for (NodeId nodeId : tree.subtreeNodes(rootId)) {
            require(tree.isAlive(nodeId), label + ": subtree must contain only alive nodes");
            require(!seen[static_cast<size_t>(nodeId)], label + ": forest components must be disjoint");
            seen[static_cast<size_t>(nodeId)] = 1;

            auto children = collectNodeIds(tree.children(nodeId));
            requireEqual(static_cast<int>(children.size()), tree.numChildren(nodeId), label + ": child count");
            requireEqual(static_cast<int>(children.empty()), static_cast<int>(tree.isLeaf(nodeId)), label + ": leaf marker");

            for (NodeId childId : children) {
                requireEqual(tree.parent(childId), nodeId, label + ": child parent consistency");
                require(tree.hasChild(nodeId, childId), label + ": hasChild consistency");
            }

            auto directProperParts = collectPixelIds(tree.properPart(nodeId));
            requireEqual(static_cast<int>(directProperParts.size()), tree.properPartCardinality(nodeId), label + ": direct proper-part count");
            requireEqual(static_cast<int>(directProperParts.size()), tree.properPartCardinality(nodeId), label + ": node proper-part count");

            std::vector<int> expectedPixelsOfCC;
            for (NodeId subtreeNodeId : tree.subtreeNodes(nodeId)) {
                for (PixelId pixelId : tree.properPart(subtreeNodeId)) {
                    expectedPixelsOfCC.push_back(pixelId);
                }
            }
            std::vector<int> actualPixelsOfCC;
            for (PixelId pixelId = 0; pixelId < tree.numPixels(); ++pixelId) {
                const NodeId smallestNodeId = tree.smallestNode(pixelId);
                if (smallestNodeId == InvalidNode) {
                    continue;
                }
                NodeId currentNodeId = smallestNodeId;
                while (currentNodeId != InvalidNode) {
                    if (currentNodeId == nodeId) {
                        actualPixelsOfCC.push_back(pixelId);
                        break;
                    }
                    const NodeId parentNodeId = tree.parent(currentNodeId);
                    if (parentNodeId == InvalidNode || parentNodeId == currentNodeId) {
                        break;
                    }
                    currentNodeId = parentNodeId;
                }
            }
            requireVectorEqual(toSorted(std::move(actualPixelsOfCC)), toSorted(std::move(expectedPixelsOfCC)), label + ": pixels of connected component");
        }
    }

    requireEqual(selfParentRoots, 1, label + ": exactly one main self-parent root");
    for (NodeId nodeId : aliveNodes) {
        require(seen[static_cast<size_t>(nodeId)], label + ": every alive node must belong to a forest component");
    }

    std::vector<int> smallestNodeMaps(tree.numPixels(), 0);
    for (NodeId nodeId : aliveNodes) {
        for (PixelId pixel : tree.properPart(nodeId)) {
            requireEqual(tree.smallestNode(pixel), nodeId, label + ": smallest-node mapping");
            smallestNodeMaps[static_cast<size_t>(pixel)] += 1;
        }
    }
    for (PixelId pixelId = 0; pixelId < tree.numPixels(); ++pixelId) {
        requireEqual(smallestNodeMaps[static_cast<size_t>(pixelId)], 1, label + ": every proper part must have exactly one smallest node");
    }

    const auto mainSubtree = collectNodeIds(tree.subtreeNodes(tree.root()));
    requireVectorEqual(toSorted(mainSubtree), toSorted(collectNodeIds(tree.breadthFirstTraversal())), label + ": BFS coverage");
    requireVectorEqual(toSorted(mainSubtree), toSorted(collectNodeIds(tree.postOrder())), label + ": post-order coverage");

    for (NodeId nodeId : mainSubtree) {
        auto descendants = collectNodeIds(tree.descendants(nodeId));
        requireEqual(static_cast<int>(descendants.size()), tree.numDescendants(nodeId), label + ": descendant count");

        auto path = collectNodeIds(tree.ancestors(nodeId));
        require(!path.empty(), label + ": path to root must not be empty");
        requireEqual(path.front(), nodeId, label + ": path to root start");
        requireEqual(path.back(), tree.root(), label + ": path to root end");
    }

    if (mainSubtree.size() <= 12) {
        for (NodeId lhs : mainSubtree) {
            for (NodeId rhs : mainSubtree) {
                const NodeId lca = tree.lowestCommonAncestor(lhs, rhs);
                require(containsNode(mainSubtree, lca), label + ": LCA must belong to main component");
                require(tree.isAncestor(lca, lhs), label + ": LCA must be ancestor of lhs");
                require(tree.isAncestor(lca, rhs), label + ": LCA must be ancestor of rhs");
            }
        }
    }
}

template <class URBG> NodeId randomElement(std::vector<NodeId> values, URBG& rng) {
    std::uniform_int_distribution<size_t> dist(0, values.size() - 1);
    return values[dist(rng)];
}

} // namespace

int main() {
    {
        auto tree = makeComponentTree(makeComponentTreeFixture(), true);
        requireTreeInvariantSnapshot(*tree, "initial max-tree forest invariants");

        tree->mergeNodeIntoParent(5);
        requireTreeInvariantSnapshot(*tree, "after merge forest invariants");

        auto editor = tree->edit();
        const NodeId reused = editor.createDetachedNode();
        requireEqual(reused, 5, "reused node id");
        editor.attach(4, reused);
        editor.movePixelToProperPart(5, 4, 5);
        editor.mergeProperParts(5, 4);
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
        std::mt19937 rng(isMaxtree ? 1337u : 4242u);

        for (int step = 0; step < 20; ++step) {
            const auto attachedNodes = collectNodeIds(tree->subtreeNodes(tree->root()));
            std::vector<int> opCodes;

            if (tree->getNumFreeNodeSlots() > 0) {
                opCodes.push_back(0); // allocate + attach
            }

            {
                bool hasSource = false;
                bool hasTarget = attachedNodes.size() > 1;
                for (NodeId nodeId : attachedNodes) {
                    if (!collectPixelIds(tree->properPart(nodeId)).empty()) {
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
                    if (nodeId != tree->root()) {
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
                    if (!collectNodeIds(tree->children(nodeId)).empty()) {
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
                    if (nodeId != tree->root() && tree->isLeaf(nodeId)) {
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
                    if (nodeId != tree->root()) {
                        prunableNodes.push_back(nodeId);
                    }
                }
                if (!prunableNodes.empty()) {
                    opCodes.push_back(6); // prune subtree
                }
            }

            require(!opCodes.empty(), "connected mutation round must always have an available operation");
            const int op = opCodes[std::uniform_int_distribution<size_t>(0, opCodes.size() - 1)(rng)];

            auto editor = tree->edit();
            switch (op) {
            case 0: {
                const NodeId newNodeId = editor.createDetachedNode();
                require(newNodeId != InvalidNode, "createDetachedNode during invariant test");
                std::vector<NodeId> nonRootNodes;
                for (NodeId nodeId : attachedNodes) {
                    if (nodeId != tree->root()) {
                        nonRootNodes.push_back(nodeId);
                    }
                }
                if (!nonRootNodes.empty()) {
                    const NodeId childId = randomElement(nonRootNodes, rng);
                    const NodeId parentId = tree->parent(childId);
                    editor.reparent(childId, newNodeId);
                    editor.attach(parentId, newNodeId);
                } else {
                    const auto rootProperParts = collectPixelIds(tree->properPart(tree->root()));
                    require(!rootProperParts.empty(), "single-node tree must own a proper part");
                    editor.movePixelToProperPart(newNodeId, tree->root(), rootProperParts.front());
                    editor.attach(tree->root(), newNodeId);
                }
                break;
            }
            case 1: {
                std::vector<NodeId> sources;
                for (NodeId nodeId : attachedNodes) {
                    if (!collectPixelIds(tree->properPart(nodeId)).empty()) {
                        sources.push_back(nodeId);
                    }
                }
                const NodeId sourceId = randomElement(sources, rng);
                std::vector<NodeId> targets = attachedNodes;
                targets.erase(std::remove(targets.begin(), targets.end(), sourceId), targets.end());
                const NodeId targetId = randomElement(targets, rng);
                auto directProperParts = collectPixelIds(tree->properPart(sourceId));
                const PixelId pixelId = randomElement(directProperParts, rng);
                editor.movePixelToProperPart(targetId, sourceId, pixelId);
                break;
            }
            case 2: {
                std::vector<NodeId> sources;
                for (NodeId nodeId : attachedNodes) {
                    if (!collectPixelIds(tree->properPart(nodeId)).empty()) {
                        sources.push_back(nodeId);
                    }
                }
                const NodeId sourceId = randomElement(sources, rng);
                std::vector<NodeId> targets = attachedNodes;
                targets.erase(std::remove(targets.begin(), targets.end(), sourceId), targets.end());
                const NodeId targetId = randomElement(targets, rng);
                editor.mergeProperParts(targetId, sourceId);
                break;
            }
            case 3: {
                std::vector<NodeId> movableNodes;
                for (NodeId nodeId : attachedNodes) {
                    if (nodeId != tree->root()) {
                        movableNodes.push_back(nodeId);
                    }
                }
                const NodeId nodeId = randomElement(movableNodes, rng);
                auto subtree = collectNodeIds(tree->subtreeNodes(nodeId));
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
                    if (!collectNodeIds(tree->children(nodeId)).empty()) {
                        sources.push_back(nodeId);
                    }
                }
                const NodeId sourceId = randomElement(sources, rng);
                auto subtree = collectNodeIds(tree->subtreeNodes(sourceId));
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
                    if (nodeId != tree->root() && tree->isLeaf(nodeId)) {
                        mergeableLeaves.push_back(nodeId);
                    }
                }
                editor.mergeNodeIntoParent(randomElement(mergeableLeaves, rng));
                break;
            }
            case 6: {
                std::vector<NodeId> prunableNodes;
                for (NodeId nodeId : attachedNodes) {
                    if (nodeId != tree->root()) {
                        prunableNodes.push_back(nodeId);
                    }
                }
                editor.pruneNode(randomElement(prunableNodes, rng));
                break;
            }
            default:
                require(false, "unexpected invariant-test operation");
            }

            bool removedEmptySubtree = true;
            while (removedEmptySubtree) {
                removedEmptySubtree = false;
                const auto currentNodes = collectNodeIds(tree->subtreeNodes(tree->root()));
                for (NodeId nodeId : currentNodes) {
                    if (nodeId == tree->root()) {
                        continue;
                    }
                    if (collectNodeIds(tree->nodeSupport(nodeId)).empty()) {
                        editor.pruneNode(nodeId);
                        removedEmptySubtree = true;
                        break;
                    }
                }
            }
            editor.commit();

            requireTreeInvariantSnapshot(*tree, (isMaxtree ? "max-tree" : "min-tree") + std::string(" connected mutation step ") + std::to_string(step));

            const auto exportedHigra = exportFlatHigraHierarchy(*tree);
            auto rebuilt = makeTreeFromHigraParent(exportedHigra.first, tree->numRows(), tree->numColumns(), isMaxtree);
            requireTreeInvariantSnapshot(*rebuilt, (isMaxtree ? "max-tree" : "min-tree") + std::string(" rebuilt mutation step ") + std::to_string(step));
            const auto reexportedHigra = exportFlatHigraHierarchy(*rebuilt);
            requireVectorEqual(reexportedHigra.first, exportedHigra.first,
                               (isMaxtree ? "max-tree" : "min-tree") + std::string(" Higra parent round-trip step ") + std::to_string(step));
            requireVectorEqual(reexportedHigra.second, exportedHigra.second,
                               (isMaxtree ? "max-tree" : "min-tree") + std::string(" Higra altitude round-trip step ") + std::to_string(step));
        }
    }

    return 0;
}
