#include "support/TestSupport.hpp"
#include "mmcfilters/filters/AttributeFilters.hpp"
#include "mmcfilters/trees/detail/TreeKernelAccess.hpp"

#include <numeric>
#include <tuple>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

namespace {

using Event = std::tuple<char, NodeId, NodeId>;

std::vector<Event> events(const MorphologicalTree& tree, NodeId root) {
    std::vector<Event> result;
    detail::TreeKernelAccess::traversePostOrder(
        tree, root,
        [&](NodeId node) { result.emplace_back('b', node, InvalidNode); },
        [&](NodeId parent, NodeId child) { result.emplace_back('m', parent, child); },
        [&](NodeId node) { result.emplace_back('e', node, InvalidNode); });
    return result;
}

void checkCallbackOrder() {
    const std::vector<NodeId> parents{4, 4, 0, 1, 4};
    const std::vector<NodeId> owner{0, 1, 2, 3, 4};
    const std::vector<int> altitude(5, 0);
    auto valued = makeValuedHierarchy<int>(parents, owner, altitude, 4, MorphologicalTreeSemantics{});
    const auto& tree = valued.topology();
    const std::vector<Event> expected{
        {'b',4,-1}, {'b',0,-1}, {'b',2,-1}, {'e',2,-1}, {'m',0,2}, {'e',0,-1}, {'m',4,0},
        {'b',1,-1}, {'b',3,-1}, {'e',3,-1}, {'m',1,3}, {'e',1,-1}, {'m',4,1}, {'e',4,-1}};
    require(events(tree,4) == expected, "callbacks preserve depth-first child order with a nonzero root");
    require(events(tree,0) == std::vector<Event>{{'b',0,-1},{'b',2,-1},{'e',2,-1},{'m',0,2},{'e',0,-1}},
            "subtree callbacks stop before visiting siblings or merging into the outer parent");
    require(events(tree,2) == std::vector<Event>{{'b',2,-1},{'e',2,-1}}, "a leaf has only pre and post callbacks");
    valued.mergeNodeIntoParent(2);
    require(events(valued.topology(),0) == std::vector<Event>{{'b',0,-1},{'e',0,-1}},
            "callbacks skip inactive slots after an edit");
}

std::vector<Event> referenceEvents(const MorphologicalTree& tree, NodeId root) {
    std::vector<Event> pending{{'b', root, InvalidNode}}, result;
    while (!pending.empty()) {
        const auto event = pending.back();
        pending.pop_back();
        result.push_back(event);
        const auto [kind, node, child] = event;
        (void)child;
        if (kind == 'b') {
            pending.emplace_back('e', node, InvalidNode);
            std::vector<NodeId> children;
            for (NodeId childNode : tree.children(node)) children.push_back(childNode);
            for (auto it = children.rbegin(); it != children.rend(); ++it) {
                pending.emplace_back('m', node, *it);
                pending.emplace_back('b', *it, InvalidNode);
            }
        }
    }
    return result;
}

void checkLongBranchingTree() {
    std::vector<NodeId> parents(8197);
    for (NodeId node = 1; node < 8197; ++node) parents[node] = node - 1;
    for (NodeId parent : {4094, 4095, 4096, 4097, 8192}) parents.push_back(parent);
    std::vector<NodeId> owner(parents.size());
    std::iota(owner.begin(), owner.end(), 0);
    const std::vector<int> altitude(parents.size(), 0);
    const auto valued = makeValuedHierarchy<int>(parents, owner, altitude, 0, MorphologicalTreeSemantics{});
    for (NodeId root : {0, 4090}) {
        require(events(valued.topology(), root) == referenceEvents(valued.topology(), root),
                "callback order remains exact when long branching paths enter and leave the fixed ancestor cache");
    }
}

void checkDeepTree(bool isMax) {
    constexpr int size = 1'000'000;
    auto image = Image<float>::create(1, size);
    for (int i=0; i<size; ++i) image->rawData()[i] = static_cast<float>(i);
    auto valued = isMax ? MorphologicalTreeFactory::createMaxTree(image,1.0)
                        : MorphologicalTreeFactory::createMinTree(image,1.0);
    const auto& tree = valued.topology();
    requireEqual(tree.numInternalNodeSlots(),size,"ramp has one node per distinct level");
    const auto [names,values] = AttributeComputation::computeAttributes<double>(valued,
        {Area,Volume,RelativeVolume,MeanGrayLevel,GrayLevelVariance,GrayLevelHeight,SubtreeHeight,DepthNode,
         IsLeafNode,IsRootNode,NumChildrenNode,NumSiblingsNode,NumDescendantsNode,NumLeafDescendantsNode,
         LeafRatioNode,BalanceNode,AvgChildHeightNode});
    const auto at = [&](NodeId node, Attribute attribute) { return values[names.linearIndex(node,attribute)]; };
    for (NodeId node : tree.aliveNodeIds()) {
        const double h = valued.nodeAltitude(node);
        const double count = isMax ? size-h : h+1;
        const double mean = isMax ? (h+size-1)/2 : h/2;
        requireEqual(at(node,Area),count,"ramp support cardinality");
        requireEqual(at(node,Volume),count*mean,"ramp intensity sum");
        requireEqual(at(node,RelativeVolume),count*count,"ramp relative volume");
        requireEqual(at(node,MeanGrayLevel),mean,"ramp mean");
        requireEqual(at(node,GrayLevelHeight),count-1,"ramp level height");
        requireEqual(at(node,SubtreeHeight),count-1,"ramp subtree height");
        requireEqual(at(node,DepthNode),size-count,"ramp depth");
        requireEqual(at(node,NumChildrenNode),count>1 ? 1.0 : 0.0,"ramp child count");
        requireEqual(at(node,NumDescendantsNode),count-1,"ramp descendant count");
        require(std::isfinite(at(node,GrayLevelVariance)),"ramp variance remains finite");
    }
    const auto reconstructed = valued.reconstructFromNodeAltitudes();
    const auto area = AttributeComputation::computeSingleAttribute<double>(valued,Area);
    const auto filtered = AttributeFilters<float>(valued).filteringByPruningMax(area.second.data(),0.0);
    const std::vector<std::uint8_t> keep(static_cast<std::size_t>(size),1);
    const auto filteredMask = AttributeFilters<float>(valued).filteringByPruningMax(std::span<const std::uint8_t>(keep));
    for (int i=0; i<size; ++i) {
        requireEqual(reconstructed->rawData()[i],image->rawData()[i],"deep reconstruction identity");
        requireEqual(filtered->rawData()[i],image->rawData()[i],"deep attribute filtering with every node preserved");
        requireEqual(filteredMask->rawData()[i],image->rawData()[i],"deep mask filtering with every node preserved");
    }
}

} // namespace

int main() {
    checkCallbackOrder();
    checkLongBranchingTree();
    checkDeepTree(true);
    checkDeepTree(false);
}
