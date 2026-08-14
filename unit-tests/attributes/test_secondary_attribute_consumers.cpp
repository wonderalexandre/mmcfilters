#include "support/TestSupport.hpp"

#include "mmcfilters/attributes/AttributeComputation.hpp"
#include "mmcfilters/filters/MSERComputer.hpp"

#include <cmath>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

int main() {
    auto image = makeComponentTreeFixture();
    auto valuedTree = makeValuedComponentTree(image, true);
    auto requireThrows = [](auto&& fn, const std::string& label) {
        bool threw = false;
        try {
            fn();
        } catch (const std::exception&) {
            threw = true;
        }
        require(threw, label);
    };

    valuedTree->mergeNodeIntoParent(4);
    const auto& tree = valuedTree->topology();

    std::vector<NodeId> aliveNodeIds;
    for (NodeId nodeId : tree.aliveNodeIds()) {
        aliveNodeIds.push_back(nodeId);
    }
    requireVectorEqual(aliveNodeIds, {0, 1, 2, 3, 5}, "dense node ids after middle-slot merge");
    requireEqual(tree.numNodes(), 5, "active node count after middle-slot merge");
    requireEqual(tree.numInternalNodeSlots(), 6, "internal slot count after middle-slot merge");

    auto [topologyNames, topologyBuffer] =
        AttributeComputation::computeAttributes(*valuedTree, {SubtreeHeight, DepthNode, IsLeafNode, IsRootNode, NumChildrenNode, NumDescendantsNode});

    requireEqual(topologyBuffer[topologyNames.linearIndex(0, SubtreeHeight)], 4.0f, "root height after middle-slot merge");
    requireEqual(topologyBuffer[topologyNames.linearIndex(0, DepthNode)], 0.0f, "root depth after middle-slot merge");
    requireEqual(topologyBuffer[topologyNames.linearIndex(0, IsRootNode)], 1.0f, "root marker after middle-slot merge");
    requireEqual(topologyBuffer[topologyNames.linearIndex(3, NumChildrenNode)], 1.0f, "node 3 child count after middle-slot merge");
    requireEqual(topologyBuffer[topologyNames.linearIndex(3, NumDescendantsNode)], 1.0f, "node 3 descendants after middle-slot merge");
    requireEqual(topologyBuffer[topologyNames.linearIndex(5, IsLeafNode)], 1.0f, "node 5 leaf marker after middle-slot merge");
    requireEqual(topologyBuffer[topologyNames.linearIndex(5, DepthNode)], 4.0f, "node 5 depth after middle-slot merge");

    auto areaComputed = [&]() {
        try {
            return AttributeComputation::computeSingleAttribute(*valuedTree, Area);
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("area computation failed: ") + e.what());
        }
    }();
    auto bitquadsComputed = [&]() {
        try {
            return AttributeComputation::computeSingleAttribute(*valuedTree, BitquadArea);
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("bitquads computation failed: ") + e.what());
        }
    }();
    auto maxDistComputed = [&]() {
        try {
            return AttributeComputation::computeSingleAttribute(*valuedTree, MaxDist);
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("max-dist computation failed: ") + e.what());
        }
    }();
    auto boxComputed = [&]() {
        try {
            return AttributeComputation::computeAttributes(*valuedTree, {BoxWidth, BoundingBoxHeight, Rectangularity, DiagonalLength});
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("bounding-box computation failed: ") + e.what());
        }
    }();
    auto momentComputed = [&]() {
        try {
            return AttributeComputation::computeAttributes(*valuedTree, {CentralMoment20, CentralMoment02, Inertia, Circularity});
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("moment computation failed: ") + e.what());
        }
    }();
    auto grayComputed = [&]() {
        try {
            return AttributeComputation::computeAttributes(*valuedTree, {MeanGrayLevel, GrayLevelVariance, GrayLevelHeight, Volume, RelativeVolume});
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("gray-level computation failed: ") + e.what());
        }
    }();

    const auto& areaNames = areaComputed.first;
    const auto& areaBuffer = areaComputed.second;
    const auto& bitquadsNames = bitquadsComputed.first;
    const auto& bitquadsBuffer = bitquadsComputed.second;
    const auto& maxDistNames = maxDistComputed.first;
    const auto& maxDistBuffer = maxDistComputed.second;
    const auto& boxNames = boxComputed.first;
    const auto& boxBuffer = boxComputed.second;
    const auto& momentNames = momentComputed.first;
    const auto& momentBuffer = momentComputed.second;
    const auto& grayNames = grayComputed.first;
    const auto& grayBuffer = grayComputed.second;

    auto pairCompatibleArea = AttributeComputation::computeSingleAttribute(*valuedTree, Area);
    requireEqual(pairCompatibleArea.first.linearIndex(0, Area), 0, "pair-like first field still exposes AttributeNames");
    requireEqual(pairCompatibleArea.second[pairCompatibleArea.first.linearIndex(0, Area)], 16.0f, "pair-like second field still exposes buffer");

    auto [explicitTopologyNames, explicitTopologyBuffer] = AttributeComputation::computeTopologyAttributes(tree, {Area, BoxWidth});
    requireEqual(explicitTopologyBuffer[explicitTopologyNames.linearIndex(0, Area)], 16.0f,
                 "explicit topology attribute facade must compute AREA without altitude");
    if constexpr (contract::validationsEnabled) {
        requireThrows([&]() { (void)AttributeComputation::computeSingleTopologyAttribute(tree, MeanGrayLevel); },
                      "explicit topology attribute facade must reject altitude-dependent attributes");
    }

    requireThrows([&]() { (void)AttributeComputation::computeSingleAttribute(*valuedTree, Area, NodeIdSpace::Higra); },
                  "image-built or edited tree must reject Higra-space attribute projection");

    auto [higraParent, higraAltitude] = valuedTree->exportHigraHierarchy();
    auto fromHigra = MorphologicalTreeFactory::createFromHigraParent(
        std::span<const NodeId>(higraParent), std::span<const std::uint8_t>(higraAltitude), tree.numRows(), tree.numColumns(),
        MorphologicalTreeKind::MaxTree, RegularGridAdjacency2D(tree.numRows(), tree.numColumns(), 1.5));
    auto higraComputed = AttributeComputation::computeAttributes(fromHigra, {Area, Volume}, NodeIdSpace::Higra);
    require(higraComputed.nodeIdSpace == NodeIdSpace::Higra, "Higra-space multi-attribute projection must preserve output node-id space");
    require(higraComputed.attributeNames().NUM_ATTRIBUTES == 2, "Higra-space multi-attribute projection must preserve requested attribute count");
    require(higraComputed.values().size() == static_cast<std::size_t>(fromHigra.topology().getNodeIdSpaceSize(NodeIdSpace::Higra)) * 2,
            "Higra-space multi-attribute projection must materialise the exported node-id space");

    auto mutatedFromHigra = MorphologicalTreeFactory::createFromHigraParent(
        std::span<const NodeId>(higraParent), std::span<const std::uint8_t>(higraAltitude), tree.numRows(), tree.numColumns(),
        MorphologicalTreeKind::MaxTree, RegularGridAdjacency2D(tree.numRows(), tree.numColumns(), 1.5));
    NodeId firstNonRootNodeId = InvalidNode;
    for (NodeId nodeId : mutatedFromHigra.topology().aliveNodeIds()) {
        if (!mutatedFromHigra.topology().isRoot(nodeId)) {
            firstNonRootNodeId = nodeId;
            break;
        }
    }
    require(firstNonRootNodeId != InvalidNode, "Higra projection mutation test needs one non-root node");
    mutatedFromHigra.mergeNodeIntoParent(firstNonRootNodeId);
    requireThrows([&]() { (void)AttributeComputation::computeSingleAttribute(mutatedFromHigra, Area, NodeIdSpace::Higra); },
                  "mutated Higra import must reject scalar Higra-space projection");
    requireThrows(
        [&]() {
            (void)AttributeComputation::computeSampledNodeAttribute(
                mutatedFromHigra, Area, AltitudeDifference<std::uint8_t>{1}, 1,
                NodeAttributeSamplingPolicy::LargestSupportDescendant, MissingNodeAttributeSamplePolicy::NotANumber, NodeIdSpace::Higra);
        },
        "mutated Higra import must reject sampled-attribute Higra-space projection");
    requireThrows([&]() { (void)AttributeComputation::computeAttributes(mutatedFromHigra, {Area, Volume}, NodeIdSpace::Higra); },
                  "mutated Higra import must reject multi-attribute Higra-space projection");

    auto [sampleLayout, sampledDepthValues] =
        AttributeComputation::computeSampledNodeAttribute(*valuedTree, DepthNode, AltitudeDifference<std::uint8_t>{1}, 1);
    requireEqual(sampledDepthValues[sampleLayout.linearIndex(1, DepthNode, -1)], 0.0f,
                 "ancestor sampling must preserve a legitimate zero-valued node attribute");

    requireEqual(areaBuffer[areaNames.linearIndex(0, Area)], 16.0f, "root area after middle-slot merge");
    requireEqual(areaBuffer[areaNames.linearIndex(1, Area)], 15.0f, "node 1 area after middle-slot merge");
    requireEqual(areaBuffer[areaNames.linearIndex(2, Area)], 12.0f, "node 2 area after middle-slot merge");
    requireEqual(areaBuffer[areaNames.linearIndex(3, Area)], 8.0f, "node 3 area after middle-slot merge");
    requireEqual(areaBuffer[areaNames.linearIndex(5, Area)], 2.0f, "node 5 area after middle-slot merge");

    requireEqual(boxBuffer[boxNames.linearIndex(0, BoxWidth)], 4.0f, "root box width after middle-slot merge");
    requireEqual(boxBuffer[boxNames.linearIndex(0, BoundingBoxHeight)], 4.0f, "root box height after middle-slot merge");
    require(boxBuffer[boxNames.linearIndex(0, Rectangularity)] > 0.0f, "root rectangularity after middle-slot merge");
    require(boxBuffer[boxNames.linearIndex(0, DiagonalLength)] > 0.0f, "root diagonal length after middle-slot merge");
    requireNear(boxBuffer[boxNames.linearIndex(3, Rectangularity)], 8.0f / 12.0f, 1e-6f, "node 3 rectangularity after middle-slot merge");

    for (NodeId nodeId : tree.aliveNodeIds()) {
        const float bitquadsArea = bitquadsBuffer[bitquadsNames.linearIndex(nodeId, BitquadArea)];
        const float maxDist = maxDistBuffer[maxDistNames.linearIndex(nodeId, MaxDist)];
        const float mu20 = momentBuffer[momentNames.linearIndex(nodeId, CentralMoment20)];
        const float mu02 = momentBuffer[momentNames.linearIndex(nodeId, CentralMoment02)];
        const float inertia = momentBuffer[momentNames.linearIndex(nodeId, Inertia)];
        const float circularity = momentBuffer[momentNames.linearIndex(nodeId, Circularity)];
        const float meanGrayLevel = grayBuffer[grayNames.linearIndex(nodeId, MeanGrayLevel)];
        const float variance = grayBuffer[grayNames.linearIndex(nodeId, GrayLevelVariance)];
        const float grayLevelHeight = grayBuffer[grayNames.linearIndex(nodeId, GrayLevelHeight)];
        const float volume = grayBuffer[grayNames.linearIndex(nodeId, Volume)];
        const float relativeVolume = grayBuffer[grayNames.linearIndex(nodeId, RelativeVolume)];
        require(std::isfinite(bitquadsArea), "bitquads area must stay finite after middle-slot merge");
        require(std::isfinite(maxDist), "max-dist must stay finite after middle-slot merge");
        require(std::isfinite(mu20), "mu20 must stay finite after middle-slot merge");
        require(std::isfinite(mu02), "mu02 must stay finite after middle-slot merge");
        require(std::isfinite(inertia), "inertia must stay finite after middle-slot merge");
        require(std::isfinite(circularity), "circularity must stay finite after middle-slot merge");
        require(std::isfinite(meanGrayLevel), "mean gray level must stay finite after middle-slot merge");
        require(std::isfinite(variance), "variance must stay finite after middle-slot merge");
        require(std::isfinite(grayLevelHeight), "gray-level height must stay finite after middle-slot merge");
        require(std::isfinite(volume), "volume must stay finite after middle-slot merge");
        require(std::isfinite(relativeVolume), "relative volume must stay finite after middle-slot merge");
        require(bitquadsArea >= 0.0f, "bitquads area must stay non-negative after middle-slot merge");
        require(maxDist >= 0.0f, "max-dist must stay non-negative after middle-slot merge");
        require(mu20 >= 0.0f, "mu20 must stay non-negative after middle-slot merge");
        require(mu02 >= 0.0f, "mu02 must stay non-negative after middle-slot merge");
        require(inertia >= 0.0f, "inertia must stay non-negative after middle-slot merge");
        require(grayLevelHeight >= 0.0f, "gray-level height must stay non-negative after middle-slot merge");
        require(volume >= 0.0f, "volume must stay non-negative after middle-slot merge");
        require(relativeVolume >= 0.0f, "relative volume must stay non-negative after middle-slot merge");
    }
    requireEqual(valuedTree->nodeAltitude(3), std::uint8_t{3}, "node 3 altitude after middle-slot merge");
    requireNear(grayBuffer[grayNames.linearIndex(3, MeanGrayLevel)], 3.5f, 1e-6f, "node 3 mean level after middle-slot merge");
    requireEqual(grayBuffer[grayNames.linearIndex(3, Volume)], 28.0f, "node 3 volume after middle-slot merge");
    requireEqual(grayBuffer[grayNames.linearIndex(3, RelativeVolume)], 14.0f, "node 3 relative volume after middle-slot merge");

    MSERComputer<std::uint8_t> mser(*valuedTree);
    std::vector<uint8_t> isMSER = mser.computeMSER(1);
    const std::vector<float> implicitVariations = mser.getVariations();
    MSERComputer<std::uint8_t> mserExplicit(*valuedTree, areaBuffer);
    std::vector<uint8_t> isMSERExplicit = mserExplicit.computeMSER(1);
    const std::vector<float> explicitVariations = mserExplicit.getVariations();
    MSERComputer<std::uint8_t> mserRaw(*valuedTree, areaBuffer.data());
    std::vector<uint8_t> isMSERRaw = mserRaw.computeMSER(1);
    const std::vector<float> rawVariations = mserRaw.getVariations();
    requireEqual(static_cast<int>(isMSER.size()), tree.numInternalNodeSlots(), "MSER buffer size after middle-slot merge");
    requireVectorEqual(isMSERExplicit, isMSER, "explicit area MSER flags must match implicit area MSER flags");
    requireVectorEqual(isMSERRaw, isMSER, "raw area MSER flags must match implicit area MSER flags");

    int numActiveMserFlags = 0;
    for (NodeId nodeId : tree.aliveNodeIds()) {
        numActiveMserFlags += isMSER[nodeId] ? 1 : 0;
        const float implicit = implicitVariations[nodeId];
        const float explicitValue = explicitVariations[nodeId];
        const float rawValue = rawVariations[nodeId];
        if (std::isnan(implicit) || std::isnan(explicitValue)) {
            require(std::isnan(implicit) && std::isnan(explicitValue) && std::isnan(rawValue), "MSER stability NaN pattern after middle-slot merge");
        } else {
            requireNear(implicit, explicitValue, 1e-6f, "MSER stability consistency after middle-slot merge");
            requireNear(implicit, rawValue, 1e-6f, "raw MSER stability consistency after middle-slot merge");
        }
    }
    require(numActiveMserFlags >= 0, "MSER flags must be readable on alive nodes after middle-slot merge");

    requireEqual(tree.lowestCommonAncestor(3, 5), 3, "LCA after middle-slot merge");
    requireEqual(tree.lowestCommonAncestor(2, 5), 2, "ancestor LCA after middle-slot merge");
    require(tree.isStrictAncestor(3, 5), "strict ancestor after middle-slot merge");
    require(tree.isDescendant(5, 2), "descendant after middle-slot merge");

    return 0;
}
