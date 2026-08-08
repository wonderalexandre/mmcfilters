#include "support/TestSupport.hpp"

#include "mmcfilters/attributes/AttributeComputation.hpp"
#include "mmcfilters/filters/MSERComputer.hpp"

#include <cmath>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

int main() {
    auto image = makeComponentTreeFixture();
    auto weighted = makeWeightedComponentTree(image, true);
    auto requireThrows = [](auto&& fn, const std::string& label) {
        bool threw = false;
        try {
            fn();
        } catch (const std::exception&) {
            threw = true;
        }
        require(threw, label);
    };

    weighted->mergeNodeIntoParent(4);
    const auto& tree = weighted->topology();

    std::vector<NodeId> aliveNodeIds;
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        aliveNodeIds.push_back(nodeId);
    }
    requireVectorEqual(aliveNodeIds, {0, 1, 2, 3, 5}, "dense node ids after middle-slot merge");
    requireEqual(tree.getNumNodes(), 5, "active node count after middle-slot merge");
    requireEqual(tree.getNumInternalNodeSlots(), 6, "internal slot count after middle-slot merge");

    auto [topologyNames, topologyBuffer] =
        AttributeComputation::computeAttributes(*weighted, {HEIGHT_NODE, DEPTH_NODE, IS_LEAF_NODE, IS_ROOT_NODE, NUM_CHILDREN_NODE, NUM_DESCENDANTS_NODE});

    requireEqual(topologyBuffer[topologyNames.linearIndex(0, HEIGHT_NODE)], 4.0f, "root height after middle-slot merge");
    requireEqual(topologyBuffer[topologyNames.linearIndex(0, DEPTH_NODE)], 0.0f, "root depth after middle-slot merge");
    requireEqual(topologyBuffer[topologyNames.linearIndex(0, IS_ROOT_NODE)], 1.0f, "root marker after middle-slot merge");
    requireEqual(topologyBuffer[topologyNames.linearIndex(3, NUM_CHILDREN_NODE)], 1.0f, "node 3 child count after middle-slot merge");
    requireEqual(topologyBuffer[topologyNames.linearIndex(3, NUM_DESCENDANTS_NODE)], 1.0f, "node 3 descendants after middle-slot merge");
    requireEqual(topologyBuffer[topologyNames.linearIndex(5, IS_LEAF_NODE)], 1.0f, "node 5 leaf marker after middle-slot merge");
    requireEqual(topologyBuffer[topologyNames.linearIndex(5, DEPTH_NODE)], 4.0f, "node 5 depth after middle-slot merge");

    auto areaComputed = [&]() {
        try {
            return AttributeComputation::computeSingleAttribute(*weighted, AREA);
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("area computation failed: ") + e.what());
        }
    }();
    auto bitquadsComputed = [&]() {
        try {
            return AttributeComputation::computeSingleAttribute(*weighted, BITQUADS_AREA);
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("bitquads computation failed: ") + e.what());
        }
    }();
    auto maxDistComputed = [&]() {
        try {
            return AttributeComputation::computeSingleAttribute(*weighted, MAX_DIST);
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("max-dist computation failed: ") + e.what());
        }
    }();
    auto boxComputed = [&]() {
        try {
            return AttributeComputation::computeAttributes(*weighted, {BOX_WIDTH, BOX_HEIGHT, RECTANGULARITY, DIAGONAL_LENGTH});
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("bounding-box computation failed: ") + e.what());
        }
    }();
    auto momentComputed = [&]() {
        try {
            return AttributeComputation::computeAttributes(*weighted, {CENTRAL_MOMENT_20, CENTRAL_MOMENT_02, INERTIA, CIRCULARITY});
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("moment computation failed: ") + e.what());
        }
    }();
    auto grayComputed = [&]() {
        try {
            return AttributeComputation::computeAttributes(*weighted, {LEVEL, MEAN_LEVEL, VARIANCE_LEVEL, GRAY_HEIGHT, VOLUME, RELATIVE_VOLUME});
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

    auto pairCompatibleArea = AttributeComputation::computeSingleAttribute(*weighted, AREA);
    requireEqual(pairCompatibleArea.first.linearIndex(0, AREA), 0, "pair-like first field still exposes AttributeNames");
    requireEqual(pairCompatibleArea.second[pairCompatibleArea.first.linearIndex(0, AREA)], 16.0f, "pair-like second field still exposes buffer");

    auto [explicitTopologyNames, explicitTopologyBuffer] = AttributeComputation::computeTopologyAttributes(tree, {AREA, BOX_WIDTH});
    requireEqual(explicitTopologyBuffer[explicitTopologyNames.linearIndex(0, AREA)], 16.0f,
                 "explicit topology attribute facade must compute AREA without altitude");
    requireThrows([&]() { (void)AttributeComputation::computeSingleTopologyAttribute(tree, LEVEL); },
                  "explicit topology attribute facade must reject altitude-dependent attributes");

    requireThrows([&]() { (void)AttributeComputation::computeSingleAttribute(*weighted, AREA, NodeIdSpace::HIGRA); },
                  "image-built or edited tree must reject Higra-space attribute projection");

    auto [higraParent, higraAltitude] = weighted->exportHigraHierarchy();
    auto fromHigra = MorphologicalTreeFactory::createFromHigraParent(
        std::span<const NodeId>(higraParent), std::span<const std::uint8_t>(higraAltitude), tree.getNumRowsOfGridDomain2D(), tree.getNumColsOfGridDomain2D(),
        MorphologicalTreeKind::MAX_TREE, RegularGridAdjacency2D(tree.getNumRowsOfGridDomain2D(), tree.getNumColsOfGridDomain2D(), 1.5));
    auto higraComputed = AttributeComputation::computeAttributes(fromHigra, {AREA, VOLUME}, NodeIdSpace::HIGRA);
    require(higraComputed.nodeIdSpace == NodeIdSpace::HIGRA, "Higra-space multi-attribute projection must preserve output node-id space");
    require(higraComputed.attributeNames().NUM_ATTRIBUTES == 2, "Higra-space multi-attribute projection must preserve requested attribute count");
    require(higraComputed.values().size() == static_cast<std::size_t>(fromHigra.topology().getNodeIdSpaceSize(NodeIdSpace::HIGRA)) * 2,
            "Higra-space multi-attribute projection must materialise the exported node-id space");

    auto mutatedFromHigra = MorphologicalTreeFactory::createFromHigraParent(
        std::span<const NodeId>(higraParent), std::span<const std::uint8_t>(higraAltitude), tree.getNumRowsOfGridDomain2D(), tree.getNumColsOfGridDomain2D(),
        MorphologicalTreeKind::MAX_TREE, RegularGridAdjacency2D(tree.getNumRowsOfGridDomain2D(), tree.getNumColsOfGridDomain2D(), 1.5));
    NodeId firstNonRootNodeId = InvalidNode;
    for (NodeId nodeId : mutatedFromHigra.topology().getAliveNodeIds()) {
        if (!mutatedFromHigra.topology().isRoot(nodeId)) {
            firstNonRootNodeId = nodeId;
            break;
        }
    }
    require(firstNonRootNodeId != InvalidNode, "Higra projection mutation test needs one non-root node");
    mutatedFromHigra.mergeNodeIntoParent(firstNonRootNodeId);
    requireThrows([&]() { (void)AttributeComputation::computeSingleAttribute(mutatedFromHigra, AREA, NodeIdSpace::HIGRA); },
                  "mutated Higra import must reject scalar Higra-space projection");
    requireThrows(
        [&]() {
            (void)AttributeComputation::computeSingleAttributeWithDelta(mutatedFromHigra, AREA, AltitudeDiff<std::uint8_t>{1}, 1, "null-padding",
                                                                        NodeIdSpace::HIGRA);
        },
        "mutated Higra import must reject delta Higra-space projection");
    requireThrows([&]() { (void)AttributeComputation::computeAttributes(mutatedFromHigra, {AREA, VOLUME}, NodeIdSpace::HIGRA); },
                  "mutated Higra import must reject multi-attribute Higra-space projection");

    auto [deltaLevelNames, deltaLevelBuffer] =
        AttributeComputation::computeSingleAttributeWithDelta(*weighted, LEVEL, AltitudeDiff<std::uint8_t>{1}, 1, "last-padding");
    requireEqual(deltaLevelBuffer[deltaLevelNames.linearIndex(1, LEVEL, -1)], 0.0f, "delta padding must preserve a legitimate zero-valued ascendant");

    requireEqual(areaBuffer[areaNames.linearIndex(0, AREA)], 16.0f, "root area after middle-slot merge");
    requireEqual(areaBuffer[areaNames.linearIndex(1, AREA)], 15.0f, "node 1 area after middle-slot merge");
    requireEqual(areaBuffer[areaNames.linearIndex(2, AREA)], 12.0f, "node 2 area after middle-slot merge");
    requireEqual(areaBuffer[areaNames.linearIndex(3, AREA)], 8.0f, "node 3 area after middle-slot merge");
    requireEqual(areaBuffer[areaNames.linearIndex(5, AREA)], 2.0f, "node 5 area after middle-slot merge");

    requireEqual(boxBuffer[boxNames.linearIndex(0, BOX_WIDTH)], 4.0f, "root box width after middle-slot merge");
    requireEqual(boxBuffer[boxNames.linearIndex(0, BOX_HEIGHT)], 4.0f, "root box height after middle-slot merge");
    require(boxBuffer[boxNames.linearIndex(0, RECTANGULARITY)] > 0.0f, "root rectangularity after middle-slot merge");
    require(boxBuffer[boxNames.linearIndex(0, DIAGONAL_LENGTH)] > 0.0f, "root diagonal length after middle-slot merge");
    requireNear(boxBuffer[boxNames.linearIndex(3, RECTANGULARITY)], 8.0f / 12.0f, 1e-6f, "node 3 rectangularity after middle-slot merge");

    for (NodeId nodeId : tree.getAliveNodeIds()) {
        const float bitquadsArea = bitquadsBuffer[bitquadsNames.linearIndex(nodeId, BITQUADS_AREA)];
        const float maxDist = maxDistBuffer[maxDistNames.linearIndex(nodeId, MAX_DIST)];
        const float mu20 = momentBuffer[momentNames.linearIndex(nodeId, CENTRAL_MOMENT_20)];
        const float mu02 = momentBuffer[momentNames.linearIndex(nodeId, CENTRAL_MOMENT_02)];
        const float inertia = momentBuffer[momentNames.linearIndex(nodeId, INERTIA)];
        const float circularity = momentBuffer[momentNames.linearIndex(nodeId, CIRCULARITY)];
        const float level = grayBuffer[grayNames.linearIndex(nodeId, LEVEL)];
        const float meanLevel = grayBuffer[grayNames.linearIndex(nodeId, MEAN_LEVEL)];
        const float variance = grayBuffer[grayNames.linearIndex(nodeId, VARIANCE_LEVEL)];
        const float grayHeight = grayBuffer[grayNames.linearIndex(nodeId, GRAY_HEIGHT)];
        const float volume = grayBuffer[grayNames.linearIndex(nodeId, VOLUME)];
        const float relativeVolume = grayBuffer[grayNames.linearIndex(nodeId, RELATIVE_VOLUME)];
        require(std::isfinite(bitquadsArea), "bitquads area must stay finite after middle-slot merge");
        require(std::isfinite(maxDist), "max-dist must stay finite after middle-slot merge");
        require(std::isfinite(mu20), "mu20 must stay finite after middle-slot merge");
        require(std::isfinite(mu02), "mu02 must stay finite after middle-slot merge");
        require(std::isfinite(inertia), "inertia must stay finite after middle-slot merge");
        require(std::isfinite(circularity), "circularity must stay finite after middle-slot merge");
        require(std::isfinite(level), "level must stay finite after middle-slot merge");
        require(std::isfinite(meanLevel), "mean level must stay finite after middle-slot merge");
        require(std::isfinite(variance), "variance must stay finite after middle-slot merge");
        require(std::isfinite(grayHeight), "gray height must stay finite after middle-slot merge");
        require(std::isfinite(volume), "volume must stay finite after middle-slot merge");
        require(std::isfinite(relativeVolume), "relative volume must stay finite after middle-slot merge");
        require(bitquadsArea >= 0.0f, "bitquads area must stay non-negative after middle-slot merge");
        require(maxDist >= 0.0f, "max-dist must stay non-negative after middle-slot merge");
        require(mu20 >= 0.0f, "mu20 must stay non-negative after middle-slot merge");
        require(mu02 >= 0.0f, "mu02 must stay non-negative after middle-slot merge");
        require(inertia >= 0.0f, "inertia must stay non-negative after middle-slot merge");
        require(grayHeight >= 0.0f, "gray height must stay non-negative after middle-slot merge");
        require(volume >= 0.0f, "volume must stay non-negative after middle-slot merge");
        require(relativeVolume >= 0.0f, "relative volume must stay non-negative after middle-slot merge");
    }
    requireEqual(grayBuffer[grayNames.linearIndex(3, LEVEL)], 3.0f, "node 3 level after middle-slot merge");
    requireNear(grayBuffer[grayNames.linearIndex(3, MEAN_LEVEL)], 3.5f, 1e-6f, "node 3 mean level after middle-slot merge");
    requireEqual(grayBuffer[grayNames.linearIndex(3, VOLUME)], 28.0f, "node 3 volume after middle-slot merge");
    requireEqual(grayBuffer[grayNames.linearIndex(3, RELATIVE_VOLUME)], 14.0f, "node 3 relative volume after middle-slot merge");

    MSERComputer<std::uint8_t> mser(*weighted);
    std::vector<uint8_t> isMSER = mser.computeMSER(1);
    const std::vector<float> implicitVariations = mser.getVariations();
    MSERComputer<std::uint8_t> mserExplicit(*weighted, areaBuffer);
    std::vector<uint8_t> isMSERExplicit = mserExplicit.computeMSER(1);
    const std::vector<float> explicitVariations = mserExplicit.getVariations();
    MSERComputer<std::uint8_t> mserRaw(*weighted, areaBuffer.data());
    std::vector<uint8_t> isMSERRaw = mserRaw.computeMSER(1);
    const std::vector<float> rawVariations = mserRaw.getVariations();
    requireEqual(static_cast<int>(isMSER.size()), tree.getNumInternalNodeSlots(), "MSER buffer size after middle-slot merge");
    requireVectorEqual(isMSERExplicit, isMSER, "explicit area MSER flags must match implicit area MSER flags");
    requireVectorEqual(isMSERRaw, isMSER, "raw area MSER flags must match implicit area MSER flags");

    int numActiveMserFlags = 0;
    for (NodeId nodeId : tree.getAliveNodeIds()) {
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

    requireEqual(tree.getLowestCommonAncestor(3, 5), 3, "LCA after middle-slot merge");
    requireEqual(tree.getLowestCommonAncestor(2, 5), 2, "ancestor LCA after middle-slot merge");
    require(tree.isStrictAncestor(3, 5), "strict ancestor after middle-slot merge");
    require(tree.isDescendant(5, 2), "descendant after middle-slot merge");

    return 0;
}
