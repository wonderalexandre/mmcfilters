#pragma once

#include "AttributeComputer.hpp"
#include "AttributeComputedIncrementally.hpp"
#include "../trees/MorphologicalTree.hpp"


namespace mmcfilters {

namespace detail {
inline NodeId topologySlotOf(const MorphologicalTree&, NodeId nodeId) noexcept {
    return nodeId;
}
}

/**
 * @brief Computes structural descriptors that depend only on the tree
 * topology.
 *
 * @details
 * This computer ignores the image-domain geometry and uses only the parent /
 * child relationships encoded in the hierarchy. The descriptors cover:
 * - depth and height;
 * - leaf/root indicators;
 * - number of children, siblings, descendants, and leaf descendants;
 * - aggregated descriptors such as leaf ratio, subtree balance, and average
 *   child height.
 *
 * Some of these descriptors depend on intermediate quantities that are also
 * public outputs, such as subtree height or number of descendants. To avoid
 * forcing callers to request those attributes explicitly, the implementation
 * allocates temporary scratch buffers whenever an intermediate quantity is
 * needed internally but is not part of the requested output set.
 *
 * The computation is incremental:
 * - node-local quantities such as depth, root/leaf flags, and child counts
 *   are initialised on entry;
 * - bottom-up aggregation then propagates heights and descendant counts;
 * - a final visit converts the accumulated values into the derived ratios and
 *   balance measures.
 */
class TreeTopologyComputer : public AttributeComputer {
public:
    /**
     * @brief Returns the full family of topology-derived descriptors.
     */
    std::vector<Attribute> attributes() const override {
        return {HEIGHT_NODE, DEPTH_NODE, IS_LEAF_NODE, IS_ROOT_NODE, NUM_CHILDREN_NODE, NUM_SIBLINGS_NODE, NUM_DESCENDANTS_NODE, NUM_LEAF_DESCENDANTS_NODE, LEAF_RATIO_NODE, BALANCE_NODE, AVG_CHILD_HEIGHT_NODE};
    }

    /**
     * @brief Computes the requested topology descriptors.
     */
    void compute(MorphologicalTree& tree, const AltitudeBuffer*, std::span<float> buffer, const AttributeNames& attrNames, std::span<const Attribute> requestedAttributes, std::span<const DependencySource>) const override {
        if (PRINT_LOG) std::cout << "\n==== AttributeComputer: Computing STRUCTURE_TREE group" << std::endl;

        bool computeHeight = std::find(requestedAttributes.begin(), requestedAttributes.end(), HEIGHT_NODE) != requestedAttributes.end();
        bool computeDepth = std::find(requestedAttributes.begin(), requestedAttributes.end(), DEPTH_NODE) != requestedAttributes.end();
        bool computeIsLeaf = std::find(requestedAttributes.begin(), requestedAttributes.end(), IS_LEAF_NODE) != requestedAttributes.end();
        bool computeIsRoot = std::find(requestedAttributes.begin(), requestedAttributes.end(), IS_ROOT_NODE) != requestedAttributes.end();
        bool computeNumChildren = std::find(requestedAttributes.begin(), requestedAttributes.end(), NUM_CHILDREN_NODE) != requestedAttributes.end();
        bool computeNumSiblings = std::find(requestedAttributes.begin(), requestedAttributes.end(), NUM_SIBLINGS_NODE) != requestedAttributes.end();
        bool computeNumDescendants = std::find(requestedAttributes.begin(), requestedAttributes.end(), NUM_DESCENDANTS_NODE) != requestedAttributes.end();
        bool computeNumLeafDescendants = std::find(requestedAttributes.begin(), requestedAttributes.end(), NUM_LEAF_DESCENDANTS_NODE) != requestedAttributes.end();
        bool computeLeafRatio = std::find(requestedAttributes.begin(), requestedAttributes.end(), LEAF_RATIO_NODE) != requestedAttributes.end();
        bool computeBalance = std::find(requestedAttributes.begin(), requestedAttributes.end(), BALANCE_NODE) != requestedAttributes.end();
        bool computeAvgChildHeight = std::find(requestedAttributes.begin(), requestedAttributes.end(), AVG_CHILD_HEIGHT_NODE) != requestedAttributes.end();

        const int numNodeSlots = tree.getNumInternalNodeSlots();

        // Reuse the public output buffer when an intermediate quantity has been
        // requested explicitly; otherwise, fall back to temporary storage so
        // the derived descriptors can still be computed internally.
        std::vector<float> heightStorage(computeHeight ? 0 : numNodeSlots, 0.0f);
        float* bufferHeight = computeHeight ? buffer.data() : heightStorage.data();
        auto indexOfHeight = [&](NodeId idx) {
            return computeHeight ? attrNames.linearIndex(idx, HEIGHT_NODE) : idx;
        };

        std::vector<float> numDescStorage(computeNumDescendants ? 0 : numNodeSlots, 0.0f);
        float* bufferNumDesc = computeNumDescendants ? buffer.data() : numDescStorage.data();
        auto indexOfNumDescendants = [&](NodeId idx) {
            return computeNumDescendants ? attrNames.linearIndex(idx, NUM_DESCENDANTS_NODE) : idx;
        };

        std::vector<float> numLeafDescStorage(computeNumLeafDescendants ? 0 : numNodeSlots, 0.0f);
        float* bufferNumLeafDesc = computeNumLeafDescendants ? buffer.data() : numLeafDescStorage.data();
        auto indexOfNumLeafDescendants = [&](NodeId idx) {
            return computeNumLeafDescendants ? attrNames.linearIndex(idx, NUM_LEAF_DESCENDANTS_NODE) : idx;
        };

        std::vector<float> depthStorage(computeDepth ? 0 : numNodeSlots, 0.0f);
        float* bufferDepth = computeDepth ? buffer.data() : depthStorage.data();
        auto indexOfDepth = [&](NodeId idx) {
            return computeDepth ? attrNames.linearIndex(idx, DEPTH_NODE) : idx;
        };

        AttributeComputedIncrementally::traversePostOrder(tree,
            tree.getRoot(),
            [&](NodeId nodeId) {
                const NodeId node = detail::topologySlotOf(tree, nodeId);
                const bool isRoot = tree.isRoot(nodeId);
                const NodeId parentNodeId = tree.getNodeParent(nodeId);
                const NodeId parent = isRoot ? InvalidNode : detail::topologySlotOf(tree, parentNodeId);
                const int numChildren = tree.getNumChildren(nodeId);
                const bool isLeaf = tree.isLeaf(nodeId);

                float parentDepth = (parent != InvalidNode) ? bufferDepth[indexOfDepth(parent)] : InvalidNode;
                bufferDepth[indexOfDepth(node)] = parent != InvalidNode ? parentDepth + 1.0f : 0.0f;
                
                bufferHeight[indexOfHeight(node)] = 0.0f;
                bufferNumDesc[indexOfNumDescendants(node)] = 0.0f;
                bufferNumLeafDesc[indexOfNumLeafDescendants(node)] = isLeaf ? 1.0f : 0.0f;

                if (computeHeight)
                    buffer[attrNames.linearIndex(node, HEIGHT_NODE)] = 0.0f;
                if (computeIsLeaf)
                    buffer[attrNames.linearIndex(node, IS_LEAF_NODE)] = isLeaf ? 1.0f : 0.0f;
                if (computeIsRoot)
                    buffer[attrNames.linearIndex(node, IS_ROOT_NODE)] = isRoot ? 1.0f : 0.0f;
                if (computeNumChildren)
                    buffer[attrNames.linearIndex(node, NUM_CHILDREN_NODE)] = static_cast<float>(numChildren);
                if (computeNumSiblings)
                    buffer[attrNames.linearIndex(node, NUM_SIBLINGS_NODE)] = isRoot ? 0.0f : static_cast<float>(tree.getNumChildren(parentNodeId) - 1);
                if (computeLeafRatio)
                    buffer[attrNames.linearIndex(node, LEAF_RATIO_NODE)] = 0.0f;
                if (computeBalance)
                    buffer[attrNames.linearIndex(node, BALANCE_NODE)] = 0.0f;
                if (computeAvgChildHeight)
                    buffer[attrNames.linearIndex(node, AVG_CHILD_HEIGHT_NODE)] = 0.0f;
            },
            [&](NodeId parentNodeId, NodeId childNodeId) {
                const NodeId parent = detail::topologySlotOf(tree, parentNodeId);
                const NodeId child = detail::topologySlotOf(tree, childNodeId);
                
                bufferNumDesc[indexOfNumDescendants(parent)] += bufferNumDesc[indexOfNumDescendants(child)] + 1.0f;
                bufferNumLeafDesc[indexOfNumLeafDescendants(parent)] += bufferNumLeafDesc[indexOfNumLeafDescendants(child)];

                float childHeight = bufferHeight[indexOfHeight(child)];
                float& parentHeight = bufferHeight[indexOfHeight(parent)];
                parentHeight = std::max(parentHeight, childHeight + 1.0f);
                const int numChildren = tree.getNumChildren(parentNodeId);

                if (computeBalance) {
                    float& minH = buffer[attrNames.linearIndex(parent, BALANCE_NODE)];
                    if (numChildren == 1)
                        minH = childHeight;
                    else
                        minH = std::min(minH, childHeight);
                }

                if (computeAvgChildHeight) {
                    float& sumH = buffer[attrNames.linearIndex(parent, AVG_CHILD_HEIGHT_NODE)];
                    if (numChildren == 1)
                        sumH = childHeight;
                    else
                        sumH += childHeight;
                }
            },
            [&](NodeId idxGlobalId) {
                const NodeId idx = detail::topologySlotOf(tree, idxGlobalId);
                
                if (computeLeafRatio) {
                    float desc = bufferNumDesc[indexOfNumDescendants(idx)];
                    float folhas = bufferNumLeafDesc[indexOfNumLeafDescendants(idx)];
                    buffer[attrNames.linearIndex(idx, LEAF_RATIO_NODE)] = desc > 0.0f ? folhas / (desc + 1.0f) : 1.0f;
                }

                if (!tree.isLeaf(idxGlobalId)) {
                    if (computeBalance) {
                        float alturaMax = bufferHeight[indexOfHeight(idx)];
                        float alturaMin = buffer[attrNames.linearIndex(idx, BALANCE_NODE)];
                        buffer[attrNames.linearIndex(idx, BALANCE_NODE)] = alturaMax - alturaMin;
                    }

                    if (computeAvgChildHeight) {
                        buffer[attrNames.linearIndex(idx, AVG_CHILD_HEIGHT_NODE)] = buffer[attrNames.linearIndex(idx, AVG_CHILD_HEIGHT_NODE)] /
                                                                                     static_cast<float>(tree.getNumChildren(idxGlobalId));
                    }
                }
            }
        );
    }
};

} // namespace mmcfilters
