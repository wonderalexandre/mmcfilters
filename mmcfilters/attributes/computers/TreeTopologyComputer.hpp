#pragma once

#include "AttributeComputerDomain.hpp"
#include "AttributeComputerFamily.hpp"
#include "../detail/AttributeKernelSupport.hpp"
#include "../../trees/detail/TreeTraversalDetail.hpp"
#include "../../trees/MorphologicalTree.hpp"

#include <algorithm>
#include <array>
#include <concepts>
#include <limits>
#include <string_view>
#include <vector>

namespace mmcfilters::attributes::computers {

namespace detail {
/**
 * @brief Returns the dense attribute-buffer slot for a tree node.
 *
 * @param nodeId Identifier of the node used by the operation.
 * @return Dense attribute-buffer slot for the node.
 */
inline NodeId topologySlotOf(const MorphologicalTree&, NodeId nodeId) noexcept { return nodeId; }
} // namespace detail

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
class TreeTopologyComputer {
  public:
    /// Family name used in dependency-plan diagnostics.
    static constexpr std::string_view familyName = "tree-topology";

    /// Stable family id used by the scheduler.
    static constexpr AttributeComputerFamily family = AttributeComputerFamily::TreeTopology;

    /// Execution domain required by the computer.
    static constexpr AttributeComputerDomain domain = AttributeComputerDomain::Topology;

    /**
     * @brief Canonical list of topology-derived descriptors produced by this computer.
     */
    inline static constexpr std::array<Attribute, 11> producedAttributes{HEIGHT_NODE,       DEPTH_NODE,        IS_LEAF_NODE,         IS_ROOT_NODE,
                                                                         NUM_CHILDREN_NODE, NUM_SIBLINGS_NODE, NUM_DESCENDANTS_NODE, NUM_LEAF_DESCENDANTS_NODE,
                                                                         LEAF_RATIO_NODE,   BALANCE_NODE,      AVG_CHILD_HEIGHT_NODE};

    /**
     * @brief Computes the requested topology descriptors.
     *
     * @details
     * The computation is topology-only and ignores altitude/dependencies. All
     * output rows are indexed by dense internal node id. Intermediate quantities
     * such as height, depth, and descendant counts are kept in scratch storage
     * when they are needed to derive a requested descriptor but are not
     * themselves requested as public columns.
     *
     * @param context Operation context or diagnostic label.
     */
    template <std::floating_point Real> static void compute(const AttributeComputeContext<Real>& context) {
        computeImpl(context.tree, context.buffer, context.attrNames, context.requestedAttributes);
    }

  private:
    /**
     * @brief Computes the requested attribute values into the output buffer.
     *
     * @param tree Tree topology used by the operation.
     * @param buffer Buffer read or written by the operation.
     * @param attrNames Layout mapping attributes to buffer columns.
     * @param requestedAttributes Requested attribute subset.
     */
    template <std::floating_point Real>
    static void computeImpl(const MorphologicalTree& tree, std::span<Real> buffer, const AttributeNames& attrNames,
                            std::span<const Attribute> requestedAttributes) {
        requireAttributeBufferShape(tree, buffer, attrNames);

        bool computeHeight = std::find(requestedAttributes.begin(), requestedAttributes.end(), HEIGHT_NODE) != requestedAttributes.end();
        bool computeDepth = std::find(requestedAttributes.begin(), requestedAttributes.end(), DEPTH_NODE) != requestedAttributes.end();
        bool computeIsLeaf = std::find(requestedAttributes.begin(), requestedAttributes.end(), IS_LEAF_NODE) != requestedAttributes.end();
        bool computeIsRoot = std::find(requestedAttributes.begin(), requestedAttributes.end(), IS_ROOT_NODE) != requestedAttributes.end();
        bool computeNumChildren = std::find(requestedAttributes.begin(), requestedAttributes.end(), NUM_CHILDREN_NODE) != requestedAttributes.end();
        bool computeNumSiblings = std::find(requestedAttributes.begin(), requestedAttributes.end(), NUM_SIBLINGS_NODE) != requestedAttributes.end();
        bool computeNumDescendants = std::find(requestedAttributes.begin(), requestedAttributes.end(), NUM_DESCENDANTS_NODE) != requestedAttributes.end();
        bool computeNumLeafDescendants =
            std::find(requestedAttributes.begin(), requestedAttributes.end(), NUM_LEAF_DESCENDANTS_NODE) != requestedAttributes.end();
        bool computeLeafRatio = std::find(requestedAttributes.begin(), requestedAttributes.end(), LEAF_RATIO_NODE) != requestedAttributes.end();
        bool computeBalance = std::find(requestedAttributes.begin(), requestedAttributes.end(), BALANCE_NODE) != requestedAttributes.end();
        bool computeAvgChildHeight = std::find(requestedAttributes.begin(), requestedAttributes.end(), AVG_CHILD_HEIGHT_NODE) != requestedAttributes.end();

        const int numNodeSlots = tree.getNumInternalNodeSlots();

        // Reuse the public output buffer when an intermediate quantity has been
        // requested explicitly; otherwise, fall back to temporary storage so
        // the derived descriptors can still be computed internally.
        std::vector<Real> heightStorage(computeHeight ? 0 : numNodeSlots, Real{0});
        Real* bufferHeight = computeHeight ? buffer.data() : heightStorage.data();
        auto indexOfHeight = [&](NodeId idx) { return computeHeight ? attrNames.linearIndex(idx, HEIGHT_NODE) : idx; };

        std::vector<Real> numDescStorage(computeNumDescendants ? 0 : numNodeSlots, Real{0});
        Real* bufferNumDesc = computeNumDescendants ? buffer.data() : numDescStorage.data();
        auto indexOfNumDescendants = [&](NodeId idx) { return computeNumDescendants ? attrNames.linearIndex(idx, NUM_DESCENDANTS_NODE) : idx; };

        std::vector<Real> numLeafDescStorage(computeNumLeafDescendants ? 0 : numNodeSlots, Real{0});
        Real* bufferNumLeafDesc = computeNumLeafDescendants ? buffer.data() : numLeafDescStorage.data();
        auto indexOfNumLeafDescendants = [&](NodeId idx) { return computeNumLeafDescendants ? attrNames.linearIndex(idx, NUM_LEAF_DESCENDANTS_NODE) : idx; };

        std::vector<Real> depthStorage(computeDepth ? 0 : numNodeSlots, Real{0});
        Real* bufferDepth = computeDepth ? buffer.data() : depthStorage.data();
        auto indexOfDepth = [&](NodeId idx) { return computeDepth ? attrNames.linearIndex(idx, DEPTH_NODE) : idx; };

        std::vector<Real> minChildHeightStorage(computeBalance ? numNodeSlots : 0, Real{0});

        ::mmcfilters::detail::traversePostOrder(
            tree, tree.getRoot(),
            [&](NodeId nodeId) {
                const NodeId node = detail::topologySlotOf(tree, nodeId);
                const bool isRoot = tree.isRoot(nodeId);
                const NodeId parentNodeId = tree.getNodeParent(nodeId);
                const NodeId parent = isRoot ? InvalidNode : detail::topologySlotOf(tree, parentNodeId);
                const int numChildren = tree.getNumChildren(nodeId);
                const bool isLeaf = tree.isLeaf(nodeId);

                Real parentDepth = (parent != InvalidNode) ? bufferDepth[indexOfDepth(parent)] : Real{0};
                bufferDepth[indexOfDepth(node)] = parent != InvalidNode ? parentDepth + Real{1} : Real{0};

                bufferHeight[indexOfHeight(node)] = Real{0};
                bufferNumDesc[indexOfNumDescendants(node)] = Real{0};
                bufferNumLeafDesc[indexOfNumLeafDescendants(node)] = isLeaf ? Real{1} : Real{0};

                if (computeHeight)
                    buffer[attrNames.linearIndex(node, HEIGHT_NODE)] = Real{0};
                if (computeIsLeaf)
                    buffer[attrNames.linearIndex(node, IS_LEAF_NODE)] = isLeaf ? Real{1} : Real{0};
                if (computeIsRoot)
                    buffer[attrNames.linearIndex(node, IS_ROOT_NODE)] = isRoot ? Real{1} : Real{0};
                if (computeNumChildren)
                    buffer[attrNames.linearIndex(node, NUM_CHILDREN_NODE)] = static_cast<Real>(numChildren);
                if (computeNumSiblings)
                    buffer[attrNames.linearIndex(node, NUM_SIBLINGS_NODE)] = isRoot ? Real{0} : static_cast<Real>(tree.getNumChildren(parentNodeId) - 1);
                if (computeLeafRatio)
                    buffer[attrNames.linearIndex(node, LEAF_RATIO_NODE)] = Real{0};
                if (computeBalance) {
                    minChildHeightStorage[static_cast<std::size_t>(node)] = std::numeric_limits<Real>::infinity();
                    buffer[attrNames.linearIndex(node, BALANCE_NODE)] = Real{0};
                }
                if (computeAvgChildHeight)
                    buffer[attrNames.linearIndex(node, AVG_CHILD_HEIGHT_NODE)] = Real{0};
            },
            [&](NodeId parentNodeId, NodeId childNodeId) {
                const NodeId parent = detail::topologySlotOf(tree, parentNodeId);
                const NodeId child = detail::topologySlotOf(tree, childNodeId);

                bufferNumDesc[indexOfNumDescendants(parent)] += bufferNumDesc[indexOfNumDescendants(child)] + Real{1};
                bufferNumLeafDesc[indexOfNumLeafDescendants(parent)] += bufferNumLeafDesc[indexOfNumLeafDescendants(child)];

                Real childHeight = bufferHeight[indexOfHeight(child)];
                Real& parentHeight = bufferHeight[indexOfHeight(parent)];
                parentHeight = std::max(parentHeight, childHeight + Real{1});
                const int numChildren = tree.getNumChildren(parentNodeId);

                if (computeBalance) {
                    Real& minChildHeight = minChildHeightStorage[static_cast<std::size_t>(parent)];
                    minChildHeight = std::min(minChildHeight, childHeight);
                }

                if (computeAvgChildHeight) {
                    Real& sumH = buffer[attrNames.linearIndex(parent, AVG_CHILD_HEIGHT_NODE)];
                    if (numChildren == 1)
                        sumH = childHeight;
                    else
                        sumH += childHeight;
                }
            },
            [&](NodeId idxGlobalId) {
                const NodeId idx = detail::topologySlotOf(tree, idxGlobalId);

                if (computeLeafRatio) {
                    Real desc = bufferNumDesc[indexOfNumDescendants(idx)];
                    Real leafCount = bufferNumLeafDesc[indexOfNumLeafDescendants(idx)];
                    buffer[attrNames.linearIndex(idx, LEAF_RATIO_NODE)] =
                        desc > Real{0} ? ::mmcfilters::attributes::numeric::safeDivide(leafCount, desc + Real{1}) : Real{1};
                }

                if (!tree.isLeaf(idxGlobalId)) {
                    if (computeBalance) {
                        const Real maxChildHeight = bufferHeight[indexOfHeight(idx)] - Real{1};
                        const Real minChildHeight = minChildHeightStorage[static_cast<std::size_t>(idx)];
                        buffer[attrNames.linearIndex(idx, BALANCE_NODE)] = maxChildHeight - minChildHeight;
                    }

                    if (computeAvgChildHeight) {
                        buffer[attrNames.linearIndex(idx, AVG_CHILD_HEIGHT_NODE)] = ::mmcfilters::attributes::numeric::safeDivide(
                            buffer[attrNames.linearIndex(idx, AVG_CHILD_HEIGHT_NODE)], static_cast<Real>(tree.getNumChildren(idxGlobalId)));
                    }
                }
            });
    }

  public:
    /**
     * @brief Materializes topology descriptors for one-pixel unit supports.
     *
     * A unit proper part is represented as a degenerate one-node tree: root and
     * leaf flags are true, depths/heights/child counts are zero, and leaf ratio
     * is one.
     *
     * @param context Operation context or diagnostic label.
     */
    template <std::floating_point Real> static void computeUnitRows(const UnitAttributeComputeContext<Real>& context) {
        const MorphologicalTree& tree = context.tree;
        std::span<const NodeId> unitProperParts = context.unitProperParts;
        std::span<Real> buffer = context.buffer;
        const AttributeNames& attrNames = context.attrNames;
        std::span<const Attribute> requestedAttributes = context.requestedAttributes;

        requireUnitAttributeBufferShape(tree, unitProperParts, buffer, attrNames);

        const bool computeHeight = requestsAttribute(requestedAttributes, HEIGHT_NODE);
        const bool computeDepth = requestsAttribute(requestedAttributes, DEPTH_NODE);
        const bool computeIsLeaf = requestsAttribute(requestedAttributes, IS_LEAF_NODE);
        const bool computeIsRoot = requestsAttribute(requestedAttributes, IS_ROOT_NODE);
        const bool computeNumChildren = requestsAttribute(requestedAttributes, NUM_CHILDREN_NODE);
        const bool computeNumSiblings = requestsAttribute(requestedAttributes, NUM_SIBLINGS_NODE);
        const bool computeNumDescendants = requestsAttribute(requestedAttributes, NUM_DESCENDANTS_NODE);
        const bool computeNumLeafDescendants = requestsAttribute(requestedAttributes, NUM_LEAF_DESCENDANTS_NODE);
        const bool computeLeafRatio = requestsAttribute(requestedAttributes, LEAF_RATIO_NODE);
        const bool computeBalance = requestsAttribute(requestedAttributes, BALANCE_NODE);
        const bool computeAvgChildHeight = requestsAttribute(requestedAttributes, AVG_CHILD_HEIGHT_NODE);

        if (!computeHeight && !computeDepth && !computeIsLeaf && !computeIsRoot && !computeNumChildren && !computeNumSiblings && !computeNumDescendants &&
            !computeNumLeafDescendants && !computeLeafRatio && !computeBalance && !computeAvgChildHeight) {
            return;
        }

        for (NodeId leafIndex = 0; leafIndex < static_cast<NodeId>(unitProperParts.size()); ++leafIndex) {
            if (computeHeight) {
                buffer[attrNames.linearIndex(leafIndex, HEIGHT_NODE)] = Real{0};
            }
            if (computeDepth) {
                buffer[attrNames.linearIndex(leafIndex, DEPTH_NODE)] = Real{0};
            }
            if (computeIsLeaf) {
                buffer[attrNames.linearIndex(leafIndex, IS_LEAF_NODE)] = Real{1};
            }
            if (computeIsRoot) {
                buffer[attrNames.linearIndex(leafIndex, IS_ROOT_NODE)] = Real{1};
            }
            if (computeNumChildren) {
                buffer[attrNames.linearIndex(leafIndex, NUM_CHILDREN_NODE)] = Real{0};
            }
            if (computeNumSiblings) {
                buffer[attrNames.linearIndex(leafIndex, NUM_SIBLINGS_NODE)] = Real{0};
            }
            if (computeNumDescendants) {
                buffer[attrNames.linearIndex(leafIndex, NUM_DESCENDANTS_NODE)] = Real{0};
            }
            if (computeNumLeafDescendants) {
                buffer[attrNames.linearIndex(leafIndex, NUM_LEAF_DESCENDANTS_NODE)] = Real{1};
            }
            if (computeLeafRatio) {
                buffer[attrNames.linearIndex(leafIndex, LEAF_RATIO_NODE)] = Real{1};
            }
            if (computeBalance) {
                buffer[attrNames.linearIndex(leafIndex, BALANCE_NODE)] = Real{0};
            }
            if (computeAvgChildHeight) {
                buffer[attrNames.linearIndex(leafIndex, AVG_CHILD_HEIGHT_NODE)] = Real{0};
            }
        }
    }
};

} // namespace mmcfilters::attributes::computers
