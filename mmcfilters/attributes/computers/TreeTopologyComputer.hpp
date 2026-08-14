#pragma once

#include "AttributeComputerDomain.hpp"
#include "AttributeComputerFamily.hpp"
#include "../detail/AttributeKernelSupport.hpp"
#include "../../trees/detail/CommittedTreeAccess.hpp"
#include "../../trees/detail/TreeTraversalDetail.hpp"
#include "../../trees/MorphologicalTree.hpp"
#include "../../utils/Contract.hpp"

#include <algorithm>
#include <array>
#include <concepts>
#include <limits>
#include <string_view>
#include <vector>

namespace mmcfilters::attributes::computers {

namespace detail {

/** @brief Requested subset of the tree-topology family. */
struct TreeTopologyRequest {
    bool height = false;             ///< Whether subtree height is requested or required.
    bool depth = false;              ///< Whether node depth is requested.
    bool isLeaf = false;             ///< Whether the leaf indicator is requested.
    bool isRoot = false;             ///< Whether the root indicator is requested.
    bool numChildren = false;        ///< Whether direct-child count is requested.
    bool numSiblings = false;        ///< Whether sibling count is requested.
    bool numDescendants = false;     ///< Whether descendant count is requested or required.
    bool numLeafDescendants = false; ///< Whether leaf-descendant count is requested or required.
    bool leafRatio = false;          ///< Whether leaf ratio is requested.
    bool balance = false;            ///< Whether subtree balance is requested.
    bool avgChildHeight = false;     ///< Whether average child height is requested.

    /** @brief Reports whether at least one topology descriptor is requested. @return True when any request flag is set. */
    [[nodiscard]] bool any() const noexcept {
        return height || depth || isLeaf || isRoot || numChildren || numSiblings || numDescendants || numLeafDescendants || leafRatio || balance ||
               avgChildHeight;
    }

    /** @brief Reports whether a height buffer is required. @return True for height-dependent requests. */
    [[nodiscard]] bool needsHeight() const noexcept { return height || balance || avgChildHeight; }
    /** @brief Reports whether descendant counters are required. @return True for descendant-dependent requests. */
    [[nodiscard]] bool needsDescendantCounts() const noexcept { return numDescendants || numLeafDescendants || leafRatio; }

    /**
     * @brief Builds the selection mask from requested scalar attributes.
     * @param requestedAttributes Requested scalar attributes.
     * @return Tree-topology selection mask.
     */
    [[nodiscard]] static TreeTopologyRequest from(std::span<const Attribute> requestedAttributes) {
        return {.height = requestsAttribute(requestedAttributes, SubtreeHeight),
                .depth = requestsAttribute(requestedAttributes, DepthNode),
                .isLeaf = requestsAttribute(requestedAttributes, IsLeafNode),
                .isRoot = requestsAttribute(requestedAttributes, IsRootNode),
                .numChildren = requestsAttribute(requestedAttributes, NumChildrenNode),
                .numSiblings = requestsAttribute(requestedAttributes, NumSiblingsNode),
                .numDescendants = requestsAttribute(requestedAttributes, NumDescendantsNode),
                .numLeafDescendants = requestsAttribute(requestedAttributes, NumLeafDescendantsNode),
                .leafRatio = requestsAttribute(requestedAttributes, LeafRatioNode),
                .balance = requestsAttribute(requestedAttributes, BalanceNode),
                .avgChildHeight = requestsAttribute(requestedAttributes, AvgChildHeightNode)};
    }
};

namespace kernel {

/**
 * @brief Computes requested topology descriptors over an established tree.
 * @param context Established tree, output layout, and output buffer.
 * @param request Topology columns to materialize.
 */
template <std::floating_point Real>
inline void computeTreeTopology(const AttributeComputeContext<Real>& context, const TreeTopologyRequest& request) {
    if (!request.any()) {
        return;
    }

    const int stride = context.attrNames.NUM_ATTRIBUTES;
    const auto offsetOf = [&](Attribute attribute) { return context.attrNames.indexMap.find(attribute)->second; };
    const int heightOffset = request.height ? offsetOf(SubtreeHeight) : 0;
    const int depthOffset = request.depth ? offsetOf(DepthNode) : 0;
    const int isLeafOffset = request.isLeaf ? offsetOf(IsLeafNode) : 0;
    const int isRootOffset = request.isRoot ? offsetOf(IsRootNode) : 0;
    const int numChildrenOffset = request.numChildren ? offsetOf(NumChildrenNode) : 0;
    const int numSiblingsOffset = request.numSiblings ? offsetOf(NumSiblingsNode) : 0;
    const int numDescendantsOffset = request.numDescendants ? offsetOf(NumDescendantsNode) : 0;
    const int numLeafDescendantsOffset = request.numLeafDescendants ? offsetOf(NumLeafDescendantsNode) : 0;
    const int leafRatioOffset = request.leafRatio ? offsetOf(LeafRatioNode) : 0;
    const int balanceOffset = request.balance ? offsetOf(BalanceNode) : 0;
    const int avgChildHeightOffset = request.avgChildHeight ? offsetOf(AvgChildHeightNode) : 0;
    const auto outputIndex = [&](NodeId node, int offset) { return static_cast<std::size_t>(node * stride + offset); };

    const int numNodeSlots = context.tree.numInternalNodeSlots();
    std::vector<Real> heightStorage(request.needsHeight() && !request.height ? static_cast<std::size_t>(numNodeSlots) : 0, Real{0});
    std::vector<Real> numDescendantStorage(request.needsDescendantCounts() && !request.numDescendants ? static_cast<std::size_t>(numNodeSlots) : 0,
                                           Real{0});
    std::vector<Real> numLeafDescendantStorage(
        request.needsDescendantCounts() && !request.numLeafDescendants ? static_cast<std::size_t>(numNodeSlots) : 0, Real{0});
    std::vector<Real> minimumChildHeight(request.balance ? static_cast<std::size_t>(numNodeSlots) : 0, Real{0});

    const auto heightAt = [&](NodeId node) -> Real& {
        return request.height ? context.buffer[outputIndex(node, heightOffset)] : heightStorage[static_cast<std::size_t>(node)];
    };
    const auto numDescendantsAt = [&](NodeId node) -> Real& {
        return request.numDescendants ? context.buffer[outputIndex(node, numDescendantsOffset)]
                                      : numDescendantStorage[static_cast<std::size_t>(node)];
    };
    const auto numLeafDescendantsAt = [&](NodeId node) -> Real& {
        return request.numLeafDescendants ? context.buffer[outputIndex(node, numLeafDescendantsOffset)]
                                          : numLeafDescendantStorage[static_cast<std::size_t>(node)];
    };

    const NodeId root = context.tree.root();
    ::mmcfilters::detail::kernel::traversePostOrder(
        context.tree, root,
        [&](NodeId node) {
            const bool isRoot = node == root;
            const NodeId parent = isRoot ? InvalidNode : ::mmcfilters::detail::CommittedTreeAccess::nodeParent(context.tree, node);
            const int numChildren = ::mmcfilters::detail::CommittedTreeAccess::numChildren(context.tree, node);
            const bool isLeaf = numChildren == 0;

            if (request.depth) {
                context.buffer[outputIndex(node, depthOffset)] =
                    isRoot ? Real{0} : context.buffer[outputIndex(parent, depthOffset)] + Real{1};
            }
            if (request.needsHeight()) {
                heightAt(node) = Real{0};
            }
            if (request.needsDescendantCounts()) {
                numDescendantsAt(node) = Real{0};
                numLeafDescendantsAt(node) = isLeaf ? Real{1} : Real{0};
            }
            if (request.isLeaf)
                context.buffer[outputIndex(node, isLeafOffset)] = isLeaf ? Real{1} : Real{0};
            if (request.isRoot)
                context.buffer[outputIndex(node, isRootOffset)] = isRoot ? Real{1} : Real{0};
            if (request.numChildren)
                context.buffer[outputIndex(node, numChildrenOffset)] = static_cast<Real>(numChildren);
            if (request.numSiblings)
                context.buffer[outputIndex(node, numSiblingsOffset)] =
                    isRoot ? Real{0} : static_cast<Real>(::mmcfilters::detail::CommittedTreeAccess::numChildren(context.tree, parent) - 1);
            if (request.leafRatio)
                context.buffer[outputIndex(node, leafRatioOffset)] = Real{0};
            if (request.balance) {
                minimumChildHeight[static_cast<std::size_t>(node)] = std::numeric_limits<Real>::infinity();
                context.buffer[outputIndex(node, balanceOffset)] = Real{0};
            }
            if (request.avgChildHeight)
                context.buffer[outputIndex(node, avgChildHeightOffset)] = Real{0};
        },
        [&](NodeId parent, NodeId child) {
            if (request.needsDescendantCounts()) {
                numDescendantsAt(parent) += numDescendantsAt(child) + Real{1};
                numLeafDescendantsAt(parent) += numLeafDescendantsAt(child);
            }
            if (request.needsHeight()) {
                const Real childHeight = heightAt(child);
                heightAt(parent) = std::max(heightAt(parent), childHeight + Real{1});
                if (request.balance) {
                    Real& minimumHeight = minimumChildHeight[static_cast<std::size_t>(parent)];
                    minimumHeight = std::min(minimumHeight, childHeight);
                }
                if (request.avgChildHeight) {
                    context.buffer[outputIndex(parent, avgChildHeightOffset)] += childHeight;
                }
            }
        },
        [&](NodeId node) {
            const int numChildren = ::mmcfilters::detail::CommittedTreeAccess::numChildren(context.tree, node);
            if (request.leafRatio) {
                const Real descendantCount = numDescendantsAt(node);
                context.buffer[outputIndex(node, leafRatioOffset)] =
                    descendantCount > Real{0}
                        ? ::mmcfilters::attributes::numeric::safeDivide(numLeafDescendantsAt(node), descendantCount + Real{1})
                        : Real{1};
            }
            if (numChildren != 0) {
                if (request.balance) {
                    context.buffer[outputIndex(node, balanceOffset)] =
                        heightAt(node) - Real{1} - minimumChildHeight[static_cast<std::size_t>(node)];
                }
                if (request.avgChildHeight) {
                    context.buffer[outputIndex(node, avgChildHeightOffset)] = ::mmcfilters::attributes::numeric::safeDivide(
                        context.buffer[outputIndex(node, avgChildHeightOffset)], static_cast<Real>(numChildren));
                }
            }
        });
}

} // namespace kernel

template <std::floating_point Real>
inline void validateTreeTopologyContext(const AttributeComputeContext<Real>& context) {
    requireAttributeBufferShape(context.tree, context.buffer, context.attrNames);
    requireRequestedAttributeColumns(context);
}
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
    inline static constexpr std::array<Attribute, 11> producedAttributes{SubtreeHeight,       DepthNode,        IsLeafNode,         IsRootNode,
                                                                         NumChildrenNode, NumSiblingsNode, NumDescendantsNode, NumLeafDescendantsNode,
                                                                         LeafRatioNode,   BalanceNode,      AvgChildHeightNode};

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
        const detail::TreeTopologyRequest request = detail::TreeTopologyRequest::from(context.requestedAttributes);
        MMCFILTERS_CONTRACT_CHECKED_ONLY(detail::validateTreeTopologyContext(context));
        detail::kernel::computeTreeTopology(context, request);
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
        std::span<const PixelId> unitPixels = context.unitPixels;
        std::span<Real> buffer = context.buffer;
        const AttributeNames& attrNames = context.attrNames;
        std::span<const Attribute> requestedAttributes = context.requestedAttributes;

        requireUnitAttributeBufferShape(tree, unitPixels, buffer, attrNames);

        const bool computeHeight = requestsAttribute(requestedAttributes, SubtreeHeight);
        const bool computeDepth = requestsAttribute(requestedAttributes, DepthNode);
        const bool computeIsLeaf = requestsAttribute(requestedAttributes, IsLeafNode);
        const bool computeIsRoot = requestsAttribute(requestedAttributes, IsRootNode);
        const bool computeNumChildren = requestsAttribute(requestedAttributes, NumChildrenNode);
        const bool computeNumSiblings = requestsAttribute(requestedAttributes, NumSiblingsNode);
        const bool computeNumDescendants = requestsAttribute(requestedAttributes, NumDescendantsNode);
        const bool computeNumLeafDescendants = requestsAttribute(requestedAttributes, NumLeafDescendantsNode);
        const bool computeLeafRatio = requestsAttribute(requestedAttributes, LeafRatioNode);
        const bool computeBalance = requestsAttribute(requestedAttributes, BalanceNode);
        const bool computeAvgChildHeight = requestsAttribute(requestedAttributes, AvgChildHeightNode);

        if (!computeHeight && !computeDepth && !computeIsLeaf && !computeIsRoot && !computeNumChildren && !computeNumSiblings && !computeNumDescendants &&
            !computeNumLeafDescendants && !computeLeafRatio && !computeBalance && !computeAvgChildHeight) {
            return;
        }

        for (NodeId leafIndex = 0; leafIndex < static_cast<NodeId>(unitPixels.size()); ++leafIndex) {
            if (computeHeight) {
                buffer[attrNames.linearIndex(leafIndex, SubtreeHeight)] = Real{0};
            }
            if (computeDepth) {
                buffer[attrNames.linearIndex(leafIndex, DepthNode)] = Real{0};
            }
            if (computeIsLeaf) {
                buffer[attrNames.linearIndex(leafIndex, IsLeafNode)] = Real{1};
            }
            if (computeIsRoot) {
                buffer[attrNames.linearIndex(leafIndex, IsRootNode)] = Real{1};
            }
            if (computeNumChildren) {
                buffer[attrNames.linearIndex(leafIndex, NumChildrenNode)] = Real{0};
            }
            if (computeNumSiblings) {
                buffer[attrNames.linearIndex(leafIndex, NumSiblingsNode)] = Real{0};
            }
            if (computeNumDescendants) {
                buffer[attrNames.linearIndex(leafIndex, NumDescendantsNode)] = Real{0};
            }
            if (computeNumLeafDescendants) {
                buffer[attrNames.linearIndex(leafIndex, NumLeafDescendantsNode)] = Real{1};
            }
            if (computeLeafRatio) {
                buffer[attrNames.linearIndex(leafIndex, LeafRatioNode)] = Real{1};
            }
            if (computeBalance) {
                buffer[attrNames.linearIndex(leafIndex, BalanceNode)] = Real{0};
            }
            if (computeAvgChildHeight) {
                buffer[attrNames.linearIndex(leafIndex, AvgChildHeightNode)] = Real{0};
            }
        }
    }
};

} // namespace mmcfilters::attributes::computers
