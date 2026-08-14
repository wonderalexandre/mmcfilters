#pragma once

#include "AttributeComputerDomain.hpp"
#include "AttributeComputerFamily.hpp"
#include "../detail/AttributeKernelSupport.hpp"
#include "../../contours/ContoursComputedIncrementally.hpp"
#include "../../trees/MorphologicalTree.hpp"
#include "../../trees/TreeAltitudeAlgorithms.hpp"
#include "../../trees/detail/CommittedTreeAccess.hpp"
#include "../../trees/detail/MorphologicalTreeConstructionContextQueries.hpp"
#include "../../trees/detail/TreeTraversalDetail.hpp"
#include "../../utils/Common.hpp"
#include "../../utils/Contract.hpp"

#include "detail/maxdist/EdtDIFT.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

namespace mmcfilters::attributes::computers::detail {

inline void requireMaxDistCapabilities(const MorphologicalTree& tree) {
    if (tree.nodeAltitudeOrder() == NodeAltitudeOrder::Unconstrained) {
        throw std::invalid_argument("MAX_DIST requires a globally monotone altitude order.");
    }
}

template <AltitudeValue T> inline bool isFiniteMaxDistAltitude(T value) noexcept {
    if constexpr (std::is_floating_point_v<T>) {
        return std::isfinite(value);
    }
    return true;
}

template <AltitudeValue T> inline void validateFiniteMaxDistAltitude(std::span<const T> altitude) {
    for (std::size_t index = 0; index < altitude.size(); ++index) {
        if (!isFiniteMaxDistAltitude(altitude[index])) {
            throw std::invalid_argument("MAX_DIST requires finite altitude values; node " + std::to_string(index) + " has a non-finite altitude.");
        }
    }
}

template <AltitudeValue T> inline void validateMaxDistInput(const MorphologicalTree& tree, std::span<const T> altitude) {
    requireMaxDistCapabilities(tree);
    if (::mmcfilters::detail::constructionAdjacency(tree) == nullptr) {
        throw std::invalid_argument("MAX_DIST requires a shared or saturated construction adjacency.");
    }
    TreeAltitudeAlgorithms::validateNodeAltitudeBufferShape(tree, altitude);
    validateFiniteMaxDistAltitude(altitude);
}

namespace kernel {

/**
 * @brief Internal MAX_DIST altitude-sweep kernel.
 *
 * @details
 * MAX_DIST is evaluated by sweeping tree nodes in typed altitude order and
 * maintaining an incremental squared Euclidean distance transform over the
 * support accumulated at the current level. Local contour additions and
 * removals are shared with `ContoursComputedIncrementally`, avoiding a second
 * dense level-image contour pass.
 *
 * The kernel returns a dense vector indexed by internal node id. The public
 * `MaxDistComputer` is responsible for projecting that vector into the
 * caller-owned attribute buffer layout.
 */
class MaxDistAttribute {
  public:
    /**
     * @brief Runs the complete MAX_DIST computation in dense node-id space.
     *
     * Floating-point altitudes must be finite because the level ordering and
     * contour tests require a total finite order.
     *
     * @param tree Tree topology.
     * @param altitude Altitude data indexed by node identifier.
     * @return Values produced by the operation.
     */
    template <std::floating_point Real, AltitudeValue T>
    [[nodiscard]] static std::vector<Real> compute(const MorphologicalTree& tree, std::span<const T> altitude) {
        std::vector<Real> maxDist(static_cast<std::size_t>(tree.numInternalNodeSlots()), Real{0});
        const GridDomain2D& domain = ::mmcfilters::detail::CommittedTreeAccess::gridDomain2D(tree);
        maxdist::EdtDIFT edtDIFT(domain.rows, domain.columns);

        const ContourDeltaStore contourDeltas = ::mmcfilters::detail::kernel::extractLocalContourDeltas(tree);
        std::vector<std::vector<PixelId>> contours(static_cast<std::size_t>(tree.numInternalNodeSlots()));
        const std::size_t totalPixels = static_cast<std::size_t>(domain.rows * domain.columns);
        std::vector<std::uint8_t> removalMark(totalPixels, 0);
        std::vector<std::uint8_t> contourAdditionMark(totalPixels, 0);

        std::vector<NodeId> sortedNodes = sortedNodesByAltitude(tree, altitude);
        std::size_t groupBegin = 0;
        while (groupBegin < sortedNodes.size()) {
            std::size_t groupEnd = groupBegin + 1;
            while (groupEnd < sortedNodes.size() && sameAltitude(altitude, sortedNodes[groupBegin], sortedNodes[groupEnd])) {
                ++groupEnd;
            }

            processLevel(tree, std::span<const NodeId>(sortedNodes.data() + static_cast<std::ptrdiff_t>(groupBegin), groupEnd - groupBegin), contourDeltas,
                         edtDIFT, contours, removalMark, contourAdditionMark, maxDist);

            groupBegin = groupEnd;
        }

        return maxDist;
    }

  private:
    /** @brief Defines the `ContourDeltaStore` alias used by the component. */
    using ContourDeltaStore = ::mmcfilters::ContoursComputedIncrementally::LocalContourDeltas;

    /**
     * @brief Marks pixels.
     *
     * @param pixels Pixel identifiers.
     * @param marks Mark buffer updated by the operation.
     */
    static void markPixels(std::span<const PixelId> pixels, std::vector<std::uint8_t>& marks) {
        for (PixelId pixelId : pixels) {
            marks[static_cast<std::size_t>(pixelId)] = 1;
        }
    }

    /**
     * @brief Clears pixel marks.
     *
     * @param pixels Pixel identifiers.
     * @param marks Mark buffer updated by the operation.
     */
    static void clearPixelMarks(std::span<const PixelId> pixels, std::vector<std::uint8_t>& marks) {
        for (PixelId pixelId : pixels) {
            marks[static_cast<std::size_t>(pixelId)] = 0;
        }
    }

    /**
     * @brief Groups nodes by sorted altitude, independently of the altitude
     * value domain.
     *
     * The vector is ordered by the level sweep required by MAX_DIST: descending
     * for max-trees and ascending for min-trees. Stable sorting preserves
     * post-order among unrelated nodes that share an altitude; the declared
     * ordered contracts already forbid equal altitudes on a parent-child arc.
     *
     * @param tree Tree topology.
     * @param altitude Altitude data indexed by node identifier.
     * @return Values produced by the operation.
     */
    template <AltitudeValue T> [[nodiscard]] static std::vector<NodeId> sortedNodesByAltitude(const MorphologicalTree& tree, std::span<const T> altitude) {
        std::vector<NodeId> nodes;
        nodes.reserve(static_cast<std::size_t>(tree.numNodes()));
        ::mmcfilters::detail::kernel::traversePostOrder(
            tree, tree.root(), [](NodeId) {}, [](NodeId, NodeId) {}, [&](NodeId nodeId) { nodes.push_back(nodeId); });

        std::stable_sort(nodes.begin(), nodes.end(), [&](NodeId lhs, NodeId rhs) {
            const T lhsAltitude = altitude[static_cast<std::size_t>(lhs)];
            const T rhsAltitude = altitude[static_cast<std::size_t>(rhs)];
            return tree.nodeAltitudeOrder() == NodeAltitudeOrder::Increasing ? lhsAltitude > rhsAltitude : lhsAltitude < rhsAltitude;
        });

        return nodes;
    }

    /**
     * @brief Tests whether two nodes belong to the same altitude group.
     *
     * @param altitude Altitude data indexed by node identifier.
     * @param lhs Left-hand operand.
     * @param rhs Right-hand operand.
     * @return True if two nodes belong to the same altitude group; otherwise false.
     */
    template <AltitudeValue T> static bool sameAltitude(std::span<const T> altitude, NodeId lhs, NodeId rhs) {
        return altitude[static_cast<std::size_t>(lhs)] == altitude[static_cast<std::size_t>(rhs)];
    }

    /**
     * @brief Processes one group of nodes that share the same altitude.
     *
     * Children contours are inherited into the parent contour unless the shared
     * contour-delta store says the pixel must be removed at this node. Proper
     * parts are then added as contour seeds when they are local additions, or
     * opened as interior pixels otherwise. Only after all nodes in the level
     * group have been materialized does EdtDIFT propagate labels; this preserves
     * the simultaneous per-level sweep.
     *
     * @param tree Tree topology.
     * @param nodes Node identifiers processed by the operation.
     * @param contourDeltas Compact local contour additions and removals.
     * @param edtDIFT Distance-transform state updated by the sweep.
     * @param contours Contour data.
     * @param removalMark Generation marks for pixels scheduled for removal.
     * @param contourAdditionMark Generation marks for contour additions.
     * @param maxDist Destination MAX_DIST values indexed by node.
     */
    template <std::floating_point Real>
    static void processLevel(const MorphologicalTree& tree, std::span<const NodeId> nodes, const ContourDeltaStore& contourDeltas, maxdist::EdtDIFT& edtDIFT,
                             std::vector<std::vector<PixelId>>& contours, std::vector<std::uint8_t>& removalMark,
                             std::vector<std::uint8_t>& contourAdditionMark,
                             std::vector<Real>& maxDist) {
        if (nodes.empty()) {
            return;
        }

        std::vector<PixelId> toRemove;
        toRemove.reserve(64);
        for (NodeId nodeId : nodes) {
            std::vector<PixelId>& nodeContour = contours[static_cast<std::size_t>(nodeId)];
            nodeContour.clear();

            const auto removals = contourDeltas.removals(nodeId);
            toRemove.clear();
            toRemove.insert(toRemove.end(), removals.begin(), removals.end());
            markPixels(removals, removalMark);

            const auto additions = contourDeltas.additions(nodeId);
            std::size_t reserveSize = additions.size();
            for (NodeId childNodeId : ::mmcfilters::detail::CommittedTreeAccess::children(tree, nodeId)) {
                reserveSize += contours[static_cast<std::size_t>(childNodeId)].size();
            }
            nodeContour.reserve(reserveSize);

            for (NodeId childNodeId : ::mmcfilters::detail::CommittedTreeAccess::children(tree, nodeId)) {
                std::vector<PixelId>& childContour = contours[static_cast<std::size_t>(childNodeId)];
                for (PixelId pixelId : childContour) {
                    if (!removalMark[static_cast<std::size_t>(pixelId)]) {
                        nodeContour.push_back(pixelId);
                    }
                }
                std::vector<PixelId>().swap(childContour);
            }
            clearPixelMarks(removals, removalMark);

            if (!toRemove.empty()) {
                edtDIFT.treeRemoval(toRemove);
            }

            markPixels(additions, contourAdditionMark);
            for (PixelId pixelId : ::mmcfilters::detail::CommittedTreeAccess::properParts(tree, nodeId)) {
                edtDIFT.addPixelToBinaryImage(pixelId);

                if (contourAdditionMark[static_cast<std::size_t>(pixelId)]) {
                    nodeContour.push_back(pixelId);
                    edtDIFT.seed(pixelId);
                } else {
                    edtDIFT.open(pixelId);
                    edtDIFT.insertNeighborsPQueue(pixelId);
                }
            }
            clearPixelMarks(additions, contourAdditionMark);
        }

        // All nodes at the same altitude must enter the binary image before the
        // distance transform propagates. This preserves the mathematical
        // per-level schedule without requiring fixed 0..255 buckets.
        edtDIFT.run();

        for (NodeId nodeId : nodes) {
            maxDist[static_cast<std::size_t>(nodeId)] = static_cast<Real>(edtDIFT.maxBedt(contours[static_cast<std::size_t>(nodeId)]));
        }
    }
};

/**
 * @brief Writes MAX_DIST after output, topology, altitude, and request domains were established.
 * @param context Established tree, altitude span, MAX_DIST column, and output buffer.
 */
template <std::floating_point Real, AltitudeValue T> inline void computeMaxDist(const AltitudeAttributeComputeContext<Real, T>& context) {
    const std::vector<Real> maxDist = MaxDistAttribute::compute<Real>(context.tree, context.altitude);
    const int stride = context.attrNames.NUM_ATTRIBUTES;
    const int offset = context.attrNames.indexMap.find(MaxDist)->second;
    for (NodeId nodeId = 0; nodeId < context.tree.numInternalNodeSlots(); ++nodeId) {
        if (::mmcfilters::detail::CommittedTreeAccess::isAlive(context.tree, nodeId)) {
            context.buffer[static_cast<std::size_t>(nodeId * stride + offset)] = maxDist[static_cast<std::size_t>(nodeId)];
        }
    }
}

} // namespace kernel

} // namespace mmcfilters::attributes::computers::detail

namespace mmcfilters::attributes::computers {

/**
 * @brief Stateless MAX_DIST scalar computer.
 *
 * @details
 * The computer exposes the standard attribute-computer protocol for MAX_DIST:
 * it receives a typed altitude-aware compute context, runs the internal
 * distance-transform sweep in dense node-id space, and writes the requested
 * scalar column into the caller-owned buffer.
 */
class MaxDistComputer {
  public:
    /// Family name used in dependency-plan diagnostics.
    static constexpr std::string_view familyName = "max-dist";

    /// Stable family id used by the scheduler.
    static constexpr AttributeComputerFamily family = AttributeComputerFamily::MaxDist;

    /// Execution domain required by the computer.
    static constexpr AttributeComputerDomain domain = AttributeComputerDomain::Altitude;

    /**
     * @brief Canonical list of scalar descriptors materialized by this computer.
     */
    inline static constexpr std::array<Attribute, 1> producedAttributes{MaxDist};

    /**
     * @brief Rejects tree kinds for which MAX_DIST is not defined.
     *
     * @param tree Tree topology.
     */
    static void requireSupportedTreeKind(const MorphologicalTree& tree) { detail::requireMaxDistCapabilities(tree); }

    /**
     * @brief Computes MAX_DIST and writes it into the requested output layout.
     *
     * @details
     * Requires max-tree or min-tree topology with adjacency metadata and a
     * dense typed altitude span. Floating-point altitude values must be finite.
     *
     * @param context Operation context or diagnostic label.
     */
    template <std::floating_point Real, AltitudeValue T> static void compute(const AltitudeAttributeComputeContext<Real, T>& context) {
        MMCFILTERS_CONTRACT_CHECKED_ONLY(requireAttributeBufferShape(context.tree, context.buffer, context.attrNames);
                                         detail::validateMaxDistInput(context.tree, context.altitude));
        if (!requestsAttribute(context.requestedAttributes, MaxDist)) {
            return;
        }

        detail::kernel::computeMaxDist(context);
    }

    /**
     * @brief Materializes MAX_DIST for one-pixel unit supports.
     *
     * A one-pixel support has no interior distance from its contour, so its
     * unit MAX_DIST value is zero.
     *
     * @param context Operation context or diagnostic label.
     */
    template <std::floating_point Real, AltitudeValue T> static void computeUnitRows(const AltitudeUnitAttributeComputeContext<Real, T>& context) {
        requireUnitAttributeBufferShape(context.tree, context.unitPixels, context.buffer, context.attrNames);
        TreeAltitudeAlgorithms::validateNodeAltitudeBufferShape(context.tree, context.altitude);
        if (!requestsAttribute(context.requestedAttributes, MaxDist)) {
            return;
        }

        for (NodeId leafIndex = 0; leafIndex < static_cast<NodeId>(context.unitPixels.size()); ++leafIndex) {
            context.buffer[context.attrNames.linearIndex(leafIndex, MaxDist)] = Real{0};
        }
    }
};

} // namespace mmcfilters::attributes::computers
