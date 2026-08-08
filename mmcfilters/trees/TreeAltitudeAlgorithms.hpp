#pragma once

#include "../utils/Altitude.hpp"
#include "../utils/Image.hpp"
#include "MorphologicalTree.hpp"
#include "detail/HigraExportLayoutDetail.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace mmcfilters {

/**
 * @brief Pure operations over a topology and an explicit altitude buffer.
 *
 * `TreeAltitudeAlgorithms` intentionally owns no state. It is the neutral home for
 * algorithms that need only `MorphologicalTree + std::span<const T>`, leaving
 * `WeightedMorphologicalTree<std::uint8_t>` as an owner/adaptor for canonical mutable state.
 */
class TreeAltitudeAlgorithms {
  public:
    /**
     * @brief Validates that an altitude buffer covers the dense internal-node domain.
     *
     * @param tree Tree topology used by the operation.
     * @param altitude Altitude data indexed by node identifier.
     */
    template <AltitudeValue T> static void validateAltitudeBufferShape(const MorphologicalTree& tree, std::span<const T> altitude) {
        if (altitude.size() != static_cast<std::size_t>(tree.getNumInternalNodeSlots())) {
            throw std::runtime_error("Altitude buffer size must match the dense internal-node domain.");
        }
    }

    /**
     * @brief Rejects non-finite floating-point altitudes while compiling to a no-op for integral types.
     *
     * @param altitude Altitude data indexed by node identifier.
     * @param index Zero-based index used by the operation.
     * @param context Operation context or diagnostic label.
     */
    template <AltitudeValue T> static void validateFiniteAltitudeValue(T altitude, std::size_t index, const char* context) {
        if constexpr (std::is_floating_point_v<T>) {
            const long double level = static_cast<long double>(altitude);
            if (!std::isfinite(level)) {
                std::ostringstream oss;
                oss << context << " requires finite floating-point altitudes; value at index " << index << " is " << level << ".";
                throw std::invalid_argument(oss.str());
            }
        }
    }

    /**
     * @brief Rejects non-finite floating-point altitudes in a contiguous input range.
     *
     * @param altitude Altitude data indexed by node identifier.
     * @param context Operation context or diagnostic label.
     */
    template <AltitudeValue T> static void validateFiniteAltitudeValues(std::span<const T> altitude, const char* context) {
        if constexpr (std::is_floating_point_v<T>) {
            for (std::size_t index = 0; index < altitude.size(); ++index) {
                validateFiniteAltitudeValue(altitude[index], index, context);
            }
        }
    }

    /**
     * @brief Rejects non-finite floating-point pixels before using an image as altitude source.
     *
     * @param image Image used by the operation.
     * @param context Operation context or diagnostic label.
     */
    template <AltitudeValue T> static void validateFiniteImageAltitudes(const ImagePtr<T>& image, const char* context) {
        if constexpr (std::is_floating_point_v<T>) {
            if (!image) {
                throw std::invalid_argument("Image altitude validation requires a non-null image.");
            }
            validateFiniteAltitudeValues(std::span<const T>(image->rawData(), static_cast<std::size_t>(image->getSize())), context);
        }
    }

    /**
     * @brief Reads one node altitude from an explicit altitude buffer.
     *
     * @param altitude Altitude data indexed by node identifier.
     * @param nodeId Identifier of the node used by the operation.
     * @return The requested node altitude from an explicit altitude buffer.
     */
    template <AltitudeValue T> [[nodiscard]] static T getAltitude(std::span<const T> altitude, NodeId nodeId) {
        if (nodeId < 0 || static_cast<std::size_t>(nodeId) >= altitude.size()) {
            throw std::invalid_argument("Altitude access requires a valid internal NodeId.");
        }
        return altitude[static_cast<std::size_t>(nodeId)];
    }

    /**
     * @brief Computes the altitude difference between one node and its parent.
     *
     * @param tree Tree topology used by the operation.
     * @param altitude Altitude data indexed by node identifier.
     * @param nodeId Identifier of the node used by the operation.
     * @return The computed altitude difference between one node and its parent.
     */
    template <AltitudeValue T> [[nodiscard]] static AltitudeDiff<T> getNodeResidue(const MorphologicalTree& tree, std::span<const T> altitude, NodeId nodeId) {
        validateAltitudeBufferShape(tree, altitude);
        if (!tree.isAlive(nodeId)) {
            throw std::invalid_argument("Node residue requires a live internal NodeId.");
        }
        const NodeId parentNodeId = tree.getNodeParent(nodeId);
        if (parentNodeId == InvalidNode || parentNodeId == nodeId) {
            return static_cast<AltitudeDiff<T>>(getAltitude(altitude, nodeId));
        }
        return static_cast<AltitudeDiff<T>>(getAltitude(altitude, nodeId)) - static_cast<AltitudeDiff<T>>(getAltitude(altitude, parentNodeId));
    }

    /**
     * @brief Converts one altitude value to `uint8_t`, rejecting values outside the output domain.
     *
     * @param altitude Altitude data indexed by node identifier.
     * @param nodeId Identifier of the node used by the operation.
     * @param context Operation context or diagnostic label.
     * @return The converted one altitude value to uint8_t, rejecting values outside the output domain.
     */
    template <AltitudeValue T> [[nodiscard]] static std::uint8_t requireUInt8AltitudeValue(T altitude, NodeId nodeId, const char* context) {
        const long double level = static_cast<long double>(altitude);
        if constexpr (std::is_floating_point_v<T>) {
            if (!std::isfinite(level)) {
                std::ostringstream oss;
                oss << context << " requires finite node altitudes in the uint8 domain [0, 255]; node " << nodeId << " has altitude " << level << ".";
                throw std::invalid_argument(oss.str());
            }
        }
        if (level < 0.0L || level > 255.0L) {
            std::ostringstream oss;
            oss << context << " requires node altitudes in the uint8 domain [0, 255]; node " << nodeId << " has altitude " << level << ".";
            throw std::invalid_argument(oss.str());
        }
        return static_cast<std::uint8_t>(altitude);
    }

    /**
     * @brief Validates all live node altitudes before materialising an `ImageUInt8`.
     *
     * @param tree Tree topology used by the operation.
     * @param altitude Altitude data indexed by node identifier.
     * @param context Operation context or diagnostic label.
     */
    template <AltitudeValue T> static void validateUInt8AltitudeDomain(const MorphologicalTree& tree, std::span<const T> altitude, const char* context) {
        validateAltitudeBufferShape(tree, altitude);
        for (NodeId nodeId : tree.getAliveNodeIds()) {
            (void)requireUInt8AltitudeValue(getAltitude(altitude, nodeId), nodeId, context);
        }
    }

    /**
     * @brief Reconstructs a typed image from topology ownership and explicit node altitudes.
     *
     * The reconstructed pixel type is the altitude type itself. This method does
     * not clamp, convert, or validate altitude values; it only checks the buffer
     * shape before indexing it through the topology.
     *
     * @param tree Tree topology used by the operation.
     * @param altitude Altitude data indexed by node identifier.
     * @param context Operation context or diagnostic label.
     * @return The reconstructed typed image from topology ownership and explicit node altitudes.
     */
    template <AltitudeValue T>
    [[nodiscard]] static ImagePtr<T> reconstructImage(const MorphologicalTree& tree, std::span<const T> altitude,
                                                      const char* context = "TreeAltitudeAlgorithms::reconstructImage") {
        (void)context;
        tree.requireNotEditing(context);
        validateAltitudeBufferShape(tree, altitude);
        ImagePtr<T> image = Image<T>::create(tree.getNumRowsOfGridDomain2D(), tree.getNumColsOfGridDomain2D());
        auto imgBuffer = image->rawData();
        for (int properPartId = 0; properPartId < tree.getNumTotalProperParts(); ++properPartId) {
            const NodeId nodeId = tree.getProperPartOwner(properPartId);
            imgBuffer[static_cast<std::size_t>(properPartId)] = getAltitude(altitude, nodeId);
        }
        return image;
    }

    /**
     * @brief Exports a live rooted topology and explicit altitudes to a compact parent/altitude representation.
     *
     * @details
     * This method owns only the structural weighted export. It reuses
     * `detail::computeExportedHigraLayout()` so the parent/altitude export has
     * the same compact id convention used by attribute-buffer projection in
     * `AttributeComputation`.
     *
     * @param tree Tree topology used by the operation.
     * @param altitude Altitude data indexed by node identifier.
     * @return The exported live rooted topology and explicit altitudes to a compact parent/altitude representation.
     */
    template <AltitudeValue T>
    [[nodiscard]] static std::pair<std::vector<NodeId>, std::vector<T>> exportHigraHierarchy(const MorphologicalTree& tree, std::span<const T> altitude) {
        tree.requireNotEditing("TreeAltitudeAlgorithms::exportHigraHierarchy");
        const detail::ExportedHigraLayout layout = detail::computeExportedHigraLayout(tree, altitude);
        const NodeId numLeaves = layout.numLeaves;
        const NodeId numVertices = layout.numVertices;

        std::vector<NodeId> parent(static_cast<std::size_t>(numVertices), InvalidNode);
        std::vector<T> exportedAltitude(static_cast<std::size_t>(numVertices), T{});

        for (NodeId oldNodeId : layout.sortedNodes) {
            const NodeId newNodeId = layout.nodeToHigra[static_cast<std::size_t>(oldNodeId)];
            exportedAltitude[static_cast<std::size_t>(newNodeId)] = getAltitude(altitude, oldNodeId);
        }

        for (NodeId leafIndex = 0; leafIndex < numLeaves; ++leafIndex) {
            const NodeId properPart = layout.properParts[static_cast<std::size_t>(leafIndex)];
            const NodeId ownerNodeId = tree.getProperPartOwner(properPart);
            if (ownerNodeId == InvalidNode || !tree.isAlive(ownerNodeId)) {
                throw std::runtime_error("Each proper part must belong to one alive node when exporting a compact Higra hierarchy.");
            }
            parent[static_cast<std::size_t>(leafIndex)] = layout.nodeToHigra[static_cast<std::size_t>(ownerNodeId)];
            exportedAltitude[static_cast<std::size_t>(leafIndex)] = getAltitude(altitude, ownerNodeId);
        }

        for (NodeId oldNodeId : layout.sortedNodes) {
            const NodeId newNodeId = layout.nodeToHigra[static_cast<std::size_t>(oldNodeId)];
            const NodeId oldParentNodeId = tree.getNodeParent(oldNodeId);
            parent[static_cast<std::size_t>(newNodeId)] =
                oldParentNodeId == oldNodeId ? newNodeId : layout.nodeToHigra[static_cast<std::size_t>(oldParentNodeId)];
        }

        return {std::move(parent), std::move(exportedAltitude)};
    }

    /**
     * @brief Validates the hierarchy's declared global altitude order.
     *
     * @param tree Tree topology used by the operation.
     * @param altitude Altitude data indexed by node identifier.
     */
    template <AltitudeValue T> static void validateMonotoneAltitude(const MorphologicalTree& tree, std::span<const T> altitude) {
        validateAltitudeBufferShape(tree, altitude);
        const AltitudeOrder altitudeOrder = tree.getAltitudeOrder();
        bool increasingFromRoot = false;
        switch (altitudeOrder) {
        case AltitudeOrder::INCREASING_FROM_ROOT:
            increasingFromRoot = true;
            break;
        case AltitudeOrder::DECREASING_FROM_ROOT:
            increasingFromRoot = false;
            break;
        case AltitudeOrder::UNCONSTRAINED:
            return;
        }

        for (NodeId nodeId : tree.getAliveNodeIds()) {
            if (tree.isRoot(nodeId)) {
                continue;
            }

            const NodeId parentNodeId = tree.getNodeParent(nodeId);
            if (parentNodeId == InvalidNode || !tree.isAlive(parentNodeId)) {
                throw std::runtime_error("Monotonic validation requires every alive non-root node to have an alive parent.");
            }

            if (increasingFromRoot) {
                if (getAltitude(altitude, parentNodeId) >= getAltitude(altitude, nodeId)) {
                    throw std::runtime_error("Hierarchy altitude buffer must be strictly increasing from parent to child.");
                }
            } else if (getAltitude(altitude, parentNodeId) <= getAltitude(altitude, nodeId)) {
                throw std::runtime_error("Hierarchy altitude buffer must be strictly decreasing from parent to child.");
            }
        }
    }
};

} // namespace mmcfilters
