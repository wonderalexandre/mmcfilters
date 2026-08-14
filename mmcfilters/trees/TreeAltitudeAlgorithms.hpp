#pragma once

#include "../utils/Altitude.hpp"
#include "../utils/Contract.hpp"
#include "../utils/Image.hpp"
#include "MorphologicalTree.hpp"
#include "detail/HigraExportLayoutDetail.hpp"

#include <algorithm>
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

namespace detail::tree_altitude {

template <class Contribution>
    requires(std::is_arithmetic_v<Contribution> && !std::is_same_v<std::remove_cv_t<Contribution>, bool>)
[[nodiscard]] std::vector<Contribution> reconstructNodeContributionValues(const MorphologicalTree& tree,
                                                                           std::span<const Contribution> nodeContributions,
                                                                           const char* context) {
    tree.requireNotEditing(context);
    MMCFILTERS_CONTRACT_REQUIRE(
        nodeContributions.size() == static_cast<std::size_t>(tree.numInternalNodeSlots()),
        throw std::invalid_argument(std::string(context) + " nodeContributions size must match the internal node slot count."));
    if constexpr (std::is_floating_point_v<Contribution> && contract::validationsEnabled) {
        for (Contribution contribution : nodeContributions) {
            if (!std::isfinite(contribution)) {
                throw std::invalid_argument(std::string(context) + " requires finite nodeContributions.");
            }
        }
    }

    std::vector<Contribution> accumulated(static_cast<std::size_t>(tree.numInternalNodeSlots()), Contribution{});
    const NodeId root = tree.root();
    accumulated[static_cast<std::size_t>(root)] = nodeContributions[static_cast<std::size_t>(root)];

    std::vector<NodeId> pending{root};
    while (!pending.empty()) {
        const NodeId nodeId = pending.back();
        pending.pop_back();
        for (NodeId childId : tree.children(nodeId)) {
            accumulated[static_cast<std::size_t>(childId)] =
                accumulated[static_cast<std::size_t>(nodeId)] + nodeContributions[static_cast<std::size_t>(childId)];
            pending.push_back(childId);
        }
    }

    std::vector<Contribution> pixels(static_cast<std::size_t>(tree.numPixels()), Contribution{});
    for (PixelId pixel = 0; pixel < tree.numPixels(); ++pixel) {
        pixels[static_cast<std::size_t>(pixel)] = accumulated[static_cast<std::size_t>(tree.smallestNode(pixel))];
    }
    return pixels;
}

} // namespace detail::tree_altitude

/**
 * @brief Pure operations over a topology and an explicit altitude buffer.
 *
 * `TreeAltitudeAlgorithms` intentionally owns no state. It is the neutral home for
 * algorithms that need only `MorphologicalTree + std::span<const T>`, leaving
 * `ValuedMorphologicalTree<std::uint8_t>` as an owner/adaptor for canonical mutable state.
 */
class TreeAltitudeAlgorithms {
  public:
    /**
     * @brief Validates that an altitude buffer covers the dense internal-node domain.
     *
     * @param tree Tree topology.
     * @param altitude Altitude data indexed by node identifier.
     */
    template <AltitudeValue T> static void validateNodeAltitudeBufferShape(const MorphologicalTree& tree, std::span<const T> altitude) {
        MMCFILTERS_CONTRACT_REQUIRE(altitude.size() == static_cast<std::size_t>(tree.numInternalNodeSlots()),
                                    throw std::runtime_error("Altitude buffer size must match the dense internal-node domain."));
    }

    /**
     * @brief Rejects non-finite floating-point altitudes while compiling to a no-op for integral types.
     *
     * @param altitude Altitude data indexed by node identifier.
     * @param index Zero-based index.
     * @param context Operation context or diagnostic label.
     */
    template <AltitudeValue T> static void validateFiniteAltitudeValue(T altitude, std::size_t index, const char* context) {
        if constexpr (std::is_floating_point_v<T>) {
            const long double level = static_cast<long double>(altitude);
            MMCFILTERS_CONTRACT_REQUIRE(std::isfinite(level), {
                std::ostringstream oss;
                oss << context << " requires finite floating-point altitudes; value at index " << index << " is " << level << ".";
                throw std::invalid_argument(oss.str());
            });
        }
    }

    /**
     * @brief Rejects non-finite floating-point altitudes in a contiguous input range.
     *
     * @param altitude Altitude data indexed by node identifier.
     * @param context Operation context or diagnostic label.
     */
    template <AltitudeValue T> static void validateFiniteAltitudeValues(std::span<const T> altitude, const char* context) {
        if constexpr (std::is_floating_point_v<T> && contract::validationsEnabled) {
            for (std::size_t index = 0; index < altitude.size(); ++index) {
                validateFiniteAltitudeValue(altitude[index], index, context);
            }
        }
    }

    /**
     * @brief Rejects non-finite floating-point pixels before using an image as altitude source.
     *
     * @param image Image.
     * @param context Operation context or diagnostic label.
     */
    template <AltitudeValue T> static void validateFiniteImageAltitudes(const ImagePtr<T>& image, const char* context) {
        if constexpr (std::is_floating_point_v<T> && contract::validationsEnabled) {
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
     * @param nodeId Dense internal node identifier.
     * @return The requested node altitude from an explicit altitude buffer.
     */
    template <AltitudeValue T> [[nodiscard]] static T nodeAltitude(std::span<const T> altitude, NodeId nodeId) {
        MMCFILTERS_CONTRACT_REQUIRE(nodeId >= 0 && static_cast<std::size_t>(nodeId) < altitude.size(),
                                    throw std::invalid_argument("Altitude access requires a valid internal NodeId."));
        return altitude[static_cast<std::size_t>(nodeId)];
    }

    /**
     * @brief Computes the altitude difference between one node and its parent.
     *
     * The reconstruction baseline is fixed at zero throughout the project.
     * Consequently, a root-like node whose parent is invalid or itself has
     * residue equal to its own altitude.
     *
     * @param tree Tree topology.
     * @param altitude Altitude data indexed by node identifier.
     * @param nodeId Dense internal node identifier.
     * @return The computed altitude difference between one node and its parent.
     */
    template <AltitudeValue T> [[nodiscard]] static AltitudeDifference<T> nodeResidue(const MorphologicalTree& tree, std::span<const T> altitude, NodeId nodeId) {
        validateNodeAltitudeBufferShape(tree, altitude);
        MMCFILTERS_CONTRACT_REQUIRE(tree.isAlive(nodeId), throw std::invalid_argument("Node residue requires a live internal NodeId."));
        const NodeId parentNodeId = tree.parent(nodeId);
        if (parentNodeId == InvalidNode || parentNodeId == nodeId) {
            return static_cast<AltitudeDifference<T>>(nodeAltitude(altitude, nodeId));
        }
        return static_cast<AltitudeDifference<T>>(nodeAltitude(altitude, nodeId)) - static_cast<AltitudeDifference<T>>(nodeAltitude(altitude, parentNodeId));
    }

    /**
     * @brief Converts one altitude value to `uint8_t`, rejecting values outside the output domain.
     *
     * @param altitude Altitude data indexed by node identifier.
     * @param nodeId Dense internal node identifier.
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
     * @param tree Tree topology.
     * @param altitude Altitude data indexed by node identifier.
     * @param context Operation context or diagnostic label.
     */
    template <AltitudeValue T> static void validateUInt8AltitudeDomain(const MorphologicalTree& tree, std::span<const T> altitude, const char* context) {
        validateNodeAltitudeBufferShape(tree, altitude);
        for (NodeId nodeId : tree.aliveNodeIds()) {
            (void)requireUInt8AltitudeValue(nodeAltitude(altitude, nodeId), nodeId, context);
        }
    }

    /**
     * @brief Reconstructs a typed image from topology storage and explicit node altitudes.
     *
     * The reconstructed pixel type is the altitude type itself. This method does
     * not clamp, convert, or validate altitude values; it only checks the buffer
     * shape before indexing it through the topology.
     *
     * @param tree Tree topology.
     * @param altitude Altitude data indexed by node identifier.
     * @param context Operation context or diagnostic label.
     * @return The reconstructed typed image from topology storage and explicit node altitudes.
     */
    template <AltitudeValue T>
    [[nodiscard]] static ImagePtr<T> reconstructFromNodeAltitudes(const MorphologicalTree& tree, std::span<const T> altitude,
                                                      const char* context = "TreeAltitudeAlgorithms::reconstructFromNodeAltitudes") {
        (void)context;
        tree.requireNotEditing(context);
        validateNodeAltitudeBufferShape(tree, altitude);
        ImagePtr<T> image = Image<T>::create(tree.numRows(), tree.numColumns());
        auto imgBuffer = image->rawData();
        for (PixelId pixel = 0; pixel < tree.numPixels(); ++pixel) {
            const NodeId nodeId = tree.smallestNode(pixel);
            imgBuffer[static_cast<std::size_t>(pixel)] = nodeAltitude(altitude, nodeId);
        }
        return image;
    }

    /**
     * @brief Reconstructs an image by summing node contributions on every root-to-node branch.
     *
     * The reconstruction baseline is fixed at zero. Therefore the accumulated
     * value at the root is its own contribution, and each non-root accumulated
     * value is the parent accumulation plus the node contribution.
     *
     * @param tree Tree topology indexing the contribution buffer.
     * @param nodeContributions Dense contribution buffer indexed by internal node id.
     * @param context Operation name used in diagnostics.
     * @return Image whose pixel value is the sum of contributions of all nodes containing that pixel.
     */
    template <class Contribution>
        requires(std::is_arithmetic_v<Contribution> && !std::is_same_v<std::remove_cv_t<Contribution>, bool>)
    [[nodiscard]] static ImagePtr<Contribution> reconstructFromNodeContributions(
        const MorphologicalTree& tree, std::span<const Contribution> nodeContributions,
        const char* context = "TreeAltitudeAlgorithms::reconstructFromNodeContributions") {
        std::vector<Contribution> pixelValues = detail::tree_altitude::reconstructNodeContributionValues(tree, nodeContributions, context);
        ImagePtr<Contribution> image = Image<Contribution>::create(tree.numRows(), tree.numColumns());
        Contribution* pixels = image->rawData();
        std::copy(pixelValues.begin(), pixelValues.end(), pixels);
        return image;
    }

    /**
     * @brief Exports a live rooted topology and explicit altitudes to a compact parent/altitude representation.
     *
     * @details
     * This method owns only the structural valuedTree export. It reuses
     * `detail::computeExportedHigraLayout()` so the parent/altitude export has
     * the same compact id convention used by attribute-buffer projection in
     * `AttributeComputation`.
     *
     * @param tree Tree topology.
     * @param altitude Altitude data indexed by node identifier.
     * @return The exported live rooted topology and explicit altitudes to a compact parent/altitude representation.
     */
    template <AltitudeValue T>
    [[nodiscard]] static std::pair<std::vector<NodeId>, std::vector<T>> exportHigraHierarchy(const MorphologicalTree& tree, std::span<const T> altitude) {
        tree.requireNotEditing("TreeAltitudeAlgorithms::exportHigraHierarchy");
        const detail::ExportedHigraLayout layout = detail::computeExportedHigraLayout(tree, altitude);
        const int numLeaves = layout.numLeaves;
        const int numVertices = layout.numVertices;

        std::vector<NodeId> parent(static_cast<std::size_t>(numVertices), InvalidNode);
        std::vector<T> exportedAltitude(static_cast<std::size_t>(numVertices), T{});

        for (NodeId oldNodeId : layout.sortedNodes) {
            const NodeId newNodeId = layout.nodeToHigra[static_cast<std::size_t>(oldNodeId)];
            exportedAltitude[static_cast<std::size_t>(newNodeId)] = nodeAltitude(altitude, oldNodeId);
        }

        for (NodeId leafIndex = 0; leafIndex < numLeaves; ++leafIndex) {
            const PixelId pixel = layout.properParts[static_cast<std::size_t>(leafIndex)];
            const NodeId smallestNodeId = tree.smallestNode(pixel);
            if (smallestNodeId == InvalidNode || !tree.isAlive(smallestNodeId)) {
                throw std::runtime_error("Each proper part must belong to one alive node when exporting a compact Higra hierarchy.");
            }
            parent[static_cast<std::size_t>(leafIndex)] = layout.nodeToHigra[static_cast<std::size_t>(smallestNodeId)];
            exportedAltitude[static_cast<std::size_t>(leafIndex)] = nodeAltitude(altitude, smallestNodeId);
        }

        for (NodeId oldNodeId : layout.sortedNodes) {
            const NodeId newNodeId = layout.nodeToHigra[static_cast<std::size_t>(oldNodeId)];
            const NodeId oldParentNodeId = tree.parent(oldNodeId);
            parent[static_cast<std::size_t>(newNodeId)] =
                oldParentNodeId == oldNodeId ? newNodeId : layout.nodeToHigra[static_cast<std::size_t>(oldParentNodeId)];
        }

        return {std::move(parent), std::move(exportedAltitude)};
    }

    /**
     * @brief Validates the hierarchy's declared global altitude order.
     *
     * @param tree Tree topology.
     * @param altitude Altitude data indexed by node identifier.
     */
    template <AltitudeValue T> static void validateMonotoneNodeAltitudes(const MorphologicalTree& tree, std::span<const T> altitude) {
        validateNodeAltitudeBufferShape(tree, altitude);
        const NodeAltitudeOrder nodeAltitudeOrder = tree.nodeAltitudeOrder();
        bool increasingFromRoot = false;
        switch (nodeAltitudeOrder) {
        case NodeAltitudeOrder::Increasing:
            increasingFromRoot = true;
            break;
        case NodeAltitudeOrder::Decreasing:
            increasingFromRoot = false;
            break;
        case NodeAltitudeOrder::Unconstrained:
            return;
        }

        for (NodeId nodeId : tree.aliveNodeIds()) {
            if (tree.isRoot(nodeId)) {
                continue;
            }

            const NodeId parentNodeId = tree.parent(nodeId);
            if (parentNodeId == InvalidNode || !tree.isAlive(parentNodeId)) {
                throw std::runtime_error("Monotonic validation requires every alive non-root node to have an alive parent.");
            }

            if (increasingFromRoot) {
                if (nodeAltitude(altitude, parentNodeId) >= nodeAltitude(altitude, nodeId)) {
                    throw std::runtime_error("Hierarchy altitude buffer must be strictly increasing from parent to child.");
                }
            } else if (nodeAltitude(altitude, parentNodeId) <= nodeAltitude(altitude, nodeId)) {
                throw std::runtime_error("Hierarchy altitude buffer must be strictly decreasing from parent to child.");
            }
        }
    }
};

} // namespace mmcfilters
