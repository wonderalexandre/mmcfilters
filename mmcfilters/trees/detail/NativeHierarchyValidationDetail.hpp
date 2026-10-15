#pragma once

#include "../NativeHierarchy.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mmcfilters::detail {

/**
 * @brief Owns native hierarchy buffers supplied by an internal producer or a validated import.
 *
 * @tparam T Altitude type stored by the hierarchy.
 */
template <AltitudeValue T> struct OwnedNativeHierarchyStorage {
    /** @brief Dense node identifier of the node parent. */
    std::vector<NodeId> parent;
    /** @brief Dense identifier of the smallest node. */
    std::vector<NodeId> smallestNodeMap;
    /** @brief Stores node altitudes. */
    std::vector<T> nodeAltitudes;
    /** @brief Dense node identifier of the root. */
    NodeId root = InvalidNode;
    /** @brief Grid domain2 d. */
    std::optional<GridDomain2D> gridDomain2D;
    /** @brief Semantics. */
    MorphologicalTreeSemantics semantics;
};

/**
 * @brief Owning, move-only native hierarchy for factory materialization.
 *
 * Internal producers move their existing buffers into this representation.
 * Public span imports are validated first and copied once into it. In both
 * cases materialization transfers the buffers into the valued-tree owner
 * without another copy. Internal producers are responsible for structural validity.
 */
template <AltitudeValue T> class OwnedNativeHierarchy {
    /** @brief Storage. */
    OwnedNativeHierarchyStorage<T> storage_;

  public:
    /**
     * @brief Takes ownership of native hierarchy storage without structural validation.
     *
     * @param storage Storage object that owns the native hierarchy arrays.
     */
    explicit OwnedNativeHierarchy(OwnedNativeHierarchyStorage<T>&& storage) : storage_(std::move(storage)) {}

    /**
     * @brief Disables copy construction.
     */
    OwnedNativeHierarchy(const OwnedNativeHierarchy&) = delete;
    /**
     * @brief Disables copy assignment.
     */
    OwnedNativeHierarchy& operator=(const OwnedNativeHierarchy&) = delete;
    /**
     * @brief Transfers ownership of the native hierarchy.
     */
    OwnedNativeHierarchy(OwnedNativeHierarchy&&) noexcept = default;
    /**
     * @brief Assigns the supplied object state.
     *
     * @return Reference to the assigned object.
     */
    OwnedNativeHierarchy& operator=(OwnedNativeHierarchy&&) noexcept = default;

    /**
     * @brief Releases ownership of the stored native hierarchy.
     *
     * @return Owned native hierarchy storage.
     */
    [[nodiscard]] OwnedNativeHierarchyStorage<T> release() && { return std::move(storage_); }
};

/**
 * @brief Complete structural validation for external native hierarchy imports.
 */
class NativeHierarchyValidation {
    /**
     * @brief Validates descriptor.
     *
     * @param hierarchy Native hierarchy data to validate or convert.
     */
    template <AltitudeValue T> static void validateDescriptor(const NativeHierarchyView<T>& hierarchy) {
        if (hierarchy.parent.size() > static_cast<std::size_t>(std::numeric_limits<NodeId>::max()) ||
            hierarchy.smallestNodeMap.size() > static_cast<std::size_t>(std::numeric_limits<NodeId>::max())) {
            throw std::invalid_argument("Native hierarchy domain exceeds the NodeId range.");
        }
        if (hierarchy.parent.empty()) {
            throw std::invalid_argument("Native hierarchy requires at least one internal node.");
        }
        if (hierarchy.smallestNodeMap.empty()) {
            throw std::invalid_argument("Native hierarchy requires at least one proper part.");
        }
        if (hierarchy.nodeAltitudes.size() != hierarchy.parent.size()) {
            throw std::invalid_argument("Native hierarchy altitude size must match the internal-node domain.");
        }

        const int numNodes = static_cast<int>(hierarchy.parent.size());
        if (hierarchy.root < 0 || hierarchy.root >= numNodes) {
            throw std::invalid_argument("Native hierarchy requires a valid root node id.");
        }
        if (hierarchy.gridDomain2D && hierarchy.smallestNodeMap.size() != hierarchy.gridDomain2D->size("Native hierarchy 2D pixel domain")) {
            throw std::invalid_argument("Native hierarchy pixel domain must match the attached 2D grid.");
        }

        switch (hierarchy.semantics.kind) {
        case MorphologicalTreeKind::Generic:
        case MorphologicalTreeKind::MaxTree:
        case MorphologicalTreeKind::MinTree:
        case MorphologicalTreeKind::TreeOfShapes:
        case MorphologicalTreeKind::ExtremumResidualTree:
            break;
        default:
            throw std::invalid_argument("Native hierarchy kind is not supported.");
        }

        validateMorphologicalTreeSemantics(hierarchy.semantics);
        const auto requireMatchingGrid = [&hierarchy](const RegularGridAdjacency2D& adjacency, const char* context) {
            if (!hierarchy.gridDomain2D) {
                throw std::invalid_argument(std::string(context) + " requires an attached 2D grid domain.");
            }
            if (adjacency.getNumRows() != hierarchy.gridDomain2D->rows || adjacency.getNumColumns() != hierarchy.gridDomain2D->columns) {
                throw std::invalid_argument(std::string(context) + " must match the attached 2D grid.");
            }
        };
        if (const auto* context = std::get_if<SharedAdjacencyContext>(&hierarchy.semantics.constructionContext)) {
            requireMatchingGrid(context->adjacency, "SharedAdjacencyContext adjacency");
        } else if (const auto* convention = std::get_if<TopographicConvention>(&hierarchy.semantics.constructionContext)) {
            if (convention->infinityPixel < 0) {
                throw std::invalid_argument("TopographicConvention infinity pixel must be non-negative.");
            }
            if (!hierarchy.gridDomain2D) {
                throw std::invalid_argument("TopographicConvention requires an attached 2D grid domain.");
            }
            const std::int64_t extension = convention->domainExtension == TopographicDomainExtension::ExteriorRing ? 1 : -1;
            const std::int64_t activeRows = 2 * static_cast<std::int64_t>(hierarchy.gridDomain2D->rows) + extension;
            const std::int64_t activeColumns = 2 * static_cast<std::int64_t>(hierarchy.gridDomain2D->columns) + extension;
            if (activeRows <= 0 || activeColumns <= 0 || static_cast<std::int64_t>(convention->infinityPixel) >= activeRows * activeColumns) {
                throw std::invalid_argument("TopographicConvention infinity pixel must belong to the active topographic domain.");
            }
            if (const auto* immersion = std::get_if<ComplementaryGridImmersion>(&convention->immersion)) {
                requireMatchingGrid(immersion->complementaryAdjacencies.minAdjacency, "Topographic minimum adjacency");
                requireMatchingGrid(immersion->complementaryAdjacencies.maxAdjacency, "Topographic maximum adjacency");
            }
        }
    }

  public:
    /**
     * @brief Validates hierarchy dimensions and semantic metadata without scanning topology buffers.
     *
     * @param hierarchy Native hierarchy descriptor.
     */
    template <AltitudeValue T> static void validateDescriptorOnly(const NativeHierarchyView<T>& hierarchy) { validateDescriptor(hierarchy); }

    /**
     * @brief Performs the complete structural validation required for a public
     * raw-buffer import.
     *
     * @param hierarchy Hierarchy data.
     */
    template <AltitudeValue T> static void validateComplete(const NativeHierarchyView<T>& hierarchy) {
        validateDescriptor(hierarchy);

        const int numNodes = static_cast<int>(hierarchy.parent.size());
        std::vector<int> remainingChildren(hierarchy.parent.size(), 0);

        for (NodeId node = 0; node < numNodes; ++node) {
            const NodeId parent = hierarchy.parent[static_cast<std::size_t>(node)];
            if (node == hierarchy.root) {
                if (parent != node) {
                    throw std::invalid_argument("Native hierarchy root must be self-parented.");
                }
                continue;
            }
            if (parent < 0 || parent >= numNodes) {
                throw std::invalid_argument("Native hierarchy parent is outside the internal-node domain.");
            }
            if (parent == node) {
                throw std::invalid_argument("Native hierarchy contains a detached self-parented non-root node.");
            }
            ++remainingChildren[static_cast<std::size_t>(parent)];
        }

        std::vector<std::uint8_t> subtreeHasSupport(hierarchy.parent.size(), std::uint8_t{0});
        for (PixelId pixel = 0; pixel < static_cast<PixelId>(hierarchy.smallestNodeMap.size()); ++pixel) {
            const NodeId smallestNode = hierarchy.smallestNodeMap[static_cast<std::size_t>(pixel)];
            if (smallestNode < 0 || smallestNode >= numNodes) {
                throw std::invalid_argument("Native hierarchy smallest node is outside the internal-node domain.");
            }
            subtreeHasSupport[static_cast<std::size_t>(smallestNode)] = 1;
        }

        std::vector<NodeId> ready;
        ready.reserve(hierarchy.parent.size());
        for (NodeId node = 0; node < numNodes; ++node) {
            if (remainingChildren[static_cast<std::size_t>(node)] == 0) {
                ready.push_back(node);
            }
        }

        std::size_t processedNodes = 0;
        for (std::size_t index = 0; index < ready.size(); ++index) {
            const NodeId node = ready[index];
            if (subtreeHasSupport[static_cast<std::size_t>(node)] == 0) {
                throw std::invalid_argument("Native hierarchy contains a node whose subtree support is empty.");
            }
            ++processedNodes;
            if (node == hierarchy.root) {
                continue;
            }

            const NodeId parent = hierarchy.parent[static_cast<std::size_t>(node)];
            subtreeHasSupport[static_cast<std::size_t>(parent)] = 1;
            int& parentRemaining = remainingChildren[static_cast<std::size_t>(parent)];
            --parentRemaining;
            if (parentRemaining == 0) {
                ready.push_back(parent);
            }
        }
        if (processedNodes != hierarchy.parent.size()) {
            throw std::invalid_argument("Native hierarchy contains a cycle or nodes disconnected from the root.");
        }
    }
};

/**
 * @brief Validates a native hierarchy and copies it into owned storage.
 *
 * @param hierarchy Native hierarchy data to validate or convert.
 * @return Owning copy of the validated native hierarchy.
 */
template <AltitudeValue T> [[nodiscard]] OwnedNativeHierarchy<T> validateAndCopyNativeHierarchy(NativeHierarchyView<T> hierarchy) {
    NativeHierarchyValidation::validateComplete(hierarchy);
    return OwnedNativeHierarchy<T>(
        OwnedNativeHierarchyStorage<T>{std::vector<NodeId>(hierarchy.parent.begin(), hierarchy.parent.end()),
                                       std::vector<NodeId>(hierarchy.smallestNodeMap.begin(), hierarchy.smallestNodeMap.end()),
                                       std::vector<T>(hierarchy.nodeAltitudes.begin(), hierarchy.nodeAltitudes.end()), hierarchy.root, hierarchy.gridDomain2D,
                                       std::move(hierarchy.semantics)});
}

/**
 * @brief Takes ownership of a trusted internal producer result without structural validation.
 *
 * @param parent Parent node indexed by internal node identifier.
 * @param smallestNodeMap Smallest node indexed by pixel identifier.
 * @param nodeAltitudes Altitude values indexed by internal node identifier.
 * @param root Root node of the operation.
 * @param gridDomain2D Two-dimensional grid metadata associated with the hierarchy.
 * @param semantics Semantic metadata associated with the hierarchy.
 * @return Owned native hierarchy.
 */
template <AltitudeValue T>
[[nodiscard]] OwnedNativeHierarchy<T> makeOwnedNativeHierarchy(std::vector<NodeId>&& parent, std::vector<NodeId>&& smallestNodeMap,
                                                             std::vector<T>&& nodeAltitudes, NodeId root, std::optional<GridDomain2D> gridDomain2D,
                                                             MorphologicalTreeSemantics semantics) {
    return OwnedNativeHierarchy<T>(OwnedNativeHierarchyStorage<T>{std::move(parent), std::move(smallestNodeMap), std::move(nodeAltitudes), root,
                                                                 gridDomain2D, std::move(semantics)});
}

} // namespace mmcfilters::detail
