#pragma once

#include "MorphologicalTreeSemantics.hpp"
#include "ProperPartDomain.hpp"
#include "../utils/Common.hpp"

#include <optional>
#include <utility>
#include <vector>

namespace mmcfilters {

/**
 * @brief Owning, move-only native hierarchy returned by builders and consumed by materialization.
 *
 * Builders move their existing buffers into this representation.
 * Higra import conversion also produces this representation after validation.
 * Materialization transfers the buffers into the valued-tree owner without
 * another copy. The constructor transfers ownership without structural
 * validation; callers remain responsible for providing valid hierarchy data.
 * This representation also supports producer scalar types outside AltitudeValue; these
 * results cannot be materialized as valued trees without an explicit conversion.
 *
 * @tparam T Scalar type stored by the producer; valued-tree materialization requires AltitudeValue.
 */
template <class T> struct NativeHierarchyOwned {
    /** @brief Parent identifier indexed by internal node. */
    std::vector<NodeId> parent;
    /** @brief Smallest node identifier indexed by pixel. */
    std::vector<NodeId> smallestNodeMap;
    /** @brief Altitude indexed by internal node. */
    std::vector<T> nodeAltitudes;
    /** @brief Identifier of the self-parented root. */
    NodeId root;
    /** @brief Optional two-dimensional pixel domain. */
    std::optional<GridDomain2D> gridDomain2D;
    /** @brief Semantic metadata, including owned construction adjacencies. */
    MorphologicalTreeSemantics semantics;

    /**
     * @brief Takes ownership of native hierarchy buffers without structural validation.
     *
     * @param parent Parent identifier indexed by internal node.
     * @param smallestNodeMap Smallest node identifier indexed by pixel.
     * @param nodeAltitudes Altitude indexed by internal node.
     * @param root Identifier of the root.
     * @param gridDomain2D Optional two-dimensional pixel domain.
     * @param semantics Semantic metadata to transfer with the buffers.
     */
    NativeHierarchyOwned(std::vector<NodeId>&& parent, std::vector<NodeId>&& smallestNodeMap, std::vector<T>&& nodeAltitudes, NodeId root,
                         std::optional<GridDomain2D> gridDomain2D, MorphologicalTreeSemantics semantics)
        : parent(std::move(parent)), smallestNodeMap(std::move(smallestNodeMap)), nodeAltitudes(std::move(nodeAltitudes)), root(root),
          gridDomain2D(gridDomain2D), semantics(std::move(semantics)) {}

    /**
     * @brief Disables copy construction.
     */
    NativeHierarchyOwned(const NativeHierarchyOwned&) = delete;
    /**
     * @brief Disables copy assignment.
     */
    NativeHierarchyOwned& operator=(const NativeHierarchyOwned&) = delete;
    /**
     * @brief Transfers ownership of the native hierarchy.
     */
    NativeHierarchyOwned(NativeHierarchyOwned&&) noexcept = default;
    /**
     * @brief Assigns the supplied object state.
     *
     * @return Reference to the assigned object.
     */
    NativeHierarchyOwned& operator=(NativeHierarchyOwned&&) noexcept = default;
};

/**
 * @brief Creates an owning native hierarchy by moving its buffers without structural validation.
 *
 * @param parent Parent node indexed by internal node identifier.
 * @param smallestNodeMap Smallest node indexed by pixel identifier.
 * @param nodeAltitudes Altitude values indexed by internal node identifier.
 * @param root Root node of the operation.
 * @param gridDomain2D Two-dimensional grid metadata associated with the hierarchy.
 * @param semantics Semantic metadata associated with the hierarchy.
 * @return Owned native hierarchy.
 */
template <class T>
[[nodiscard]] NativeHierarchyOwned<T> makeNativeHierarchyOwned(std::vector<NodeId>&& parent, std::vector<NodeId>&& smallestNodeMap,
                                                             std::vector<T>&& nodeAltitudes, NodeId root, std::optional<GridDomain2D> gridDomain2D,
                                                             MorphologicalTreeSemantics semantics) {
    return NativeHierarchyOwned<T>(std::move(parent), std::move(smallestNodeMap), std::move(nodeAltitudes), root, gridDomain2D, std::move(semantics));
}

} // namespace mmcfilters
