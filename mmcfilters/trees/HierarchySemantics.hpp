#pragma once

#include "../utils/RegularGridAdjacency2D.hpp"

#include <optional>
#include <stdexcept>
#include <utility>
#include <variant>

namespace mmcfilters {

/**
 * @brief Descriptive family label used for diagnostics and introspection.
 *
 * Algorithms should query explicit hierarchy capabilities such as
 * `AltitudeOrder` and `AdjacencyMode` instead of dispatching on this label.
 */
enum class MorphologicalTreeKind { MAX_TREE = 0, MIN_TREE = 1, TREE_OF_SHAPES = 2, SELF_DUAL_RESIDUAL_TREE = 3, GENERIC = 4 };

/**
 * @brief Global ordering constraint of node altitudes along parent/child arcs.
 */
enum class AltitudeOrder {
    /// Every alive parent-child arc satisfies `altitude(parent) < altitude(child)`.
    INCREASING_FROM_ROOT,
    /// Every alive parent-child arc satisfies `altitude(parent) > altitude(child)`.
    DECREASING_FROM_ROOT,
    /// No global parent-child altitude order is declared.
    UNCONSTRAINED
};

/**
 * @brief Shape of the adjacency context attached to a hierarchy.
 */
enum class AdjacencyMode { NONE, UNIFORM, DIRECTIONAL };

/** @brief Marks a hierarchy with no attached adjacency context. */
struct NoAdjacency {};

/** @brief One adjacency relation applies uniformly to the hierarchy. */
struct UniformGridAdjacency2D {
    /// Immutable regular-grid relation shared by every hierarchy arc.
    RegularGridAdjacency2D relation;
};

/**
 * @brief Separate adjacency relations apply to decreasing and increasing arcs.
 */
struct DirectionalGridAdjacency2D {
    /// Relation used when following decreasing-altitude arcs.
    RegularGridAdjacency2D decreasing;
    /// Relation used when following increasing-altitude arcs.
    RegularGridAdjacency2D increasing;
};

using AdjacencyContext = std::variant<NoAdjacency, UniformGridAdjacency2D, DirectionalGridAdjacency2D>;

/**
 * @brief Generic semantic capabilities associated with a morphological hierarchy.
 *
 * The topology and proper-part ownership remain the primary tree contract.
 * These capabilities describe optional interpretation rules without changing
 * the structural representation or specializing it for a concrete builder.
 */
struct HierarchySemantics {
    /// Descriptive family label; algorithms must dispatch on capabilities instead.
    MorphologicalTreeKind descriptiveKind = MorphologicalTreeKind::GENERIC;
    /// Global altitude-order capability of the hierarchy.
    AltitudeOrder altitudeOrder = AltitudeOrder::UNCONSTRAINED;
    /// Optional uniform or directional adjacency context.
    AdjacencyContext adjacency = NoAdjacency{};

    /**
     * @brief Returns the shape of the attached adjacency context.
     *
     * @return The shape of the attached adjacency context.
     */
    [[nodiscard]] AdjacencyMode adjacencyMode() const noexcept {
        if (std::holds_alternative<UniformGridAdjacency2D>(adjacency)) {
            return AdjacencyMode::UNIFORM;
        }
        if (std::holds_alternative<DirectionalGridAdjacency2D>(adjacency)) {
            return AdjacencyMode::DIRECTIONAL;
        }
        return AdjacencyMode::NONE;
    }
};

/**
 * @brief Maps a descriptive family label to its default global altitude order.
 *
 * @param kind Morphological-tree family.
 * @return The mapped descriptive family label to its default global altitude order.
 */
inline AltitudeOrder defaultAltitudeOrder(MorphologicalTreeKind kind) noexcept {
    switch (kind) {
    case MorphologicalTreeKind::MAX_TREE:
        return AltitudeOrder::INCREASING_FROM_ROOT;
    case MorphologicalTreeKind::MIN_TREE:
        return AltitudeOrder::DECREASING_FROM_ROOT;
    case MorphologicalTreeKind::GENERIC:
    case MorphologicalTreeKind::SELF_DUAL_RESIDUAL_TREE:
    case MorphologicalTreeKind::TREE_OF_SHAPES:
        return AltitudeOrder::UNCONSTRAINED;
    }
    return AltitudeOrder::UNCONSTRAINED;
}

/**
 * @brief Constructs semantic capabilities from a family label and adjacency.
 *
 * Producer and import adapters use this helper to select the conventional
 * altitude order for a known family while keeping the tree representation
 * independent of the producer.
 *
 * @param kind Morphological-tree family.
 * @param adjacency Adjacency relation used by the operation.
 * @param directionalAdjacency Whether the hierarchy uses directional adjacency.
 * @return The constructed semantic capabilities from a family label and adjacency.
 */
inline HierarchySemantics makeHierarchySemantics(MorphologicalTreeKind kind, std::optional<RegularGridAdjacency2D> adjacency = std::nullopt,
                                                 std::optional<DirectionalGridAdjacency2D> directionalAdjacency = std::nullopt) {
    if (adjacency && directionalAdjacency) {
        throw std::invalid_argument("A hierarchy cannot use uniform and directional adjacency contexts simultaneously.");
    }

    AdjacencyContext context = NoAdjacency{};
    if (directionalAdjacency) {
        context = std::move(*directionalAdjacency);
    } else if (adjacency) {
        context = UniformGridAdjacency2D{std::move(*adjacency)};
    }

    return HierarchySemantics{kind, defaultAltitudeOrder(kind), std::move(context)};
}

} // namespace mmcfilters
