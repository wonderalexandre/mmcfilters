#pragma once

#include "MorphologicalTreeSemantics.hpp"
#include "ProperPartDomain.hpp"
#include "../utils/Altitude.hpp"
#include "../utils/Common.hpp"

#include <optional>
#include <span>

namespace mmcfilters {

/**
 * @brief Typed non-owning output contract shared by native hierarchy producers.
 *
 * Builders may use any internal algorithm. At the factory boundary they expose
 * the same independent node, pixel, and altitude domains plus optional
 * interpretation capabilities. Materialization copies these synchronous views
 * into the owning valued tree.
 */
template <AltitudeValue T> struct NativeHierarchyView {
    /// Parent id for every internal-node slot.
    std::span<const NodeId> parent;
    /// Smallest node containing each pixel.
    std::span<const NodeId> smallestNodeMap;
    /// Altitude value for every internal-node slot.
    std::span<const T> nodeAltitudes;
    /// Root node id.
    NodeId root = InvalidNode;
    /// Optional row/column interpretation of pixel ids.
    std::optional<GridDomain2D> gridDomain2D;
    /// Altitude-order and adjacency capabilities of the hierarchy.
    MorphologicalTreeSemantics semantics;
};

} // namespace mmcfilters
