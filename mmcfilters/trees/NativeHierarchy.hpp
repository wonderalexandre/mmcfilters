#pragma once

#include "HierarchySemantics.hpp"
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
 * the same independent node, proper-part, and altitude domains plus optional
 * interpretation capabilities. Materialization copies these synchronous views
 * into the owning weighted tree.
 */
template <AltitudeValue T> struct NativeHierarchyView {
    /// Parent id for every internal-node slot.
    std::span<const NodeId> nodeParent;
    /// Direct owning node for every proper part.
    std::span<const NodeId> properPartOwner;
    /// Altitude value for every internal-node slot.
    std::span<const T> altitude;
    /// Root node id.
    NodeId root = InvalidNode;
    /// Optional row/column interpretation of proper-part ids.
    std::optional<GridDomain2D> gridDomain2D;
    /// Altitude-order and adjacency capabilities of the hierarchy.
    HierarchySemantics semantics;
};

} // namespace mmcfilters
