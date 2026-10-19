#pragma once

#include "MorphologicalTreeSemantics.hpp"
#include "ProperPartDomain.hpp"
#include "../utils/Altitude.hpp"
#include "../utils/Common.hpp"

#include <optional>
#include <span>

namespace mmcfilters {

/**
 * @brief View of borrowed native hierarchy buffers for internal conversion and validation.
 *
 * The three spans borrow node-parent, pixel-owner, and node-altitude buffers.
 * Those buffers must remain valid while a consumer uses the view. Domain and
 * semantic metadata are stored by value. This descriptor does not provide a
 * public import entry point; external hierarchy import uses the Higra format.
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
