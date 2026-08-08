#pragma once

#include "../HierarchySemantics.hpp"
#include "HigraImportLayoutDetail.hpp"
#include "NativeHierarchyValidationDetail.hpp"

#include <cmath>
#include <optional>
#include <sstream>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace mmcfilters::detail {

/**
 * @brief Result of adapting a compact leaf-first Higra hierarchy.
 *
 * `internalNodeOffset` records the affine mapping from dense internal node ids
 * back to the imported external id domain. It is interoperability provenance,
 * not morphological-tree semantics.
 */
template <AltitudeValue T> struct HigraHierarchyImport {
    /** @brief Stores the hierarchy. */
    ValidatedNativeHierarchy<T> hierarchy;
    /** @brief Stores the internal node offset. */
    NodeId internalNodeOffset = InvalidNode;
};

/**
 * @brief Converts Higra parent/altitude arrays into the generic native format.
 *
 * The adapter owns all leaf-first id policy. `MorphologicalTree` subsequently
 * materializes only independent node-parent and proper-part-owner domains.
 *
 * @param parent Parent node used by the operation.
 * @param higraAltitude Altitude or level represented by `higraAltitude`.
 * @param rows Number of rows in the domain.
 * @param cols Number of columns in the domain.
 * @param kind Morphological-tree family.
 * @param adjacency Adjacency relation used by the operation.
 * @return The converted Higra parent/altitude arrays into the generic native format.
 */
template <AltitudeValue T>
[[nodiscard]] HigraHierarchyImport<T> adaptHigraHierarchy(std::span<const NodeId> parent, std::span<const T> higraAltitude, int rows, int cols,
                                                          MorphologicalTreeKind kind, std::optional<RegularGridAdjacency2D> adjacency) {
    if (kind == MorphologicalTreeKind::SELF_DUAL_RESIDUAL_TREE) {
        throw std::invalid_argument(
            "Higra import does not accept SELF_DUAL_RESIDUAL_TREE because residual-tree proper parts are not represented by the Higra leaf layout.");
    }
    if (!adjacency && kind != MorphologicalTreeKind::TREE_OF_SHAPES) {
        throw std::invalid_argument("Higra import of max/min trees requires an explicit adjacency relation.");
    }

    const GridDomain2D gridDomain{rows, cols};
    HierarchySemantics semantics = makeHierarchySemantics(kind, std::move(adjacency));
    const NodeId numProperParts = static_cast<NodeId>(gridDomain.size("Higra import 2D proper-part domain"));
    if (const auto* uniformAdjacency = std::get_if<UniformGridAdjacency2D>(&semantics.adjacency)) {
        if (uniformAdjacency->relation.getNumRows() * uniformAdjacency->relation.getNumCols() != numProperParts) {
            throw std::invalid_argument("Higra leaf count must match image domain size.");
        }
    }

    const HigraImportLayout layout(numProperParts, parent.size());
    if (higraAltitude.size() != static_cast<std::size_t>(layout.numVertices())) {
        throw std::invalid_argument("Higra altitude buffer size must match the preserved imported Higra hierarchy.");
    }
    const auto validateFiniteAltitude = [&](NodeId externalNode) {
        if constexpr (std::is_floating_point_v<T>) {
            const long double level = static_cast<long double>(higraAltitude[static_cast<std::size_t>(externalNode)]);
            if (!std::isfinite(level)) {
                std::ostringstream message;
                message << "WeightedMorphologicalTree Higra altitude input requires finite floating-point altitudes; value at index " << externalNode << " is "
                        << level << ".";
                throw std::invalid_argument(message.str());
            }
        }
    };

    std::vector<NodeId> properPartOwner(static_cast<std::size_t>(numProperParts), InvalidNode);
    for (NodeId properPart = 0; properPart < numProperParts; ++properPart) {
        validateFiniteAltitude(properPart);
        properPartOwner[static_cast<std::size_t>(properPart)] =
            layout.internalNodeId(parent[static_cast<std::size_t>(properPart)], "Each Higra leaf must point to an internal node.");
    }

    const NodeId numNodes = layout.numInternalNodes();
    std::vector<NodeId> nodeParent(static_cast<std::size_t>(numNodes), InvalidNode);
    std::vector<T> altitude(static_cast<std::size_t>(numNodes));
    NodeId root = InvalidNode;
    for (NodeId node = 0; node < numNodes; ++node) {
        const NodeId higraNode = layout.higraNodeId(node);
        const NodeId parentHigraNode = parent[static_cast<std::size_t>(higraNode)];
        validateFiniteAltitude(higraNode);
        altitude[static_cast<std::size_t>(node)] = higraAltitude[static_cast<std::size_t>(higraNode)];

        if (HigraImportLayout::isRootLink(higraNode, parentHigraNode)) {
            if (root != InvalidNode) {
                throw std::invalid_argument("A Higra hierarchy must encode exactly one root.");
            }
            nodeParent[static_cast<std::size_t>(node)] = node;
            root = node;
            continue;
        }

        nodeParent[static_cast<std::size_t>(node)] =
            layout.internalNodeId(parentHigraNode, "Each Higra internal node must point to another internal node or to itself.");
    }
    if (root == InvalidNode) {
        throw std::invalid_argument("A Higra hierarchy must encode exactly one root.");
    }

    const NativeHierarchyView<T> view{nodeParent, properPartOwner, altitude, root, gridDomain, semantics};
    NativeTopologyProof proof = NativeHierarchyValidation::validateComplete(view);
    return {makeValidatedNativeHierarchy<T>(std::move(nodeParent), std::move(properPartOwner), std::move(altitude), root, gridDomain, std::move(semantics),
                                            std::move(proof)),
            layout.numLeaves()};
}

} // namespace mmcfilters::detail
