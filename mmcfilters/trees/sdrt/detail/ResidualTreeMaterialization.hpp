#pragma once

/**
 * @file ResidualTreeMaterialization.hpp
 * @brief Validation and native-buffer materialization for residual-tree events.
 */

#include "ResidualTreeEventAssembler.hpp"
#include "../../HierarchySemantics.hpp"
#include "../../detail/NativeHierarchyValidationDetail.hpp"
#include "../../../utils/Image.hpp"
#include "../../../utils/RegularGridAdjacency2D.hpp"

#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcfilters::sdrt::detail {

/** @brief Validated residual-tree buffers ready to be retained by a builder. */
template <AltitudeValue T> struct MaterializedResidualTree {
    int rows = 0;                                      ///< Source image rows.
    int cols = 0;                                      ///< Source image columns.
    NodeId root = InvalidNode;                         ///< Dense residual-tree root.
    std::vector<NodeId> nodeParent;                    ///< Parent indexed by residual node.
    std::vector<NodeId> properPartOwner;               ///< Direct owner indexed by pixel.
    std::vector<T> altitude;                           ///< Altitude indexed by residual node.
    mmcfilters::detail::NativeTopologyProof topologyProof; ///< Established native topology proof.
};

/**
 * @brief Validates assembler output and returns native residual-tree buffers.
 * @param image Image reconstructed by the output hierarchy.
 * @param adjacency Adjacency semantics assigned to the hierarchy.
 * @param output Finalized event-assembler output.
 * @return Validated materialized buffers.
 */
template <AltitudeValue T>
[[nodiscard]] MaterializedResidualTree<T> materializeResidualTree(const ImagePtr<T>& image, const RegularGridAdjacency2D& adjacency,
                                                                  typename ResidualTreeEventAssembler<T>::Output output) {
    const NativeHierarchyView<T> hierarchy{output.nodeParent,
                                           output.properPartOwner,
                                           output.altitude,
                                           0,
                                           GridDomain2D{image->getNumRows(), image->getNumCols()},
                                           makeHierarchySemantics(MorphologicalTreeKind::SELF_DUAL_RESIDUAL_TREE, adjacency)};
    auto topologyProof = mmcfilters::detail::NativeHierarchyValidation::validateWithEstablishedSupport(
        hierarchy, std::move(output.subtreeSupportProof), [&image, &output](NodeId properPart, NodeId owner) {
            if (output.altitude[static_cast<std::size_t>(owner)] != (*image)[properPart]) {
                throw std::runtime_error("Min/max residual assembler reconstruction differs from the input image.");
            }
        });
    return MaterializedResidualTree<T>{image->getNumRows(), image->getNumCols(), NodeId{0}, std::move(output.nodeParent),
                                       std::move(output.properPartOwner), std::move(output.altitude), std::move(topologyProof)};
}

} // namespace mmcfilters::sdrt::detail
