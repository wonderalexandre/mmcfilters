#pragma once

#include "../HierarchySemantics.hpp"
#include "ComponentTreeUnionFind.hpp"
#include "NativeHierarchyValidationDetail.hpp"
#include "../../utils/Image.hpp"

#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcfilters::detail {

namespace kernel {

/**
 * @brief Polarity selected by the concrete image component-tree producer.
 */
enum class ComponentTreePolarity { MAX_TREE, MIN_TREE };

/**
 * @brief Concrete max-tree/min-tree producer over a regular 2D image.
 *
 * The producer preserves the historical union-find ordering and translates
 * its pixel-parent representation directly into the generic native hierarchy
 * buffers consumed by `MorphologicalTreeFactory`.
 */
template <AltitudeValue T> class ComponentTreeProducer {
    /** @brief Stores the polarity. */
    ComponentTreePolarity polarity_;
    /** @brief Stores the adjacency. */
    RegularGridAdjacency2D adjacency_;

  public:
    /**
     * @brief Creates a component-tree producer using an explicit adjacency.
     *
     * @param polarity Polarity that determines whether maxima or minima are processed.
     * @param adjacency Adjacency relation used by the operation.
     */
    ComponentTreeProducer(ComponentTreePolarity polarity, RegularGridAdjacency2D adjacency) : polarity_(polarity), adjacency_(std::move(adjacency)) {}

    /**
     * @brief Builds proven owning topology, ownership, and altitude buffers.
     *
     * Internal node ids are assigned in exactly the historical ordered-pixel
     * materialization order. The altitude representative is the smallest
     * row-major direct proper part owned by each node, matching the previous
     * linked-list based inference even for equal floating-point values with
     * different object representations.
     *
     * @param image Image used by the operation.
     * @return The resulting proven owning topology, ownership, and altitude buffers.
     */
    [[nodiscard]] ValidatedNativeHierarchy<T> build(const ImagePtr<T>& image) {
        assert(polarity_ == ComponentTreePolarity::MAX_TREE || polarity_ == ComponentTreePolarity::MIN_TREE);
        const bool isMaxTree = polarity_ == ComponentTreePolarity::MAX_TREE;
        const MorphologicalTreeKind descriptiveKind = isMaxTree ? MorphologicalTreeKind::MAX_TREE : MorphologicalTreeKind::MIN_TREE;

        const GridDomain2D gridDomain{image->getNumRows(), image->getNumCols()};
        RegularGridAdjacency2D adjacency = std::move(adjacency_);
        assert(adjacency.getNumRows() == gridDomain.rows && adjacency.getNumCols() == gridDomain.cols);
        ComponentTreeUnionFind unionFind(&adjacency, isMaxTree);
        auto [pixelParent, orderedPixels, numBuiltNodes] = unionFind.template build<T>(image);

        const int numPixels = image->getSize();
        const auto* imageValues = image->rawData();
        std::vector<NodeId> nodeParent;
        std::vector<T> altitude;
        std::vector<NodeId> representativeProperPart;
        nodeParent.reserve(static_cast<std::size_t>(numBuiltNodes));
        altitude.reserve(static_cast<std::size_t>(numBuiltNodes));
        representativeProperPart.reserve(static_cast<std::size_t>(numBuiltNodes));
        std::vector<NodeId> properPartOwner(static_cast<std::size_t>(numPixels), InvalidNode);

        TopologicalNativeHierarchyRecorder proofRecorder(static_cast<std::size_t>(numBuiltNodes), static_cast<std::size_t>(numPixels), NodeId{0});

        NodeId root = InvalidNode;
        for (int index = 0; index < numPixels; ++index) {
            const NodeId properPart = orderedPixels[static_cast<std::size_t>(index)];
            NodeId owner = InvalidNode;

            if (properPart == pixelParent[static_cast<std::size_t>(properPart)]) {
                assert(root == InvalidNode);
                owner = static_cast<NodeId>(nodeParent.size());
                root = owner;
                nodeParent.push_back(owner);
                altitude.push_back(imageValues[static_cast<std::size_t>(properPart)]);
                representativeProperPart.push_back(properPart);
                proofRecorder.recordSupportedNode(owner, owner);
            } else {
                const NodeId parentProperPart = pixelParent[static_cast<std::size_t>(properPart)];
                const NodeId parentOwner = properPartOwner[static_cast<std::size_t>(parentProperPart)];
                if (imageValues[static_cast<std::size_t>(properPart)] != imageValues[static_cast<std::size_t>(parentProperPart)]) {
                    owner = static_cast<NodeId>(nodeParent.size());
                    nodeParent.push_back(parentOwner);
                    altitude.push_back(imageValues[static_cast<std::size_t>(properPart)]);
                    representativeProperPart.push_back(properPart);
                    proofRecorder.recordSupportedNode(owner, parentOwner);
                } else {
                    owner = parentOwner;
                }
            }

            properPartOwner[static_cast<std::size_t>(properPart)] = owner;
            if (properPart < representativeProperPart[static_cast<std::size_t>(owner)]) {
                representativeProperPart[static_cast<std::size_t>(owner)] = properPart;
                altitude[static_cast<std::size_t>(owner)] = imageValues[static_cast<std::size_t>(properPart)];
            }
            proofRecorder.recordProperPartOwner(owner);
        }

        assert(root == 0 && nodeParent.size() == static_cast<std::size_t>(numBuiltNodes));

        return makeValidatedNativeHierarchy<T>(std::move(nodeParent), std::move(properPartOwner), std::move(altitude), root, gridDomain,
                                               makeHierarchySemantics(descriptiveKind, std::move(adjacency)), std::move(proofRecorder).finish());
    }
};

} // namespace kernel

} // namespace mmcfilters::detail
