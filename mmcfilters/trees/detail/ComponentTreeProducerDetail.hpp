#pragma once

#include "../MorphologicalTreeSemantics.hpp"
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
enum class ComponentTreePolarity { MaxTree, MinTree };

/**
 * @brief Concrete max-tree/min-tree producer over a regular 2D image.
 *
 * The producer preserves the historical union-find ordering and translates
 * its pixel-parent representation directly into the generic native hierarchy
 * buffers consumed by `MorphologicalTreeFactory`.
 */
template <AltitudeValue T> class ComponentTreeProducer {
    /** @brief Polarity. */
    ComponentTreePolarity polarity_;
    /** @brief Adjacency. */
    RegularGridAdjacency2D adjacency_;

  public:
    /**
     * @brief Creates a component-tree producer using an explicit adjacency.
     *
     * @param polarity Polarity that determines whether maxima or minima are processed.
     * @param adjacency Adjacency relation.
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
     * @param image Image.
     * @return The resulting proven owning topology, ownership, and altitude buffers.
     */
    [[nodiscard]] ValidatedNativeHierarchy<T> build(const ImagePtr<T>& image) {
        assert(polarity_ == ComponentTreePolarity::MaxTree || polarity_ == ComponentTreePolarity::MinTree);
        const bool isMaxTree = polarity_ == ComponentTreePolarity::MaxTree;
        const MorphologicalTreeKind kind = isMaxTree ? MorphologicalTreeKind::MaxTree : MorphologicalTreeKind::MinTree;

        const GridDomain2D gridDomain{image->getNumRows(), image->getNumColumns()};
        RegularGridAdjacency2D adjacency = std::move(adjacency_);
        assert(adjacency.getNumRows() == gridDomain.rows && adjacency.getNumColumns() == gridDomain.columns);
        ComponentTreeUnionFind unionFind(&adjacency, isMaxTree);
        auto [pixelParent, orderedPixels, numBuiltNodes] = unionFind.template build<T>(image);

        const int numPixels = image->getSize();
        const auto* imageValues = image->rawData();
        std::vector<NodeId> nodeParent;
        std::vector<T> altitude;
        std::vector<PixelId> representativePixel;
        nodeParent.reserve(static_cast<std::size_t>(numBuiltNodes));
        altitude.reserve(static_cast<std::size_t>(numBuiltNodes));
        representativePixel.reserve(static_cast<std::size_t>(numBuiltNodes));
        std::vector<NodeId> smallestNodeMap(static_cast<std::size_t>(numPixels), InvalidNode);

        TopologicalNativeHierarchyRecorder proofRecorder(static_cast<std::size_t>(numBuiltNodes), static_cast<std::size_t>(numPixels), NodeId{0});

        NodeId root = InvalidNode;
        for (int index = 0; index < numPixels; ++index) {
            const PixelId pixel = orderedPixels[static_cast<std::size_t>(index)];
            NodeId smallestNodeId = InvalidNode;

            if (pixel == pixelParent[static_cast<std::size_t>(pixel)]) {
                assert(root == InvalidNode);
                smallestNodeId = static_cast<NodeId>(nodeParent.size());
                root = smallestNodeId;
                nodeParent.push_back(smallestNodeId);
                altitude.push_back(imageValues[static_cast<std::size_t>(pixel)]);
                representativePixel.push_back(pixel);
                proofRecorder.recordSupportedNode(smallestNodeId, smallestNodeId);
            } else {
                const PixelId parentPixel = pixelParent[static_cast<std::size_t>(pixel)];
                const NodeId parentOwner = smallestNodeMap[static_cast<std::size_t>(parentPixel)];
                if (imageValues[static_cast<std::size_t>(pixel)] != imageValues[static_cast<std::size_t>(parentPixel)]) {
                    smallestNodeId = static_cast<NodeId>(nodeParent.size());
                    nodeParent.push_back(parentOwner);
                    altitude.push_back(imageValues[static_cast<std::size_t>(pixel)]);
                    representativePixel.push_back(pixel);
                    proofRecorder.recordSupportedNode(smallestNodeId, parentOwner);
                } else {
                    smallestNodeId = parentOwner;
                }
            }

            smallestNodeMap[static_cast<std::size_t>(pixel)] = smallestNodeId;
            if (pixel < representativePixel[static_cast<std::size_t>(smallestNodeId)]) {
                representativePixel[static_cast<std::size_t>(smallestNodeId)] = pixel;
                altitude[static_cast<std::size_t>(smallestNodeId)] = imageValues[static_cast<std::size_t>(pixel)];
            }
            proofRecorder.recordSmallestNode(smallestNodeId);
        }

        assert(root == 0 && nodeParent.size() == static_cast<std::size_t>(numBuiltNodes));

        return makeValidatedNativeHierarchy<T>(std::move(nodeParent), std::move(smallestNodeMap), std::move(altitude), root, gridDomain,
                                               makeMorphologicalTreeSemantics(kind, SharedAdjacencyContext{std::move(adjacency)}),
                                               std::move(proofRecorder).finish());
    }
};

} // namespace kernel

} // namespace mmcfilters::detail
