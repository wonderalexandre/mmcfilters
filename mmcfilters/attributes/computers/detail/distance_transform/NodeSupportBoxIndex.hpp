#pragma once

#include "../../../../contours/detail/NodeSupportIndex.hpp"

#include <algorithm>
#include <limits>
#include <vector>

namespace mmcfilters::attributes::computers::detail::distance_transform {

using ::mmcfilters::contours::detail::SupportInterval;

/**
 * @brief Smallest axis-aligned box containing a non-empty node support.
 */
struct SupportBox2D {
    int rowMin = std::numeric_limits<int>::max();
    int rowMax = -1;
    int columnMin = std::numeric_limits<int>::max();
    int columnMax = -1;

    /**
     * @brief Tests whether no support pixel has been accumulated.
     */
    [[nodiscard]] bool empty() const noexcept { return rowMax < rowMin || columnMax < columnMin; }

    /**
     * @brief Returns the box height.
     */
    [[nodiscard]] int rows() const noexcept { return empty() ? 0 : rowMax - rowMin + 1; }

    /**
     * @brief Returns the box width.
     */
    [[nodiscard]] int columns() const noexcept { return empty() ? 0 : columnMax - columnMin + 1; }

    /**
     * @brief Returns the number of pixels in the box.
     */
    [[nodiscard]] std::size_t area() const noexcept { return static_cast<std::size_t>(rows()) * static_cast<std::size_t>(columns()); }

    /**
     * @brief Expands the box to contain a row-major global pixel.
     */
    void include(PixelId pixel, int globalColumns) noexcept {
        const int row = pixel / globalColumns;
        const int column = pixel % globalColumns;
        rowMin = std::min(rowMin, row);
        rowMax = std::max(rowMax, row);
        columnMin = std::min(columnMin, column);
        columnMax = std::max(columnMax, column);
    }

    /**
     * @brief Expands the box to contain another support box.
     */
    void include(const SupportBox2D& other) noexcept {
        if (other.empty()) {
            return;
        }
        rowMin = std::min(rowMin, other.rowMin);
        rowMax = std::max(rowMax, other.rowMax);
        columnMin = std::min(columnMin, other.columnMin);
        columnMax = std::max(columnMax, other.columnMax);
    }

    /**
     * @brief Maps a contained global pixel to the translated local domain.
     */
    [[nodiscard]] PixelId localPixel(PixelId globalPixel, int globalColumns) const noexcept {
        const int globalRow = globalPixel / globalColumns;
        const int globalColumn = globalPixel % globalColumns;
        return (globalRow - rowMin) * columns() + (globalColumn - columnMin);
    }
};

/**
 * @brief Support index with bounding boxes for distance-transform domains.
 *
 * Pixel-contour consumers use the support index alone. Boxes are computed only
 * for distance transforms, which need translated rectangular work domains.
 */
class NodeSupportBoxIndex : public ::mmcfilters::contours::detail::NodeSupportIndex {
  public:
    /** @brief Builds the support intervals and their boxes in O(P + N). */
    explicit NodeSupportBoxIndex(const MorphologicalTree& tree)
        : NodeSupportIndex(tree), supportBoxByNode_(static_cast<std::size_t>(tree.numInternalNodeSlots())) {
        const int imageColumns = tree.requireGridDomain2D("NodeSupportBoxIndex").columns;
        for (NodeId node : tree.postOrder()) {
            auto& supportBox = supportBoxByNode_[static_cast<std::size_t>(node)];
            for (PixelId pixel : ::mmcfilters::detail::CommittedTreeAccess::properParts(tree, node)) {
                supportBox.include(pixel, imageColumns);
            }
            if (node != tree.root()) {
                supportBoxByNode_[static_cast<std::size_t>(::mmcfilters::detail::CommittedTreeAccess::nodeParent(tree, node))].include(supportBox);
            }
        }
    }

    /** @brief Returns the support box after checking the tree and node. */
    [[nodiscard]] const SupportBox2D& supportBox(NodeId node) const {
        (void)supportInterval(node);
        return establishedSupportBox(node);
    }

    /** @brief Returns the support box for a caller-established live node. */
    [[nodiscard]] const SupportBox2D& establishedSupportBox(NodeId node) const noexcept {
        return supportBoxByNode_[static_cast<std::size_t>(node)];
    }

  private:
    std::vector<SupportBox2D> supportBoxByNode_;
};

} // namespace mmcfilters::attributes::computers::detail::distance_transform
