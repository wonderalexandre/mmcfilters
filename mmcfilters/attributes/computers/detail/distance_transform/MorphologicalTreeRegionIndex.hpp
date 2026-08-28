#pragma once

#include "../../../../trees/MorphologicalTree.hpp"
#include "../../../../trees/detail/CommittedTreeAccess.hpp"
#include "../../../../utils/Common.hpp"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace mmcfilters::attributes::computers::detail::distance_transform {

/**
 * @brief Inclusive axis-aligned bounds of one non-empty node support.
 */
struct RegionBox2D {
    int rowMin = std::numeric_limits<int>::max();
    int rowMax = -1;
    int columnMin = std::numeric_limits<int>::max();
    int columnMax = -1;

    /**
     * @brief Tests whether no support pixel has been accumulated.
     */
    [[nodiscard]] bool empty() const noexcept { return rowMax < rowMin || columnMax < columnMin; }

    /**
     * @brief Returns the active-domain height.
     */
    [[nodiscard]] int rows() const noexcept { return empty() ? 0 : rowMax - rowMin + 1; }

    /**
     * @brief Returns the active-domain width.
     */
    [[nodiscard]] int columns() const noexcept { return empty() ? 0 : columnMax - columnMin + 1; }

    /**
     * @brief Returns the number of cells in the active rectangular domain.
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
    void include(const RegionBox2D& other) noexcept {
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
 * @brief Half-open interval in the shared pre-order support storage.
 */
struct RegionSupportInterval {
    std::size_t begin = 0;
    std::size_t end = 0;

    /**
     * @brief Returns the support cardinality represented by the interval.
     */
    [[nodiscard]] std::size_t size() const noexcept { return end - begin; }
};

/**
 * @brief Dense immutable support index for one committed morphological tree.
 *
 * Proper parts are appended once in pre-order. Consequently, all proper parts
 * in a node subtree occupy one contiguous interval. Support cardinalities give
 * the interval ends, while support boxes are accumulated bottom-up. Dead node
 * slots retain empty intervals and boxes.
 *
 * The referenced tree must outlive the index. Returned spans borrow the
 * index-owned pixel storage and remain valid only while the index exists and
 * the tree stays at the captured mutation version.
 */
class MorphologicalTreeRegionIndex {
  public:
    MorphologicalTreeRegionIndex(const MorphologicalTreeRegionIndex&) = delete;
    MorphologicalTreeRegionIndex& operator=(const MorphologicalTreeRegionIndex&) = delete;
    MorphologicalTreeRegionIndex(MorphologicalTreeRegionIndex&&) = delete;
    MorphologicalTreeRegionIndex& operator=(MorphologicalTreeRegionIndex&&) = delete;

    /**
     * @brief Builds the O(P + N) support index for a committed 2D tree.
     */
    explicit MorphologicalTreeRegionIndex(const MorphologicalTree& tree)
        : tree_(tree), mutationVersion_(tree.getMutationVersion()), domain_(requireDomain(tree)),
          intervals_(static_cast<std::size_t>(tree.numInternalNodeSlots())), boxes_(static_cast<std::size_t>(tree.numInternalNodeSlots())) {
        supportPixels_.reserve(static_cast<std::size_t>(tree_.numPixels()));
        std::vector<NodeId> preOrderNodes;
        preOrderNodes.reserve(static_cast<std::size_t>(tree_.numNodes()));
        std::vector<std::size_t> supportCardinalities(static_cast<std::size_t>(tree_.numInternalNodeSlots()), std::size_t{0});

        for (NodeId node : ::mmcfilters::detail::CommittedTreeAccess::subtree(tree_, tree_.root())) {
            const std::size_t nodeIndex = static_cast<std::size_t>(node);
            preOrderNodes.push_back(node);
            intervals_[nodeIndex].begin = supportPixels_.size();

            for (PixelId pixel : ::mmcfilters::detail::CommittedTreeAccess::properParts(tree_, node)) {
                supportPixels_.push_back(pixel);
                ++supportCardinalities[nodeIndex];
                boxes_[nodeIndex].include(pixel, domain_.columns);
            }
        }

        for (auto iterator = preOrderNodes.rbegin(); iterator != preOrderNodes.rend(); ++iterator) {
            const NodeId node = *iterator;
            const std::size_t nodeIndex = static_cast<std::size_t>(node);
            const std::size_t supportCardinality = supportCardinalities[nodeIndex];
            intervals_[nodeIndex].end = intervals_[nodeIndex].begin + supportCardinality;
            if (node != tree_.root()) {
                const NodeId parent = ::mmcfilters::detail::CommittedTreeAccess::nodeParent(tree_, node);
                const std::size_t parentIndex = static_cast<std::size_t>(parent);
                supportCardinalities[parentIndex] += supportCardinality;
                boxes_[parentIndex].include(boxes_[nodeIndex]);
            }
        }
    }

    /**
     * @brief Returns the contiguous support pixels of one live node.
     */
    [[nodiscard]] std::span<const PixelId> support(NodeId node) const {
        const RegionSupportInterval interval = supportInterval(node);
        return std::span<const PixelId>(supportPixels_).subspan(interval.begin, interval.size());
    }

    /**
     * @brief Returns the half-open support interval of one live node.
     */
    [[nodiscard]] RegionSupportInterval supportInterval(NodeId node) const {
        requireLiveNode(node);
        return intervals_[static_cast<std::size_t>(node)];
    }

    /**
     * @brief Returns the minimal non-empty support box of one live node.
     */
    [[nodiscard]] const RegionBox2D& boundingBox(NodeId node) const {
        requireLiveNode(node);
        return boxes_[static_cast<std::size_t>(node)];
    }

    /**
     * @brief Returns support for a caller-established live node.
     */
    [[nodiscard]] std::span<const PixelId> establishedSupport(NodeId node) const noexcept {
        const RegionSupportInterval interval = establishedSupportInterval(node);
        return std::span<const PixelId>(supportPixels_).subspan(interval.begin, interval.size());
    }

    /**
     * @brief Returns the support interval of a caller-established live node.
     */
    [[nodiscard]] RegionSupportInterval establishedSupportInterval(NodeId node) const noexcept {
        return intervals_[static_cast<std::size_t>(node)];
    }

    /**
     * @brief Returns the support box of a caller-established live node.
     */
    [[nodiscard]] const RegionBox2D& establishedBoundingBox(NodeId node) const noexcept {
        return boxes_[static_cast<std::size_t>(node)];
    }

    /**
     * @brief Returns the number of proper-part pixels stored exactly once.
     */
    [[nodiscard]] std::size_t numIndexedPixels() const {
        tree_.requireMutationVersion(mutationVersion_, "MorphologicalTreeRegionIndex");
        return supportPixels_.size();
    }

    /**
     * @brief Rejects use with a different tree, mutation version, or domain.
     */
    void requireCompatibleTree(const MorphologicalTree& tree) const {
        if (&tree != &tree_) {
            throw std::invalid_argument("Morphological-tree region index belongs to a different tree instance.");
        }
        if (tree_.getMutationVersion() != mutationVersion_) {
            throw std::logic_error("Morphological-tree region index belongs to an older tree mutation version.");
        }
        const GridDomain2D domain = tree_.requireGridDomain2D("MorphologicalTreeRegionIndex");
        if (domain.rows != domain_.rows || domain.columns != domain_.columns) {
            throw std::logic_error("Morphological-tree region index no longer matches the captured 2D domain.");
        }
    }

  private:
    [[nodiscard]] static GridDomain2D requireDomain(const MorphologicalTree& tree) {
        tree.requireNotEditing("MorphologicalTreeRegionIndex");
        if (!tree.hasGridDomain2D()) {
            throw std::invalid_argument("Morphological-tree region indexing requires a regular 2D pixel domain.");
        }
        const GridDomain2D domain = *tree.gridDomain2D();
        if (domain.rows <= 0 || domain.columns <= 0 || domain.columns > std::numeric_limits<int>::max() / domain.rows ||
            domain.rows * domain.columns != tree.numPixels()) {
            throw std::invalid_argument("Morphological-tree region indexing requires a consistent non-empty 2D pixel domain.");
        }
        return domain;
    }

    void requireLiveNode(NodeId node) const {
        tree_.requireMutationVersion(mutationVersion_, "MorphologicalTreeRegionIndex");
        if (!tree_.isNode(node) || !::mmcfilters::detail::CommittedTreeAccess::isAlive(tree_, node)) {
            throw std::out_of_range("Morphological-tree region index received a non-live node id.");
        }
    }

    const MorphologicalTree& tree_;
    std::size_t mutationVersion_ = 0;
    GridDomain2D domain_{};
    std::vector<PixelId> supportPixels_;
    std::vector<RegionSupportInterval> intervals_;
    std::vector<RegionBox2D> boxes_;
};

} // namespace mmcfilters::attributes::computers::detail::distance_transform
