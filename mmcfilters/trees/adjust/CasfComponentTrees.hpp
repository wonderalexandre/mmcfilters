#pragma once

#include "DualMinMaxTreeIncrementalFilter.hpp"
#include "../MorphologicalTreeFactory.hpp"

#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcfilters::adjust {

/**
 * @brief Increasing attribute used to select CASF pruning candidates.
 *
 * The CASF only requires the selected attribute to be increasing on the
 * component-tree order. Area is available for every image domain. Bounding-box
 * attributes are computed from the image grid coordinates and expose three
 * scalar projections of the same cached box state.
 */
enum class CasfComponentTreesAttribute {
    AREA,
    BOUNDING_BOX_WIDTH,
    BOUNDING_BOX_HEIGHT,
    BOUNDING_BOX_DIAGONAL
};

/**
 * @brief Connected alternating sequential filter on paired component trees.
 *
 * This class implements a connected alternating sequential filter (CASF) on a
 * fixed image-domain adjacency. Instead of rebuilding both component trees
 * after every threshold, it owns a mutable min-tree / max-tree pair and updates
 * one tree incrementally whenever rooted subtrees are pruned from the dual one.
 *
 * Supported workflow:
 *
 * - build the initial min-tree and max-tree from the same input image;
 * - choose one increasing attribute used to select pruning candidates
 *   (`AREA`, `BOUNDING_BOX_WIDTH`, `BOUNDING_BOX_HEIGHT`, or
 *   `BOUNDING_BOX_DIAGONAL`);
 * - run `filter(thresholds)` on one or more threshold sequences applied
 *   successively to the current state;
 * - inspect or export the current dynamic tree state through `minTree`,
 *   `maxTree`, `exportMinTree`, and `exportMaxTree`.
 *
 * For each threshold, the standard two CASF half-steps are performed:
 *
 * - select maximal max-tree nodes whose increasing attribute is below the
 *   threshold, update the min-tree incrementally, then prune those max-tree
 *   subtrees;
 * - select maximal min-tree nodes under the same rule, update the max-tree
 *   incrementally, then prune those min-tree subtrees.
 *
 * Internal state owned by this class:
 *
 * - the adjacency relation of the image domain;
 * - the mutable weighted min-tree and max-tree;
 * - one incremental attribute computer per tree and their output buffers;
 * - the dual-tree adjustment helper used to propagate pruning from one tree to
 *   the other;
 * - temporary buffers reused for breadth-first maximal candidate selection.
 *
 * Update semantics:
 *
 * - calling `filter(...)` with a non-empty threshold sequence mutates the
 *   internal trees;
 * - later `filter(...)` calls continue from the current filtered state;
 * - an empty threshold sequence leaves the state unchanged and simply
 *   reconstructs the current image from the min-tree.
 *
 * Attribute buffers are initialized once and then kept synchronized through
 * `DualMinMaxTreeIncrementalFilter`; they are not globally recomputed at every
 * threshold.
 *
 * Ownership note:
 *
 * CASF intentionally remains a mutable `WeightedMorphologicalTree<T>` owner.
 * `WeightedTreeView<T>` is suited to read-only kernels; here each threshold step
 * changes tree topology, proper-part ownership, altitude state, and dynamic
 * attribute buffers together.
 */
template<AltitudeValue T>
class CasfComponentTrees {
private:
    using tree_t = WeightedMorphologicalTree<T>;
    using image_ptr_t = ImagePtr<T>;
    using attribute_computer_t = DynamicTreeAttributeComputer<T>;

    AdjacencyRelation adjacency_;
    tree_t minTree_;
    tree_t maxTree_;
    std::unique_ptr<attribute_computer_t> minAttributeComputer_;
    std::unique_ptr<attribute_computer_t> maxAttributeComputer_;
    std::unique_ptr<DualMinMaxTreeIncrementalFilter<T>> adjust_;
    std::vector<double> minAttributeBuffer_;
    std::vector<double> maxAttributeBuffer_;
    std::vector<NodeId> pruneCandidateQueue_;
    std::vector<NodeId> selectedPruneCandidates_;
    CasfComponentTreesAttribute attribute_ = CasfComponentTreesAttribute::AREA;

    /**
     * @brief Creates the incremental attribute computer matching the selected CASF attribute.
     * @details Area uses a scalar area computer. Bounding-box attributes reuse
     * the row-major image embedding already stored by `WeightedMorphologicalTree<T>`.
     */
    static std::unique_ptr<attribute_computer_t> makeAttributeComputer(CasfComponentTreesAttribute attribute) {
        switch (attribute) {
            case CasfComponentTreesAttribute::AREA:
                return std::make_unique<DynamicAreaAttributeComputer<T>>();
            case CasfComponentTreesAttribute::BOUNDING_BOX_WIDTH:
                return std::make_unique<DynamicBoundingBoxAttributeComputer<T>>(BoundingBoxMeasure::WIDTH);
            case CasfComponentTreesAttribute::BOUNDING_BOX_HEIGHT:
                return std::make_unique<DynamicBoundingBoxAttributeComputer<T>>(BoundingBoxMeasure::HEIGHT);
            case CasfComponentTreesAttribute::BOUNDING_BOX_DIAGONAL:
                return std::make_unique<DynamicBoundingBoxAttributeComputer<T>>(BoundingBoxMeasure::DIAGONAL_LENGTH);
        }
        throw std::runtime_error("Unknown CASF component-tree attribute.");
    }

    /**
     * @brief Selects maximal non-root pruning candidates whose attribute is below or equal to a threshold.
     *
     * The traversal follows the same semantics used by the Higra CASF: descend
     * while `attribute(node) > threshold`; otherwise select the current non-root
     * node and stop descending below it. The root is never returned because
     * pruning the root would not define a valid half-step. Internal queue/output
     * buffers are reused across calls to avoid allocation in the threshold loop.
     */
    const std::vector<NodeId>& selectPruneCandidates(const tree_t& tree, const std::vector<double>& attribute, double threshold) {
        pruneCandidateQueue_.clear();
        selectedPruneCandidates_.clear();

        const MorphologicalTree& topology = tree.topology();
        const NodeId root = topology.getRoot();
        if (root == InvalidNode || !topology.isAlive(root)) {
            return selectedPruneCandidates_;
        }

        pruneCandidateQueue_.push_back(root);
        size_t head = 0;

        while (head < pruneCandidateQueue_.size()) {
            const NodeId nodeId = pruneCandidateQueue_[head++];
            if (!topology.isAlive(nodeId)) {
                continue;
            }

            const double nodeAttribute = attribute[static_cast<size_t>(nodeId)];
            if (!topology.isRoot(nodeId) && nodeAttribute <= threshold) {
                selectedPruneCandidates_.push_back(nodeId);
                continue;
            }

            for (NodeId childId : topology.getChildren(nodeId)) {
                if (topology.isAlive(childId)) {
                    pruneCandidateQueue_.push_back(childId);
                }
            }
        }

        return selectedPruneCandidates_;
    }

    /**
     * @brief Applies one threshold step of the CASF on the current dynamic state.
     *
     * The step first applies the extensive half-step (max-tree pruning with
     * min-tree update) and then the anti-extensive half-step (min-tree pruning
     * with max-tree update). This order matches the original component-tree CASF
     * schedule.
     */
    void applyFilterStep(double threshold) {
        const std::vector<NodeId> maxCandidates = selectPruneCandidates(maxTree_, maxAttributeBuffer_, threshold);
        adjust_->pruneMaxTreeAndUpdateMinTree(maxCandidates);

        const std::vector<NodeId> minCandidates = selectPruneCandidates(minTree_, minAttributeBuffer_, threshold);
        adjust_->pruneMinTreeAndUpdateMaxTree(minCandidates);
    }

public:
    /**
     * @brief Initializes the CASF state from the input image and the chosen attribute.
     * @param image Input gray-level image.
     * @param attribute Increasing attribute used to select pruning candidates.
     * @param radius Radius used to build the image-domain adjacency relation.
     */
    CasfComponentTrees(image_ptr_t image, CasfComponentTreesAttribute attribute = CasfComponentTreesAttribute::AREA, double radius = 1.5)
        : adjacency_(image ? image->getNumRows() : 0, image ? image->getNumCols() : 0, radius),
          minTree_(MorphologicalTreeFactory::createMinTree(image, radius)),
          maxTree_(MorphologicalTreeFactory::createMaxTree(image, radius)),
          minAttributeComputer_(makeAttributeComputer(attribute)),
          maxAttributeComputer_(makeAttributeComputer(attribute)),
          attribute_(attribute) {
        if (!image || image->getNumRows() <= 0 || image->getNumCols() <= 0 || image->getSize() <= 0) {
            throw std::invalid_argument("CasfComponentTrees requires a non-empty image.");
        }

        adjust_ = std::make_unique<DualMinMaxTreeIncrementalFilter<T>>(&minTree_, &maxTree_, adjacency_);
        adjust_->setAttributeComputer(*minAttributeComputer_, *maxAttributeComputer_, minAttributeBuffer_, maxAttributeBuffer_);
        minAttributeComputer_->computeAttribute(minTree_, minAttributeBuffer_);
        maxAttributeComputer_->computeAttribute(maxTree_, maxAttributeBuffer_);

        const size_t maxNodes = static_cast<size_t>(
            std::max(minTree_.topology().getNumInternalNodeSlots(), maxTree_.topology().getNumInternalNodeSlots()));
        pruneCandidateQueue_.reserve(maxNodes);
        selectedPruneCandidates_.reserve(maxNodes);
    }

    /**
     * @brief Runs the CASF on the threshold sequence and returns the filtered image.
     *
     * Thresholds are applied to the current internal state. The object is
     * therefore stateful: successive non-empty calls continue filtering the
     * result of the previous ones.
     */
    [[nodiscard]] image_ptr_t filter(const std::vector<double>& thresholds) {
        for (double threshold : thresholds) {
            applyFilterStep(threshold);
        }
        return minTree_.reconstructionImage();
    }

    /**
     * @brief Returns the current min-tree state.
     */
    [[nodiscard]] const tree_t& minTree() const noexcept {
        return minTree_;
    }

    /**
     * @brief Returns the current max-tree state.
     */
    [[nodiscard]] const tree_t& maxTree() const noexcept {
        return maxTree_;
    }

    /**
     * @brief Returns the increasing attribute configured for this CASF instance.
     */
    [[nodiscard]] CasfComponentTreesAttribute attribute() const noexcept {
        return attribute_;
    }

    /**
     * @brief Exports the current min-tree as a compact static parent/altitude pair.
     * @details The pair follows the local `WeightedMorphologicalTree<T>`
     * export convention and can be compared with Higra-style static hierarchy
     * outputs in tests and benchmarks.
     */
    [[nodiscard]] std::pair<std::vector<NodeId>, std::vector<T>> exportMinTree() const {
        return minTree_.exportHigraHierarchy();
    }

    /**
     * @brief Exports the current max-tree as a compact static parent/altitude pair.
     * @details The pair follows the local `WeightedMorphologicalTree<T>`
     * export convention and can be compared with Higra-style static hierarchy
     * outputs in tests and benchmarks.
     */
    [[nodiscard]] std::pair<std::vector<NodeId>, std::vector<T>> exportMaxTree() const {
        return maxTree_.exportHigraHierarchy();
    }
};

} // namespace mmcfilters::adjust
