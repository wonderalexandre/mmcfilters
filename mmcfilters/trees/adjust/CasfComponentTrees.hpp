#pragma once

#include "DualMinMaxTreeIncrementalFilter.hpp"

#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcfilters::adjust {

/**
 * @brief Increasing attribute used to select CASF pruning candidates.
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
 * `CasfComponentTrees` owns a mutable min-tree/max-tree pair built from the
 * same input image. For each threshold, it performs the standard two half-steps:
 *
 * - select maximal max-tree nodes whose increasing attribute is below the
 *   threshold, update the min-tree incrementally, then prune those max-tree
 *   subtrees;
 * - select maximal min-tree nodes under the same rule, update the max-tree
 *   incrementally, then prune those min-tree subtrees.
 *
 * The object is stateful: consecutive calls to `filter(...)` continue from the
 * current filtered trees. Attribute buffers are initialized once and then kept
 * synchronized through `DualMinMaxTreeIncrementalFilter`; they are not globally
 * recomputed at every threshold.
 */
class CasfComponentTrees {
private:
    AdjacencyRelation adjacency_;
    WeightedMorphologicalTree minTree_;
    WeightedMorphologicalTree maxTree_;
    std::unique_ptr<DynamicTreeAttributeComputer> minAttributeComputer_;
    std::unique_ptr<DynamicTreeAttributeComputer> maxAttributeComputer_;
    std::unique_ptr<DualMinMaxTreeIncrementalFilter<AltitudeType>> adjust_;
    std::vector<double> minAttributeBuffer_;
    std::vector<double> maxAttributeBuffer_;
    std::vector<NodeId> pruneCandidateQueue_;
    std::vector<NodeId> selectedPruneCandidates_;
    CasfComponentTreesAttribute attribute_ = CasfComponentTreesAttribute::AREA;

    static std::unique_ptr<DynamicTreeAttributeComputer> makeAttributeComputer(CasfComponentTreesAttribute attribute) {
        switch (attribute) {
            case CasfComponentTreesAttribute::AREA:
                return std::make_unique<DynamicAreaAttributeComputer>();
            case CasfComponentTreesAttribute::BOUNDING_BOX_WIDTH:
                return std::make_unique<DynamicBoundingBoxAttributeComputer>(BoundingBoxMeasure::WIDTH);
            case CasfComponentTreesAttribute::BOUNDING_BOX_HEIGHT:
                return std::make_unique<DynamicBoundingBoxAttributeComputer>(BoundingBoxMeasure::HEIGHT);
            case CasfComponentTreesAttribute::BOUNDING_BOX_DIAGONAL:
                return std::make_unique<DynamicBoundingBoxAttributeComputer>(BoundingBoxMeasure::DIAGONAL_LENGTH);
        }
        throw std::runtime_error("Unknown CASF component-tree attribute.");
    }

    /**
     * @brief Selects maximal non-root nodes satisfying the threshold.
     *
     * The traversal descends while the current node does not satisfy the
     * criterion. Once a non-root node is selected, its descendants are skipped,
     * preserving the maximal-candidate semantics used by the Higra CASF.
     */
    const std::vector<NodeId>& selectPruneCandidates(const WeightedMorphologicalTree& tree, const std::vector<double>& attribute, double threshold) {
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
     * @brief Applies one extensive and one anti-extensive half-step.
     */
    void applyFilterStep(double threshold) {
        const std::vector<NodeId> maxCandidates = selectPruneCandidates(maxTree_, maxAttributeBuffer_, threshold);
        adjust_->pruneMaxTreeAndUpdateMinTree(maxCandidates);

        const std::vector<NodeId> minCandidates = selectPruneCandidates(minTree_, minAttributeBuffer_, threshold);
        adjust_->pruneMinTreeAndUpdateMaxTree(minCandidates);
    }

public:
    /**
     * @brief Builds the initial min-tree/max-tree pair and bootstraps attributes.
     */
    CasfComponentTrees(ImageUInt8Ptr image, CasfComponentTreesAttribute attribute = CasfComponentTreesAttribute::AREA, double radius = 1.5)
        : adjacency_(image ? image->getNumRows() : 0, image ? image->getNumCols() : 0, radius),
          minTree_(WeightedMorphologicalTree::createComponentTree(image, false, radius)),
          maxTree_(WeightedMorphologicalTree::createComponentTree(image, true, radius)),
          minAttributeComputer_(makeAttributeComputer(attribute)),
          maxAttributeComputer_(makeAttributeComputer(attribute)),
          attribute_(attribute) {
        if (!image || image->getNumRows() <= 0 || image->getNumCols() <= 0 || image->getSize() <= 0) {
            throw std::invalid_argument("CasfComponentTrees requires a non-empty image.");
        }

        adjust_ = std::make_unique<DualMinMaxTreeIncrementalFilter<AltitudeType>>(&minTree_, &maxTree_, adjacency_);
        adjust_->setAttributeComputer(*minAttributeComputer_, *maxAttributeComputer_, minAttributeBuffer_, maxAttributeBuffer_);
        minAttributeComputer_->computeAttribute(minTree_, minAttributeBuffer_);
        maxAttributeComputer_->computeAttribute(maxTree_, maxAttributeBuffer_);

        const size_t maxNodes = static_cast<size_t>(
            std::max(minTree_.topology().getNumInternalNodeSlots(), maxTree_.topology().getNumInternalNodeSlots()));
        pruneCandidateQueue_.reserve(maxNodes);
        selectedPruneCandidates_.reserve(maxNodes);
    }

    /**
     * @brief Applies all thresholds to the current state and returns the reconstructed image.
     */
    ImageUInt8Ptr filter(const std::vector<double>& thresholds) {
        for (double threshold : thresholds) {
            applyFilterStep(threshold);
        }
        return minTree_.reconstructionImage();
    }

    const WeightedMorphologicalTree& minTree() const noexcept {
        return minTree_;
    }

    const WeightedMorphologicalTree& maxTree() const noexcept {
        return maxTree_;
    }

    CasfComponentTreesAttribute attribute() const noexcept {
        return attribute_;
    }

    /**
     * @brief Exports the current min-tree as a compact Higra-style parent/altitude pair.
     */
    std::pair<std::vector<NodeId>, std::vector<AltitudeType>> exportMinTree() const {
        return minTree_.exportHigraHierarchy();
    }

    /**
     * @brief Exports the current max-tree as a compact Higra-style parent/altitude pair.
     */
    std::pair<std::vector<NodeId>, std::vector<AltitudeType>> exportMaxTree() const {
        return maxTree_.exportHigraHierarchy();
    }
};

using ComponentTreeCasf = CasfComponentTrees;
using ComponentTreeCasfAttribute = CasfComponentTreesAttribute;

} // namespace mmcfilters::adjust
