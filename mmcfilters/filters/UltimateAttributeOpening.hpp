#pragma once

#include "../utils/Image.hpp"
#include "../utils/Common.hpp"
#include "../trees/WeightedMorphologicalTree.hpp"
#include "../trees/WeightedTreeView.hpp"
#include "DepthStableRegionComputer.hpp"
#include "MSERComputer.hpp"

#include <cmath>
#include <concepts>
#include <memory>
#include <stack>
#include <stdexcept>
#include <string>
#include <vector>

namespace mmcfilters {

/**
 * @brief Computes an Ultimate Attribute Opening by accumulating maximal contrasts.
 *
 * @details
 * The object consumes an increasing attribute buffer indexed by dense internal
 * `NodeId` and computes, for each image point, the maximum altitude contrast and
 * the associated attribute index selected by the UAO traversal. The result is
 * materialised through image-producing accessors after `execute(...)` or
 * `executeWithMSER(...)`.
 *
 * Like other weighted-tree filter objects, it captures the tree mutation version
 * at construction time and rejects public operations after topology mutation.
 *
 * @tparam T Altitude type used by the weighted tree or weighted view.
 * @tparam Real Attribute-buffer floating-point type.
 */
template <AltitudeValue T, std::floating_point Real = float> class UltimateAttributeOpening {
  public:
    /// Floating-point type used for the input attribute buffer and criteria.
    using attribute_value_type = Real;

  protected:
    /// @cond INTERNAL
    /** @brief Defines the weighted altitude-view type used by the filter. */
    using AltitudeView = WeightedTreeView<T>;

    /** @brief Stores the largest criterion evaluated by the last execution. */
    Real maxCriterion = Real{0};
    /** @brief References the monotonically increasing attribute criteria. */
    const Real* const attrs_increasing;
    /** @brief Owns a copied attribute buffer when ownership is required. */
    std::shared_ptr<Real[]> ownedAttrsIncreasing_;
    /** @brief References the weighted-tree owner when one was supplied. */
    const WeightedMorphologicalTree<T>* weighted_ = nullptr;
    /** @brief Stores the non-owning weighted-tree view supplied at construction. */
    AltitudeView view_;
    /** @brief References the tree topology processed by the filter. */
    const MorphologicalTree& tree;
    /** @brief Stores the topology mutation version captured at construction. */
    std::size_t treeMutationVersion_ = 0;
    /** @brief Stores the maximum contrast selected for each tree node. */
    std::vector<T> maxContrastLUT;
    /** @brief Stores the selected attribute index for each tree node. */
    std::vector<int> associatedIndexLUT;
    /** @brief Marks the nodes selected for filtering by the current criterion. */
    std::vector<uint8_t> selectedForFiltering;

    /**
     * @brief Returns the active weighted-tree view.
     *
     * @return The active weighted-tree view.
     */
    AltitudeView view() const { return weighted_ != nullptr ? weighted_->asView() : view_; }

    /**
     * @brief Validates stable tree.
     *
     * @param context Operation name used in diagnostics.
     */
    void requireStableTree(const char* context) const { tree.requireMutationVersion(treeMutationVersion_, context); }

    /**
     * @brief Returns the altitude of a tree node.
     *
     * @param view Weighted-tree view used by the operation.
     * @param nodeId Identifier of the node used by the operation.
     * @return Altitude associated with the node.
     */
    static T altitudeOf(const AltitudeView& view, NodeId nodeId) { return view.getAltitude(nodeId); }

    /**
     * @brief Validates attribute pointer.
     *
     * @param attr Attribute requested by the operation.
     * @param context Operation name used in diagnostics.
     */
    static void requireAttributePointer(const Real* attr, const char* context) {
        if (attr == nullptr) {
            throw std::invalid_argument(std::string(context) + " requires a non-null attribute buffer.");
        }
    }

    /**
     * @brief Validates attribute buffer.
     *
     * @param tree Tree topology used by the operation.
     * @param attr Attribute requested by the operation.
     * @param context Operation name used in diagnostics.
     * @return Pointer to the resulting object.
     */
    static const Real* requireAttributeBuffer(const MorphologicalTree& tree, const std::vector<Real>& attr, const char* context) {
        if (attr.size() != static_cast<std::size_t>(tree.getNumInternalNodeSlots())) {
            throw std::invalid_argument(std::string(context) + " attribute size must match the internal node slot count.");
        }
        return attr.data();
    }

    /**
     * @brief Computes altitude difference.
     *
     * @param lhs Left-hand operand.
     * @param rhs Right-hand operand.
     * @return Computed altitude difference.
     */
    static T absoluteAltitudeDifference(AltitudeDiff<T> lhs, AltitudeDiff<T> rhs) {
        const long double difference = std::abs(static_cast<long double>(lhs) - static_cast<long double>(rhs));
        // UAO stores the contrast in the same type as the altitude by API
        // decision. For signed integral altitudes, contrasts that do not fit in
        // T can lose information here; avoiding that requires a separate
        // contrast type, which this API intentionally does not introduce.
        return static_cast<T>(difference);
    }

    /**
     * @brief Computes uao.
     *
     * @param view Weighted-tree view used by the operation.
     * @param currentNodeId Identifier represented by `currentNodeId`.
     * @param altitudeNodeNotInNR Altitude or level represented by `altitudeNodeNotInNR`.
     * @param qPropag Propagation queue used while updating the tree response.
     * @param isCalculateResidue Flag controlling is calculate residue.
     */
    void computeUAO(const AltitudeView& view, NodeId currentNodeId, AltitudeDiff<T> altitudeNodeNotInNR, bool qPropag, bool isCalculateResidue) {
        const NodeId parentNodeId = tree.getNodeParent(currentNodeId);
        const AltitudeDiff<T> altitudeNodeInNR = static_cast<AltitudeDiff<T>>(altitudeOf(view, currentNodeId));
        bool flagPropag = false;
        T contrast = T{};

        if (this->isSelectedForPruning(currentNodeId)) {
            altitudeNodeNotInNR = static_cast<AltitudeDiff<T>>(altitudeOf(view, parentNodeId));
            if (this->attrs_increasing[currentNodeId] <= this->maxCriterion) {
                isCalculateResidue = hasNodeSelectedInPrimitive(currentNodeId);
            }
        }

        if (this->attrs_increasing[currentNodeId] <= this->maxCriterion) {
            if (isCalculateResidue) {
                contrast = absoluteAltitudeDifference(altitudeNodeInNR, altitudeNodeNotInNR);
            }

            if (this->maxContrastLUT[parentNodeId] >= contrast) {
                this->maxContrastLUT[currentNodeId] = this->maxContrastLUT[parentNodeId];
                this->associatedIndexLUT[currentNodeId] = this->associatedIndexLUT[parentNodeId];
            } else {
                this->maxContrastLUT[currentNodeId] = contrast;
                this->associatedIndexLUT[currentNodeId] =
                    !qPropag ? static_cast<int>(this->attrs_increasing[currentNodeId] + 1) : this->associatedIndexLUT[parentNodeId];
                flagPropag = true;
            }
        }

        for (NodeId childNodeId : tree.getChildren(currentNodeId)) {
            this->computeUAO(view, childNodeId, altitudeNodeNotInNR, flagPropag, isCalculateResidue);
        }
    }

    /**
     * @brief Implements ultimate attribute opening for the selected node mask.
     *
     * @param maxCriterion Largest attribute criterion evaluated by the filter.
     * @param selectedForFiltering Per-node flags indicating which nodes are filtered.
     */
    void executeImpl(Real maxCriterion, const std::vector<uint8_t>& selectedForFiltering) {
        const AltitudeView altitudeView = view();
        this->maxCriterion = maxCriterion;
        this->selectedForFiltering = selectedForFiltering;

        for (NodeId id : tree.getAliveNodeIds()) {
            maxContrastLUT[id] = T{};
            associatedIndexLUT[id] = 0;
        }

        const NodeId rootNodeId = tree.getRoot();
        const AltitudeDiff<T> level = static_cast<AltitudeDiff<T>>(altitudeOf(altitudeView, rootNodeId));
        for (NodeId childNodeId : tree.getChildren(rootNodeId)) {
            computeUAO(altitudeView, childNodeId, level, false, false);
        }
    }

    /**
     * @brief Tests whether selected for pruning holds.
     *
     * @param currentNodeId Identifier represented by `currentNodeId`.
     * @return True when selected for pruning; otherwise false.
     */
    bool isSelectedForPruning(NodeId currentNodeId) const {
        const NodeId parentNodeId = tree.getNodeParent(currentNodeId);
        if (parentNodeId == InvalidNode) {
            return false;
        }
        return this->attrs_increasing[currentNodeId] != this->attrs_increasing[parentNodeId];
    }

    /**
     * @brief Tests whether node selected in primitive holds.
     *
     * @param currentNodeId Identifier represented by `currentNodeId`.
     * @return True when node selected in primitive; otherwise false.
     */
    bool hasNodeSelectedInPrimitive(NodeId currentNodeId) const {
        std::stack<NodeId> stack;
        stack.push(currentNodeId);
        while (!stack.empty()) {
            const NodeId nodeId = stack.top();
            stack.pop();
            if (selectedForFiltering[nodeId]) {
                return true;
            }

            for (NodeId childNodeId : tree.getChildren(nodeId)) {
                if (this->attrs_increasing[childNodeId] == this->attrs_increasing[nodeId]) {
                    stack.push(childNodeId);
                }
            }
        }
        return false;
    }
    /// @endcond

  public:
    /**
     * @brief Creates a UAO computation over a non-owning weighted view.
     *
     * @param view Weighted tree view whose topology and altitude define the
     * reconstruction domain.
     * @param attrs_increasing Shared increasing-attribute buffer indexed by
     * dense internal `NodeId`.
     */
    UltimateAttributeOpening(const AltitudeView& view, const std::shared_ptr<Real[]>& attrs_increasing)
        : UltimateAttributeOpening(view, attrs_increasing.get()) {
        this->ownedAttrsIncreasing_ = attrs_increasing;
    }

    /**
     * @brief Creates a UAO computation over a non-owning weighted view.
     *
     * @param view Weighted tree view whose topology and altitude define the
     * reconstruction domain.
     * @param attrs_increasing Increasing-attribute values indexed by dense
     * internal `NodeId`.
     * @throws std::invalid_argument If `attrs_increasing` does not match the
     * internal node slot count.
     */
    UltimateAttributeOpening(const AltitudeView& view, const std::vector<Real>& attrs_increasing)
        : UltimateAttributeOpening(view, requireAttributeBuffer(view.topology(), attrs_increasing, "UltimateAttributeOpening")) {}

    /**
     * @brief Creates a UAO computation over a non-owning weighted view.
     *
     * @param view Weighted tree view whose topology and altitude define the
     * reconstruction domain.
     * @param attrs_increasing Non-null increasing-attribute buffer indexed by
     * dense internal `NodeId`.
     * @throws std::invalid_argument If `attrs_increasing` is null or if the view
     * topology is stale.
     */
    UltimateAttributeOpening(const AltitudeView& view, const Real* attrs_increasing)
        : attrs_increasing(attrs_increasing), view_(view), tree(view_.topology()), treeMutationVersion_(tree.getMutationVersion()),
          maxContrastLUT(this->tree.getNumInternalNodeSlots()), associatedIndexLUT(this->tree.getNumInternalNodeSlots()) {
        view_.requireTopologyUnchanged("UltimateAttributeOpening");
        requireAttributePointer(attrs_increasing, "UltimateAttributeOpening");
        this->selectedForFiltering.assign(this->tree.getNumInternalNodeSlots(), true);
    }

    /**
     * @brief Creates a UAO computation over a borrowed weighted tree.
     *
     * The weighted tree is borrowed; it must outlive this object.
     *
     * @param weighted Weighted tree whose topology and altitude define the
     * reconstruction domain.
     * @param attrs_increasing Shared increasing-attribute buffer indexed by
     * dense internal `NodeId`.
     */
    UltimateAttributeOpening(const WeightedMorphologicalTree<T>& weighted, const std::shared_ptr<Real[]>& attrs_increasing)
        : UltimateAttributeOpening(weighted, attrs_increasing.get()) {
        this->ownedAttrsIncreasing_ = attrs_increasing;
    }

    /**
     * @brief Creates a UAO computation over a borrowed weighted tree.
     *
     * The weighted tree is borrowed; it must outlive this object.
     *
     * @param weighted Weighted tree whose topology and altitude define the
     * reconstruction domain.
     * @param attrs_increasing Increasing-attribute values indexed by dense
     * internal `NodeId`.
     * @throws std::invalid_argument If `attrs_increasing` does not match the
     * internal node slot count.
     */
    UltimateAttributeOpening(const WeightedMorphologicalTree<T>& weighted, const std::vector<Real>& attrs_increasing)
        : UltimateAttributeOpening(weighted.asView(), attrs_increasing) {
        weighted_ = &weighted;
    }

    /**
     * @brief Creates a UAO computation over a borrowed weighted tree.
     *
     * The weighted tree is borrowed; it must outlive this object.
     *
     * @param weighted Weighted tree whose topology and altitude define the
     * reconstruction domain.
     * @param attrs_increasing Non-null increasing-attribute buffer indexed by
     * dense internal `NodeId`.
     */
    UltimateAttributeOpening(const WeightedMorphologicalTree<T>& weighted, const Real* attrs_increasing)
        : UltimateAttributeOpening(weighted.asView(), attrs_increasing) {
        weighted_ = &weighted;
    }

    /**
     * @brief Destroys the ultimate-attribute-opening evaluator.
     */
    ~UltimateAttributeOpening() = default;

  public:
    /**
     * @brief Executes UAO using all internal tree nodes as selectable candidates.
     *
     * @param maxCriterion Maximum increasing-attribute threshold considered by
     * the UAO traversal.
     * @throws std::logic_error If the tree topology changed after construction.
     */
    void execute(Real maxCriterion) {
        requireStableTree("UltimateAttributeOpening::execute");
        std::vector<uint8_t> tmp(this->tree.getNumInternalNodeSlots(), true);
        executeImpl(maxCriterion, tmp);
    }

    /**
     * @brief Executes UAO with an explicit node-selection mask.
     *
     * @param maxCriterion Maximum increasing-attribute threshold considered by
     * the UAO traversal.
     * @param selectedForFiltering Dense internal-node mask marking selectable
     * primitive nodes.
     * @throws std::invalid_argument If `selectedForFiltering` does not match the
     * internal node slot count.
     * @throws std::logic_error If the tree topology changed after construction.
     */
    void execute(Real maxCriterion, const std::vector<uint8_t>& selectedForFiltering) {
        requireStableTree("UltimateAttributeOpening::execute");
        if (selectedForFiltering.size() != static_cast<std::size_t>(this->tree.getNumInternalNodeSlots())) {
            throw std::invalid_argument("UltimateAttributeOpening::execute selectedForFiltering size must match the internal node slot count.");
        }
        executeImpl(maxCriterion, selectedForFiltering);
    }

    /**
     * @brief Executes UAO with an MSER-derived node-selection mask.
     *
     * @param maxCriterion Maximum increasing-attribute threshold considered by
     * the UAO traversal.
     * @param deltaMSER Altitude delta used to compute the MSER stability mask.
     * @throws std::logic_error If this object was constructed from a view rather
     * than a weighted tree owner, or if the tree topology changed after
     * construction.
     */
    void executeWithMSER(Real maxCriterion, AltitudeDiff<T> deltaMSER) {
        requireStableTree("UltimateAttributeOpening::executeWithMSER");
        if (weighted_ == nullptr) {
            throw std::logic_error(
                "UltimateAttributeOpening::executeWithMSER requires a WeightedMorphologicalTree owner because MSER uses the tree-owned altitude.");
        }
        MSERComputer<T, Real> mser(*weighted_);
        executeImpl(maxCriterion, mser.computeMSER(deltaMSER));
    }

    /**
     * @brief Executes UAO with a depth-stability node-selection mask.
     *
     * @param maxCriterion Maximum increasing-attribute threshold considered by
     * the UAO traversal.
     * @param depthDelta Positive number of tree edges used to build the stability
     * window. Altitude is not read by the stability selection.
     * @throws std::invalid_argument If `depthDelta` is not positive.
     * @throws std::logic_error If the tree topology changed after construction.
     */
    void executeWithDepthStability(Real maxCriterion, int depthDelta) {
        requireStableTree("UltimateAttributeOpening::executeWithDepthStability");
        DepthStableRegionComputer<Real> stabilityComputer(this->tree);
        executeImpl(maxCriterion, stabilityComputer.computeByDepth(depthDelta));
    }

    /**
     * @brief Returns the per-pixel maximum UAO contrast image.
     *
     * @return Image on the original image domain using altitude type `T`.
     * @throws std::logic_error If the tree topology changed after construction.
     */
    [[nodiscard]] ImagePtr<T> getMaxContrastImage() const {
        requireStableTree("UltimateAttributeOpening::getMaxContrastImage");
        const int size = this->tree.getNumColsOfGridDomain2D() * this->tree.getNumRowsOfGridDomain2D();
        ImagePtr<T> imgOut = Image<T>::create(this->tree.getNumRowsOfGridDomain2D(), this->tree.getNumColsOfGridDomain2D());
        auto out = imgOut->rawData();

        for (int pidx = 0; pidx < size; pidx++) {
            out[pidx] = this->maxContrastLUT[tree.getProperPartOwner(pidx)];
        }
        return imgOut;
    }

    /**
     * @brief Returns the per-pixel associated attribute-index image.
     *
     * @return Signed integer image on the original image domain.
     * @throws std::logic_error If the tree topology changed after construction.
     */
    [[nodiscard]] ImageInt32Ptr getAssociatedImage() const {
        requireStableTree("UltimateAttributeOpening::getAssociatedImage");
        const int size = this->tree.getNumColsOfGridDomain2D() * this->tree.getNumRowsOfGridDomain2D();
        ImageInt32Ptr imgOut = ImageInt32::create(this->tree.getNumRowsOfGridDomain2D(), this->tree.getNumColsOfGridDomain2D());
        auto out = imgOut->rawData();

        for (int pidx = 0; pidx < size; pidx++) {
            out[pidx] = this->associatedIndexLUT[tree.getProperPartOwner(pidx)];
        }
        return imgOut;
    }

    /**
     * @brief Returns a color rendering of the associated-index image.
     *
     * @return RGB-like color image produced from `getAssociatedImage()`.
     * @throws std::logic_error If the tree topology changed after construction.
     */
    [[nodiscard]] ImageUInt8Ptr getAssociatedColorImage() const {
        requireStableTree("UltimateAttributeOpening::getAssociatedColorImage");
        return ImageUtils::createRandomColor(this->getAssociatedImage()->rawData(), this->tree.getNumRowsOfGridDomain2D(),
                                             this->tree.getNumColsOfGridDomain2D());
    }
};

} // namespace mmcfilters
