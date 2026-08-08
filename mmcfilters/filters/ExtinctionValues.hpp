#pragma once

#include "../trees/TreeAltitudeAlgorithms.hpp"
#include "../trees/WeightedMorphologicalTree.hpp"
#include "../trees/WeightedTreeView.hpp"
#include "../trees/detail/HierarchyCapabilityValidation.hpp"
#include "../trees/saliency/HierarchySaliencyMapValidation.hpp"
#include "../trees/saliency/HierarchySaliencyMap.hpp"
#include "../trees/saliency/HierarchicalWatershedSaliency.hpp"
#include "../utils/Image.hpp"
#include "../utils/Common.hpp"
#include "../contours/ContoursComputedIncrementally.hpp"

#include <algorithm>
#include <cmath>
#include <concepts>
#include <limits>
#include <memory>
#include <span>
#include <stack>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mmcfilters {

/**
 * @brief Record describing one regional extremum and its extinction value.
 */
template <std::floating_point Real = float> struct RegionalExtremaNode {
    /// Leaf node that represents the regional extremum in a max-tree/min-tree.
    NodeId leaf;

    /// Highest node retained before the extremum merges with a stronger branch.
    NodeId cutoffNode;

    /// Attribute value at `cutoffNode`, or `numeric_limits<Real>::max()` for
    /// the dominant extremum that survives until the root.
    Real extinction;

    /**
     * @brief Builds one extinction-value record.
     *
     * @param leaf Leaf node representing the regional extremum.
     * @param cutoffNode Node where the extremum stops being dominant.
     * @param extinction Attribute value associated with the cutoff.
     */
    RegionalExtremaNode(NodeId leaf, NodeId cutoffNode, Real extinction) : leaf(leaf), cutoffNode(cutoffNode), extinction(extinction) {}
};

/**
 * @brief Explicit extinction-extrema selection policy.
 *
 * The same policy object is consumed by image reconstruction (`filtering`) and
 * by contour visualization (`contourMap`). It only controls which regional
 * extrema are retained; output-specific scoring is configured separately.
 */
template <std::floating_point Real = float> class ExtinctionSelectionPolicy {
    /** @brief Enumerates the supported mode values. */
    enum class Mode { TopK, MinimumExtinction };

    /** @brief Stores the mode. */
    Mode mode_ = Mode::TopK;
    /** @brief Stores the extrema to keep. */
    int extremaToKeep_ = 0;
    /** @brief Stores the threshold. */
    Real threshold_ = Real{0};

  public:
    /**
     * @brief Select the strongest extrema by decreasing extinction ranking.
     *
     * @param extremaToKeep Maximum number of regional extrema retained.
     * Validation is deferred to the operation that consumes the policy.
     * @return Selection policy configured for top-k extinction filtering.
     */
    [[nodiscard]] static ExtinctionSelectionPolicy byTopK(int extremaToKeep) noexcept {
        ExtinctionSelectionPolicy policy;
        policy.mode_ = Mode::TopK;
        policy.extremaToKeep_ = extremaToKeep;
        return policy;
    }

    /**
     * @brief Select every extremum whose extinction value is at least `threshold`.
     *
     * @param threshold Minimum extinction value retained. Validation is deferred
     * to the operation that consumes the policy.
     * @return Selection policy configured for threshold extinction filtering.
     */
    [[nodiscard]] static ExtinctionSelectionPolicy byThreshold(Real threshold) noexcept {
        ExtinctionSelectionPolicy policy;
        policy.mode_ = Mode::MinimumExtinction;
        policy.threshold_ = threshold;
        return policy;
    }

    /**
     * @brief Returns true when the policy selects the first ranked extrema.
     *
     * @return True for top-k policies; false for threshold policies.
     */
    [[nodiscard]] bool selectsTopK() const noexcept { return mode_ == Mode::TopK; }

    /**
     * @brief Maximum number of ranked extrema retained by a top-k policy.
     *
     * @return Stored top-k count. The value is meaningful only when
     * `selectsTopK()` returns true.
     */
    [[nodiscard]] int extremaToKeep() const noexcept { return extremaToKeep_; }

    /**
     * @brief Minimum accepted extinction value for a threshold policy.
     *
     * @return Stored extinction threshold. The value is meaningful only when
     * `selectsTopK()` returns false.
     */
    [[nodiscard]] Real minimumExtinction() const noexcept { return threshold_; }
};

/**
 * @brief Score convention used when drawing extinction contour maps.
 */
enum class ExtinctionContourScorePolicy {
    /// Use dense rank scores over the selected extrema.
    RankScore,
    /// Use the raw extinction value associated with each selected extremum.
    ExtinctionValue
};

/**
 * @brief Computes and stores extinction values for regional extrema.
 *
 * @details
 * `ExtinctionValues` implements the classical leaf-extrema extinction ranking
 * for hierarchies that declare a globally monotone altitude order. Standard
 * max-tree and min-tree producers provide this capability. In this setting,
 * the regional extrema processed by the algorithm are the tree leaves. The supplied scalar
 * attribute is indexed by dense internal `NodeId`, must have one value for every
 * internal node slot, and is interpreted so that larger values represent
 * stronger extrema. Results are sorted in decreasing extinction order and can be
 * consumed either as records, a filtered reconstruction, or a contour saliency
 * map.
 *
 * The strongest extremum has no stronger merge point. Its extinction value is
 * represented by the explicit finite sentinel `numeric_limits<Real>::max()`.
 *
 * Standard tree-of-shapes and self-dual residual-tree producers declare
 * `AltitudeOrder::UNCONSTRAINED` and are rejected because their complete
 * regional-extrema set is not generally equivalent to `tree.getLeaves()`.
 * Acceptance is based on the altitude-order capability, not the descriptive
 * tree kind.
 *
 * The object records the tree mutation version at construction time. Public
 * operations reject use after the underlying topology changes.
 *
 * @par Primary algorithmic reference
 * Alexandre Gonçalves Silva and Roberto de Alencar Lotufo, "Efficient
 * computation of new extinction values from extended component tree," Pattern
 * Recognition Letters, 32(1):79-90, 2011.
 * [DOI 10.1016/j.patrec.2010.07.019](https://doi.org/10.1016/j.patrec.2010.07.019).
 * The component-tree branch traversal follows Algorithm 1; this class receives
 * an already constructed tree and attribute buffer instead of incrementally
 * computing them during tree construction.
 *
 * @tparam T Altitude type used by the weighted tree or weighted view.
 * @tparam Real Attribute-buffer floating-point type.
 */
template <AltitudeValue T, std::floating_point Real = float> class ExtinctionValues {
  protected:
    /// @cond INTERNAL
    /** @brief Defines the weighted altitude-view type used by the computation. */
    using AltitudeView = WeightedTreeView<T>;

    /** @brief Stores the regional extrema ordered by the extinction computation. */
    std::vector<RegionalExtremaNode<Real>> regionalExtremaNodes_;
    /** @brief Stores the extinction value indexed by tree node. */
    std::vector<Real> extinctionValueAttribute_;
    /** @brief Stores the non-owning weighted-tree view supplied at construction. */
    AltitudeView view_;
    /** @brief References the weighted-tree owner when one was supplied. */
    const WeightedMorphologicalTree<T>* weighted_ = nullptr;
    /** @brief References the tree topology processed by the computation. */
    const MorphologicalTree& tree;
    /** @brief Stores the topology mutation version captured at construction. */
    std::size_t treeMutationVersion_ = 0;

    /** @brief Associates one regional-extremum record with its deterministic rank. */
    struct SelectedExtremum {
        /** @brief References the selected regional-extremum record. */
        const RegionalExtremaNode<Real>* record = nullptr;
        /** @brief Stores the ranking score assigned to the selected extremum. */
        int rankScore = 0;
    };

    /**
     * @brief Returns the freshest weighted view used for reconstruction.
     *
     * Instances built from a `WeightedMorphologicalTree` keep a pointer to that
     * owner so later altitude-buffer replacements are observed. Instances built
     * from an external `WeightedTreeView` keep the view snapshot supplied at
     * construction time.
     *
     * @return Weighted view over the same topology captured by this object.
     */
    [[nodiscard]] AltitudeView view() const { return weighted_ != nullptr ? weighted_->asView() : view_; }

    /**
     * @brief Rejects operations after the borrowed tree topology changes.
     *
     * Extinction records store dense `NodeId` indexes and are valid only for the
     * topology version captured by the constructor.
     *
     * @param context Name of the public/internal operation reported in errors.
     * @throws std::logic_error If the tree mutation version no longer matches.
     */
    void requireStableTree(const char* context) const { tree.requireMutationVersion(treeMutationVersion_, context); }

    /**
     * @brief Validates raw attribute pointer arguments.
     *
     * @param attr Attribute buffer passed to an extinction constructor.
     * @param context Name used in the exception message.
     * @throws std::invalid_argument If `attr` is null.
     */
    static void requireAttributePointer(const Real* attr, const char* context) {
        if (attr == nullptr) {
            throw std::invalid_argument(std::string(context) + " requires a non-null attribute buffer.");
        }
    }

    /**
     * @brief Validates vector attribute buffers and returns their raw pointer.
     *
     * @param tree Topology that defines the required dense internal-node size.
     * @param attr Attribute vector indexed by dense internal `NodeId`.
     * @param context Name used in the exception message.
     * @return Non-null pointer to the vector storage.
     * @throws std::invalid_argument If `attr.size()` is not the internal node
     * slot count of `tree`.
     */
    static const Real* requireAttributeBuffer(const MorphologicalTree& tree, const std::vector<Real>& attr, const char* context) {
        if (attr.size() != static_cast<std::size_t>(tree.getNumInternalNodeSlots())) {
            throw std::invalid_argument(std::string(context) + " attribute size must match the internal node slot count.");
        }
        return attr.data();
    }

    /**
     * @brief Reads one node altitude through the active weighted view.
     *
     * This wrapper keeps reconstruction code independent from whether the object
     * was built from a weighted owner or from an external view.
     *
     * @param view Weighted view containing the altitude buffer.
     * @param nodeId Dense internal node id to read.
     * @return Altitude value associated with `nodeId`.
     */
    static T altitudeOf(const AltitudeView& view, NodeId nodeId) { return view.getAltitude(nodeId); }

    /**
     * @brief Collects regional extrema candidates for component-tree extinction.
     *
     * For a hierarchy with globally monotone altitude order, regional extrema
     * are represented by tree leaves. Standard tree-of-shapes and self-dual
     * residual-tree producers are rejected before this helper is used because
     * they declare unconstrained altitude order and their extrema are not
     * generally equivalent to `getLeaves()`.
     *
     * @return Leaf node ids in topology order.
     */
    std::vector<NodeId> collectComponentTreeExtrema() const { return this->tree.getLeaves(); }

    /**
     * @brief Sentinel used for the dominant extremum.
     *
     * The dominant branch never meets a strictly stronger branch before the root.
     * A finite maximum sentinel preserves the ranking while avoiding infinities in
     * the public record list.
     *
     * @return `numeric_limits<Real>::max()`.
     */
    static Real dominantExtremumSentinel() noexcept { return std::numeric_limits<Real>::max(); }

    /**
     * @brief Validates top-k selection counts.
     *
     * @param extremaToKeep Number of strongest extrema requested.
     * @param context Name used in the exception message.
     * @throws std::invalid_argument If `extremaToKeep` is negative.
     */
    static void requireNonNegativeExtremaToKeep(int extremaToKeep, const char* context) {
        if (extremaToKeep < 0) {
            throw std::invalid_argument(std::string(context) + " requires a non-negative extremaToKeep value.");
        }
    }

    /**
     * @brief Validates extinction threshold selection values.
     *
     * @param threshold Minimum accepted extinction value.
     * @param context Name used in the exception message.
     * @throws std::invalid_argument If `threshold` is NaN or infinite.
     */
    static void requireFiniteThreshold(Real threshold, const char* context) {
        if (!std::isfinite(threshold)) {
            throw std::invalid_argument(std::string(context) + " requires a finite extinction threshold.");
        }
    }

    /**
     * @brief Validates the selected extrema policy before image operations.
     *
     * Top-k policies require non-negative counts. Threshold policies require a
     * finite threshold. Validation happens once at the selection boundary so
     * `filtering` and `contourMap` share the same acceptance rules.
     *
     * @param policy Selection policy supplied by the caller.
     * @param context Name used in the exception message.
     * @throws std::invalid_argument If the policy contains an invalid value.
     */
    static void requireValidSelectionPolicy(const ExtinctionSelectionPolicy<Real>& policy, const char* context) {
        if (policy.selectsTopK()) {
            requireNonNegativeExtremaToKeep(policy.extremaToKeep(), context);
        } else {
            requireFiniteThreshold(policy.minimumExtinction(), context);
        }
    }

    /**
     * @brief Validates extinction values before using them as formal valuations.
     *
     * Formal edge saliency maps require finite, non-negative hierarchy
     * valuations. The node id is accepted for diagnostic context and future
     * extension; the current message is intentionally compact.
     *
     * @param value Extinction value to validate.
     * @param nodeId Node associated with the value.
     * @param context Name used in the exception message.
     * @throws std::invalid_argument If `value` is not finite or is negative.
     */
    static void requireFormalExtinctionValue(Real value, NodeId nodeId, const char* context) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(std::string(context) + " requires finite extinction values.");
        }
        if (value < Real{0}) {
            throw std::invalid_argument(std::string(context) + " requires non-negative extinction values.");
        }
        (void)nodeId;
    }

    /**
     * @brief Extends regional-extremum extinction records to all hierarchy nodes.
     *
     * Each regional-extremum leaf receives its extinction value. Every ancestor
     * then receives the maximum extinction value among selected descendant
     * extrema, computed by a post-order propagation. The resulting dense
     * node-indexed attribute is non-decreasing toward the root and can be
     * validated/projected by `HierarchySaliencyMap`. This extension is an
     * intermediate hierarchical-watershed quantity; the persistence of a merge
     * is subsequently the minimum of its two child extinctions.
     */
    void buildExtinctionValueAttribute() {
        extinctionValueAttribute_.assign(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), Real{0});
        for (const RegionalExtremaNode<Real>& record : regionalExtremaNodes_) {
            Real& localValue = extinctionValueAttribute_[static_cast<std::size_t>(record.leaf)];
            if (localValue < record.extinction) {
                localValue = record.extinction;
            }
        }

        for (NodeId nodeId : tree.getPostOrderNodes()) {
            Real& nodeValue = extinctionValueAttribute_[static_cast<std::size_t>(nodeId)];
            for (NodeId childId : tree.getChildren(nodeId)) {
                const Real childValue = extinctionValueAttribute_[static_cast<std::size_t>(childId)];
                if (nodeValue < childValue) {
                    nodeValue = childValue;
                }
            }
        }
    }

    /**
     * @brief Checks the cached extinction hierarchy valuation.
     *
     * This validates both the regional-extremum record values and the dense
     * hierarchy attribute generated by `buildExtinctionValueAttribute`.
     *
     * @param context Name used in exception messages.
     * @throws std::invalid_argument If a record is invalid or if the dense
     * valuation violates the required non-negative hierarchy order.
     */
    void validateExtinctionValueAttribute(const char* context) const {
        for (const RegionalExtremaNode<Real>& record : regionalExtremaNodes_) {
            requireFormalExtinctionValue(record.extinction, record.leaf, context);
        }
        HierarchySaliencyMapValidation::validateHierarchyValuation(tree, std::span<const Real>(extinctionValueAttribute_),
                                                                   HierarchyValuationPolicy::AllowLevelCollapse,
                                                                   HierarchyValuationRangePolicy::RequireNonNegative, context);
    }

    /**
     * @brief Selects regional-extremum records according to a caller policy.
     *
     * The returned vector keeps the global decreasing extinction order. For a
     * top-k policy, it contains the first `k` sorted records. For a threshold
     * policy, it contains every record satisfying `extinction >= threshold`.
     * `rankScore` is assigned after selection, so contour rank scores are dense
     * over the selected set.
     *
     * @param policy Selection policy supplied to `filtering` or `contourMap`.
     * @param context Name used in exception messages.
     * @return Selected records plus their dense rank score.
     * @throws std::invalid_argument If `policy` is invalid.
     */
    [[nodiscard]] std::vector<SelectedExtremum> selectExtrema(const ExtinctionSelectionPolicy<Real>& policy, const char* context) const {
        requireValidSelectionPolicy(policy, context);

        std::vector<const RegionalExtremaNode<Real>*> records;
        if (policy.selectsTopK()) {
            const int extremaToSelect = std::min(policy.extremaToKeep(), static_cast<int>(regionalExtremaNodes_.size()));
            records.reserve(static_cast<std::size_t>(extremaToSelect));
            for (int i = 0; i < extremaToSelect; ++i) {
                records.push_back(&regionalExtremaNodes_[static_cast<std::size_t>(i)]);
            }
        } else {
            records.reserve(regionalExtremaNodes_.size());
            for (const RegionalExtremaNode<Real>& record : regionalExtremaNodes_) {
                if (record.extinction >= policy.minimumExtinction()) {
                    records.push_back(&record);
                }
            }
        }

        std::vector<SelectedExtremum> selected;
        selected.reserve(records.size());
        const int selectedCount = static_cast<int>(records.size());
        for (int i = 0; i < selectedCount; ++i) {
            selected.push_back(SelectedExtremum{records[static_cast<std::size_t>(i)], selectedCount - i});
        }
        return selected;
    }

    /**
     * @brief Reconstructs an image from a dense keep/remove node criterion.
     *
     * The criterion marks regional-extremum leaves that must be preserved. This
     * helper propagates preserved marks to ancestors, then traverses the tree
     * top-down. Preserved child branches are explored; unpreserved child subtrees
     * are collapsed to the current node altitude.
     *
     * @param criterion Dense internal-node boolean buffer; non-zero means keep.
     * @param context Name used in exception messages.
     * @return Reconstructed image on the original image domain.
     * @throws std::logic_error If the topology changed after construction.
     */
    [[nodiscard]] ImagePtr<T> filteringFromExtremumCriterion(std::vector<uint8_t> criterion, const char* context) const {
        requireStableTree(context);
        const AltitudeView altitudeView = view();
        for (NodeId nodeId : tree.getPostOrderNodes()) {
            if (!tree.isRoot(nodeId) && criterion[static_cast<std::size_t>(nodeId)]) {
                criterion[static_cast<std::size_t>(tree.getNodeParent(nodeId))] = true;
            }
        }

        ImagePtr<T> imgOutputPtr = Image<T>::create(tree.getNumRowsOfGridDomain2D(), tree.getNumColsOfGridDomain2D(), T{});
        auto imgOutput = imgOutputPtr->rawData();
        std::stack<NodeId> stack;
        stack.push(tree.getRoot());
        while (!stack.empty()) {
            const NodeId nodeId = stack.top();
            stack.pop();
            const T level = altitudeOf(altitudeView, nodeId);
            for (int pixel : tree.getProperParts(nodeId)) {
                imgOutput[pixel] = level;
            }
            for (NodeId childNodeId : tree.getChildren(nodeId)) {
                if (criterion[static_cast<std::size_t>(childNodeId)]) {
                    stack.push(childNodeId);
                } else {
                    for (NodeId subtreeNodeId : tree.getNodeSubtree(childNodeId)) {
                        for (int pixel : tree.getProperParts(subtreeNodeId)) {
                            imgOutput[pixel] = level;
                        }
                    }
                }
            }
        }
        return imgOutputPtr;
    }

    /**
     * @brief Reconstructs an image from selected regional-extremum records.
     *
     * This adapter converts selected records to the dense leaf criterion consumed
     * by `filteringFromExtremumCriterion`.
     *
     * @param selected Regional extrema chosen by `selectExtrema`.
     * @param context Name used in exception messages.
     * @return Reconstructed image on the original image domain.
     */
    [[nodiscard]] ImagePtr<T> filteringFromSelectedExtrema(const std::vector<SelectedExtremum>& selected, const char* context) const {
        std::vector<uint8_t> criterion(tree.getNumInternalNodeSlots(), false);
        for (const SelectedExtremum& item : selected) {
            criterion[static_cast<std::size_t>(item.record->leaf)] = true;
        }
        return filteringFromExtremumCriterion(std::move(criterion), context);
    }

    /**
     * @brief Computes extinction records from the supplied increasing attribute.
     *
     * This method implements the component-tree extinction algorithm described in:
     *
     * Alexandre Gonçalves Silva and Roberto de Alencar Lotufo,
     * "Efficient computation of new extinction values from extended component
     * tree," Pattern Recognition Letters, 32(1):79-90, 2011,
     * [DOI 10.1016/j.patrec.2010.07.019](https://doi.org/10.1016/j.patrec.2010.07.019).
     *
     * The implementation processes component-tree leaves as regional extrema,
     * climbs each leaf branch until it reaches the first merge with a stronger
     * or tie-resolved already visited branch, stores the corresponding cutoff
     * node and extinction value, and then sorts records in decreasing extinction
     * order. The dominant extremum, which survives to the root, receives the
     * finite sentinel returned by `dominantExtremumSentinel()`.
     *
     * @param attr Non-null dense attribute buffer indexed by internal `NodeId`.
     * Larger values must represent stronger extrema.
     * @throws std::invalid_argument If `attr` is null.
     */
    void initialize(const Real* attr) {
        requireAttributePointer(attr, "ExtinctionValues");
        std::vector<NodeId> leaves = collectComponentTreeExtrema();
        regionalExtremaNodes_.reserve(leaves.size());
        std::vector<uint8_t> visited(this->tree.getNumInternalNodeSlots(), false);
        for (NodeId leafNodeId : leaves) {
            Real extinction = dominantExtremumSentinel();
            NodeId cutoffNodeId = leafNodeId;
            NodeId parentNodeId = this->tree.getNodeParent(cutoffNodeId);
            bool flag = true;
            while (flag && !this->tree.isRoot(cutoffNodeId)) {
                if (this->tree.getNumChildren(parentNodeId) > 1) {
                    for (NodeId sonNodeId : this->tree.getChildren(parentNodeId)) {
                        if (flag) {
                            if (visited[sonNodeId] && sonNodeId != cutoffNodeId && attr[sonNodeId] == attr[cutoffNodeId]) {
                                flag = false;
                            } else if (sonNodeId != cutoffNodeId && attr[sonNodeId] > attr[cutoffNodeId]) {
                                flag = false;
                            }
                            visited[sonNodeId] = true;
                        }
                    }
                }
                if (flag) {
                    cutoffNodeId = parentNodeId;
                    parentNodeId = this->tree.getNodeParent(cutoffNodeId);
                }
            }
            if (!this->tree.isRoot(cutoffNodeId)) {
                extinction = attr[cutoffNodeId];
            }
            regionalExtremaNodes_.emplace_back(leafNodeId, cutoffNodeId, extinction);
        }

        std::sort(regionalExtremaNodes_.begin(), regionalExtremaNodes_.end(), [](const auto& a, const auto& b) {
            if (a.extinction != b.extinction) {
                return a.extinction > b.extinction;
            }
            if (a.cutoffNode != b.cutoffNode) {
                return a.cutoffNode < b.cutoffNode;
            }
            return a.leaf < b.leaf;
        });
        buildExtinctionValueAttribute();
    }
    /// @endcond

  public:
    /// Scalar type used for input attributes and extinction values.
    using value_type = Real;

    /**
     * @brief Computes extinction values from a weighted view and shared attribute buffer.
     *
     * @param view Weighted tree view whose topology and altitude define the
     * reconstruction domain.
     * @param attr Shared buffer with one scalar attribute value per internal
     * node slot.
     * @throws std::invalid_argument If `attr` is null.
     * @throws std::logic_error If `view` is stale.
     */
    ExtinctionValues(const AltitudeView& view, const std::shared_ptr<Real[]>& attr) : ExtinctionValues(view, attr.get()) {}

    /**
     * @brief Computes extinction values from a weighted view and vector attribute buffer.
     *
     * @param view Weighted tree view whose topology and altitude define the
     * reconstruction domain.
     * @param attr Vector with one scalar attribute value per internal node slot.
     * @throws std::invalid_argument If `attr` does not match the internal node
     * slot count of `view.topology()`.
     * @throws std::logic_error If `view` is stale.
     */
    ExtinctionValues(const AltitudeView& view, const std::vector<Real>& attr)
        : ExtinctionValues(view, requireAttributeBuffer(view.topology(), attr, "ExtinctionValues")) {}

    /**
     * @brief Computes extinction values from a weighted view and raw attribute buffer.
     *
     * @param view Weighted tree view whose topology and altitude define the
     * reconstruction domain.
     * @param attr Non-null buffer indexed by dense internal `NodeId`.
     * @throws std::invalid_argument If `attr` is null.
     * @throws std::logic_error If `view` is stale.
     */
    ExtinctionValues(const AltitudeView& view, const Real* attr) : view_(view), tree(view_.topology()), treeMutationVersion_(tree.getMutationVersion()) {
        view_.requireTopologyUnchanged("ExtinctionValues");
        detail::validateGlobalMonotoneAltitudeOrder(this->tree, "ExtinctionValues");
        initialize(attr);
    }

    /**
     * @brief Computes extinction values from a weighted tree and shared attribute buffer.
     *
     * The weighted tree is borrowed; it must outlive this object.
     *
     * @param weighted Weighted tree whose topology and altitude define the
     * reconstruction domain.
     * @param attr Shared buffer with one scalar attribute value per internal
     * node slot.
     * @throws std::invalid_argument If `attr` is null.
     */
    ExtinctionValues(const WeightedMorphologicalTree<T>& weighted, const std::shared_ptr<Real[]>& attr) : ExtinctionValues(weighted.asView(), attr.get()) {
        weighted_ = &weighted;
    }

    /**
     * @brief Computes extinction values from a weighted tree and vector attribute buffer.
     *
     * The weighted tree is borrowed; it must outlive this object.
     *
     * @param weighted Weighted tree whose topology and altitude define the
     * reconstruction domain.
     * @param attr Vector with one scalar attribute value per internal node slot.
     * @throws std::invalid_argument If `attr` does not match the internal node
     * slot count of the tree.
     */
    ExtinctionValues(const WeightedMorphologicalTree<T>& weighted, const std::vector<Real>& attr) : ExtinctionValues(weighted.asView(), attr) {
        weighted_ = &weighted;
    }

    /**
     * @brief Computes extinction values from a weighted tree and raw attribute buffer.
     *
     * The weighted tree is borrowed; it must outlive this object.
     *
     * @param weighted Weighted tree whose topology and altitude define the
     * reconstruction domain.
     * @param attr Non-null buffer indexed by dense internal `NodeId`.
     * @throws std::invalid_argument If `attr` is null.
     */
    ExtinctionValues(const WeightedMorphologicalTree<T>& weighted, const Real* attr) : ExtinctionValues(weighted.asView(), attr) { weighted_ = &weighted; }

    /**
     * @brief Builds a contour-valued image from selected extinction events.
     *
     * `contourMap` is an image-domain visualization, not the formal edge-indexed
     * saliency map of a hierarchy. The selection policy chooses which regional
     * extrema are kept, and the score policy chooses the value written on each
     * retained cutoff-node contour. If several selected extrema share the same
     * cutoff node, the strongest score for that node is kept.
     *
     * @param selection Selection policy shared with `filtering`.
     * @param scorePolicy Value convention for retained contours.
     * @return Floating-point image on the original image domain.
     * @throws std::invalid_argument If `selection` is invalid.
     * @throws std::logic_error If the tree topology changed after construction.
     */
    [[nodiscard]] ImagePtr<Real> contourMap(const ExtinctionSelectionPolicy<Real>& selection, ExtinctionContourScorePolicy scorePolicy) const {
        constexpr const char* context = "ExtinctionValues::contourMap";
        requireStableTree(context);
        const std::vector<SelectedExtremum> selected = selectExtrema(selection, context);

        std::vector<uint8_t> keep(tree.getNumInternalNodeSlots(), false);
        std::vector<Real> scoreByNode(tree.getNumInternalNodeSlots(), Real{0});
        std::vector<NodeId> keptNodes;
        for (const SelectedExtremum& item : selected) {
            const NodeId cutoffNode = item.record->cutoffNode;
            Real score = Real{0};
            switch (scorePolicy) {
            case ExtinctionContourScorePolicy::RankScore:
                score = static_cast<Real>(item.rankScore);
                break;
            case ExtinctionContourScorePolicy::ExtinctionValue:
                score = item.record->extinction;
                break;
            default:
                throw std::invalid_argument(std::string(context) + " received an unknown extinction contour score policy.");
            }
            if (!keep[static_cast<std::size_t>(cutoffNode)]) {
                keptNodes.push_back(cutoffNode);
                keep[static_cast<std::size_t>(cutoffNode)] = true;
                scoreByNode[static_cast<std::size_t>(cutoffNode)] = score;
            } else if (scoreByNode[static_cast<std::size_t>(cutoffNode)] < score) {
                scoreByNode[static_cast<std::size_t>(cutoffNode)] = score;
            }
        }

        ImagePtr<Real> imgOutputPtr = Image<Real>::create(tree.getNumRowsOfGridDomain2D(), tree.getNumColsOfGridDomain2D(), Real{0});
        auto contourOutput = imgOutputPtr->rawData();

        auto contours = ContoursComputedIncrementally::extractCompactContours(tree);
        for (NodeId node : keptNodes) {
            for (int p : contours.getContour(node)) {
                contourOutput[p] = scoreByNode[static_cast<std::size_t>(node)];
            }
        }

        return imgOutputPtr;
    }

    /**
     * @brief Returns extinction values extended from extrema to every hierarchy node.
     *
     * The extinction records associate one value with each leaf extremum and store
     * the cutoff node where that extremum stops being dominant. This method follows
     * the standard extinction-attribute extension used for hierarchy saliency: each
     * extremum leaf receives its extinction value, and every non-leaf node receives
     * the maximum extinction value among the extrema contained in its subtree.
     *
     * The resulting dense node attribute is computed during initialization,
     * non-decreasing toward the root. It can be passed to
     * `HierarchySaliencyMap::computeSaliencyEdgeMap` when the monotone node
     * projection is explicitly desired. The Cousty hierarchical-watershed path
     * is `computeFormalSaliencyEdgeMap`, which additionally computes merge
     * persistences. Raw extinction values are
     * required to be finite and non-negative; use
     * `computeRankedExtinctionValueAttribute` when a compact integer level scale
     * is preferable to the raw values and the dominant-extremum sentinel.
     *
     * @return Cached dense node extinction attribute.
     * @throws std::logic_error If the tree topology changed after construction.
     * @throws std::invalid_argument If cached values are not a valid non-negative
     * hierarchy valuation.
     */
    [[nodiscard]] const std::vector<Real>& getExtinctionValueAttribute() const {
        constexpr const char* context = "ExtinctionValues::getExtinctionValueAttribute";
        requireStableTree(context);
        validateExtinctionValueAttribute(context);
        return extinctionValueAttribute_;
    }

    /**
     * @brief Builds a dense integer extinction attribute from extinction levels.
     *
     * Raw extinction values can include the dominant-extremum sentinel
     * `numeric_limits<Real>::max()`, which is useful for ordering records but
     * awkward as a display scale. This helper preserves the induced hierarchy
     * order while replacing distinct live-node values by dense ranks `0..k-1`.
     * Equal extinction levels remain equal, so intentional level collapse is
     * preserved.
     *
     * @return Dense integer hierarchy valuation indexed by internal `NodeId`.
     * @throws std::logic_error If the tree topology changed after construction.
     * @throws std::invalid_argument If the raw extinction valuation is invalid.
     */
    [[nodiscard]] std::vector<int> computeRankedExtinctionValueAttribute() const {
        const std::vector<Real>& valuation = getExtinctionValueAttribute();
        return HierarchySaliencyMapValidation::rankHierarchyValuation(tree, std::span<const Real>(valuation), HierarchyValuationPolicy::AllowLevelCollapse);
    }

    /**
     * @brief Computes the formal hierarchical-watershed extinction saliency map.
     *
     * This is the Section 8.1 persistence path of Cousty et al.: a Kruskal BPTAO
     * is represented by an MST, each binary merge receives the minimum of the
     * maximum descendant extinctions of its two children, and the resulting
     * persistence-weighted MST is converted to its full-graph QFZ saliency map.
     * The return type is an edge map, not the image-domain contour visualization
     * returned by `contourMap`.
     *
     * @par Primary reference
     * Jean Cousty, Laurent Najman, Yukiko Kenmochi, and Silvio Guimarães,
     * "Hierarchical segmentations with graphs: quasi-flat zones, minimum spanning
     * trees, and saliency maps," Journal of Mathematical Imaging and Vision,
     * 60(4):479-502, 2018,
     * [DOI 10.1007/s10851-017-0768-7](https://doi.org/10.1007/s10851-017-0768-7),
     * Section 8.1.
     *
     * @param adjacency Explicit adjacency relation used to enumerate image-domain
     * graph edges.
     * @return Edge-indexed saliency map with raw extinction values.
     * @throws std::logic_error If the tree topology changed after construction.
     * @throws std::invalid_argument If the raw extinction valuation is invalid.
     */
    [[nodiscard]] EdgeSaliencyMap<Real> computeFormalSaliencyEdgeMap(const RegularGridAdjacency2D& adjacency) const {
        const std::vector<Real>& valuation = getExtinctionValueAttribute();
        return HierarchicalWatershedSaliency::compute(view(), std::span<const Real>(valuation), adjacency);
    }

    /**
     * @brief Computes the formal extinction saliency edge map using stored adjacency.
     *
     * The weighted tree must carry a construction adjacency relation. Use the
     * overload receiving `RegularGridAdjacency2D` when the tree was imported without
     * stored adjacency metadata.
     *
     * @return Edge-indexed saliency map with raw extinction values.
     * @throws std::logic_error If the tree topology changed after construction or
     * no stored adjacency is available.
     * @throws std::invalid_argument If the raw extinction valuation is invalid.
     */
    [[nodiscard]] EdgeSaliencyMap<Real> computeFormalSaliencyEdgeMap() const {
        return computeFormalSaliencyEdgeMap(HierarchySaliencyMap::requireProjectionAdjacency(tree, "ExtinctionValues::computeFormalSaliencyEdgeMap"));
    }

    /**
     * @brief Computes a ranked formal extinction saliency edge map.
     *
     * This variant first computes the persistence-based hierarchical watershed,
     * then ranks only distinct values that occur on its final graph edges. The
     * dominant-extremum sentinel therefore cannot create an unused rank.
     *
     * @param adjacency Explicit adjacency relation used to enumerate image-domain
     * graph edges.
     * @return Edge-indexed saliency map with dense integer extinction ranks.
     * @throws std::logic_error If the tree topology changed after construction.
     * @throws std::invalid_argument If the raw extinction valuation is invalid.
     */
    [[nodiscard]] EdgeSaliencyMap<int> computeRankedFormalSaliencyEdgeMap(const RegularGridAdjacency2D& adjacency) const {
        const std::vector<Real>& valuation = getExtinctionValueAttribute();
        return HierarchicalWatershedSaliency::computeRanked(view(), std::span<const Real>(valuation), adjacency);
    }

    /**
     * @brief Computes a ranked formal extinction saliency edge map using stored adjacency.
     *
     * The weighted tree must carry a construction adjacency relation. Use the
     * overload receiving `RegularGridAdjacency2D` when the tree was imported without
     * stored adjacency metadata.
     *
     * @return Edge-indexed saliency map with dense integer extinction ranks.
     * @throws std::logic_error If the tree topology changed after construction or
     * no stored adjacency is available.
     * @throws std::invalid_argument If the raw extinction valuation is invalid.
     */
    [[nodiscard]] EdgeSaliencyMap<int> computeRankedFormalSaliencyEdgeMap() const {
        return computeRankedFormalSaliencyEdgeMap(
            HierarchySaliencyMap::requireProjectionAdjacency(tree, "ExtinctionValues::computeRankedFormalSaliencyEdgeMap"));
    }

    /**
     * @brief Projects the max-descendant extinction attribute directly by LCA.
     *
     * This method preserves the pre-correction behavior of
     * `computeFormalSaliencyEdgeMap`. It induces a valid monotone hierarchy, but
     * it is not the hierarchical-watershed persistence construction.
     *
     * @param adjacency Explicit adjacency relation used to enumerate image-domain
     * graph edges.
     * @return Edge-indexed map obtained by LCA projection of the max-descendant
     * extinction valuation.
     * @throws std::logic_error If the tree topology changed after construction.
     * @throws std::invalid_argument If the valuation or projection graph is invalid.
     */
    [[nodiscard]] EdgeSaliencyMap<Real> computeMonotoneExtinctionProjection(const RegularGridAdjacency2D& adjacency) const {
        const std::vector<Real>& valuation = getExtinctionValueAttribute();
        return HierarchySaliencyMap::computeSaliencyEdgeMap(tree, adjacency, std::span<const Real>(valuation), HierarchyValuationPolicy::AllowLevelCollapse);
    }

    /**
     * @brief Stored-adjacency overload of `computeMonotoneExtinctionProjection`.
     *
     * @return Edge-indexed map obtained by monotone LCA projection.
     * @throws std::logic_error If the topology changed or no unambiguous stored
     * adjacency is available.
     * @throws std::invalid_argument If the valuation or projection graph is invalid.
     */
    [[nodiscard]] EdgeSaliencyMap<Real> computeMonotoneExtinctionProjection() const {
        return computeMonotoneExtinctionProjection(
            HierarchySaliencyMap::requireProjectionAdjacency(tree, "ExtinctionValues::computeMonotoneExtinctionProjection"));
    }

    /**
     * @brief Computes canonical effective-edge ranks for the monotone projection.
     *
     * @param adjacency Explicit adjacency relation used to enumerate image-domain
     * graph edges.
     * @return Edge-indexed monotone projection with dense effective-edge ranks.
     * @throws std::logic_error If the tree topology changed after construction.
     * @throws std::invalid_argument If the valuation or projection graph is invalid.
     */
    [[nodiscard]] EdgeSaliencyMap<int> computeRankedMonotoneExtinctionProjection(const RegularGridAdjacency2D& adjacency) const {
        const std::vector<Real>& valuation = getExtinctionValueAttribute();
        return HierarchySaliencyMap::computeCanonicalRankedSaliencyEdgeMap(tree, adjacency, std::span<const Real>(valuation),
                                                                           HierarchyValuationPolicy::AllowLevelCollapse);
    }

    /**
     * @brief Stored-adjacency overload of `computeRankedMonotoneExtinctionProjection`.
     *
     * @return Edge-indexed monotone projection with dense effective-edge ranks.
     * @throws std::logic_error If the topology changed or no unambiguous stored
     * adjacency is available.
     * @throws std::invalid_argument If the valuation or projection graph is invalid.
     */
    [[nodiscard]] EdgeSaliencyMap<int> computeRankedMonotoneExtinctionProjection() const {
        return computeRankedMonotoneExtinctionProjection(
            HierarchySaliencyMap::requireProjectionAdjacency(tree, "ExtinctionValues::computeRankedMonotoneExtinctionProjection"));
    }

    /**
     * @brief Reconstructs an image from selected regional extrema.
     *
     * @param selection Selection policy used to retain extrema by rank or by
     * extinction threshold.
     * @return Image on the original image domain using altitude type `T`.
     * @throws std::invalid_argument If `selection` is invalid.
     * @throws std::logic_error If the tree topology changed after construction.
     */
    [[nodiscard]] ImagePtr<T> filtering(const ExtinctionSelectionPolicy<Real>& selection) const {
        constexpr const char* context = "ExtinctionValues::filtering";
        const std::vector<SelectedExtremum> selected = selectExtrema(selection, context);
        return filteringFromSelectedExtrema(selected, context);
    }

    /**
     * @brief Returns regional-extremum records sorted by decreasing extinction.
     *
     * @return Immutable record vector kept by this object.
     * @throws std::logic_error If the tree topology changed after construction.
     */
    [[nodiscard]] const std::vector<RegionalExtremaNode<Real>>& getRegionalExtrema() const {
        requireStableTree("ExtinctionValues::getRegionalExtrema");
        return regionalExtremaNodes_;
    }
};

} // namespace mmcfilters
