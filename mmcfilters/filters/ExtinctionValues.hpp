#pragma once

#include "../trees/TreeAltitudeAlgorithms.hpp"
#include "../trees/WeightedMorphologicalTree.hpp"
#include "../trees/WeightedTreeView.hpp"
#include "../trees/detail/TreeKindValidation.hpp"
#include "../utils/Image.hpp"
#include "../utils/Common.hpp"
#include "../contours/ContoursComputedIncrementally.hpp"

#include <algorithm>
#include <concepts>
#include <limits>
#include <memory>
#include <stack>
#include <stdexcept>
#include <string>
#include <vector>

namespace mmcfilters {

/**
 * @brief Record describing one regional extremum and its extinction value.
 */
template <std::floating_point Real = float>
struct RegionalExtremaNode {
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
    RegionalExtremaNode(NodeId leaf, NodeId cutoffNode, Real extinction)
        : leaf(leaf), cutoffNode(cutoffNode), extinction(extinction) {}
};

/**
 * @brief Computes and stores extinction values for regional extrema.
 *
 * @details
 * `ExtinctionValues` implements the classical leaf-extrema extinction ranking
 * for max-trees and min-trees. In this component-tree setting, the regional
 * extrema processed by the algorithm are the tree leaves. The supplied scalar
 * attribute is indexed by dense internal `NodeId`, must have one value for every
 * internal node slot, and is interpreted so that larger values represent
 * stronger extrema. Results are sorted in decreasing extinction order and can be
 * consumed either as records, a filtered reconstruction, or a contour saliency
 * map.
 *
 * The strongest extremum has no stronger merge point. Its extinction value is
 * represented by the explicit finite sentinel `numeric_limits<Real>::max()`.
 *
 * Trees of shapes and self-dual residual trees are intentionally rejected by the
 * public constructors because their complete regional-extrema set is not
 * generally equivalent to `tree.getLeaves()`.
 *
 * The object records the tree mutation version at construction time. Public
 * operations reject use after the underlying topology changes.
 *
 * @tparam T Altitude type used by the weighted tree or weighted view.
 * @tparam Real Attribute-buffer floating-point type.
 */
template<AltitudeValue T, std::floating_point Real = float>
class ExtinctionValues {
protected:
    /// @cond INTERNAL
    using AltitudeView = WeightedTreeView<T>;

    std::vector<RegionalExtremaNode<Real>> regionalExtremaNodes;
    AltitudeView view_;
    const WeightedMorphologicalTree<T>* weighted_ = nullptr;
    const MorphologicalTree& tree;
    std::size_t treeMutationVersion_ = 0;

    AltitudeView view() const {
        return weighted_ != nullptr ? weighted_->asView() : view_;
    }

    void requireStableTree(const char* context) const {
        tree.requireMutationVersion(treeMutationVersion_, context);
    }

    static void requireAttributePointer(const Real* attr, const char* context) {
        if (attr == nullptr) {
            throw std::invalid_argument(std::string(context) + " requires a non-null attribute buffer.");
        }
    }

    static const Real* requireAttributeBuffer(const MorphologicalTree& tree, const std::vector<Real>& attr, const char* context) {
        if (attr.size() != static_cast<std::size_t>(tree.getNumInternalNodeSlots())) {
            throw std::invalid_argument(std::string(context) + " attribute size must match the internal node slot count.");
        }
        return attr.data();
    }

    static T altitudeOf(const AltitudeView& view, NodeId nodeId) {
        return view.getAltitude(nodeId);
    }

    std::vector<NodeId> collectComponentTreeExtrema() const {
        return this->tree.getLeaves();
    }

    static Real dominantExtremumSentinel() noexcept {
        return std::numeric_limits<Real>::max();
    }

    static void requireNonNegativeExtremaToKeep(int extremaToKeep, const char* context) {
        if (extremaToKeep < 0) {
            throw std::invalid_argument(std::string(context) + " requires a non-negative extremaToKeep value.");
        }
    }

    void initialize(const Real* attr) {
        requireAttributePointer(attr, "ExtinctionValues");
        std::vector<NodeId> leaves = collectComponentTreeExtrema();
        regionalExtremaNodes.reserve(leaves.size());
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
            regionalExtremaNodes.emplace_back(leafNodeId, cutoffNodeId, extinction);
        }

        std::sort(regionalExtremaNodes.begin(), regionalExtremaNodes.end(), [](const auto& a, const auto& b) {
            if (a.extinction != b.extinction) {
                return a.extinction > b.extinction;
            }
            if (a.cutoffNode != b.cutoffNode) {
                return a.cutoffNode < b.cutoffNode;
            }
            return a.leaf < b.leaf;
        });
    }
    /// @endcond

public:
    /**
     * @brief Computes extinction values from a weighted view and shared attribute buffer.
     *
     * @param view Weighted tree view whose topology and altitude define the
     * reconstruction domain.
     * @param attr Shared buffer with one scalar attribute value per internal
     * node slot.
     */
    ExtinctionValues(const AltitudeView& view, const std::shared_ptr<Real[]>& attr)
        : ExtinctionValues(view, attr.get()) {}

    /**
     * @brief Computes extinction values from a weighted view and vector attribute buffer.
     *
     * @throws std::invalid_argument If `attr` does not match the internal node
     * slot count of `view.topology()`.
     */
    ExtinctionValues(const AltitudeView& view, const std::vector<Real>& attr)
        : ExtinctionValues(view, requireAttributeBuffer(view.topology(), attr, "ExtinctionValues")) {}

    /**
     * @brief Computes extinction values from a weighted view and raw attribute buffer.
     *
     * @param view Weighted tree view whose topology and altitude define the
     * reconstruction domain.
     * @param attr Non-null buffer indexed by dense internal `NodeId`.
     * @throws std::invalid_argument If `attr` is null or if the view topology is stale.
     */
    ExtinctionValues(const AltitudeView& view, const Real* attr)
        : view_(view),
          tree(view_.topology()),
          treeMutationVersion_(tree.getMutationVersion()) {
        view_.requireTopologyUnchanged("ExtinctionValues");
        detail::validateComponentTreeKind(this->tree, "ExtinctionValues");
        initialize(attr);
    }

    /**
     * @brief Computes extinction values from a weighted tree and shared attribute buffer.
     *
     * The weighted tree is borrowed; it must outlive this object.
     */
    ExtinctionValues(const WeightedMorphologicalTree<T>& weighted, const std::shared_ptr<Real[]>& attr)
        : ExtinctionValues(weighted.asView(), attr.get()) {
        weighted_ = &weighted;
    }

    /**
     * @brief Computes extinction values from a weighted tree and vector attribute buffer.
     *
     * The weighted tree is borrowed; it must outlive this object.
     *
     * @throws std::invalid_argument If `attr` does not match the internal node
     * slot count of the tree.
     */
    ExtinctionValues(const WeightedMorphologicalTree<T>& weighted, const std::vector<Real>& attr)
        : ExtinctionValues(weighted.asView(), attr) {
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
     */
    ExtinctionValues(const WeightedMorphologicalTree<T>& weighted, const Real* attr)
        : ExtinctionValues(weighted.asView(), attr) {
        weighted_ = &weighted;
    }

    /**
     * @brief Builds a contour saliency image from the strongest extrema.
     *
     * @param extremaToKeep Maximum number of extrema retained from the
     * descending extinction ranking.
     * @param unweighted When true, contours receive rank-like scores; when
     * false, contours receive their extinction values.
     * @return Floating-point image on the original image domain.
     * @throws std::logic_error If the tree topology changed after construction.
     */
    [[nodiscard]] ImagePtr<Real> saliencyMap(int extremaToKeep, bool unweighted = true) {
        requireStableTree("ExtinctionValues::saliencyMap");
        requireNonNegativeExtremaToKeep(extremaToKeep, "ExtinctionValues::saliencyMap");
        std::vector<uint8_t> keep(tree.getNumInternalNodeSlots(), false);
        std::vector<Real> extinctionByNode(tree.getNumInternalNodeSlots(), Real{0});
        std::vector<NodeId> keptNodes;
        const int leafToKeep = std::min(extremaToKeep, static_cast<int>(regionalExtremaNodes.size()));
        for (int i = 0; i < leafToKeep; ++i) {
            const NodeId cutoffNode = this->regionalExtremaNodes[i].cutoffNode;
            if (!keep[cutoffNode]) {
                keptNodes.push_back(cutoffNode);
            }
            keep[cutoffNode] = true;
            extinctionByNode[cutoffNode] = unweighted ? static_cast<Real>(leafToKeep - i) : this->regionalExtremaNodes[i].extinction;
        }

        ImagePtr<Real> imgOutputPtr = Image<Real>::create(tree.getNumRowsOfImage(), tree.getNumColsOfImage(), Real{0});
        auto saliencyOutput = imgOutputPtr->rawData();

        auto contours = ContoursComputedIncrementally::extractCompactContours(tree);
        for (NodeId node : keptNodes) {
            for (int p : contours.getContour(node)) {
                saliencyOutput[p] = extinctionByNode[node];
            }
        }

        return imgOutputPtr;
    }

    /**
     * @brief Reconstructs an image by keeping the strongest extrema.
     *
     * @param extremaToKeep Maximum number of leaf extrema retained from the
     * descending extinction ranking.
     * @return Image on the original image domain using altitude type `T`.
     * @throws std::logic_error If the tree topology changed after construction.
     */
    [[nodiscard]] ImagePtr<T> filtering(int extremaToKeep) {
        requireStableTree("ExtinctionValues::filtering");
        requireNonNegativeExtremaToKeep(extremaToKeep, "ExtinctionValues::filtering");
        const AltitudeView altitudeView = view();
        std::vector<uint8_t> criterion(tree.getNumInternalNodeSlots(), false);
        const int leafToKeep = std::min(extremaToKeep, static_cast<int>(regionalExtremaNodes.size()));
        for (int i = 0; i < leafToKeep; i++) {
            criterion[regionalExtremaNodes[i].leaf] = true;
        }
        for (NodeId nodeId : tree.getPostOrderNodes()) {
            if (!tree.isRoot(nodeId) && criterion[nodeId]) {
                criterion[tree.getNodeParent(nodeId)] = true;
            }
        }

        ImagePtr<T> imgOutputPtr = Image<T>::create(tree.getNumRowsOfImage(), tree.getNumColsOfImage(), T{});
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
                if (criterion[childNodeId]) {
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
     * @brief Returns the extinction records sorted by decreasing extinction.
     *
     * @return Mutable record vector kept by this object.
     * @throws std::logic_error If the tree topology changed after construction.
     */
    std::vector<RegionalExtremaNode<Real>>& getExtinctionValues() {
        requireStableTree("ExtinctionValues::getExtinctionValues");
        return regionalExtremaNodes;
    }
};


} // namespace mmcfilters
