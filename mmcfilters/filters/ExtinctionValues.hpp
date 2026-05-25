#pragma once

#include "../trees/TreeAltitudeAlgorithms.hpp"
#include "../trees/WeightedMorphologicalTree.hpp"
#include "../trees/WeightedTreeView.hpp"
#include "../utils/Image.hpp"
#include "../utils/Common.hpp"
#include "../contours/ContoursComputedIncrementally.hpp"

#include <algorithm>
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
struct RegionalExtremaNode {
    /// Leaf node that represents the regional extremum.
    NodeId leaf;

    /// Highest node retained before the extremum merges with a stronger branch.
    NodeId cutoffNode;

    /// Attribute value at `cutoffNode`, or infinity-like sentinel at the root.
    float extinction;

    /**
     * @brief Builds one extinction-value record.
     *
     * @param leaf Leaf node representing the regional extremum.
     * @param cutoffNode Node where the extremum stops being dominant.
     * @param extinction Attribute value associated with the cutoff.
     */
    RegionalExtremaNode(NodeId leaf, NodeId cutoffNode, float extinction)
        : leaf(leaf), cutoffNode(cutoffNode), extinction(extinction) {}
};

/**
 * @brief Computes and stores extinction values for regional extrema.
 *
 * @details
 * `ExtinctionValues` ranks tree leaves by the persistence of a supplied scalar
 * attribute. The attribute buffer is indexed by dense internal `NodeId` and must
 * have one value for every internal node slot of the tree. Results are sorted in
 * decreasing extinction order and can be consumed either as records, a filtered
 * reconstruction, or a contour saliency map.
 *
 * The object records the tree mutation version at construction time. Public
 * operations reject use after the underlying topology changes.
 *
 * @tparam T Altitude type used by the weighted tree or weighted view.
 */
template<AltitudeValue T>
class ExtinctionValues {
protected:
    /// @cond INTERNAL
    using AltitudeView = WeightedTreeView<T>;

    std::vector<RegionalExtremaNode> regionalExtremaNodes;
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

    static void requireAttributePointer(const float* attr, const char* context) {
        if (attr == nullptr) {
            throw std::invalid_argument(std::string(context) + " requires a non-null attribute buffer.");
        }
    }

    static const float* requireAttributeBuffer(const MorphologicalTree& tree, const std::vector<float>& attr, const char* context) {
        if (attr.size() != static_cast<std::size_t>(tree.getNumInternalNodeSlots())) {
            throw std::invalid_argument(std::string(context) + " attribute size must match the internal node slot count.");
        }
        return attr.data();
    }

    static T altitudeOf(const AltitudeView& view, NodeId nodeId) {
        return view.getAltitude(nodeId);
    }

    void initialize(const float* attr) {
        requireAttributePointer(attr, "ExtinctionValues");
        std::vector<NodeId> leaves = this->tree.getLeaves();
        regionalExtremaNodes.reserve(leaves.size());
        std::vector<uint8_t> visited(this->tree.getNumInternalNodeSlots(), false);
        for (NodeId leafNodeId : leaves) {
            float extinction = std::numeric_limits<float>::max();
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
            return a.extinction > b.extinction;
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
    ExtinctionValues(const AltitudeView& view, const std::shared_ptr<float[]>& attr)
        : ExtinctionValues(view, attr.get()) {}

    /**
     * @brief Computes extinction values from a weighted view and vector attribute buffer.
     *
     * @throws std::invalid_argument If `attr` does not match the internal node
     * slot count of `view.topology()`.
     */
    ExtinctionValues(const AltitudeView& view, const std::vector<float>& attr)
        : ExtinctionValues(view, requireAttributeBuffer(view.topology(), attr, "ExtinctionValues")) {}

    /**
     * @brief Computes extinction values from a weighted view and raw attribute buffer.
     *
     * @param view Weighted tree view whose topology and altitude define the
     * reconstruction domain.
     * @param attr Non-null buffer indexed by dense internal `NodeId`.
     * @throws std::invalid_argument If `attr` is null or if the view topology is stale.
     */
    ExtinctionValues(const AltitudeView& view, const float* attr)
        : view_(view),
          tree(view_.topology()),
          treeMutationVersion_(tree.getMutationVersion()) {
        view_.requireTopologyUnchanged("ExtinctionValues");
        initialize(attr);
    }

    /**
     * @brief Computes extinction values from a weighted tree and shared attribute buffer.
     *
     * The weighted tree is borrowed; it must outlive this object.
     */
    ExtinctionValues(const WeightedMorphologicalTree<T>& weighted, const std::shared_ptr<float[]>& attr)
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
    ExtinctionValues(const WeightedMorphologicalTree<T>& weighted, const std::vector<float>& attr)
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
    ExtinctionValues(const WeightedMorphologicalTree<T>& weighted, const float* attr)
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
     * @return Float image on the original image domain.
     * @throws std::logic_error If the tree topology changed after construction.
     */
    [[nodiscard]] ImageFloatPtr saliencyMap(int extremaToKeep, bool unweighted = true) {
        requireStableTree("ExtinctionValues::saliencyMap");
        std::vector<uint8_t> keep(tree.getNumInternalNodeSlots(), false);
        std::vector<float> extinctionByNode(tree.getNumInternalNodeSlots(), 0.0f);
        std::vector<NodeId> keptNodes;
        const int leafToKeep = std::min(extremaToKeep, static_cast<int>(regionalExtremaNodes.size()));
        for (int i = 0; i < leafToKeep; ++i) {
            const NodeId cutoffNode = this->regionalExtremaNodes[i].cutoffNode;
            if (!keep[cutoffNode]) {
                keptNodes.push_back(cutoffNode);
            }
            keep[cutoffNode] = true;
            extinctionByNode[cutoffNode] = unweighted ? static_cast<float>(leafToKeep - i) : this->regionalExtremaNodes[i].extinction;
        }

        ImageFloatPtr imgOutputPtr = ImageFloat::create(tree.getNumRowsOfImage(), tree.getNumColsOfImage(), 0);
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
    std::vector<RegionalExtremaNode>& getExtinctionValues() {
        requireStableTree("ExtinctionValues::getExtinctionValues");
        return regionalExtremaNodes;
    }
};


} // namespace mmcfilters
