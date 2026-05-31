#pragma once

#include "../attributes/AttributeComputation.hpp"
#include "../trees/detail/TreeStabilityNeighborhood.hpp"
#include "../trees/MorphologicalTree.hpp"
#include "../utils/Common.hpp"
#include "detail/VariationMeasure.hpp"

#include <concepts>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcfilters {

/**
 * @brief Detects stable regions from a topological depth window.
 *
 * @details
 * This class is the self-dual/topological counterpart of `MSERComputer`.
 * `depthDelta` is a number of tree edges, not an altitude distance:
 *
 * - the ascendant of a node is obtained by climbing exactly `depthDelta`
 *   parent links;
 * - the descendant is selected among nodes exactly `depthDelta` child links
 *   below the node;
 * - when several descendants exist, the largest-area one is selected, with ties
 *   resolved by the smallest `NodeId`;
 * - altitude is never read.
 *
 * The variation value uses the same area-variation form as MSER:
 *
 * `variation(node) = (attr(asc(node)) - attr(desc(node))) / attr(node)`
 *
 * and selected nodes are strict local minima of this variation. The default
 * attribute is `AREA`, computed from topology only.
 */
template<std::floating_point Real = float>
class DepthStableRegionComputer {
public:
    /// Floating-point type used to store variation scores.
    using variation_value_type = Real;

private:
    const MorphologicalTree& tree;
    const Real* attrView_ = nullptr;
    std::vector<Real> ownedAttr_;
    Real maxVariation = Real{10};
    Real minAttr = Real{0};
    Real maxAttr;
    int num = 0;
    std::vector<Real> variation;
    std::vector<NodeId> ascendants;
    std::vector<NodeId> descendants;

    static void validateOwnedAttributeSize(const MorphologicalTree& tree, const std::vector<Real>& attr) {
        if (attr.size() != static_cast<std::size_t>(tree.getNumInternalNodeSlots())) {
            throw std::invalid_argument("DepthStableRegionComputer attribute size must match the internal node slot count.");
        }
    }

    DepthStableRegionComputer(const MorphologicalTree& tree, const Real* attr, std::vector<Real> ownedAttr)
        : tree(tree),
          attrView_(nullptr),
          ownedAttr_(std::move(ownedAttr)),
          maxAttr(static_cast<Real>(tree.getNumColsOfImage() * tree.getNumRowsOfImage())) {
        this->attrView_ = this->ownedAttr_.empty() ? attr : this->ownedAttr_.data();
    }

public:
    /**
     * @brief Creates a detector backed by an owned increasing-attribute buffer.
     */
    DepthStableRegionComputer(const MorphologicalTree& tree, std::vector<Real> attr)
        : DepthStableRegionComputer(tree, nullptr, [&]() {
            validateOwnedAttributeSize(tree, attr);
            return std::move(attr);
        }()) {}

    /**
     * @brief Creates a detector backed by a non-owning attribute buffer.
     */
    DepthStableRegionComputer(const MorphologicalTree& tree, const Real* attr)
        : DepthStableRegionComputer(tree, attr, {}) {
        if (attr == nullptr) {
            throw std::invalid_argument("DepthStableRegionComputer requires a non-null attribute buffer for the raw-pointer constructor.");
        }
    }

    /**
     * @brief Creates a detector that lazily computes topology-only `AREA`.
     */
    explicit DepthStableRegionComputer(const MorphologicalTree& tree)
        : DepthStableRegionComputer(tree, nullptr, {}) {}

    /**
     * @brief Computes the stable-region indicator vector for a positive depth delta.
     */
    [[nodiscard]] std::vector<uint8_t> computeByDepth(int depthDelta) {
        detail::StabilityNeighborhood neighborhood =
            detail::computeDepthStabilityNeighborhood(this->tree, depthDelta);
        this->ascendants = std::move(neighborhood.ascendants);
        this->descendants = std::move(neighborhood.descendants);

        auto attrAt = [this](NodeId node) -> Real {
            return this->getAttribute(node);
        };
        this->variation =
            detail::computeVariationsFromNeighborhood<Real>(
                this->tree,
                this->ascendants,
                this->descendants,
                attrAt);
        return detail::selectStrictVariationMinima<Real>(
            this->tree,
            this->variation,
            this->ascendants,
            this->descendants,
            attrAt,
            this->maxVariation,
            this->minAttr,
            this->maxAttr,
            this->num);
    }

    /**
     * @brief Returns the variation score currently associated with a node.
     */
    [[nodiscard]] Real getVariation(NodeId node) {
        detail::validateStabilityNeighborhoodShape(
            this->tree,
            this->ascendants,
            this->descendants,
            "DepthStableRegionComputer::getVariation");
        if (this->variation.size() != static_cast<std::size_t>(this->tree.getNumInternalNodeSlots())) {
            throw std::logic_error("DepthStableRegionComputer::getVariation requires computeByDepth to run first.");
        }
        auto attrAt = [this](NodeId nodeId) -> Real {
            return this->getAttribute(nodeId);
        };
        return detail::computeVariationValue<Real>(
            node,
            this->ascendants,
            this->descendants,
            attrAt);
    }

    /**
     * @brief Returns the increasing attribute used by the criterion.
     */
    [[nodiscard]] Real getAttribute(NodeId node) {
        if (this->attrView_ == nullptr) {
            auto area = AttributeComputation::computeSingleTopologyAttribute<Real>(this->tree, AREA);
            this->ownedAttr_ = std::move(area.second);
            this->attrView_ = this->ownedAttr_.data();
        }
        return this->attrView_[static_cast<std::size_t>(node)];
    }

    /**
     * @brief Returns the node with minimum variation in the depth window.
     */
    [[nodiscard]] NodeId nodeWithMinimumVariationInWindow(NodeId node) const {
        detail::validateStabilityNeighborhoodShape(
            this->tree,
            this->ascendants,
            this->descendants,
            "DepthStableRegionComputer::nodeWithMinimumVariationInWindow");
        if (this->variation.size() != static_cast<std::size_t>(this->tree.getNumInternalNodeSlots())) {
            throw std::logic_error("DepthStableRegionComputer::nodeWithMinimumVariationInWindow requires computeByDepth to run first.");
        }
        return detail::nodeWithMinimumVariationInWindow<Real>(
            node,
            this->variation,
            this->ascendants,
            this->descendants);
    }

    /**
     * @brief Returns the ascendant used in the current depth window.
     */
    [[nodiscard]] NodeId ascendantInStabilityWindow(NodeId node) const { return this->ascendants[static_cast<std::size_t>(node)]; }

    /**
     * @brief Returns the descendant used in the current depth window.
     */
    [[nodiscard]] NodeId descendantInStabilityWindow(NodeId node) const { return this->descendants[static_cast<std::size_t>(node)]; }

    /**
     * @brief Returns the current variation array, indexed by node slot.
     */
    std::vector<Real>& getVariations() { return variation; }

    /**
     * @brief Returns the number of nodes selected in the last run.
     */
    [[nodiscard]] int getNumNodes() const { return num; }

    /**
     * @brief Sets the maximum accepted variation value.
     */
    void setMaxVariation(Real value) { this->maxVariation = value; }

    /**
     * @brief Sets the lower bound of the accepted attribute interval.
     */
    void setMinAttribute(Real value) { this->minAttr = value; }

    /**
     * @brief Sets the upper bound of the accepted attribute interval.
     */
    void setMaxAttribute(Real value) { this->maxAttr = value; }
};

} // namespace mmcfilters
