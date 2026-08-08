#pragma once

#include "../attributes/AttributeComputation.hpp"
#include "../trees/detail/TreeStabilityNeighborhood.hpp"
#include "../trees/MorphologicalTree.hpp"
#include "../utils/Common.hpp"
#include "detail/VariationMeasure.hpp"

#include <concepts>
#include <limits>
#include <stdexcept>
#include <string>
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
template <std::floating_point Real = float> class DepthStableRegionComputer {
  public:
    /// Floating-point type used to store variation scores.
    using variation_value_type = Real;

  private:
    /** @brief References the tree used by the component. */
    const MorphologicalTree& tree;
    /** @brief Stores the external attr. */
    const Real* externalAttr_ = nullptr;
    /** @brief Stores the owned attr. */
    std::vector<Real> ownedAttr_;
    /** @brief Stores the max variation. */
    Real maxVariation = Real{10};
    /** @brief Stores the min attr. */
    Real minAttr = Real{0};
    /** @brief Stores the max attr. */
    Real maxAttr;
    /** @brief Stores the num. */
    int num = 0;
    /** @brief Indicates whether a successful computation result is available. */
    bool hasComputed_ = false;
    /** @brief Stores the variation. */
    std::vector<Real> variation;
    /** @brief Stores the ascendants. */
    std::vector<NodeId> ascendants;
    /** @brief Stores the descendants. */
    std::vector<NodeId> descendants;

    /**
     * @brief Returns the active owned or borrowed attribute buffer.
     *
     * @return The active buffer, or `nullptr` when lazy `AREA` computation is pending.
     */
    [[nodiscard]] const Real* attributeData() const noexcept { return ownedAttr_.empty() ? externalAttr_ : ownedAttr_.data(); }

    /**
     * @brief Rejects access to results before the first successful computation.
     *
     * @param context Public operation requesting computed state.
     */
    void requireComputed(const char* context) const {
        if (!hasComputed_) {
            throw std::logic_error(std::string(context) + " requires computeByDepth to run first.");
        }
    }

    /**
     * @brief Validates owned attribute size.
     *
     * @param tree Tree topology used by the operation.
     * @param attr Attribute requested by the operation.
     */
    static void validateOwnedAttributeSize(const MorphologicalTree& tree, const std::vector<Real>& attr) {
        if (attr.size() != static_cast<std::size_t>(tree.getNumInternalNodeSlots())) {
            throw std::invalid_argument("DepthStableRegionComputer attribute size must match the internal node slot count.");
        }
    }

    /**
     * @brief Creates a depth-stability evaluator over a tree and attribute buffer.
     *
     * @param tree Tree topology used by the operation.
     * @param attr Attribute requested by the operation.
     * @param ownedAttr Owned attribute storage used when no external buffer is supplied.
     */
    DepthStableRegionComputer(const MorphologicalTree& tree, const Real* attr, std::vector<Real> ownedAttr)
        : tree(tree), externalAttr_(attr), ownedAttr_(std::move(ownedAttr)),
          maxAttr(static_cast<Real>(tree.getNumColsOfGridDomain2D() * tree.getNumRowsOfGridDomain2D())) {}

  public:
    /**
     * @brief Creates a detector backed by an owned increasing-attribute buffer.
     *
     * @param tree Tree topology used by the operation.
     * @param attr Attribute requested by the operation.
     */
    DepthStableRegionComputer(const MorphologicalTree& tree, std::vector<Real> attr)
        : DepthStableRegionComputer(tree, nullptr, [&]() {
              validateOwnedAttributeSize(tree, attr);
              return std::move(attr);
          }()) {}

    /**
     * @brief Creates a detector backed by a non-owning attribute buffer.
     *
     * @param tree Tree topology used by the operation.
     * @param attr Attribute requested by the operation.
     */
    DepthStableRegionComputer(const MorphologicalTree& tree, const Real* attr) : DepthStableRegionComputer(tree, attr, {}) {
        if (attr == nullptr) {
            throw std::invalid_argument("DepthStableRegionComputer requires a non-null attribute buffer for the raw-pointer constructor.");
        }
    }

    /**
     * @brief Creates a detector that lazily computes topology-only `AREA`.
     *
     * @param tree Tree topology used by the operation.
     */
    explicit DepthStableRegionComputer(const MorphologicalTree& tree) : DepthStableRegionComputer(tree, nullptr, {}) {}

    /** @brief Copies the evaluator while preserving owned-buffer safety. */
    DepthStableRegionComputer(const DepthStableRegionComputer&) = default;
    /** @brief Moves the evaluator while preserving its referenced tree. */
    DepthStableRegionComputer(DepthStableRegionComputer&&) noexcept = default;
    /** @brief Assignment is unavailable because the evaluator stores a tree reference. */
    DepthStableRegionComputer& operator=(const DepthStableRegionComputer&) = delete;
    /** @brief Move assignment is unavailable because the evaluator stores a tree reference. */
    DepthStableRegionComputer& operator=(DepthStableRegionComputer&&) = delete;

    /**
     * @brief Computes the stable-region indicator vector for a positive depth delta.
     *
     * @param depthDelta Topological depth offset used by the operation.
     * @return The computed stable-region indicator vector for a positive depth delta.
     */
    [[nodiscard]] std::vector<uint8_t> computeByDepth(int depthDelta) {
        detail::StabilityNeighborhood neighborhood = detail::computeDepthStabilityNeighborhood(this->tree, depthDelta);
        this->ascendants = std::move(neighborhood.ascendants);
        this->descendants = std::move(neighborhood.descendants);

        auto attrAt = [this](NodeId node) -> Real { return this->getAttribute(node); };
        this->variation = detail::computeVariationsFromNeighborhood<Real>(this->tree, this->ascendants, this->descendants, attrAt);
        std::vector<uint8_t> selected = detail::selectStrictVariationMinima<Real>(this->tree, this->variation, this->ascendants, this->descendants, attrAt,
                                                                                  this->maxVariation, this->minAttr, this->maxAttr, this->num);
        this->hasComputed_ = true;
        return selected;
    }

    /**
     * @brief Returns the variation score currently associated with a node.
     *
     * @param node Node identifier used by the operation.
     * @return The variation score currently associated with a node.
     */
    [[nodiscard]] Real getVariation(NodeId node) {
        requireComputed("DepthStableRegionComputer::getVariation");
        detail::validateStabilityNeighborhoodShape(this->tree, this->ascendants, this->descendants, "DepthStableRegionComputer::getVariation");
        auto attrAt = [this](NodeId nodeId) -> Real { return this->getAttribute(nodeId); };
        return detail::computeVariationValue<Real>(node, this->ascendants, this->descendants, attrAt);
    }

    /**
     * @brief Returns the increasing attribute used by the criterion.
     *
     * @param node Node identifier used by the operation.
     * @return The increasing attribute used by the criterion.
     */
    [[nodiscard]] Real getAttribute(NodeId node) {
        if (attributeData() == nullptr) {
            auto area = AttributeComputation::computeSingleTopologyAttribute<Real>(this->tree, AREA);
            this->ownedAttr_ = std::move(area.second);
        }
        const Real* data = attributeData();
        if (data == nullptr) {
            throw std::logic_error("DepthStableRegionComputer attribute storage is unavailable.");
        }
        return data[static_cast<std::size_t>(node)];
    }

    /**
     * @brief Returns the node with minimum variation in the depth window.
     *
     * @param node Node identifier used by the operation.
     * @return The node with minimum variation in the depth window.
     */
    [[nodiscard]] NodeId nodeWithMinimumVariationInWindow(NodeId node) const {
        requireComputed("DepthStableRegionComputer::nodeWithMinimumVariationInWindow");
        detail::validateStabilityNeighborhoodShape(this->tree, this->ascendants, this->descendants,
                                                   "DepthStableRegionComputer::nodeWithMinimumVariationInWindow");
        return detail::nodeWithMinimumVariationInWindow<Real>(node, this->variation, this->ascendants, this->descendants);
    }

    /**
     * @brief Returns the ascendant used in the current depth window.
     *
     * @param node Node identifier used by the operation.
     * @return The ascendant used in the current depth window.
     */
    [[nodiscard]] NodeId ascendantInStabilityWindow(NodeId node) const {
        requireComputed("DepthStableRegionComputer::ascendantInStabilityWindow");
        return this->ascendants[static_cast<std::size_t>(node)];
    }

    /**
     * @brief Returns the descendant used in the current depth window.
     *
     * @param node Node identifier used by the operation.
     * @return The descendant used in the current depth window.
     */
    [[nodiscard]] NodeId descendantInStabilityWindow(NodeId node) const {
        requireComputed("DepthStableRegionComputer::descendantInStabilityWindow");
        return this->descendants[static_cast<std::size_t>(node)];
    }

    /**
     * @brief Returns the current variation array, indexed by node slot.
     *
     * @return The current variation array, indexed by node slot.
     */
    [[nodiscard]] std::vector<Real>& getVariations() {
        requireComputed("DepthStableRegionComputer::getVariations");
        return variation;
    }

    /**
     * @brief Returns the number of nodes selected in the last run.
     *
     * @return The number of nodes selected in the last run.
     */
    [[nodiscard]] int getNumNodes() const {
        requireComputed("DepthStableRegionComputer::getNumNodes");
        return num;
    }

    /**
     * @brief Sets the maximum accepted variation value.
     *
     * @param value Value used by the operation.
     */
    void setMaxVariation(Real value) { this->maxVariation = value; }

    /**
     * @brief Sets the lower bound of the accepted attribute interval.
     *
     * @param value Value used by the operation.
     */
    void setMinAttribute(Real value) { this->minAttr = value; }

    /**
     * @brief Sets the upper bound of the accepted attribute interval.
     *
     * @param value Value used by the operation.
     */
    void setMaxAttribute(Real value) { this->maxAttr = value; }
};

} // namespace mmcfilters
