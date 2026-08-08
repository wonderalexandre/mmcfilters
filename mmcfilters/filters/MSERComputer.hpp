#pragma once

#include "../attributes/AttributeComputation.hpp"
#include "../trees/detail/TreeStabilityNeighborhood.hpp"
#include "../trees/detail/HierarchyCapabilityValidation.hpp"
#include "../trees/MorphologicalTree.hpp"
#include "../trees/TreeAltitudeAlgorithms.hpp"
#include "../trees/WeightedMorphologicalTree.hpp"
#include "../utils/Common.hpp"
#include "detail/VariationMeasure.hpp"

#include <cassert>
#include <cmath>
#include <concepts>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mmcfilters {

/**
 * @brief Detects MSER-like nodes from a monotone increasing attribute defined
 * on the hierarchy.
 *
 * @details
 * The class implements the classical MSER stability criterion in tree form.
 * It is intentionally tied to `WeightedMorphologicalTree<T>`: the delta
 * neighbourhood is computed from the altitude buffer owned by the weighted
 * tree, not from an arbitrary external altitude array. Given a delta value,
 * each node is paired with an ascendant and a descendant located approximately
 * `delta` units away in altitude space. The node variation is then defined
 * from a monotone increasing attribute as:
 *
 * `variation(node) = (attr(asc(node)) - attr(desc(node))) / attr(node)`
 *
 * A node is marked as MSER-like when:
 * @li both delta neighbours exist;
 * @li its variation is a strict local minimum with respect to those neighbours;
 * @li the variation lies below `maxVariation`;
 * @li the attribute value is within the user-specified `[minAttr, maxAttr]`
 *     acceptance interval.
 *
 * If no external attribute buffer is supplied, the class lazily computes
 * `AREA` and uses it as the increasing attribute. The object may therefore act
 * either as a lightweight view over an externally owned buffer or as a small
 * owner of the fallback `AREA` buffer.
 *
 * @note MSER requires a globally monotone altitude capability on the topology.
 * Standard max-tree and min-tree producers provide that capability. Standard
 * tree-of-shapes and self-dual residual-tree producers publish
 * `AltitudeOrder::UNCONSTRAINED` and are therefore rejected. The descriptive
 * tree kind itself is not used as the acceptance gate.
 */
template <AltitudeValue T, std::floating_point Real = float> class MSERComputer {
  public:
    /// Floating-point type used to store variation scores.
    using variation_value_type = Real;

  private:
    /** @brief Stores the weighted. */
    const WeightedMorphologicalTree<T>& weighted_;
    /** @brief References the tree used by the component. */
    const MorphologicalTree& tree;
    /** @brief Stores the altitude. */
    const std::vector<T>& altitude_;
    /** @brief Stores the external attr mser. */
    const Real* externalAttrMser_;
    /** @brief Stores the owned attr mser. */
    std::vector<Real> ownedAttrMser_;
    /** @brief Stores the max variation. */
    Real maxVariation;
    /** @brief Stores the min attr. */
    Real minAttr;
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
    [[nodiscard]] const Real* attributeData() const noexcept { return ownedAttrMser_.empty() ? externalAttrMser_ : ownedAttrMser_.data(); }

    /**
     * @brief Rejects access to results before the first successful computation.
     *
     * @param context Public operation requesting computed state.
     */
    void requireComputed(const char* context) const {
        if (!hasComputed_) {
            throw std::logic_error(std::string(context) + " requires computeMSER to run first.");
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
            throw std::invalid_argument("MSERComputer attribute size must match the internal node slot count.");
        }
    }

    /**
     * @brief Creates an MSER evaluator over a weighted tree and attribute buffer.
     *
     * @param weighted Weighted tree used by the operation.
     * @param attr_increasing Attribute information represented by `attr_increasing`.
     * @param ownedAttr Owned attribute storage used when no external buffer is supplied.
     */
    MSERComputer(const WeightedMorphologicalTree<T>& weighted, const Real* attr_increasing, std::vector<Real> ownedAttr)
        : weighted_(weighted), tree(weighted.topology()), altitude_(weighted.getAltitudeBuffer()), externalAttrMser_(attr_increasing),
          ownedAttrMser_(std::move(ownedAttr)), maxVariation(Real{10}), minAttr(Real{0}),
          maxAttr(static_cast<Real>(this->tree.getNumColsOfGridDomain2D() * this->tree.getNumRowsOfGridDomain2D())) {
        TreeAltitudeAlgorithms::validateAltitudeBufferShape(this->tree, std::span<const T>(this->altitude_));
        detail::validateGlobalMonotoneAltitudeOrder(this->tree, "MSERComputer");
    }

  public:
    /**
     * @brief Creates an MSER detector backed by an owned attribute buffer.
     *
     * @param weighted Weighted tree used by the operation.
     * @param attr_increasing Attribute information represented by `attr_increasing`.
     */
    MSERComputer(const WeightedMorphologicalTree<T>& weighted, std::vector<Real> attr_increasing)
        : MSERComputer(weighted, nullptr, [&]() {
              validateOwnedAttributeSize(weighted.topology(), attr_increasing);
              return std::move(attr_increasing);
          }()) {}

    /**
     * @brief Creates an MSER detector backed by a non-owning attribute view.
     *
     * The raw pointer must reference one value per internal node slot.
     *
     * @param weighted Weighted tree used by the operation.
     * @param attr_increasing Attribute information represented by `attr_increasing`.
     */
    MSERComputer(const WeightedMorphologicalTree<T>& weighted, const Real* attr_increasing) : MSERComputer(weighted, attr_increasing, {}) {
        if (attr_increasing == nullptr) {
            throw std::invalid_argument("MSERComputer requires a non-null attribute buffer for the raw-pointer constructor.");
        }
    }

    /**
     * @brief Creates an MSER detector that lazily falls back to `AREA`.
     *
     * @param weighted Weighted tree used by the operation.
     */
    MSERComputer(const WeightedMorphologicalTree<T>& weighted) : MSERComputer(weighted, nullptr, {}) {}

    /** @brief Copies the evaluator while preserving owned-buffer safety. */
    MSERComputer(const MSERComputer&) = default;
    /** @brief Moves the evaluator while preserving its referenced tree. */
    MSERComputer(MSERComputer&&) noexcept = default;
    /** @brief Assignment is unavailable because the evaluator stores tree references. */
    MSERComputer& operator=(const MSERComputer&) = delete;
    /** @brief Move assignment is unavailable because the evaluator stores tree references. */
    MSERComputer& operator=(MSERComputer&&) = delete;

    /**
     * @brief Destroys the MSER evaluator.
     */
    ~MSERComputer() = default;

    /**
     * @brief Computes the MSER indicator vector for the given delta.
     *
     * @param delta Delta offset or radius used by the operation.
     * @return A dense boolean-like vector indexed by the tree's internal node
     * slots, with `true` at the nodes selected as MSER-like regions.
     *
     */
    [[nodiscard]] std::vector<uint8_t> computeMSER(AltitudeDiff<T> delta) {
        detail::StabilityNeighborhood neighborhood = detail::computeAltitudeStabilityNeighborhood(tree, std::span<const T>(altitude_), delta);
        this->ascendants = std::move(neighborhood.ascendants);
        this->descendants = std::move(neighborhood.descendants);

        auto attrAt = [this](NodeId node) -> Real { return this->getAttrMSER(node); };
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
        requireComputed("MSERComputer::getVariation");
        detail::validateStabilityNeighborhoodShape(this->tree, this->ascendants, this->descendants, "MSERComputer::getVariation");
        auto attrAt = [this](NodeId nodeId) -> Real { return this->getAttrMSER(nodeId); };
        return detail::computeVariationValue<Real>(node, this->ascendants, this->descendants, attrAt);
    }

    /**
     * @brief Returns the attribute used by the MSER criterion, lazily
     * computing `AREA` when no external buffer has been provided.
     *
     * @param node Node identifier used by the operation.
     * @return The attribute used by the MSER criterion, lazily computing AREA when no external buffer has been provided.
     */
    [[nodiscard]] Real getAttrMSER(NodeId node) {
        if (attributeData() == nullptr) {
            auto area = AttributeComputation::computeSingleAttribute<Real>(weighted_, AREA);
            ownedAttrMser_ = std::move(area.second);
        }
        const Real* data = attributeData();
        if (data == nullptr) {
            throw std::logic_error("MSERComputer attribute storage is unavailable.");
        }
        return data[node];
    }

    /**
     * @brief Returns the node with minimum variation among the current node and
     * its delta-linked neighbours.
     *
     * @param node Node identifier used by the operation.
     * @return The node with minimum variation among the current node and its delta-linked neighbours.
     */
    [[nodiscard]] NodeId nodeWithMinimumVariationInWindow(NodeId node) {
        requireComputed("MSERComputer::nodeWithMinimumVariationInWindow");
        detail::validateStabilityNeighborhoodShape(this->tree, this->ascendants, this->descendants, "MSERComputer::nodeWithMinimumVariationInWindow");
        return detail::nodeWithMinimumVariationInWindow<Real>(node, this->variation, this->ascendants, this->descendants);
    }

    /**
     * @brief Returns the ascendant used in the current stability window.
     *
     * @param node Node identifier used by the operation.
     * @return The ascendant used in the current stability window.
     */
    [[nodiscard]] NodeId ascendantInStabilityWindow(NodeId node) const {
        requireComputed("MSERComputer::ascendantInStabilityWindow");
        return this->ascendants[node];
    }

    /**
     * @brief Returns the descendant used in the current stability window.
     *
     * @param node Node identifier used by the operation.
     * @return The descendant used in the current stability window.
     */
    [[nodiscard]] NodeId descendantInStabilityWindow(NodeId node) const {
        requireComputed("MSERComputer::descendantInStabilityWindow");
        return this->descendants[node];
    }

    /**
     * @brief Returns the current variation array, indexed by node slot.
     *
     * @return The current variation array, indexed by node slot.
     */
    [[nodiscard]] std::vector<Real>& getVariations() {
        requireComputed("MSERComputer::getVariations");
        return this->variation;
    }

    /**
     * @brief Returns the number of nodes selected as MSER-like in the last run.
     *
     * @return The number of nodes selected as MSER-like in the last run.
     */
    [[nodiscard]] int getNumNodes() {
        requireComputed("MSERComputer::getNumNodes");
        return this->num;
    }

    /**
     * @brief Sets the maximum accepted variation value.
     *
     * @param maxVariation Maximum accepted stability variation.
     */
    void setMaxVariation(Real maxVariation) { this->maxVariation = maxVariation; }
    /**
     * @brief Sets the lower bound of the accepted attribute interval.
     *
     * @param minAttr Minimum accepted attribute value.
     */
    void setMinAttribute(Real minAttr) { this->minAttr = minAttr; }
    /**
     * @brief Sets the upper bound of the accepted attribute interval.
     *
     * @param maxAttr Maximum accepted attribute value.
     */
    void setMaxAttribute(Real maxAttr) { this->maxAttr = maxAttr; }
};

} // namespace mmcfilters
