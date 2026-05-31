#pragma once

#include "../attributes/AttributeComputation.hpp"
#include "../trees/detail/TreeStabilityNeighborhood.hpp"
#include "../trees/detail/TreeKindValidation.hpp"
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
 * @note MSER only makes semantic sense for max-trees and min-trees whose
 * altitude matches the weighted tree topology. Public constructors therefore
 * accept `WeightedMorphologicalTree<T>` only and reject self-dual tree kinds.
 */
template<AltitudeValue T, std::floating_point Real = float>
class MSERComputer {
public:
	/// Floating-point type used to store variation scores.
	using variation_value_type = Real;

private:
	const WeightedMorphologicalTree<T>& weighted_;
	const MorphologicalTree& tree;
	const std::vector<T>& altitude_;
	const Real* attrMserView_;
	std::vector<Real> ownedAttrMser_;
	Real maxVariation;
	Real minAttr;
	Real maxAttr;
	int num;
	std::vector<Real> variation;
	std::vector<NodeId> ascendants;
	std::vector<NodeId> descendants;

	static void validateOwnedAttributeSize(const MorphologicalTree& tree, const std::vector<Real>& attr) {
		if (attr.size() != static_cast<std::size_t>(tree.getNumInternalNodeSlots())) {
			throw std::invalid_argument("MSERComputer attribute size must match the internal node slot count.");
		}
	}

	MSERComputer(const WeightedMorphologicalTree<T>& weighted, const Real* attr_increasing, std::vector<Real> ownedAttr)
		: weighted_(weighted),
		  tree(weighted.topology()),
		  altitude_(weighted.getAltitudeBuffer()),
		  attrMserView_(nullptr),
		  ownedAttrMser_(std::move(ownedAttr)),
		  maxVariation(Real{10}),
		  minAttr(Real{0}),
		  maxAttr(static_cast<Real>(this->tree.getNumColsOfImage() * this->tree.getNumRowsOfImage())) {
		TreeAltitudeAlgorithms::validateAltitudeBufferShape(this->tree, std::span<const T>(this->altitude_));
		detail::validateComponentTreeKind(this->tree, "MSERComputer");
		this->attrMserView_ = this->ownedAttrMser_.empty() ? attr_increasing : this->ownedAttrMser_.data();
	}

public:
	/**
	 * @brief Creates an MSER detector backed by an owned attribute buffer.
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
	 */
	MSERComputer(const WeightedMorphologicalTree<T>& weighted, const Real* attr_increasing)
		: MSERComputer(weighted, attr_increasing, {}) {
		if (attr_increasing == nullptr) {
			throw std::invalid_argument("MSERComputer requires a non-null attribute buffer for the raw-pointer constructor.");
		}
	}

	/**
	 * @brief Creates an MSER detector that lazily falls back to `AREA`.
	 */
	MSERComputer(const WeightedMorphologicalTree<T>& weighted)
		: MSERComputer(weighted, nullptr, {}) { }

	~MSERComputer() = default;

	/**
	 * @brief Computes the MSER indicator vector for the given delta.
	 *
	 * @return A dense boolean-like vector indexed by the tree's internal node
	 * slots, with `true` at the nodes selected as MSER-like regions.
	 */
	[[nodiscard]] std::vector<uint8_t> computeMSER(AltitudeDiff<T> delta){
		detail::StabilityNeighborhood neighborhood =
			detail::computeAltitudeStabilityNeighborhood(
				tree,
				std::span<const T>(altitude_),
				delta);
		this->ascendants = std::move(neighborhood.ascendants);
		this->descendants = std::move(neighborhood.descendants);

		auto attrAt = [this](NodeId node) -> Real {
			return this->getAttrMSER(node);
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
	[[nodiscard]] Real getVariation(NodeId node){
		detail::validateStabilityNeighborhoodShape(
			this->tree,
			this->ascendants,
			this->descendants,
			"MSERComputer::getVariation");
		if (this->variation.size() != static_cast<std::size_t>(this->tree.getNumInternalNodeSlots())) {
			throw std::logic_error("MSERComputer::getVariation requires computeMSER to run first.");
		}
		auto attrAt = [this](NodeId nodeId) -> Real {
			return this->getAttrMSER(nodeId);
		};
		return detail::computeVariationValue<Real>(
			node,
			this->ascendants,
			this->descendants,
			attrAt);
	}

	/**
	 * @brief Returns the attribute used by the MSER criterion, lazily
	 * computing `AREA` when no external buffer has been provided.
	 */
	[[nodiscard]] Real getAttrMSER(NodeId node){
		if(attrMserView_ == nullptr) {
			auto area = AttributeComputation::computeSingleAttribute<Real>(weighted_, AREA);
			ownedAttrMser_ = std::move(area.second);
			attrMserView_ = ownedAttrMser_.data();
		}
		if(attrMserView_ == nullptr)
			return Real{0};
		else
			return this->attrMserView_[node];
	}

	/**
	 * @brief Returns the node with minimum variation among the current node and
	 * its delta-linked neighbours.
	 */
	[[nodiscard]] NodeId nodeWithMinimumVariationInWindow(NodeId node){
		detail::validateStabilityNeighborhoodShape(
			this->tree,
			this->ascendants,
			this->descendants,
			"MSERComputer::nodeWithMinimumVariationInWindow");
		if (this->variation.size() != static_cast<std::size_t>(this->tree.getNumInternalNodeSlots())) {
			throw std::logic_error("MSERComputer::nodeWithMinimumVariationInWindow requires computeMSER to run first.");
		}
		return detail::nodeWithMinimumVariationInWindow<Real>(
			node,
			this->variation,
			this->ascendants,
			this->descendants);
	}

	/**
	 * @brief Returns the ascendant used in the current stability window.
	 */
	[[nodiscard]] NodeId ascendantInStabilityWindow(NodeId node) const { return this->ascendants[node];}
	/**
	 * @brief Returns the descendant used in the current stability window.
	 */
	[[nodiscard]] NodeId descendantInStabilityWindow(NodeId node) const { return descendants[node];}
	/**
	 * @brief Returns the current variation array, indexed by node slot.
	 */
	std::vector<Real>& getVariations() { return variation; }
	/**
	 * @brief Returns the number of nodes selected as MSER-like in the last run.
	 */
	[[nodiscard]] int getNumNodes() {return  num;}
	/**
	 * @brief Sets the maximum accepted variation value.
	 */
	void setMaxVariation(Real maxVariation) { this->maxVariation = maxVariation; }
	/**
	 * @brief Sets the lower bound of the accepted attribute interval.
	 */
	void setMinAttribute(Real minAttr) { this->minAttr = minAttr; }
	/**
	 * @brief Sets the upper bound of the accepted attribute interval.
	 */
	void setMaxAttribute(Real maxAttr) { this->maxAttr = maxAttr; }
};


} // namespace mmcfilters
