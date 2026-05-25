#pragma once

#include "../attributes/AttributeComputation.hpp"
#include "../trees/detail/TreeAltitudeDeltaNeighborhood.hpp"
#include "../trees/MorphologicalTree.hpp"
#include "../trees/TreeAltitudeAlgorithms.hpp"
#include "../trees/WeightedMorphologicalTree.hpp"
#include "../utils/Common.hpp"

#include <cassert>
#include <cmath>
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
 * `delta` units away in altitude space. The node stability is then defined
 * from a monotone increasing attribute as:
 *
 * `stability(node) = (attr(asc(node)) - attr(desc(node))) / attr(node)`
 *
 * A node is marked as MSER-like when:
 * - both delta neighbours exist;
 * - its stability is a strict local minimum with respect to those neighbours;
 * - the stability lies below `maxVariation`;
 * - the attribute value is within the user-specified `[minAttr, maxAttr]`
 *   acceptance interval.
 *
 * If no external attribute buffer is supplied, the class lazily computes
 * `AREA` and uses it as the increasing attribute. The object may therefore act
 * either as a lightweight view over an externally owned buffer or as a small
 * owner of the fallback `AREA` buffer.
 *
 * @note MSER only makes semantic sense when the altitude matches the weighted
 * tree topology. Public constructors therefore accept
 * `WeightedMorphologicalTree<T>` only.
 */
template<AltitudeValue T>
class ComputerMSER {
private:
	const WeightedMorphologicalTree<T>& weighted_;
	const MorphologicalTree& tree;
	const std::vector<T>& altitude_;
	const float* attrMserView_;
	std::vector<float> ownedAttrMser_;
	float maxVariation;
	float minAttr;
	float maxAttr;
	int num;
	std::vector<float> stability;
	std::vector<NodeId> ascendants;
	std::vector<NodeId> descendants;

	static void validateOwnedAttributeSize(const MorphologicalTree& tree, const std::vector<float>& attr) {
		if (attr.size() != static_cast<std::size_t>(tree.getNumInternalNodeSlots())) {
			throw std::invalid_argument("ComputerMSER attribute size must match the internal node slot count.");
		}
	}

	ComputerMSER(const WeightedMorphologicalTree<T>& weighted, const float* attr_increasing, std::vector<float> ownedAttr)
		: weighted_(weighted),
		  tree(weighted.topology()),
		  altitude_(weighted.getAltitudeBuffer()),
		  attrMserView_(nullptr),
		  ownedAttrMser_(std::move(ownedAttr)),
		  maxVariation(10.0),
		  minAttr(0),
		  maxAttr(this->tree.getNumColsOfImage() * this->tree.getNumRowsOfImage()) {
		TreeAltitudeAlgorithms::validateAltitudeBufferShape(this->tree, std::span<const T>(this->altitude_));
		this->attrMserView_ = this->ownedAttrMser_.empty() ? attr_increasing : this->ownedAttrMser_.data();
	}

public:
	/**
	 * @brief Creates an MSER detector backed by an owned attribute buffer.
	 */
	ComputerMSER(const WeightedMorphologicalTree<T>& weighted, std::vector<float> attr_increasing)
		: ComputerMSER(weighted, nullptr, [&]() {
			validateOwnedAttributeSize(weighted.topology(), attr_increasing);
			return std::move(attr_increasing);
		}()) {}

	/**
	 * @brief Creates an MSER detector backed by a non-owning attribute view.
	 *
	 * The raw pointer must reference one value per internal node slot.
	 */
	ComputerMSER(const WeightedMorphologicalTree<T>& weighted, const float* attr_increasing)
		: ComputerMSER(weighted, attr_increasing, {}) {
		if (attr_increasing == nullptr) {
			throw std::invalid_argument("ComputerMSER requires a non-null attribute buffer for the raw-pointer constructor.");
		}
	}

	/**
	 * @brief Creates an MSER detector that lazily falls back to `AREA`.
	 */
	ComputerMSER(const WeightedMorphologicalTree<T>& weighted)
		: ComputerMSER(weighted, nullptr, {}) { }

		~ComputerMSER() = default;

	/**
	 * @brief Computes the MSER indicator vector for the given delta.
	 *
	 * @return A dense boolean-like vector indexed by the tree's internal node
	 * slots, with `true` at the nodes selected as MSER-like regions.
	 */
	[[nodiscard]] std::vector<uint8_t> computeMSER(AltitudeDiff<T> delta){
		std::pair<std::vector<NodeId>, std::vector<NodeId>> ascDesc =
			detail::computeAscendantsAndDescendantsByAltitude(
				tree,
				std::span<const T>(altitude_),
				delta);
		this->ascendants = std::move(ascDesc.first);
		this->descendants = std::move(ascDesc.second);
		this->stability.assign(tree.getNumInternalNodeSlots(), std::numeric_limits<float>::quiet_NaN());

		for (NodeId nodeId : tree.getAliveNodeIds()) {
			const NodeId node = nodeId;
			if(this->ascendants[node] != InvalidNode && this->descendants[node] != InvalidNode){
				this->stability[node] = this->getStability(node);
			}
		}

		this->num = 0;
		double maxStabilityDesc, maxStabilityAsc;
		std::vector<uint8_t> mser(this->tree.getNumInternalNodeSlots(), false);
		for (NodeId nodeId : tree.getAliveNodeIds()) {
			const NodeId node = nodeId;
			if(!std::isnan(this->stability[node]) && !std::isnan(this->stability[this->ascendants[node]]) && !std::isnan(this->stability[this->descendants[node]])){
				maxStabilityDesc = this->stability[this->descendants[node]];
				maxStabilityAsc = this->stability[this->ascendants[node]];
				if(this->stability[node] < maxStabilityDesc && this->stability[node] < maxStabilityAsc){
					if(stability[node] < this->maxVariation && this->getAttrMSER(node) >= this->minAttr && this->getAttrMSER(node) <= this->maxAttr){
						mser[node] = true;
						this->num++;
					}
				}
			}
		}
		return mser;
	}


	/**
	 * @brief Returns the stability score currently associated with a node.
	 */
	[[nodiscard]] double getStability(NodeId node){
		return (this->getAttrMSER(this->ascendants[node]) - this->getAttrMSER(this->descendants[node])) / this->getAttrMSER(node)  ;
	}

	/**
	 * @brief Returns the attribute used by the MSER criterion, lazily
	 * computing `AREA` when no external buffer has been provided.
	 */
	[[nodiscard]] float getAttrMSER(NodeId node){
			if(attrMserView_ == nullptr) {
				auto area = AttributeComputation::computeSingleAttribute(weighted_, AREA);
				ownedAttrMser_ = std::move(area.second);
				attrMserView_ = ownedAttrMser_.data();
			}
		if(attrMserView_ == nullptr)
			return 0.0f;
		else
			return this->attrMserView_[node];
	}

	/**
	 * @brief Returns the most stable node among the current node and its
	 * delta-linked neighbours.
	 */
	[[nodiscard]] NodeId getNodeInPathWithMaxStability(NodeId node){
		NodeId nodeAsc = this->ascendants[node];
		NodeId nodeDes = this->descendants[node];

        if(stability[node] <= stability[nodeDes] && stability[node] <= stability[nodeAsc]) {
            return node;
        }else if (stability[nodeDes] <= stability[nodeAsc]) {
            return nodeDes;
        }else {
            return nodeAsc;
        }

	}


	/**
	 * @brief Returns the ascendant used in the current stability window.
	 */
	[[nodiscard]] NodeId ascendantWithMaxStability(NodeId node) const { return this->ascendants[node];}
	/**
	 * @brief Returns the descendant used in the current stability window.
	 */
	[[nodiscard]] NodeId descendantWithMaxStability(NodeId node) const { return descendants[node];}
	/**
	 * @brief Returns the current stability array, indexed by node slot.
	 */
	std::vector<float>& getStabilities() { return stability; }
	/**
	 * @brief Returns the number of nodes selected as MSER-like in the last run.
	 */
	[[nodiscard]] int getNumNodes() {return  num;}
	/**
	 * @brief Sets the maximum accepted stability value.
	 */
	void setMaxVariation(float maxVariation) { this->maxVariation = maxVariation; }
	/**
	 * @brief Sets the lower bound of the accepted attribute interval.
	 */
	void setMinAttribute(float minAttr) { this->minAttr = minAttr; }
	/**
	 * @brief Sets the upper bound of the accepted attribute interval.
	 */
	void setMaxAttribute(float maxAttr) { this->maxAttr = maxAttr; }
};


} // namespace mmcfilters
