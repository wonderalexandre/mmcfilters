#pragma once

#include "../attributes/AttributeComputedIncrementally.hpp"
#include "../trees/MorphologicalTree.hpp"
#include "../trees/WeightedMorphologicalTree.hpp"
#include "../utils/Common.hpp"

#include <cassert>
#include <stdexcept>

namespace mmcfilters {

/**
 * @brief Detects MSER-like nodes from a monotone increasing attribute defined
 * on the hierarchy.
 *
 * @details
 * The class implements the classical MSER stability criterion in tree form.
 * Given a node altitude buffer and a delta value, each node is paired with an
 * ascendant and a descendant located approximately `delta` units away in
 * altitude space. The node stability is then defined from a monotone
 * increasing attribute as:
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
 * @note The caller is responsible for providing explicit node altitudes and an
 * attribute that is meaningful for MSER-like selection. The implementation uses
 * altitude for the delta neighbourhood and the attribute for the stability
 * score.
 */
class ComputerMSER {
private:
	const MorphologicalTree& tree;
	const AltitudeBuffer* altitude_;
	const float* attrMserView_;
	std::vector<float> ownedAttrMser_;
	float maxVariation;
	float minAttr;
	float maxAttr;
	int num;
	std::vector<float> stability;
	std::vector<NodeId> ascendants;
	std::vector<NodeId> descendants;

	
public:
	/**
	 * @brief Creates an MSER detector backed by an owned attribute buffer.
	 */
	ComputerMSER(const MorphologicalTree& tree, const AltitudeBuffer* altitude, std::vector<float> attr_increasing)
		: tree(tree),
		  altitude_(&WeightedMorphologicalTree::requireAltitudeBuffer(altitude)),
		  attrMserView_(nullptr),
		  ownedAttrMser_(std::move(attr_increasing)),
		  maxVariation(10.0),
		  minAttr(0),
		  maxAttr(this->tree.getNumColsOfImage() * this->tree.getNumRowsOfImage()) {
		WeightedMorphologicalTree::validateAltitudeBufferShape(this->tree, this->altitude_);
		this->attrMserView_ = this->ownedAttrMser_.data();
	}

	/**
	 * @brief Creates an MSER detector backed by a non-owning attribute view.
	 */
	ComputerMSER(const MorphologicalTree& tree, const AltitudeBuffer* altitude, const float* attr_increasing)
		: tree(tree),
		  altitude_(&WeightedMorphologicalTree::requireAltitudeBuffer(altitude)),
		  attrMserView_(attr_increasing),
		  ownedAttrMser_(),
		  maxVariation(10.0),
		  minAttr(0),
		  maxAttr(this->tree.getNumColsOfImage() * this->tree.getNumRowsOfImage()) {
		WeightedMorphologicalTree::validateAltitudeBufferShape(this->tree, this->altitude_);
	}

	ComputerMSER(const WeightedMorphologicalTree& weighted, std::vector<float> attr_increasing)
		: ComputerMSER(weighted.tree_, &weighted.altitude_, std::move(attr_increasing)) {}

	ComputerMSER(const WeightedMorphologicalTree& weighted, const float* attr_increasing)
		: ComputerMSER(weighted.tree_, &weighted.altitude_, attr_increasing) {}
	
	/**
	 * @brief Creates an MSER detector that lazily falls back to `AREA`.
	 */
	ComputerMSER(const MorphologicalTree& tree, const AltitudeBuffer* altitude)
		: ComputerMSER(tree, altitude, static_cast<const float*>(nullptr)) { }

	ComputerMSER(const WeightedMorphologicalTree& weighted)
		: ComputerMSER(weighted.tree_, &weighted.altitude_, static_cast<const float*>(nullptr)) { }

	ComputerMSER(const MorphologicalTree& tree) = delete;
	ComputerMSER(const MorphologicalTree& tree, std::vector<float> attr_increasing) = delete;
	ComputerMSER(const MorphologicalTree& tree, const float* attr_increasing) = delete;

	~ComputerMSER(){}

	/**
	 * @brief Computes the MSER indicator vector for the given delta.
	 *
	 * @return A dense boolean-like vector indexed by the tree's internal node
	 * slots, with `true` at the nodes selected as MSER-like regions.
	 */
	std::vector<uint8_t> computeMSER(int delta){
		std::pair<std::vector<NodeId>, std::vector<NodeId>> ascDesc =
			WeightedMorphologicalTree::computeAscendantsAndDescendants(tree, altitude_, delta);
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
	double getStability(NodeId node){
		return (this->getAttrMSER(this->ascendants[node]) - this->getAttrMSER(this->descendants[node])) / this->getAttrMSER(node)  ;
	}

	/**
	 * @brief Returns the attribute used by the MSER criterion, lazily
	 * computing `AREA` when no external buffer has been provided.
	 */
	float getAttrMSER(NodeId node){
			if(attrMserView_ == nullptr) {
				auto area = AttributeComputedIncrementally::computeSingleAttribute(tree, altitude_, AREA);
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
	NodeId getNodeInPathWithMaxStability(NodeId node){
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
	NodeId ascendantWithMaxStability(NodeId node) const { return this->ascendants[node];}
	/**
	 * @brief Returns the descendant used in the current stability window.
	 */
	NodeId descendantWithMaxStability(NodeId node) const { return descendants[node];}
	/**
	 * @brief Returns the current stability array, indexed by node slot.
	 */
	std::vector<float>& getStabilities() { return stability; }
	/**
	 * @brief Returns the number of nodes selected as MSER-like in the last run.
	 */
	int getNumNodes() {return  num;}
	/**
	 * @brief Sets the maximum accepted stability value.
	 */
	void setMaxVariation(float maxVariation) { this->maxVariation = maxVariation; }
	/**
	 * @brief Sets the lower bound of the accepted attribute interval.
	 */
	void setMinAttribute(int minAttr) { this->minAttr = minAttr; }
	/**
	 * @brief Sets the upper bound of the accepted attribute interval.
	 */
	void setMaxAttribute(int maxAttr) { this->maxAttr = maxAttr; }
};

} // namespace mmcfilters
