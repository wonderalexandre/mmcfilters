#pragma once

#include "AttributeComputer.hpp"
#include "ComputerAttributeBasedBitQuads.hpp"
#include "../trees/MorphologicalTree.hpp"


namespace mmcfilters {

/**
 * @brief Computes the family of descriptors derived from bit-quad counts.
 *
 * @details
 * Bit-quad descriptors are obtained by analysing local binary patterns around
 * the contour representation of each component. The low-level extraction is
 * delegated to `ComputerAttributeBasedBitQuads`, which traverses the tree and
 * accumulates the bit-quad counts required to derive Euler number, number of
 * holes, discrete and continuous perimeters, circularity, and average shape
 * measures.
 *
 * This wrapper is intentionally thin: it triggers the low-level computation
 * once, then materialises only the requested descriptors into the attribute
 * buffer.
 *
 * @note These descriptors require a tree with valid connectivity metadata. If
 * the tree has been imported without an adjacency relation, the delegated
 * implementation is expected to reject requests that depend on that missing
 * information.
 */
class BitquadsComputer : public AttributeComputer {
	public:
		/**
		 * @brief Returns the full family of bit-quad descriptors.
		 */
		std::vector<Attribute> attributes() const override {
			return {BITQUADS_AREA,
					BITQUADS_NUMBER_EULER,
					BITQUADS_NUMBER_HOLES,
					BITQUADS_PERIMETER,
					BITQUADS_PERIMETER_CONTINUOUS,
					BITQUADS_CIRCULARITY,
					BITQUADS_PERIMETER_AVERAGE,
					BITQUADS_LENGTH_AVERAGE,
					BITQUADS_WIDTH_AVERAGE};
		}

		/**
		 * @brief Computes the requested bit-quad descriptors for each live node.
		 */
		void compute(const MorphologicalTree& tree, const AltitudeBuffer*, std::span<float> buffer, const AttributeNames& attrNames, std::span<const Attribute> requestedAttributes, std::span<const DependencySource>) const override {
			if(PRINT_LOG) std::cout << "\n==== AttributeComputer: Computing BITQUADS group" << std::endl;
			auto indexOf = [&](int idx, Attribute attr) {
				return attrNames.linearIndex(idx, attr);
			};
			bool computeArea = std::find(requestedAttributes.begin(), requestedAttributes.end(), BITQUADS_AREA) != requestedAttributes.end();
			bool computeNumberEuler = std::find(requestedAttributes.begin(), requestedAttributes.end(), BITQUADS_NUMBER_EULER) != requestedAttributes.end();
			bool computeNumberHoles = std::find(requestedAttributes.begin(), requestedAttributes.end(), BITQUADS_NUMBER_HOLES) != requestedAttributes.end();
			bool computePerimeter = std::find(requestedAttributes.begin(), requestedAttributes.end(), BITQUADS_PERIMETER) != requestedAttributes.end();
			bool computePerimeterCont = std::find(requestedAttributes.begin(), requestedAttributes.end(), BITQUADS_PERIMETER_CONTINUOUS) != requestedAttributes.end();
			bool computeCircularity = std::find(requestedAttributes.begin(), requestedAttributes.end(), BITQUADS_CIRCULARITY) != requestedAttributes.end();
			bool computePerimeterAverage = std::find(requestedAttributes.begin(), requestedAttributes.end(), BITQUADS_PERIMETER_AVERAGE) != requestedAttributes.end();
			bool computeLengthAverage = std::find(requestedAttributes.begin(), requestedAttributes.end(), BITQUADS_LENGTH_AVERAGE) != requestedAttributes.end();
			bool computeWithAverage = std::find(requestedAttributes.begin(), requestedAttributes.end(), BITQUADS_WIDTH_AVERAGE) != requestedAttributes.end();


			ComputerAttributeBasedBitQuads computerBitQuads(tree);
			std::vector<AttributeBasedBitQuads> attr = computerBitQuads.getAttributes();
			for (NodeId nodeId : tree.getAliveNodeIds()) {
				const NodeId node = nodeId;
				if(computeArea)
					buffer[indexOf(node, BITQUADS_AREA)] = attr[node].getAreaDuda();
				if(computeNumberEuler)
					buffer[indexOf(node, BITQUADS_NUMBER_EULER)] = attr[node].getNumberEuler();
				if(computeNumberHoles)
					buffer[indexOf(node, BITQUADS_NUMBER_HOLES)] = attr[node].getNumberHoles();
				if(computePerimeter)
					buffer[indexOf(node, BITQUADS_PERIMETER)] = attr[node].getPerimeter();
				if(computePerimeterCont)
					buffer[indexOf(node, BITQUADS_PERIMETER_CONTINUOUS)] = attr[node].getPerimeterContinuous();
				if(computeCircularity)
					buffer[indexOf(node, BITQUADS_CIRCULARITY)] = attr[node].getCircularity();
				if(computePerimeterAverage)
					buffer[indexOf(node, BITQUADS_PERIMETER_AVERAGE)] = attr[node].getPerimeterAverage();
				if(computeLengthAverage)
					buffer[indexOf(node, BITQUADS_LENGTH_AVERAGE)] = attr[node].getLengthAverage();
				if(computeWithAverage)
					buffer[indexOf(node, BITQUADS_WIDTH_AVERAGE)] = attr[node].getWidthAverage();
			}
			
		}

		void computeUnitAttributes(
			const MorphologicalTree& tree,
			const AltitudeBuffer*,
			std::span<const NodeId> unitProperParts,
			std::span<float> buffer,
			const AttributeNames& attrNames,
			std::span<const Attribute> requestedAttributes) const override
		{
			requireUnitAttributeBufferShape(tree, unitProperParts, buffer, attrNames);
			const AdjacencyRelation* adjacency = tree.getAdjacencyRelation();
			if (adjacency == nullptr) {
				throw std::invalid_argument("BitQuads attributes require an adjacency relation.");
			}

			const bool is4Connectivity = adjacency->is4connectivity();
			const float nan = std::numeric_limits<float>::quiet_NaN();
			const float inf = std::numeric_limits<float>::infinity();

			auto unitValue = [&](Attribute attribute) -> float {
				if (is4Connectivity) {
					switch (attribute) {
						case BITQUADS_AREA: return 0.0f;
						case BITQUADS_NUMBER_EULER: return 1.0f;
						case BITQUADS_NUMBER_HOLES: return 0.0f;
						case BITQUADS_PERIMETER: return 0.0f;
						case BITQUADS_PERIMETER_CONTINUOUS: return 0.0f;
						case BITQUADS_CIRCULARITY: return nan;
						case BITQUADS_PERIMETER_AVERAGE: return 0.0f;
						case BITQUADS_LENGTH_AVERAGE: return 0.0f;
						case BITQUADS_WIDTH_AVERAGE: return nan;
						default: break;
					}
				}

				switch (attribute) {
					case BITQUADS_AREA: return 1.0f;
					case BITQUADS_NUMBER_EULER: return 0.0f;
					case BITQUADS_NUMBER_HOLES: return 1.0f;
					case BITQUADS_PERIMETER: return 4.0f;
					case BITQUADS_PERIMETER_CONTINUOUS: return 8.0f / 3.0f;
					case BITQUADS_CIRCULARITY: return 9.0f * std::numbers::pi_v<float> / 16.0f;
					case BITQUADS_PERIMETER_AVERAGE: return inf;
					case BITQUADS_LENGTH_AVERAGE: return inf;
					case BITQUADS_WIDTH_AVERAGE: return nan;
					default: break;
				}
				throw std::runtime_error("Unsupported BitQuads unit attribute.");
			};

			for (Attribute attribute : attributes()) {
				if (!requestsAttribute(requestedAttributes, attribute)) {
					continue;
				}
				const float value = unitValue(attribute);
				for (NodeId leafIndex = 0; leafIndex < static_cast<NodeId>(unitProperParts.size()); ++leafIndex) {
					buffer[attrNames.linearIndex(leafIndex, attribute)] = value;
				}
			}
		}
};

} // namespace mmcfilters
