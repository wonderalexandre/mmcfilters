#pragma once

#include "AttributeRegistry.hpp"
#include "AttributeTypes.hpp"
#include "../utils/Contract.hpp"

#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace mmcfilters {
/**
 * @brief Canonical expansion of each public attribute group.
 *
 * @details
 * The map is used by public convenience APIs and by the internal request
 * expansion helpers. Its contents therefore define the library-wide meaning of
 * each `AttributeGroup`.
 */
inline const std::unordered_map<AttributeGroup, std::vector<Attribute>>& ATTRIBUTE_GROUPS = attributes::registry::attributeGroups();

/**
 * @brief Layout for node attributes sampled at signed ancestor/current/descendant offsets.
 *
 * @details
 * The layout is node-major: for a given node, all requested sample positions
 * are stored contiguously. `NUM_ATTRIBUTES` therefore counts the number of
 * `(attribute, sampleOffset)` entries associated with each node.
 */
class NodeAttributeSampleLayout {
  public:
    /// Maps each `(attribute, sampleOffset)` key to its per-node offset.
    std::unordered_map<NodeAttributeSampleKey, int> indexMap;

    /// Number of sampled entries stored for each node.
    const int NUM_ATTRIBUTES;

    /**
     * @brief Constructs a layout from an already validated offset map.
     *
     * The constructor takes ownership of `map` and derives the per-node stride
     * from its final size.
     *
     * @param map Mapping used to initialize the layout.
     */
    NodeAttributeSampleLayout(std::unordered_map<NodeAttributeSampleKey, int>&& map)
        : indexMap(std::move(map)), NUM_ATTRIBUTES(static_cast<int>(indexMap.size())) {}

    /**
     * @brief Builds the canonical dense layout for offsets in `[-samplingRadius, samplingRadius]`.
     *
     * @param samplingRadius Maximum ancestor/descendant sample offset.
     * @param attributes Attributes requested by the operation.
     * @return The resulting canonical dense sampled layout.
     */
    [[nodiscard]] static NodeAttributeSampleLayout create(int samplingRadius, const std::vector<Attribute>& attributes) {
        std::unordered_map<NodeAttributeSampleKey, int> map;
        int offset = 0;
        for (int sampleOffset = -samplingRadius; sampleOffset <= samplingRadius; ++sampleOffset) {
            for (std::size_t i = 0; i < attributes.size(); ++i) {
                const NodeAttributeSampleKey key{attributes[i], sampleOffset};
                const auto [_, inserted] = map.emplace(key, offset++);
                MMCFILTERS_CONTRACT_REQUIRE(
                    inserted, throw std::invalid_argument("NodeAttributeSampleLayout::create received duplicate attribute " + toString(key)));
            }
        }
        return NodeAttributeSampleLayout(std::move(map));
    }

    /**
     * @brief Returns the per-node offset associated with `(attribute, sampleOffset)`.
     *
     * @param attribute Attribute requested by the operation.
     * @param sampleOffset Signed sample coordinate.
     * @return The per-node offset associated with the key.
     */
    [[nodiscard]] int getIndex(Attribute attribute, int sampleOffset) const {
        return getIndex(NodeAttributeSampleKey{attribute, sampleOffset});
    }

    /**
     * @brief Returns the per-node offset associated with a composite sample key.
     *
     * @param sampleKey Composite sampled-attribute key.
     * @return The per-node offset associated with the key.
     */
    [[nodiscard]] int getIndex(NodeAttributeSampleKey sampleKey) const { return indexMap.at(sampleKey); }

    /**
     * @brief Returns the linear index in the flat buffer for a given node and
     * sampled attribute.
     *
     * @param nodeIndex Index.
     * @param attribute Attribute requested by the operation.
     * @param sampleOffset Signed sample coordinate.
     * @return The linear index in the flat buffer.
     */
    [[nodiscard]] int linearIndex(int nodeIndex, Attribute attribute, int sampleOffset) const {
        return nodeIndex * NUM_ATTRIBUTES + getIndex(attribute, sampleOffset);
    }

    /**
     * @brief Convenience overload taking a `NodeAttributeSampleKey`.
     *
     * @param nodeIndex Index.
     * @param sampleKey Composite sampled-attribute key.
     * @return Value returned by the convenience overload.
     */
    [[nodiscard]] int linearIndex(int nodeIndex, NodeAttributeSampleKey sampleKey) const {
        return linearIndex(nodeIndex, sampleKey.attribute, sampleKey.sampleOffset);
    }

    /**
     * @brief Returns a human-readable label for a composite key.
     *
     * @param sampleKey Composite sampled-attribute key.
     * @return A human-readable label for a composite key.
     */
    [[nodiscard]] static std::string toString(NodeAttributeSampleKey sampleKey) {
        return toString(sampleKey.attribute, sampleKey.sampleOffset);
    }

    /**
     * @brief Returns a serialized label for an attribute at a sample offset.
     *
     * @param attribute Attribute requested by the operation.
     * @param sampleOffset Signed sample coordinate.
     * @return The canonical serialized label.
     */
    [[nodiscard]] static std::string toString(Attribute attribute, int sampleOffset) {
        std::string name(attributes::registry::name(attribute));

        if (sampleOffset < 0) {
            name += "_ANCESTOR_" + std::to_string(-sampleOffset);
        } else if (sampleOffset > 0) {
            name += "_DESCENDANT_" + std::to_string(sampleOffset);
        }

        return name;
    }
};

/**
 * @brief Layout object that maps scalar attributes to flat-buffer offsets.
 *
 * @details
 * `AttributeNames` describes the in-memory arrangement of a node-attribute
 * buffer. The layout is node-major:
 * - all requested attributes for node `0` come first;
 * - then all requested attributes for node `1`;
 * - and so on.
 *
 * `NUM_ATTRIBUTES` is therefore the per-node stride needed to jump from the
 * attributes of one node to the next. The actual number of floats in a buffer
 * is `numNodes * NUM_ATTRIBUTES`.
 *
 * The class is deliberately lightweight because it is shared across most of
 * the incremental attribute pipeline, dependency maps, and Python bindings.
 */
class AttributeNames {
  public:
    /// Maps each requested attribute to its per-node offset.
    std::unordered_map<Attribute, int> indexMap;

    /// Number of scalar attributes stored for each node.
    const int NUM_ATTRIBUTES;

    /**
     * @brief Constructs a layout from an already validated offset map.
     *
     * The constructor takes ownership of `map` and derives the per-node stride
     * from its final size.
     *
     * @param map Mapping used to initialize the layout.
     */
    AttributeNames(std::unordered_map<Attribute, int>&& map) : indexMap(std::move(map)), NUM_ATTRIBUTES(static_cast<int>(indexMap.size())) {}

    /**
     * @brief Builds a dense layout from an explicit attribute list.
     *
     * Duplicate attributes are rejected because a repeated key would make the
     * flat-buffer stride disagree with the lookup table.
     *
     * @param attributes Attributes requested by the operation.
     * @return The resulting dense layout from an explicit attribute list.
     */
    [[nodiscard]] static AttributeNames fromList(const std::vector<Attribute>& attributes) {
        std::unordered_map<Attribute, int> map;
        int i = 0;
        for (auto attr : attributes) {
            const auto [_, inserted] = map.emplace(attr, i++);
            MMCFILTERS_CONTRACT_REQUIRE(inserted,
                                        throw std::invalid_argument("AttributeNames::fromList received duplicate attribute " + toString(attr)));
        }
        return AttributeNames(std::move(map));
    }

    /**
     * @brief Builds a layout from a public attribute group.
     *
     * @param group Public attribute group.
     * @return The resulting layout from a public attribute group.
     */
    [[nodiscard]] static AttributeNames fromGroup(AttributeGroup group) {
        auto it = ATTRIBUTE_GROUPS.find(group);
        MMCFILTERS_CONTRACT_REQUIRE(it != ATTRIBUTE_GROUPS.end(), throw std::invalid_argument("Unknown attribute group."));
        return fromList(it->second);
    }

    /**
     * @brief Returns the per-node offset associated with `attr`.
     *
     * @param attr Attribute requested by the operation.
     * @return The per-node offset associated with attr.
     */
    [[nodiscard]] int getIndex(Attribute attr) const { return indexMap.at(attr); }

    /**
     * @brief Returns whether this layout contains `attr`.
     *
     * @param attr Attribute requested by the operation.
     * @return Whether this layout contains attr.
     */
    [[nodiscard]] bool contains(Attribute attr) const noexcept { return indexMap.contains(attr); }

    /**
     * @brief Returns the per-node offset associated with `attr`.
     *
     * @param attr Attribute requested by the operation.
     * @return The per-node offset associated with attr.
     */
    [[nodiscard]] int offset(Attribute attr) const { return getIndex(attr); }

    /**
     * @brief Returns the linear index in the flat buffer for `(node, attr)`.
     *
     * @param nodeIndex Index.
     * @param attr Attribute requested by the operation.
     * @return The linear index in the flat buffer for (node, attr).
     */
    [[nodiscard]] int linearIndex(int nodeIndex, Attribute attr) const { return nodeIndex * NUM_ATTRIBUTES + getIndex(attr); }

    /**
     * @brief Returns the stable symbolic name associated with `attr`.
     *
     * @param attr Attribute requested by the operation.
     * @return The stable symbolic name associated with attr.
     */
    [[nodiscard]] static std::string toString(Attribute attr) { return std::string(attributes::registry::name(attr)); }

    /**
     * @brief Returns a user-facing description of an attribute.
     *
     * @details
     * The descriptions are intentionally richer than `toString(...)` and are
     * meant for documentation, notebooks, and Python introspection helpers.
     *
     * @param attr Attribute requested by the operation.
     * @return A user-facing description of an attribute.
     */
    static std::string describe(Attribute attr) { return std::string(attributes::registry::description(attr)); }

    /**
     * @brief Parses a stable symbolic attribute name.
     *
     * @param str Symbolic attribute name to parse.
     * @return The matching attribute when `str` names a known descriptor, or
     * `std::nullopt` otherwise.
     *
     */
    [[nodiscard]] static std::optional<Attribute> parse(const std::string& str) { return attributes::registry::parse(str); }
};

} // namespace mmcfilters
