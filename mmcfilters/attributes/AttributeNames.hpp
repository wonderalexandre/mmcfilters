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
 * @brief Layout object for buffers that store several `(attribute, delta)`
 * combinations per node.
 *
 * @details
 * The layout is node-major: for a given node, all requested delta variants are
 * stored contiguously. `NUM_ATTRIBUTES` therefore counts the number of
 * `(attribute, delta)` entries associated with each node, not the number of
 * live nodes.
 *
 * The class is primarily used by delta-based APIs such as MSER-like stability
 * descriptors, where each node needs access to values sampled at ancestor and
 * descendant offsets around the current position.
 */
class AttributeNamesWithDelta {
  public:
    /// Maps each `(attribute, delta)` key to its per-node offset.
    std::unordered_map<AttributeKey, int> indexMap;

    /// Number of delta-augmented entries stored for each node.
    const int NUM_ATTRIBUTES;

    /**
     * @brief Constructs a layout from an already validated offset map.
     *
     * The constructor takes ownership of `map` and derives the per-node stride
     * from its final size.
     *
     * @param map Mapping used to initialize the layout.
     */
    AttributeNamesWithDelta(std::unordered_map<AttributeKey, int>&& map) : indexMap(std::move(map)), NUM_ATTRIBUTES(static_cast<int>(indexMap.size())) {}

    /**
     * @brief Builds the canonical dense layout for all deltas in `[-delta, delta]`.
     *
     * @param delta Delta offset or radius used by the operation.
     * @param attributes Attributes requested by the operation.
     * @return The resulting canonical dense layout for all deltas in [-delta, delta].
     */
    [[nodiscard]] static AttributeNamesWithDelta create(int delta, const std::vector<Attribute>& attributes) {
        std::unordered_map<AttributeKey, int> map;
        int offset = 0;
        for (int d = -delta; d <= delta; ++d) {
            for (std::size_t i = 0; i < attributes.size(); ++i) {
                const AttributeKey key{attributes[i], d};
                const auto [_, inserted] = map.emplace(key, offset++);
                MMCFILTERS_CONTRACT_REQUIRE(
                    inserted, throw std::invalid_argument("AttributeNamesWithDelta::create received duplicate attribute " + toString(key)));
            }
        }
        return AttributeNamesWithDelta(std::move(map));
    }

    /**
     * @brief Returns the per-node offset associated with `(attr, delta)`.
     *
     * @param attr Attribute requested by the operation.
     * @param delta Delta offset or radius used by the operation.
     * @return The per-node offset associated with (attr, delta).
     */
    [[nodiscard]] int getIndex(Attribute attr, int delta) const { return getIndex(AttributeKey{attr, delta}); }

    /**
     * @brief Returns the per-node offset associated with a composite delta key.
     *
     * @param attrKey Attribute information represented by `attrKey`.
     * @return The per-node offset associated with a composite delta key.
     */
    [[nodiscard]] int getIndex(AttributeKey attrKey) const { return indexMap.at(attrKey); }

    /**
     * @brief Returns the linear index in the flat buffer for a given node and
     * delta-augmented attribute.
     *
     * @param nodeIndex Index represented by `nodeIndex`.
     * @param attr Attribute requested by the operation.
     * @param delta Delta offset or radius used by the operation.
     * @return The linear index in the flat buffer for a given node and delta-augmented attribute.
     */
    [[nodiscard]] int linearIndex(int nodeIndex, Attribute attr, int delta) const { return nodeIndex * NUM_ATTRIBUTES + getIndex(attr, delta); }

    /**
     * @brief Convenience overload taking an `AttributeKey`.
     *
     * @param nodeIndex Index represented by `nodeIndex`.
     * @param attrKey Attribute information represented by `attrKey`.
     * @return Value returned by the convenience overload.
     */
    [[nodiscard]] int linearIndex(int nodeIndex, AttributeKey attrKey) const { return linearIndex(nodeIndex, attrKey.attr, attrKey.delta); }

    /**
     * @brief Returns a human-readable label for a composite key.
     *
     * @param attrKey Attribute information represented by `attrKey`.
     * @return A human-readable label for a composite key.
     */
    [[nodiscard]] static std::string toString(AttributeKey attrKey) { return toString(attrKey.attr, attrKey.delta); }

    /**
     * @brief Returns a human-readable label for an attribute at a given delta.
     *
     * @param attr Attribute requested by the operation.
     * @param delta Delta offset or radius used by the operation.
     * @return A human-readable label for an attribute at a given delta.
     */
    [[nodiscard]] static std::string toString(Attribute attr, int delta) {
        std::string name(attributes::registry::name(attr));

        if (delta < 0) {
            name += "_ASC_" + std::to_string(-delta);
        } else if (delta > 0) {
            name += "_DESC_" + std::to_string(delta);
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
     * @param nodeIndex Index represented by `nodeIndex`.
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
