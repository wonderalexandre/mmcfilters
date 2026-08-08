#pragma once

#include <functional>
#include <type_traits>
#include <variant>

namespace mmcfilters {

/**
 * @brief Enumeration of scalar node attributes supported by the library.
 *
 * @details
 * Each enumerator identifies a descriptor that can be materialised per live
 * node of a `MorphologicalTree`. The values are used in three roles:
 * - as the public request space exposed by the API;
 * - as keys in dependency maps and attribute layouts;
 * - as dispatch keys inside the attribute pipeline and topology backend.
 *
 * The enumeration is intentionally flat: grouped requests are expressed
 * separately through `AttributeGroup`.
 */
enum class Attribute {
    AREA,

    VOLUME,
    RELATIVE_VOLUME,
    LEVEL,
    GRAY_HEIGHT,
    MEAN_LEVEL,
    VARIANCE_LEVEL,

    BOX_WIDTH,
    BOX_HEIGHT,
    DIAGONAL_LENGTH,
    RECTANGULARITY,
    RATIO_WH,
    BOX_COL_MIN,
    BOX_COL_MAX,
    BOX_ROW_MIN,
    BOX_ROW_MAX,

    CENTRAL_MOMENT_20,
    CENTRAL_MOMENT_02,
    CENTRAL_MOMENT_11,
    CENTRAL_MOMENT_30,
    CENTRAL_MOMENT_03,
    CENTRAL_MOMENT_21,
    CENTRAL_MOMENT_12,

    HU_MOMENT_1,
    HU_MOMENT_2,
    HU_MOMENT_3,
    HU_MOMENT_4,
    HU_MOMENT_5,
    HU_MOMENT_6,
    HU_MOMENT_7,

    INERTIA,
    COMPACTNESS,
    ECCENTRICITY,
    LENGTH_MAJOR_AXIS,
    LENGTH_MINOR_AXIS,
    AXIS_ORIENTATION,
    CIRCULARITY,

    BITQUADS_AREA,
    BITQUADS_NUMBER_EULER,
    BITQUADS_NUMBER_HOLES,
    BITQUADS_PERIMETER,
    BITQUADS_PERIMETER_CONTINUOUS,
    BITQUADS_CIRCULARITY,
    BITQUADS_PERIMETER_AVERAGE,
    BITQUADS_LENGTH_AVERAGE,
    BITQUADS_WIDTH_AVERAGE,

    HEIGHT_NODE,
    DEPTH_NODE,
    IS_LEAF_NODE,
    IS_ROOT_NODE,
    NUM_CHILDREN_NODE,
    NUM_SIBLINGS_NODE,
    NUM_DESCENDANTS_NODE,
    NUM_LEAF_DESCENDANTS_NODE,
    LEAF_RATIO_NODE,
    BALANCE_NODE,

    MAX_DIST,

    AVG_CHILD_HEIGHT_NODE,

    CONTOUR_PIXELS,
    CONTOUR_PERIMETER,
    CONTOUR_SIDE_NORTH,
    CONTOUR_SIDE_WEST,
    CONTOUR_SIDE_EAST,
    CONTOUR_SIDE_SOUTH
};

/**
 * @brief Enumeration of attribute groups used for bulk requests.
 */
enum class AttributeGroup {
    ALL,
    GRAY_LEVEL,
    SHAPE,
    MOMENTS,
    BOUNDARY,
    TREE_TOPOLOGY,
};

using AttributeOrGroup = std::variant<Attribute, AttributeGroup>;
using enum Attribute;

/**
 * @brief Composite key used to index delta-augmented attribute layouts.
 */
struct AttributeKey {
    /// Attribute component of the composite lookup key.
    Attribute attr;

    /// Relative ancestor/descendant offset associated with `attr`.
    int delta = 0;

    /**
     * @brief Builds a composite key for `a` at delta offset `d`.
     *
     * A zero delta denotes the current node. Negative deltas refer to ancestor
     * samples and positive deltas refer to descendant samples in delta-aware
     * layouts.
     *
     * @param a Attribute component of the composite key.
     * @param d Delta component of the composite key.
     */
    AttributeKey(Attribute a, int d = 0) : attr(a), delta(d) {}

    /**
     * @brief Returns true when both the attribute and delta offset match.
     *
     * @param other Object to compare with or transfer from.
     * @return True when both the attribute and delta offset match.
     */
    bool operator==(const AttributeKey& other) const { return attr == other.attr && delta == other.delta; }
};

} // namespace mmcfilters

namespace std {
/** @brief Provides standard-library hashing for `mmcfilters::AttributeGroup`. */
template <> struct hash<mmcfilters::AttributeGroup> {
    /**
     * @brief Applies the function-call operator.
     *
     * @param group Attribute group whose metadata is requested.
     * @return Hash value for the supplied key.
     */
    std::size_t operator()(const mmcfilters::AttributeGroup& group) const noexcept { return static_cast<std::size_t>(group); }
};

/** @brief Provides standard-library hashing for `mmcfilters::Attribute`. */
template <> struct hash<mmcfilters::Attribute> {
    /**
     * @brief Applies the function-call operator.
     *
     * @param attr Attribute requested by the operation.
     * @return Hash value for the supplied key.
     */
    std::size_t operator()(const mmcfilters::Attribute& attr) const noexcept { return static_cast<std::size_t>(attr); }
};

/** @brief Provides standard-library hashing for `mmcfilters::AttributeOrGroup`. */
template <> struct hash<mmcfilters::AttributeOrGroup> {
    /**
     * @brief Applies the function-call operator.
     *
     * @param attr Attribute requested by the operation.
     * @return Hash value for the supplied key.
     */
    std::size_t operator()(const mmcfilters::AttributeOrGroup& attr) const {
        return std::visit([](auto&& a) -> std::size_t { return std::hash<std::decay_t<decltype(a)>>{}(a); }, attr);
    }
};

/** @brief Provides standard-library hashing for `mmcfilters::AttributeKey`. */
template <> struct hash<mmcfilters::AttributeKey> {
    /**
     * @brief Applies the function-call operator.
     *
     * @param k Key whose hash value is computed.
     * @return Hash value for the supplied key.
     */
    std::size_t operator()(const mmcfilters::AttributeKey& k) const { return std::hash<int>()(static_cast<int>(k.attr)) ^ (std::hash<int>()(k.delta) << 1); }
};
} // namespace std
