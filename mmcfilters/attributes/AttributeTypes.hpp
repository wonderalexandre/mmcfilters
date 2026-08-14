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
    Area,

    Volume,
    RelativeVolume,
    GrayLevelHeight,
    MeanGrayLevel,
    GrayLevelVariance,

    BoxWidth,
    BoundingBoxHeight,
    DiagonalLength,
    Rectangularity,
    RatioWh,
    BoxColumnMin,
    BoxColumnMax,
    BoxRowMin,
    BoxRowMax,

    CentralMoment20,
    CentralMoment02,
    CentralMoment11,
    CentralMoment30,
    CentralMoment03,
    CentralMoment21,
    CentralMoment12,

    HuMoment1,
    HuMoment2,
    HuMoment3,
    HuMoment4,
    HuMoment5,
    HuMoment6,
    HuMoment7,

    Inertia,
    Compactness,
    Eccentricity,
    LengthMajorAxis,
    LengthMinorAxis,
    AxisOrientation,
    Circularity,

    BitquadArea,
    BitquadNumberEuler,
    BitquadNumberHoles,
    BitquadPerimeter,
    BitquadPerimeterContinuous,
    BitquadCircularity,
    BitquadPerimeterAverage,
    BitquadLengthAverage,
    BitquadWidthAverage,

    SubtreeHeight,
    DepthNode,
    IsLeafNode,
    IsRootNode,
    NumChildrenNode,
    NumSiblingsNode,
    NumDescendantsNode,
    NumLeafDescendantsNode,
    LeafRatioNode,
    BalanceNode,

    MaxDist,

    AvgChildHeightNode,

    ContourPixels,
    ContourPerimeter,
    ContourSideNorth,
    ContourSideWest,
    ContourSideEast,
    ContourSideSouth
};

/**
 * @brief Enumeration of attribute groups used for bulk requests.
 */
enum class AttributeGroup {
    All,
    GrayLevel,
    Shape,
    Moments,
    Boundary,
    TreeTopology,
};

using AttributeOrGroup = std::variant<Attribute, AttributeGroup>;
using enum Attribute;

/**
 * @brief Policy used to select one representative descendant sample.
 */
enum class NodeAttributeSamplingPolicy {
    LargestSupportDescendant,
};

/**
 * @brief Policy used when a requested node-attribute sample is unavailable.
 */
enum class MissingNodeAttributeSamplePolicy {
    RepeatNearest,
    NotANumber,
    Zero,
};

/**
 * @brief Composite key used to index sampled node-attribute layouts.
 */
struct NodeAttributeSampleKey {
    /// Attribute component of the composite lookup key.
    Attribute attribute;

    /// Signed ancestor/current/representative-descendant sample coordinate.
    int sampleOffset = 0;

    /**
     * @brief Builds a composite key for `value` at `offset`.
     *
     * A zero offset denotes the current node. Negative offsets refer to
     * ancestor samples and positive offsets refer to representative-descendant
     * samples.
     *
     * @param value Attribute component of the composite key.
     * @param offset Sample coordinate of the composite key.
     */
    NodeAttributeSampleKey(Attribute value, int offset = 0) : attribute(value), sampleOffset(offset) {}

    /**
     * @brief Returns true when both the attribute and sample offset match.
     *
     * @param other Object to compare with or transfer from.
     * @return True when both key components match.
     */
    bool operator==(const NodeAttributeSampleKey& other) const {
        return attribute == other.attribute && sampleOffset == other.sampleOffset;
    }
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

/** @brief Provides standard-library hashing for `mmcfilters::NodeAttributeSampleKey`. */
template <> struct hash<mmcfilters::NodeAttributeSampleKey> {
    /**
     * @brief Applies the function-call operator.
     *
     * @param key Key whose hash value is computed.
     * @return Hash value for the supplied key.
     */
    std::size_t operator()(const mmcfilters::NodeAttributeSampleKey& key) const {
        return std::hash<int>()(static_cast<int>(key.attribute)) ^ (std::hash<int>()(key.sampleOffset) << 1);
    }
};
} // namespace std
