#pragma once

#include "AttributeNumericPolicy.hpp"
#include "../../trees/TreeAltitudeAlgorithms.hpp"
#include "../../trees/ValuedMorphologicalTree.hpp"
#include "../../trees/detail/CommittedTreeAccess.hpp"
#include "../../utils/Common.hpp"
#include "../../utils/Contract.hpp"
#include "../AttributeNames.hpp"

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <initializer_list>
#include <span>
#include <stdexcept>

namespace mmcfilters {

/**
 * @brief Non-owning view over a typed dependency attribute buffer.
 *
 * @details
 * Attribute kernels receive upstream dependencies through lightweight views
 * rather than ownership-bearing objects. A `DependencySourceT<Real>` bundles
 * the dependency layout and a borrowed pointer to values of the selected real
 * type. Rows are always indexed in dense internal `NodeId` space unless a
 * caller explicitly documents a different projection.
 */
template <std::floating_point Real> struct DependencySourceT {
    /// Layout describing the dependency buffer; borrowed and non-owning.
    const AttributeNames* attrNames = nullptr;

    /// Flat dependency values indexed through `attrNames`; borrowed and non-owning.
    const Real* buffer = nullptr;
};

/**
 * @brief Named resolver over borrowed dependency attribute buffers.
 *
 * @details
 * `DependencyResolver` lets kernels ask for semantic dependencies by
 * `Attribute`, which makes the contract explicit while preserving the same
 * borrowed-buffer representation. It does not own or copy dependency data; it
 * only validates that the borrowed sources are present and expose the requested
 * columns.
 */
template <std::floating_point Real> class DependencyResolver {
    /** @brief Sources buffer. */
    std::span<const DependencySourceT<Real>> sources_;

    /**
     * @brief Validates valid source.
     *
     * @param source Source value or object.
     */
    static void requireValidSource(const DependencySourceT<Real>& source) {
        if (source.attrNames == nullptr || source.buffer == nullptr) {
            throw std::invalid_argument("Requested attribute computation received an invalid dependency source.");
        }
    }

  public:
    /**
     * @brief Creates a resolver over borrowed dependency sources.
     *
     * The span and every referenced `AttributeNames`/buffer must outlive the
     * resolver and any compute context that owns it.
     *
     * @param sources Input.
     */
    explicit DependencyResolver(std::span<const DependencySourceT<Real>> sources) noexcept : sources_(sources) {}

    /**
     * @brief Returns the first source containing `requiredAttribute`.
     *
     * @param requiredAttribute Attribute information.
     * @return The first source containing requiredAttribute.
     *
     * @throws std::invalid_argument if any inspected source is invalid or no
     * source contains the requested attribute.
     *
     */
    [[nodiscard]] const DependencySourceT<Real>& require(Attribute requiredAttribute) const {
        for (const DependencySourceT<Real>& source : sources_) {
            requireValidSource(source);
            if (source.attrNames->contains(requiredAttribute)) {
                return source;
            }
        }
        throw std::invalid_argument("Requested attribute computation is missing a required dependency attribute.");
    }

    /**
     * @brief Returns one source that contains all requested attributes.
     *
     * @details
     * This is used by computers whose formulas consume a coherent dependency
     * block, such as several central moments stored in the same buffer layout.
     *
     * @param requiredAttributes Attribute information.
     * @return One source that contains all requested attributes.
     *
     * @throws std::invalid_argument if the request is empty, if an inspected
     * source is invalid, or if no single source contains every requested
     * attribute.
     *
     */
    [[nodiscard]] const DependencySourceT<Real>& requireAll(std::initializer_list<Attribute> requiredAttributes) const {
        if (requiredAttributes.size() == 0) {
            throw std::invalid_argument("Requested attribute computation must name at least one required dependency attribute.");
        }

        for (const DependencySourceT<Real>& source : sources_) {
            requireValidSource(source);
            const bool hasAll =
                std::all_of(requiredAttributes.begin(), requiredAttributes.end(), [&](Attribute attribute) { return source.attrNames->contains(attribute); });
            if (hasAll) {
                return source;
            }
        }
        throw std::invalid_argument("Requested attribute computation is missing a dependency source with all required attributes.");
    }
};

/**
 * @brief Non-owning topology attribute computation context.
 *
 * @details
 * This is the public adapter boundary for topology/support computer families.
 * It intentionally borrows every object: the tree, output buffer, attribute
 * layout, request list, and dependency views must outlive the context.
 *
 * `buffer` is interpreted through `attrNames` in dense internal-node order.
 * `requestedAttributes` is the subset the caller expects this family to fill;
 * hidden dependency sources can be accessed by semantic name through
 * `dependencies`.
 */
template <std::floating_point Real> struct AttributeComputeContext {
    /// Topology whose dense internal node ids index `buffer`.
    const MorphologicalTree& tree;

    /// Caller-owned output buffer in `tree` internal-node space.
    std::span<Real> buffer;

    /// Column layout for `buffer`.
    const AttributeNames& attrNames;

    /// Attribute subset to materialize in this call.
    std::span<const Attribute> requestedAttributes;

    /// Borrowed dependency views, kept for APIs that need raw span access.
    std::span<const DependencySourceT<Real>> dependencySources;

    /// Semantic dependency lookup over `dependencySources`.
    DependencyResolver<Real> dependencies;

    /**
     * @brief Binds borrowed topology, output, request, and dependency views.
     *
     * @param tree_ Tree topology.
     * @param buffer_ Buffer read or written by the operation.
     * @param attrNames_ Layout that maps attributes to buffer columns.
     * @param requestedAttributes_ Attributes requested for materialization.
     * @param dependencySources_ Input.
     */
    AttributeComputeContext(const MorphologicalTree& tree_, std::span<Real> buffer_, const AttributeNames& attrNames_,
                            std::span<const Attribute> requestedAttributes_, std::span<const DependencySourceT<Real>> dependencySources_ = {}) noexcept
        : tree(tree_), buffer(buffer_), attrNames(attrNames_), requestedAttributes(requestedAttributes_), dependencySources(dependencySources_),
          dependencies(dependencySources_) {}
};

/**
 * @brief Non-owning altitude-aware attribute computation context.
 *
 * @details
 * Altitude-dependent kernels receive the same topology/output contract as
 * `AttributeComputeContext` plus a typed altitude span indexed by dense
 * internal node id. The altitude span is borrowed and must outlive the context.
 */
template <std::floating_point Real, AltitudeValue T> struct AltitudeAttributeComputeContext {
    /// Topology whose dense internal node ids index `altitude` and `buffer`.
    const MorphologicalTree& tree;

    /// Borrowed altitude values indexed by internal node id.
    std::span<const T> altitude;

    /// Caller-owned output buffer in `tree` internal-node space.
    std::span<Real> buffer;

    /// Column layout for `buffer`.
    const AttributeNames& attrNames;

    /// Attribute subset to materialize in this call.
    std::span<const Attribute> requestedAttributes;

    /// Borrowed dependency views, kept for APIs that need raw span access.
    std::span<const DependencySourceT<Real>> dependencySources;

    /// Semantic dependency lookup over `dependencySources`.
    DependencyResolver<Real> dependencies;

    /**
     * @brief Binds borrowed topology, altitude, output, request, and dependencies.
     *
     * @param tree_ Tree topology.
     * @param altitude_ Altitude data indexed by node identifier.
     * @param buffer_ Buffer read or written by the operation.
     * @param attrNames_ Layout that maps attributes to buffer columns.
     * @param requestedAttributes_ Attributes requested for materialization.
     * @param dependencySources_ Input.
     */
    AltitudeAttributeComputeContext(const MorphologicalTree& tree_, std::span<const T> altitude_, std::span<Real> buffer_, const AttributeNames& attrNames_,
                                    std::span<const Attribute> requestedAttributes_, std::span<const DependencySourceT<Real>> dependencySources_ = {}) noexcept
        : tree(tree_), altitude(altitude_), buffer(buffer_), attrNames(attrNames_), requestedAttributes(requestedAttributes_),
          dependencySources(dependencySources_), dependencies(dependencySources_) {}
};

/**
 * @brief Non-owning unit-row computation context for topology/support families.
 *
 * @details
 * A compact exported Higra layout contains one unit row per pixel leaf.
 * `unitPixels` lists pixel identifiers in output row order; `attrNames`
 * describes the columns of each row in `buffer`.
 */
template <std::floating_point Real> struct UnitAttributeComputeContext {
    /// Topology supplying the smallest node of each pixel in `unitPixels`.
    const MorphologicalTree& tree;

    /// Pixel identifiers in output row order.
    std::span<const PixelId> unitPixels;

    /// Caller-owned output buffer in pixel row order.
    std::span<Real> buffer;

    /// Column layout for `buffer`.
    const AttributeNames& attrNames;

    /// Attribute subset to materialize in this call.
    std::span<const Attribute> requestedAttributes;

    /**
     * @brief Binds borrowed topology, unit rows, output, and request.
     *
     * @param tree_ Tree topology.
     * @param unitPixels_ Pixel identifiers in output row order.
     * @param buffer_ Buffer read or written by the operation.
     * @param attrNames_ Layout that maps attributes to buffer columns.
     * @param requestedAttributes_ Attributes requested for materialization.
     */
    UnitAttributeComputeContext(const MorphologicalTree& tree_, std::span<const PixelId> unitPixels_, std::span<Real> buffer_,
                                const AttributeNames& attrNames_, std::span<const Attribute> requestedAttributes_) noexcept
        : tree(tree_), unitPixels(unitPixels_), buffer(buffer_), attrNames(attrNames_), requestedAttributes(requestedAttributes_) {}
};

/**
 * @brief Non-owning unit-row context for altitude-aware families.
 *
 * @details
 * This is the unit-row counterpart of
 * `AltitudeAttributeComputeContext<Real, T>`. The altitude span is indexed in
 * dense internal node-id space; each unit row reads the altitude of its
 * pixel's smallest node.
 */
template <std::floating_point Real, AltitudeValue T> struct AltitudeUnitAttributeComputeContext {
    /// Topology supplying the smallest node of each pixel in `unitPixels`.
    const MorphologicalTree& tree;

    /// Borrowed altitude values indexed by internal node id.
    std::span<const T> altitude;

    /// Pixel identifiers in output row order.
    std::span<const PixelId> unitPixels;

    /// Caller-owned output buffer in pixel row order.
    std::span<Real> buffer;

    /// Column layout for `buffer`.
    const AttributeNames& attrNames;

    /// Attribute subset to materialize in this call.
    std::span<const Attribute> requestedAttributes;

    /**
     * @brief Binds borrowed topology, altitude, unit rows, output, and request.
     *
     * @param tree_ Tree topology.
     * @param altitude_ Altitude data indexed by node identifier.
     * @param unitPixels_ Pixel identifiers in output row order.
     * @param buffer_ Buffer read or written by the operation.
     * @param attrNames_ Layout that maps attributes to buffer columns.
     * @param requestedAttributes_ Attributes requested for materialization.
     */
    AltitudeUnitAttributeComputeContext(const MorphologicalTree& tree_, std::span<const T> altitude_, std::span<const PixelId> unitPixels_,
                                        std::span<Real> buffer_, const AttributeNames& attrNames_, std::span<const Attribute> requestedAttributes_) noexcept
        : tree(tree_), altitude(altitude_), unitPixels(unitPixels_), buffer(buffer_), attrNames(attrNames_),
          requestedAttributes(requestedAttributes_) {}
};

template <AltitudeValue T> inline T unitAltitude(const MorphologicalTree& tree, std::span<const T> altitudeView, PixelId pixel);

namespace detail::kernel {
template <AltitudeValue T> inline T unitAltitude(const MorphologicalTree& tree, std::span<const T> altitudeView, PixelId pixel);
}

/**
 * @brief Validates a dense internal-node output buffer.
 *
 * @details
 * All ordinary attribute computers write one row per internal node slot,
 * including dead slots. Live-node semantics are enforced by traversal code;
 * this helper only checks the flat storage shape.
 *
 * @param tree Tree topology.
 * @param buffer Buffer read or written by the operation.
 * @param attrNames Layout that maps attributes to buffer columns.
 */
template <std::floating_point Real>
inline void requireAttributeBufferShape(const MorphologicalTree& tree, std::span<Real> buffer, const AttributeNames& attrNames) {
    const std::size_t expectedSize = static_cast<std::size_t>(tree.numInternalNodeSlots()) * static_cast<std::size_t>(attrNames.NUM_ATTRIBUTES);
    if (buffer.size() != expectedSize) {
        throw std::invalid_argument("Attribute output buffer size must match the dense internal-node domain and requested attributes.");
    }
}

/** @brief Validates that every requested scalar has an output column. */
template <std::floating_point Real>
inline void requireRequestedAttributeColumns(const AttributeComputeContext<Real>& context) {
    for (Attribute attribute : context.requestedAttributes) {
        if (!context.attrNames.contains(attribute)) {
            throw std::invalid_argument("Requested attribute computation requires a matching output column.");
        }
    }
}

template <std::floating_point Real, AltitudeValue T>
inline void requireRequestedAttributeColumns(const AltitudeAttributeComputeContext<Real, T>& context) {
    for (Attribute attribute : context.requestedAttributes) {
        if (!context.attrNames.contains(attribute)) {
            throw std::invalid_argument("Requested attribute computation requires a matching output column.");
        }
    }
}

/** @brief Finds a valid dependency source after the caller established dependency storage. */
template <std::floating_point Real>
inline const DependencySourceT<Real>* findDependencySource(std::span<const DependencySourceT<Real>> sources, Attribute attribute) {
    for (const DependencySourceT<Real>& source : sources) {
        if (source.attrNames->contains(attribute)) {
            return &source;
        }
    }
    return nullptr;
}

/**
 * @brief Tests whether `requestedAttributes` asks a kernel to fill `attribute`.
 *
 * @param requestedAttributes Attributes requested for materialization.
 * @param attribute Attribute requested by the operation.
 * @return True if requestedAttributes asks a kernel to fill attribute; otherwise false.
 */
inline bool requestsAttribute(std::span<const Attribute> requestedAttributes, Attribute attribute) {
    return std::find(requestedAttributes.begin(), requestedAttributes.end(), attribute) != requestedAttributes.end();
}

/**
 * @brief Validates a unit-attribute output buffer with one row per pixel.
 *
 * @details
 * Unit rows are used when projecting internal tree-node attributes to compact
 * exported Higra layouts. Each row corresponds to one pixel listed in
 * `unitPixels`, not to an internal tree node.
 *
 * @param tree Tree topology.
 * @param unitPixels Pixel identifiers in output row order.
 * @param buffer Buffer read or written by the operation.
 * @param attrNames Layout that maps attributes to buffer columns.
 */
template <std::floating_point Real>
inline void requireUnitAttributeBufferShape(const MorphologicalTree& tree, std::span<const PixelId> unitPixels, std::span<Real> buffer,
                                            const AttributeNames& attrNames) {
    const std::size_t expectedSize = unitPixels.size() * static_cast<std::size_t>(attrNames.NUM_ATTRIBUTES);
    MMCFILTERS_CONTRACT_REQUIRE(buffer.size() == expectedSize,
                                throw std::invalid_argument("Unit-attribute buffer size must match the exported leaf domain and requested attributes."));
    if constexpr (contract::validationsEnabled) {
        for (const PixelId pixel : unitPixels) {
            if (!tree.isPixel(pixel)) {
                throw std::invalid_argument("Unit-attribute computation requires valid pixel ids.");
            }
        }
    }
}

namespace detail::kernel {

/**
 * @brief Reads the altitude of an established pixel's smallest node.
 * @param tree Established tree topology.
 * @param altitudeView Established altitude span.
 * @param pixel Established pixel identifier.
 * @return Altitude of the smallest node.
 */
template <AltitudeValue T> inline T unitAltitude(const MorphologicalTree& tree, std::span<const T> altitudeView, PixelId pixel) {
    const NodeId smallestNodeId = detail::CommittedTreeAccess::smallestNodeMap(tree, pixel);
    return altitudeView[static_cast<std::size_t>(smallestNodeId)];
}

} // namespace detail::kernel

/**
 * @brief Reads the altitude of a pixel's smallest node.
 *
 * @param tree Tree topology.
 * @param altitudeView Altitude values indexed by internal node slot.
 * @param pixel Pixel identifier.
 * @return Altitude of the pixel's smallest node.
 */
template <AltitudeValue T> inline T unitAltitude(const MorphologicalTree& tree, std::span<const T> altitudeView, PixelId pixel) {
    TreeAltitudeAlgorithms::validateNodeAltitudeBufferShape(tree, altitudeView);
    MMCFILTERS_CONTRACT_REQUIRE(tree.isPixel(pixel),
                                throw std::invalid_argument("Unit altitude computation requires a valid proper-part id."));
    const NodeId smallestNodeId = detail::CommittedTreeAccess::smallestNodeMap(tree, pixel);
    if (smallestNodeId == InvalidNode || !detail::CommittedTreeAccess::isAlive(tree, smallestNodeId)) {
        throw std::runtime_error("Unit-attribute computation requires every proper part to have an alive smallest node.");
    }
    return altitudeView[static_cast<std::size_t>(smallestNodeId)];
}

} // namespace mmcfilters
