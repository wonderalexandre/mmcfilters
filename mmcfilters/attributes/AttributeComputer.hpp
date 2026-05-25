#pragma once

#include <algorithm>
#include <cstddef>
#include <initializer_list>
#include <span>
#include <stdexcept>
#include <vector>

#include "../utils/Common.hpp"
#include "../trees/TreeAltitudeAlgorithms.hpp"
#include "../trees/WeightedMorphologicalTree.hpp"
#include "AttributeNames.hpp"


namespace mmcfilters {

/**
 * @brief Non-owning view over a dependency attribute buffer.
 *
 * @details
 * Attribute computers receive their upstream dependencies through lightweight
 * views rather than ownership-bearing objects. A `DependencySource` therefore
 * bundles:
 * - the attribute layout describing how to interpret the flat buffer;
 * - the raw float buffer holding the dependency values.
 *
 * The view is valid only for the duration of the computation that requested
 * it. Ownership remains with the attribute pipeline.
 */
struct DependencySource {
    /// Layout describing the dependency buffer; borrowed and non-owning.
    const AttributeNames* attrNames = nullptr;

    /// Flat dependency values indexed through `attrNames`; borrowed and non-owning.
    const float* buffer = nullptr;
};

/**
 * @brief Returns a validated dependency source for a required attribute.
 *
 * @details
 * Attribute computers and typed kernels both consume dependency buffers. This
 * helper centralises the common checks so every caller gets the same handling
 * for missing, null, or layout-incompatible dependency sources.
 */
inline const DependencySource& requireDependencySourceForAttribute(
    std::span<const DependencySource> dependencySources,
    size_t sourceIndex,
    Attribute requiredAttribute)
{
    if (sourceIndex >= dependencySources.size()) {
        throw std::invalid_argument("Requested attribute computation is missing a required dependency source.");
    }
    const DependencySource& source = dependencySources[sourceIndex];
    if (source.attrNames == nullptr || source.buffer == nullptr) {
        throw std::invalid_argument("Requested attribute computation received an invalid dependency source.");
    }
    if (!source.attrNames->contains(requiredAttribute)) {
        throw std::invalid_argument("Requested attribute computation received a dependency source without the required attribute.");
    }
    return source;
}

/**
 * @brief Returns a dependency source that contains all required attributes.
 */
inline const DependencySource& requireDependencySourceForAttributes(
    std::span<const DependencySource> dependencySources,
    size_t sourceIndex,
    std::initializer_list<Attribute> requiredAttributes)
{
    if (requiredAttributes.size() == 0) {
        throw std::invalid_argument("Requested attribute computation must name at least one required dependency attribute.");
    }

    const DependencySource& source = requireDependencySourceForAttribute(
        dependencySources,
        sourceIndex,
        *requiredAttributes.begin());
    for (const Attribute requiredAttribute : requiredAttributes) {
        if (!source.attrNames->contains(requiredAttribute)) {
            throw std::invalid_argument("Requested attribute computation received a dependency source without a required attribute.");
        }
    }
    return source;
}

/**
 * @brief Optional altitude view used inside the virtual attribute-computer API.
 *
 * An empty value means the request is genuinely topology-only. A present value
 * is a borrowed span over the caller-owned canonical `uint8_t` altitude buffer.
 * Keeping this as a view avoids forcing internal dispatch to traffic in nullable
 * owning-buffer pointers.
 */
using AttributeAltitudeView = OptionalAltitudeSpan<std::uint8_t>;

inline AttributeAltitudeView makeAttributeAltitudeView(
    const AltitudeBuffer<std::uint8_t>& altitude) noexcept
{
    return AltitudeSpan<std::uint8_t>(altitude);
}

/**
 * @brief Advanced interface for canonical attribute-computer implementations.
 *
 * @details
 * `AttributeComputer` is the advanced implementation contract used by concrete
 * attribute-family computers and implementation-level tests that need to run a
 * concrete computer directly. Ordinary attribute users should prefer
 * `AttributeComputation`, whose public facade accepts
 * `WeightedMorphologicalTree<T>` and `WeightedTreeView<T>` for typed altitude
 * computation and exposes explicit topology-only entry points for requests that
 * do not depend on altitude.
 *
 * Each concrete computer:
 * - declares the attributes it can naturally produce together;
 * - writes its output into a flat float buffer indexed through
 *   `AttributeNames`.
 *
 * The interface is intentionally generic enough to cover several patterns:
 * - single-attribute computers such as `AREA`;
 * - families of tightly related descriptors such as bounding-box attributes;
 * - grouped outputs whose internal computation shares intermediate state.
 *
 * Callers may request only a subset of the attributes supported by a computer.
 * Implementations are therefore expected to respect
 * `requestedAttributes` and materialise only what has been asked for, even if
 * the same traversal computes several related quantities.
 *
 * This interface intentionally remains tied to the canonical
 * `MorphologicalTree + AttributeAltitudeView` computer contract. New typed
 * public computation should be added through the typed pipeline or kernels in
 * `detail`, not by adding overloads to the ordinary public facade.
 */
class AttributeComputer {
public:
    virtual ~AttributeComputer() = default;

    /**
     * @brief Computes the attributes declared by `attributes()` into `buffer`.
     * @details This overload forwards to the fully parameterised version using
     * the computer's natural attribute set.
     */
    virtual void compute(const MorphologicalTree& tree, AttributeAltitudeView altitude, std::span<float> buffer, const AttributeNames& attrNames, std::span<const DependencySource> dependencySources = {}) const{
        compute(tree, altitude, buffer, attrNames, this->attributes(), dependencySources);
    }

    /**
     * @brief Topology-only convenience wrapper.
     */
    void compute(const MorphologicalTree& tree, std::span<float> buffer, const AttributeNames& attrNames, std::span<const DependencySource> dependencySources = {}) const{
        compute(tree, AttributeAltitudeView{}, buffer, attrNames, dependencySources);
    }

    /**
     * @brief Owner convenience wrapper using the tree's canonical altitude.
     */
    void compute(const WeightedMorphologicalTree<std::uint8_t>& tree, std::span<float> buffer, const AttributeNames& attrNames, std::span<const DependencySource> dependencySources = {}) const{
        compute(tree.topology(), makeAttributeAltitudeView(tree.getAltitudeBuffer()), buffer, attrNames, dependencySources);
    }

    /**
     * @brief Computes the requested attributes into `buffer`.
     * @param tree Tree on which the attribute is evaluated.
     * @param altitude Optional altitude data indexed by dense internal `NodeId`.
     * Empty views are accepted only by topology-only computers.
     * @param buffer Output buffer indexed through `attrNames`.
     * @param attrNames Attribute layout associated with `buffer`.
     * @param requestedAttributes Subset of the computer's natural attributes to
     * materialise in this invocation.
     * @param dependencySources Precomputed dependency buffers keyed by
     * attribute name and already validated by the caller.
     */
    virtual void compute(const MorphologicalTree& tree, AttributeAltitudeView altitude, std::span<float> buffer, const AttributeNames& attrNames, std::span<const Attribute> requestedAttributes, std::span<const DependencySource> dependencySources = {}) const = 0;

    /**
     * @brief Computes attribute values for unit components represented by
     * individual proper parts.
     *
     * @details
     * Compact Higra exports contain leaves/proper parts before internal nodes.
     * This hook lets each concrete computer define the value of its attributes
     * on a one-pixel component instead of relying on a generic placeholder.
     * The output buffer is indexed by the position of each proper part in
     * `unitProperParts`, not by the global proper-part id.
     */
    virtual void computeUnitAttributes(const MorphologicalTree& tree, AttributeAltitudeView altitude, std::span<const NodeId> unitProperParts, std::span<float> buffer, const AttributeNames& attrNames, std::span<const Attribute> requestedAttributes) const = 0;

    /**
     * @brief Topology-only convenience wrapper.
     */
    void compute(const MorphologicalTree& tree, std::span<float> buffer, const AttributeNames& attrNames, std::span<const Attribute> requestedAttributes, std::span<const DependencySource> dependencySources = {}) const{
        compute(tree, AttributeAltitudeView{}, buffer, attrNames, requestedAttributes, dependencySources);
    }

    /**
     * @brief Owner convenience wrapper using the tree's canonical altitude.
     */
    void compute(const WeightedMorphologicalTree<std::uint8_t>& tree, std::span<float> buffer, const AttributeNames& attrNames, std::span<const Attribute> requestedAttributes, std::span<const DependencySource> dependencySources = {}) const{
        compute(tree.topology(), makeAttributeAltitudeView(tree.getAltitudeBuffer()), buffer, attrNames, requestedAttributes, dependencySources);
    }

    /**
     * @brief Returns the attribute set naturally produced by this computer.
     */
    [[nodiscard]] virtual std::vector<Attribute> attributes() const = 0;

protected:

    /**
     * @brief Validates a dense internal-node output buffer.
     */
    static void requireAttributeBufferShape(const MorphologicalTree& tree, std::span<float> buffer, const AttributeNames& attrNames){
        const size_t expectedSize = static_cast<size_t>(tree.getNumInternalNodeSlots()) * static_cast<size_t>(attrNames.NUM_ATTRIBUTES);
        if (buffer.size() != expectedSize) {
            throw std::invalid_argument("Attribute output buffer size must match the dense internal-node domain and requested attributes.");
        }
    }

    /**
     * @brief Returns a validated dependency source for a required attribute.
     */
    static const DependencySource& requireDependencySource(std::span<const DependencySource> dependencySources, size_t sourceIndex, Attribute requiredAttribute){
        return requireDependencySourceForAttribute(dependencySources, sourceIndex, requiredAttribute);
    }

    /**
     * @brief Tests whether `requestedAttributes` asks the computer to fill
     * `attribute`.
     */
    static bool requestsAttribute(std::span<const Attribute> requestedAttributes, Attribute attribute) {
        return std::find(requestedAttributes.begin(), requestedAttributes.end(), attribute) != requestedAttributes.end();
    }

    /**
     * @brief Validates a proper-part-indexed unit-attribute output buffer.
     */
    static void requireUnitAttributeBufferShape(const MorphologicalTree& tree, std::span<const NodeId> unitProperParts, std::span<float> buffer, const AttributeNames& attrNames) {
        const size_t expectedSize = unitProperParts.size() * static_cast<size_t>(attrNames.NUM_ATTRIBUTES);
        if (buffer.size() != expectedSize) {
            throw std::invalid_argument("Unit-attribute buffer size must match the exported leaf domain and requested attributes.");
        }
        for (const NodeId properPart : unitProperParts) {
            if (!tree.isProperPart(properPart)) {
                throw std::invalid_argument("Unit-attribute computation requires valid proper-part ids.");
            }
        }
    }

    /**
     * @brief Reads the altitude assigned to the one-pixel component represented
     * by `properPart`.
     */
    static std::uint8_t unitAltitude(const MorphologicalTree& tree, AttributeAltitudeView altitude, NodeId properPart) {
        const std::span<const std::uint8_t> altitudeView = TreeAltitudeAlgorithms::requireAltitudeSpan(altitude);
        TreeAltitudeAlgorithms::validateAltitudeBufferShape(tree, altitudeView);
        if (!tree.isProperPart(properPart)) {
            throw std::invalid_argument("Unit altitude computation requires a valid proper-part id.");
        }
        const NodeId ownerNodeId = tree.getProperPartOwner(properPart);
        if (ownerNodeId == InvalidNode || !tree.isAlive(ownerNodeId)) {
            throw std::runtime_error("Unit-attribute computation requires every proper part to have an alive owner.");
        }
        return TreeAltitudeAlgorithms::getAltitude(altitudeView, ownerNodeId);
    }
};

} // namespace mmcfilters
