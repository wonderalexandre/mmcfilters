#pragma once

#include <span>

#include "../utils/Common.hpp"
#include "../trees/WeightedMorphologicalTree.hpp"
#include "AttributeNames.hpp"


namespace mmcfilters {

class AttributeNames;

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
 * it. Ownership remains with the incremental attribute pipeline.
 */
struct DependencySource {
    const AttributeNames* attrNames = nullptr;
    const float* buffer = nullptr;
};

/**
 * @brief Common interface for attribute computers operating on `MorphologicalTree`.
 *
 * @details
 * An `AttributeComputer` is the unit of work used by the incremental
 * attribute pipeline. Each computer:
 * - declares the attributes it can naturally produce together;
 * - declares the dependencies required from upstream computers;
 * - writes its output into a flat float buffer indexed through
 *   `AttributeNames`.
 *
 * The interface is intentionally generic enough to cover several patterns:
 * - single-attribute computers such as `AREA`;
 * - families of tightly related descriptors such as bounding-box attributes;
 * - grouped outputs whose internal computation shares intermediate state.
 *
 * The public pipeline may request only a subset of the attributes supported by
 * a computer. Implementations are therefore expected to respect
 * `requestedAttributes` and materialise only what has been asked for, even if
 * the same traversal computes several related quantities.
 *
 * See `AttributeComputedIncrementally.hpp` for the higher-level architectural
 * overview of how computers, layouts, the factory, and the dependency cache
 * fit together.
 */
class AttributeComputer {
public:
    virtual ~AttributeComputer() = default;

    /**
     * @brief Computes the attributes declared by `attributes()` into `buffer`.
     * @details This overload forwards to the fully parameterised version using
     * the computer's natural attribute set.
     */
    virtual void compute(
        MorphologicalTree& tree,
        const AltitudeBuffer* altitude,
        std::span<float> buffer,
        const AttributeNames& attrNames,
        std::span<const DependencySource> dependencySources = {}) const
    {
        compute(tree, altitude, buffer, attrNames, this->attributes(), dependencySources);
    }

    /**
     * @brief Compatibility wrapper for callers still passing `MorphologicalTree`.
     */
    void compute(
        MorphologicalTree& tree,
        std::span<float> buffer,
        const AttributeNames& attrNames,
        std::span<const DependencySource> dependencySources = {}) const
    {
        compute(tree, nullptr, buffer, attrNames, dependencySources);
    }

    /**
     * @brief Compatibility wrapper for callers passing `WeightedMorphologicalTree`.
     */
    void compute(
        WeightedMorphologicalTree& tree,
        std::span<float> buffer,
        const AttributeNames& attrNames,
        std::span<const DependencySource> dependencySources = {}) const
    {
        compute(tree.tree, &tree.altitude, buffer, attrNames, dependencySources);
    }

    /**
     * @brief Computes the requested attributes into `buffer`.
     * @param tree Tree on which the attribute is evaluated.
     * @param buffer Output buffer indexed through `attrNames`.
     * @param attrNames Attribute layout associated with `buffer`.
     * @param requestedAttributes Subset of the computer's natural attributes to
     * materialise in this invocation.
     * @param dependencySources Precomputed dependency buffers keyed by
     * attribute name and already validated by the caller.
     */
    virtual void compute(
        MorphologicalTree& tree,
        const AltitudeBuffer* altitude,
        std::span<float> buffer,
        const AttributeNames& attrNames,
        std::span<const Attribute> requestedAttributes,
        std::span<const DependencySource> dependencySources = {}) const = 0;

    /**
     * @brief Compatibility wrapper for callers still passing `MorphologicalTree`.
     */
    void compute(
        MorphologicalTree& tree,
        std::span<float> buffer,
        const AttributeNames& attrNames,
        std::span<const Attribute> requestedAttributes,
        std::span<const DependencySource> dependencySources = {}) const
    {
        compute(tree, nullptr, buffer, attrNames, requestedAttributes, dependencySources);
    }

    /**
     * @brief Compatibility wrapper for callers passing `WeightedMorphologicalTree`.
     */
    void compute(
        WeightedMorphologicalTree& tree,
        std::span<float> buffer,
        const AttributeNames& attrNames,
        std::span<const Attribute> requestedAttributes,
        std::span<const DependencySource> dependencySources = {}) const
    {
        compute(tree.tree, &tree.altitude, buffer, attrNames, requestedAttributes, dependencySources);
    }

    /**
     * @brief Returns the attribute set naturally produced by this computer.
     */
    virtual std::vector<Attribute> attributes() const = 0;

    /**
     * @brief Returns the dependency attributes required by this computer.
     *
     * @details
     * Dependencies may be declared either as individual attributes or as
     * attribute groups when the implementation expects a coherent family to be
     * materialised together.
     */
    virtual std::vector<AttributeOrGroup> requiredAttributes() const { return {}; }
};

} // namespace mmcfilters
