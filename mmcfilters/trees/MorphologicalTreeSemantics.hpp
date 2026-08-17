#pragma once

/**
 * @file MorphologicalTreeSemantics.hpp
 * @brief Immutable scientific metadata attached to a morphological tree.
 */

#include "../utils/Common.hpp"
#include "../utils/RegularGridAdjacency2D.hpp"

#include <utility>
#include <stdexcept>
#include <variant>

namespace mmcfilters {

/** @brief Declared construction family of one morphological tree. */
enum class MorphologicalTreeKind {
    Generic,
    MaxTree,
    MinTree,
    TreeOfShapes,
    UnrestrictedResidualTree,
    SaturatedResidualTree
};

/** @brief Global ordering constraint of node altitudes along parent-child arcs. */
enum class NodeAltitudeOrder {
    /// Every alive parent-child arc satisfies `nodeAltitude(parent) < nodeAltitude(child)`.
    Increasing,
    /// Every alive parent-child arc satisfies `nodeAltitude(parent) > nodeAltitude(child)`.
    Decreasing,
    /// No global parent-child altitude order is declared.
    Unconstrained
};

/** @brief Explicitly records that no construction context is available. */
struct NoConstructionContext {};

/** @brief Records one adjacency shared by both construction polarities. */
struct SharedAdjacencyContext {
    RegularGridAdjacency2D adjacency; ///< Shared construction adjacency.
};

/** @brief Records the shared adjacency and infinity pixel of a saturated residual construction. */
struct SaturatedResidualContext {
    RegularGridAdjacency2D adjacency; ///< Shared construction adjacency.
    PixelId infinityPixel = 0;        ///< Declared infinity pixel in the active domain.
};

/** @brief Ordered pair of complementary minimum and maximum adjacencies. */
struct ComplementaryAdjacencies {
    RegularGridAdjacency2D minAdjacency; ///< Adjacency used for minimum-oriented processing.
    RegularGridAdjacency2D maxAdjacency; ///< Adjacency used for maximum-oriented processing.
};

/** @brief Selects the self-dual span-valued tree-of-shapes immersion. */
struct SelfDualSpanImmersion {};

/** @brief Canonical ordered pairing of complementary 4- and 8-connectivities. */
enum class ComplementaryPairing { Min4Max8, Min8Max4 };

/**
 * @brief Selects a canonical complementary-grid immersion without binding a domain.
 *
 * The pairing is resolved into a `ComplementaryGridImmersion` carrying concrete
 * adjacencies once the source domain is known. Only the resolved form is
 * retained by a construction result.
 */
struct CanonicalComplementaryGridImmersion {
    ComplementaryPairing pairing = ComplementaryPairing::Min4Max8; ///< Canonical minimum/maximum pairing.
};

/** @brief Selects the optimized complementary-grid immersion. */
struct ComplementaryGridImmersion {
    ComplementaryAdjacencies complementaryAdjacencies; ///< Ordered minimum/maximum adjacency pair.
};

/** @brief Tagged tree-of-shapes immersion family. */
using TreeOfShapesImmersion = std::variant<CanonicalComplementaryGridImmersion, ComplementaryGridImmersion, SelfDualSpanImmersion>;

/** @brief Relation between the source pixel domain and active topographic domain. */
enum class TopographicDomainExtension { ExteriorRing, None };

/**
 * @brief Scale in which tree-of-shapes node altitudes are published.
 *
 * `UInt8` publishes source gray levels unchanged and is available only for
 * complementary-grid immersions, whose construction levels never leave the
 * source level set. `ExactDoubled` publishes doubled units, where a source
 * level `v` is represented by `2 * v` and odd values represent the half levels
 * that a self-dual span immersion can introduce.
 */
enum class TopographicAltitudeEncoding { UInt8, ExactDoubled };

/** @brief Complete discrete convention retained by a tree-of-shapes result. */
struct TopographicConvention {
    TreeOfShapesImmersion immersion = CanonicalComplementaryGridImmersion{}; ///< Selected topographic immersion.
    TopographicDomainExtension domainExtension = TopographicDomainExtension::ExteriorRing; ///< Active-domain extension convention.
    PixelId infinityPixel = 0; ///< Declared infinity pixel in the active topographic domain.
    TopographicAltitudeEncoding altitudeEncoding = TopographicAltitudeEncoding::UInt8; ///< Published node-altitude scale.
};

/**
 * @brief Builds a self-dual span convention in one call.
 *
 * The immersion, extension, and encoding of a self-dual span construction are
 * coupled: the exterior ring carries the boundary reference level, which admits
 * only `TopographicAltitudeEncoding::ExactDoubled`. Gathering them here keeps
 * call sites short and the combination coherent.
 *
 * @param domainExtension Active-domain extension convention.
 * @param infinityPixel Declared infinity pixel in the active topographic domain.
 * @param altitudeEncoding Published node-altitude scale.
 * @return A self-dual span convention carrying the requested fields.
 */
inline TopographicConvention selfDualSpanConvention(TopographicDomainExtension domainExtension = TopographicDomainExtension::ExteriorRing,
                                                    PixelId infinityPixel = 0,
                                                    TopographicAltitudeEncoding altitudeEncoding = TopographicAltitudeEncoding::ExactDoubled) {
    return TopographicConvention{SelfDualSpanImmersion{}, domainExtension, infinityPixel, altitudeEncoding};
}

/** @brief Tagged construction context retained by morphological-tree semantics. */
using MorphologicalTreeConstructionContext =
    std::variant<NoConstructionContext, SharedAdjacencyContext, SaturatedResidualContext, TopographicConvention>;

/**
 * @brief Immutable scientific interpretation attached to one morphological tree.
 *
 * The metadata do not alter topology. Algorithms query the typed construction
 * context they require and do not infer capabilities from `kind`.
 */
struct MorphologicalTreeSemantics {
    MorphologicalTreeKind kind = MorphologicalTreeKind::Generic; ///< Declared construction/result kind.
    NodeAltitudeOrder nodeAltitudeOrder = NodeAltitudeOrder::Unconstrained; ///< Global parent-child altitude order.
    MorphologicalTreeConstructionContext constructionContext = NoConstructionContext{}; ///< Typed retained construction context.
};

/** @brief Returns the conventional node-altitude order of a declared tree kind. @param kind Declared tree kind. @return Conventional altitude order. */
inline NodeAltitudeOrder defaultNodeAltitudeOrder(MorphologicalTreeKind kind) noexcept {
    switch (kind) {
    case MorphologicalTreeKind::MaxTree:
        return NodeAltitudeOrder::Increasing;
    case MorphologicalTreeKind::MinTree:
        return NodeAltitudeOrder::Decreasing;
    case MorphologicalTreeKind::Generic:
    case MorphologicalTreeKind::TreeOfShapes:
    case MorphologicalTreeKind::UnrestrictedResidualTree:
    case MorphologicalTreeKind::SaturatedResidualTree:
        return NodeAltitudeOrder::Unconstrained;
    }
    return NodeAltitudeOrder::Unconstrained;
}

/** @brief Constructs semantics using the conventional altitude order of `kind`. @param kind Declared tree kind. @param constructionContext Retained typed context. @return Coherent semantics descriptor. */
inline MorphologicalTreeSemantics makeMorphologicalTreeSemantics(
    MorphologicalTreeKind kind,
    MorphologicalTreeConstructionContext constructionContext = NoConstructionContext{}) {
    return MorphologicalTreeSemantics{kind, defaultNodeAltitudeOrder(kind), std::move(constructionContext)};
}

/**
 * @brief Validates the scientific compatibility of a semantics descriptor.
 *
 * `NoConstructionContext` is valid for every kind because it explicitly
 * records unavailable provenance. Every concrete context is restricted to the
 * construction families whose parameters it represents.
 * @param semantics Descriptor to validate.
 */
inline void validateMorphologicalTreeSemantics(const MorphologicalTreeSemantics& semantics) {
    if (semantics.kind != MorphologicalTreeKind::Generic && semantics.nodeAltitudeOrder != defaultNodeAltitudeOrder(semantics.kind)) {
        throw std::invalid_argument("Morphological-tree kind and node-altitude order are incompatible.");
    }
    if (std::holds_alternative<NoConstructionContext>(semantics.constructionContext)) {
        return;
    }
    if (std::holds_alternative<SharedAdjacencyContext>(semantics.constructionContext)) {
        if (semantics.kind == MorphologicalTreeKind::Generic || semantics.kind == MorphologicalTreeKind::MaxTree ||
            semantics.kind == MorphologicalTreeKind::MinTree || semantics.kind == MorphologicalTreeKind::UnrestrictedResidualTree) {
            return;
        }
        throw std::invalid_argument("SharedAdjacencyContext is incompatible with the declared morphological-tree kind.");
    }
    if (std::holds_alternative<SaturatedResidualContext>(semantics.constructionContext)) {
        if (semantics.kind == MorphologicalTreeKind::SaturatedResidualTree) {
            return;
        }
        throw std::invalid_argument("SaturatedResidualContext requires SaturatedResidualTree kind.");
    }
    if (std::holds_alternative<TopographicConvention>(semantics.constructionContext)) {
        if (semantics.kind == MorphologicalTreeKind::TreeOfShapes) {
            return;
        }
        throw std::invalid_argument("TopographicConvention requires TreeOfShapes kind.");
    }
    throw std::invalid_argument("Unsupported morphological-tree construction context.");
}

} // namespace mmcfilters
