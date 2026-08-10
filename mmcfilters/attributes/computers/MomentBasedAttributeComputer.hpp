#pragma once

#include "AttributeComputerDomain.hpp"
#include "AttributeComputerFamily.hpp"
#include "../detail/AttributeKernelSupport.hpp"
#include "../../trees/detail/CommittedTreeAccess.hpp"
#include "../../trees/detail/TreeTraversalDetail.hpp"
#include "../../trees/MorphologicalTree.hpp"
#include "../../utils/Contract.hpp"

#include <array>
#include <concepts>
#include <cmath>
#include <limits>
#include <numbers>
#include <string_view>
#include <vector>

namespace mmcfilters::attributes::computers {

namespace detail {

/** @brief Selection mask for central-moment output columns. */
struct CentralMomentsRequest {
    bool mu20 = false; ///< Whether the central moment mu20 is requested.
    bool mu02 = false; ///< Whether the central moment mu02 is requested.
    bool mu11 = false; ///< Whether the central moment mu11 is requested.
    bool mu30 = false; ///< Whether the central moment mu30 is requested.
    bool mu03 = false; ///< Whether the central moment mu03 is requested.
    bool mu21 = false; ///< Whether the central moment mu21 is requested.
    bool mu12 = false; ///< Whether the central moment mu12 is requested.

    /** @brief Reports whether at least one central moment is requested. @return True when any request flag is set. */
    [[nodiscard]] bool any() const noexcept { return mu20 || mu02 || mu11 || mu30 || mu03 || mu21 || mu12; }

    /**
     * @brief Builds the selection mask from requested scalar attributes.
     * @param requestedAttributes Requested scalar attributes.
     * @return Central-moment selection mask.
     */
    [[nodiscard]] static CentralMomentsRequest from(std::span<const Attribute> requestedAttributes) {
        return {.mu20 = requestsAttribute(requestedAttributes, CENTRAL_MOMENT_20),
                .mu02 = requestsAttribute(requestedAttributes, CENTRAL_MOMENT_02),
                .mu11 = requestsAttribute(requestedAttributes, CENTRAL_MOMENT_11),
                .mu30 = requestsAttribute(requestedAttributes, CENTRAL_MOMENT_30),
                .mu03 = requestsAttribute(requestedAttributes, CENTRAL_MOMENT_03),
                .mu21 = requestsAttribute(requestedAttributes, CENTRAL_MOMENT_21),
                .mu12 = requestsAttribute(requestedAttributes, CENTRAL_MOMENT_12)};
    }
};

/** @brief Selection mask for Hu-moment output columns. */
struct HuMomentsRequest {
    bool hu1 = false; ///< Whether the first Hu invariant is requested.
    bool hu2 = false; ///< Whether the second Hu invariant is requested.
    bool hu3 = false; ///< Whether the third Hu invariant is requested.
    bool hu4 = false; ///< Whether the fourth Hu invariant is requested.
    bool hu5 = false; ///< Whether the fifth Hu invariant is requested.
    bool hu6 = false; ///< Whether the sixth Hu invariant is requested.
    bool hu7 = false; ///< Whether the seventh Hu invariant is requested.

    /** @brief Reports whether at least one Hu moment is requested. @return True when any request flag is set. */
    [[nodiscard]] bool any() const noexcept { return hu1 || hu2 || hu3 || hu4 || hu5 || hu6 || hu7; }

    /**
     * @brief Builds the selection mask from requested scalar attributes.
     * @param requestedAttributes Requested scalar attributes.
     * @return Hu-moment selection mask.
     */
    [[nodiscard]] static HuMomentsRequest from(std::span<const Attribute> requestedAttributes) {
        return {.hu1 = requestsAttribute(requestedAttributes, HU_MOMENT_1),
                .hu2 = requestsAttribute(requestedAttributes, HU_MOMENT_2),
                .hu3 = requestsAttribute(requestedAttributes, HU_MOMENT_3),
                .hu4 = requestsAttribute(requestedAttributes, HU_MOMENT_4),
                .hu5 = requestsAttribute(requestedAttributes, HU_MOMENT_5),
                .hu6 = requestsAttribute(requestedAttributes, HU_MOMENT_6),
                .hu7 = requestsAttribute(requestedAttributes, HU_MOMENT_7)};
    }
};

/** @brief Selection mask for attributes derived from second-order moments. */
struct MomentDerivedRequest {
    bool inertia = false;         ///< Whether normalized inertia is requested.
    bool compactness = false;     ///< Whether compactness is requested.
    bool eccentricity = false;    ///< Whether eccentricity is requested.
    bool majorAxis = false;       ///< Whether major-axis length is requested.
    bool minorAxis = false;       ///< Whether minor-axis length is requested.
    bool axisOrientation = false; ///< Whether principal-axis orientation is requested.
    bool circularity = false;     ///< Whether moment-based circularity is requested.

    /** @brief Reports whether at least one derived attribute is requested. @return True when any request flag is set. */
    [[nodiscard]] bool any() const noexcept {
        return inertia || compactness || eccentricity || majorAxis || minorAxis || axisOrientation || circularity;
    }

    /**
     * @brief Builds the selection mask from requested scalar attributes.
     * @param requestedAttributes Requested scalar attributes.
     * @return Moment-derived selection mask.
     */
    [[nodiscard]] static MomentDerivedRequest from(std::span<const Attribute> requestedAttributes) {
        return {.inertia = requestsAttribute(requestedAttributes, INERTIA),
                .compactness = requestsAttribute(requestedAttributes, COMPACTNESS),
                .eccentricity = requestsAttribute(requestedAttributes, ECCENTRICITY),
                .majorAxis = requestsAttribute(requestedAttributes, LENGTH_MAJOR_AXIS),
                .minorAxis = requestsAttribute(requestedAttributes, LENGTH_MINOR_AXIS),
                .axisOrientation = requestsAttribute(requestedAttributes, AXIS_ORIENTATION),
                .circularity = requestsAttribute(requestedAttributes, CIRCULARITY)};
    }
};

template <std::floating_point Real> inline constexpr Real maximumFiniteEccentricity() noexcept { return static_cast<Real>(1.0e6); }

namespace kernel {

/**
 * @brief Computes requested central moments over established tree supports.
 * @param context Established tree, output layout, and output buffer.
 * @param request Central-moment columns to materialize.
 */
template <std::floating_point Real>
inline void computeCentralMoments(const AttributeComputeContext<Real>& context, const CentralMomentsRequest& request) {
    if (!request.any()) {
        return;
    }

    struct RawMoments {
        Real m00 = Real{0};
        Real m10 = Real{0};
        Real m01 = Real{0};
        Real m20 = Real{0};
        Real m02 = Real{0};
        Real m11 = Real{0};
        Real m30 = Real{0};
        Real m03 = Real{0};
        Real m21 = Real{0};
        Real m12 = Real{0};

        void add(const RawMoments& other) {
            m00 += other.m00;
            m10 += other.m10;
            m01 += other.m01;
            m20 += other.m20;
            m02 += other.m02;
            m11 += other.m11;
            m30 += other.m30;
            m03 += other.m03;
            m21 += other.m21;
            m12 += other.m12;
        }
    };

    const int stride = context.attrNames.NUM_ATTRIBUTES;
    const auto offsetOf = [&](Attribute attribute) { return context.attrNames.indexMap.find(attribute)->second; };
    const int mu20Offset = request.mu20 ? offsetOf(CENTRAL_MOMENT_20) : 0;
    const int mu02Offset = request.mu02 ? offsetOf(CENTRAL_MOMENT_02) : 0;
    const int mu11Offset = request.mu11 ? offsetOf(CENTRAL_MOMENT_11) : 0;
    const int mu30Offset = request.mu30 ? offsetOf(CENTRAL_MOMENT_30) : 0;
    const int mu03Offset = request.mu03 ? offsetOf(CENTRAL_MOMENT_03) : 0;
    const int mu21Offset = request.mu21 ? offsetOf(CENTRAL_MOMENT_21) : 0;
    const int mu12Offset = request.mu12 ? offsetOf(CENTRAL_MOMENT_12) : 0;
    const auto outputIndex = [&](NodeId node, int offset) { return static_cast<std::size_t>(node * stride + offset); };

    const GridDomain2D& domain = ::mmcfilters::detail::CommittedTreeAccess::gridDomain2D(context.tree);
    std::vector<RawMoments> raw(static_cast<std::size_t>(context.tree.getNumInternalNodeSlots()));
    ::mmcfilters::detail::kernel::traversePostOrder(
        context.tree, context.tree.getRoot(),
        [&](NodeId node) {
            RawMoments& moments = raw[static_cast<std::size_t>(node)];
            moments = RawMoments{};
            for (int properPart : ::mmcfilters::detail::CommittedTreeAccess::properParts(context.tree, node)) {
                const Real x = static_cast<Real>(properPart % domain.cols);
                const Real y = static_cast<Real>(properPart / domain.cols);
                const Real x2 = x * x;
                const Real y2 = y * y;
                moments.m00 += Real{1};
                moments.m10 += x;
                moments.m01 += y;
                moments.m20 += x2;
                moments.m02 += y2;
                moments.m11 += x * y;
                moments.m30 += x2 * x;
                moments.m03 += y2 * y;
                moments.m21 += x2 * y;
                moments.m12 += x * y2;
            }
        },
        [&](NodeId parent, NodeId child) { raw[static_cast<std::size_t>(parent)].add(raw[static_cast<std::size_t>(child)]); },
        [&](NodeId node) {
            const RawMoments& moments = raw[static_cast<std::size_t>(node)];
            if (moments.m00 <= Real{0}) {
                return;
            }
            const Real cx = moments.m10 / moments.m00;
            const Real cy = moments.m01 / moments.m00;
            if (request.mu20)
                context.buffer[outputIndex(node, mu20Offset)] = moments.m20 - Real{2} * cx * moments.m10 + cx * cx * moments.m00;
            if (request.mu02)
                context.buffer[outputIndex(node, mu02Offset)] = moments.m02 - Real{2} * cy * moments.m01 + cy * cy * moments.m00;
            if (request.mu11)
                context.buffer[outputIndex(node, mu11Offset)] = moments.m11 - cx * moments.m01 - cy * moments.m10 + cx * cy * moments.m00;
            if (request.mu30)
                context.buffer[outputIndex(node, mu30Offset)] =
                    moments.m30 - Real{3} * cx * moments.m20 + Real{3} * cx * cx * moments.m10 - cx * cx * cx * moments.m00;
            if (request.mu03)
                context.buffer[outputIndex(node, mu03Offset)] =
                    moments.m03 - Real{3} * cy * moments.m02 + Real{3} * cy * cy * moments.m01 - cy * cy * cy * moments.m00;
            if (request.mu21)
                context.buffer[outputIndex(node, mu21Offset)] = moments.m21 - Real{2} * cx * moments.m11 - cy * moments.m20 +
                                                                Real{2} * cx * cy * moments.m10 + cx * cx * moments.m01 -
                                                                cx * cx * cy * moments.m00;
            if (request.mu12)
                context.buffer[outputIndex(node, mu12Offset)] = moments.m12 - Real{2} * cy * moments.m11 - cx * moments.m02 +
                                                                Real{2} * cx * cy * moments.m01 + cy * cy * moments.m10 -
                                                                cx * cy * cy * moments.m00;
        });
}

/**
 * @brief Computes requested Hu invariants from established area and central-moment dependencies.
 * @param context Established tree, output layout, and output buffer.
 * @param request Hu-invariant columns to materialize.
 * @param momentDependency Established central-moment dependency.
 * @param areaDependency Established area dependency.
 */
template <std::floating_point Real>
inline void computeHuMoments(const AttributeComputeContext<Real>& context, const HuMomentsRequest& request,
                             const DependencySourceT<Real>& momentDependency, const DependencySourceT<Real>& areaDependency) {
    if (!request.any()) {
        return;
    }

    const int stride = context.attrNames.NUM_ATTRIBUTES;
    const auto offsetOf = [&](Attribute attribute) { return context.attrNames.indexMap.find(attribute)->second; };
    const int hu1Offset = request.hu1 ? offsetOf(HU_MOMENT_1) : 0;
    const int hu2Offset = request.hu2 ? offsetOf(HU_MOMENT_2) : 0;
    const int hu3Offset = request.hu3 ? offsetOf(HU_MOMENT_3) : 0;
    const int hu4Offset = request.hu4 ? offsetOf(HU_MOMENT_4) : 0;
    const int hu5Offset = request.hu5 ? offsetOf(HU_MOMENT_5) : 0;
    const int hu6Offset = request.hu6 ? offsetOf(HU_MOMENT_6) : 0;
    const int hu7Offset = request.hu7 ? offsetOf(HU_MOMENT_7) : 0;
    const auto outputIndex = [&](NodeId node, int offset) { return static_cast<std::size_t>(node * stride + offset); };

    const int momentStride = momentDependency.attrNames->NUM_ATTRIBUTES;
    const auto momentOffset = [&](Attribute attribute) { return momentDependency.attrNames->indexMap.find(attribute)->second; };
    const auto momentIndex = [&](NodeId node, Attribute attribute) {
        return static_cast<std::size_t>(node * momentStride + momentOffset(attribute));
    };
    const int areaStride = areaDependency.attrNames->NUM_ATTRIBUTES;
    const int areaOffset = areaDependency.attrNames->indexMap.find(AREA)->second;
    const auto areaIndex = [&](NodeId node) { return static_cast<std::size_t>(node * areaStride + areaOffset); };
    const auto normMoment = [](Real area, Real moment, int p, int q) {
        return ::mmcfilters::attributes::numeric::safeDivide(moment, std::pow(area, static_cast<Real>(p + q + 2) / Real{2}));
    };

    ::mmcfilters::detail::kernel::traversePostOrder(
        context.tree, context.tree.getRoot(), [](NodeId) {}, [](NodeId, NodeId) {},
        [&](NodeId node) {
            const Real area = areaDependency.buffer[areaIndex(node)];
            if (area <= Real{0}) {
                return;
            }
            const Real eta20 = normMoment(area, momentDependency.buffer[momentIndex(node, CENTRAL_MOMENT_20)], 2, 0);
            const Real eta02 = normMoment(area, momentDependency.buffer[momentIndex(node, CENTRAL_MOMENT_02)], 0, 2);
            const Real eta11 = normMoment(area, momentDependency.buffer[momentIndex(node, CENTRAL_MOMENT_11)], 1, 1);
            const Real eta30 = normMoment(area, momentDependency.buffer[momentIndex(node, CENTRAL_MOMENT_30)], 3, 0);
            const Real eta03 = normMoment(area, momentDependency.buffer[momentIndex(node, CENTRAL_MOMENT_03)], 0, 3);
            const Real eta21 = normMoment(area, momentDependency.buffer[momentIndex(node, CENTRAL_MOMENT_21)], 2, 1);
            const Real eta12 = normMoment(area, momentDependency.buffer[momentIndex(node, CENTRAL_MOMENT_12)], 1, 2);

            if (request.hu1)
                context.buffer[outputIndex(node, hu1Offset)] = eta20 + eta02;
            if (request.hu2)
                context.buffer[outputIndex(node, hu2Offset)] = std::pow(eta20 - eta02, 2) + Real{4} * std::pow(eta11, 2);
            if (request.hu3)
                context.buffer[outputIndex(node, hu3Offset)] = std::pow(eta30 - Real{3} * eta12, 2) + std::pow(Real{3} * eta21 - eta03, 2);
            if (request.hu4)
                context.buffer[outputIndex(node, hu4Offset)] = std::pow(eta30 + eta12, 2) + std::pow(eta21 + eta03, 2);
            if (request.hu5)
                context.buffer[outputIndex(node, hu5Offset)] =
                    (eta30 - Real{3} * eta12) * (eta30 + eta12) * (std::pow(eta30 + eta12, 2) - Real{3} * std::pow(eta21 + eta03, 2)) +
                    (Real{3} * eta21 - eta03) * (eta21 + eta03) * (Real{3} * std::pow(eta30 + eta12, 2) - std::pow(eta21 + eta03, 2));
            if (request.hu6)
                context.buffer[outputIndex(node, hu6Offset)] =
                    (eta20 - eta02) * (std::pow(eta30 + eta12, 2) - std::pow(eta21 + eta03, 2)) +
                    Real{4} * eta11 * (eta30 + eta12) * (eta21 + eta03);
            if (request.hu7)
                context.buffer[outputIndex(node, hu7Offset)] =
                    (Real{3} * eta21 - eta03) * (eta30 + eta12) * (std::pow(eta30 + eta12, 2) - Real{3} * std::pow(eta21 + eta03, 2)) -
                    (eta30 - Real{3} * eta12) * (eta21 + eta03) * (Real{3} * std::pow(eta30 + eta12, 2) - std::pow(eta21 + eta03, 2));
        });
}

/**
 * @brief Computes requested shape descriptors from established moment dependencies.
 * @param context Established tree, output layout, and output buffer.
 * @param request Derived columns to materialize.
 * @param momentDependency Established second-order central moments.
 * @param areaDependency Established area dependency.
 */
template <std::floating_point Real>
inline void computeMomentDerived(const AttributeComputeContext<Real>& context, const MomentDerivedRequest& request,
                                 const DependencySourceT<Real>& momentDependency, const DependencySourceT<Real>& areaDependency) {
    if (!request.any()) {
        return;
    }

    const int stride = context.attrNames.NUM_ATTRIBUTES;
    const auto offsetOf = [&](Attribute attribute) { return context.attrNames.indexMap.find(attribute)->second; };
    const int inertiaOffset = request.inertia ? offsetOf(INERTIA) : 0;
    const int compactnessOffset = request.compactness ? offsetOf(COMPACTNESS) : 0;
    const int eccentricityOffset = request.eccentricity ? offsetOf(ECCENTRICITY) : 0;
    const int majorAxisOffset = request.majorAxis ? offsetOf(LENGTH_MAJOR_AXIS) : 0;
    const int minorAxisOffset = request.minorAxis ? offsetOf(LENGTH_MINOR_AXIS) : 0;
    const int orientationOffset = request.axisOrientation ? offsetOf(AXIS_ORIENTATION) : 0;
    const int circularityOffset = request.circularity ? offsetOf(CIRCULARITY) : 0;
    const auto outputIndex = [&](NodeId node, int offset) { return static_cast<std::size_t>(node * stride + offset); };

    const int momentStride = momentDependency.attrNames->NUM_ATTRIBUTES;
    const int mu20Offset = momentDependency.attrNames->indexMap.find(CENTRAL_MOMENT_20)->second;
    const int mu02Offset = momentDependency.attrNames->indexMap.find(CENTRAL_MOMENT_02)->second;
    const int mu11Offset = momentDependency.attrNames->indexMap.find(CENTRAL_MOMENT_11)->second;
    const auto momentIndex = [&](NodeId node, int offset) { return static_cast<std::size_t>(node * momentStride + offset); };
    const int areaStride = areaDependency.attrNames->NUM_ATTRIBUTES;
    const int areaOffset = areaDependency.attrNames->indexMap.find(AREA)->second;
    const auto areaIndex = [&](NodeId node) { return static_cast<std::size_t>(node * areaStride + areaOffset); };

    ::mmcfilters::detail::kernel::traversePostOrder(
        context.tree, context.tree.getRoot(), [](NodeId) {}, [](NodeId, NodeId) {},
        [&](NodeId node) {
            const Real mu20 = momentDependency.buffer[momentIndex(node, mu20Offset)];
            const Real mu02 = momentDependency.buffer[momentIndex(node, mu02Offset)];
            const Real mu11 = momentDependency.buffer[momentIndex(node, mu11Offset)];
            const Real area = areaDependency.buffer[areaIndex(node)];
            const Real discriminant = std::pow(mu20 - mu02, 2) + Real{4} * std::pow(mu11, 2);
            const Real sqrtDiscriminant = ::mmcfilters::attributes::numeric::safeSqrt(discriminant);
            const Real lambda1 = mu20 + mu02 + sqrtDiscriminant;
            const Real lambda2 = mu20 + mu02 - sqrtDiscriminant;

            if (request.majorAxis)
                context.buffer[outputIndex(node, majorAxisOffset)] =
                    ::mmcfilters::attributes::numeric::safeSqrt(::mmcfilters::attributes::numeric::safeDivide(Real{2} * lambda1, area));
            if (request.minorAxis)
                context.buffer[outputIndex(node, minorAxisOffset)] =
                    ::mmcfilters::attributes::numeric::safeSqrt(::mmcfilters::attributes::numeric::safeDivide(Real{2} * lambda2, area));
            if (request.eccentricity) {
                const Real eps = std::numeric_limits<Real>::epsilon();
                if (lambda1 <= eps && std::abs(lambda2) <= eps) {
                    context.buffer[outputIndex(node, eccentricityOffset)] = Real{1};
                } else if (lambda2 <= eps) {
                    context.buffer[outputIndex(node, eccentricityOffset)] = maximumFiniteEccentricity<Real>();
                } else {
                    context.buffer[outputIndex(node, eccentricityOffset)] = ::mmcfilters::attributes::numeric::clampUpper(
                        ::mmcfilters::attributes::numeric::safeDivide(lambda1, lambda2, maximumFiniteEccentricity<Real>()),
                        maximumFiniteEccentricity<Real>());
                }
            }
            if (request.compactness)
                context.buffer[outputIndex(node, compactnessOffset)] = (Real{1} / (Real{2} * std::numbers::pi_v<Real>)) *
                                                                       ::mmcfilters::attributes::numeric::safeDivide(area, mu20 + mu02);
            if (request.axisOrientation) {
                if (mu20 != mu02 || mu11 != Real{0}) {
                    const Real radians = Real{0.5} * std::atan2(Real{2} * mu11, mu20 - mu02);
                    const Real degrees = radians * (Real{180} / std::numbers::pi_v<Real>);
                    context.buffer[outputIndex(node, orientationOffset)] = std::fmod(std::abs(degrees), Real{360});
                } else {
                    context.buffer[outputIndex(node, orientationOffset)] = Real{0};
                }
            }
            if (request.inertia) {
                const Real areaSquared = area * area;
                context.buffer[outputIndex(node, inertiaOffset)] = ::mmcfilters::attributes::numeric::safeDivide(mu20, areaSquared) +
                                                                  ::mmcfilters::attributes::numeric::safeDivide(mu02, areaSquared);
            }
            if (request.circularity) {
                const Real eps = std::numeric_limits<Real>::epsilon();
                if (lambda1 <= eps && std::abs(lambda2) <= eps) {
                    context.buffer[outputIndex(node, circularityOffset)] = Real{1};
                } else if (lambda1 <= eps || lambda2 <= eps) {
                    context.buffer[outputIndex(node, circularityOffset)] = Real{0};
                } else {
                    context.buffer[outputIndex(node, circularityOffset)] = ::mmcfilters::attributes::numeric::safeDivide(lambda2, lambda1);
                }
            }
        });
}

} // namespace kernel

template <std::floating_point Real> inline void validateCentralMomentsContext(const AttributeComputeContext<Real>& context) {
    requireAttributeBufferShape(context.tree, context.buffer, context.attrNames);
    requireRequestedAttributeColumns(context);
    static_cast<void>(context.tree.requireGridDomain2D("CentralMomentsComputer"));
}

template <std::floating_point Real> inline void validateMomentDependencyContext(const AttributeComputeContext<Real>& context) {
    requireAttributeBufferShape(context.tree, context.buffer, context.attrNames);
    requireRequestedAttributeColumns(context);
}
} // namespace detail

/**
 * @brief Computes geometric central moments up to third order.
 *
 * @details
 * The computer exposes the classical central moments required by the rest of
 * the moment-based descriptors:
 * - second order: `mu20`, `mu02`, `mu11`;
 * - third order: `mu30`, `mu03`, `mu21`, `mu12`.
 *
 * The implementation is organised in two passes over the hierarchy:
 * 1. accumulate raw geometric moments over each subtree support;
 * 2. convert those raw moments to central moments at the node centroid.
 *
 * Raw moments are additive across disjoint supports, so child values can be
 * merged directly into their parent without the centroid-shift error that
 * occurs when summing already-centralised child moments.
 */
class CentralMomentsComputer {
  public:
    /// Family name used in dependency-plan diagnostics.
    static constexpr std::string_view familyName = "central-moments";

    /// Stable family id used by the scheduler.
    static constexpr AttributeComputerFamily family = AttributeComputerFamily::CentralMoments;

    /// Execution domain required by the computer.
    static constexpr AttributeComputerDomain domain = AttributeComputerDomain::Topology;

    /**
     * @brief Canonical list of central moments produced together.
     */
    inline static constexpr std::array<Attribute, 7> producedAttributes{CENTRAL_MOMENT_20, CENTRAL_MOMENT_02, CENTRAL_MOMENT_11, CENTRAL_MOMENT_30,
                                                                        CENTRAL_MOMENT_03, CENTRAL_MOMENT_21, CENTRAL_MOMENT_12};

    /**
     * @brief Computes the requested central moments.
     *
     * @details
     * The output buffer is indexed by dense internal node id and interpreted by
     * `context.attrNames`. Coordinates are row-major image coordinates
     * interpreted as `(x = col, y = row)`. The method accumulates raw moments
     * bottom-up and converts them to central moments after each subtree support
     * has been assembled.
     *
     * @param context Operation context or diagnostic label.
     */
    template <std::floating_point Real> static void compute(const AttributeComputeContext<Real>& context) {
        const detail::CentralMomentsRequest request = detail::CentralMomentsRequest::from(context.requestedAttributes);
        MMCFILTERS_CONTRACT_CHECKED_ONLY(detail::validateCentralMomentsContext(context));
        detail::kernel::computeCentralMoments(context, request);
    }

  public:
    /**
     * @brief Materializes central moments for one-pixel unit supports.
     *
     * All central moments of a one-pixel support are zero.
     *
     * @param context Operation context or diagnostic label.
     */
    template <std::floating_point Real> static void computeUnitRows(const UnitAttributeComputeContext<Real>& context) {
        requireUnitAttributeBufferShape(context.tree, context.unitProperParts, context.buffer, context.attrNames);
        constexpr std::array<Attribute, 7> zeroAttributes{CENTRAL_MOMENT_20, CENTRAL_MOMENT_02, CENTRAL_MOMENT_11, CENTRAL_MOMENT_30,
                                                          CENTRAL_MOMENT_03, CENTRAL_MOMENT_21, CENTRAL_MOMENT_12};
        for (const Attribute attribute : zeroAttributes) {
            if (!requestsAttribute(context.requestedAttributes, attribute)) {
                continue;
            }
            for (NodeId leafIndex = 0; leafIndex < static_cast<NodeId>(context.unitProperParts.size()); ++leafIndex) {
                context.buffer[context.attrNames.linearIndex(leafIndex, attribute)] = Real{0};
            }
        }
    }
};

/**
 * @brief Computes the seven Hu invariant moments.
 *
 * @details
 * Hu moments are derived from normalised central moments and provide
 * descriptors that are invariant to translation, scale, and rotation. This
 * computer expects the central moments and `AREA` to have been materialised by
 * upstream computers, then evaluates the seven classical Hu combinations for
 * each live node independently.
 *
 * @note The normalisation step is undefined for degenerate zero-area supports.
 * In practice such supports should not occur for live nodes in a valid tree.
 */
class HuMomentsComputer {
  public:
    /// Family name used in dependency-plan diagnostics.
    static constexpr std::string_view familyName = "hu-moments";

    /// Stable family id used by the scheduler.
    static constexpr AttributeComputerFamily family = AttributeComputerFamily::HuMoments;

    /// Execution domain required by the computer.
    static constexpr AttributeComputerDomain domain = AttributeComputerDomain::Topology;

    /**
     * @brief Canonical list of Hu moments produced together.
     */
    inline static constexpr std::array<Attribute, 7> producedAttributes{HU_MOMENT_1, HU_MOMENT_2, HU_MOMENT_3, HU_MOMENT_4,
                                                                        HU_MOMENT_5, HU_MOMENT_6, HU_MOMENT_7};

    /**
     * @brief Computes the requested Hu invariant moments.
     *
     * @details
     * Requires one dependency source containing all central moments and another
     * containing `AREA`. Sources are resolved semantically through
     * `context.dependencies`, so their order in the dependency span is not part
     * of the contract.
     *
     * @param context Operation context or diagnostic label.
     */
    template <std::floating_point Real> static void compute(const AttributeComputeContext<Real>& context) {
        const detail::HuMomentsRequest request = detail::HuMomentsRequest::from(context.requestedAttributes);
        MMCFILTERS_CONTRACT_CHECKED_ONLY(detail::validateMomentDependencyContext(context));
        if (!request.any()) {
            return;
        }

        const DependencySourceT<Real>* momentDependency = nullptr;
        const DependencySourceT<Real>* areaDependency = nullptr;
        if constexpr (contract::validationsEnabled) {
            momentDependency = &context.dependencies.requireAll({CENTRAL_MOMENT_20, CENTRAL_MOMENT_02, CENTRAL_MOMENT_11, CENTRAL_MOMENT_30,
                                                                 CENTRAL_MOMENT_03, CENTRAL_MOMENT_21, CENTRAL_MOMENT_12});
            areaDependency = &context.dependencies.require(AREA);
        } else {
            momentDependency = ::mmcfilters::findDependencySource(context.dependencySources, CENTRAL_MOMENT_20);
            areaDependency = ::mmcfilters::findDependencySource(context.dependencySources, AREA);
        }
        detail::kernel::computeHuMoments(context, request, *momentDependency, *areaDependency);
    }

  public:
    /**
     * @brief Materializes Hu moments for one-pixel unit supports.
     *
     * The degenerate one-pixel support uses zero for all Hu invariants in the
     * exported unit path.
     *
     * @param context Operation context or diagnostic label.
     */
    template <std::floating_point Real> static void computeUnitRows(const UnitAttributeComputeContext<Real>& context) {
        requireUnitAttributeBufferShape(context.tree, context.unitProperParts, context.buffer, context.attrNames);
        constexpr std::array<Attribute, 7> zeroAttributes{HU_MOMENT_1, HU_MOMENT_2, HU_MOMENT_3, HU_MOMENT_4, HU_MOMENT_5, HU_MOMENT_6, HU_MOMENT_7};
        for (const Attribute attribute : zeroAttributes) {
            if (!requestsAttribute(context.requestedAttributes, attribute)) {
                continue;
            }
            for (NodeId leafIndex = 0; leafIndex < static_cast<NodeId>(context.unitProperParts.size()); ++leafIndex) {
                context.buffer[context.attrNames.linearIndex(leafIndex, attribute)] = Real{0};
            }
        }
    }
};

/**
 * @brief Computes higher-level shape descriptors derived from second-order
 * central moments.
 *
 * @details
 * The descriptors produced here summarise the geometry of the equivalent
 * second-moment ellipse associated with each node support. They include:
 * - lengths of the major and minor axes;
 * - eccentricity and circularity-like ratios derived from the eigenvalues of
 *   the inertia matrix;
 * - compactness and inertia;
 * - axis orientation in degrees.
 *
 * The computer consumes `mu20`, `mu02`, `mu11`, and `AREA`. From those values
 * it builds the 2x2 second-moment matrix, extracts its eigenvalues, and uses
 * them to derive the final descriptors.
 *
 * @note All computations are local once the dependencies have been produced.
 * No additional tree traversal state is required beyond iterating over the
 * live nodes in post-order.
 */
class MomentBasedAttributeComputer {
  public:
    /// Family name used in dependency-plan diagnostics.
    static constexpr std::string_view familyName = "moment-derived";

    /// Stable family id used by the scheduler.
    static constexpr AttributeComputerFamily family = AttributeComputerFamily::MomentDerived;

    /// Execution domain required by the computer.
    static constexpr AttributeComputerDomain domain = AttributeComputerDomain::Topology;

    /**
     * @brief Largest finite eccentricity emitted for degenerate one-dimensional supports.
     *
     * @return Largest finite eccentricity emitted for degenerate one-dimensional supports.
     */
    template <std::floating_point Real> static constexpr Real maxFiniteEccentricity() noexcept {
        return detail::maximumFiniteEccentricity<Real>();
    }

    /**
     * @brief Canonical list of moment-derived descriptors produced by this computer.
     */
    inline static constexpr std::array<Attribute, 7> producedAttributes{INERTIA,           COMPACTNESS,      ECCENTRICITY, LENGTH_MAJOR_AXIS,
                                                                        LENGTH_MINOR_AXIS, AXIS_ORIENTATION, CIRCULARITY};

    /**
     * @brief Computes the requested moment-derived descriptors.
     *
     * @details
     * Requires dependencies for `CENTRAL_MOMENT_20`, `CENTRAL_MOMENT_02`,
     * `CENTRAL_MOMENT_11`, and `AREA`. The dependencies are resolved by name,
     * then each live node is evaluated independently from the second-order
     * moment matrix.
     *
     * @param context Operation context or diagnostic label.
     */
    template <std::floating_point Real> static void compute(const AttributeComputeContext<Real>& context) {
        const detail::MomentDerivedRequest request = detail::MomentDerivedRequest::from(context.requestedAttributes);
        MMCFILTERS_CONTRACT_CHECKED_ONLY(detail::validateMomentDependencyContext(context));
        if (!request.any()) {
            return;
        }

        const DependencySourceT<Real>* momentDependency = nullptr;
        const DependencySourceT<Real>* areaDependency = nullptr;
        if constexpr (contract::validationsEnabled) {
            momentDependency = &context.dependencies.requireAll({CENTRAL_MOMENT_20, CENTRAL_MOMENT_02, CENTRAL_MOMENT_11});
            areaDependency = &context.dependencies.require(AREA);
        } else {
            momentDependency = ::mmcfilters::findDependencySource(context.dependencySources, CENTRAL_MOMENT_20);
            areaDependency = ::mmcfilters::findDependencySource(context.dependencySources, AREA);
        }
        detail::kernel::computeMomentDerived(context, request, *momentDependency, *areaDependency);
    }

  public:
    /**
     * @brief Materializes moment-derived descriptors for one-pixel unit supports.
     *
     * Eccentricity and circularity are defined as `1` for the degenerate unit
     * support; other descriptors are zero.
     *
     * @param context Operation context or diagnostic label.
     */
    template <std::floating_point Real> static void computeUnitRows(const UnitAttributeComputeContext<Real>& context) {
        requireUnitAttributeBufferShape(context.tree, context.unitProperParts, context.buffer, context.attrNames);
        constexpr std::array<Attribute, 7> zeroAttributes{COMPACTNESS,      ECCENTRICITY, LENGTH_MAJOR_AXIS, LENGTH_MINOR_AXIS,
                                                          AXIS_ORIENTATION, INERTIA,      CIRCULARITY};
        for (const Attribute attribute : zeroAttributes) {
            if (!requestsAttribute(context.requestedAttributes, attribute)) {
                continue;
            }
            for (NodeId leafIndex = 0; leafIndex < static_cast<NodeId>(context.unitProperParts.size()); ++leafIndex) {
                context.buffer[context.attrNames.linearIndex(leafIndex, attribute)] =
                    (attribute == ECCENTRICITY || attribute == CIRCULARITY) ? Real{1} : Real{0};
            }
        }
    }
};

} // namespace mmcfilters::attributes::computers
