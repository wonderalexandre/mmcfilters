#pragma once

#include "AttributeComputerDomain.hpp"
#include "AttributeComputerFamily.hpp"
#include "../detail/AttributeKernelSupport.hpp"
#include "../../trees/detail/TreeTraversalDetail.hpp"
#include "../../trees/MorphologicalTree.hpp"
#include "../../utils/Image.hpp"

#include <algorithm>
#include <array>
#include <concepts>
#include <cmath>
#include <limits>
#include <numbers>
#include <string_view>
#include <vector>

namespace mmcfilters::attributes::computers {

namespace detail {
/**
 * @brief Returns the dense attribute-buffer slot for a tree node.
 *
 * @param nodeId Identifier of the node used by the operation.
 * @return Dense attribute-buffer slot for the node.
 */
inline NodeId momentSlotOf(const MorphologicalTree&, NodeId nodeId) noexcept { return nodeId; }
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
        computeImpl(context.tree, context.buffer, context.attrNames, context.requestedAttributes);
    }

  private:
    /**
     * @brief Computes the requested attribute values into the output buffer.
     *
     * @param tree Tree topology used by the operation.
     * @param buffer Buffer read or written by the operation.
     * @param attrNames Layout mapping attributes to buffer columns.
     * @param requested Requested attribute subset.
     */
    template <std::floating_point Real>
    static void computeImpl(const MorphologicalTree& tree, std::span<Real> buffer, const AttributeNames& attrNames, std::span<const Attribute> requested) {
        requireAttributeBufferShape(tree, buffer, attrNames);

        bool computeMu20 = std::find(requested.begin(), requested.end(), CENTRAL_MOMENT_20) != requested.end();
        bool computeMu02 = std::find(requested.begin(), requested.end(), CENTRAL_MOMENT_02) != requested.end();
        bool computeMu11 = std::find(requested.begin(), requested.end(), CENTRAL_MOMENT_11) != requested.end();
        bool computeMu30 = std::find(requested.begin(), requested.end(), CENTRAL_MOMENT_30) != requested.end();
        bool computeMu03 = std::find(requested.begin(), requested.end(), CENTRAL_MOMENT_03) != requested.end();
        bool computeMu21 = std::find(requested.begin(), requested.end(), CENTRAL_MOMENT_21) != requested.end();
        bool computeMu12 = std::find(requested.begin(), requested.end(), CENTRAL_MOMENT_12) != requested.end();
        if (!computeMu20 && !computeMu02 && !computeMu11 && !computeMu30 && !computeMu03 && !computeMu21 && !computeMu12) {
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

        const int numCols = tree.getNumColsOfGridDomain2D();
        std::vector<RawMoments> raw(static_cast<std::size_t>(tree.getNumInternalNodeSlots()));
        auto indexOf = [&](NodeId idx, Attribute attr) { return attrNames.linearIndex(idx, attr); };

        ::mmcfilters::detail::traversePostOrder(
            tree, tree.getRoot(),
            [&](NodeId nodeId) {
                const NodeId node = detail::momentSlotOf(tree, nodeId);
                raw[static_cast<std::size_t>(node)] = RawMoments{};
                for (int p : tree.getProperParts(nodeId)) {
                    auto [py, px] = ImageUtils::to2D(p, numCols);
                    const Real x = static_cast<Real>(px);
                    const Real y = static_cast<Real>(py);
                    const Real x2 = x * x;
                    const Real y2 = y * y;
                    RawMoments& moments = raw[static_cast<std::size_t>(node)];
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
            [&](NodeId parentNodeId, NodeId childNodeId) {
                const NodeId parent = detail::momentSlotOf(tree, parentNodeId);
                const NodeId child = detail::momentSlotOf(tree, childNodeId);
                raw[static_cast<std::size_t>(parent)].add(raw[static_cast<std::size_t>(child)]);
            },
            [&](NodeId nodeId) {
                const NodeId node = detail::momentSlotOf(tree, nodeId);
                const RawMoments& moments = raw[static_cast<std::size_t>(node)];
                if (moments.m00 <= Real{0}) {
                    return;
                }

                const Real cx = moments.m10 / moments.m00;
                const Real cy = moments.m01 / moments.m00;

                if (computeMu20) {
                    const Real mu20 = moments.m20 - Real{2} * cx * moments.m10 + cx * cx * moments.m00;
                    buffer[indexOf(node, CENTRAL_MOMENT_20)] = mu20;
                }
                if (computeMu02) {
                    const Real mu02 = moments.m02 - Real{2} * cy * moments.m01 + cy * cy * moments.m00;
                    buffer[indexOf(node, CENTRAL_MOMENT_02)] = mu02;
                }
                if (computeMu11) {
                    const Real mu11 = moments.m11 - cx * moments.m01 - cy * moments.m10 + cx * cy * moments.m00;
                    buffer[indexOf(node, CENTRAL_MOMENT_11)] = mu11;
                }
                if (computeMu30) {
                    const Real mu30 = moments.m30 - Real{3} * cx * moments.m20 + Real{3} * cx * cx * moments.m10 - cx * cx * cx * moments.m00;
                    buffer[indexOf(node, CENTRAL_MOMENT_30)] = mu30;
                }
                if (computeMu03) {
                    const Real mu03 = moments.m03 - Real{3} * cy * moments.m02 + Real{3} * cy * cy * moments.m01 - cy * cy * cy * moments.m00;
                    buffer[indexOf(node, CENTRAL_MOMENT_03)] = mu03;
                }
                if (computeMu21) {
                    const Real mu21 = moments.m21 - Real{2} * cx * moments.m11 - cy * moments.m20 + Real{2} * cx * cy * moments.m10 + cx * cx * moments.m01 -
                                      cx * cx * cy * moments.m00;
                    buffer[indexOf(node, CENTRAL_MOMENT_21)] = mu21;
                }
                if (computeMu12) {
                    const Real mu12 = moments.m12 - Real{2} * cy * moments.m11 - cx * moments.m02 + Real{2} * cx * cy * moments.m01 + cy * cy * moments.m10 -
                                      cx * cy * cy * moments.m00;
                    buffer[indexOf(node, CENTRAL_MOMENT_12)] = mu12;
                }
            });
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
        computeImpl(context.tree, context.buffer, context.attrNames, context.requestedAttributes, context.dependencySources);
    }

  private:
    /**
     * @brief Computes the requested attribute values into the output buffer.
     *
     * @param tree Tree topology used by the operation.
     * @param buffer Buffer read or written by the operation.
     * @param attrNames Layout mapping attributes to buffer columns.
     * @param requested Requested attribute subset.
     * @param dependencySources Available dependency-attribute sources.
     */
    template <std::floating_point Real>
    static void computeImpl(const MorphologicalTree& tree, std::span<Real> buffer, const AttributeNames& attrNames, std::span<const Attribute> requested,
                            std::span<const DependencySourceT<Real>> dependencySources) {
        requireAttributeBufferShape(tree, buffer, attrNames);

        auto indexOf = [&](NodeId idx, Attribute attr) { return attrNames.linearIndex(idx, attr); };
        const DependencyResolver<Real> dependencies{dependencySources};
        const auto& muDependency = dependencies.requireAll({
            CENTRAL_MOMENT_20,
            CENTRAL_MOMENT_02,
            CENTRAL_MOMENT_11,
            CENTRAL_MOMENT_30,
            CENTRAL_MOMENT_03,
            CENTRAL_MOMENT_21,
            CENTRAL_MOMENT_12,
        });
        const auto& areaDependency = dependencies.require(AREA);
        auto indexOfMu = [&](NodeId idx, Attribute attr) { return muDependency.attrNames->linearIndex(idx, attr); };
        auto indexOfArea = [&](NodeId idx) { return areaDependency.attrNames->linearIndex(idx, AREA); };

        auto normMoment = [](Real area, Real moment, int p, int q) {
            return ::mmcfilters::attributes::numeric::safeDivide(moment, std::pow(area, static_cast<Real>(p + q + 2) / Real{2}));
        };

        bool computeHu1 = std::find(requested.begin(), requested.end(), HU_MOMENT_1) != requested.end();
        bool computeHu2 = std::find(requested.begin(), requested.end(), HU_MOMENT_2) != requested.end();
        bool computeHu3 = std::find(requested.begin(), requested.end(), HU_MOMENT_3) != requested.end();
        bool computeHu4 = std::find(requested.begin(), requested.end(), HU_MOMENT_4) != requested.end();
        bool computeHu5 = std::find(requested.begin(), requested.end(), HU_MOMENT_5) != requested.end();
        bool computeHu6 = std::find(requested.begin(), requested.end(), HU_MOMENT_6) != requested.end();
        bool computeHu7 = std::find(requested.begin(), requested.end(), HU_MOMENT_7) != requested.end();

        ::mmcfilters::detail::traversePostOrder(
            tree, tree.getRoot(), [](NodeId) {}, [](NodeId, NodeId) {},
            [&](NodeId idxGlobalId) {
                const NodeId idx = detail::momentSlotOf(tree, idxGlobalId);
                Real mu20 = muDependency.buffer[indexOfMu(idx, CENTRAL_MOMENT_20)];
                Real mu02 = muDependency.buffer[indexOfMu(idx, CENTRAL_MOMENT_02)];
                Real mu11 = muDependency.buffer[indexOfMu(idx, CENTRAL_MOMENT_11)];
                Real mu30 = muDependency.buffer[indexOfMu(idx, CENTRAL_MOMENT_30)];
                Real mu03 = muDependency.buffer[indexOfMu(idx, CENTRAL_MOMENT_03)];
                Real mu21 = muDependency.buffer[indexOfMu(idx, CENTRAL_MOMENT_21)];
                Real mu12 = muDependency.buffer[indexOfMu(idx, CENTRAL_MOMENT_12)];
                Real area = areaDependency.buffer[indexOfArea(idx)];
                if (area <= Real{0}) {
                    return;
                }

                // Normalise the central moments before evaluating the Hu invariants.
                Real eta20 = normMoment(area, mu20, 2, 0);
                Real eta02 = normMoment(area, mu02, 0, 2);
                Real eta11 = normMoment(area, mu11, 1, 1);
                Real eta30 = normMoment(area, mu30, 3, 0);
                Real eta03 = normMoment(area, mu03, 0, 3);
                Real eta21 = normMoment(area, mu21, 2, 1);
                Real eta12 = normMoment(area, mu12, 1, 2);

                // Evaluate the seven classical Hu combinations.
                if (computeHu1)
                    buffer[indexOf(idx, HU_MOMENT_1)] = eta20 + eta02; // The first Hu moment corresponds to a normalised inertia term.
                if (computeHu2)
                    buffer[indexOf(idx, HU_MOMENT_2)] = std::pow(eta20 - eta02, 2) + Real{4} * std::pow(eta11, 2);
                if (computeHu3)
                    buffer[indexOf(idx, HU_MOMENT_3)] = std::pow(eta30 - Real{3} * eta12, 2) + std::pow(Real{3} * eta21 - eta03, 2);
                if (computeHu4)
                    buffer[indexOf(idx, HU_MOMENT_4)] = std::pow(eta30 + eta12, 2) + std::pow(eta21 + eta03, 2);
                if (computeHu5)
                    buffer[indexOf(idx, HU_MOMENT_5)] =
                        (eta30 - Real{3} * eta12) * (eta30 + eta12) * (std::pow(eta30 + eta12, 2) - Real{3} * std::pow(eta21 + eta03, 2)) +
                        (Real{3} * eta21 - eta03) * (eta21 + eta03) * (Real{3} * std::pow(eta30 + eta12, 2) - std::pow(eta21 + eta03, 2));
                if (computeHu6)
                    buffer[indexOf(idx, HU_MOMENT_6)] =
                        (eta20 - eta02) * (std::pow(eta30 + eta12, 2) - std::pow(eta21 + eta03, 2)) + Real{4} * eta11 * (eta30 + eta12) * (eta21 + eta03);
                if (computeHu7)
                    buffer[indexOf(idx, HU_MOMENT_7)] =
                        (Real{3} * eta21 - eta03) * (eta30 + eta12) * (std::pow(eta30 + eta12, 2) - Real{3} * std::pow(eta21 + eta03, 2)) -
                        (eta30 - Real{3} * eta12) * (eta21 + eta03) * (Real{3} * std::pow(eta30 + eta12, 2) - std::pow(eta21 + eta03, 2));
            });
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
    template <std::floating_point Real> static constexpr Real maxFiniteEccentricity() noexcept { return static_cast<Real>(1.0e6); }

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
        computeImpl(context.tree, context.buffer, context.attrNames, context.requestedAttributes, context.dependencySources);
    }

  private:
    /**
     * @brief Computes the requested attribute values into the output buffer.
     *
     * @param tree Tree topology used by the operation.
     * @param buffer Buffer read or written by the operation.
     * @param attrNames Layout mapping attributes to buffer columns.
     * @param requestedAttributes Requested attribute subset.
     * @param dependencySources Available dependency-attribute sources.
     */
    template <std::floating_point Real>
    static void computeImpl(const MorphologicalTree& tree, std::span<Real> buffer, const AttributeNames& attrNames,
                            std::span<const Attribute> requestedAttributes, std::span<const DependencySourceT<Real>> dependencySources) {
        requireAttributeBufferShape(tree, buffer, attrNames);

        auto indexOfMajorAxis = [&](int idx) { return attrNames.linearIndex(idx, LENGTH_MAJOR_AXIS); };
        auto indexOfMinorAxis = [&](int idx) { return attrNames.linearIndex(idx, LENGTH_MINOR_AXIS); };
        auto indexOfEccentricity = [&](int idx) { return attrNames.linearIndex(idx, ECCENTRICITY); };
        auto indexOfCompactness = [&](int idx) { return attrNames.linearIndex(idx, COMPACTNESS); };
        auto indexOfAxisOrientation = [&](int idx) { return attrNames.linearIndex(idx, AXIS_ORIENTATION); };
        auto indexOfInertia = [&](int idx) { return attrNames.linearIndex(idx, INERTIA); };
        auto indexOfCircularity = [&](int idx) { return attrNames.linearIndex(idx, CIRCULARITY); };

        bool computeMajorAxis = std::find(requestedAttributes.begin(), requestedAttributes.end(), LENGTH_MAJOR_AXIS) != requestedAttributes.end();
        bool computeMinorAxis = std::find(requestedAttributes.begin(), requestedAttributes.end(), LENGTH_MINOR_AXIS) != requestedAttributes.end();
        bool computeEccentricity = std::find(requestedAttributes.begin(), requestedAttributes.end(), ECCENTRICITY) != requestedAttributes.end();
        bool computeCompactness = std::find(requestedAttributes.begin(), requestedAttributes.end(), COMPACTNESS) != requestedAttributes.end();
        bool computeAxisOrientation = std::find(requestedAttributes.begin(), requestedAttributes.end(), AXIS_ORIENTATION) != requestedAttributes.end();
        bool computeInertia = std::find(requestedAttributes.begin(), requestedAttributes.end(), INERTIA) != requestedAttributes.end();
        bool computeCircularity = std::find(requestedAttributes.begin(), requestedAttributes.end(), CIRCULARITY) != requestedAttributes.end();

        const DependencyResolver<Real> dependencies{dependencySources};
        const auto& momentDependency = dependencies.requireAll({
            CENTRAL_MOMENT_20,
            CENTRAL_MOMENT_02,
            CENTRAL_MOMENT_11,
        });
        const auto& areaDependency = dependencies.require(AREA);
        auto indexMu20 = [&](int idx) { return momentDependency.attrNames->linearIndex(idx, CENTRAL_MOMENT_20); };
        auto indexMu02 = [&](int idx) { return momentDependency.attrNames->linearIndex(idx, CENTRAL_MOMENT_02); };
        auto indexMu11 = [&](int idx) { return momentDependency.attrNames->linearIndex(idx, CENTRAL_MOMENT_11); };
        auto indexArea = [&](int idx) { return areaDependency.attrNames->linearIndex(idx, AREA); };

        ::mmcfilters::detail::traversePostOrder(
            tree, tree.getRoot(), [&](NodeId) {}, [&](NodeId, NodeId) {},
            [&](NodeId idxGlobalId) {
                const NodeId idx = detail::momentSlotOf(tree, idxGlobalId);
                Real mu20 = momentDependency.buffer[indexMu20(idx)];
                Real mu02 = momentDependency.buffer[indexMu02(idx)];
                Real mu11 = momentDependency.buffer[indexMu11(idx)];
                Real area = areaDependency.buffer[indexArea(idx)];

                const Real discriminant = std::pow(mu20 - mu02, 2) + Real{4} * std::pow(mu11, 2);
                const Real sqrtDiscriminant = ::mmcfilters::attributes::numeric::safeSqrt(discriminant);
                const Real lambda1 = mu20 + mu02 + sqrtDiscriminant; // Largest eigenvalue of the inertia matrix.
                const Real lambda2 = mu20 + mu02 - sqrtDiscriminant; // Smallest eigenvalue of the inertia matrix.

                if (computeMajorAxis) {
                    buffer[indexOfMajorAxis(idx)] =
                        ::mmcfilters::attributes::numeric::safeSqrt(::mmcfilters::attributes::numeric::safeDivide(Real{2} * lambda1, area));
                }
                if (computeMinorAxis) {
                    buffer[indexOfMinorAxis(idx)] =
                        ::mmcfilters::attributes::numeric::safeSqrt(::mmcfilters::attributes::numeric::safeDivide(Real{2} * lambda2, area));
                }
                if (computeEccentricity) {
                    const Real eps = std::numeric_limits<Real>::epsilon();
                    if (lambda1 <= eps && std::abs(lambda2) <= eps) {
                        buffer[indexOfEccentricity(idx)] = Real{1};
                    } else if (lambda2 <= eps) {
                        buffer[indexOfEccentricity(idx)] = maxFiniteEccentricity<Real>();
                    } else {
                        buffer[indexOfEccentricity(idx)] = ::mmcfilters::attributes::numeric::clampUpper(
                            ::mmcfilters::attributes::numeric::safeDivide(lambda1, lambda2, maxFiniteEccentricity<Real>()), maxFiniteEccentricity<Real>());
                    }
                }
                if (computeCompactness) {
                    Real denom = mu20 + mu02;
                    buffer[indexOfCompactness(idx)] =
                        (Real{1} / (Real{2} * std::numbers::pi_v<Real>)) * ::mmcfilters::attributes::numeric::safeDivide(area, denom);
                }
                if (computeAxisOrientation) {
                    // Guard the orientation computation against the fully isotropic case.
                    if (mu20 != mu02 || mu11 != 0) {
                        Real radians = Real{0.5} * std::atan2(Real{2} * mu11, mu20 - mu02);            // Orientation in radians.
                        Real degrees = radians * (Real{180} / std::numbers::pi_v<Real>);               // Converted to degrees for the public API.
                        buffer[indexOfAxisOrientation(idx)] = std::fmod(std::abs(degrees), Real{360}); // Orientation stored in [0, 360].
                    } else {
                        buffer[indexOfAxisOrientation(idx)] = Real{0}; // Degenerate isotropic support: the principal direction is undefined.
                    }
                }
                if (computeInertia) {
                    const Real areaSquared = area * area;
                    Real normMu20 = ::mmcfilters::attributes::numeric::safeDivide(mu20, areaSquared);
                    Real normMu02 = ::mmcfilters::attributes::numeric::safeDivide(mu02, areaSquared);
                    buffer[indexOfInertia(idx)] = normMu20 + normMu02;
                }
                if (computeCircularity) {
                    const Real eps = std::numeric_limits<Real>::epsilon();
                    if (lambda1 <= eps && std::abs(lambda2) <= eps) {
                        buffer[indexOfCircularity(idx)] = Real{1};
                    } else if (lambda1 <= eps || lambda2 <= eps) {
                        buffer[indexOfCircularity(idx)] = Real{0};
                    } else {
                        buffer[indexOfCircularity(idx)] = ::mmcfilters::attributes::numeric::safeDivide(lambda2, lambda1);
                    }
                }
            });
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
