#pragma once

#include "AttributeComputerDomain.hpp"
#include "AttributeComputerFamily.hpp"
#include "../detail/AttributeKernelSupport.hpp"
#include "../../contours/ContourTraceComputation.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <limits>
#include <numbers>
#include <span>
#include <string_view>

namespace mmcfilters::attributes::computers {

namespace detail {

/** @brief Discrete moments in coordinates relative to one boundary pixel. */
struct FilledShapeMoments {
    long double area = 0;       ///< Zeroth moment.
    long double column = 0;     ///< First column moment.
    long double row = 0;        ///< First row moment.
    long double column2 = 0;    ///< Second column moment.
    long double row2 = 0;       ///< Second row moment.
    long double columnRow = 0;  ///< Mixed column-row moment.

    /**
     * @brief Adds signed pixel-center moments enclosed by one boundary.
     *
     * Vertical edges delimit horizontal pixel runs. The discrete primitives
     * k, k(k-1)/2, and k(k-1)(2k-1)/6 sum 1, x, and x^2 along each run.
     * East edges add these sums; west edges subtract them. Internal boundaries
     * therefore subtract hole moments. Horizontal edges contribute zero.
     *
     * @param trace Borrowed ordered trace.
     * @param boundary Boundary to integrate.
     * @param columns Image width.
     * @param originColumn Column of the local coordinate origin.
     * @param originRow Row of the local coordinate origin.
     * @param secondOrder Whether second moments are needed.
     */
    void addBoundary(ContourTraceView trace, const ContourBoundary& boundary, int columns,
                     int originColumn, int originRow, bool secondOrder) {
        for (ContourEdge edge : trace.boundaryEdges(boundary)) {
            if (edge.side != ContourSide::East && edge.side != ContourSide::West) {
                continue;
            }
            const bool east = edge.side == ContourSide::East;
            const long double k = static_cast<long double>(edge.pixel % columns - originColumn) + (east ? 1 : 0);
            const long double y = edge.pixel / columns - originRow;
            const long double signedWidth = east ? k : -k;
            const long double columnSum = signedWidth * (k - 1) / 2;
            area += signedWidth;
            column += columnSum;
            row += signedWidth * y;
            if (secondOrder) {
                column2 += columnSum * (2 * k - 1) / 3;
                row2 += signedWidth * y * y;
                columnRow += columnSum * y;
            }
        }
    }
};

} // namespace detail

/**
 * @brief Shape descriptors of the region enclosed by the external boundary.
 *
 * Holes are filled and each enclosed pixel has unit mass at its center.
 * Computation consumes one incremental trace traversal without storing filled
 * masks or retaining traces. Every live node must have exactly one external
 * boundary. Altitude is used only to select connectivity for lower/upper shapes.
 */
class FilledShapeAttributeComputer {
  public:
    static constexpr std::string_view familyName = "filled-shape"; ///< Name used in dependency diagnostics.
    static constexpr AttributeComputerFamily family = AttributeComputerFamily::FilledShape; ///< Scheduler family.
    static constexpr AttributeComputerDomain domain = AttributeComputerDomain::Topology; ///< Support-based computation.

    /// Canonical output order of the family.
    inline static constexpr std::array<Attribute, 12> producedAttributes{
        FilledArea, FilledCentroidRow, FilledCentroidColumn, FilledLengthMajorAxis, FilledLengthMinorAxis,
        FilledAxisOrientation, FilledEccentricity, FilledInertia, HoleAreaFraction, FilledCentroidDisplacementNormalized,
        FilledCompactness, FilledCircularity};

    /** @brief Computes requested filled-region descriptors. @param context Tree, requests, and destination columns. */
    template <std::floating_point Real> static void compute(const AttributeComputeContext<Real>& context) {
        MMCFILTERS_CONTRACT_CHECKED_ONLY(validateContext(context));
        if (!context.requestedAttributes.empty()) {
            computeTraces(context, ContourTraceComputation(context.tree));
        }
    }

    /**
     * @brief Computes descriptors using altitudes to resolve directional connectivity.
     * @param context Valued tree, requests, and destination columns.
     */
    template <std::floating_point Real, AltitudeValue T> static void compute(const AltitudeAttributeComputeContext<Real, T>& context) {
        const AttributeComputeContext<Real> topologyContext{context.tree, context.buffer, context.attrNames, context.requestedAttributes};
        MMCFILTERS_CONTRACT_CHECKED_ONLY(validateContext(topologyContext));
        if (!context.requestedAttributes.empty()) {
            const ValuedMorphologicalTreeView<T> view(context.tree, context.altitude);
            computeTraces(topologyContext, ContourTraceComputation(view));
        }
    }

    /**
     * @brief Writes unit-support values in the supplied pixel order.
     *
     * Area, eccentricity, and circularity are one; the centroid is the pixel coordinate.
     * All remaining descriptors are zero.
     * @param context Pixel order, requests, and destination columns.
     */
    template <std::floating_point Real> static void computeUnitRows(const UnitAttributeComputeContext<Real>& context) {
        requireUnitAttributeBufferShape(context.tree, context.unitPixels, context.buffer, context.attrNames);
        const int columns = context.tree.requireGridDomain2D("FilledShapeAttributeComputer").columns;
        const auto offsets = outputOffsets(context.attrNames, context.requestedAttributes);
        for (NodeId row = 0; row < static_cast<NodeId>(context.unitPixels.size()); ++row) {
            const PixelId pixel = context.unitPixels[static_cast<std::size_t>(row)];
            const std::array<long double, 12> values{1, static_cast<long double>(pixel / columns),
                static_cast<long double>(pixel % columns), 0, 0, 0, 1, 0, 0, 0, 0, 1};
            writeRow(context.buffer, context.attrNames.NUM_ATTRIBUTES, row, offsets, values);
        }
    }

  private:
    /** @brief Validates the borrowed output layout. @param context Tree, requests, and destination columns. */
    template <std::floating_point Real> static void validateContext(const AttributeComputeContext<Real>& context) {
        requireAttributeBufferShape(context.tree, context.buffer, context.attrNames);
        requireRequestedAttributeColumns(context);
    }

    /**
     * @brief Resolves requested column offsets once per computation.
     * @param names Output layout.
     * @param requested Scalars to materialize.
     * @return Column offsets in family order, or -1 for unrequested scalars.
     */
    static std::array<int, 12> outputOffsets(const AttributeNames& names, std::span<const Attribute> requested) {
        std::array<int, 12> offsets;
        offsets.fill(-1);
        for (Attribute attribute : requested) {
            const auto it = std::find(producedAttributes.begin(), producedAttributes.end(), attribute);
            if (it == producedAttributes.end()) {
                throw std::invalid_argument("FilledShapeAttributeComputer received an attribute from another family.");
            }
            offsets[static_cast<std::size_t>(it - producedAttributes.begin())] = names.linearIndex(0, attribute);
        }
        return offsets;
    }

    /**
     * @brief Copies requested values into one output row.
     * @param buffer Destination buffer.
     * @param stride Columns per row.
     * @param node Output row identifier.
     * @param offsets Requested columns in family order.
     * @param values Computed descriptors in family order.
     */
    template <std::floating_point Real>
    static void writeRow(std::span<Real> buffer, int stride, NodeId node, const std::array<int, 12>& offsets,
                         const std::array<long double, 12>& values) {
        const std::size_t row = static_cast<std::size_t>(node) * static_cast<std::size_t>(stride);
        for (std::size_t i = 0; i < offsets.size(); ++i) {
            if (offsets[i] >= 0) {
                buffer[row + static_cast<std::size_t>(offsets[i])] = static_cast<Real>(values[i]);
            }
        }
    }

    /**
     * @brief Integrates boundaries and materializes each live node immediately.
     * @param context Tree, requests, and destination columns.
     * @param traces Incremental trace computation with resolved connectivity.
     */
    template <std::floating_point Real>
    static void computeTraces(const AttributeComputeContext<Real>& context, const ContourTraceComputation& traces) {
        const auto offsets = outputOffsets(context.attrNames, context.requestedAttributes);
        const bool secondOrder = offsets[10] >= 0 || offsets[11] >= 0 ||
            std::any_of(offsets.begin() + 3, offsets.begin() + 8, [](int offset) { return offset >= 0; });
        const bool compareSupport = offsets[8] >= 0 || offsets[9] >= 0;
        const int columns = context.tree.numColumns();
        traces.forEachTrace([&](NodeId node, ContourTraceView trace) {
            const ContourBoundary external = trace.externalBoundary();
            const PixelId origin = (*trace.boundaryEdges(external).begin()).pixel;
            const int originColumn = origin % columns;
            const int originRow = origin / columns;
            detail::FilledShapeMoments filled;
            filled.addBoundary(trace, external, columns, originColumn, originRow, secondOrder);
            const long double area = filled.area;
            const long double centroidColumn = filled.column / area;
            const long double centroidRow = filled.row / area;
            std::array<long double, 12> values{area, centroidRow + originRow, centroidColumn + originColumn, 0, 0, 0, 1, 0, 0, 0, 0, 1};
            if (secondOrder) {
                const long double mu20 = std::max(0.L, filled.column2 - filled.column * centroidColumn);
                const long double mu02 = std::max(0.L, filled.row2 - filled.row * centroidRow);
                const long double mu11 = filled.columnRow - filled.column * centroidRow;
                const long double discriminant = std::hypot(mu20 - mu02, 2 * mu11);
                const long double lambda1 = mu20 + mu02 + discriminant;
                const long double lambda2 = std::max(0.L, mu20 + mu02 - discriminant);
                values[3] = std::sqrt(2 * lambda1 / area);
                values[4] = std::sqrt(2 * lambda2 / area);
                if (mu20 != mu02 || mu11 != 0) {
                    values[5] = std::abs(std::atan2(2 * mu11, mu20 - mu02) * 90 / std::numbers::pi_v<long double>);
                }
                const long double epsilon = std::numeric_limits<Real>::epsilon();
                values[6] = lambda1 <= epsilon ? 1 : (lambda2 <= epsilon ? 1e6L : std::min(lambda1 / lambda2, 1e6L));
                values[7] = (mu20 + mu02) / (area * area);
                values[10] = mu20 + mu02 > epsilon ? area / (2 * std::numbers::pi_v<long double> * (mu20 + mu02)) : 0;
                values[11] = lambda1 <= epsilon ? 1 : (lambda2 <= epsilon ? 0 : lambda2 / lambda1);
            }
            if (compareSupport) {
                auto support = filled;
                for (const ContourBoundary& boundary : trace.boundaries()) {
                    if (boundary.kind == ContourBoundaryKind::Internal) {
                        support.addBoundary(trace, boundary, columns, originColumn, originRow, false);
                    }
                }
                values[8] = std::clamp(1 - support.area / area, 0.L, 1.L);
                values[9] = std::hypot(centroidColumn - support.column / support.area,
                                       centroidRow - support.row / support.area) / std::sqrt(area);
            }
            writeRow(context.buffer, context.attrNames.NUM_ATTRIBUTES, node, offsets, values);
        });
    }
};

} // namespace mmcfilters::attributes::computers
