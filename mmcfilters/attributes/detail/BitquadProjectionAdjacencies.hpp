#pragma once

#include "../../trees/MorphologicalTree.hpp"

#include <optional>

namespace mmcfilters::detail {

/**
 * @brief Materializes the two adjacencies used by the current bitquad scalar projection.
 *
 * This is a transition-only algorithm query. Complementary immersions retain
 * their minimum/maximum adjacencies explicitly. The self-dual span immersion
 * preserves the existing 4/4 scalar-projection behavior without representing
 * those derived relations as construction metadata.
 */
inline std::optional<ComplementaryAdjacencies> currentBitquadProjectionAdjacencies(const MorphologicalTree& tree) {
    if (const auto* adjacencies = tree.complementaryAdjacencies()) {
        return *adjacencies;
    }
    const TopographicConvention* convention = tree.topographicConvention();
    if (convention == nullptr || !std::holds_alternative<SelfDualSpanImmersion>(convention->immersion) || !tree.hasGridDomain2D()) {
        return std::nullopt;
    }
    const GridDomain2D& domain = *tree.gridDomain2D();
    return ComplementaryAdjacencies{
        RegularGridAdjacency2D(domain.rows, domain.columns, 1.0), RegularGridAdjacency2D(domain.rows, domain.columns, 1.0)};
}

} // namespace mmcfilters::detail
