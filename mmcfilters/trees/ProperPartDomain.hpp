#pragma once

#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>

namespace mmcfilters {

/**
 * @brief Shape metadata optionally attached to the pixel domain.
 *
 * A connected-subset tree model only requires a finite pixel domain. A regular
 * 2D grid is therefore an optional capability used by image reconstruction and
 * geometric algorithms, not part of the core topology contract.
 */
struct GridDomain2D {
    /// Number of grid rows.
    int rows = 0;
    /// Number of grid columns.
    int columns = 0;

    /**
     * @brief Creates an empty domain.
     */
    constexpr GridDomain2D() noexcept = default;

    /**
     * @brief Creates a domain with the supplied row and column counts.
     *
     * @param numRows Number of rows in the domain.
     * @param numColumns Number of columns in the domain.
     */
    constexpr GridDomain2D(int numRows, int numColumns) noexcept : rows(numRows), columns(numColumns) {}

    /**
     * @brief Returns the number of grid elements after validating the shape.
     *
     * @param context Operation context or diagnostic label.
     * @return The number of grid elements after validating the shape.
     */
    [[nodiscard]] std::size_t size(const char* context) const {
        if (rows <= 0 || columns <= 0) {
            throw std::invalid_argument(std::string(context) + " requires a non-empty 2D grid.");
        }
        const std::size_t gridSize = static_cast<std::size_t>(rows) * static_cast<std::size_t>(columns);
        if (gridSize > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            throw std::invalid_argument(std::string(context) + " exceeds the supported pixel-id range.");
        }
        return gridSize;
    }
};

} // namespace mmcfilters
