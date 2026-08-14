#pragma once

#include "Image.hpp"

#include <memory>

namespace mmcfilters::detail {

/**
 * @brief Image allocation after positive, non-overflowing dimensions were established.
 *
 * Internal reconstruction kernels use this entry point when their dimensions
 * come from an already validated tree grid domain. Shared allocation delegates
 * to the ordinary owning factory so every instantiated pixel type uses the same
 * portable shared-ownership control block.
 */
class CommittedImageAccess {
  public:
    /** @brief Allocates an image value from established dimensions. @param rows Established row count. @param columns Established column count. @return Owned image value. */
    template <typename PixelType> [[nodiscard]] static Image<PixelType> createValue(int rows, int columns) {
        return Image<PixelType>(rows, columns, typename Image<PixelType>::EstablishedDimensionsTag{});
    }

    /** @brief Allocates a shared image from established dimensions. @param rows Established row count. @param columns Established column count. @return Shared owned image. */
    template <typename PixelType> [[nodiscard]] static std::shared_ptr<Image<PixelType>> create(int rows, int columns) {
        return Image<PixelType>::create(rows, columns);
    }
};

} // namespace mmcfilters::detail
