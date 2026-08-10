#pragma once

#include "Image.hpp"

#include <memory>

namespace mmcfilters::detail {

/**
 * @brief Image allocation after positive, non-overflowing dimensions were established.
 *
 * Internal reconstruction kernels use this entry point when their dimensions
 * come from an already validated tree grid domain.  Public image construction
 * remains defensive and unchanged.
 */
class CommittedImageAccess {
  public:
    /** @brief Allocates an image value from established dimensions. @param rows Established row count. @param cols Established column count. @return Owned image value. */
    template <typename PixelType> [[nodiscard]] static Image<PixelType> createValue(int rows, int cols) {
        return Image<PixelType>(rows, cols, typename Image<PixelType>::EstablishedDimensionsTag{});
    }

    /** @brief Allocates a shared image from established dimensions. @param rows Established row count. @param cols Established column count. @return Shared owned image. */
    template <typename PixelType> [[nodiscard]] static std::shared_ptr<Image<PixelType>> create(int rows, int cols) {
        return std::shared_ptr<Image<PixelType>>(new Image<PixelType>(rows, cols, typename Image<PixelType>::EstablishedDimensionsTag{}));
    }
};

} // namespace mmcfilters::detail
