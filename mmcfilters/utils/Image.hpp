#pragma once

/**
 * @file Image.hpp
 * @brief Contiguous row-major image container and coordinate helpers.
 */

#include <memory>
#include <cstdint>
#include <cstddef>
#include <algorithm>
#include <utility>
#include <cstdlib>
#include <limits>
#include <stdexcept>

namespace mmcfilters {

/**
 * @brief Generic row-major 2D image with contiguous storage and shared ownership.
 *
 * `Image<PixelType>` stores a rectangular image in a single contiguous
 * row-major buffer. The class exposes simple factory helpers for owned and
 * non-owned memory, deep copy support, linear access, and a few convenience
 * operations used throughout the library.
 *
 * Ownership semantics:
 *
 * - `create(rows, cols)` allocates and owns a fresh buffer;
 * - `create(rows, cols, initValue)` allocates, owns, and fills the buffer;
 * - `fromExternal(rawPtr, rows, cols)` wraps external memory without taking ownership;
 * - `fromRaw(rawPtr, rows, cols)` wraps external memory and assumes ownership.
 *
 * @tparam PixelType Pixel scalar type, such as `uint8_t`, `int32_t`, or `float`.
 */
template <typename PixelType> class Image {
  private:
    /** @brief Stores the num rows. */
    int numRows;
    /** @brief Stores the num cols. */
    int numCols;
    /** @brief Stores the data. */
    std::shared_ptr<PixelType[]> data;
    /** @brief Defines the `Ptr` alias used by the component. */
    using Ptr = std::shared_ptr<Image<PixelType>>;

    /**
     * @brief Checks and converts size.
     *
     * @param rows Number of rows in the domain.
     * @param cols Number of columns in the domain.
     * @return Validated number of pixels in the image.
     */
    static std::size_t checkedSize(int rows, int cols) {
        if (rows <= 0 || cols <= 0) {
            throw std::invalid_argument("Image dimensions must be positive.");
        }
        const auto rowCount = static_cast<std::size_t>(rows);
        const auto colCount = static_cast<std::size_t>(cols);
        if (rowCount > static_cast<std::size_t>(std::numeric_limits<int>::max()) / colCount) {
            throw std::overflow_error("Image dimensions exceed the supported int-indexed size.");
        }
        return rowCount * colCount;
    }

  public:
    /// Pixel scalar type stored by this image.
    using Type = PixelType;

    /**
     * @brief Allocates an owned row-major image buffer.
     *
     * @param rows Number of image rows.
     * @param cols Number of image columns.
     * @throws std::invalid_argument If either dimension is not positive.
     * @throws std::overflow_error If `rows * cols` exceeds the supported size.
     */
    Image(int rows, int cols) : numRows(rows), numCols(cols), data(new PixelType[checkedSize(rows, cols)], std::default_delete<PixelType[]>()) {}

    /**
     * @brief Creates an owned image with uninitialised pixel values.
     *
     * @param rows Number of image rows.
     * @param cols Number of image columns.
     * @return Shared owning image wrapper.
     * @throws std::invalid_argument If either dimension is not positive.
     * @throws std::overflow_error If `rows * cols` exceeds the supported size.
     */
    [[nodiscard]] static Ptr create(int rows, int cols) { return std::make_shared<Image>(rows, cols); }

    /**
     * @brief Creates an owned image and fills every pixel with `initValue`.
     *
     * @param rows Number of image rows.
     * @param cols Number of image columns.
     * @param initValue Value assigned to every pixel.
     * @return Shared owning image wrapper.
     * @throws std::invalid_argument If either dimension is not positive.
     * @throws std::overflow_error If `rows * cols` exceeds the supported size.
     */
    [[nodiscard]] static Ptr create(int rows, int cols, PixelType initValue) {
        auto img = create(rows, cols);
        img->fill(initValue);
        return img;
    }

    /**
     * @brief Wraps external row-major memory without taking ownership.
     *
     * The caller remains responsible for keeping `rawPtr` alive while the image
     * wrapper is used.
     *
     * @param rawPtr Pointer to `rows * cols` row-major pixels.
     * @param rows Number of image rows.
     * @param cols Number of image columns.
     * @return Shared non-owning image wrapper.
     * @throws std::invalid_argument If `rawPtr` is null or dimensions are invalid.
     * @throws std::overflow_error If `rows * cols` exceeds the supported size.
     */
    [[nodiscard]] static Ptr fromExternal(PixelType* rawPtr, int rows, int cols) {
        if (rawPtr == nullptr) {
            throw std::invalid_argument("Image::fromExternal requires a non-null raw pointer.");
        }
        auto img = create(rows, cols);
        img->data = std::shared_ptr<PixelType[]>(rawPtr, [](PixelType*) {
            // Empty deleter: the wrapper does not release external memory.
        });
        return img;
    }

    /**
     * @brief Wraps row-major memory and takes ownership through `delete[]`.
     *
     * @param rawPtr Pointer allocated with `new PixelType[]`.
     * @param rows Number of image rows.
     * @param cols Number of image columns.
     * @return Shared owning image wrapper.
     * @throws std::invalid_argument If `rawPtr` is null or dimensions are invalid.
     * @throws std::overflow_error If `rows * cols` exceeds the supported size.
     */
    [[nodiscard]] static Ptr fromRaw(PixelType* rawPtr, int rows, int cols) {
        if (rawPtr == nullptr) {
            throw std::invalid_argument("Image::fromRaw requires a non-null raw pointer.");
        }
        auto img = create(rows, cols);
        img->data = std::shared_ptr<PixelType[]>(rawPtr, std::default_delete<PixelType[]>());
        return img;
    }

    /**
     * @brief Fills every pixel with `value`.
     *
     * Complexity: O(rows * cols).
     *
     * @param value Value used by the operation.
     */
    void fill(PixelType value) { std::fill_n(data.get(), numRows * numCols, value); }

    /**
     * @brief Returns true when shape and pixel values match `other`.
     *
     * @param other Non-null image pointer to compare against.
     * @return True only when dimensions and every row-major pixel match.
     *
     * Complexity: O(rows * cols).
     */
    bool isEqual(const Ptr& other) const {
        if (numRows != other->numRows || numCols != other->numCols)
            return false;
        int n = numRows * numCols;
        for (int i = 0; i < n; ++i) {
            if (data[i] != (*other)[i])
                return false;
        }
        return true;
    }

    /**
     * @brief Returns a deep copy with its own owned buffer.
     *
     * @return New image with the same dimensions and pixel values.
     *
     * Complexity: O(rows * cols).
     */
    [[nodiscard]] Ptr clone() const {
        auto newImg = create(numRows, numCols);
        std::copy(data.get(), data.get() + (numRows * numCols), newImg->data.get());
        return newImg;
    }

    /**
     * @brief Returns the shared pointer that owns or wraps the raw buffer.
     *
     * The returned pointer shares the same ownership/lifetime policy as this
     * image: either owning storage, externally wrapped storage, or adopted raw
     * storage.
     *
     * @return The shared pointer that owns or wraps the raw buffer.
     */
    std::shared_ptr<PixelType[]> rawDataPtr() { return data; }

    /**
     * @brief Returns a mutable pointer to the contiguous row-major buffer.
     *
     * The pointer remains valid while the image object and its shared buffer are
     * alive and no external owner invalidates externally wrapped memory.
     *
     * @return A mutable pointer to the contiguous row-major buffer.
     */
    PixelType* rawData() { return data.get(); }

    /**
     * @brief Returns the number of image rows.
     *
     * @return The number of image rows.
     */
    int getNumRows() const { return numRows; }

    /**
     * @brief Returns the number of image columns.
     *
     * @return The number of image columns.
     */
    int getNumCols() const { return numCols; }

    /**
     * @brief Returns the total number of pixels.
     *
     * @return The total number of pixels.
     */
    int getSize() const { return numRows * numCols; }

    /**
     * @brief Returns mutable linear access to pixel `index`.
     *
     * `index` is interpreted in row-major order as `row * getNumCols() + col`.
     * The operator does not perform bounds checking.
     *
     * @param index Zero-based index used by the operation.
     * @return Mutable linear access to pixel index.
     */
    PixelType& operator[](int index) { return data[index]; }

    /**
     * @brief Returns immutable linear access to pixel `index`.
     *
     * `index` is interpreted in row-major order as `row * getNumCols() + col`.
     * The operator does not perform bounds checking.
     *
     * @param index Zero-based index used by the operation.
     * @return Immutable linear access to pixel index.
     */
    const PixelType& operator[](int index) const { return data[index]; }
};

/// 8-bit unsigned image container.
using ImageUInt8 = Image<uint8_t>;
/// 32-bit signed integer image container.
using ImageInt32 = Image<int32_t>;
/// Single-precision floating-point image container.
using ImageFloat = Image<float>;

/// Shared pointer to an 8-bit unsigned image.
using ImageUInt8Ptr = std::shared_ptr<ImageUInt8>;
/// Shared pointer to a 32-bit signed integer image.
using ImageInt32Ptr = std::shared_ptr<ImageInt32>;
/// Shared pointer to a single-precision floating-point image.
using ImageFloatPtr = std::shared_ptr<ImageFloat>;

/**
 * @brief Shared pointer alias for an image with arbitrary pixel type.
 */
template <typename PixelType> using ImagePtr = std::shared_ptr<Image<PixelType>>;

/**
 * @brief Utility functions for basic image conversion and manipulation.
 *
 * Groups helpers for converting between 1D/2D coordinates and for generating
 * colour visualisations from integer labels.
 */
class ImageUtils {
  public:
    /**
     * @brief Converts `(row, col)` to a row-major linear index.
     *
     * @param row Zero-based row coordinate.
     * @param col Zero-based column coordinate.
     * @param numCols Number of columns in the image domain.
     * @return Linear index `row * numCols + col`.
     */
    inline static int to1D(int row, int col, int numCols) noexcept { return row * numCols + col; }

    /**
     * @brief Converts a row-major linear index to `(row, col)`.
     *
     * @param index Row-major linear index.
     * @param numCols Number of columns in the image domain.
     * @return Pair `(row, col)`.
     */
    inline static std::pair<int, int> to2D(int index, int numCols) noexcept {
        int row = index / numCols;
        int col = index - row * numCols;
        return {row, col};
    }

    /**
     * @brief Creates a random-colour visualisation from an integer-labelled image.
     *
     * The output is a row-major `uint8_t` image with three adjacent columns per
     * input pixel, storing RGB triples as `(R, G, B)`. Labels are used as array
     * indices and therefore must be non-negative.
     *
     * @param img Pointer to `numRowsOfImage * numColsOfImage` integer labels.
     * @param numRowsOfImage Number of input rows.
     * @param numColsOfImage Number of input columns.
     * @return RGB visualisation encoded as an image with `numColsOfImage * 3`
     * columns.
     */
    [[nodiscard]] static ImageUInt8Ptr createRandomColor(int* img, int numRowsOfImage, int numColsOfImage) {
        int max = 0;
        int sizeImage = numColsOfImage * numRowsOfImage;
        for (int i = 0; i < sizeImage; i++) {
            if (img[i] > max)
                max = img[i];
        }

        std::unique_ptr<int[]> r(new int[max + 1]);
        std::unique_ptr<int[]> g(new int[max + 1]);
        std::unique_ptr<int[]> b(new int[max + 1]);
        r[0] = 0;
        g[0] = 0;
        r[0] = 0;
        for (int i = 1; i <= max; i++) {
            r[i] = rand() % 256;
            g[i] = rand() % 256;
            b[i] = rand() % 256;
        }

        int sizeOutput = sizeImage * 3; // [(R,G,B), (R,G,B), ...]
        ImageUInt8Ptr outImage = ImageUInt8::create(numRowsOfImage, numColsOfImage * 3);

        auto output = outImage->rawData();
        // Initialise with zero.
        std::fill_n(output, sizeOutput, 0);

        for (int pidx = 0; pidx < sizeImage; pidx++) {
            int cpidx = pidx * 3; // (coloured) for 3 channels
            output[cpidx] = r[img[pidx]];
            output[cpidx + 1] = g[img[pidx]];
            output[cpidx + 2] = b[img[pidx]];
        }
        return outImage;
    }
};

} // namespace mmcfilters
