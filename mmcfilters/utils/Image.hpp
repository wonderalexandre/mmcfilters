#pragma once

#include <memory>
#include <cstdint>
#include <algorithm>
#include <utility>
#include <cstdlib>

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
template <typename PixelType>
class Image {
    private:
        int numRows;
        int numCols;
        std::shared_ptr<PixelType[]> data;
        using Ptr = std::shared_ptr<Image<PixelType>>;
        
    public:
    using Type = PixelType;

    Image(int rows, int cols): numRows(rows), numCols(cols), data(new PixelType[rows * cols], std::default_delete<PixelType[]>()) {}

    static Ptr create(int rows, int cols) {
        return std::make_shared<Image>(rows, cols);
    }

    static Ptr create(int rows, int cols, PixelType initValue) {
        auto img = create(rows, cols);
        img->fill(initValue);
        return img;
    }

    static Ptr fromExternal(PixelType* rawPtr, int rows, int cols) {
        auto img = create(rows, cols);
        img->data = std::shared_ptr<PixelType[]>(rawPtr, [](PixelType*) {
            // deleter vazio: não libera o ponteiro
        });
        return img;
    }

    static Ptr fromRaw(PixelType* rawPtr, int rows, int cols) {
        auto img = create(rows, cols);
        img->data = std::shared_ptr<PixelType[]>(rawPtr, std::default_delete<PixelType[]>());
        return img;
    }

    
    void fill(PixelType value) {
        std::fill_n(data.get(), numRows * numCols, value);
    }

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
    
    Ptr clone() const {
        auto newImg = create(numRows, numCols);
        std::copy(data.get(), data.get() + (numRows * numCols), newImg->data.get());
        return newImg;
    }

    std::shared_ptr<PixelType[]> rawDataPtr(){ return data; }
    PixelType* rawData() { return data.get(); }
    int getNumRows() const { return numRows; }
    int getNumCols() const { return numCols; }
    int getSize() const { return numRows * numCols; }
    PixelType& operator[](int index) { return data[index]; }
    const PixelType& operator[](int index) const { return data[index]; }


};

// Aliases
using ImageUInt8 = Image<uint8_t>;
using ImageInt32 = Image<int32_t>;
using ImageFloat = Image<float>;

using ImageUInt8Ptr = std::shared_ptr<ImageUInt8>;
using ImageInt32Ptr = std::shared_ptr<ImageInt32>;
using ImageFloatPtr = std::shared_ptr<ImageFloat>;

template <typename PixelType>
using ImagePtr = std::shared_ptr<Image<PixelType>>;


/**
 * @brief Funções utilitárias para conversões e manipulação básica de imagens.
 *
 * Agrupa helpers relacionados à transformação entre coordenadas 1D/2D e à
 * geração de imagens coloridas a partir de rótulos inteiros, facilitando
 * depuração e visualização de resultados intermediários.
 */
class ImageUtils{
public:
    /**
     * @brief Converts `(row, col)` to a row-major linear index.
     */
    inline static int to1D(int row, int col, int numCols) noexcept{
        return row * numCols + col;
    }

    /**
     * @brief Converts a row-major linear index to `(row, col)`.
     */
    inline static std::pair<int, int> to2D(int index, int numCols) noexcept {
        int row = index / numCols;
        int col = index - row * numCols;
        return {row, col};
    }

    /**
     * @brief Creates a random-colour visualisation from an integer-labelled image.
     */
    static ImageUInt8Ptr createRandomColor(int* img, int numRowsOfImage, int numColsOfImage){
        int max = 0;
        int sizeImage = numColsOfImage * numRowsOfImage;
        for (int i = 0; i < sizeImage; i++){
            if (img[i] > max)
                max = img[i];
        }

        std::unique_ptr<int[]> r(new int[max + 1]);
        std::unique_ptr<int[]> g(new int[max + 1]);
        std::unique_ptr<int[]> b(new int[max + 1]);
        r[0] = 0;
        g[0] = 0;
        r[0] = 0;
        for (int i = 1; i <= max; i++){
            r[i] = rand() % 256;
            g[i] = rand() % 256;
            b[i] = rand() % 256;
        }
        
        int sizeOutput = sizeImage * 3; // [(R,G,B), (R,G,B), ...]
        ImageUInt8Ptr outImage = ImageUInt8::create(numRowsOfImage, numColsOfImage * 3);
        
        auto output = outImage->rawData();
            // Inicializa com zero
        std::fill_n(output, sizeOutput, 0);

        for (int pidx = 0; pidx < sizeImage; pidx++){
            int cpidx = pidx * 3; // (coloured) for 3 channels
            output[cpidx]     = r[img[pidx]];
            output[cpidx + 1] = g[img[pidx]];
            output[cpidx + 2] = b[img[pidx]];
        }
        return outImage;
    }


};

}
