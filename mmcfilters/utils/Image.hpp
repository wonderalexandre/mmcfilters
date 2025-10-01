#pragma once

#include <memory> // Ponteiros inteligentes
#include <cstdint> // Tipos fixos
#include <algorithm> // Algoritmos STL
#include <utility> // Utilidades diversas
#include <cstdlib>  // para rand()

namespace mmcfilters {

/**
 * @brief Classe de imagem genérica 2D com armazenamento contíguo e controle de vida via std::shared_ptr.
 *
 * A `Image<PixelType>` representa uma imagem 2D em ordem row-major, encapsulando
 * largura, altura e um buffer contíguo gerenciado por `std::shared_ptr<PixelType[]>`.
 * Fornece utilitários para criação/copiar/preencher e acesso indexado em 1D.
 *
 * ## Semântica de propriedade do buffer
 * - `create(rows, cols)`: aloca um novo buffer e o gerencia (deleter padrão).
 * - `create(rows, cols, initValue)`: idem, preenchendo com o valor inicial.
 * - `fromExternal(rawPtr, rows, cols)`: **não** assume a propriedade
 *   (o deleter é vazio). Útil quando o ciclo de vida do ponteiro bruto é externo.
 * - `fromRaw(rawPtr, rows, cols)`: **assume** a propriedade do array
 *   (deleter padrão de array). Use quando a instância deve gerenciar a memória.
 *
 * ## Layout de memória
 * - Acesso linear (row-major): índice `i = row * numCols + col`.
 * - Operador `operator[](int)` fornece acesso por índice linear.
 *
 *
 * ## Exemplo de uso
 * @code
 * using ImgU8 = Image<uint8_t>;
 * auto img = ImgU8::create(480, 640, 0);     // aloca e zera
 * img->fill(255);                            // preenche com 255
 * int idx = ImageUtils::to1D(10, 20, img->getNumCols());
 * (*img)[idx] = 128;                         // acesso linear
 * auto clone = img->clone();                  // deep copy
 * bool eq = img->isEqual(clone);              // true
 * @endcode
 *
 * @tparam PixelType tipo do pixel armazenado (ex.: uint8_t, int32_t, float).
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
    // Converte (row, col) para índice 1D (row-major)
    inline static int to1D(int row, int col, int numCols) noexcept{
        return row * numCols + col;
    }

    // Converte índice 1D para (row, col) (row-major)
    inline static std::pair<int, int> to2D(int index, int numCols) noexcept {
        int row = index / numCols;
        int col = index - row * numCols;  // evita operador % => int col = index % numCols;
        return {row, col};
    }

    // Cria uma imagem colorida aleatória a partir de uma imagem em escala de cinza: [(R,G,B), (R,G,B), ...]
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