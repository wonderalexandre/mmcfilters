#ifndef COMMONS_HPP  
#define COMMONS_HPP  


#define NDEBUG  // Remove os asserts do código
#include <cassert>
#include <cstdint>
#include <list>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <memory>
#include <limits>
#include <algorithm>

#define PRINT_LOG 0 
#define PRINT_DEBUG 0 




template <typename PixelType>
class Image {
    private:
        int numRows;
        int numCols;
        std::shared_ptr<PixelType[]> data;
        using Ptr = std::shared_ptr<Image<PixelType>>;

    public:
    
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

template <typename T>
using ImagePtr = std::shared_ptr<Image<T>>;

#endif 
