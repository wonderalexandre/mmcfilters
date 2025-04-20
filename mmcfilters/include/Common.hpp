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

#define PRINT_LOG 1 
#define PRINT_DEBUG 0 

// Forward declaration 
class MorphologicalTree;
class NodeMT;


using NodeMTPtr = std::shared_ptr<NodeMT>;
using MorphologicalTreePtr = std::shared_ptr<MorphologicalTree>; 

using PixelValueType = uint8_t;

struct Image {
    int numRows;
    int numCols;
    std::shared_ptr<PixelValueType[]> data;


    Image(int rows, int cols) : numRows(rows), numCols(cols), data(new PixelValueType[rows * cols], std::default_delete<PixelValueType[]>()) {}

    static std::shared_ptr<Image> create(int rows, int cols) {
        return std::make_shared<Image>(rows, cols);
    }

    static std::shared_ptr<Image> create(int rows, int cols, PixelValueType initValue) {
        auto img = std::make_shared<Image>(rows, cols);
        img->fill(initValue);
        return img;
    }

    static std::shared_ptr<Image> fromExternal(PixelValueType* rawPtr, int rows, int cols) {
        auto img = std::make_shared<Image>(rows, cols);
        img->data = std::shared_ptr<PixelValueType[]>(rawPtr, [](PixelValueType*) {
            // deleter vazio: não libera o ponteiro
        });
        return img;
    }

    static std::shared_ptr<Image> fromRaw(PixelValueType* rawPtr, int rows, int cols) {
        auto img = std::make_shared<Image>(rows, cols);
        img->data = std::shared_ptr<PixelValueType[]>(rawPtr, std::default_delete<PixelValueType[]>());
        return img;
    }

    PixelValueType* rawData() { return data.get(); }

    PixelValueType& operator[](int index) { return data[index]; }

    const PixelValueType& operator[](int index) const { return data[index]; }

    void fill(PixelValueType value) {
        std::fill(data.get(), data.get() + numRows * numCols, value);
    }

    bool isEqual(std::shared_ptr<Image> other) const {
        if (numRows != other->numRows || numCols != other->numCols)
            return false;
    
        int n = numRows * numCols;
        for (int i = 0; i < n; ++i) {
            if (data[i] != (*other)[i])
                return false;
        }
    
        return true;
    }
};

using ImagePtr = std::shared_ptr<Image>;

#endif 
