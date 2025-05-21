#include <utility>
#include <algorithm>
#include "../include/Common.hpp"

#ifndef IMAGE_UTILS_H
#define IMAGE_UTILS_H

class ImageUtils{
    public:

    
        // Converte (row, col) para índice 1D (row-major)
        static int to1D(int row, int col, int numCols) {
            return row * numCols + col;
        }

        // Converte índice 1D para (row, col) (row-major)
        static std::pair<int, int> to2D(int index, int numCols) {
            int row = index / numCols;
            int col = index % numCols;
            return std::make_pair(row, col);
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

#endif