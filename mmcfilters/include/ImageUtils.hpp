#include <utility>
#include <algorithm>


#ifndef IMAGE_UTILS_H
#define IMAGE_UTILS_H

class ImageUtils{
    public:
        
        static int to1D(int x, int y, int numCols) {
            return y * numCols + x;
        }

        static std::pair<int, int> to2D(int index, int numCols) {
            return std::make_pair(index % numCols, index / numCols);
        }


        static int* createRandomColor(int* img, int numColsOfImage, int numRowsOfImage){
            int max = 0;
            int sizeImage = numColsOfImage * numRowsOfImage;
            for (int i = 0; i < sizeImage; i++){
                if (img[i] > max)
                    max = img[i];
            }

            int r[max + 1];
            int g[max + 1];
            int b[max + 1];
            r[0] = 0;
            g[0] = 0;
            r[0] = 0;
            for (int i = 1; i <= max; i++){
                r[i] = rand() % 256;
                g[i] = rand() % 256;
                b[i] = rand() % 256;
            }

            int* output = new int[sizeImage * 3];
            for (int pidx = 0; pidx < (sizeImage * 3); pidx++){
                output[pidx] = 0;
            }

            for (int pidx = 0; pidx < sizeImage; pidx++){
                int cpidx = pidx * 3; // (coloured) for 3 channels
                output[cpidx]     = r[img[pidx]];
                output[cpidx + 1] = g[img[pidx]];
                output[cpidx + 2] = b[img[pidx]];
            }
            return output;
        }

        
};

#endif