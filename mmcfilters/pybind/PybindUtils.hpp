
#include "../include/Common.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

#ifndef PYBIND_UTILS_H
#define PYBIND_UTILS_H

class PybindUtils{
    public:
        static py::array_t<PixelType> toNumpy(PixelType* data, int size) {
            // Cria um capsule que sabe como liberar o ponteiro
            py::capsule free_when_done(data, [](void* f) {
                delete[] static_cast<PixelType*>(f);
            });
        
            // Cria o py::array com o capsule responsável por liberar a memória
            return py::array_t<PixelType>(
                { size },               // shape (tamanho do vetor)
                { sizeof(PixelType) },        // strides (distância entre elementos)
                data,                   // ponteiro para os dados
                free_when_done          // capsule que cuida da liberação
            );
        }
    

        static py::array_t<int> toNumpyInt(int* data, int size) {
            // Cria capsule com função de destruição
            py::capsule free_when_done(data, [](void* f) {
                delete[] static_cast<int*>(f);
            });
        
            // Cria o array NumPy com os dados e o capsule
            return py::array_t<int>(
                { size },                // shape (1D)
                { sizeof(int) },       // strides
                data,                    // ponteiro para os dados
                free_when_done           // capsule que cuida da liberação
            );
        }

        static py::array_t<float> toNumpyFloat(float* data, int size) {
            // Cria capsule com função de destruição
            py::capsule free_when_done(data, [](void* f) {
                delete[] static_cast<float*>(f);
            });
        
            // Cria o array NumPy com os dados e o capsule
            return py::array_t<float>(
                { size },                // shape (1D)
                { sizeof(float) },       // strides
                data,                    // ponteiro para os dados
                free_when_done           // capsule que cuida da liberação
            );
        }
        




};

#endif
