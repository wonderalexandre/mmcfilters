
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

#ifndef PYBIND_UTILS_H
#define PYBIND_UTILS_H

class PybindUtils{
    public:
        static py::array_t<int> toNumpy(int* data, int size) {
            // Cria um capsule que sabe como liberar o ponteiro
            py::capsule free_when_done(data, [](void* f) {
                delete[] static_cast<int*>(f);
            });
        
            // Cria o py::array com o capsule responsável por liberar a memória
            return py::array_t<int>(
                { size },               // shape (tamanho do vetor)
                { sizeof(int) },        // strides (distância entre elementos)
                data,                   // ponteiro para os dados
                free_when_done          // capsule que cuida da liberação
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
