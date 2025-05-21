
#include "../include/Common.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

#ifndef PYBIND_UTILS_H
#define PYBIND_UTILS_H

class PybindUtils{
    public:
        
        static py::array_t<uint8_t> toNumpy(ImageUInt8Ptr image) {
            
            std::shared_ptr<uint8_t[]> buffer = image->rawDataPtr();
            int n = image->getSize();
            std::shared_ptr<uint8_t[]> bufferCopy = buffer;

            py::capsule free_when_done(new std::shared_ptr<uint8_t[]>(bufferCopy), [](void* ptr) {
                // Converte de volta e destrói corretamente
                delete reinterpret_cast<std::shared_ptr<uint8_t[]>*>(ptr);
            });
            
            py::array_t<uint8_t> numpy = py::array(py::buffer_info(
                buffer.get(),
                sizeof(uint8_t),
                py::format_descriptor<uint8_t>::value,
                1,
                { n },
                { sizeof(uint8_t) }
            ), free_when_done);
            
            return numpy;

/*
            uint8_t* data = image->rawData();
            int size = image->getNumRows() * image->getNumCols();
            // Cria um capsule que sabe como liberar o ponteiro
            py::capsule free_when_done(data, [](void* f) {
                delete[] static_cast<uint8_t*>(f);
            });
        
            // Cria o py::array com o capsule responsável por liberar a memória
            return py::array_t<uint8_t>(
                { size },               // shape (tamanho do vetor)
                { sizeof(uint8_t) },        // strides (distância entre elementos)
                data,                   // ponteiro para os dados
                free_when_done          // capsule que cuida da liberação
            );
            */
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

        static py::array_t<float> toNumpyShared_ptr(std::shared_ptr<float[]> buffer, int n){
            std::shared_ptr<float[]> bufferCopy = buffer;

            py::capsule free_when_done(new std::shared_ptr<float[]>(bufferCopy), [](void* ptr) {
                // Converte de volta e destrói corretamente
                delete reinterpret_cast<std::shared_ptr<float[]>*>(ptr);
            });
            
            py::array_t<float> numpy = py::array(py::buffer_info(
                buffer.get(),
                sizeof(float),
                py::format_descriptor<float>::value,
                1,
                { n },
                { sizeof(float) }
            ), free_when_done);
            
            return numpy;
        }
        

        static std::shared_ptr<float[]> toShared_ptr(py::array_t<float>& arr) {
            // Cria um capsule que sabe como liberar o ponteiro
            return std::shared_ptr<float[]>(
                static_cast<float*>(arr.request().ptr),
                [obj = py::object(arr)](float*) mutable { obj.dec_ref(); }
            );
        }


};

#endif
