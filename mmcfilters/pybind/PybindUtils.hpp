
#include "../include/Common.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

#ifndef PYBIND_UTILS_H
#define PYBIND_UTILS_H

class PybindUtils{
    public:
        
        template <typename PixelType>
        static py::array_t<PixelType> toNumpy(ImagePtr<PixelType> image) {
            
            std::shared_ptr<PixelType[]> buffer = image->rawDataPtr();
            int n = image->getSize();
            std::shared_ptr<PixelType[]> bufferCopy = buffer;

            py::capsule free_when_done(new std::shared_ptr<PixelType[]>(bufferCopy), [](void* ptr) {
                // Converte de volta e destrói corretamente
                delete reinterpret_cast<std::shared_ptr<PixelType[]>*>(ptr);
            });
            
            py::array_t<PixelType> numpy = py::array(py::buffer_info(
                buffer.get(),
                sizeof(PixelType),
                py::format_descriptor<PixelType>::value,
                1,
                { n },
                { sizeof(PixelType) }
            ), free_when_done);
            
            return numpy;
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
        
        /*
        static std::shared_ptr<float[]> toShared_ptr(py::array_t<float>& arr) {
            // Cria um capsule que sabe como liberar o ponteiro
            return std::shared_ptr<float[]>(
                static_cast<float*>(arr.request().ptr),
                [obj = py::object(arr)](float*) mutable { obj.dec_ref(); }
            );
        }*/

        static std::shared_ptr<float[]> toShared_ptr(py::array_t<float>& arr) {
            // Captura o objeto Python no deleter — isso garante que o buffer não será liberado prematuramente
            return std::shared_ptr<float[]>(
                static_cast<float*>(arr.request().ptr),
                [obj = py::object(arr)](float*) mutable {
                    //  manter o py::object vivo
                }
            );
        }


};

#endif
