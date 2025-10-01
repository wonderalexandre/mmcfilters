#pragma once

#include "../mmcfilters/utils/Common.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <array>

namespace mmcfilters {

namespace py = pybind11;

/**
 * @brief Funções auxiliares para converter entre estruturas C++ e NumPy.
 */
class PybindUtils{
    public:
   
        template <typename PixelType>
        static py::array_t<PixelType> toNumpy(ImagePtr<PixelType> image) {
            int numCols = image->getNumCols();
            int numRows = image->getNumRows();

            std::shared_ptr<PixelType[]> buffer = image->rawDataPtr();
            std::shared_ptr<PixelType[]> bufferCopy = buffer;

            py::capsule free_when_done(new std::shared_ptr<PixelType[]>(bufferCopy), [](void* ptr) {
                delete reinterpret_cast<std::shared_ptr<PixelType[]>*>(ptr);
            });

            // 2D shape: (numRows, numCols), row-major strides
            const py::ssize_t itemsize = sizeof(PixelType);
            const std::array<py::ssize_t, 2> shape   = { static_cast<py::ssize_t>(numRows), static_cast<py::ssize_t>(numCols) };
            const std::array<py::ssize_t, 2> strides = { static_cast<py::ssize_t>(numCols) * itemsize, itemsize };

            py::array_t<PixelType> numpy(
                shape,
                strides,
                buffer.get(),
                free_when_done
            );

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

} // namespace mmcfilters
