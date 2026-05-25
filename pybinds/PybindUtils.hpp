#pragma once

#include "../mmcfilters/utils/Image.hpp"
#include "../mmcfilters/utils/Common.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <array>
#include <sstream>
#include <string_view>
#include <vector>

namespace mmcfilters {

namespace py = pybind11;

/**
 * @brief Helper functions for converting between C++ buffers/images and NumPy arrays.
 */
class PybindUtils{
    public:

        static void require1DArray(const py::buffer_info& bufferInfo, py::ssize_t expectedSize, std::string_view argumentName) {
            if (bufferInfo.ndim != 1) {
                std::ostringstream oss;
                oss << argumentName << " must be a 1D array";
                throw std::invalid_argument(oss.str());
            }
            if (bufferInfo.shape[0] != expectedSize) {
                std::ostringstream oss;
                oss << argumentName << " must have length " << expectedSize
                    << ", got " << bufferInfo.shape[0];
                throw std::invalid_argument(oss.str());
            }
        }

        template <class T>
        static void requireVectorSize(const std::vector<T>& values, std::size_t expectedSize, std::string_view argumentName) {
            if (values.size() != expectedSize) {
                std::ostringstream oss;
                oss << argumentName << " must have length " << expectedSize
                    << ", got " << values.size();
                throw std::invalid_argument(oss.str());
            }
        }
   
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
	        // Create a capsule that owns destruction of the raw buffer.
	        py::capsule free_when_done(data, [](void* f) {
	            delete[] static_cast<int*>(f);
	        });

	        // Create a NumPy view backed by the capsule-owned data.
	        return py::array_t<int>(
	            { size },                // shape (1D)
	            { sizeof(int) },       // strides
	            data,                    // data pointer
	            free_when_done           // capsule that releases the buffer
	        );
	    }

	    static py::array_t<float> toNumpyFloat(float* data, int size) {
	        // Create a capsule that owns destruction of the raw buffer.
	        py::capsule free_when_done(data, [](void* f) {
	            delete[] static_cast<float*>(f);
	        });

	        // Create a NumPy view backed by the capsule-owned data.
	        return py::array_t<float>(
	            { size },                // shape (1D)
	            { sizeof(float) },       // strides
	            data,                    // data pointer
	            free_when_done           // capsule that releases the buffer
	        );
	    }

        static py::array_t<float> toNumpyShared_ptr(std::shared_ptr<float[]> buffer, int n){
	        std::shared_ptr<float[]> bufferCopy = buffer;

	        py::capsule free_when_done(new std::shared_ptr<float[]>(bufferCopy), [](void* ptr) {
	            // Convert back and destroy the shared_ptr holder.
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

        static py::array_t<float> toNumpyOwned(std::vector<float>&& buffer, int n) {
            auto* owned = new std::vector<float>(std::move(buffer));
            py::capsule free_when_done(owned, [](void* ptr) {
                delete reinterpret_cast<std::vector<float>*>(ptr);
            });

            return py::array_t<float>(
                {n},
                {static_cast<py::ssize_t>(sizeof(float))},
                owned->data(),
                free_when_done
            );
        }

        static py::array_t<float> toNumpyOwned2D(std::vector<float>&& buffer, int rows, int cols) {
            auto* owned = new std::vector<float>(std::move(buffer));
            py::capsule free_when_done(owned, [](void* ptr) {
                delete reinterpret_cast<std::vector<float>*>(ptr);
            });

            return py::array_t<float>(
                {rows, cols},
                {static_cast<py::ssize_t>(sizeof(float) * cols), static_cast<py::ssize_t>(sizeof(float))},
                owned->data(),
                free_when_done
            );
        }
        
	    template <typename Array>
	    static std::shared_ptr<float[]> toShared_ptr(const Array& arr) {
	        // Capture the Python object in the deleter so the buffer stays alive.
	        return std::shared_ptr<float[]>(
	            static_cast<float*>(arr.request().ptr),
	            [obj = py::object(arr)](float*) mutable {
	                // Keep the py::object alive until the shared_ptr is destroyed.
	            }
	        );
	    }


};

} // namespace mmcfilters
