#pragma once

#include "../mmcfilters/trees/MorphologicalTree.hpp"
#include "../mmcfilters/utils/Image.hpp"
#include "../mmcfilters/utils/Common.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <array>
#include <concepts>
#include <stdexcept>
#include <string>
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
        enum class FloatingDType {
            Float32,
            Float64,
        };

        static FloatingDType parseFloatingDType(py::object dtype, std::string_view argumentName = "dtype") {
            if (dtype.is_none()) {
                return FloatingDType::Float32;
            }

            py::object numpy = py::module_::import("numpy");
            py::object normalized = numpy.attr("dtype")(dtype);
            const std::string name = py::str(normalized.attr("name")).cast<std::string>();
            if (name == "float32") {
                return FloatingDType::Float32;
            }
            if (name == "float64") {
                return FloatingDType::Float64;
            }
            throw std::invalid_argument(std::string(argumentName) + " must be np.float32 or np.float64");
        }

        static FloatingDType parseFloatingArrayDType(const py::array& array, std::string_view argumentName) {
            py::object numpy = py::module_::import("numpy");
            py::object normalized = numpy.attr("dtype")(array.dtype());
            const std::string name = py::str(normalized.attr("name")).cast<std::string>();
            if (name == "float32") {
                return FloatingDType::Float32;
            }
            if (name == "float64") {
                return FloatingDType::Float64;
            }
            throw std::invalid_argument(std::string(argumentName) + " must be a 1D np.float32 or np.float64 array.");
        }

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

        template <std::floating_point Real>
        static py::array_t<Real, py::array::c_style> requireNodeAttributeArray(py::array attr, const MorphologicalTree& tree, std::string_view argumentName = "attr") {
            const py::buffer_info info = attr.request();
            require1DArray(info, tree.getNumInternalNodeSlots(), argumentName);
            if (info.strides[0] != static_cast<py::ssize_t>(sizeof(Real))) {
                throw std::invalid_argument(std::string(argumentName) + " must be C-contiguous.");
            }
            return py::reinterpret_borrow<py::array_t<Real, py::array::c_style>>(attr);
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
    

        template <typename Real>
        static py::array_t<Real> toNumpyOwned(std::vector<Real>&& buffer, int n) {
            auto* owned = new std::vector<Real>(std::move(buffer));
            py::capsule free_when_done(owned, [](void* ptr) {
                delete reinterpret_cast<std::vector<Real>*>(ptr);
            });

            return py::array_t<Real>(
                {n},
                {static_cast<py::ssize_t>(sizeof(Real))},
                owned->data(),
                free_when_done
            );
        }

        template <typename Real>
        static py::array_t<Real> toNumpyOwned2D(std::vector<Real>&& buffer, int rows, int cols) {
            auto* owned = new std::vector<Real>(std::move(buffer));
            py::capsule free_when_done(owned, [](void* ptr) {
                delete reinterpret_cast<std::vector<Real>*>(ptr);
            });

            return py::array_t<Real>(
                {rows, cols},
                {static_cast<py::ssize_t>(sizeof(Real) * cols), static_cast<py::ssize_t>(sizeof(Real))},
                owned->data(),
                free_when_done
            );
        }
        
        template <typename Scalar, typename Array>
        static std::shared_ptr<Scalar[]> toSharedPtr(const Array& arr) {
	        // Capture the Python object in the deleter so the buffer stays alive.
	        return std::shared_ptr<Scalar[]>(
	            static_cast<Scalar*>(arr.request().ptr),
	            [obj = py::object(arr)](Scalar*) mutable {
	                // Keep the py::object alive until the shared_ptr is destroyed.
	            }
	        );
	    }

	    template <typename Array>
	    static std::shared_ptr<float[]> toShared_ptr(const Array& arr) {
	        return toSharedPtr<float>(arr);
	    }


};

} // namespace mmcfilters
